package algorithms.imageProcessing.scaleSpace;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.NearestPoints;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.features.CornerRegion;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.util.CornerArray;
import algorithms.util.PairIntArray;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * The code is implemented from interpreting several papers by the authors
 * Farzin Mokhtarian and Alan Mackworth.
 *
 * They prescribe a method for detecting features and corners that is scale
 * invariant, rotation and translation invariant and does not create
 * artifacts as a side effect.
 *
 * The method finds edges in an image and then calculates the curvature of
 * the edges to find "inflection" points along the curve.  Those points of
 * inflection as a function of scale parameter t are then findable in
 * another image that may have rotation or translation, for example, using
 * a search method such as A* to find the best matching features in t space.
 * The process of creating the scale based curves is repeated for increasing
 * sigma until no points are found with curvature = 0.
 *
 * The method uses the recursive and separable properties of Gaussians where
 * possible.  (NOTE, not finished implementing the recursive portion).
 *
 * @author nichole
 */
public class CurvatureScaleSpaceCornerDetector extends
    AbstractCurvatureScaleSpaceMapper {

    /**
     * corners detected in the image.  the false corners have been removed
     */
    protected PairIntArray corners = new PairIntArray();

    /**
     * this is not ready for use yet.  when implemented it should hold
     * a sublist of corners that are better to use for matching the same
     * in other images.
     * TODO: populate this with edgeCornerRegionMap values
     */
    protected PairIntArray cornersForMatching = new PairIntArray();

    protected Map<Integer, List<CornerRegion>> edgeCornerRegionMap = new
        HashMap<Integer, List<CornerRegion>>();
    
    protected float factorIncreaseForCurvatureMinimum = 1.f;
        
    protected boolean doStoreCornerRegions = true;
            
    public CurvatureScaleSpaceCornerDetector(final ImageExt input) {

        super(input);
    }
    
    public CurvatureScaleSpaceCornerDetector(final GreyscaleImage input) {

        super(input);
    }

    public CurvatureScaleSpaceCornerDetector(final ImageExt input,
        List<PairIntArray> theEdges) {

        super(input, theEdges);
    }
    
    public void doNotStoreCornerRegions() {
        doStoreCornerRegions = false;
    }
    
    public void increaseFactorForCurvatureMinimum(float factor) {
        factorIncreaseForCurvatureMinimum = factor;
    }
    
    public void resetFactorForCurvatureMinimum() {
        factorIncreaseForCurvatureMinimum = 1.f;
    }
    
    public void findCorners() {
        
        // extract edges and junction maps:
        initialize();
        
        /*
        Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > maps =
            findCornersInScaleSpaceMaps(edges, corners);
        */
        // changing to use hough transforms to clean the lines of corners 
        // created by line rendering
        Map<Integer,Set<PairIntWithIndex>> emptyJunctionsMap 
            = new HashMap<Integer,Set<PairIntWithIndex>>();
                
        CSSCornerMaker cornerMaker = new CSSCornerMaker(img.getWidth(), img.getHeight());
        cornerMaker.doNotStoreCornerRegions();
        cornerMaker.disableJaggedLineCorrections();
        List<CornerArray> cornerList =
            cornerMaker.findCornersInScaleSpaceMaps(edges, emptyJunctionsMap);
        
        /*
        filter out small curvature:
        see line 64 of comments in CornerRegion.java.
        for curvature smaller than 0.2 won't see changes in slope in the
             neighboring 2 points on either side.
        */
        Map<PairInt, Float> cornerMap = new HashMap<PairInt, Float>();
        for (int i = 0; i < cornerList.size(); ++i) {
            CornerArray ca = cornerList.get(i);
            for (int idx = 0; idx < ca.getN(); ++idx) {
                float curvature = ca.getCurvature(idx);
                
                if (Math.abs(curvature) > 0.05) {// should not set above 0.07...
                                        
                    cornerMap.put(new PairInt(Math.round(ca.getX(idx)),
                        Math.round(ca.getY(idx))), 
                        Float.valueOf(ca.getCurvature(idx)));
                }
            }
        }
        
        Set<PairInt> edgeJunctions = new HashSet<PairInt>();
        for (Entry<Integer, Set<Integer>> entry : junctionMap.entrySet()) {
            int pixIdx = entry.getKey().intValue();
            int col = img.getCol(pixIdx);
            int row = img.getRow(pixIdx);
            PairInt p = new PairInt(col, row);
            edgeJunctions.add(p);
        }
        
        if (this.filterProducts.getHoughLines() != null) {
            //remove intra-line corners from corner map
            cornerMaker.useHoughTransformationToFilterCornersForOrdered(edges, 
                cornerMap, edgeJunctions, filterProducts.getHoughLines(),
                img.getWidth(), img.getHeight());
        } else {
            //remove intra-line corners from corner map
            cornerMaker.useHoughTransformationToFilterCornersForOrdered(edges, 
                cornerMap, edgeJunctions, img.getWidth(), img.getHeight());
        }
        
        corners = new PairIntArray();
        for (Entry<PairInt, Float> entry : cornerMap.entrySet()) {
            int x = entry.getKey().getX();
            int y = entry.getKey().getY();
            corners.add(x, y);
        }
        
    }
    
    public PairIntArray getCorners() {
        return corners;
    }

    public PairIntArray getCornersForMatching() {
        return cornersForMatching;
    }

    /**
     * <em>note, this is not ready for use yet.</em>
     * it's meant to be a sublist of
     * the variable "corners" selected to be better for matching the same
     * corners in other images.
     * @return
     */
    public PairIntArray getCornersForMatchingInOriginalReferenceFrame() {

        PairIntArray co = new PairIntArray();
        /*for (int i = 0; i < cornersForMatching.getN(); i++) {
            int x = cornersForMatching.getX(i);
            int y = cornersForMatching.getY(i);
            x += this.trimmedXOffset;
            y += this.trimmedYOffset;
            co.add(x, y);
        }*/

        return co;
    }

    public PairIntArray getCornersInOriginalReferenceFrame() {
        
PairIntArray cp = corners.copy();
MultiArrayMergeSort.sortByYThenX(cp);
  
        PairIntArray co = new PairIntArray();
        for (int i = 0; i < corners.getN(); i++) {
            int x = corners.getX(i);
            int y = corners.getY(i);
            x += this.trimmedXOffset;
            y += this.trimmedYOffset;
            co.add(x, y);
        }

        return co;
    }

    public List<PairIntArray> getEdgesInOriginalReferenceFrame() {

        List<PairIntArray> output = new ArrayList<PairIntArray>();

        for (int i = 0; i < edges.size(); i++) {

            PairIntArray ce = new PairIntArray();

            PairIntArray edge = edges.get(i);

            for (int j = 0; j < edge.getN(); j++) {
                int x = edge.getX(j);
                int y = edge.getY(j);
                x += this.trimmedXOffset;
                y += this.trimmedYOffset;
                ce.add(x, y);
            }

            output.add(ce);
        }

        return output;
    }

    private List<PairIntArray> copy(List<PairIntArray> edges) {
        List<PairIntArray> copied = new ArrayList<PairIntArray>();
        for (PairIntArray edge : edges) {
            copied.add(edge.copy());
        }
        return copied;
    }

    public Map<Integer, List<CornerRegion>> getEdgeCornerRegionMap() {
        return edgeCornerRegionMap;
    }
    
    public Set<CornerRegion> getEdgeCornerRegions() {
        
        Set<CornerRegion> set = new HashSet<CornerRegion>();
        
        for (Entry<Integer, List<CornerRegion>> entry : edgeCornerRegionMap.entrySet()) {
            set.addAll(entry.getValue());
        }
        
        return set;
    }
    
     public Set<CornerRegion> getEdgeCornerRegionsInOriginalReferenceFrame() {
        
        Set<CornerRegion> set = new HashSet<CornerRegion>();
        
        for (Entry<Integer, List<CornerRegion>> entry : edgeCornerRegionMap.entrySet()) {
            for (CornerRegion cr : entry.getValue()) {
                CornerRegion crCopy = cr.copy();
                for (int i = 0; i < cr.getX().length; ++i) {
                    int x = cr.getX()[i] + this.trimmedXOffset;
                    int y = cr.getY()[i] + this.trimmedYOffset;
                    crCopy.set(i, cr.getK()[i], x, y);
                }
                set.add(crCopy);
            }
        }
        
        return set;
    }
     
    public Set<CornerRegion> getEdgeCornerRegionsInOriginalReferenceFrame(
        boolean removeAmbiguousPeaks) {
        
        if (!removeAmbiguousPeaks) {
            return getEdgeCornerRegionsInOriginalReferenceFrame();
        }

        Set<CornerRegion> set = getEdgeCornerRegions(removeAmbiguousPeaks);
       
        Set<CornerRegion> edited = new HashSet<CornerRegion>();
        
        for (CornerRegion cr : set) {
            CornerRegion crCopy = cr.copy();
            for (int i = 0; i < crCopy.getX().length; ++i) {
                int x = crCopy.getX()[i] + this.trimmedXOffset;
                int y = crCopy.getY()[i] + this.trimmedYOffset;
                crCopy.set(i, crCopy.getK()[i], x, y);
            }
            edited.add(crCopy);
        }
        
        return edited;
    }
    
    public Set<CornerRegion> getEdgeCornerRegions(boolean removeAmbiguousPeaks) {
        
        if (!removeAmbiguousPeaks) {
            return getEdgeCornerRegions();
        }
/*        
        Map<Integer, Set<Integer> > theJunctionMap = new HashMap<Integer, Set<Integer>>();

        Map<Integer, PairInt> theJunctionLocationMap = new HashMap<Integer, PairInt>();
        
        EdgeExtractorWithJunctions.findJunctions(edges, 
            theJunctionMap, theJunctionLocationMap, img.getWidth(), img.getHeight());

        this.junctionLocationMap = theJunctionLocationMap;
        this.junctionMap = theJunctionMap;
  */      
        Set<CornerRegion> set = new HashSet<CornerRegion>();
       
        for (Entry<Integer, List<CornerRegion>> entry : edgeCornerRegionMap.entrySet()) {
            set.addAll(entry.getValue());
        }
        
        return set;
    }

    /**
     * search for the point (x, y) within the junction maps, and if found,
     * construct best defining CornerRegion for the point.  The best defining
     * region for the curvature neighbors appears to be the lowest intensity 
     * pixels when the edge surrounds a void or dark patch.
     * This method may need to be revised or used with segmentation or 
     * object recognition (in which case it probably belongs in a class
     * to specialize handling of the junction points for the EdgeExtractorWithJunctions).
     * 
     * Note that the CornerRegions are not necessarily in the same frame
     * as the original image (they do not have trimmedXOffset and trimmedYOffset
     * added to them).
     * 
     * NOTE: there is a flaw here that this instance does not retain the
     * scale space curves or the curvature of all edges currently,
     * so 2 dummy values for the neighboring curvatures are returned
     * (a flag is set in the CornerRegion to indicate that).
     * 
     * @param x
     * @param y
     * @param k not used, just stored in the returned CornerRegion
     * @return 
     */
    protected List<CornerRegion> searchForCornerRegionWithinJunctions(int x, 
        int y, float k) {
        
        Integer pixIndex = Integer.valueOf(img.getIndex(x, y));

        // if the corner is a junction, these are it's neighbors:
        Set<Integer> junctionAdjPixIndexes = junctionMap.get(pixIndex);

        PairInt edgeLocation = junctionLocationMap.get(pixIndex);
                
        if ((junctionAdjPixIndexes == null || junctionAdjPixIndexes.isEmpty())
            && (edgeLocation == null)) {
            return null;
        }
        
        if (junctionAdjPixIndexes == null || junctionAdjPixIndexes.isEmpty()) {
            // because edgeLocation is not null, this pixel is in a junction.
            // iterate over each of the 8 neighbors and if it's in the
            // junctionLocationMap, add it to junctionAdjPixIndexes for
            // making CornerRegions
            if (junctionAdjPixIndexes == null) {
                junctionAdjPixIndexes = new HashSet<Integer>();
            }
            int[] dx8 = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
            int[] dy8 = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
            for (int i = 0; i < dx8.length; ++i) {
                int x1 = x + dx8[i];
                int y1 = y + dy8[i];
                if (x1 < 0 || y1 < 0 || (x1 > (img.getWidth() - 1)) || 
                (y1 > (img.getHeight() - 1))) {
                    continue;
                }
                Integer pixIndex1 = Integer.valueOf(img.getIndex(x1, y1));
                if (junctionLocationMap.containsKey(pixIndex1)) {
                    junctionAdjPixIndexes.add(pixIndex1);
                }
            }
        }
        
        List<CornerRegion> output = new ArrayList<CornerRegion>();
        
        List<Integer> list = new ArrayList<Integer>(junctionAdjPixIndexes);
        
        for (int i = 0; i < list.size(); ++i) {
            
            Integer pixIndex1 = list.get(i);
            
            PairInt edgeLocation1 = junctionLocationMap.get(pixIndex1);
            
            PairIntArray edge1 = edges.get(edgeLocation1.getX());
            
            //TODO: there's a bug in the junctionMap state!!!!  needs to be refactored
            int x1 = edge1.getX(edgeLocation1.getY());
            int y1 = edge1.getY(edgeLocation1.getY());
            
            for (int j = (i + 1); j < list.size(); ++j) {
                
                Integer pixIndex2 = list.get(j);
                                
                PairInt edgeLocation2 = junctionLocationMap.get(pixIndex2);
                
                PairIntArray edge2 = edges.get(edgeLocation2.getX());
                
                //TODO: there's a bug in the junctionMap state!!!!
                int x2 = edge2.getX(edgeLocation2.getY());
                int y2 = edge2.getY(edgeLocation2.getY());
            
                // making only 3 pixel regions instead of 5 for now
                
                CornerRegion cr = new CornerRegion(edgeLocation.getX(), 3, 1);
                cr.setFlagThatNeighborsHoldDummyValues();

                // add the image offsets
                cr.set(0, 0.1f*k, x1, y1);
                cr.set(1, k, x, y);
                cr.set(2, 0.1f*k, x2, y2);
                
                output.add(cr);
            }
        }
        
        return output;
    }
    
    /**
     * search for the point (x, y) within the junction maps and return true
     * if found.
     * 
     * @param x
     * @param y
     * @return 
     */
    protected boolean isInAJunction(int x, int y) {
        
        Integer pixIndex = Integer.valueOf(img.getIndex(x, y));

        Set<Integer> junctionAdjPixIndexes = junctionMap.get(pixIndex);

        if (junctionAdjPixIndexes != null && !junctionAdjPixIndexes.isEmpty()) {
            return true;
        }
        
        return junctionLocationMap.containsKey(pixIndex);
    }
    
    private void includeJunctionsInCorners() {

        int nTot = corners.getN();

        int[] x = new int[nTot];
        int[] y = new int[nTot];

        System.arraycopy(corners.getX(), 0, x, 0, corners.getN());
        System.arraycopy(corners.getY(), 0, y, 0, corners.getN());

        Set<PairInt> corners2 = new HashSet<PairInt>();
        for (int i = 0; i < corners.getN(); i++) {
            int x2 = corners.getX(i);
            int y2 = corners.getY(i);
            corners2.add(new PairInt(x2, y2));
        }

        //TODO: should be '3'?
        float sepDist = 4;

        final NearestPoints nearestPoints = new NearestPoints(x, y);

        for (Entry<Integer, Set<Integer>> entry : junctionMap.entrySet()) {

            int pixIdx = entry.getKey().intValue();

            int xP = img.getCol(pixIdx);
            int yP = img.getRow(pixIdx);

            Set<PairInt> overlapping = nearestPoints.findNeighbors(xP, yP,
                sepDist);

            for (PairInt p : overlapping) {
                corners2.remove(p);
            }

            corners2.add(new PairInt(xP, yP));
        }

        corners = new PairIntArray();

        for (PairInt p : corners2) {
            corners.add(p.getX(), p.getY());
        }

    }
    
    /**
     * extract from junction maps, a map w/ key = edge index,
     * values = junctions on edge as (x,y) and edge curve index.
     * @return 
     */
    private Map<Integer, Set<PairIntWithIndex>> getJunctionCoordinates() {
        
        Map<Integer, Set<PairIntWithIndex>> edgeCoordsMap = new HashMap<Integer, Set<PairIntWithIndex>>();
                
        for (Entry<Integer, Set<Integer>> entry : junctionMap.entrySet()) {

            int pixIdx = entry.getKey().intValue();

            int xP = img.getCol(pixIdx);
            int yP = img.getRow(pixIdx);
            
            PairInt edgeLocation = junctionLocationMap.get(Integer.valueOf(pixIdx));
            
            Integer edgeIndex = Integer.valueOf(edgeLocation.getX());
            
            PairIntWithIndex p = new PairIntWithIndex(xP, yP, edgeLocation.getY());
            
            Set<PairIntWithIndex> set = edgeCoordsMap.get(edgeIndex);
            if (set == null) {
                set = new HashSet<PairIntWithIndex>();
                edgeCoordsMap.put(edgeIndex, set);
            }
            
            set.add(p);
        }
        
        return edgeCoordsMap;
    }
}
