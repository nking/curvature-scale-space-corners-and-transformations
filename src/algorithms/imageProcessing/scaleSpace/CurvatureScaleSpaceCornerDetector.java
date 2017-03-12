package algorithms.imageProcessing.scaleSpace;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.NearestPoints;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
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
                
        CSSCornerMaker cornerMaker = new CSSCornerMaker(img.getWidth(), img.getHeight());
        cornerMaker.doNotStoreCornerRegions();
        List<CornerArray> cornerList =
            cornerMaker.findCornersInScaleSpaceMaps(edges);
        
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
    
}
