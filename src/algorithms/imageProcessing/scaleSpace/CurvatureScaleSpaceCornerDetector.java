package algorithms.imageProcessing.scaleSpace;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.NearestPoints;
import algorithms.imageProcessing.CannyEdgeFilter;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.features.CornerRegion;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.MiscDebug;
import algorithms.util.PairIntArray;
import algorithms.util.PairInt;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

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
    
    /**
     * corners populated if extractSkyline is true
     */
    protected PairIntArray skylineCorners = new PairIntArray();

    protected boolean enableJaggedLineCorrections = true;
    
    protected float factorIncreaseForCurvatureMinimum = 1.f;
    
    protected boolean performWholeImageCorners = true;
    
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
    
    /**
     * set the edge detector to create edges that are better for outdoor
     * conditions and calculate corners only for the skyline.  
     * Note that the skyline extraction is currently
     * a long running process.
     */
    void calculateSkylineCornersOnly() {
        
        useOutdoorModeAndExtractSkyline();
        
        performWholeImageCorners = false;
    }
  
    public void enableJaggedLineCorrections() {
        enableJaggedLineCorrections = true;
    }
    public void disableJaggedLineCorrections() {
        enableJaggedLineCorrections = false;
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

        if (extractSkyline && !skylineEdges.isEmpty()) {
            calculateSkylineCorners();
        }

        if (!performWholeImageCorners) {
            return;
        }
             
        // not re-using return maps for now, but they are available here
        // while refactoring the public method returns and signatures
        Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > maps =
            findCornersInScaleSpaceMaps(edges, useOutdoorMode, corners);        
    }
    
    @Override
    protected void reinitializeSpecialization() {
        corners = new PairIntArray();
    }

    /**
     * find corners iteratively until approximately the number desired are
     * found.  This method is useful for creating corners in stereo projection
     * images - it helps to adjust the image intensity levels so that
     * similar edges can be formed in both images.
     * Presumably, if calibration of the images were possible, findCorners()
     * should alone provide a stable result (where calibration of the images
     * are steps such as bias removal, flat field corrections, intensity
     * corrections using standard candles, and corrections to make the
     * point spread functions similar for the images).
     *
     * @param approxNumberOfCornersDesired
     * @param filterOnlyAboveThisNumberOfCorners if the default number of corners
     * produced is this larger or larger, the method will iteratively
     * increase the lower intensity filter until approxNumberOfCornersDesired
     * are produced, else if the default number of corners is less
     * than useOnlyAboveThisNumberOfCorners, the method will not filter
     * the image further.
     */
    public void findCornersIteratively(int approxNumberOfCornersDesired,
        int filterOnlyAboveThisNumberOfCorners) {

        //TODO: this method needs to be refactored
        
        if (!performWholeImageCorners) {
            throw new IllegalStateException(
                "performWholeImageCorners is currently set to false");
        }

        float lowerThresholdStepSize = 1.0f;
        
        int nCorners = corners.getN();

        float lowThreshold = (useOutdoorMode) ?
            CannyEdgeFilter.defaultOutdoorLowThreshold :
            CannyEdgeFilter.defaultLowThreshold;
        
        if ((nCorners > 0) && (nCorners < filterOnlyAboveThisNumberOfCorners)) {
            return;
        } else if ((nCorners > 0) && (nCorners < approxNumberOfCornersDesired)) {
            return;
        } else if (state.ordinal() < CurvatureScaleSpaceMapperState.INITIALIZED.ordinal()) {
            findCorners();
        }

        nCorners = corners.getN();

        if (nCorners < filterOnlyAboveThisNumberOfCorners) {
            return;
        }

        //TODO: this could be improved to assert a minimum presence of corners
        // at boundaries and throughout the image compared to the first round.
        // In other words, would not want to remove all corners for an important
        // part of the image intersection with another image.

        List<PairIntArray> prevEdges = copy(this.edges);
        PairIntArray prevCorners = corners.copy();

        int nIter = 0;
            
        while ((nCorners > 0) && (nCorners > approxNumberOfCornersDesired)) {

            log.info("nCorners=" + nCorners);

            prevEdges = copy(this.edges);
            prevCorners = this.corners.copy();

            //TODO: needs adjustments:
            float additionalBlurSigma = 0;
            if (nIter == 0) {
                //additionalBlurSigma = 1;
                /*if (nCorners > 1000) {
                    lowerThresholdStepSize  = 1.5f;
                } else {*/
                    lowerThresholdStepSize = 0.5f;
                //}
            }
            
            lowThreshold += lowerThresholdStepSize;
            
            reinitialize(lowThreshold, additionalBlurSigma);

            findCornersInScaleSpaceMaps(edges, useOutdoorMode, corners);

            includeJunctionsInCorners();
            
            nCorners = corners.getN();
            
            nIter++;
        }

        if (Math.abs(corners.getN() - approxNumberOfCornersDesired) >
            Math.abs(prevCorners.getN() - approxNumberOfCornersDesired)) {
            
            this.edges = prevEdges;
            this.corners = prevCorners;
        }
        
    }

    /**
     * Find corners in the edges by creating scale space maps for each edge that
     * range from a maximum value determined in getMaxSIGMAForECSS() and
     * determine the corners in the highest sigma of those maps and refine the
     * corner locations in the smaller sigma maps.  The corners are found using
     * the curvature minima and maxima points in the curve.  A lower threshold
     * is determined and used during the maxima finding and minima comparison.
     * Each corner candidate is larger than one of the adjacent minima by
     * a factor such as 2 or 3.
     * The results are set in the instance variable corners.
     * The returned variable is the set of scale space maps which might be
     * useful for other purposes, but are no longer needed for the corner
     * determination.
     *
     * @param theEdges
     * @param doUseOutdoorMode
     * @param outputCorners
     * @return scale space maps for each edge
     */
    protected Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> >
    findCornersInScaleSpaceMaps(final List<PairIntArray> theEdges, 
        final boolean doUseOutdoorMode, final PairIntArray outputCorners) {

        CSSCornerMaker cornerMaker = new CSSCornerMaker(this.img.getWidth(),
            this.img.getHeight());
        
        if (enableJaggedLineCorrections) {
            cornerMaker.enableJaggedLineCorrections();
        }
        
        if (!doStoreCornerRegions) {
            cornerMaker.doNotStoreCornerRegions();
        }
        
        cornerMaker.increaseFactorForCurvatureMinimum(
            factorIncreaseForCurvatureMinimum);
        
        //map w/ key = edge index,
        //   values = junctions on edge as (x,y) and edge curve indexs
        Map<Integer, Set<PairIntWithIndex>> junctionCoords = getJunctionCoordinates();
        
        Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > sMap =
            cornerMaker.findCornersInScaleSpaceMaps(theEdges, junctionCoords,
                doUseOutdoorMode, outputCorners); 
        
        if (doStoreCornerRegions) {
            edgeCornerRegionMap.clear();
            Map<Integer, List<CornerRegion>> tmp = cornerMaker.getEdgeCornerRegionMap();
            edgeCornerRegionMap.putAll(tmp);
        }
    
        //TODO: place in the replacement for aspects:
        /*if (true) {
            try {
                Image imgCp = this.img.copyToColorGreyscale();
                ImageIOHelper.addAlternatingColorCurvesToImage(edges, imgCp, 3);
                MiscDebug.writeImage(imgCp, "_dbg_edges_");
                imgCp = imgCp = this.img.copyToColorGreyscale();
                for (Entry<Integer, List<CornerRegion>> entry : edgeCornerRegionMap.entrySet()) {
                    for (CornerRegion cr : entry.getValue()) {
                        int x = cr.getX()[cr.getKMaxIdx()];
                        int y = cr.getY()[cr.getKMaxIdx()];
                        ImageIOHelper.addPointToImage(x, y, imgCp, 2, 255, 0, 0);
                    }
                }
                MiscDebug.writeImage(imgCp, "_dbg_cornerregions_");
                imgCp = this.img.copyToColorGreyscale();
                for (int ii = 0; ii < outputCorners.getN(); ++ii) {
                    int x = outputCorners.getX(ii);
                    int y = outputCorners.getY(ii);
                    ImageIOHelper.addPointToImage(x, y, imgCp, 2, 255, 0, 0);
                }
                MiscDebug.writeImage(imgCp, "_dbg_corners_");
            } catch (IOException ex) {
                Logger.getLogger(this.class.getName()).log(Level.SEVERE, null, ex);
            }
        }*/
        
        return sMap;
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

    public PairIntArray getSkylineCornersInOriginalReferenceFrame() {

        PairIntArray co = new PairIntArray();

        for (int i = 0; i < skylineCorners.getN(); i++) {
            int x = skylineCorners.getX(i);
            int y = skylineCorners.getY(i);
            x += this.trimmedXOffset;
            y += this.trimmedYOffset;
            co.add(x, y);
        }

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

    public List<PairIntArray> getSkylineEdgesInOriginalReferenceFrame() {

        List<PairIntArray> output = new ArrayList<PairIntArray>();

        for (int i = 0; i < skylineEdges.size(); i++) {

            PairIntArray ce = new PairIntArray();

            PairIntArray edge = skylineEdges.get(i);

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

    //NOT READY FOR USE YET
    private void calculateSkylineCorners() {

        if (!extractSkyline) {
            return;
        }

        PairIntArray theSkylineCorners = new PairIntArray();
    
        enableJaggedLineCorrections = false;
        
        Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > maps =
            findCornersInScaleSpaceMaps(skylineEdges, false, theSkylineCorners);

        // TODO: improve this for resolution
        
        // smoothing does not make fewer points, so using min curvature factor
        
        if (filterProducts == null) {
            throw new IllegalStateException("Error in use of method.  filterProducts should not be null.");
        }
        
        float nGoalCorners = (filterProducts.getGradientXY().getNPixels() 
            < 500000) ? 50.f : 100.f;
        
        // want about 50 points
        float factor = (float)theSkylineCorners.getN()/nGoalCorners;
        
        if (factor > 1.2) {
            
            //TODO: consider saving the space curves to make recalc faster
            
            increaseFactorForCurvatureMinimum(4.0f * factor);
            
            PairIntArray theSkylineCorners2 = new PairIntArray();
            
            maps = findCornersInScaleSpaceMaps(skylineEdges, false, 
                theSkylineCorners2);
            
            log.info("before curvature factor change, number of skyline corners=" 
                + theSkylineCorners.getN());
            log.info("after=" + theSkylineCorners2.getN());
            
            theSkylineCorners = theSkylineCorners2;
            
            resetFactorForCurvatureMinimum();
        }
        
        if (theSkylineCorners.getN() > 0) {
            skylineCorners.addAll(theSkylineCorners);
        }
        
        enableJaggedLineCorrections = true;
        
        log.info("number of skyline corners=" + theSkylineCorners.getN());
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
