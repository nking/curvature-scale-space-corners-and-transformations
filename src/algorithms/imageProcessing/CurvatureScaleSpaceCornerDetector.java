package algorithms.imageProcessing;

import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;

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
 * possible.
 * 
 * @author nichole
 */
public class CurvatureScaleSpaceCornerDetector extends 
    AbstractCurvatureScaleSpaceMapper {
    
    protected PairIntArray corners = new PairIntArray();
    
    protected boolean doNotUseNoisyEdgeCorners = false;    
    
    public CurvatureScaleSpaceCornerDetector(final GreyscaleImage input) {
                
        super(input);
    }
    
    /**
     * constructor with option needed for images that are line drawings or
     * are solid blocks of color.
     * @param input
     * @param doUseLineMode 
     */
    public CurvatureScaleSpaceCornerDetector(final GreyscaleImage input,
        boolean doUseLineMode) {
                        
        super(input, true);
    }
    
    public CurvatureScaleSpaceCornerDetector(final GreyscaleImage input, 
        List<PairIntArray> theEdges) {
        
        super(input, theEdges);
    }
     
    /**
     * when constructing the corner list, exclude corners from the edges which
     * are very noisy.  NOTE: This feature is not yet implemented, so has no effect
     * on the results.
     */
    public void doNotUseNoisyEdgeCorners() {
        doNotUseNoisyEdgeCorners = true;
    }
    
    public void findCorners() {
         
        initialize();
               
        // not re-using return maps for now, but they are available here
        // while refactoring the public method returns and signatures
        Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > maps = 
            findCornersInScaleSpaceMaps();
        
    }

    /**
     * maxSigma is defined by the ECSS algorithm in:
     * 2006, "Performance evaluation of corner detectors using consistency and 
     * accuracy measures" by Farzin Mokhtarian and Farahnaz Mohanna in
     * Computer Vision and Image Understanding, vol 102, pp 81-94.
     * @param nPoints
     * @return 
     */
    private SIGMA getMaxSIGMAForECSS(int nPoints) {
        
        //ECSS:
        //    < 200, sigma = 2
        //    <= 600, sigma=3
        //    else    sigma=4
        if (nPoints < 200) {
            return SIGMA.TWO;
        } else if (nPoints < 601) {
            return SIGMA.THREE;
        } else {
            return SIGMA.FOUR;
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
     * @return scale space maps for each edge
     */
    protected Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > 
    findCornersInScaleSpaceMaps() {
        
        Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > scaleSpaceMaps 
            = new HashMap<PairIntArray, Map<SIGMA, ScaleSpaceCurve> >();
      
        // perform the analysis on the edgesScaleSpaceMaps
        for (int i = 0; i < edges.size(); i++) {
  
            if (doNotUseNoisyEdgeCorners) {
                if (highChangeEdges.contains(Integer.valueOf(i))) {
                    continue;
                }
            }
            
            final PairIntArray edge = edges.get(i);
            
            final Map<SIGMA, ScaleSpaceCurve> map = 
                createLowUpperThresholdScaleSpaceMaps(edge);
            
            scaleSpaceMaps.put(edge, map);

            PairIntArray edgeCorners = findCornersInScaleSpaceMap(edge, map, i);
            
            log.log(Level.FINE, 
                "{0}) number of corners adding ={1} for edge={2}", 
                new Object[]{Integer.valueOf(i), 
                    Integer.valueOf(edgeCorners.getN()), 
                    Integer.valueOf(i)}); 
                   
//edgeCorners = removeFalseCorners(edgeCorners);
            
            //store xc and yc for the edge
            for (int ii = 0; ii < edgeCorners.getN(); ii++) {
                corners.add(edgeCorners.getX(ii), edgeCorners.getY(ii));
            }
        }
        
        return scaleSpaceMaps;
    }

    /**
     * find the corners in the given scale space map.  
     * The corners are found using the curvature minima and maxima points in 
     * the curve.  A lower threshold is determined and used during the maxima 
     * finding and minima comparison.  Each corner candidate is larger than one 
     * of the adjacent minima by a factor such as 2 or 3.
     * 
     * The returned variable is the set of (x,y) points for the candidate
     * corners found.
     * 
     * @param scaleSpace scale space map for an edge
     * @param edgeNumber the edgeNumber of the scaleSpace.  it's passed for
     * debugging purposes.
     * @param correctForJaggedLines
     * @return 
     */
    protected PairFloatArray findCornersInScaleSpaceMap(
        final ScaleSpaceCurve scaleSpace, int edgeNumber, 
        boolean correctForJaggedLines) {
        
        float[] k = Arrays.copyOf(scaleSpace.getK(), scaleSpace.getK().length);
        
        float[] outputLowThreshold = new float[1];
        
        List<Integer> minimaAndMaximaIndexes = findMinimaAndMaximaInCurvature(
            k, outputLowThreshold);

        List<Integer> maxCandidateCornerIndexes = findCandidateCornerIndexes(
            k, minimaAndMaximaIndexes, outputLowThreshold[0]);
    
        PairFloatArray xy = new PairFloatArray(maxCandidateCornerIndexes.size());
        
        int nRemoved = 0;
        
        if (correctForJaggedLines) {
            List<Integer> remove = new ArrayList<Integer>();
            for (int ii = 0; ii < maxCandidateCornerIndexes.size(); ii++) {
                int idx = maxCandidateCornerIndexes.get(ii);
                if ((idx > 1) && (idx < (scaleSpace.getSize() - 2))) {
                    
                    boolean isDueToJaggedLine = isDueToJaggedLine(idx, scaleSpace);
                    if (!isDueToJaggedLine) {
                        xy.add(scaleSpace.getX(idx), scaleSpace.getY(idx));
                    } else {
                        remove.add(Integer.valueOf(ii));
                    }
                } else {
                    xy.add(scaleSpace.getX(idx), scaleSpace.getY(idx));
                }
            }
            nRemoved += remove.size();
            for (int ii = (remove.size() - 1); ii > -1; ii--) {
                int idx = remove.get(ii).intValue();
                maxCandidateCornerIndexes.remove(idx);
            }
        } else {
            for (int ii = 0; ii < maxCandidateCornerIndexes.size(); ii++) {
                int idx = maxCandidateCornerIndexes.get(ii);
                xy.add(scaleSpace.getX(idx), scaleSpace.getY(idx));
            }
        }

        log.info("NREMOVED=" + nRemoved);
        
        return xy;
    }
    
    /**
     * determine the corners in the highest sigma of those maps and refine the
     * corner locations in the smaller sigma maps.
     * 
     * @param edge
     * @param scaleSpaceCurves
     * @param edgeNumber
     * @return 
     */
    private PairIntArray findCornersInScaleSpaceMap(final PairIntArray edge, 
        final Map<SIGMA, ScaleSpaceCurve> scaleSpaceCurves, int edgeNumber) {
       
        SIGMA maxSIGMA = getMaxSIGMAForECSS(edge.getN());         
           
        ScaleSpaceCurve maxScaleSpace = scaleSpaceCurves.get(maxSIGMA);
           
        if (maxScaleSpace == null) {
            throw new IllegalStateException(
                "could not find the scale space for max sigma");
        }
        
        if (maxScaleSpace.getK() == null) {
            return new PairIntArray(0);
        }
 
        PairFloatArray candidateCornersXY = 
            findCornersInScaleSpaceMap(maxScaleSpace, edgeNumber, true);
        
        SIGMA sigma = SIGMA.divideBySQRT2(maxSIGMA);  

        SIGMA prevSigma = maxSIGMA;
                
        //find the corners in the higher res scale space curves
        while (sigma != null) {

            ScaleSpaceCurve scaleSpace = scaleSpaceCurves.get(sigma);

            refinePrimaryCoordinates(candidateCornersXY.getX(), 
                candidateCornersXY.getY(), scaleSpace, prevSigma, 
                edgeNumber);
            
            prevSigma = sigma;

            sigma = SIGMA.divideBySQRT2(sigma);
        }
        
        log.log(Level.FINE, "number of corners adding ={0}", 
            Integer.valueOf(candidateCornersXY.getN()));
        
        PairIntArray edgeCorners = new PairIntArray(candidateCornersXY.getN());
                                 
        //store xc and yc for the edge
        for (int ii = 0; ii < candidateCornersXY.getN(); ii++) {
            
            edgeCorners.add((int)candidateCornersXY.getX(ii), 
                (int)candidateCornersXY.getY(ii));
        }
        
        return edgeCorners;
    }
    
    /**
     * Construct scale space images (that is X(t, sigma), y(t, sigma), and
     * k(t, sigma) where t is the intervals of spacing along the curve
     * and is valued from 0 to 1 and k is the curvature).
     * 
     * The range of the scale space maps is from sigma = 0.5 up to a maximum
     * value determined in getMaxSIGMAForECSS.
     * 
     * The results are returned as a map keyed by sigma.
     */
    private Map<SIGMA, ScaleSpaceCurve> createLowUpperThresholdScaleSpaceMaps(
        final PairIntArray edge) {
        
        ScaleSpaceCurvature scaleSpaceHelper = new ScaleSpaceCurvature();
        
        Map<SIGMA, ScaleSpaceCurve> map = new HashMap<SIGMA, ScaleSpaceCurve>();
                    
        SIGMA sigma = SIGMA.ZEROPOINTFIVE;

        //TODO: take a look at the highest resolution curvature array
        //   and see if can derive the same upper limit by S/N estimate
        //   and assuming want to smooth to 2 o3 3 times the noise
        SIGMA maxSIGMA = getMaxSIGMAForECSS(edge.getN());
            
        // this increases by a factor of sqrt(2)
        float resultingSigma = SIGMA.getValue(sigma);

        ScaleSpaceCurve lastCurve = null;

        while (sigma.compareTo(maxSIGMA) < 1) {

            ScaleSpaceCurve curve;
                
            if (lastCurve == null) {
                curve = scaleSpaceHelper.computeCurvature(edge, sigma, 
                    resultingSigma);
            } else {
                curve = scaleSpaceHelper.computeCurvature(
                    lastCurve.getXYCurve(), sigma, resultingSigma);
            }

            map.put(sigma, curve);

            sigma = SIGMA.increaseToFactorBySQRT2(resultingSigma);
                    
            resultingSigma *= Math.sqrt(2);
              
            lastCurve = curve;
                        
        }
        
        if (!map.containsKey(maxSIGMA)) {
            
            ScaleSpaceCurve curve = scaleSpaceHelper.computeCurvature(edge, 
                maxSIGMA, SIGMA.getValue(maxSIGMA));
            
            map.put(maxSIGMA, curve);
        }
        
        return map;
    }
   
    /**
     * given histogram h, find the first peak and the subsequent first minima
     * after it.
     * 
     * @param h
     * @return 
     */
    protected int[] findFirstPeakAndMinimum(HistogramHolder h) {
        
        float[] xh = h.getXHist();
        float[] yh = h.getYHistFloat();
        
        int firstMinIdx = -1;
        int yFirstPeakIdx = -1;
        float lastY = yh[0];
        boolean incr = true;
        for (int i = 1; i < xh.length; i++) {
            if (yFirstPeakIdx > -1) {
                if ((yh[i] < lastY) && incr) {
                    incr = false;
                } else if ((yh[i] > lastY) && !incr) {
                    firstMinIdx = i - 1;
                    break;
                }
            } else {
                if ((yh[i] < lastY) && incr) {
                    incr = false;
                    yFirstPeakIdx = i - 1;
                }
            }
            lastY = yh[i];
        }
        firstMinIdx++;
        
        return new int[]{yFirstPeakIdx, firstMinIdx};
    }
    
    /**
     * find the minima and maxima of the curvature k and return the low threshold
     * used in the in-out variable outputLowThreshold and a list of the
     * indexes for the minima as negative values and the maxima as positive 
     * values.
     * @param k
     * @param outputLowThreshold
     * @return 
     */
    protected List<Integer> findMinimaAndMaximaInCurvature(float[] k,
        float[] outputLowThreshold) {
        
        if ((k == null) || (k.length == 0)) {
            return new ArrayList<Integer>();
        }
        
        for (int ii = 1; ii < k.length; ii++) {
            if (k[ii] < 0.0f) {
                k[ii] *= -1.0f;
            }
        }
        
        float kMax = MiscMath.findMax(k);
        
        // determine float lowThresh
        HistogramHolder h = Histogram.calculateSturgesHistogramRemoveZeroTail(
            k, Errors.populateYErrorsBySqrt(k));
        
        int[] firstPeakAndMinIdx = findFirstPeakAndMinimum(h);
       
        if (firstPeakAndMinIdx[1] > 3) {
            firstPeakAndMinIdx[1] = firstPeakAndMinIdx[0] + 1;
        }
        if (firstPeakAndMinIdx[1] >= (h.getYHist().length >> 1)) {
            firstPeakAndMinIdx[1] = 1;
        }
        
        if (kMax/(h.getXHist()[firstPeakAndMinIdx[1]]) > 10000) {
            
            h = Histogram.calculateSturgesHistogram(
                0, kMax, k, Errors.populateYErrorsBySqrt(k));
            
            firstPeakAndMinIdx = findFirstPeakAndMinimum(h);
        }
        
        // sum intensity <= firstMinIdx and then after to compare
        long sum0 = 0;
        long sum1 = 0;
        for (int i = 0; i < h.getXHist().length; i++) {
            if (i <= firstPeakAndMinIdx[1]) {
                sum0 += h.getYHist()[i];
            } else {
                sum1 += h.getYHist()[i];
            }
        }
        
        float divSum = (float)sum0/(float)sum1;
        outputLowThreshold[0] = 0;
        if (divSum > 10) {
            outputLowThreshold[0] = h.getXHist()[firstPeakAndMinIdx[0]];
        } else if ((firstPeakAndMinIdx[1] > 0) && (sum1 > 0)) {
            outputLowThreshold[0] = (h.getXHist()[firstPeakAndMinIdx[1]] 
                + h.getXHist()[firstPeakAndMinIdx[1] - 1])/2;
        } else if ((firstPeakAndMinIdx[1] > 0) && (sum1 == 0)) {
            firstPeakAndMinIdx[1]--;
        }
                
        log.fine("lowThresh=" + outputLowThreshold[0] 
            + " sum0=" + sum0 + " sum1=" + sum1 + " divsum=" + divSum
            + " firstPeakAndMinIdx[0]=" + firstPeakAndMinIdx[0] 
            + " firstPeakAndMinIdx[1]=" + firstPeakAndMinIdx[1]);

        /*
        storing the minima and maxima in the same array list.
        the minima have -1*index within k
        and the maxima keep their positive values of the index within k.
        */
        List<Integer> minMaxIndexes = new ArrayList<Integer>();
        
        float lastK = k[0];
        int lastMinIdx = -1;
        boolean incr = true;
        for (int ii = 1; ii < k.length; ii++) {
            
            if ((k[ii] < lastK) && incr) {
                if (k[ii - 1] > outputLowThreshold[0]) {
                    minMaxIndexes.add(Integer.valueOf(ii - 1));
                }
                incr = false;
            } else if ((k[ii] > lastK) && !incr) {
                // values below outputLowThreshold[0] are handled by 
                // callers.  TODO: redesign the caller and this method
                // to not need to understand peculiarities of the data.
                minMaxIndexes.add(Integer.valueOf(-1*(ii - 1)));
                incr = true;
            }

            lastK = k[ii];
        }
        
        if (incr) {
            // add the last point
             minMaxIndexes.add(Integer.valueOf(k.length - 1));
        }
        
        return minMaxIndexes; 
    }
  
    /**
     * given the curvature array and a list of the indexes of the
     * minima and maxima of the curvature array, find the candidate corner
     * indexes with respect to the k array.
     * 
     * @param k
     * @param minMaxIndexes
     * @param lowThreshold
     * @return 
     */
    protected List<Integer> findCandidateCornerIndexes(float[] k, 
        List<Integer> minMaxIndexes, float lowThreshold) {
        
        // find peaks where k[ii] is > factorAboveMin* adjacent local minima 
        
        float factorAboveMin = 3.0f; // 10 misses some corners

        List<Integer> cornerCandidates = new ArrayList<Integer>();

        // choose candidates from minMaxIndexes that are 
        //     >= factorAboveMin one adjacent min
        for (int ii = 0; ii < minMaxIndexes.size(); ii++) {

            int idx = minMaxIndexes.get(ii).intValue();

            if (idx > -1) {
                // this is maximum
                
                boolean found = false;
                
                // compare to preceding minimum
                for (int iii = (ii - 1); iii > -1; iii--) {
                    int idx2 = minMaxIndexes.get(iii).intValue();
                    if (idx2 < 0) {
                        float compare = k[-1*idx2];
                        if (compare < lowThreshold) {
                            // avoids divide by very small number sometimes
                            compare = lowThreshold;
                        }
                        if (k[idx] >= factorAboveMin * compare) {
                            cornerCandidates.add(Integer.valueOf(idx));
                            found = true;
                        }
                        break;
                    }
                }
                if (found) {
                    continue;
                }

                //compare to proceeding minimum
                for (int iii = (ii + 1); iii < minMaxIndexes.size(); iii++) {
                    int idx2 = minMaxIndexes.get(iii).intValue();
                    if (idx2 < 0) {
                        float compare = k[-1*idx2];
                        if (compare < lowThreshold) {
                            // avoids divide by very small number sometimes
                            compare = lowThreshold;
                        }
                        if (k[idx] >= factorAboveMin * compare) {
                            cornerCandidates.add(Integer.valueOf(idx));
                        }
                        break;
                    }
                }
            }
        }
        
        return cornerCandidates;
    }
   
    /**
     * refine the primary coordinates given in xc and yc using the corner
     * candidate results of the given scale space map.  Matches to a corner
     * are found using a rough separation limit determined by the previous
     * scale map's sigma. 
     * 
     * @param xc
     * @param yc
     * @param scaleSpace
     * @param previousSigma
     * @param edgeNumber included for debugging
     */
    private void refinePrimaryCoordinates(float[] xc, float[] yc, 
        ScaleSpaceCurve scaleSpace, SIGMA previousSigma,
        int edgeNumber) {
        
        if (scaleSpace == null || scaleSpace.getK() == null) {
            //TODO: follow up on NPE here
            return;
        }
        
        PairFloatArray xy2 = 
            findCornersInScaleSpaceMap(scaleSpace, edgeNumber, false);
       
        // roughly estimating maxSep as the ~FWZI of the gaussian
        //TODO: this may need to be altered to a smaller value
        float maxSepSq = 2*Gaussian1D.estimateHWZI(previousSigma, 0.01f);
        maxSepSq *= maxSepSq;
            
        // revise the points in {xc, yc} to the closest in {xc2, yc2}
        boolean[] matchedNew = new boolean[xy2.getN()];
        for (int j = 0; j < xc.length; j++) {
            float x = xc[j];
            float y = yc[j];
            float minSepSq = Float.MAX_VALUE;
            int minSepIdx = -1;
            for (int jj = 0; jj < xy2.getN(); jj++) {
                if (matchedNew[jj]) {
                    continue;
                }
                float x2 = xy2.getX(jj);
                float y2 = xy2.getY(jj);
                float dx = x2 - x;
                float dy = y2 - y;
                float sepSq = (dx*dx) + (dy*dy);
                if (sepSq < minSepSq) {
                    minSepSq = sepSq;
                    minSepIdx = jj;
                }
            }
            
            if (minSepIdx > -1) {
                if (minSepSq < maxSepSq) {
                    matchedNew[minSepIdx] = true;
                    xc[j] = xy2.getX(minSepIdx);
                    yc[j] = xy2.getY(minSepIdx);
                }
            }
        }

    }

    public PairIntArray getCorners() {
        return corners;
    }
    
    public PairIntArray getCornersInOriginalReferenceFrame() {
        
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

    private PairIntArray removeFalseCorners(PairIntArray edgeCorners) {
        
        if (edgeCorners.getN() < 3) {
            return edgeCorners;
        } 
        
        float tolerance = 0.5f;//0.2f;
        
        PairIntArray output = new PairIntArray();
        
        /*
        -- for each point in edge,
           -- calculate the slope to the next corner.
           -- for the points where the slope is the same within some tolerance,
              remove all points in between.
        */
        float[] slopes = new float[edgeCorners.getN() - 1];
        
        for (int i = 0; i < (edgeCorners.getN() - 1); i++) {
            
            float diffX = edgeCorners.getX(i + 1) - edgeCorners.getX(i);
            float diffY = edgeCorners.getY(i + 1) - edgeCorners.getY(i);
            
            slopes[i] = (diffX == 0) ? 1 : (diffY/diffX);
            
            log.fine("(" + edgeCorners.getX(i) + ", " + edgeCorners.getY(i) +
                ")  slope=" + slopes[i]); 
        }
        
        boolean[] remove = new boolean[slopes.length];
        
        int startSameIdx = -1;
        float lastSlope = Integer.MAX_VALUE;
        for (int i = 0; i < slopes.length; i++) {
            
            float slope = slopes[i];
            
            boolean withinTol = false;
            
            if (Math.abs(slope - lastSlope) < tolerance) {
                if (startSameIdx == -1) {
                    startSameIdx = i;
                }
                withinTol = true;
            } 
            
            if (!withinTol || (i == (slopes.length - 1))) {
                
                if (startSameIdx > -1) {
                    
                    //TODO: may need to estimate mean slope of removable section
                    // and stdev and only remove if stdev is small
                    
                    // mark the points in between for removal
                    // startSameIdx __  i
                    //   0           1  2
                    int nRemovable = i - startSameIdx - 1;
                    if (nRemovable > 0) {
                        for (int j = (startSameIdx + 1); j < i; j++) {
                            remove[j] = true;
                        }
                    }
                    
                    startSameIdx = -1;
                }
            }
            
            lastSlope = slope;
        }
        
        for (int i = 0; i < remove.length; i++) {
            if (!remove[i]) {
                output.add(edgeCorners.getX(i), edgeCorners.getY(i));
            }
        }
        
        output.add(edgeCorners.getX(edgeCorners.getN() - 1), 
            edgeCorners.getY(edgeCorners.getN() - 1));
        
        return output;
    }

    private boolean isDueToJaggedLine(int idx, ScaleSpaceCurve scaleSpace) {
        
        // range check for idx larger and smaller than image boundaries
        // by 2 pixels has been performed by invoker
        
        float prev2X = scaleSpace.getX(idx - 2);
        float prev2Y = scaleSpace.getY(idx - 2);
        float prev1X = scaleSpace.getX(idx - 1);
        float prev1Y = scaleSpace.getY(idx - 1);
        float x = scaleSpace.getX(idx);
        float y = scaleSpace.getY(idx);
        float next1X = scaleSpace.getX(idx + 1);
        float next1Y = scaleSpace.getY(idx + 1);
        float next2X = scaleSpace.getX(idx + 2);
        float next2Y = scaleSpace.getY(idx + 2);
        float k = scaleSpace.getK(idx);
        
        if ((x > 80) && (x < 90) && (y > 140) && (y < 150)) {
            int z = 1;
        }
        
        // check for vertical lines,
        if (prev2X == prev1X) {
            if (prev2Y == (prev1Y + 1)) {
                /* check for vertical lines, increasing in Y with idx
                0
                1   {1}
                2   {2}
                  @  @  @
                 {4}    {4}
                 {5}    {5}
                */
                if (next2X == next1X) {
                    if ((next2X == (prev1X + 1)) ||(next2X == prev1X) ||
                        (next2X == (prev1X - 1))) {
                        if ((x == prev1X) || (x == (prev1X - 1))
                            || (x == (prev1X + 1))) {

                            if ((y == (prev1Y - 1)) && (next1Y == (y - 1)) &&
                                (next2Y == (next1Y - 1))) {

                                // skip this corner
                                return true;
                            }
                        }
                    }
                }
            } else if (prev2Y == (prev1Y - 1)) {
                /* check for vertical lines, decreasing in Y with idx
                0
                1 {5}  {5}
                2 {4}  {4}
                  @  @  @
                    {2}
                    {1}
                */
                if (next2X == next1X) {
                    if ((next2X == (prev1X + 1)) || (next2X == prev1X) ||
                        (next2X == (prev1X - 1))) {
                        if ((x == prev1X) || (x == (prev1X - 1))
                            || (x == (prev1X + 1))) {

                            if ((y == (prev1Y + 1)) && (next1Y == (y + 1)) &&
                                (next2Y == (next1Y + 1))) {

                                // skip this corner
                                return true;
                            }
                        }
                    }
                }
            } 
        } else if (prev2Y == prev1Y) {
            if (prev2X == (prev1X + 1)) {
                /* check for horizontal lines, decreasing in X with idx
                0
                1 {5} {4} @
                2         @ {2} {1}
                3 {5} {4} @
                */
                if (next2Y == next1Y) {
                    if ((next2Y == (prev1Y + 1)) || (next2Y == prev1Y) ||
                        (next2Y == (prev1Y - 1)) ) {
                        if ((y == prev1Y) || (y == (prev1Y - 1))
                            || (y == (prev1Y + 1))) {

                            if ((x == (prev1X - 1)) && (next1X == (x - 1)) 
                                && (next2X == (next1X - 1))) {
                                
                                /*
                                check one more pixel on the "next" end
                                */
                                if ((idx + 3) < scaleSpace.getSize()) {
                                    float next3X = scaleSpace.getX(idx + 3);
                                    float next3Y = scaleSpace.getY(idx + 3);
                                    if ((next3Y == next2Y) && (next3X == (next2X - 1))) {
                                        /*
                                        check one more pixel on the "prev" end
                                        */
                                        if (idx > 2) {
                                            float prev3X = scaleSpace.getX(idx - 3);
                                            float prev3Y = scaleSpace.getY(idx - 3);
                                            if ((prev3Y == prev2Y) && 
                                                (prev3X == (prev2X + 1))) {
                                                // skip this corner
                                                return true;
                                            }
                                        } else {
                                            return true;
                                        }
                                    }
                                } else {
                                    /*
                                    check one more pixel on the "prev" end
                                    */
                                    if (idx > 2) {
                                        float prev3X = scaleSpace.getX(idx - 3);
                                        float prev3Y = scaleSpace.getY(idx - 3);
                                        if ((prev3Y == prev2Y)
                                            && (prev3X == (prev2X + 1))) {
                                            // skip this corner
                                            return true;
                                        }
                                    } else {
                                        return true;
                                    }
                                    //TODO: make same checks for 3rd pixel in other blocks
                                }
                            }
                        }
                    }
                }
            } else if (prev2X == (prev1X + 1)) {
                /* check for horizontal lines, increasing in X with idx
                0
                1         @ {4} {5}
                2 {1} {2} @
                3         @ {4} {5}
                */
                if (next2Y == next1Y) {
                    if ((next2Y == (prev1Y + 1)) || (next2Y == prev1Y) ||
                        (next2Y == (prev1Y - 1)) ) {
                        if ((y == prev1Y) || (y == (prev1Y - 1))
                            || (y == (prev1Y + 1))) {

                            if ((x == (prev1X + 1)) && (next1X == (x + 1)) 
                                && (next2X == (next1X + 1))) {

                                // skip this corner
                                return true;
                            }
                        }
                    }
                }
            }
        }
        
        // look for the prev and next pairs to be a line
        float dPrevX = prev2X - prev1X;
        float dPrevY = prev2Y - prev1Y;
        float dNextX = next2X - next1X;
        float dNextY = next2Y - next1Y;
        float prevSlope = (dPrevX == 0) ? 1 : dPrevY/dPrevX;
        float nextSlope = (dNextX == 0) ? 1 : dNextY/dNextX;
        
        float dPrevNextX = prev1X - next1X;  
        float dPrevNextY = prev1Y - next1Y;
        
        if (prevSlope == nextSlope) {
            if (prevSlope == -1) {
                if ((dPrevNextX == -2) && (dPrevNextY == 1)) {
                    return true;
                } else if ((dPrevNextX == -1) && (dPrevNextY == 2)) {
                    return true;
                } else if ((dPrevNextX == 1) && (dPrevNextY == -2)) {
                    return true;
                } else if ((dPrevNextX == 2) && (dPrevNextY == -1)) {
                    return true;
                } else if ((dPrevNextX == 2) && (dPrevNextY == 0)) {
                    return true;
                }
            } else if (prevSlope == 1) {
                if ((dPrevNextX == 2) && (dPrevNextY == 0)) {
                    return true;
                } else if ((dPrevNextX == 2) && (dPrevNextY == -1)) {
                    return true;
                }
            } else if (((prevSlope == 0) || (prevSlope == -0)) 
                && ((dPrevNextX == 2) || (dPrevNextX == -2)) && (dPrevNextY == -2)
            ) {
                return true;
            } 
        }
        
        return false;
    }

    private boolean isAGoodMatchingCandidate(int idx, ScaleSpaceCurve scaleSpace) {
       
        float prev2X = scaleSpace.getX(idx - 2);
        float prev2Y = scaleSpace.getY(idx - 2);
        float prev1X = scaleSpace.getX(idx - 1);
        float prev1Y = scaleSpace.getY(idx - 1);
        float x = scaleSpace.getX(idx);
        float y = scaleSpace.getY(idx);
        float next1X = scaleSpace.getX(idx + 1);
        float next1Y = scaleSpace.getY(idx + 1);
        float next2X = scaleSpace.getX(idx + 2);
        float next2Y = scaleSpace.getY(idx + 2);
        float dPrevX = prev2X - prev1X;
        float dPrevY = prev2Y - prev1Y;
        float dNextX = next2X - next1X;
        float dNextY = next2Y - next1Y;
        float prevSlope = (dPrevX == 0) ? 1 : dPrevY / dPrevX;
        float nextSlope = (dNextX == 0) ? 1 : dNextY / dNextX;

        float dPrevNextX = prev1X - next1X;
        float dPrevNextY = prev1Y - next1Y;

        if (Math.abs(Math.abs(prevSlope) - Math.abs(nextSlope)) == 1) {
            if ((prevSlope == 0) || (nextSlope == 0)) {
                return true;
            }
        }
        
        return false;
    }

}
