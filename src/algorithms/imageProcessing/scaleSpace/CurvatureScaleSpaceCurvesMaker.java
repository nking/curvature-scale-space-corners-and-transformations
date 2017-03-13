package algorithms.imageProcessing.scaleSpace;

import algorithms.imageProcessing.SIGMA;
import algorithms.util.PairIntArray;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class CurvatureScaleSpaceCurvesMaker {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * Create for an edge, 
     * X(t,sigma), Y(t,sigma), k(t, sigma) for sigma=1 until sigma=the level at 
     * which there are no more values of k that equal 0, that is until no 
     * more inflection points are found on the extremely smoothed curve.  
     * 
     * Note that the method re-uses iterative convolution, so each interval
     * of sigma is 2^(1/8) times the previous convolution result using a kernel
     * of sigma = 2^(1/8) each time.  The error in determining the peak of
     * a contour in the resulting scale space curve should be less than 10%.
     * TODO: calculate error in peak determination....
     * 
     * Note that the method follows the Mokhtarian and Mackworth papers 
     * referenced in the other classes.
     * 
     * @param edge
     * @param sigmaPowerFactor power to use for base 2 in creating the next 
     * larger kernel convolution.  a coarse grained set of intervals is (1./8.) 
     * while a finer grained set of intervals is (1./32.)
     * @param sigmaStart the starting sigma of the gaussian convolution operations
     * @param sigmaEnd the last sigma of the gaussian convolution operations, 
     * exclusive.
     * @return 
     */
    public Map<Float, ScaleSpaceCurve> createScaleSpaceMetricsForEdge(
    PairIntArray edge, float sigmaPowerFactor, SIGMA sigmaStart, SIGMA sigmaEnd) {
        
        ScaleSpaceCurvature scaleSpaceHelper = new ScaleSpaceCurvature();
            
        Map<Float, ScaleSpaceCurve> scaleSpaceMap = new HashMap<Float,
            ScaleSpaceCurve>();
        
        float sigma = SIGMA.getValue(sigmaStart);

        float resultingSigma = sigma;

        boolean hasInflectionPoints = true;

        ScaleSpaceCurve lastCurve = null;

        /*TODO: implement the recursive methods for gaussian derivatives for
        these using their separable properties.
        https://hal.inria.fr/inria-00074778/document
        http://homepage.tudelft.nl/e3q6n/publications/1995/SP95TYLV/SP95TYLV.pdf
        */
        
        while (hasInflectionPoints
            && (resultingSigma < SIGMA.getValue(sigmaEnd))) {

            ScaleSpaceCurve curve;
            
//log.info("trimmedXOffset=" + trimmedXOffset + " trimmedYOffset=" + trimmedYOffset);

            if (lastCurve == null) {
                curve = scaleSpaceHelper.computeCurvature(edge, sigma, 
                    resultingSigma);
            } else {
                curve = scaleSpaceHelper.computeCurvature(
                    lastCurve.getXYCurve(), sigma, resultingSigma);
            }

            scaleSpaceMap.put(sigma, curve);

            hasInflectionPoints = (curve.getKIsZeroIdxSize() > 0);
            
            //log.info("sigma=" + sigma + " nZeros=" + curve.getKIsZeroIdxSize());

            if (hasInflectionPoints) {
                
                sigma = resultingSigma;

                resultingSigma *= sigmaPowerFactor;
            }

            lastCurve = curve;
        }
        
        return scaleSpaceMap;
    }
 
    /**
     * create a scale space curve image from the given sigma scale space inflection
     * points.
     * The points for which "t" and sigma are extracted are the inflection points,
     * that is, the zero-crossings in "t" versus curvature, that is, where 
     * the curvature is zero.
     * (The number of inflection points decreases as sigma increases.)
     * 
     * <pre>
     *        |    *
     * sigma  |   * *     **
     *        |   * *     **
     *        ----------------------
     *          scale free axis t
     * </pre>
     * @param scaleSpaceMap
     * @param edgeNumber a number placed in the ScaleSpaceCurveImage to help
     * identify the edge for debugging.
     * @param edgeLength the length of the edge.  this is the number of points
     * that went into the scale free length axis.
     * @return 
     */
    public ScaleSpaceCurveImage convertScaleSpaceMapToSparseImage(
        Map<Float, ScaleSpaceCurve> scaleSpaceMap, int edgeNumber,
        int edgeLength) {
        
        SortedMap<Float, ScaleSpaceCurve> sortedMap = 
            new TreeMap<Float, ScaleSpaceCurve>(
                new CurvatureScaleSpaceImageMaker.DescendingScaleSpaceComparator());
        
        sortedMap.putAll(scaleSpaceMap);
        
        Float maxSigma = sortedMap.isEmpty() ? null : sortedMap.firstKey();
        
        if ((maxSigma != null) && 
            sortedMap.get(maxSigma).getKIsZeroIdxSize() == 0) {
            
            sortedMap.remove(maxSigma);
            
            maxSigma = sortedMap.isEmpty() ? null : sortedMap.firstKey();
        }
        
        Iterator<Map.Entry<Float, ScaleSpaceCurve> > iter = 
            sortedMap.entrySet().iterator();
        
        int sImageRowIdx = 0;
        
        ScaleSpaceCurveImage spaceImage = new ScaleSpaceCurveImage(
            sortedMap.size());
        
        spaceImage.setEdgeNumber(edgeNumber);
        
        spaceImage.setEdgeSize(edgeLength);
        
        while (iter.hasNext()) {
            
            Map.Entry<Float, ScaleSpaceCurve> entry = iter.next();
            
            float sigma = entry.getKey().floatValue();
            
            ScaleSpaceCurve scaleSpaceCurve = entry.getValue();
            
            int nPoints = scaleSpaceCurve.getSize();
                            
            int nz = scaleSpaceCurve.getKIsZeroIdxSize();
                
            float[] row = new float[nz];
                        
            for (int i = 0; i < nz; i++) {

                int idx = scaleSpaceCurve.getKIsZeroIdx()[i];

                float t = (float)idx/(float)nPoints;
                
                row[i] = t;                
            }
            
            spaceImage.setRow(sImageRowIdx, row);
            
            spaceImage.setSigma(sImageRowIdx, sigma);
            
            spaceImage.setXYCoords(sImageRowIdx, scaleSpaceCurve.getKIsZeroX(),
                scaleSpaceCurve.getKIsZeroY());
            
            sImageRowIdx++;
        }
        
        return spaceImage;
    }
    
}
