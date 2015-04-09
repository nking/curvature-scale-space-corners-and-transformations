package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 given an image, creates scale space maps to find inflection points
    and creates scale space contours for the closed curves from those
    points.
 * 
 * Based upon the algorithm contained in
 * <pre>
 * IEEE 'TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE, VOL. PAMI-8, 
 * NO. 1. JANUARY 1986.  "Scale-Based Description and Recognition of Planar 
 * Curves and Two-Dimensional Shapes" by FARZIN MOKHTARIAN AND ALAN MACKWORTH
 * </pre>
 * 
 * @author nichole
 */
public final class CurvatureScaleSpaceImageMaker extends 
    AbstractCurvatureScaleSpaceMapper {
    
    /**
     * array of closed curve edges which are a subset of the list, 'edges'.
     * NOTE, the closedCurves coordinates here are referenced by the return 
     * value of createInflectionContours();
     */
    protected List<PairIntArray> closedCurves = 
        new ArrayList<PairIntArray>();
    
    public CurvatureScaleSpaceImageMaker(ImageExt input) {
        
        super(input);        
    }
    
    public CurvatureScaleSpaceImageMaker(ImageExt input, 
        List<PairIntArray> theEdges) {
        
        super(input, theEdges);
    }
    
    @Override
    protected void initialize() {
        
        super.initialize();
        
        createListOfClosedCurves();
    } 
    
    @Override
    protected void reinitializeSpecialization() {
    }

    /**
     * Create for each edge in the instance variable edges, 
     * X(t,sigma), Y(t,sigma), k(t, sigma) for sigma=1 until sigma=the level at 
     * which there are no more values of k that equal 0, that is until no 
     * more inflection points are found on the extremely smoothed curve.  
     * 
     
     */
    protected Map<PairIntArray, Map<Float, ScaleSpaceCurve> > 
        createScaleSpaceMetricsForInflectionPoints() {

        /*
        SIGMA=0:
           X(t,sigma), Y(t,sigma), k(t, sigma) and t 
              where t is the indexes normalized to the range 0 to 1.
        */

        Map<PairIntArray, Map<Float, ScaleSpaceCurve> > scaleSpaceMaps
            = new HashMap<PairIntArray, Map<Float, ScaleSpaceCurve> >();
                
        for (int i = 0; i < edges.size(); i++) {
            
            PairIntArray edge = edges.get(i);
            
            Map<Float, ScaleSpaceCurve> scaleSpaceMap = 
                createScaleSpaceMetricsForEdge2(edge);
    
           scaleSpaceMaps.put(edge, scaleSpaceMap);
                
        }
        
        return scaleSpaceMaps;
    }
        
    /**
     * NOTE: prefer use of createScaleSpaceMetricsForForEdge2() for now.
     * See notes below.
     * 
     * Create for an edge, 
     * X(t,sigma), Y(t,sigma), k(t, sigma) for sigma=1 until sigma=the level at 
     * which there are no more values of k that equal 0, that is until no 
     * more inflection points are found on the extremely smoothed curve.  
     * 
     * Note that the method re-uses iterative convolution, so each interval
     * of sigma is sqrt(2) times the previous.
     * Scale space images created to look for zero crossings in the curvature
     * are: sigma=1, sigma=sqrt(2), sigma=2, sigma=2*sqrt(2), sigma=4,
     * sigma=4*sqrt(2), sigma=8, sigma=8*sqrt(2), sigma=16, sigma=16*sqrt(2),
     * sigma=32, sigma=32*sqrt(2), sigma=64, sigma=64*sqrt(2), sigma=128.
     * Sometimes, the peaks for a contour do not close for the last non-zero
     * crossings convolution or for other peaks underneath.  The intervals
     * where the last peak occurs are skipped over. 
     * TODO: An efficient means of backtracking once a contour disappears could
     * be added to a method like this with an active contour finder in
     * order to minimize the number of convolutions (and the size of the kernels
     * of convolution) used.
     * For now, one should prefer the method createScaleSpaceMetricsForForEdge2
     * which uses a smaller interval for the iterative kernel sigma factor
     * so the peaks as a single point are present in these data.
     * (It uses sigma=1.189207115002721 instead of sqrt(2), so spans
     * sigma=1 to 4 in 7 steps rather than the 3 here.)
     * 
     * NOTE: one could imagine using these results with the mapper tailored
     * to tolerate an error of up to sqrt(2) in the peak (scale transformations)
     * with the understanding that lower peaks in the transformation 
     * solution should help reduce the error.  The cost function might need
     * to be altered for something like that.
     * 
     * @param edge
     * @return 
     */
    protected Map<SIGMA, ScaleSpaceCurve> createScaleSpaceMetricsForEdge(
    PairIntArray edge) {
        
        /*
        SIGMA=0:
           X(t,sigma), Y(t,sigma), k(t, sigma) and t 
              where t is the indexes normalized to the range 0 to 1.
        */
        
        ScaleSpaceCurvature scaleSpaceHelper = new ScaleSpaceCurvature();
            
        Map<SIGMA, ScaleSpaceCurve> scaleSpaceMap = new HashMap<SIGMA,
            ScaleSpaceCurve>();
        
        SIGMA sigma = SIGMA.ONE;

        // this increases by a factor of sqrt(2)
        float resultingSigma = SIGMA.getValue(sigma);

        boolean hasInflectionPoints = true;

        ScaleSpaceCurve lastCurve = null;

        while (hasInflectionPoints && (sigma != null) 
            && (resultingSigma < SIGMA.getValue(SIGMA.TWOHUNDREDANDFIFTYSIX))) {

            ScaleSpaceCurve curve;

            if (lastCurve == null) {
                curve = scaleSpaceHelper.computeCurvature(edge, sigma, 
                    resultingSigma);
            } else {
                curve = scaleSpaceHelper.computeCurvature(
                    lastCurve.getXYCurve(), sigma, resultingSigma);
            }

            scaleSpaceMap.put(sigma, curve);

            hasInflectionPoints = (curve.getKIsZeroIdxSize() > 0);
            
            log.fine("sigma=" + sigma + " nZeros=" + curve.getKIsZeroIdxSize());

            if (hasInflectionPoints) {
                //sigma = SIGMA.increaseToFactorBy2(resultingSigma);
                //resultingSigma *= 2;
                sigma = SIGMA.increaseToFactorBySQRT2(resultingSigma);

                resultingSigma *= Math.sqrt(2);
                
            }

            lastCurve = curve;
        }
        
        return scaleSpaceMap;
    }
    
    public ScaleSpaceCurveImage convertScaleSpaceMapToSparseImage(
        Map<Float, ScaleSpaceCurve> scaleSpaceMap, int edgeNumber,
        int edgeLength) {
        
        initialize();
        
        /*       |    *
          sigma  |   * *     **
                 |   * *     **
                 ----------------------
                   scale free axis t
        */
    
        SortedMap<Float, ScaleSpaceCurve> sortedMap = 
            new TreeMap<Float, ScaleSpaceCurve>(
                new DescendingScaleSpaceComparator());
        
        sortedMap.putAll(scaleSpaceMap);
        
        Float maxSigma = sortedMap.isEmpty() ? null : sortedMap.firstKey();
        
        if ((maxSigma != null) && 
            sortedMap.get(maxSigma).getKIsZeroIdxSize() == 0) {
            
            sortedMap.remove(maxSigma);
            
            maxSigma = sortedMap.isEmpty() ? null : sortedMap.firstKey();
        }
        
        Iterator<Entry<Float, ScaleSpaceCurve> > iter = 
            sortedMap.entrySet().iterator();
        
        int rowIdx = 0;
        
        ScaleSpaceCurveImage spaceImage = new ScaleSpaceCurveImage(
            sortedMap.size());
        
        spaceImage.setEdgeNumber(edgeNumber);
        
        spaceImage.setEdgeSize(edgeLength);
        
        while (iter.hasNext()) {
            
            Entry<Float, ScaleSpaceCurve> entry = iter.next();
            
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
            
            spaceImage.setRow(rowIdx, row);
            
            spaceImage.setSigma(rowIdx, sigma);
            
            spaceImage.setXYCoords(rowIdx, scaleSpaceCurve.getKIsZeroX(),
                scaleSpaceCurve.getKIsZeroY());
            
            rowIdx++;
        }
        
        return spaceImage;
    }
    
    /**
     * Create for an edge, 
     * X(t,sigma), Y(t,sigma), k(t, sigma) for sigma=1 until sigma=the level at 
     * which there are no more values of k that equal 0, that is until no 
     * more inflection points are found on the extremely smoothed curve.  
     * 
     * Note that the method re-uses iterative convolution, so each interval
     * of sigma is 2^(1/8) times the previous convolution result using a kernel
     * of sigma = 2^(1/8) each time.  The error in determining the peak of
     * a contour in the resulting scale space curve should be < 10%.
     * TODO: calculate error in peak determination....
     * 
     * @param edge
     * @return 
     */
    protected Map<Float, ScaleSpaceCurve> createScaleSpaceMetricsForEdge2(
    PairIntArray edge) {
        
        // if use 2^(1/8) as a sigma factor should result in an error less than 10%
        // in determing the peak of a contour.  smaller factors have smaller
        // errors than that.
        float factor = (float)Math.pow(2, 1./32.);
        
        /*
        SIGMA=0:
           X(t,sigma), Y(t,sigma), k(t, sigma) and t 
              where t is the indexes normalized to the range 0 to 1.
        */
        
        ScaleSpaceCurvature scaleSpaceHelper = new ScaleSpaceCurvature();
            
        Map<Float, ScaleSpaceCurve> scaleSpaceMap = new HashMap<Float,
            ScaleSpaceCurve>();
        
        float sigma = SIGMA.getValue(SIGMA.ONE);

        float resultingSigma = sigma;

        boolean hasInflectionPoints = true;

        ScaleSpaceCurve lastCurve = null;

        while (hasInflectionPoints
            && (resultingSigma < SIGMA.getValue(SIGMA.TWOHUNDREDANDFIFTYSIX))) {

            ScaleSpaceCurve curve;

            if (lastCurve == null) {
                curve = scaleSpaceHelper.computeCurvature(edge, sigma, 
                    resultingSigma);
            } else {
                curve = scaleSpaceHelper.computeCurvature(
                    lastCurve.getXYCurve(), sigma, resultingSigma);
            }

            scaleSpaceMap.put(sigma, curve);

            hasInflectionPoints = (curve.getKIsZeroIdxSize() > 0);
            
            log.fine("sigma=" + sigma + " nZeros=" + curve.getKIsZeroIdxSize());

            if (hasInflectionPoints) {
                
                sigma = resultingSigma;

                resultingSigma *= factor;
            }

            lastCurve = curve;
        }
        
        return scaleSpaceMap;
    }

    protected void createListOfClosedCurves() {
        
        // parse edges to find only the closed curves 
        
        if (state.ordinal() < CurvatureScaleSpaceMapperState.INITIALIZED.ordinal()) {
            initialize();
        }
        
        closedCurves.clear();
        
        for (PairIntArray edge : edges) {
            if (edge instanceof PairIntArrayWithColor) {
                if (((PairIntArrayWithColor)edge).getColor() == 1) {
                    closedCurves.add(edge);
                }
            }
        }
    }
    
    public List<PairIntArray> getClosedCurves() {
        return closedCurves;
    }
    
    private static class DescendingScaleSpaceComparator implements 
        Comparator<Float> {
        
        @Override
        public int compare(Float o1, Float o2) {
            return o2.compareTo(o1);
        }
    }
}
