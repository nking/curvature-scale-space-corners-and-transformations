package algorithms.imageProcessing.scaleSpace;

import algorithms.compGeometry.PerimeterFinder2;
import algorithms.connected.ConnectedPointsFinder;
import algorithms.imageProcessing.CannyEdgeFilterAdaptive;
import algorithms.imageProcessing.EdgeFilterProducts;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.PostLineThinnerCorrections;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.features.mser.MSEREdges;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

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
public final class CurvatureScaleSpaceImageMaker {
    
    protected CurvatureScaleSpaceMapperState state = 
        CurvatureScaleSpaceMapperState.UNINITIALIZED;
    
    /**
     * array of closed curve clockwise ordered edges;
     */
    protected List<PairIntArray> closedCurves = 
        new ArrayList<PairIntArray>();
    
    protected boolean useLineDrawingMode = false;
    
    protected EdgeFilterProducts filterProducts = null;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public CurvatureScaleSpaceImageMaker(ImageExt input) {
        
        initialize(input);
    }
    
    public CurvatureScaleSpaceImageMaker(ImageExt input, boolean useLineDrawingMode) {
        
        this.useLineDrawingMode = useLineDrawingMode;
        
        initialize(input);
    }
    
    protected void initialize(ImageExt img) {
        
        if (state.ordinal() < 
            CurvatureScaleSpaceMapperState.INITIALIZED.ordinal()) {
            
            MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
            
            int w = img.getWidth();
            int h = img.getHeight();
            
            // extract closed curve edges
            
            if (useLineDrawingMode) {
                
                GreyscaleImage img2 = img.copyToGreyscale2();
                
                CannyEdgeFilterAdaptive filter = new CannyEdgeFilterAdaptive();
                filter.setToUseLineDrawingMode();
            filter.setToDebug();
                filter.applyFilter(img2);
        
                filterProducts = filter.getFilterProducts();
        
                GreyscaleImage gXY = filterProducts.getGradientXY();
                //PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
                //pltc.extremeCornerRemover(gXY);
                
                TIntSet nzs = new TIntHashSet();
                for (int i = 0; i < gXY.getNPixels(); ++i) {
                    if (gXY.getValue(i) > 0) {
                        nzs.add(i);
                    }
                }

                ConnectedPointsFinder cFinder = new ConnectedPointsFinder(gXY.getWidth(), 
                    gXY.getHeight());
                cFinder.setToUse8Neighbors();
                cFinder.findConnectedPointGroups(nzs);

                PerimeterFinder2 finder2 = new PerimeterFinder2();

                for (int i = 0; i < cFinder.getNumberOfGroups(); ++i) {
                    TIntSet set = cFinder.getXY(i);
                    if (isOnImageBounds(set, w, h)) {
                        continue;
                    }
                    PairIntArray ordered = null;
                    try {
                        ordered = finder2.orderTheBoundary(set, gXY.getWidth(), 
                            gXY.getHeight());
                    } catch (Exception e) {
                        continue;
                    }
                    PairIntArrayWithColor closedCurve 
                        = new PairIntArrayWithColor(ordered);
                    closedCurve.setAsClosedCurve();

                    closedCurves.add(closedCurve);
                }
                
            } else {
                MSEREdges mserEdges = new MSEREdges(img);
                mserEdges.setToLowerContrast();
                mserEdges.mergeAndExtractEdges();
                List<TIntSet> edgeSets = mserEdges.getEdges();
            
                PerimeterFinder2 finder2 = new PerimeterFinder2();

                for (int i = 0; i < edgeSets.size(); ++i) {

                    TIntSet set = edgeSets.get(i);

                    if (isOnImageBounds(set, w, h)) {
                        continue;
                    }
                    
                    PairIntArray ordered = null;
                    try {
                        ordered = finder2.orderTheBoundary(
                            set, w, h);
                    } catch (Exception e) {
                        continue;
                    }
                    
                    // NOTE: ordered is clockwise, but other code expects 
                    //  CCW, so reversing for now
                    ordered.reverse();
                    
                    PairIntArrayWithColor closedCurve = new PairIntArrayWithColor(
                        ordered);
                    closedCurve.setAsClosedCurve();

                    closedCurves.add(closedCurve);
                }
            }
            
            state = CurvatureScaleSpaceMapperState.EDGES_EXTRACTED;

            state = CurvatureScaleSpaceMapperState.INITIALIZED;
        }
                
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
                
        for (int i = 0; i < closedCurves.size(); i++) {
            
            PairIntArray edge = closedCurves.get(i);
            
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
                
        /*       |    *
          sigma  |   * *     **
                 |   * *     **
                 ----------------------
                   scale free axis t
        */
    
        CurvatureScaleSpaceCurvesMaker csscMaker = new CurvatureScaleSpaceCurvesMaker();
        
        ScaleSpaceCurveImage spaceImage = 
            csscMaker.convertScaleSpaceMapToSparseImage(
            scaleSpaceMap, edgeNumber, edgeLength);
        
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
 
        CurvatureScaleSpaceCurvesMaker csscMaker = new CurvatureScaleSpaceCurvesMaker();
        
        // if use 2^(1/8) as a sigma factor should result in an error less than 10%
        // in determing the peak of a contour.  smaller factors have smaller
        // errors than that.
        float factor = (float)Math.pow(2, 1./32.);
        
        /*
        SIGMA=0:
           X(t,sigma), Y(t,sigma), k(t, sigma) and t 
              where t is the indexes normalized to the range 0 to 1.
        */
        
        Map<Float, ScaleSpaceCurve> scaleSpaceMap = 
            csscMaker.createScaleSpaceMetricsForEdge(edge, factor,
                SIGMA.ONE, SIGMA.TWOHUNDREDANDFIFTYSIX);
        
        return scaleSpaceMap;
    }
    
    public List<PairIntArray> getClosedCurves() {
        return closedCurves;
    }

    private boolean isOnImageBounds(TIntSet set, int w, int h) {

        TIntIterator iter = set.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int y = pixIdx/w;
            int x = pixIdx - (y * w);
            if (x == 0 || y == 0 || (x == (w - 1)) || (y == (h - 1))) {
                return true;
            }
        }
        
        return false;
    }
    
    public static class DescendingScaleSpaceComparator implements 
        Comparator<Float> {
        
        @Override
        public int compare(Float o1, Float o2) {
            return o2.compareTo(o1);
        }
    }
}
