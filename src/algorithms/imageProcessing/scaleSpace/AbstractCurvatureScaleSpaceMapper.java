package algorithms.imageProcessing.scaleSpace;

import algorithms.imageProcessing.CannyEdgeFilter;
import algorithms.imageProcessing.CannyEdgeFilterSettings;
import algorithms.imageProcessing.ClosedCurveAndJunctionFinder;
import algorithms.imageProcessing.EdgeExtractorWithJunctions;
import algorithms.imageProcessing.EdgeFilterProducts;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.IEdgeExtractor;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.SkylineExtractor;
import algorithms.misc.HistogramHolder;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayComparator;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public abstract class AbstractCurvatureScaleSpaceMapper {

    protected CurvatureScaleSpaceMapperState state = 
        CurvatureScaleSpaceMapperState.UNINITIALIZED;
    
    protected GreyscaleImage img;
    
    protected final ImageExt originalImg;
    
    /**
     * edges extracted from image.  if an instance of PairIntArrayWithColor
     * is present, that holds a color field in which a value of '1' means
     * the curve is closed.
     */
    protected List<PairIntArray> edges = new ArrayList<PairIntArray>();
    
    /**
     * if extractSkline is true, this is populated with the sky line
     * edge(s).
     */
    protected final List<PairIntArray> skylineEdges = new ArrayList<PairIntArray>();
        
    protected boolean doNotNormalizeByHistogram = false;
    
    protected boolean doNotFindJunctions = false;
    
    protected boolean useLineDrawingMode = false;
        
    protected boolean useLowestHighIntensityCutoff = false;
    
    protected boolean useLowHighIntensityCutoff = false;
    
    protected boolean extractSkyline = false;
    
    protected final int trimmedXOffset;
    
    protected final int trimmedYOffset;
    
    protected boolean useOutdoorMode = false;
    
    protected HistogramHolder imgHistogram = null;
    
    protected EdgeFilterProducts filterProducts = null;
    
    /**
     * map with key = center of junction pixel coordinates; 
     * value = set of adjacent pixels when there are more than the preceding 
     * and next.
     */
    protected Map<Integer, Set<Integer>> junctionMap = new HashMap<Integer, Set<Integer>>();

    /**
     * map with key = pixel coordinates of all pixels involved in junctions;
     * value = PairInt holding index of edge that pixel is located in and
     * holding the index within that edge of the pixel.
     * for example, a pixel located in edges(0) at offset=100
     * would have PairInt(0, 100) as a value.
     */
    protected Map<Integer, PairInt> junctionLocationMap = new HashMap<Integer, PairInt>();
     
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * constructor w/ input image which is operated on.  the same instance
     * input is modified by this class.
     * 
     * @param input 
     */
    public AbstractCurvatureScaleSpaceMapper(ImageExt input) {
        
        img = input.copyToGreyscale();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        originalImg = (ImageExt)input.copyImage();
            
        int[] offsetXY = imageProcessor.shrinkImageToFirstNonZeros(img);
        
        trimmedXOffset = offsetXY[0];
        
        trimmedYOffset = offsetXY[1];
    }
    
    public AbstractCurvatureScaleSpaceMapper(GreyscaleImage input) {
        
        img = input.copyImage();
        
        ImageProcessor imageProcessor = new ImageProcessor();

        originalImg = ImageIOHelper.convertImage(input);
            
        int[] offsetXY = imageProcessor.shrinkImageToFirstNonZeros(img);
        
        trimmedXOffset = offsetXY[0];
        
        trimmedYOffset = offsetXY[1];
    }
    
    /**
     * constructor with input image and the already extracted edges.
     * The input image is needed only for debugging purposes and 
     * may be removed as an argument after testing is complete.
     * @param input
     * @param theEdges 
     */
    public AbstractCurvatureScaleSpaceMapper(ImageExt input, 
        List<PairIntArray> theEdges) {
        
        img = input.copyToGreyscale();
        
        ImageProcessor ImageProcessor = new ImageProcessor();
        
        originalImg = (ImageExt)input.copyImage();
        
        int[] offsetXY = ImageProcessor.shrinkImageToFirstNonZeros(img);
        
        trimmedXOffset = offsetXY[0];
        
        trimmedYOffset = offsetXY[1];
        
        this.edges = new ArrayList<PairIntArray>(theEdges);
        
        state = CurvatureScaleSpaceMapperState.INITIALIZED;
    }

    /**
     * apply histogram normalization before processing.  For some images, this
     * will increase the contrast of fainter features.
     */
    public void doNotPerformHistogramEqualization() {
        this.doNotNormalizeByHistogram = true;
    }
    
    public void doNotFindJunctions() {
        doNotFindJunctions = true;
    }
    
    public void useLowestHighIntensityCutoff() {
        useLowestHighIntensityCutoff = true;
    }
    
    public void useLowHighIntensityCutoff() {
        useLowHighIntensityCutoff = true;
    }
    
    public void useLineDrawingMode() {
        useLineDrawingMode = true;
    }
    
    public void useOutdoorMode() {
        
        useLowestHighIntensityCutoff();
        
        useOutdoorMode = true;
    }
    
    /**
     * set the edge detector to create edges that are better for outdoor
     * conditions and also extract the skyline from the intermediate
     * image products.  Note that the skyline extraction is currently
     * a long running process.
     */
    public void useOutdoorModeAndExtractSkyline() {
        
        useLowestHighIntensityCutoff();
        
        useOutdoorMode = true;
        
        extractSkyline = true;
    }
    
    protected void initialize() {
        
        if (state.ordinal() < 
            CurvatureScaleSpaceMapperState.INITIALIZED.ordinal()) {
            
            // (1) apply an edge filter
            applyEdgeFilter();
            
            if (extractSkyline && skylineEdges.isEmpty()) {
                
                List<PairIntArray> skyEdges = extractSkyline();
                                    
                skylineEdges.addAll(skyEdges);
            }
            
            // (2) extract edges and create junction maps
            extractEdges();
            
            //TODO: note that there may be a need to search for closed
            //      curves in the EdgeContourExtractor instead of here
            //      in order to create shapes instead of creating
            //      lines preferentially.
            markTheClosedCurves();
            
            state = CurvatureScaleSpaceMapperState.INITIALIZED;
        }
    }
    
    protected abstract void reinitializeSpecialization();
    
    protected void reinitialize(float cannyLowThreshold, float additionalBlurSigma) {
        
        if (state.ordinal() < 
            CurvatureScaleSpaceMapperState.INITIALIZED.ordinal()) {
            
            throw new IllegalStateException(
                "initialize() should have been invoked before reinitialize(...)");
            
        } else {
                       
            GreyscaleImage input = filterProducts.getGradientXY().copyImage();
            GreyscaleImage gTheta = filterProducts.getTheta();
            HistogramHolder hist = imgHistogram;
        
            reinitializeSpecialization();
            
            edges.clear();
            
            if (additionalBlurSigma > 0) {
                
                // (1) apply an edge filter
                
                CannyEdgeFilter filter = new CannyEdgeFilter();

                CannyEdgeFilterSettings settings = getCannyEdgeFilterSettings();
        
                filter.setSetters(settings);
                
                filter.setAdditionalImageBlur(additionalBlurSigma);
                
                filter.overrideLowThreshold(cannyLowThreshold);
                
                img = originalImg.copyToGreyscale();
                
                filter.applyFilter(img);
                
                filterProducts = filter.getEdgeFilterProducts();
        
                imgHistogram = filter.getImgHistogram();
                
                img = filterProducts.getGradientXY();
                
            } else {
            
                // (1) re-apply an edge filter

                CannyEdgeFilter filter = new CannyEdgeFilter();
                CannyEdgeFilterSettings settings = getCannyEdgeFilterSettings();

                filter.setSetters(settings);

                filter.overrideLowThreshold(cannyLowThreshold);

                filter.reApply2LayerFilter(input, gTheta, hist);

                img = input;
                
            }
            
            // (2) extract edges
            extractEdges();
            
            //TODO: note that there may be a need to search for closed
            //      curves in the EdgeContourExtractor instead of here
            //      in order to create shapes instead of creating
            //      lines preferentially.
            // (3) look for t-junctions and closed curves
            markTheClosedCurves();
            
            state = CurvatureScaleSpaceMapperState.INITIALIZED;
        }
    }

    protected void applyEdgeFilter() {
        
        CannyEdgeFilter filter = new CannyEdgeFilter();

        CannyEdgeFilterSettings settings = getCannyEdgeFilterSettings();
        
        filter.setSetters(settings);
                
        filter.applyFilter(img);
        
        filterProducts = filter.getEdgeFilterProducts();
        
        imgHistogram = filter.getImgHistogram();
                                
        state = CurvatureScaleSpaceMapperState.EDGE_FILTERED;
    }
    
    /**
     * use the output from the canny edge filter and the theta image
     * it produced to find the skyline in the image, extract it as
     * an edge(s) and make several pixels adjustment to align with
     * the output of the canny edge filter (which should be a binary
     * image with edges as 0's).   NOTE that this is currently a long
     * running process.
     * 
     */
    protected List<PairIntArray> extractSkyline() {
                
        try {            
            
            CannyEdgeFilterSettings settings = getCannyEdgeFilterSettings();
            
            SkylineExtractor skylineExtractor = new SkylineExtractor();
            
            PairIntArray outputSkyCentroid = new PairIntArray();
            GreyscaleImage out = skylineExtractor.createSkyline(
                filterProducts.getTheta(), filterProducts.getGradientXY(),
                this.originalImg, settings, outputSkyCentroid);
             
            List<PairIntArray> skyEdges = skylineExtractor.getSkylineEdges();
            
            if (skyEdges == null) {
                return new ArrayList<PairIntArray>();
            }
            
            Collections.sort(skyEdges, new PairIntArrayComparator());

            //reverse the list so the edges with largest numbers of points are
            // at smaller indexes
            int n = skyEdges.size();
            if (n > 1) {
                int end = n >> 1;
                for (int i = 0; i < end; i++) {
                    int idx2 = n - i - 1;                
                    PairIntArray swap = skyEdges.get(i);
                    skyEdges.set(i, skyEdges.get(idx2));
                    skyEdges.set(idx2, swap);
                }
            }
            
            return skyEdges;
            
        } catch(IOException e) {
            log.severe(e.getMessage());
            
        } catch(NoSuchAlgorithmException e) {
            log.severe(e.getMessage());
        }
        
        return new ArrayList<PairIntArray>();
    }
    
    public CannyEdgeFilterSettings getCannyEdgeFilterSettings() {
    
        CannyEdgeFilterSettings settings = new CannyEdgeFilterSettings();
        
        if (useOutdoorMode) {
            settings.setUseOutdoorMode();
        }
        
        if (doNotNormalizeByHistogram) {
            settings.setDoNotNormalizeByHistogram();
        }
        if (useLineDrawingMode) {
            settings.setUseLineDrawingMode();
        }
        
        if (useLowestHighIntensityCutoff) {
            settings.setOverrideHighThreshold(1.0f);
        } else if (useLowHighIntensityCutoff) {
            settings.setOverrideHighThreshold(2.0f);
        }
        
        return settings;
    }
    
    protected void extractEdges() {
        
        IEdgeExtractor contourExtractor;
        
        contourExtractor = new EdgeExtractorWithJunctions(img);
        
        //TODO: for closed curves, might consider always ordering them the
        // same here (counter clock wise)
        
        List<PairIntArray> tmpEdges = contourExtractor.findEdges();
       
        edges.clear();
        edges.addAll(tmpEdges);
        
        if (!doNotFindJunctions && (contourExtractor instanceof 
            EdgeExtractorWithJunctions)) {

            junctionMap.clear();
            junctionLocationMap.clear();
            
            Map<Integer, Set<Integer>> jm
                = ((EdgeExtractorWithJunctions) contourExtractor)
                .getJunctionMap();

            if (!jm.isEmpty()) {

                junctionMap.putAll(jm);

                junctionLocationMap.putAll(
                    ((EdgeExtractorWithJunctions) contourExtractor)
                    .getLocatorForJunctionAssociatedPoints()
                );
            }
        }
        
        state = CurvatureScaleSpaceMapperState.EDGES_EXTRACTED;
        
        log.fine("edges extracted");
    }

    protected void markTheClosedCurves() {
        
        ClosedCurveAndJunctionFinder ccjFinder = 
            new ClosedCurveAndJunctionFinder();
        
        ccjFinder.findClosedCurves(edges);
       
    }
  
    public List<PairIntArray> getEdges() {
        return edges;
    }
    
    public List<PairIntArray> getSkylineEdges() {
        return skylineEdges;
    }

    protected void addCurveToImage(PairIntArray edge, GreyscaleImage input, 
        int nExtraForDot, int value) {
        
        for (int i = 0; i < edge.getN(); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            for (int dx = -1 * nExtraForDot; dx < (nExtraForDot + 1); dx++) {
                float xx = x + dx;
                if ((xx > -1) && (xx < (input.getWidth() - 1))) {
                    for (int dy = -1 * nExtraForDot; dy < (nExtraForDot + 1); 
                        dy++) {
                        
                        float yy = y + dy;
                        if ((yy > -1) && (yy < (input.getHeight() - 1))) {
                            input.setValue((int) xx, (int) yy, value);
                        }
                    }
                }
            }
        }
    }
    
    public boolean getInitialized() {
        
        return (state.ordinal() >= 
            CurvatureScaleSpaceMapperState.INITIALIZED.ordinal());
    }
    
    public GreyscaleImage getImage() {
        return img;
    }

    public Image getOriginalImage() {
        return originalImg;
    }
    
    public int getTrimmedXOffset() {
        return trimmedXOffset;
    }
    
    public int getTrimmedYOffset() {
        return trimmedYOffset;
    }
   
    public EdgeFilterProducts getEdgeFilterProducts() {
        return filterProducts;
    }
    
    public Map<Integer, Set<Integer>> getJunctionMap() {        
        return junctionMap;
    }
    
    /*
    The making of a curvature scale space image is in
    "Scale-Based Description and Recognition of Planar Curves and Two-Dimensional
    Shapes" by FARZIN MOKHTARIAN AND ALAN MACKWORTH
    IEEE 'TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE,
    VOL. PAMI-8, NO. 1. JANUARY 1986
    https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0CCIQFjAA&url=https%3A%2F%2Fwww.cs.ubc.ca%2F~mack%2FPublications%2FIEEE-PAMI86.pdf&ei=jiIFVJGNLIa0igLv74DgDw&usg=AFQjCNHj7v2JaUDqSFkQZSNOSpqBbfbOAQ&sig2=L08nOsKD1Mw_XJX-EPmY-w&bvm=bv.74115972,d.cGE
    planar curve:
    f_curve = {x(t), y(t)}
    t = linear function of the path length bounded by values [0, 1], that is,
    one can make this by scaling the range os indexes for x and y
    for a curve to values between 0 and 1.
    If f_curve is closed, x(t) and y(t) are periodic functions.
    The curvature, k, is the the change of the angle of the tangent line at
    point P on arc s with respect to the arc length s.
    #  /
    # /
    #/      /|
    P   ds / | dy
    #/      /__|
    # /        dx
    #  /
    /
    / theta
    ._________
    ds^2 = dx^2 + dy^2
    ds = sqrt(dx^2 + dy^2) = sqrt(1 + (dy/dx)^2)*dx = sqrt((dx/dy)^2 + 1)*dy
    k = dTheta/ds = 1/rho
    where rho is the radius of the circle of curvature at point P
    dTheta   dTheta   dx
    ------ = ------ * --
    ds       dx      ds
    theta = tan^-1 (dy / dx)
    d                 d/dx (dy/dx)      d^2y/dx^2
    dTheta/dx = -- arctan(dy/dx) = ------------- = -------------
    dx                 1 + (dy/dx)^2   1 + (dy/dx)^2
    dx    1             1
    -- = ------ = -------------------
    ds   ds/dx    sqrt(1 + (dy/dx)^2)
    and use y' = (dy/dx)
    and use y" = (d^2y/dx^2)
    dTheta          y"                  1
    k =  ------ = --------------- * -------------------
    ds     (1 + (dy/dx)^2)   sqrt(1 + (dy/dx)^2)
    d^2y/dx^2
    = ---------------------  for planar curves
    (1 + (dy/dx)^2)^(1.5)
     * the sign of k is + if y" is + and is - if y" is -. the absolute value
    might be used instead though.
    NOTE that if dy/dx doesn’t exist at a point, such as where
    the tangent line is parallel to the y-axis,
    one can invert the y/x relationships in k to x/y
    (d^2x/dy^2)
    k = ---------------------
    (1 + (dx/dy)^2)^(1.5)
     * Need to express k in terms of a function of t, the parameteric form of k
    dTheta   dTheta   dt     1     dTheta
    k = ------ = ------ * -- = ----- * ------
    ds       dt      ds   ds/dt     dt
    where (ds/dt)^2 = (dx/dt)^2 + (dy/dt)^2
    dy   dy/dt
    tan(theta) = -- = -----
    dx   dx/dt
    d
    --(tan(theta)) = sec^2(theta) * (dTheta/dt)
    dt
    (d^2y/dt^2)   (dy/dt)*(d^2x/dt^2)
    = ----------- - -------------------
    (dx/dt)           (dx/dt)^2
    (d^2y/dt^2)*(dx/dt) - (dy/dt)*(d^2x/dt^2)
    = -----------------------------------------
    (dx/dt)^2
    dTheta        1         (d^2y/dt^2)*(dx/dt) - (dy/dt)*(d^2x/dt^2)
    so ------ = ------------ * -----------------------------------------
    dt      sec^2(theta)                   (dx/dt)^2
    1           (d^2y/dt^2)*(dx/dt) - (dy/dt)*(d^2x/dt^2)
    = ---------------- * -----------------------------------------
    1 + tan^2(theta)                 (dx/dt)^2
    1           (d^2y/dt^2)*(dx/dt) - (dy/dt)*(d^2x/dt^2)
    = ---------------- * -----------------------------------------
    1 + (dy/dt)^2                  (dx/dt)^2
    ---------
    (dx/dt)^2
    (d^2y/dt^2)*(dx/dt) - (dy/dt)*(d^2x/dt^2)
    = ------------------------------------------
    (dx/dt)^2 + (dy/dt)^2
     * now can return to
    1      dTheta
    k_geodesic = ----- *  ------
    ds/dt     dt
    (d^2y/dt^2)*(dx/dt) - (dy/dt)*(d^2x/dt^2)
    = ---------------------------------------------------------
    (((dx/dt)^2 + (dy/dt)^2)^(0.5)) * ((dx/dt)^2 + (dy/dt)^2)
    (d^2y/dt^2)*(dx/dt) - (dy/dt)*(d^2x/dt^2)
    = -----------------------------------------
    ((dx/dt)^2 + (dy/dt)^2)^(1.5)
    REWRITE in terms of code:
    X_dot(t,o~) * Y_dot_dot(t,o~) - Y_dot(t,o~) * X_dot_dot(t,o~)
    k(t,o~) = ----------------------------------------------------------------
    (X^2(t,o~) + Y^2(t,o~))^1.5
    where o~ denotes the width of the Gaussian
    convolve X and Y w/ one dimensional gaussian kernel each:
    X(t, o~) = Integ(x(v) * exp(-(v)^2/2o~^2) * dv)
    Y(t, o~) = Integ(y(v) * exp(-(v)^2/2o~^2) * dv)
    Integ denotes the integral evaluated from -infinity to +infinity.
    First Deriv:
    X_dot(t,o~) = Integ(x(v) * (-2*(v)) * exp(-(v)^2/2o~^2) * dv)
    Y_dot(t,o~) = Integ(y(v) * (-2*(v)) * exp(-(v)^2/2o~^2) * dv)
    Second Deriv:
    X_dot_dot(t,o~) = Integ(x(v) * (-2 + 4 * (v)^2)) * exp(-(v)^2/2o~^2) * dv)
    Y_dot_dot(t,o~) = Integ(y(v) * (-2 + 4 * (v)^2)) * exp(-(v)^2/2o~^2) * dv)
    The curvture of a straight line is zero.
    Points where k = 0 are called the points of inflection.
     */

}
