package algorithms.imageProcessing;

import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayComparator;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
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
    protected List<PairIntArray> skylineEdges = new ArrayList<PairIntArray>();
    
    protected PairIntArray tCorners = new PairIntArray();
    
    protected boolean doNotNormalizeByHistogram = false;
    
    protected boolean useLineDrawingMode = false;
    
    protected boolean doNotStraightenEdges = false;
    
    protected boolean useLowestHighIntensityCutoff = false;
    
    protected boolean useLowHighIntensityCutoff = false;
    
    protected boolean extractSkyline = false;
    
    protected final int trimmedXOffset;
    
    protected final int trimmedYOffset;
    
    protected boolean useOutdoorMode = false;
    
    protected HistogramHolder imgHistogram = null;
    
    protected GreyscaleImage gradientXY = null;
    
    protected GreyscaleImage theta = null;
     
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * constructor w/ input image which is operated on.  the same instance
     * input is modified by this class.
     * 
     * @param input 
     */
    public AbstractCurvatureScaleSpaceMapper(ImageExt input) {
        
        img = input.copyToGreyscale();
        
        ImageProcessor ImageProcessor = new ImageProcessor();
        
        originalImg = (ImageExt)input.copyImage();
            
        int[] offsetXY = ImageProcessor.shrinkImageToFirstNonZeros(img);
        
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
            
            if (extractSkyline) {
                extractSkyline();
            }
            
            // (2) extract edges
            extractEdges();
            
            //TODO: note that there may be a need to search for closed
            //      curves in the EdgeContourExtractor instead of here
            //      in order to create shapes instead of creating
            //      lines preferentially.
            // (3) look for t-junctions and closed curves
            markTheClosedCurvesAndTCorners();
            
            state = CurvatureScaleSpaceMapperState.INITIALIZED;
        }
    }
    
    protected abstract void reinitializeSpecialization();
    
    protected void reinitialize(float cannyLowThreshold) {
        
        if (state.ordinal() < 
            CurvatureScaleSpaceMapperState.INITIALIZED.ordinal()) {
            
            throw new IllegalStateException(
                "initialize() should have been invoked before reinitialize(...)");
            
        } else {
            
            GreyscaleImage input = gradientXY.copyImage();
            GreyscaleImage gTheta = theta;
            HistogramHolder hist = imgHistogram;
        
            reinitializeSpecialization();
            
            edges.clear();
            
            // (1) apply an edge filter
            
            CannyEdgeFilter filter = new CannyEdgeFilter();
            CannyEdgeFilterSettings settings = getCannyEdgeFilterSettings();
        
            filter.setSetters(settings);
            
            filter.overrideLowThreshold(cannyLowThreshold);
                        
            filter.reApply2LayerFilter(input, gTheta, hist);
            
            img = input;
            
            // (2) extract edges
            extractEdges();
            
            //TODO: note that there may be a need to search for closed
            //      curves in the EdgeContourExtractor instead of here
            //      in order to create shapes instead of creating
            //      lines preferentially.
            // (3) look for t-junctions and closed curves
            markTheClosedCurvesAndTCorners();
            
            state = CurvatureScaleSpaceMapperState.INITIALIZED;
        }
    }

    protected void applyEdgeFilter() {
        
        CannyEdgeFilter filter = new CannyEdgeFilter();

        CannyEdgeFilterSettings settings = 
            getCannyEdgeFilterSettings();
        
        filter.setSetters(settings);
                
        filter.applyFilter(img);
        
        gradientXY = filter.getGradientXY();
        
        theta = filter.getTheta();
        
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
    protected void extractSkyline() {
        
        PairIntArray outputSkyCentroid = new PairIntArray();
        
        try {            
            
            CannyEdgeFilterSettings settings = getCannyEdgeFilterSettings();
            
            SkylineExtractor skylineExtractor = new SkylineExtractor();
            
            GreyscaleImage out = skylineExtractor.createSkyline(theta, 
                gradientXY, this.originalImg, settings, outputSkyCentroid);
            
            EdgeExtractor contourExtractor = new EdgeExtractor(out);
            
            if ((img.getWidth() > 300) && (img.getHeight() > 300)) {
                contourExtractor.overrideEdgeSizeLowerLimit(50);
            }
            
            contourExtractor.useRepeatConnectAndTrim();

            List<PairIntArray> skyEdges = contourExtractor.findEdges();
            
            if (skyEdges.isEmpty()) {
                return;
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

            skylineEdges.addAll(skyEdges);
            
        } catch(IOException e) {
            log.severe(e.getMessage());
            
        } catch(NoSuchAlgorithmException e) {
            log.severe(e.getMessage());
        }
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
        
        EdgeExtractor contourExtractor;
        
        if (doNotStraightenEdges) { 
            contourExtractor = new EdgeExtractor(img);
        } else {
            contourExtractor = new EdgeExtractor(img, gradientXY);
        }
        
        List<PairIntArray> tmpEdges = contourExtractor.findEdges();
       
        edges.clear();
        edges.addAll(tmpEdges);
        
        GreyscaleImage output = img.createWithDimensions();
        for (PairIntArray edge : edges) {
            addCurveToImage(edge, output, 0, 255);
        }
        
        img = output;
        
        state = CurvatureScaleSpaceMapperState.EDGES_EXTRACTED;
        
        log.fine("edges extracted");
    }

    protected void markTheClosedCurvesAndTCorners() {
        
        ClosedCurveAndJunctionFinder ccjFinder = 
            new ClosedCurveAndJunctionFinder();
        
        ccjFinder.findClosedCurves(edges);
        
        /*
              ---@---
                 |
                 |
        */
        // Fill small gaps in edge contours. When the gap forms a T-junction, 
        // mark it as a T-corner.
        // May have already filled these in in the contourextractor.
        // TODO: revisit this and consider tracking t-junctions in the
        // contourExtractor when filling a gap
        // see http://www.diva-portal.org/smash/get/diva2:457189/FULLTEXT01.pdf
        //  pg 26
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
    
    private void applyBestXShift(List<PairIntArray> skyEdges, GreyscaleImage img, 
        int deltaX, int absoluteMaxDelta) {
        
        applyBestXOrYShift(skyEdges, img, deltaX, absoluteMaxDelta, true);
    }
    private void applyBestYShift(List<PairIntArray> skyEdges, GreyscaleImage img, 
        int deltaY, int absoluteMaxDelta) {
        
        applyBestXOrYShift(skyEdges, img, deltaY, absoluteMaxDelta, false);
    }
    private void applyBestXOrYShift(List<PairIntArray> skyEdges, GreyscaleImage img, 
        int delta, int absoluteMaxDelta, boolean shiftX) {
        
        int bestShift = 0;
        int bestSum = Integer.MIN_VALUE;
        
        // 0 to 12, dy=+1
        // 0 to -12, dy=-1
        
        int start = 0;
        int stop = 0;
        if (delta > 0) {
            stop = absoluteMaxDelta;
        } else {
            //-12 to 0
            start = -1*absoluteMaxDelta;
            delta *= -1;
        }
        
        for (int shift = start; shift < stop; shift += delta) {
            
            int sum = 0;
            
            for (PairIntArray edge : skyEdges) {
                for (int i = 0; i < edge.getN(); i++) {
                    int x = edge.getX(i);
                    int y = edge.getY(i);
                    
                    int xs = x;
                    int ys = y;
                    
                    if (shiftX) {
                        xs += shift;
                    } else {
                        ys += shift;
                    }
                    
                    int v = img.getValue(xs, ys);
                    
                    if (v > 0) {
                        sum++;
                    }
                }
            }
            
            if (sum > bestSum) {
                bestSum = sum;
                bestShift = shift;
            }
        }
        
        if (bestShift != Integer.MIN_VALUE) {
            for (PairIntArray edge : skyEdges) {
                for (int i = 0; i < edge.getN(); i++) {
                    int x = edge.getX(i);
                    int y = edge.getY(i);
                    
                    if (shiftX) {
                        x += bestShift;
                    } else {
                        y += bestShift;
                    }
                    
                    edge.set(i, x, y);
                }
            }
        }
    }

    private void applyShift(List<PairIntArray> skyEdges, int shift, 
        boolean shiftX) {
        
        for (PairIntArray edge : skyEdges) {
            for (int i = 0; i < edge.getN(); i++) {
                int x = edge.getX(i);
                int y = edge.getY(i);

                if (shiftX) {
                    x += shift;
                } else {
                    y += shift;
                }

                edge.set(i, x, y);
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
   
    public GreyscaleImage getTheta() {
        return theta;
    }
    
    public GreyscaleImage getGradientXY() {
        return gradientXY;
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
