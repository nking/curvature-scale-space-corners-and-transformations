package algorithms.imageProcessing.scaleSpace;

import algorithms.compGeometry.PerimeterFinder2;
import algorithms.connected.ConnectedPointsFinder;
import algorithms.imageProcessing.CannyEdgeFilterAdaptive;
import algorithms.imageProcessing.ClosedCurveAndJunctionFinder;
import algorithms.imageProcessing.EdgeExtractorSimple;
import algorithms.imageProcessing.EdgeFilterProducts;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.PostLineThinnerCorrections;
import algorithms.imageProcessing.SpurRemover;
import algorithms.imageProcessing.features.mser.MSEREdges;
import algorithms.misc.MiscDebug;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.List;
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
     * is present whether it is a closed curve or not can be checked.
     */
    protected List<PairIntArray> edges = new ArrayList<PairIntArray>();
            
    protected boolean useLineDrawingMode = false;
                        
    protected EdgeFilterProducts filterProducts = null;
    
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
    }
    
    public AbstractCurvatureScaleSpaceMapper(GreyscaleImage input) {
        
        img = input.copyImage();
        
        ImageProcessor imageProcessor = new ImageProcessor();

        originalImg = ImageIOHelper.convertImage(input);
        
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
        
        this.edges = new ArrayList<PairIntArray>(theEdges);
        
        state = CurvatureScaleSpaceMapperState.INITIALIZED;
    }

    public void useLineDrawingMode() {
        useLineDrawingMode = true;
    }
    
    protected void initialize() {
        
        if (state.ordinal() < 
            CurvatureScaleSpaceMapperState.INITIALIZED.ordinal()) {
            
            // extract edges
            extractEdges();
                        
            state = CurvatureScaleSpaceMapperState.INITIALIZED;
        }
    }
    
    protected void extractEdgesLineMode() {
        
        CannyEdgeFilterAdaptive filter = new CannyEdgeFilterAdaptive();

        if (useLineDrawingMode) {
            filter.setToUseLineDrawingMode();
        }
        //filter.setToDebug();
        
        filter.applyFilter(img);
        
        filterProducts = filter.getFilterProducts();
        
        GreyscaleImage gXY = filterProducts.getGradientXY();
        PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
        pltc.extremeStaircaseRemover(gXY);
        
        //GreyscaleImage tmp = gXY.copyImage();
        //tmp.multiply(255.f);
        //MiscDebug.writeImage(tmp, "_GXY_" + MiscDebug.getCurrentTimeFormatted());
        
        edges.clear();
        
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
            
            edges.add(closedCurve);
        }
                
        state = CurvatureScaleSpaceMapperState.EDGE_FILTERED;
    }
    
    protected void extractEdges() {
        if (useLineDrawingMode) {
            extractEdgesLineMode();
        } else {
            extractEdgesNotLineMode();
        }
    }
    
    protected void extractEdgesNotLineMode() {
       
        edges.clear();
        
        int w = originalImg.getWidth();
        int h = originalImg.getHeight();

        MSEREdges mserEdges = new MSEREdges(this.originalImg);
        mserEdges.mergeAndExtractEdges();
        this.filterProducts = mserEdges.getEdgeFilterProducts();
        List<TIntSet> edgeSets = mserEdges.getEdges();
        
        PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
        
        PerimeterFinder2 finder2 = new PerimeterFinder2();
        for (int i = 0; i < edgeSets.size(); ++i) {
            TIntSet set = edgeSets.get(i); 
            
            pltc.extremeStaircaseRemover(set, w, h);
            
            PairIntArray ordered = null;
            try {
                ordered = finder2.orderTheBoundary(
                    set, w, h);
            } catch (Exception e) {
                log.severe(e.getMessage());
                continue;
            }

            PairIntArrayWithColor closedCurve = new PairIntArrayWithColor(
                ordered);
            closedCurve.setAsClosedCurve();

            edges.add(closedCurve);
        }
                       
        state = CurvatureScaleSpaceMapperState.EDGES_EXTRACTED;
        
        log.fine("edges extracted");
    }
  
    public List<PairIntArray> getEdges() {
        return edges;
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
   
    public EdgeFilterProducts getEdgeFilterProducts() {
        return filterProducts;
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
