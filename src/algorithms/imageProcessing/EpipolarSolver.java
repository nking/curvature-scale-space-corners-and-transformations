package algorithms.imageProcessing;

import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.ejml.simple.SimpleMatrix;

/**
 * class encapsulating the steps from scale calculation to matching corners
 * to making correspondence lists to solving for epipolar projection.
 * 
 * @author nichole
 */
public class EpipolarSolver {
    
    private final ImageExt img1;
    private final ImageExt img2;
    
    private GreyscaleImage gsImg1 = null;
    private GreyscaleImage gsImg2 = null;
        
    private Set<CornerRegion> cornerRegions1 = null;
    private Set<CornerRegion> cornerRegions2 = null;

    private enum State {
        INITIALIZED, SOLVED, NO_SOLUTION
    }
    private State state = null;
    
    private final boolean doDetermineScale;
    
    private final boolean debug;
    
    private final String debugTagPrefix;
    
    private TransformationParameters params = null;
            
    private float scaleTol = 0.2f;
    
    private float rotationInRadiansTol = (float)(20. * Math.PI/180.);
    
    //TODO: revise this...
    private int transXYTol = 20;
        
    private PairFloatArray solutionLeftXY = null;
    private PairFloatArray solutionRightXY = null;
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public EpipolarSolver(ImageExt image1, ImageExt image2) {
        img1 = image1;
        img2 = image2;
        doDetermineScale = true;
        debug = false;
        debugTagPrefix = "";
        state = State.INITIALIZED;
    }
    
    public EpipolarSolver(ImageExt image1, ImageExt image2, 
        String debugTagPrefix) {
        img1 = image1;
        img2 = image2;
        doDetermineScale = true;
        debug = true;
        this.debugTagPrefix = debugTagPrefix;
        state = State.INITIALIZED;
    }
    
    /**
     * constructor accepting transformation parameters.  Note, for best results,
     * the standard deviations within parameters should be populated because they
     * are used as tolerances in matching.
     * @param image1
     * @param image2
     * @param parameters 
     */
    public EpipolarSolver(ImageExt image1, ImageExt image2, 
        TransformationParameters parameters) {
        img1 = image1;
        img2 = image2;
        doDetermineScale = false;
        params = parameters;
        debug = false;
        debugTagPrefix = "";
        state = State.INITIALIZED;
    }
    
    /**
     * constructor accepting transformation parameters and a debugging tag for
     * image names.  Note, for best results, the standard deviations within 
     * parameters should be populated because they are used as tolerances in 
     * matching.
     * @param image1
     * @param image2
     * @param parameters
     * @param debugTagPrefix 
     */
    public EpipolarSolver(ImageExt image1, ImageExt image2, 
        TransformationParameters parameters, String debugTagPrefix) {
        img1 = image1;
        img2 = image2;
        doDetermineScale = false;
        params = parameters;
        debug = true;
        this.debugTagPrefix = debugTagPrefix;
        state = State.INITIALIZED;
    }
    
    public StereoProjectionTransformerFit solve() throws IOException, 
        NoSuchAlgorithmException {
        
        if (!state.equals(State.INITIALIZED)) {
            throw new IllegalStateException("solve() has already been invoked");
        }
        
        FeatureMatcherWrapper wrapper = null;
        
        if (params != null) {
            if (debug) {
                wrapper = new FeatureMatcherWrapper(img1, img2, params, debugTagPrefix);
            } else {
                wrapper = new FeatureMatcherWrapper(img1, img2, params);
            }
        } else {
            if (debug) {
                wrapper = new FeatureMatcherWrapper(img1, img2, debugTagPrefix);
            } else {
                wrapper = new FeatureMatcherWrapper(img1, img2);
            }
        }
        
        CorrespondenceList cl = wrapper.matchFeatures();
        
        if (cl == null) {
            state = State.NO_SOLUTION;
            return null;
        }
        
        List<PairInt> points1 = cl.getPoints1();
        List<PairInt> points2 = cl.getPoints2();
        
        int n = points1.size();
        
        PairFloatArray matchedLeftXY = new PairFloatArray(n);
        PairFloatArray matchedRightXY = new PairFloatArray(n);
        PairFloatArray outputLeftXY = new PairFloatArray(n);
        PairFloatArray outputRightXY = new PairFloatArray(n);
        
        for (int i = 0; i < n; ++i) {
            matchedLeftXY.add(points1.get(i).getX(), points1.get(i).getY());
            matchedRightXY.add(points2.get(i).getX(), points2.get(i).getY());
        }
        
        RANSACSolver solver = new RANSACSolver();
        
        StereoProjectionTransformerFit fit = solver.calculateEpipolarProjection(
            matchedLeftXY, matchedRightXY,
            outputLeftXY, outputRightXY);
        
        if (fit != null) {
            
            //TODO: check that stdev is reasonable
            
            state = State.SOLVED;
            
            this.solutionLeftXY = outputLeftXY;
            
            this.solutionRightXY = outputRightXY;
            
            if (debug) {
                plotFit(fit);
            }
            
            return fit;
            
        } else {
                        
            state = State.NO_SOLUTION;
            
            return null;
        }
    }
    
    private void plotFit(StereoProjectionTransformerFit fit) {
       
        Image img1Cp = img1.copyImage();
        Image img2Cp = img2.copyImage();
        
        SimpleMatrix input1 =
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(solutionLeftXY);

        SimpleMatrix input2 =
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(solutionRightXY);

        for (int ii = 0; ii < input1.numCols(); ii++) {
            double x = input1.get(0, ii);
            double y = input1.get(1, ii);
            ImageIOHelper.addPointToImage((float) x, (float) y, img1Cp, 3,
                255, 0, 0);
        }
        for (int ii = 0; ii < input2.numCols(); ii++) {
            double x2 = input2.get(0, ii);
            double y2 = input2.get(1, ii);
            ImageIOHelper.addPointToImage((float) x2, (float) y2, img2Cp, 3,
                255, 0, 0);
        }

        StereoProjectionTransformer spTransformer = new
            StereoProjectionTransformer();

        for (int ii = 0; ii < input2.numCols(); ii++) {
            
            int[] rgb = ImageIOHelper.getNextRGB(ii);
                        
            SimpleMatrix epipolarLinesInLeft = fit.getFundamentalMatrix().transpose().mult(input2);
            
            PairIntArray leftLine = spTransformer.getEpipolarLine(
                epipolarLinesInLeft, img1Cp.getWidth(), img1Cp.getHeight(), ii);
            
            ImageIOHelper.addCurveToImage(leftLine, img1, 0, rgb[0], rgb[1], rgb[2]);
        }

        for (int ii = 0; ii < input1.numCols(); ii++) {
            
            int[] rgb = ImageIOHelper.getNextRGB(ii);
            
            SimpleMatrix epipolarLinesInRight = fit.getFundamentalMatrix().mult(input1);
            
            PairIntArray rightLine = spTransformer.getEpipolarLine(
                epipolarLinesInRight, img2.getWidth(), img2.getHeight(), ii);
            
            ImageIOHelper.addCurveToImage(rightLine, img2, 0, 
                rgb[0], rgb[1], rgb[2]);
        }

        String dirPath;
        try {
            dirPath = ResourceFinder.findDirectory("bin");
            
            ImageIOHelper.writeOutputImage(
                dirPath + "/tmp_m_1_" + debugTagPrefix + ".png", img1);
        
            ImageIOHelper.writeOutputImage(
                dirPath + "/tmp_m_2_" + debugTagPrefix + ".png", img2);
        
        } catch (IOException ex) {
            Logger.getLogger(EpipolarSolver.class.getName()).log(Level.SEVERE, null, ex);
        }        
    }
        
}
