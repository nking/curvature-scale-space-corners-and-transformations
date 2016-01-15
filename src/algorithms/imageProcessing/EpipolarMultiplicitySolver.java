package algorithms.imageProcessing;

import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.ejml.simple.SimpleMatrix;

/**
 * class encapsulating the steps from generating a degenerate correspondenc list
 * to solving for the best epipolar projection.
 
 * @author nichole
 */
public class EpipolarMultiplicitySolver {

    private final ImageExt img1;
    private final ImageExt img2;

    private static enum State {
        INITIALIZED, SOLVED, NO_SOLUTION
    }
    private State state = null;

    private PairIntArray solutionLeftXY = null;
    private PairIntArray solutionRightXY = null;

    private final FeatureMatcherSettings featureSettings;

    private Logger log = Logger.getLogger(this.getClass().getName());

    public EpipolarMultiplicitySolver(ImageExt image1, ImageExt image2,
        FeatureMatcherSettings settings) {

        img1 = image1;
        img2 = image2;

        featureSettings = settings.copy();

        state = State.INITIALIZED;
    }
    
    public EpipolarTransformationFit solve() throws IOException,
        NoSuchAlgorithmException {

        if (!state.equals(State.INITIALIZED)) {
            throw new IllegalStateException("solve() has already been invoked");
        }

        //NOTE: will change this to the non euclidean solver soon
        NonEuclideanSegmentFeatureMatcher2 wrapper
            = new NonEuclideanSegmentFeatureMatcher2(img1, img2, featureSettings);

        boolean solved = wrapper.match();

        if (!solved) {
            state = State.NO_SOLUTION;
            return null;
        }
         
        if (wrapper.getSolutionMatched2Multiplicity().size() < 7) {
            state = State.NO_SOLUTION;
            return null;
        }
        
        List<PairInt> matchedLeftXY = new ArrayList<PairInt>();       
        List<List<PairInt>> matchedRightXYs = new ArrayList<List<PairInt>>();
        for (int i = 0; i < wrapper.getSolutionMatched1().size(); ++i) {
            List<PairInt> list = wrapper.getSolutionMatched2Multiplicity().get(i);
            if (list.size() > 0) {
                matchedRightXYs.add(list);
                matchedLeftXY.add(wrapper.getSolutionMatched1().get(i));
            }
        }
        assert(matchedLeftXY.size() == matchedRightXYs.size());
             
        PairIntArray outputLeftXY = new PairIntArray();
        PairIntArray outputRightXY = new PairIntArray();
     
        RANSACMultiplicitySolver solver = new RANSACMultiplicitySolver();
        
        EpipolarTransformationFit fit = solver.calculateEpipolarProjection(
            matchedLeftXY, matchedRightXYs, outputLeftXY, outputRightXY);
        
        if (fit != null) {

            //TODO: check that stdev is reasonable
            state = State.SOLVED;

            this.solutionLeftXY = outputLeftXY;

            this.solutionRightXY = outputRightXY;

            if (featureSettings.debug()) {
                plotFit(fit);
            }

            return fit;

        } else {

            state = State.NO_SOLUTION;

            return null;
        }
    }

    private void plotFit(EpipolarTransformationFit fit) {

        if (this.img1 == null || this.img2 == null) {
            return;
        }
        
        Image img1Cp = img1.copyImage();
        Image img2Cp = img2.copyImage();

        SimpleMatrix input1
            = StereoProjectionTransformer.rewriteInto3ColumnMatrix(solutionLeftXY);

        SimpleMatrix input2
            = StereoProjectionTransformer.rewriteInto3ColumnMatrix(solutionRightXY);

        StereoProjectionTransformer spTransformer = new StereoProjectionTransformer();

        for (int ii = 0; ii < input2.numCols(); ii++) {

            int[] rgb = ImageIOHelper.getNextRGB(ii);

            SimpleMatrix epipolarLinesInLeft = fit.getFundamentalMatrix().transpose().mult(input2);

            PairIntArray leftLine = spTransformer.getEpipolarLine(
                epipolarLinesInLeft, img1Cp.getWidth(), img1Cp.getHeight(), ii);

            ImageIOHelper.addCurveToImage(leftLine, img1Cp, 0, rgb[0], rgb[1], rgb[2]);
        }

        for (int ii = 0; ii < input1.numCols(); ii++) {

            int[] rgb = ImageIOHelper.getNextRGB(ii);

            SimpleMatrix epipolarLinesInRight = fit.getFundamentalMatrix().mult(input1);

            PairIntArray rightLine = spTransformer.getEpipolarLine(
                epipolarLinesInRight, img2Cp.getWidth(), img2Cp.getHeight(), ii);

            ImageIOHelper.addCurveToImage(rightLine, img2Cp, 0,
                rgb[0], rgb[1], rgb[2]);
        }

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

        String dirPath;
        try {
            dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.writeOutputImage(
                dirPath + "/tmp_m_1_" + featureSettings.getDebugTag() + ".png", img1Cp);

            ImageIOHelper.writeOutputImage(
                dirPath + "/tmp_m_2_" + featureSettings.getDebugTag() + ".png", img2Cp);

        } catch (IOException ex) {
            Logger.getLogger(EpipolarMultiplicitySolver.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}
