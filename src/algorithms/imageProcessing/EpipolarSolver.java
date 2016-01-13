package algorithms.imageProcessing;

import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.ejml.simple.SimpleMatrix;

/**
 * class encapsulating the steps from scale calculation to matching corners to
 * making correspondence lists to solving for epipolar projection.
 *
 * @author nichole
 */
public class EpipolarSolver {

    private final ImageExt img1;
    private final ImageExt img2;

    private enum State {
        INITIALIZED, SOLVED, NO_SOLUTION
    }
    private State state = null;

    private PairFloatArray solutionLeftXY = null;
    private PairFloatArray solutionRightXY = null;

    private final FeatureMatcherSettings featureSettings;

    private Logger log = Logger.getLogger(this.getClass().getName());

    public EpipolarSolver(ImageExt image1, ImageExt image2,
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
        EuclideanSegmentFeatureMatcher2 wrapper
            = new EuclideanSegmentFeatureMatcher2(img1, img2, featureSettings);

        boolean solved = wrapper.match();

        if (!solved) {
            state = State.NO_SOLUTION;
            return null;
        }

        List<PairInt> points1, points2;

        if (wrapper.getSolutionMatched1().size() < 7) {

            EuclideanSegmentFeatureMatcher wrapper2
                = new EuclideanSegmentFeatureMatcher(img1, img2, featureSettings);

            CorrespondenceList cl = wrapper2.matchFeatures();

            if (cl == null) {
                state = State.NO_SOLUTION;
                return null;
            }

            log.info("params from scale calc: scale=" + cl.getScale()
                + " rot(deg)=" + cl.getRotationInDegrees()
                + " tx=" + cl.getTranslationX() + " ty=" + cl.getTranslationY());

            points1 = cl.getPoints1();
            points2 = cl.getPoints2();

        } else {

            points1 = wrapper.getSolutionMatched1();
            points2 = wrapper.getSolutionMatched2();
            TransformationParameters params = wrapper.getSolutionTransformationParameters();

            log.info("params from scale calc: scale=" + params.getScale()
                + " rot(deg)=" + params.getRotationInDegrees()
                + " tx=" + params.getTranslationX() + " ty=" + params.getTranslationY());
        }

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

        EpipolarTransformationFit fit = solver.calculateEpipolarProjection(
            matchedLeftXY, matchedRightXY, outputLeftXY, outputRightXY);

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
        
        Image img1Cp = img1.copyImage();
        Image img2Cp = img2.copyImage();

        SimpleMatrix input1
            = StereoProjectionTransformer.rewriteInto3ColumnMatrix(solutionLeftXY);

        SimpleMatrix input2
            = StereoProjectionTransformer.rewriteInto3ColumnMatrix(solutionRightXY);

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

        StereoProjectionTransformer spTransformer = new StereoProjectionTransformer();

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

        for (int ii = 0; ii < input1.numCols(); ii++) {
            double x = input1.get(0, ii);
            double y = input1.get(1, ii);
            ImageIOHelper.addPointToImage((float) x, (float) y, img1, 3,
                255, 0, 0);
        }
        for (int ii = 0; ii < input2.numCols(); ii++) {
            double x2 = input2.get(0, ii);
            double y2 = input2.get(1, ii);
            ImageIOHelper.addPointToImage((float) x2, (float) y2, img2, 3,
                255, 0, 0);
        }

        String dirPath;
        try {
            dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.writeOutputImage(
                dirPath + "/tmp_m_1_" + featureSettings.getDebugTag() + ".png", img1);

            ImageIOHelper.writeOutputImage(
                dirPath + "/tmp_m_2_" + featureSettings.getDebugTag() + ".png", img2);

        } catch (IOException ex) {
            Logger.getLogger(EpipolarSolver.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}
