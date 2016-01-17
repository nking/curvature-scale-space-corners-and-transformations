package algorithms.imageProcessing.features;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.transform.EuclideanTransformationFit;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * class encapsulating the steps from scale calculation to matching corners to
 * making correspondence lists to solving for epipolar projection.
 *
 * @author nichole
 */
public class EuclideanSolver {

    private final ImageExt img1;
    private final ImageExt img2;

    private enum State {
        INITIALIZED, SOLVED, NO_SOLUTION
    }
    private State state = null;

    private PairIntArray solutionLeftXY = null;
    private PairIntArray solutionRightXY = null;

    private final FeatureMatcherSettings featureSettings;

    private Logger log = Logger.getLogger(this.getClass().getName());

    public EuclideanSolver(ImageExt image1, ImageExt image2,
        FeatureMatcherSettings settings) {

        img1 = image1;
        img2 = image2;

        featureSettings = settings.copy();

        state = State.INITIALIZED;
    }

    public EuclideanTransformationFit solve() throws IOException,
        NoSuchAlgorithmException {

        if (!state.equals(State.INITIALIZED)) {
            throw new IllegalStateException("solve() has already been invoked");
        }

        NonEuclideanSegmentFeatureMatcher matcher = 
            new NonEuclideanSegmentFeatureMatcher(img1, img2, featureSettings);
        
        boolean solved = matcher.match();

        if (!solved || (matcher.getSolutionMatched1().size() < 7)) {
            state = State.NO_SOLUTION;
            return null;
        }

        List<PairInt> points1 = matcher.getSolutionMatched1();
        List<PairInt> points2 = matcher.getSolutionMatched2();

        int n = points1.size();

        PairIntArray matchedLeftXY = new PairIntArray(n);
        PairIntArray matchedRightXY = new PairIntArray(n);
        PairIntArray outputLeftXY = new PairIntArray(n);
        PairIntArray outputRightXY = new PairIntArray(n);

        for (int i = 0; i < n; ++i) {
            matchedLeftXY.add(points1.get(i).getX(), points1.get(i).getY());
            matchedRightXY.add(points2.get(i).getX(), points2.get(i).getY());
        }
        
        if (featureSettings.debug()) {
            plotFit(matchedLeftXY, matchedRightXY, "ransac_input");
        }

        RANSACEuclideanSolver solver = new RANSACEuclideanSolver();
        EuclideanTransformationFit fit = solver.calculateEuclideanTransformation(
            matchedLeftXY, matchedRightXY, outputLeftXY, outputRightXY);
        
        if (fit != null) {

            //TODO: check that stdev is reasonable
            state = State.SOLVED;

            this.solutionLeftXY = outputLeftXY;

            this.solutionRightXY = outputRightXY;

            if (featureSettings.debug()) {
                plotFit(outputLeftXY, outputRightXY, "tmp_m");
            }

            return fit;

        } else {

            state = State.NO_SOLUTION;

            return null;
        }
    }

    private void plotFit(PairIntArray xy1, PairIntArray xy2, String fileSuffix) {
        
        Image img1Cp = img1.copyImage();
        Image img2Cp = img2.copyImage();

        for (int ii = 0; ii < xy1.getN(); ii++) {
            double x = xy1.getX(ii);
            double y = xy1.getY(ii);
            ImageIOHelper.addPointToImage((float) x, (float) y, img1Cp, 3,
                255, 0, 0);
        }
        for (int ii = 0; ii < xy2.getN(); ii++) {
            double x = xy2.getX(ii);
            double y = xy2.getY(ii);
            ImageIOHelper.addPointToImage((float) x, (float) y, img2Cp, 3,
                255, 0, 0);
        }

        String dirPath;
        try {
            dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.writeOutputImage(
                dirPath + "/" + fileSuffix + "_1_" + featureSettings.getDebugTag() + ".png", img1Cp);

            ImageIOHelper.writeOutputImage(
                dirPath + "/" + fileSuffix + "_2_" + featureSettings.getDebugTag() + ".png", img2Cp);

        } catch (IOException ex) {
            Logger.getLogger(EuclideanSolver.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}
