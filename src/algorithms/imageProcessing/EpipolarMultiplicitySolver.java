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
   Class accepting a map of points from image1 matched to possibly multiple
   points in image 2, and determining the best epipolar solution for the
   true matches using a modified RANSAC algorithm.
   The algorithm follows "Generalized RANSAC framework for relaxed correspondence problems"
   by Zhang and Kosecka in choosing uniformly randomly from each point1
   as a single entity then if chosen, choosing randomly from the point2's
   it is possibly matched to.
 
 * @author nichole
 */
public class EpipolarMultiplicitySolver {

    private enum State {
        INITIALIZED, SOLVED, NO_SOLUTION
    }
    private State state = null;

    private PairFloatArray solutionLeftXY = null;
    private PairFloatArray solutionRightXY = null;
    
    private final List<PairInt> inputMatches1;
    
    private final List<List<PairInt>> inputMultipleMatches2;

    private Logger log = Logger.getLogger(this.getClass().getName());

    private GreyscaleImage debugImg1 = null;
    
    private GreyscaleImage debugImg2 = null;
    
    private String debugTag = "";
    
    public EpipolarMultiplicitySolver(List<PairInt> inputMatches1,
        List<List<PairInt>> inputMultipleMatches2) {
                
        this.inputMatches1 = inputMatches1;
        
        this.inputMultipleMatches2 = inputMultipleMatches2;
        
        state = State.INITIALIZED;
    }
    
    public void setImagesForDebuggingPlots(GreyscaleImage img1, 
        GreyscaleImage img2, String debugTag) {
        
        this.debugImg1 = img1;
        this.debugImg2 = img2;
        this.debugTag = debugTag;
    }

    public EpipolarTransformationFit solve() throws IOException,
        NoSuchAlgorithmException {

        if (!state.equals(State.INITIALIZED)) {
            throw new IllegalStateException("solve() has already been invoked");
        }

        List<PairInt> points1, points2;

        throw new UnsupportedOperationException("not yet implemented");

        /*
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
solver.debugImg1 = img1;
solver.debugImg2 = img2;
solver.debugTag = debugTag;
        
        EpipolarTransformationFit fit = solver.calculateEpipolarProjection(
            matchedLeftXY, matchedRightXY, outputLeftXY, outputRightXY);

        if (fit != null) {

            //TODO: check that stdev is reasonable
            state = State.SOLVED;

            this.solutionLeftXY = outputLeftXY;

            this.solutionRightXY = outputRightXY;

            plotFit(fit);

            return fit;

        } else {

            state = State.NO_SOLUTION;

            return null;
        }
            */
    }

    private void plotFit(EpipolarTransformationFit fit) {

        if (debugImg1 == null || debugImg2 == null) {
            return;
        }
        
        Image img1Cp = debugImg1.copyToColorGreyscale();
        Image img2Cp = debugImg2.copyToColorGreyscale();

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
                dirPath + "/tmp_m_1_" + debugTag + ".png", img1Cp);

            ImageIOHelper.writeOutputImage(
                dirPath + "/tmp_m_2_" + debugTag + ".png", img2Cp);

        } catch (IOException ex) {
            Logger.getLogger(EpipolarMultiplicitySolver.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}
