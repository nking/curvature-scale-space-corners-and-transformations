package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.ResourceFinder;
import algorithms.util.PairFloatArray;
import algorithms.util.PairFloatArrayUnmodifiable;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.RangeInt;
import java.awt.Color;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.List;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;
import org.ejml.simple.*;
import static org.junit.Assert.fail;

/**
 *
 * @author nichole
 */
public class PointMatcher3Test extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public PointMatcher3Test() {
    }

    /*
    for more datasets:
    http://www.robots.ox.ac.uk/~vgg/data/data-mview.html
    */

    public void testPreSearches() throws Exception {

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1437014214666L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);

        //image size range 250 to 10000
        //  these are altered below
        int imageWidth = 650;
        int imageHeight = 400;

        float scale = 1.0f;

        boolean spacer = true;
        
        float maxOfMinDiffRots = Float.MIN_VALUE;
        float maxOfBestDiffRots = Float.MIN_VALUE;
        
        for (int nn = 0; nn < 1; ++nn) { // repeat number of tests
        for (int rotType = 0; rotType < 5; ++rotType) {
            for (int nTest = 10; nTest < 12 /*20*/; ++nTest) { // this increases nPoints

                PointMatcher pointMatcher = new PointMatcher();

                PairIntArray unmatchedLeftXY = new PairIntArray();

                float rotInDegrees = (int)(sr.nextFloat() * 20);

                if (rotType == 1) {
                    rotInDegrees = sr.nextBoolean() ? (90 + rotInDegrees) :
                        (90 - rotInDegrees);
                    //image size range 250 to 10000
                    imageWidth = Math.abs(5000) + 250;
                    imageHeight = Math.abs(5000) + 250;
                    if ((imageWidth & 1) == 1) {
                        imageWidth++;
                    }
                    if ((imageHeight & 1) == 1) {
                        imageHeight++;
                    }
                } else if (rotType == 2) {
                    rotInDegrees = sr.nextBoolean() ? (180 + rotInDegrees) :
                        (180 - rotInDegrees);
                    //image size range 250 to 10000
                    imageWidth = Math.abs(250) + 250;
                    imageHeight = Math.abs(250) + 250;
                    if ((imageWidth & 1) == 1) {
                        imageWidth++;
                    }
                    if ((imageHeight & 1) == 1) {
                        imageHeight++;
                    }
                } else if (rotType == 3) {
                    rotInDegrees = sr.nextBoolean() ? (270 + rotInDegrees) :
                        (270 - rotInDegrees);
                    //image size range 250 to 10000
                    imageWidth = Math.abs(1000) + 250;
                    imageHeight = Math.abs(1000) + 250;
                    if ((imageWidth & 1) == 1) {
                        imageWidth++;
                    }
                    if ((imageHeight & 1) == 1) {
                        imageHeight++;
                    }
                } else if (rotType == 4) {
                    rotInDegrees = (360 - rotInDegrees);
                    //image size range 250 to 10000
                    imageWidth = Math.abs(500) + 250;
                    imageHeight = Math.abs(500) + 250;
                    if ((imageWidth & 1) == 1) {
                        imageWidth++;
                    }
                    if ((imageHeight & 1) == 1) {
                        imageHeight++;
                    }
                }

                int transX = (int)(0.25f * sr.nextFloat() * (1 + sr.nextInt(imageWidth)));
                int transY = (int)(0.05f * sr.nextFloat() * imageHeight);
                if (sr.nextBoolean()) {
                    transX *= -1;
                }
                if (sr.nextBoolean()) {
                    transY *= -1;
                }

                TransformationParameters params = new TransformationParameters();
                params.setRotationInDegrees(rotInDegrees);
                params.setScale(scale);
                params.setTranslationX(transX);
                params.setTranslationY(transY);
                params.setOriginX(0);
                params.setOriginY(0);

                int nPoints = (nTest + 1) * 7;

                log.info("test for nPoints=" + nPoints 
                    + "\nparams=" + params.toString());

                for (int i = 0; i < nPoints; ++i) {
                    int x = (imageWidth/4) + sr.nextInt(imageWidth/4);
                    int y = (imageHeight/4) + sr.nextInt(imageHeight/4);
                    unmatchedLeftXY.add(x, y);
                }

                // ===  transform the right points  ======

                Transformer transformer = new Transformer();

                PairIntArray unmatchedRightXY = transformer.applyTransformation(
                    params, unmatchedLeftXY);

                // consider adding small deviations to the right

                // add small number of points in left that are not in right
                //   and vice versa
                int nAdd = (int)(sr.nextFloat()*nPoints/10);
                if (nAdd == 0) {
                    nAdd = 1;
                }
                for (int i = 0; i < nAdd; ++i) {
                    int x = sr.nextInt(imageWidth);
                    int y = sr.nextInt(imageHeight);
                    unmatchedLeftXY.add(x, y);
                }
                for (int i = 0; i < nAdd; ++i) {
                    int x = sr.nextInt(imageWidth);
                    int y = sr.nextInt(imageHeight);
                    unmatchedRightXY.add(x, y);
                }

                // TODO: consider scrambling the order of the right points

                // --- TODO: in the difference between the left and right regions,
                //     need to generate points in the right

                float setsFractionOfImage = 1.0f;

                int nMaxMatchable = nPoints;

                int nExpected = nMaxMatchable;
                int nEps = (int)Math.round(Math.sqrt(nMaxMatchable)/2.);
                   
                boolean useGreedyMatching = true;
                
                // ------- assert that preSearch0 gets the answer within <> degrees of rotation -------
                TransformationPointFit[] fits =
                    pointMatcher.preSearch0(
                    unmatchedLeftXY, unmatchedRightXY, scale,
                    useGreedyMatching, setsFractionOfImage);
                
                assert(fits != null);
                
                log.info("preSearch0 FITs:");
                
                float minRotationDiff = Float.MAX_VALUE;

                for (int i = 0; i < fits.length; ++i) {

                    TransformationPointFit fit = fits[i];
                    
                    if (fit == null) {
                        continue;
                    }
                    TransformationParameters fitParams = fit.getParameters();
                    int diffN = Math.abs(nExpected - fit.getNumberOfMatchedPoints());
                    float diffRotDeg = getAngleDifference(
                        fitParams.getRotationInDegrees(), rotInDegrees);
                    float diffScale = Math.abs(fitParams.getScale() - scale);
                    float diffTransX = Math.abs(fitParams.getTranslationX() - transX);
                    float diffTransY = Math.abs(fitParams.getTranslationY() - transY);

                    String space = spacer ? "    " : "";
                    log.info(space + "diff =" +
                        String.format(
                        " dRotDeg=%f, dScale=%f, dTransX=%f, dTransY=%f  nEps=%d dNPoints=%d meanDiffModel=%f",
                        diffRotDeg, diffScale, diffTransX, 
                        diffTransY, nEps, diffN, (float)fit.getMeanDistFromModel()));
                    
                    if (diffRotDeg < 0) {
                        diffRotDeg *= -1;
                    }
                    if (diffRotDeg < minRotationDiff) {
                        minRotationDiff = diffRotDeg;
                    }
                    
                    if (i == 0) {
                        if (diffRotDeg > maxOfBestDiffRots) {
                            maxOfBestDiffRots = diffRotDeg;
                        }
                    }
                }
                spacer = !spacer;
                log.info("ps0: " + nPoints 
                    + " points min difference from expected rotation=" + minRotationDiff);
                //assertTrue(minRotationDiff <= 10);
                
                if (minRotationDiff > maxOfMinDiffRots) {
                    maxOfMinDiffRots = minRotationDiff;
                }
                       
                /*try simplex, else grid search within 90 of best fit
                    
                if (true) {
                    continue;
                }
                */
                
                // ------ assert preSearch1 ---------
                TransformationPointFit[] fits2 = pointMatcher.preSearch1(
                    fits, unmatchedLeftXY, unmatchedRightXY,
                    useGreedyMatching, setsFractionOfImage);
                
                assert(fits2 != null);
                
                log.info("preSearch1 FITs:");
                
                minRotationDiff = Float.MAX_VALUE;

                TransformationPointFit fit = fits2[0];

                if (fit == null) {
                    continue;
                }
                TransformationParameters fitParams = fit.getParameters();
                int diffN = Math.abs(nExpected - fit.getNumberOfMatchedPoints());
                float diffRotDeg = getAngleDifference(
                    fitParams.getRotationInDegrees(), rotInDegrees);
                float diffScale = Math.abs(fitParams.getScale() - scale);
                float diffTransX = Math.abs(fitParams.getTranslationX() - transX);
                float diffTransY = Math.abs(fitParams.getTranslationY() - transY);

                String space = spacer ? "    *" : "*";
                log.info(space + "diff =" +
                    String.format(
                    " dRotDeg=%f, dScale=%f, dTransX=%f, dTransY=%f  nEps=%d dNPoints=%d meanDiffModel=%f",
                    diffRotDeg, diffScale, diffTransX, 
                    diffTransY, nEps, diffN, (float)fit.getMeanDistFromModel()));

                if (diffRotDeg < 0) {
                    diffRotDeg *= -1;
                }
                spacer = !spacer;
                //assertTrue(minRotationDiff <= 10);
            }
        }
        }
        
        log.info("Largest difference in rotation from expected=" + maxOfMinDiffRots
            + "  max of best fits diff from expected rotation=" + maxOfBestDiffRots);
    }

    public void estCalculateTranslationFromGridThenDownhillSimplex()
        throws Exception {

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1436162255215L;
        //seed = 1436200575218L;
        //seed = 1436226579115L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);

        //image size range 250 to 10000
        //  these are altered below
        int imageWidth = 650;
        int imageHeight = 400;

        // ---- random testing for stereo imaging ----

        float scale = 1.0f;

        int halfRange = 4;

        for (int nRuns = 0; nRuns < 10; ++nRuns) { // this increases the number of tests
            for (int rotType = 0; rotType < 4; ++rotType) {
                for (int nTest = 0; nTest < 20; ++nTest) { // this increases nPoints

                    PointMatcher pointMatcher = new PointMatcher();

                    PairIntArray unmatchedLeftXY = new PairIntArray();

                    float rotInDegrees = (int)(sr.nextFloat() * 20);

                    if (rotType == 1) {
                        rotInDegrees = sr.nextBoolean() ? (270 + rotInDegrees) :
                            (270 - rotInDegrees);
                        halfRange = 10;
                        //image size range 250 to 10000
                        imageWidth = Math.abs(250) + 250;
                        imageHeight = Math.abs(250) + 250;
                        if ((imageWidth & 1) == 1) {
                            imageWidth++;
                        }
                        if ((imageHeight & 1) == 1) {
                            imageHeight++;
                        }
                    } else if (rotType == 2) {
                        rotInDegrees = sr.nextBoolean() ? (180 + rotInDegrees) :
                            (180 - rotInDegrees);
                        halfRange = 20;
                        //image size range 250 to 10000
                        imageWidth = Math.abs(1000) + 250;
                        imageHeight = Math.abs(1000) + 250;
                        if ((imageWidth & 1) == 1) {
                            imageWidth++;
                        }
                        if ((imageHeight & 1) == 1) {
                            imageHeight++;
                        }
                    } else if (rotType == 3) {
                        rotInDegrees = sr.nextBoolean() ? (90 + rotInDegrees) :
                            (90 - rotInDegrees);
                        halfRange = 100;
                        //image size range 250 to 10000
                        imageWidth = Math.abs(2000) + 250;
                        imageHeight = Math.abs(2000) + 250;
                        if ((imageWidth & 1) == 1) {
                            imageWidth++;
                        }
                        if ((imageHeight & 1) == 1) {
                            imageHeight++;
                        }
                    } else if (rotType == 4) {
                        rotInDegrees = sr.nextBoolean() ? (45 + rotInDegrees) :
                            (45 - rotInDegrees);
                        halfRange = 115;
                        //image size range 250 to 10000
                        imageWidth = Math.abs(5000) + 250;
                        imageHeight = Math.abs(5000) + 250;
                        if ((imageWidth & 1) == 1) {
                            imageWidth++;
                        }
                        if ((imageHeight & 1) == 1) {
                            imageHeight++;
                        }
                    }

                    int deltaTransX = 15;
                    int deltaTransY = 15;

                    int transX = (int)(0.25f * sr.nextFloat() * (1 + sr.nextInt(imageWidth)));
                    int transY = (int)(0.05f * sr.nextFloat() * imageHeight);
                    if (sr.nextBoolean()) {
                        transX *= -1;
                    }
                    if (sr.nextBoolean()) {
                        transY *= -1;
                    }

                    TransformationParameters params = new TransformationParameters();
                    params.setRotationInDegrees(rotInDegrees);
                    params.setScale(scale);
                    params.setTranslationX(transX);
                    params.setTranslationY(transY);

                    int nPoints = (nTest + 1) * 7;

                    log.info("\ntest for nPoints=" + nPoints + "\nparams=" + params.toString()
                        + "\ndeltaTransX=" + deltaTransX + " deltaTransY=" + deltaTransY
                        + " translation range=" + (2*halfRange));

                    for (int i = 0; i < nPoints; ++i) {
                        int x = (imageWidth/4) + sr.nextInt(imageWidth/4);
                        int y = (imageHeight/4) + sr.nextInt(imageHeight/4);
                        unmatchedLeftXY.add(x, y);
                    }

                    // ===  transform the right points  ======

                    Transformer transformer = new Transformer();

                    PairIntArray unmatchedRightXY = transformer.applyTransformation(
                        params, unmatchedLeftXY, imageWidth >> 1, imageHeight >> 1);

                    // consider adding small deviations to the right

                    // add small number of points in left that are not in right
                    //   and vice versa
                    int nAdd = (int)(sr.nextFloat()*nPoints/10);
                    if (nAdd == 0) {
                        nAdd = 1;
                    }
                    for (int i = 0; i < nAdd; ++i) {
                        int x = sr.nextInt(imageWidth);
                        int y = sr.nextInt(imageHeight);
                        unmatchedLeftXY.add(x, y);
                    }
                    for (int i = 0; i < nAdd; ++i) {
                        int x = sr.nextInt(imageWidth);
                        int y = sr.nextInt(imageHeight);
                        unmatchedRightXY.add(x, y);
                    }

                    // change the order of points to make sure not biased
                    // by order
                    unmatchedRightXY.reverse();

                    // tolerance factor set to 0.5 of cell size makes perfect
                    // matching easier.  a larger tolerance may lead to a
                    // small number of false matches
                    float tolTransX = pointMatcher.getTolFactor() * deltaTransX;
                    float tolTransY = pointMatcher.getTolFactor() * deltaTransY;
                    if (tolTransX < 1) {
                        tolTransX = 1;
                    }
                    if (tolTransY < 1) {
                        tolTransY = 1;
                    }

                    boolean setsAreMatched = false;
                    float setsFractionOfImage = 1.0f;

                    int nMaxMatchable = nPoints;

                    RangeInt transXStartStop = new RangeInt(
                        (int)(params.getTranslationX() - halfRange),
                        (int)(params.getTranslationX() + halfRange));
                    RangeInt transYStartStop = new RangeInt(
                        (int)(params.getTranslationY() - halfRange),
                        (int)(params.getTranslationY() + halfRange));

                    PairFloatArrayUnmodifiable scaledRotatedLeft =
                        pointMatcher.scaleAndRotate(unmatchedLeftXY,
                            params.getRotationInRadians(), params.getScale(),
                             imageWidth >> 1, imageHeight >> 1);

                    // quick test of the evaluation function for perfect solution
                    // this method uses "optimal" matching
                    TransformationPointFit checkFit =
                        //pointMatcher.evaluateFitForUnMatchedOptimal(
                        pointMatcher.evaluateFitForUnMatchedGreedy(
                        scaledRotatedLeft, transX, transY, tolTransX, tolTransY,
                        unmatchedRightXY,
                        params.getScale(), params.getRotationInRadians());

                    /*
                    one can decrease the tolerance in translation error when
                    matching, or increase the tolerance in numbers when comparing
                    matched amounts.
                    */
                    int nExpected = nMaxMatchable;
                    double densX = ((double)nPoints/(double)imageWidth);
                    double densY = ((double)nPoints/(double)imageHeight);
                    int nEps = (int)Math.round(Math.sqrt(nMaxMatchable)/2.);

                    /*if ((densX > 0.15) && (densY > 0.15)) {
                        nEps = (int)Math.round(Math.sqrt(nMaxMatchable)/2.);
                    }*/

                    log.info("point density  n/width=" + densX + " n/height=" + densY);

                    log.info("check of evalFit, checkFit=" + checkFit.toString());
                    /*
                    assertTrue(
                        Math.abs(checkFit.getParameters().getRotationInRadians()
                        - params.getRotationInRadians()) < 0.1);
                    assertTrue(
                        Math.abs(checkFit.getParameters().getScale()
                        - params.getScale()) < 0.1);
                    assertTrue(
                        Math.abs(checkFit.getParameters().getTranslationX()
                        - params.getTranslationX()) < 0.1);
                    assertTrue(
                        Math.abs(checkFit.getParameters().getTranslationY()
                        - params.getTranslationY()) < 0.1);
                    assertTrue(
                        Math.abs(checkFit.getNumberOfMatchedPoints() -
                        nMaxMatchable) <= nEps);
                    assertTrue(checkFit.getMeanDistFromModel() < 1);
                    assertTrue(checkFit.getStDevFromMean() < 0.5);
                    */

                    int nIter = 0;
                    int nMaxIter = 10;

                    int dsMaxIter = pointMatcher.getDsNMaxIter();
                    double dens = (double)nPoints/(double)imageWidth;

                    boolean converged = false;

                    while ((nIter == 0) || (!converged && (nIter < nMaxIter))) {

                        if (nIter > 0) {
                            break;
                        }

                        TransformationPointFit fit =
                            pointMatcher.calculateTranslationFromGridThenDownhillSimplex(
                                scaledRotatedLeft,  unmatchedLeftXY, unmatchedRightXY,
                                imageWidth, imageHeight, imageWidth, imageHeight,
                                params.getRotationInRadians(), params.getScale(),
                                transXStartStop, deltaTransX,
                                transYStartStop, deltaTransY,
                                tolTransX, tolTransY, setsAreMatched,
                                setsFractionOfImage);

         //assert(fit != null);
                        TransformationParameters fitParams = fit.getParameters();
                        int diffN = Math.abs(nExpected - fit.getNumberOfMatchedPoints());
                        float diffRotDeg = getAngleDifference(
                            fitParams.getRotationInDegrees(), rotInDegrees);
                        float diffScale = Math.abs(fitParams.getScale() - scale);
                        float diffTransX = Math.abs(fitParams.getTranslationX() - transX);
                        float diffTransY = Math.abs(fitParams.getTranslationY() - transY);

                        double epsTrans = 3;
                        double epsRot = 5;

                        log.info("FINAL FIT=" + fit.toString());
                        log.info("diff result and expected =" +
                            String.format(
                            " dNPoints=%d, dRotDeg=%f, dScale=%f, dTransX=%f, dTransY=%f  nEps=%d  tEps=%f",
                            diffN, diffRotDeg, diffScale, diffTransX, diffTransY, nEps, (float)epsTrans)
                        );

                        if ((diffN <= nEps) && (diffRotDeg <= epsRot) &&
                            (diffScale < 0.2) && (diffTransX <= epsTrans) &&
                            (diffTransY <= epsTrans)) {
                            converged = true;
                        }

                        nIter++;
                    }

                    if (nIter > 1) {
                        log.info("needed change for translation delta for"
                            + " dens=" + dens + " nPoints=" + nPoints
                            + " w=" + imageWidth + " h=" + imageHeight
                            + " deltaTransX=" + deltaTransX
                            + " deltaTransY=" + deltaTransY
                            + " dsMaxIter=" + dsMaxIter
                        );
                    } else {
                        log.info("density=" + dens + " nPoints=" + nPoints
                            + " w=" + imageWidth + " h=" + imageHeight
                            + " deltaTransX=" + deltaTransX
                            + " deltaTransY=" + deltaTransY
                            + " dsMaxIter=" + dsMaxIter
                        );
                    }

        //         assertTrue(converged);

                    /* For having rotation and scale correct already,
                    a transXDelta and transYDelta of values 2 leads
                    to convergence when using a downhill simplex with max
                    number of iterations of 50.
                    Increasing the transXDelta and transYDelta to values of
                    5 also converges.

                    Increasing the transXDelta and transYDelta to values of
                    7 results in residuals in translation of a little more than
                    3 so is beginning to be large.

                    Increasing the transXDelta and transYDelta to values of
                    10 results in residuals in translation of about 5, so that
                    is probably not always going to lead to the correct solution.

                    */
                }
            }
        }
    }

    public void estCalculateEuclideanTransformation()
        throws Exception {

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1436749125398L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);

        //image size range 250 to 10000
        //  these are altered below
        int imageWidth = 650;
        int imageHeight = 400;

        // ---- random testing for stereo imaging ----

        float scale = 1.0f;

        for (int nRuns = 0; nRuns < 1; ++nRuns) { // this increases the number of tests
            for (int rotType = 0; rotType < 2; ++rotType) {
                for (int nTest = 0; nTest < 10/*20*/; ++nTest) { // this increases nPoints

                    PointMatcher pointMatcher = new PointMatcher();

                    PairIntArray unmatchedLeftXY = new PairIntArray();

                    float rotInDegrees = (int)(sr.nextFloat() * 20);

                    if (rotType == 1) {

                        rotInDegrees = 360 - rotInDegrees;
                        imageWidth = Math.abs(500) + 250;
                        imageHeight = Math.abs(500) + 250;
                        if ((imageWidth & 1) == 1) {
                            imageWidth++;
                        }
                        if ((imageHeight & 1) == 1) {
                            imageHeight++;
                        }
                    } else if (rotType == 2) {

                        rotInDegrees = 360 - rotInDegrees;
                        imageWidth = Math.abs(5000) + 250;
                        imageHeight = Math.abs(5000) + 250;
                        if ((imageWidth & 1) == 1) {
                            imageWidth++;
                        }
                        if ((imageHeight & 1) == 1) {
                            imageHeight++;
                        }
                    }

                    int transX = (int)(0.25f * sr.nextFloat() * (1 + sr.nextInt(imageWidth)));
                    int transY = (int)(0.05f * sr.nextFloat() * imageHeight);
                    if (sr.nextBoolean()) {
                        transX *= -1;
                    }
                    if (sr.nextBoolean()) {
                        transY *= -1;
                    }

                    TransformationParameters params = new TransformationParameters();
                    params.setRotationInDegrees(rotInDegrees);
                    params.setScale(scale);
                    params.setTranslationX(transX);
                    params.setTranslationY(transY);

                    int nPoints = (nTest + 1) * 7;

                    log.info("\ntest for nPoints=" + nPoints + " nTest=" + nTest
                        + " rotType=" + rotType + " nRuns=" + nRuns
                        + "\nparams=" + params.toString()
                    );

                    for (int i = 0; i < nPoints; ++i) {
                        int x = (imageWidth/4) + sr.nextInt(imageWidth/4);
                        int y = (imageHeight/4) + sr.nextInt(imageHeight/4);
                        unmatchedLeftXY.add(x, y);
                    }

                    // ===  transform the right points  ======

                    Transformer transformer = new Transformer();

                    PairIntArray unmatchedRightXY = transformer.applyTransformation(
                        params, unmatchedLeftXY, imageWidth >> 1, imageHeight >> 1);

                    // consider adding small deviations to the right

                    // add small number of points in left that are not in right
                    //   and vice versa
                    int nAdd = (int)(sr.nextFloat()*nPoints/10);
                    if (nAdd == 0) {
                        nAdd = 1;
                    }
                    for (int i = 0; i < nAdd; ++i) {
                        int x = sr.nextInt(imageWidth);
                        int y = sr.nextInt(imageHeight);
                        unmatchedLeftXY.add(x, y);
                    }
                    for (int i = 0; i < nAdd; ++i) {
                        int x = sr.nextInt(imageWidth);
                        int y = sr.nextInt(imageHeight);
                        unmatchedRightXY.add(x, y);
                    }

                    // change the order of points to make sure not biased
                    // by order
                    unmatchedRightXY.reverse();

                    float setsFractionOfImage = 1.0f;

                    int nMaxMatchable = nPoints;

                    int nExpected = nMaxMatchable;
                    int nEps = (int)Math.round(Math.sqrt(nMaxMatchable)/2.);

                    TransformationPointFit fit =
                        pointMatcher.calculateEuclideanTransformation(
                        unmatchedLeftXY, unmatchedRightXY,
                        setsFractionOfImage);

                    assert(fit != null);

                    log.info("FIT=" + fit.toString());

                    boolean converged = pointMatcher.hasConverged(fit, nMaxMatchable);

                    TransformationParameters fitParams = fit.getParameters();
                    int diffN = Math.abs(nExpected - fit.getNumberOfMatchedPoints());
                    float diffRotDeg = getAngleDifference(
                        fitParams.getRotationInDegrees(), rotInDegrees);
                    float diffScale = Math.abs(fitParams.getScale() - scale);
                    float diffTransX = Math.abs(fitParams.getTranslationX() - transX);
                    float diffTransY = Math.abs(fitParams.getTranslationY() - transY);

                    double epsTrans = 3;
                    double epsRot = 3;

                    log.info("FINAL FIT=" + fit.toString());
                    log.info("diff result and expected =" +
                        String.format(
                        " dNPoints=%d, dRotDeg=%f, dScale=%f, dTransX=%f, dTransY=%f  nEps=%d  tEps=%f",
                        diffN, diffRotDeg, diffScale, diffTransX, diffTransY, nEps, (float)epsTrans)
                    );

                    if ((diffN <= nEps) && (diffRotDeg <= epsRot) &&
                        (diffScale < 0.2) && (diffTransX <= epsTrans) &&
                        (diffTransY <= epsTrans)) {
                        converged = true;
                    }

                    assertTrue(converged);
                }
            }
        }
    }

    private void debugPlot(PairIntArray set0,
        PairIntArray set1, String label) {

        int minX0 = MiscMath.findMin(set0.getX(), set0.getN());
        int maxX0 = MiscMath.findMax(set0.getX(), set0.getN());
        int minY0 = MiscMath.findMin(set0.getY(), set0.getN());
        int maxY0 = MiscMath.findMax(set0.getY(), set0.getN());

        int minX1 = MiscMath.findMin(set1.getX(), set1.getN());
        int maxX1 = MiscMath.findMax(set1.getX(), set1.getN());
        int minY1 = MiscMath.findMin(set1.getY(), set1.getN());
        int maxY1 = MiscMath.findMax(set1.getY(), set1.getN());

        float minX = Math.min(minX0, minX1);
        float maxX = Math.max(maxX0, maxX1);

        float minY = Math.min(minY0, minY1);
        float maxY = Math.max(maxY0, maxY1);

        try {

            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(
                minX, maxX, minY, maxY);

            plotter.addPlot(set0.getX(), set0.getY(), set1.getX(), set1.getY(), label);

            plotter.writeFile(MiscDebug.getCurrentTimeFormatted());

        } catch (IOException e) {

        }
    }

    private void debugPlot(PairFloatArrayUnmodifiable set0,
        PairIntArray set1, String label) {

        float minX0 = MiscMath.findMin(set0.getX(), set0.getN());
        float maxX0 = MiscMath.findMax(set0.getX(), set0.getN());
        float minY0 = MiscMath.findMin(set0.getY(), set0.getN());
        float maxY0 = MiscMath.findMax(set0.getY(), set0.getN());

        int minX1 = MiscMath.findMin(set1.getX(), set1.getN());
        int maxX1 = MiscMath.findMax(set1.getX(), set1.getN());
        int minY1 = MiscMath.findMin(set1.getY(), set1.getN());
        int maxY1 = MiscMath.findMax(set1.getY(), set1.getN());

        float minX = Math.min(minX0, minX1);
        float maxX = Math.max(maxX0, maxX1);

        float minY = Math.min(minY0, minY1);
        float maxY = Math.max(maxY0, maxY1);

        try {

            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(
                minX, maxX, minY, maxY);

            plotter.addPlot(set0.getX(), set0.getY(), set1.getX(), set1.getY(),
                label);

            plotter.writeFile(MiscDebug.getCurrentTimeFormatted());

        } catch (IOException e) {

        }
    }

    protected class DensityTranslationResults {
        final double density;
        final int nPoints;
        final int imageWidth;
        final int imageHeight;
        final double translationDelta;
        public DensityTranslationResults(double theDensity, int numberOfPoints,
            int theImageWidth, int theImageHeight, double theTranslationDelta) {
            density = theDensity;
            nPoints = numberOfPoints;
            imageWidth = theImageWidth;
            imageHeight = theImageHeight;
            translationDelta = theTranslationDelta;
        }
    }

    /*
    https://github.com/jesolem/PCV
    adapted from code licensed under  BSD license (2-clause "Simplified BSD License").
    private SimpleMatrix calculateCameraMatrixFromFundamentalMatrix(SimpleMatrix fm) {
        SimpleMatrix u = fm.svd().getU();
        SimpleMatrix leftE = u.extractVector(true, 2);
        leftE = leftE.divide(u.get(2, 2));

        double[][] skewSymmetric = new double[3][];
        skewSymmetric[0] = new double[]{0, -leftE.get(0, 2), leftE.get(0, 1)};
        skewSymmetric[1] = new double[]{leftE.get(0, 2), 0, -leftE.get(0, 0)};
        skewSymmetric[2] = new double[]{-leftE.get(0, 1), leftE.get(0, 0), 0};

        SimpleMatrix sMatrix = new SimpleMatrix(skewSymmetric);
        double v = sMatrix.dot( fm.transpose() );


        //essential matrix from F: E = transpose(K2)*F*K1
        SimpleMatrix k = DataForTests.getVenturiCameraIntrinsics()
            .mult(DataForTests.getVenturiRotationMatrixForImage001());
        SimpleMatrix k2 = DataForTests.getVenturiCameraIntrinsics()
            .mult(DataForTests.getVenturiRotationMatrixForImage010());

        SimpleMatrix essentialMatrix = k2.transpose().mult(fm).mult(k);

        int z = 1;
    }*/

    public void estSkyline() throws Exception {

        String[] fileNames = new String[] {
            "brown_lowe_2003_image1.jpg",
            //"brown_lowe_2003_image1_rot.jpg",
            //"brown_lowe_2003_image2.jpg",
            /*"venturi_mountain_j6_0001.png",
            //"venturi_mountain_j6_0010.png",
            "seattle.jpg",
            "arches.jpg",
            "stinson_beach.jpg",
            "cloudy_san_jose.jpg",
            "stonehenge.jpg",
            "norwegian_mtn_range.jpg",
            "halfdome.jpg",
            "costa_rica.jpg",
            "new-mexico-sunrise_w725_h490.jpg",
            "arizona-sunrise-1342919937GHz.jpg",*/
            //"sky_with_rainbow.jpg",
            //"sky_with_rainbow2.jpg",
            //"patagonia_snowy_foreground.jpg",
            //"mt_rainier_snowy_field.jpg"
            //"klein_matterhorn_snowy_foreground.jpg"
            //"30.jpg",
            //"arches_sun_01.jpg",
            //"stlouis_arch.jpg",
            //"contrail.jpg"
        };

        for (String fileName : fileNames) {

            log.info("fileName=" + fileName);

            // revisit infl points.  is there a threshold removing points?
            String filePath1 = ResourceFinder.findFileInTestResources(fileName);
            ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
            int image1Width = img1.getWidth();
            int image1Height = img1.getHeight();

      /*
      ImageProcessor ip = new ImageProcessor();
      ip.testFilter(img1);
      if (true) {
          return;
      }
      */

            List<PairIntArray> edges1 = null;
            PairIntArray points1 = null;

            int nPreferredCorners = 100;
            int nCrit = 500;

            SkylineExtractor.setDebugName(fileName);

            CurvatureScaleSpaceCornerDetector detector = new
                CurvatureScaleSpaceCornerDetector(img1);

            detector.useOutdoorModeAndExtractSkyline();
            //detector.findCornersIteratively(nPreferredCorners, nCrit);
            SkylineExtractor.setDebugName(fileName);
            detector.findCorners();
            /*edges1 = detector.getEdgesInOriginalReferenceFrame();
            points1 = detector.getCornersInOriginalReferenceFrame();

            Image image1 = ImageIOHelper.readImageAsGrayScale(filePath1);

            for (PairIntArray edge : edges1) {
                ImageIOHelper.addCurveToImage(edge, image1, 2,
                    Color.YELLOW.getRed(), Color.YELLOW.getGreen(),
                    Color.YELLOW.getBlue());
            }

            ImageIOHelper.addCurveToImage(points1, image1, 1, 255, 0, 0);

            String dirPath = ResourceFinder.findDirectory("bin");
            String outFilePath = dirPath + "/tmp1_edges_and_corners_" +
                fileName + ".png";

            ImageIOHelper.writeOutputImage(outFilePath, image1);*/
        }
    }

    private void examineIterativeCorners() throws Exception {


        String fileName1 = "brown_lowe_2003_image1.jpg";
        String fileName2 = "brown_lowe_2003_image2.jpg";

        /*
        String fileName1 = "venturi_mountain_j6_0001.png";
        String fileName2 = "venturi_mountain_j6_0010.png";
        */

        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        int image1Width = img1.getWidth();
        int image1Height = img1.getHeight();

        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
        int image2Width = img2.getWidth();
        int image2Height = img2.getHeight();

        //List<PairIntArray> edges1 = null;
        List<PairIntArray> edgesSkyline1 = null;
        //List<PairIntArray> edges2 = null;
        List<PairIntArray> edgesSkyline2 = null;
        PairIntArray points1 = null;
        PairIntArray points2 = null;

        int nPreferredCorners = 100;
        int nCrit = 500;

        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img1);
        //detector.useOutdoorMode();
        //detector.useOutdoorModeAndExtractSkyline();
        detector.calculateSkylineCornersOnly();
        //detector.findCornersIteratively(nPreferredCorners, nCrit);
        detector.findCorners();
        //edges1 = detector.getEdgesInOriginalReferenceFrame();
        points1 = detector.getSkylineCornersInOriginalReferenceFrame();
        edgesSkyline1 = detector.getSkylineEdgesInOriginalReferenceFrame();

        detector = new CurvatureScaleSpaceCornerDetector(img2);
        //detector.useOutdoorMode();
        //detector.useOutdoorModeAndExtractSkyline();
        detector.calculateSkylineCornersOnly();
        //detector.findCornersIteratively(nPreferredCorners, nCrit);
        detector.findCorners();
        //edges2 = detector.getEdgesInOriginalReferenceFrame();
        points2 = detector.getSkylineCornersInOriginalReferenceFrame();
        edgesSkyline2 = detector.getSkylineEdgesInOriginalReferenceFrame();

        Image image1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        String dirPath = ResourceFinder.findDirectory("bin");
        Image image2 = ImageIOHelper.readImageAsGrayScale(filePath2);

        // ===== use partitioned matching =====
        PairIntArray outputMatchedScene = new PairIntArray();
        PairIntArray outputMatchedModel = new PairIntArray();

        PointMatcher pointMatcher = new PointMatcher();
        //pointMatcher.setCostToNumMatchedAndDiffFromModel();
        //pointMatcher.setCostToDiffFromModel();

        boolean useLargestToleranceForOutput = true;
        float setsFractionOfImage = 1.0f;

        TransformationPointFit trFit =
            pointMatcher.performMatching(points1, points2,
                outputMatchedScene, outputMatchedModel,
                useLargestToleranceForOutput, setsFractionOfImage);

        log.info("rough euclidean from skyline alone=" + trFit.toString());

/*
the skyline is hopefully useful in getting the transformation solution
into the local search region

now should try rough euclidean match using the transformation from
skyline as a start of solution with whole image corners and prefer the
later if better.

then make a larger output matched scene and model set of matched points to
determine the epipolar projection

(for brown & lee, the skyline-only region of overlap is small so
difficult to make a precise epipolar projection solution.  it's
better to have points of correspondence spread over as much of the
images as possible)
*/

        /*
        pointMatcher.performMatching(points1, points2,
        image1Width >> 1, image1Height >> 1,
            image2Width >> 1, image2Height >> 1,
            outputMatchedScene, outputMatchedModel, 1.0f);
        */

        if (outputMatchedScene.getN() < 7) {
            // no epipolar solution. need at least 7 points
            return;
        }

        PairFloatArray finalOutputMatchedScene = new PairFloatArray();
        PairFloatArray finalOutputMatchedModel = new PairFloatArray();

        RANSACSolver ransacSolver = new RANSACSolver();

        StereoProjectionTransformerFit sFit = ransacSolver
            .calculateEpipolarProjection(
                StereoProjectionTransformer.rewriteInto3ColumnMatrix(outputMatchedScene),
                StereoProjectionTransformer.rewriteInto3ColumnMatrix(outputMatchedModel),
                finalOutputMatchedScene, finalOutputMatchedModel);

        overplotEpipolarLines(sFit.getFundamentalMatrix(),
            outputMatchedScene.toPairFloatArray(), outputMatchedModel.toPairFloatArray(),
            ImageIOHelper.readImage(filePath1),
            ImageIOHelper.readImage(filePath2),
            image1Width,
            img1.getHeight(), img2.getWidth(), img2.getHeight());

        /*
        StereoProjectionTransformer st = new StereoProjectionTransformer();
        SimpleMatrix fm =
            st.calculateEpipolarProjectionForPerfectlyMatched(
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(outputMatchedScene),
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(outputMatchedModel));

        overplotEpipolarLines(fm,
            outputMatchedScene.toPairFloatArray(), outputMatchedModel.toPairFloatArray(),
            ImageIOHelper.readImage(filePath1),
            ImageIOHelper.readImage(filePath2),
            image1Width,
            img1.getHeight(), img2.getWidth(), img2.getHeight());
        */

        System.out.println("test done");
    }

    //@Test
    public void est155() throws Exception {

        PairIntArray scene = new PairIntArray();
        PairIntArray model = new PairIntArray();
        DataForTests.readBrownAndLoweMatches(scene, model);

        String fileName1 = "brown_lowe_2003_image1.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        int xSceneCentroid = img1.getWidth() >> 1;
        int ySceneCentroid = img1.getHeight() >> 1;
        String fileName2 = "brown_lowe_2003_image2.jpg";
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        int xModelCentroid = img2.getWidth() >> 1;
        int yModelCentroid = img2.getHeight() >> 1;

        StereoProjectionTransformer st = new StereoProjectionTransformer();

        SimpleMatrix left =
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(scene);
        SimpleMatrix right =
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(model);

        SimpleMatrix fm =
            st.calculateEpipolarProjectionForPerfectlyMatched(left, right);

        overplotEpipolarLines(fm, scene.toPairFloatArray(),
            model.toPairFloatArray(), img1, img2, img1.getWidth(),
            img1.getHeight(), img2.getWidth(), img2.getHeight());

    }

    //@Test
    public void est156() throws Exception {

        PairIntArray scene = new PairIntArray();
        PairIntArray model = new PairIntArray();
        DataForTests.readBrownAndLoweCorners(scene, model);

        String fileName1 = "brown_lowe_2003_image1.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        int xSceneCentroid = img1.getWidth() >> 1;
        int ySceneCentroid = img1.getHeight() >> 1;
        String fileName2 = "brown_lowe_2003_image2.jpg";
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        int xModelCentroid = img2.getWidth() >> 1;
        int yModelCentroid = img2.getHeight() >> 1;

        StereoProjectionTransformer st = new StereoProjectionTransformer();

        PairIntArray outputMatchedScene = new PairIntArray();
        PairIntArray outputMatchedModel = new PairIntArray();

        StereoProjectionTransformerFit fit =
            st.calculateEpipolarProjectionForUnmatched(scene, model,
            xSceneCentroid, ySceneCentroid,
            xModelCentroid, yModelCentroid,
            outputMatchedScene, outputMatchedModel);

        overplotEpipolarLines(fit.getFundamentalMatrix(),
            scene.toPairFloatArray(), model.toPairFloatArray(),
            img1, img2, img1.getWidth(),
            img1.getHeight(), img2.getWidth(), img2.getHeight());
    }

    private void overplotEpipolarLines(SimpleMatrix fm, PairFloatArray set1,
        PairFloatArray set2, Image img1, Image img2, int image1Width,
        int image1Height, int image2Width, int image2Height) throws IOException {

        String flNum = "";

        overplotEpipolarLines(fm, set1, set2, img1, img2, image1Width,
            image1Height, image2Width, image2Height, flNum);
    }

    private void overplotEpipolarLines(SimpleMatrix fm, PairFloatArray set1,
        PairFloatArray set2, Image img1, Image img2, int image1Width,
        int image1Height, int image2Width, int image2Height, String outfileNumber)
        throws IOException {

        SimpleMatrix input1 =
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(set1);

        SimpleMatrix input2 =
            StereoProjectionTransformer.rewriteInto3ColumnMatrix(set2);

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

        StereoProjectionTransformer spTransformer = new
            StereoProjectionTransformer();

        Color clr = null;
        for (int ii = 0; ii < input2.numCols(); ii++) {
            clr = getColor(clr);
            SimpleMatrix epipolarLinesInLeft = fm.transpose().mult(input2);
            PairIntArray leftLine = spTransformer.getEpipolarLine(
                epipolarLinesInLeft, image1Width, image1Height, ii);
            ImageIOHelper.addCurveToImage(leftLine, img1, 0,
                clr.getRed(), clr.getGreen(), clr.getBlue());
        }

        clr = null;
        for (int ii = 0; ii < input1.numCols(); ii++) {
            clr = getColor(clr);
            SimpleMatrix epipolarLinesInRight = fm.mult(input1);
            PairIntArray rightLine = spTransformer.getEpipolarLine(
                epipolarLinesInRight, img2.getWidth(), img2.getHeight(), ii);
            ImageIOHelper.addCurveToImage(rightLine, img2, 0,
                clr.getRed(), clr.getGreen(), clr.getBlue());
        }

        String dirPath = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(
            dirPath + "/tmp_m_1_" + outfileNumber + ".png", img1);
        ImageIOHelper.writeOutputImage(
            dirPath + "/tmp_m_2_" + outfileNumber + ".png", img2);
    }

    public static void main(String[] args) {

        try {
            PointMatcher3Test test = new PointMatcher3Test();

            test.testPreSearches();
            //test.testCalculateTranslationFromGridThenDownhillSimplex();
            //test.testCalculateTransformationWithGridSearch();

            //test.testCalculateEuclideanTransformation();
            //test.testPreSearch();

            /*
            tests for :
            -- for same set w/ projection
            -- for same set w/ projection and noise
            tests for scales which are close to 1 and less than 2
            */

        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            fail(e.getMessage());
        }
    }

    private Color getColor(Color clr) {
        if ((clr == null) || clr.equals(Color.MAGENTA)) {
            return Color.BLUE;
        }
        if (clr.equals(Color.BLUE)) {
            return Color.PINK;
        } else if (clr.equals(Color.PINK)) {
            return Color.GREEN;
        } else if (clr.equals(Color.GREEN)) {
            return Color.RED;
        } else if (clr.equals(Color.RED)) {
            return Color.CYAN;
        } else if (clr.equals(Color.CYAN)) {
            return Color.MAGENTA;
        } else if (clr.equals(Color.MAGENTA)) {
            return Color.LIGHT_GRAY;
        } else {
            return Color.ORANGE;
        }
    }

    private float getAngleDifference(float rotDegrees0, float rotDegrees1) {
         /*
         I  |  0
        ---------
         II | III
        */
        int q0 = 0;
        if (rotDegrees0 >= 270) {
            q0 = 3;
        } else if (rotDegrees0 >= 180) {
            q0 = 2;
        } else if (rotDegrees0 >= 90) {
            q0 = 1;
        }
        int q1 = 0;
        if (rotDegrees1 >= 270) {
            q1 = 3;
        } else if (rotDegrees1 >= 180) {
            q1 = 2;
        } else if (rotDegrees1 >= 90) {
            q1 = 1;
        }

        /*
         I  |  0
        ---------
         II | III
        */
        float angleDiff = -1;
        if (q0 == 0){
            if (q1 == 0) {
                if (rotDegrees0 > rotDegrees1) {
                    angleDiff = rotDegrees0 - rotDegrees1;
                } else {
                    angleDiff = rotDegrees1 - rotDegrees0;
                }
            } else if (q1 == 1) {
                angleDiff = (rotDegrees1 - rotDegrees0);
            } else if (q1 == 2) {
                float diff = rotDegrees1 - rotDegrees0;
                if (diff > 180) {
                    diff = 360 - diff;
                }
                angleDiff = diff;
            } else {
                angleDiff = Math.abs(360 - rotDegrees1 + rotDegrees0);
            }
        } else if (q0 == 1) {
            /*
             I  |  0
             ---------
             II | III
             */
            if (q1 == 0) {
                angleDiff = (rotDegrees1 - rotDegrees0);
            } else if (q1 == 1) {
                if (rotDegrees0 > rotDegrees1) {
                    angleDiff = rotDegrees0 - rotDegrees1;
                } else {
                    angleDiff = rotDegrees1 - rotDegrees0;
                }
            } else if (q1 == 2) {
                angleDiff = (rotDegrees1 - rotDegrees0);
            } else {
                float diff = rotDegrees1 - rotDegrees0;
                if (diff > 180) {
                    diff = 360 - diff;
                }
                angleDiff = diff;
            }
        } else if (q0 == 2) {
            /*
             I  |  0
             ---------
             II | III
             */
            if (q1 == 0) {
                float diff = rotDegrees1 - rotDegrees0;
                if (diff > 180) {
                    diff = 360 - diff;
                }
                angleDiff = diff;
            } else if (q1 == 1) {
                angleDiff = (rotDegrees0 - rotDegrees1);
            } else if (q1 == 2) {
                if (rotDegrees0 > rotDegrees1) {
                    angleDiff = rotDegrees0 - rotDegrees1;
                } else {
                    angleDiff = rotDegrees1 - rotDegrees0;
                }
            } else {
                angleDiff = (rotDegrees1 - rotDegrees0);
            }
        } else if (q0 == 3) {
            /*
             I  |  0
             ---------
             II | III
             */
            if (q1 == 0) {
                angleDiff = (360 - rotDegrees0 + rotDegrees1);
            } else if (q1 == 1) {
                float diff = (rotDegrees0 - rotDegrees1);
                if (diff > 180) {
                    diff = 360 - diff;
                }
                angleDiff = diff;
            } else if (q1 == 2) {
                angleDiff = (rotDegrees0 - rotDegrees1);
            } else {
                if (rotDegrees0 > rotDegrees1) {
                    angleDiff = rotDegrees0 - rotDegrees1;
                } else {
                    angleDiff = rotDegrees1 - rotDegrees0;
                }
            }
        }

        if (angleDiff > 359) {
            angleDiff = angleDiff - 360;
        }

        return angleDiff;
    }

    private void check(float rotationInDegrees, float scale,
        float transX, float transY,
        int image1CentroidX, int image1CentroidY, float rotRangeInDegrees) {

        Transformer tr = new Transformer();

        for (int i = 0; i < 3; ++i) {

            TransformationParameters params0 = new TransformationParameters();
            params0.setRotationInDegrees(rotationInDegrees - rotRangeInDegrees);
            params0.setScale(scale);
            params0.setTranslationX(transX);
            params0.setTranslationY(transY);

            TransformationParameters params = new TransformationParameters();
            params.setRotationInDegrees(rotationInDegrees);
            params.setScale(scale);
            params.setTranslationX(transX);
            params.setTranslationY(transY);

            TransformationParameters params1 = new TransformationParameters();
            params1.setRotationInDegrees(rotationInDegrees + rotRangeInDegrees);
            params1.setScale(scale);
            params1.setTranslationX(transX);
            params1.setTranslationY(transY);


            double xPt = 0;
            double yPt = 0;
            double[] transformedXY0 = tr.applyTransformation(params0,
                image1CentroidX, image1CentroidY, xPt, yPt);
            double[] transformedXY = tr.applyTransformation(params,
                image1CentroidX, image1CentroidY, xPt, yPt);
            double[] transformedXY1 = tr.applyTransformation(params1,
                image1CentroidX, image1CentroidY, xPt, yPt);

            log.info(
               String.format(
                   "  EXPECTED DIFFS SAMPLE for rotRange=%f pt at (0,0) dx's={%d,%d} dy's={%d,%d}",
                   rotRangeInDegrees,
                   (int)Math.round(Math.abs(transformedXY0[0] - transformedXY[0])),
                   (int)Math.round(Math.abs(transformedXY1[0] - transformedXY[0])),
                   (int)Math.round(Math.abs(transformedXY0[1] - transformedXY[1])),
                   (int)Math.round(Math.abs(transformedXY1[1] - transformedXY[1]))
                   ));


            xPt = image1CentroidX;
            yPt = image1CentroidY;
            transformedXY0 = tr.applyTransformation(params0,
                image1CentroidX, image1CentroidY, xPt, yPt);
            transformedXY = tr.applyTransformation(params,
                image1CentroidX, image1CentroidY, xPt, yPt);
            transformedXY1 = tr.applyTransformation(params1,
                image1CentroidX, image1CentroidY, xPt, yPt);

            log.info(
               String.format(
                   "  EXPECTED DIFFS SAMPLE for rotRange=%f pt at middle of image dx's={%d,%d} dy's={%d,%d}",
                   rotRangeInDegrees,
                   (int)Math.round(Math.abs(transformedXY0[0] - transformedXY[0])),
                   (int)Math.round(Math.abs(transformedXY1[0] - transformedXY[0])),
                   (int)Math.round(Math.abs(transformedXY0[1] - transformedXY[1])),
                   (int)Math.round(Math.abs(transformedXY1[1] - transformedXY[1]))
                   ));
        }
    }
}
