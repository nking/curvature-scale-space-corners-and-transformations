package algorithms.imageProcessing;

import static algorithms.imageProcessing.EdgeMatcher.minTolerance;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

/**
 * class to match the edges extracted from two images.  The use of corners and
 * PointMatcher should be preferred over this.
 
 * @author nichole
 */
public final class EdgeMatcher extends AbstractPointMatcher {

    private final Logger log = Logger.getLogger(this.getClass().getName());

    protected static int minTolerance = 5;

    //TODO: this has to be a high number for sets with projection.
    // the solution is sensitive to this value.
    private final float generalTolerance = 8;

    public static float toleranceGridFactor = 4.f;

    protected boolean debug = true;

    /**
     * the maximum number of iterations that the refinement of translation
     * will use in the downhill simplex.  The default value is 50.
     */
    private int dsNMaxIter = 50;
    protected void setDsNMaxIter(int n) {
        dsNMaxIter = n;
    }
    protected int getDsNMaxIter() {
        return dsNMaxIter;
    }
    float nEpsFactor = 2.0f;
    protected void setNEpsFactor(float f) {
        nEpsFactor = f;
    }
    protected float getNEpsFactor() {
        return nEpsFactor;
    }
    protected float cellFactor = 1.25f;
    protected float tolFactor = 0.5f;
    protected void setCellFactor(float f) {
        cellFactor = f;
    }
    protected void setTolFactor(float f) {
        tolFactor = f;
    }
    protected float getCellFactor() {
        return cellFactor;
    }
    protected float getTolFactor() {
        return tolFactor;
    }
  
    // ======= code that needs testing and revision the most
    /**
     * refine the transformation params to make a better match of edges1 to
     * edges2 where the points within edges in both sets are not necessarily
     * 1 to 1 matches (that is, the input is not expected to be matched
     * already).
     *
     * TODO: improve transformEdges to find translation for all edges
     * via a search method rather than trying all pairs of points.
     *
     * @param edges1
     * @param edges2
     * @param params
     * @param centroidX1
     * @param centroidY1
     * @param centroidX2
     * @param centroidY2
     * @return
     */
    public TransformationParameters refineTransformation(PairIntArray[] edges1,
        PairIntArray[] edges2, final TransformationParameters params,
        final int centroidX1, final int centroidY1,
        final int centroidX2, final int centroidY2) {

        if (edges1 == null || edges1.length == 0) {
            throw new IllegalArgumentException("edges1 cannot be null or empty");
        }
        if (edges2 == null || edges2.length == 0) {
            throw new IllegalArgumentException("edges2 cannot be null or empty");
        }

        //TODO: set this empirically from tests
        double convergence = 0;

        double r = params.getRotationInRadians();
        double s = params.getScale();

        double rMin = r - (10 * Math.PI/180);
        double rMax = r + (10 * Math.PI/180);
        double sMin = s - 1.5;
        double sMax = s + 1.5;

        // the positive offsets can be found w/ reflection?
        // TODO: needs testing for starter points.  these are supplying the
        // "grid search" portion of exploring more than local space
        double[] drs = new double[] {
            -5.0 * Math.PI/180.,
            -2.5 * Math.PI/180.,
            -1.0 * Math.PI/180.,
            1.0 * Math.PI/180.,
            2.5 * Math.PI/180.,
            5.0 * Math.PI/180.
        };
        double[] dss = new double[] {
            -1.0, -0.1, -0.05, 0.05 /*, 0.05, 0.1, 1.0*/
        };
        if (r == 0) {
             drs = new double[]{0};
        }
        if (s == 1) {
            dss = new double[]{0};
            sMin = 1;
        }
        if (sMin < 1) {
            sMin = 1;
        }
        if (rMin < 0) {
            rMin = 0;
        }

        int n = (1 + dss.length) * (1 + drs.length);

        TransformationPointFit[] fits = new TransformationPointFit[n];

        int count = 0;
        for (int i = 0; i <= dss.length; i++) {

            double scale = (i == 0) ? s : s + dss[i - 1];

            for (int j = 0; j <= drs.length; j++) {

                double rotation = (j == 0) ? r : r + drs[j - 1];

                fits[count] = calculateTranslationAndTransformEdges(
                    rotation, scale, edges1, edges2,
                    centroidX1, centroidY1);

                if (fits[count] != null) {
                    count++;
                }
            }
        }

        if (count < n) {
            fits = Arrays.copyOf(fits, count);
        }

        float alpha = 1;   // > 0
        float gamma = 2;   // > 1
        float beta = 0.5f;
        float tau = 0.5f;

        boolean go = true;

        int nMaxIter = 100;
        int nIter = 0;

        int bestFitIdx = 0;
        int worstFitIdx = fits.length - 1;

        int lastNMatches = Integer.MIN_VALUE;
        double lastAvgDistModel = Double.MAX_VALUE;
        int nIterSameMin = 0;

        while (go && (nIter < nMaxIter)) {

            if (fits.length == 0) {
                break;
            }

            sortByDescendingMatches(fits, 0, (fits.length - 1));

            for (int i = (fits.length - 1); i > -1; --i) {
                if (fits[i] != null) {
                    worstFitIdx = i;
                    break;
                }
            }
            if (fits[bestFitIdx] == null) {
                break;
            }

/*if (fits.length > 0) {
    log.info("best fit:  n=" + fits[bestFitIdx].getNumberOfMatchedPoints()
    + " dm=" + fits[bestFitIdx].getMeanDistFromModel()
    + " params:\n" + fits[bestFitIdx].getParameters().toString());
}*/

            if ((lastNMatches == fits[bestFitIdx].getNumberOfMatchedPoints()) &&
                (Math.abs(lastAvgDistModel -
                fits[bestFitIdx].getMeanDistFromModel()) < 0.01)) {

                nIterSameMin++;
                /*if (nIterSameMin >= 5) {
                    break;
                }*/

            } else {
                nIterSameMin = 0;
            }

            lastNMatches = fits[bestFitIdx].getNumberOfMatchedPoints();
            lastAvgDistModel = fits[bestFitIdx].getMeanDistFromModel();

            // determine center for all points excepting the worse fit
            double rSum = 0.0;
            double sSum = 0.0;
            int c = 0;
            for (int i = 0; i < (fits.length - 1); i++) {
                if (fits[i] != null) {
                    rSum += fits[i].getRotationInRadians();
                    sSum += fits[i].getScale();
                    c++;
                }
            }
            r = rSum / (double)c;
            s = sSum / (double)c;

            // "Reflection"
            double rReflect = r + (alpha *
                (r - fits[worstFitIdx].getRotationInRadians()));
            double sReflect = s + (alpha *
                (s - fits[worstFitIdx].getScale()));

            TransformationPointFit fitReflected =
                calculateTranslationAndTransformEdges(
                rReflect, sReflect, edges1, edges2,
                centroidX1, centroidY1);

            //TODO: consider putting back in a check for bounds

            int comp0 = compare(fits[bestFitIdx], fitReflected);
            int compLast = compare(fits[worstFitIdx], fitReflected);

            if ((comp0 < 1) && (compLast == 1)) {

                // replace last with f_refl
                fits[worstFitIdx] = fitReflected;

            } else if (comp0 == 1) {

                // reflected is better than best fit, so "expand"
                // "Expansion"
                double rExpansion = r + (gamma * (rReflect - r));
                double sExpansion = s + (gamma * (sReflect - s));

                TransformationPointFit fitExpansion =
                    calculateTranslationAndTransformEdges(
                    rExpansion, sExpansion, edges1, edges2,
                    centroidX1, centroidY1);

                int compR = compare(fitReflected, fitExpansion);

                if (compR == 1) {

                    // expansion fit is better than reflected fit
                    fits[worstFitIdx] = fitExpansion;

                } else {

                    fits[worstFitIdx] = fitReflected;
                }

            } else if (compLast < 1) {

                // reflected fit is worse than the worst (last) fit, so contract
                // "Contraction"
                double rContraction = r + (beta *
                    (fits[worstFitIdx].getRotationInRadians() - r));
                double sContraction = s + (beta *
                    (fits[worstFitIdx].getScale() - s));

                TransformationPointFit fitContraction =
                    calculateTranslationAndTransformEdges(
                    rContraction, sContraction, edges1, edges2,
                    centroidX1, centroidY1);

                int compC = compare(fits[worstFitIdx], fitContraction);

                if (compC > -1) {

                    fits[worstFitIdx] = fitContraction;

                } else {

                    // "Reduction"
                    for (int i = 1; i < fits.length; ++i) {

                        if (fits[i] == null) {
                            /*TODO: consider setting this
                            fits[i] = new TransformationPointFit(
                                new TransformationParameters(),
                                0, Double.MAX_VALUE, Double.MAX_VALUE,
                                Double.MAX_VALUE
                            );
                            */
                            continue;
                        }

                        float rReduction
                            = (fits[bestFitIdx].getRotationInRadians()
                            + (tau
                            * (fits[i].getRotationInRadians()
                            - fits[bestFitIdx].getRotationInRadians())));

                        float sReduction
                            = (fits[bestFitIdx].getScale()
                            + (tau
                            * (fits[i].getScale()
                            - fits[bestFitIdx].getScale())));

                        //NOTE: there's a possibility of a null fit.
                        //  instead of re-writing the fits array, will
                        //  assign a fake infinitely bad fit which will
                        //  fall to the bottom of the list after the next
                        //  sort.
                        TransformationPointFit fit =
                            calculateTranslationAndTransformEdges(
                            rReduction, sReduction, edges1, edges2,
                            centroidX1, centroidY1);

                        fits[i] = fit;
                    }
                }
            }

            log.finest("best fit so far: nMatches="
                + fits[bestFitIdx].getNumberOfMatchedPoints()
                + " diff from model=" + fits[bestFitIdx].getMeanDistFromModel()
            );

            nIter++;

            if ((fits[bestFitIdx].getNumberOfMatchedPoints() == convergence)
                && (fits[bestFitIdx].getMeanDistFromModel() == 0)) {
                go = false;
            /*} else if ((r > rMax) || (r < rMin)) {
                go = false;*/
            } else if ((s > sMax) || (s < sMin)) {
                go = false;
            }
        }

        // additional step that's helpful if not enough iterations are used,
        // is to test the summed transX, transY which represent the center
        // of the simplex against the best fit
        TransformationPointFit fitAvg = calculateTranslationAndTransformEdges(
            r, s, edges1, edges2, centroidX1, centroidY1);

        int comp = compare(fits[bestFitIdx], fitAvg);
        if (comp == 1) {
            fits[bestFitIdx] = fitAvg;
        }

        // if rotation > 2PI, subtract 2PI
        if ((fits[bestFitIdx] != null) &&
            (fits[bestFitIdx].getParameters().getRotationInRadians()
            > 2.*Math.PI)) {

            float rot = fits[bestFitIdx].getParameters().getRotationInRadians();
            while (rot >= 2*Math.PI) {
                rot -= 2*Math.PI;
            }
            fits[bestFitIdx].getParameters().setRotationInRadians(rot);
        }

        return fits[bestFitIdx].getParameters();
    }

    /**
     * TODO: improve transformEdges to find translation for all edges
     * via a search method rather than trying all pairs of points.
     *
     * Given edges1 and edges2 which we already know are matched edges due to
     * contour matching or other means, and given the rotation and scale,
     * determine the translation between the edges and return the fit.
     *
     * @param rotInRad
     * @param scl
     * @param edges1
     * @param edges2
     * @param centroidX1
     * @param centroidY1
     * @return
     */
    private TransformationPointFit calculateTranslationAndTransformEdges(
        double rotInRad, double scl,
        PairIntArray[] edges1, PairIntArray[] edges2,
        int centroidX1, int centroidY1) {

        if (edges1 == null || edges1.length == 0) {
            throw new IllegalArgumentException("edges1 cannot be null or empty");
        }
        if (edges2 == null || edges2.length == 0) {
            throw new IllegalArgumentException("edges2 cannot be null or empty");
        }
        if ((edges1.length != edges2.length)) {
            throw new IllegalArgumentException(
            "edges1 and edges2 must be the same length");
        }

        /*
        edges1 and edges2 are matched edges, but should be considered clouds
        of points rather than point to point matches within the edges.

        For each paired edge in edges1 and edges2, determine the implied
        translation by their centroids.
        These are combined in weighted average where the
        weight is the length of the edge.

        Then a small area surrounding the averaged translation is searched
        to find the best fit.
        */

        float s = (float)scl;
        float scaleTimesCosine = (float)(s * Math.cos(rotInRad));
        float scaleTimesSine = (float)(s * Math.sin(rotInRad));

        //TODO: revisit this:
        float tolTransX = 2.f * centroidX1 * 0.02f;
        float tolTransY = 2.f * centroidY1 * 0.02f;
        if (tolTransX < minTolerance) {
            tolTransX = minTolerance;
        }
        if (tolTransY < minTolerance) {
            tolTransY = minTolerance;
        }

        int nTotal = 0;
        float[] weights = new float[edges1.length];
        for (int i = 0; i < edges1.length; i++) {
            weights[i] = edges1[i].getN();
            nTotal += edges1[i].getN();
        }
        for (int i = 0; i < edges1.length; i++) {
            weights[i] /= (float)nTotal;
        }

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        float translationX = 0;
        float translationY = 0;
        for (int i = 0; i < edges1.length; i++) {

            PairIntArray edge1 = edges1[i];
            PairIntArray edge2 = edges2[i];

            double[] xycen1 = curveHelper.calculateXYCentroids(edge1);

            double[] xycen2 = curveHelper.calculateXYCentroids(edge2);

            double srX1 = centroidX1*s + (
                ((xycen1[0] - centroidX1) * scaleTimesCosine) +
                ((xycen1[1] - centroidY1) * scaleTimesSine));

            double srY1 = centroidY1*s + (
                (-(xycen1[0] - centroidX1) * scaleTimesSine) +
                ((xycen1[1] - centroidY1) * scaleTimesCosine));

            double tx = xycen2[0] - srX1;

            double ty = xycen2[1] - srY1;

            translationX += weights[i] * tx;

            translationY += weights[i] * ty;
        }

        TransformationPointFit bestFit = refineTranslationWithDownhillSimplex(
            edges1, edges2, translationX, translationY, tolTransX, tolTransY,
            20.f, 20.f, s, (float)rotInRad, centroidX1, centroidY1);

        return bestFit;
    }

    /**
     * sort the fits by descending number of matches.
     * @param fits
     * @param idxLo
     * @param idxHi, upper index, inclusive
     */
    void sortByDescendingMatches(TransformationPointFit[] fits, int idxLo,
        int idxHi) {

        if (idxLo < idxHi) {

            int idxMid = partition(fits, idxLo, idxHi);

            sortByDescendingMatches(fits, idxLo, idxMid - 1);

            sortByDescendingMatches(fits, idxMid + 1, idxHi);
        }
    }

    private int partition(TransformationPointFit[] fits, int idxLo, int idxHi) {

        TransformationPointFit x = fits[idxHi];

        int store = idxLo - 1;

        for (int i = idxLo; i < idxHi; ++i) {
            if (fitIsBetter(x, fits[i])) {
                store++;
                TransformationPointFit swap = fits[store];
                fits[store] = fits[i];
                fits[i] = swap;
            }
        }

        store++;
        TransformationPointFit swap = fits[store];
        fits[store] = fits[idxHi];
        fits[idxHi] = swap;

        return store;
    }

    private TransformationPointFit refineTranslationWithDownhillSimplex(
        final PairIntArray[] edges1, final PairIntArray[] edges2,
        float transX, float transY, float tolTransX, float tolTransY,
        float plusMinusTransX, float plusMinusTransY,
        final float scale, final float rotationRadians,
        int centroidX1, int centroidY1) {

        int n1 = 0;
        for (PairIntArray edge : edges1) {
            n1 += edge.getN();
        }
        int n2 = 0;
        for (PairIntArray edge : edges2) {
            n2 += edge.getN();
        }
        int nMaxMatchable = (n1 < n2) ? n1 : n2;

        if (nMaxMatchable == 0) {
            return null;
        }

        //TODO: revise this:
        double eps = Math.log(nMaxMatchable)/Math.log(10);

        float txMin = transX - plusMinusTransX;
        float txMax = transX + plusMinusTransX;
        float tyMin = transY - plusMinusTransY;
        float tyMax = transY + plusMinusTransY;

        // the positive offsets can be found w/ reflection
        float[] dtx = new float[]{
            -plusMinusTransX, -0.5f*plusMinusTransX,
            -0.25f*plusMinusTransX,
            -0.125f*plusMinusTransX, -0.0625f*plusMinusTransX
        };
        float[] dty = new float[]{
            -plusMinusTransY, -0.5f*plusMinusTransY,
            -0.25f*plusMinusTransY,
            -0.125f*plusMinusTransY, -0.0625f*plusMinusTransY
        };
        int n = (1 + dtx.length) * (1 + dty.length);

        TransformationPointFit[] fits = new TransformationPointFit[n];

        int count = 0;
        for (int i = 0; i <= dtx.length; i++) {

            float tx = (i == 0) ? transX : (transX + dtx[i - 1]);

            for (int j = 0; j <= dty.length; j++) {

                float ty = (i == 0) ? transY : (transY + dty[i - 1]);

                fits[count] = evaluateFit(edges1, edges2,
                    tx, ty, tolTransX, tolTransY,
                    scale, rotationRadians, centroidX1, centroidY1);

                if (fits[count] != null) {
                    count++;
                }
            }
        }

        if (count < n) {
            fits = Arrays.copyOf(fits, count);
        }

        float alpha = 1;   // > 0
        float gamma = 2;   // > 1
        float beta = 0.5f;
        float tau = 0.5f;

        boolean go = true;

        int nMaxIter = 100;
        int nIter = 0;

        int bestFitIdx = 0;
        int worstFitIdx = fits.length - 1;

        int lastNMatches = Integer.MIN_VALUE;
        double lastAvgDistModel = Double.MAX_VALUE;
        int nIterSameMin = 0;

        while (go && (nIter < nMaxIter)) {

            if (fits.length == 0) {
                break;
            }

            sortByDescendingMatches(fits, 0, (fits.length - 1));

            for (int i = (fits.length - 1); i > -1; --i) {
                if (fits[i] != null) {
                    worstFitIdx = i;
                    break;
                }
            }
            if (fits[bestFitIdx] == null) {
                break;
            }

            if ((lastNMatches == fits[bestFitIdx].getNumberOfMatchedPoints())
                && (Math.abs(lastAvgDistModel
                    - fits[bestFitIdx].getMeanDistFromModel()) < 0.01)) {

                nIterSameMin++;
                /*if (nIterSameMin >= 10) {
                    break;
                }*/

            } else {
                nIterSameMin = 0;
            }

            lastNMatches = fits[bestFitIdx].getNumberOfMatchedPoints();
            lastAvgDistModel = fits[bestFitIdx].getMeanDistFromModel();

            // determine center for all points excepting the worse fit
            float txSum = 0.0f;
            float tySum = 0.0f;
            int c = 0;
            for (int i = 0; i < (fits.length - 1); ++i) {
                if (fits[i] != null) {
                    txSum += fits[i].getTranslationX();
                    tySum += fits[i].getTranslationY();
                    ++c;
                }
            }
            transX = txSum / (float)c;
            transY = tySum / (float)c;

            // "Reflection"
            float txReflect = transX + (alpha
                * (transX - fits[worstFitIdx].getTranslationX()));
            float tyReflect = transY + (alpha
                * (transY - fits[worstFitIdx].getTranslationY()));

            TransformationPointFit fitReflected = evaluateFit(
                edges1, edges2, txReflect, tyReflect,
                tolTransX, tolTransY,
                scale, rotationRadians, centroidX1, centroidY1);

            //TODO: consider putting back in a check for bounds

            int comp0 = compare(fits[bestFitIdx], fitReflected);
            int compLast = compare(fits[worstFitIdx], fitReflected);

            if ((comp0 < 1) && (compLast == 1)) {

                // replace last with f_refl
                fits[worstFitIdx] = fitReflected;

            } else if (comp0 == 1) {

                // reflected is better than best fit, so "expand"
                // "Expansion"
                float txExpansion = transX + (gamma * (txReflect - transX));
                float tyExpansion = transY + (gamma * (tyReflect - transY));

                TransformationPointFit fitExpansion =
                    evaluateFit(edges1, edges2,
                        txExpansion, tyExpansion, tolTransX, tolTransY,
                        scale, rotationRadians, centroidX1, centroidY1);

                int compR = compare(fitReflected, fitExpansion);

                if (compR == 1) {

                    // expansion fit is better than reflected fit
                    fits[worstFitIdx] = fitExpansion;

                } else {

                    fits[worstFitIdx] = fitReflected;
                }

            } else if (compLast < 1) {

                // reflected fit is worse than the worst (last) fit, so contract
                // "Contraction"
                float txContraction = transX + (beta
                    * (fits[worstFitIdx].getTranslationX() - transX));
                float tyContraction = transY + (beta
                    * (fits[worstFitIdx].getTranslationY() - transY));

                TransformationPointFit fitContraction
                    = evaluateFit(edges1, edges2,
                        txContraction, tyContraction,
                        tolTransX, tolTransY,
                        scale, rotationRadians, centroidX1, centroidY1);

                int compC = compare(fits[worstFitIdx], fitContraction);

                if (compC > -1) {

                    fits[worstFitIdx] = fitContraction;

                } else {

                    // "Reduction"
                    for (int i = 1; i < fits.length; i++) {

                        if (fits[i] == null) {
                            /* TODO: consider setting this
                            fits[i] = new TransformationPointFit(
                                new TransformationParameters(),
                                0, Double.MAX_VALUE, Double.MAX_VALUE,
                                Double.MAX_VALUE
                            );
                            */
                            continue;
                        }

                        float txReduction
                            = (fits[bestFitIdx].getTranslationX()
                            + (tau
                            * (fits[i].getTranslationX()
                            - fits[bestFitIdx].getTranslationX())));

                        float tyReduction
                            = (fits[bestFitIdx].getTranslationY()
                            + (tau
                            * (fits[i].getTranslationY()
                            - fits[bestFitIdx].getTranslationY())));

                        //NOTE: there's a possibility of a null fit.
                        //  instead of re-writing the fits array, will
                        //  assign a fake infinitely bad fit which will
                        //  fall to the bottom of the list after the next
                        //  sort.
                        TransformationPointFit fit
                            = evaluateFit(edges1, edges2,
                                txReduction, tyReduction,
                                tolTransX, tolTransY,
                                scale, rotationRadians, centroidX1, centroidY1);

                        fits[i] = fit;
                    }
                }
            }

            log.finest("best fit so far: nMatches="
                + fits[bestFitIdx].getNumberOfMatchedPoints()
                + " diff from model=" + fits[bestFitIdx].getMeanDistFromModel()
            );

            nIter++;

            if ((fits[bestFitIdx].getNumberOfMatchedPoints() == nMaxMatchable)
                && (fits[bestFitIdx].getMeanDistFromModel() <  eps)) {
                go = false;
            } /*else if ((transX > txMax) || (transX < txMin)) {
                go = false;
            } else if ((transY > tyMax) || (transY < tyMin)) {
                go = false;
            }*/
        }

        // additional step that's helpful if not enough iterations are used,
        // is to test the summed transX, transY which represent the center
        // of the simplex against the best fit
        TransformationPointFit fitAvg = evaluateFit(edges1, edges2,
            transX, transY, tolTransX, tolTransY,
            scale, rotationRadians, centroidX1, centroidY1);

        int comp = compare(fits[bestFitIdx], fitAvg);
        if (comp == 1) {
            fits[bestFitIdx] = fitAvg;
        }

        return fits[bestFitIdx];
    }

    private TransformationPointFit evaluateFit(
        PairIntArray[] edges1, PairIntArray[] edges2,
        float translationX, float translationY,
        float tolTransX, float tolTransY, float scale,
        float rotationRadians, int centroidX1, int centroidY1) {

        if (edges1 == null || edges1.length == 0) {
            throw new IllegalArgumentException("edges1 cannot be null or empty");
        }
        if (edges2 == null || edges2.length == 0) {
            throw new IllegalArgumentException("edges2 cannot be null or empty");
        }

        int n1 = 0;
        for (PairIntArray edge : edges1) {
            n1 += edge.getN();
        }
        int n2 = 0;
        for (PairIntArray edge : edges2) {
            n2 += edge.getN();
        }
        int nMaxMatchable = (n1 < n2) ? n1 : n2;

        if (nMaxMatchable == 0) {
            return null;
        }

        List<Double> residuals = new ArrayList<Double>();

        int nTotal = 0;

        for (int i = 0; i < edges1.length; i++) {

            PairIntArray edge1 = edges1[i];

            PairIntArray edge2 = edges2[i];

            nTotal += edge1.getN();

            calculateTranslationResidualsForUnmatchedMultiplicity(
                edge1, edge2, translationX, translationY, tolTransX, tolTransY,
                scale, rotationRadians,
                centroidX1, centroidY1, residuals);
        }

        double avg = 0;
        for (Double diff : residuals) {
            avg += diff.doubleValue();
        }
        avg /= (double)residuals.size();

        double stdDev = 0;
        for (Double diff : residuals) {
            double d = diff.doubleValue() - avg;
            stdDev += (d * d);
        }
        stdDev = Math.sqrt(stdDev/((double)residuals.size() - 1));

        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians(rotationRadians);
        params.setScale(scale);
        params.setTranslationX(translationX);
        params.setTranslationY(translationY);

        TransformationPointFit fit = new TransformationPointFit(params,
            residuals.size(), avg, stdDev,
            tolTransX, tolTransY);

        return fit;
    }

    
    /**
     * calculate the residuals between edge1 points and edge2 points with the
     * assumption that there may not be one to one point matches, but instead
     * will be general point matches.  For that reason, a greedy search for
     * nearest neighbor with the possibility of multiple matches to a neighbor
     * is used.
     *
     * @param edge1
     * @param edge2
     * @param transX
     * @param transY
     * @param tolTransX
     * @param tolTransY
     * @param scale
     * @param rotationRadians
     * @param outputResiduals
     */
    private void calculateTranslationResidualsForUnmatchedMultiplicity(
        PairIntArray edge1, PairIntArray edge2,
        float transX, float transY, float tolTransX, float tolTransY,
        final float scale, final float rotationRadians,
        int centroidX1, int centroidY1, List<Double> outputResiduals) {

        if (edge1 == null || edge1.getN() == 0) {
            throw new IllegalArgumentException(
            "edge1 cannot be null or empty");
        }
        if (edge2 == null || edge2.getN() == 0) {
            throw new IllegalArgumentException(
            "edge2 cannot be null or empty");
        }

        int nMaxMatchable = (edge1.getN() < edge2.getN()) ?
            edge1.getN() : edge2.getN();

        if (nMaxMatchable == 0) {
            return;
        }

        float scaleTimesCosine = (float)(scale * Math.cos(rotationRadians));
        float scaleTimesSine = (float)(scale * Math.sin(rotationRadians));

        for (int i = 0; i < edge1.getN(); i++) {

            int x1 = edge1.getX(i);
            int y1 = edge1.getY(i);

            float transformedX = (centroidX1*scale + (
                ((x1 - centroidX1) * scaleTimesCosine) +
                ((y1 - centroidY1) * scaleTimesSine))) + transX;

            float transformedY = (centroidY1*scale + (
                (-(x1 - centroidX1) * scaleTimesSine) +
                ((y1 - centroidY1) * scaleTimesCosine))) + transY;

            double minDiff = Double.MAX_VALUE;

            for (int j = 0; j < edge2.getN(); j++) {

                int x2 = edge2.getX(j);
                int y2 = edge2.getY(j);

                float dx = x2 - transformedX;
                float dy = y2 - transformedY;

                if ((Math.abs(dx) > tolTransX) || (Math.abs(dy) > tolTransY)) {
                    continue;
                }

                float diff = (float)Math.sqrt(dx*dx + dy*dy);

                if (diff < minDiff) {
                    minDiff = diff;
                }
            }

            if (minDiff < Double.MAX_VALUE) {
                outputResiduals.add(Double.valueOf(minDiff));
            }
        }
    }

}
