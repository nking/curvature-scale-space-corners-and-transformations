package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.features.FeatureComparisonStat;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.matrix.MatrixUtil;
import algorithms.util.PairFloatArray;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import no.uib.cipr.matrix.DenseMatrix;

/**
 *
 * @author nichole
 */
public class Distances {

    private Logger log = Logger.getLogger(this.getClass().getName());

    private static double eps = 1e-12;

    /*
    Consider the problem of computing the distance from a point (p,q) 
    to a general planar polynomial f(x,y)=0. The first order approximation of the 
    function f(x,y) about the point (p,q).
    
    
    
     */
    private EpipolarFeatureTransformationFit combineErrors(EpipolarTransformationFit distanceErrors, EpipolarFeatureTransformationFit featureErrors) {

        /*
        in order to have the distance errors count as much as the SSD errors,
        need to scale them up or SSD down by a factor.

        will use the descriptor size and the average of the maximum SSD that
        could be calculated and the minimum and make a factor for the
        distances of SSDFactor/tolerance.
         */
        // sum square diffs / size = (d0*d0) + (d3*d3).../n
        double maxSSD = 255. * 255.;

        double distScaleFactor = (maxSSD / 2.) / distanceErrors.getTolerance();

        Map<Integer, Double> indexSSDErrorsMap = new HashMap<Integer, Double>();
        Map<Integer, FeatureComparisonStat> indexFeatureMap
                = new HashMap<Integer, FeatureComparisonStat>();
        for (int i = 0; i < featureErrors.getInlierIndexes().size(); ++i) {
            Integer index = featureErrors.getInlierIndexes().get(i);
            indexSSDErrorsMap.put(index, featureErrors.getErrors().get(i));
            indexFeatureMap.put(index,
                    featureErrors.getFeatureComparisonStats().get(i));
        }

        List<Integer> outputInliers = new ArrayList<Integer>();
        List<Double> outputDistances = new ArrayList<Double>();
        List<FeatureComparisonStat> fcs = new ArrayList<FeatureComparisonStat>();

        for (int i = 0; i < distanceErrors.getInlierIndexes().size(); ++i) {
            Integer index = distanceErrors.getInlierIndexes().get(i);
            Double ssd = indexSSDErrorsMap.get(index);
            if (ssd != null) {

                outputInliers.add(index);

                Double dist = distanceErrors.getErrors().get(i);
                double cost = dist.doubleValue() * ssd.doubleValue();
                outputDistances.add(Double.valueOf(cost));

                fcs.add(indexFeatureMap.get(index));
            }
        }
        double costTerm2 = 1. / (double) outputDistances.size();
        for (int i = 0; i < outputDistances.size(); ++i) {
            double err = outputDistances.get(i).doubleValue() * costTerm2 * costTerm2
                    * distScaleFactor;
            outputDistances.set(i, Double.valueOf(err));
        }

        EpipolarFeatureTransformationFit fit
                = new EpipolarFeatureTransformationFit(
                        distanceErrors.getFundamentalMatrix(),
                        outputInliers, fcs,
                        distanceErrors.getErrorType(), outputDistances,
                        distanceErrors.getTolerance());

        return fit;
    }

    /**
     * evaluate fit for already matched point lists
     *
     * @param fm
     * @param leftPoints points from left image in matrix of size 3 X nPoints.
     * row 0 is x, row 1 is y, row 2 is all 1's
     * @param rightPoints points from right image in matrix of size 3 X nPoints.
     * row 0 is x, row 1 is y, row 2 is all 1's
     * @param tolerance
     * @return
     */
    public EpipolarTransformationFit calculateEpipolarDistanceError(
            DenseMatrix fm, DenseMatrix leftPoints, DenseMatrix rightPoints,
            double tolerance) {

        if (fm == null) {
            throw new IllegalArgumentException("fm cannot be null");
        }
        if (leftPoints == null) {
            throw new IllegalArgumentException("leftPoints cannot be null");
        }
        if (rightPoints == null) {
            throw new IllegalArgumentException("rightPoints cannot be null");
        }
        int nRows = leftPoints.numRows();
        if (nRows != rightPoints.numRows()) {
            throw new IllegalArgumentException("matrices must have same number of rows");
        }

        //2D point (x,y) and line (a, b, c): dist=(a*x + b*y + c)/sqrt(a^2 + b^2)
        PairFloatArray distances = calculateDistancesFromEpipolar(fm,
                leftPoints, rightPoints);

        List<Double> errors = new ArrayList<Double>();

        List<Integer> inlierIndexes = new ArrayList<Integer>();

        for (int i = 0; i < distances.getN(); ++i) {

            float leftPtD = distances.getX(i);

            float rightPtD = distances.getY(i);

            float dist = (float) Math.sqrt(leftPtD * leftPtD + rightPtD * rightPtD);

            if (dist > tolerance) {
                continue;
            }

            inlierIndexes.add(Integer.valueOf(i));

            errors.add(Double.valueOf(dist));
        }

        EpipolarTransformationFit fit = null;

        if (errors.size() > 0) {
            fit = new EpipolarTransformationFit(fm, inlierIndexes,
                    ErrorType.DIST_TO_EPIPOLAR_LINE, errors, tolerance);
        } else {
            fit = new EpipolarTransformationFit(fm, new ArrayList<Integer>(),
                    ErrorType.DIST_TO_EPIPOLAR_LINE, new ArrayList<Double>(), tolerance);
        }

        fit.setNMaxMatchable(leftPoints.numColumns());

        fit.calculateErrorStatistics();

        return fit;
    }

    /**
     * find the distance of the given points from their respective projective
     * epipolar lines.
     *
     * @param fm
     * @param matchedLeftPoints points from left image in matrix of size 3 X
     * nPoints. row 0 is x, row 1 is y, row 2 is all 1's
     * @param matchedRightPoints points from right image in matrix of size 3 X
     * nPoints. row 0 is x, row 1 is y, row 2 is all 1's
     * @return
     */
    public PairFloatArray calculateDistancesFromEpipolar(
            DenseMatrix fm, DenseMatrix matchedLeftPoints,
            DenseMatrix matchedRightPoints) {

        if (fm == null) {
            throw new IllegalArgumentException("fm cannot be null");
        }
        if (matchedLeftPoints == null) {
            throw new IllegalArgumentException("matchedLeftPoints cannot be null");
        }
        if (matchedRightPoints == null) {
            throw new IllegalArgumentException("rightPoints cannot be null");
        }
        int nRows = matchedLeftPoints.numRows();
        if (nRows != matchedRightPoints.numRows()) {
            throw new IllegalArgumentException("matrices must have same number of rows");
        }

        /*
        u_2^T * F * u_1 = 0  where u are the x,y points in images _1 and _2
        u_1 = (x_1, y_1, 1)^T
        u_2 = (x_2, y_2, 1)^T
         */
        int n = matchedLeftPoints.numColumns();

        PairFloatArray distances = new PairFloatArray(n);

        // F * u_1
        DenseMatrix rightEpipolarLines = MatrixUtil.multiply(fm, matchedLeftPoints);

        // F^T * u_2
        DenseMatrix leftEpipolarLines = MatrixUtil.multiply(fm.transpose(),
                matchedRightPoints);

        float[] output = new float[2];

        for (int i = 0; i < matchedLeftPoints.numColumns(); i++) {

            calculatePerpDistFromLines(matchedLeftPoints, matchedRightPoints,
                    rightEpipolarLines, leftEpipolarLines,
                    i, i, output);

            distances.add(output[0], output[1]);
        }

        return distances;
    }

    /**
     *
     * @param u1 points from left image in matrix of size 3 X nPoints. row 0 is
     * x, row 1 is y, row 2 is all 1's
     * @param u2 points from right image in matrix of size 3 X nPoints. row 0 is
     * x, row 1 is y, row 2 is all 1's
     * @param fu1
     * @param invFu2
     * @param leftIdx
     * @param rightIdx
     * @param output an output variable to hold as element 0, the distance of
     * the right image point at rightIdx from from it's left epipolar line
     * projected into the right image. element 1 holds the reverse for the left
     * image point at leftIdx.
     */
    public void calculatePerpDistFromLines(DenseMatrix u1, DenseMatrix u2,
            DenseMatrix fu1, DenseMatrix invFu2,
            int leftIdx, int rightIdx, float[] output) {

        double a = fu1.get(0, leftIdx);
        double b = fu1.get(1, leftIdx);
        double c = fu1.get(2, leftIdx);

        double aplusb = Math.sqrt((a * a) + (b * b));

        //dist = (a*x + b*y + c)/sqrt(a^2 + b^2)
        double x = u2.get(0, rightIdx);
        double y = u2.get(1, rightIdx);

        double d = (a * x + b * y + c) / aplusb;

        // find the reverse distance by projection:
        double aRev = invFu2.get(0, rightIdx);
        double bRev = invFu2.get(1, rightIdx);
        double cRev = invFu2.get(2, rightIdx);

        double xL = u1.get(0, leftIdx);
        double yL = u1.get(1, leftIdx);

        double dRev = (aRev * xL + bRev * yL + cRev)
                / Math.sqrt((aRev * aRev + bRev * bRev));

        output[0] = (float) dRev;
        output[1] = (float) d;
    }

    public EpipolarTransformationFit calculateError(DenseMatrix fm,
            DenseMatrix x1, DenseMatrix x2, ErrorType errorType, double tolerance) {

        //if (errorType.equals(ErrorType.SAMPSONS)) {
        //    return calculateSampsonsError(fm, x1, x2, tolerance);
        //} else {
        return calculateEpipolarDistanceError(fm, x1, x2, tolerance);
    }
}
