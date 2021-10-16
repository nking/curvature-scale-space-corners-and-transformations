package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.features.FeatureComparisonStat;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath;
import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import algorithms.util.PairFloatArray;
import gnu.trove.list.TIntList;
import java.util.ArrayList;
import java.util.Arrays;
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
     * evaluate fit for the distances of each correspondence point to the
     * epipolar line for it projected from it's pair in the other image.
      Luong et al. 1993, "On determining the fundamental matrix?: analysis of 
     different methods and experimental results"
     *
     * @param fm
     * @param leftPoints points from left image in matrix of size 3 X nPoints.
     * row 0 is x, row 1 is y, row 2 is all 1's
     * @param rightPoints points from right image in matrix of size 3 X nPoints.
     * row 0 is x, row 1 is y, row 2 is all 1's
     * @param threshhold
     * @return
     */
    public EpipolarTransformationFit calculateEpipolarDistanceError(
            DenseMatrix fm, DenseMatrix leftPoints, DenseMatrix rightPoints,
            double threshhold) {

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
        double[][] distances = calculateDistancesFromEpipolar(fm,
            leftPoints, rightPoints);

        double[] combinedDist = combineDistances(distances);
                        
        List<Double> errors = new ArrayList<Double>();

        List<Integer> inlierIndexes = new ArrayList<Integer>();

        for (int i = 0; i < combinedDist.length; ++i) {

            if (combinedDist[i] > threshhold) {
                continue;
            }

            inlierIndexes.add(Integer.valueOf(i));

            errors.add(Double.valueOf(combinedDist[i]));
        }

        EpipolarTransformationFit fit = null;

        if (errors.size() > 0) {
            fit = new EpipolarTransformationFit(fm, inlierIndexes,
                    ErrorType.DIST_TO_EPIPOLAR_LINE, errors, threshhold);
        } else {
            fit = new EpipolarTransformationFit(fm, new ArrayList<Integer>(),
                    ErrorType.DIST_TO_EPIPOLAR_LINE, new ArrayList<Double>(), threshhold);
        }

        fit.setNMaxMatchable(leftPoints.numColumns());

        fit.calculateErrorStatistics();

        return fit;
    }
    
    /**
     * For use upon data that have been unit normal standardized.
     * evaluate fit for already matched point lists
     *
     * @param fm
     * @param leftPoints points from left image in matrix of size 3 X nPoints.
     * row 0 is x, row 1 is y, row 2 is all 1's
     * @param rightPoints points from right image in matrix of size 3 X nPoints.
     * row 0 is x, row 1 is y, row 2 is all 1's
     * @param threshhold
     * @return
     */
    public EpipolarTransformationFit calculateSampsonsError(
            DenseMatrix fm, DenseMatrix leftPoints, DenseMatrix rightPoints,
            double threshhold) {

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

        double[] distSquared = calculateEpipolarSampsonsDistanceSquared(fm, 
            leftPoints, rightPoints);

        List<Double> errors = new ArrayList<Double>();

        List<Integer> inlierIndexes = new ArrayList<Integer>();

        double d;
        for (int i = 0; i < distSquared.length; ++i) {

            d = Math.sqrt(distSquared[i]);
            
            if (d > threshhold) {
                continue;
            }

            inlierIndexes.add(Integer.valueOf(i));

            errors.add(d);
        }

        EpipolarTransformationFit fit = null;

        if (errors.size() > 0) {
            fit = new EpipolarTransformationFit(fm, inlierIndexes,
                    ErrorType.SAMPSONS, errors, threshhold);
        } else {
            fit = new EpipolarTransformationFit(fm, new ArrayList<Integer>(),
                    ErrorType.SAMPSONS, new ArrayList<Double>(), threshhold);
        }

        fit.setNMaxMatchable(leftPoints.numColumns());

        fit.calculateErrorStatistics();

        return fit;
    }
    
    /**
     * given the fundamental matrix solution and the correspondence pairs,
     * use the Luong et al. 1993 distance for each correspondence pair as
     * the closest distance of a point from the projected epipolar line of
     * its pair.   The standard deviation of the errors,
     * and the chi-squared statistic factor are used to remove
     * outliers and return the results.
      Luong et al. 1993, "On determining the fundamental matrix?: analysis of 
      different methods and experimental results".
      NOTE: there may be other references for this method.
     *
     * @param fm
     * @param leftPoints points from left image in matrix of size 3 X nPoints.
     * row 0 is x, row 1 is y, row 2 is all 1's
     * @param rightPoints points from right image in matrix of size 3 X nPoints.
     * row 0 is x, row 1 is y, row 2 is all 1's
     * @param chiSquaredStatFactor
     * @return
     */
    public EpipolarTransformationFit calculateEpipolarDistanceError2(
            DenseMatrix fm, DenseMatrix leftPoints, DenseMatrix rightPoints,
            double chiSquaredStatFactor) {

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
        double[][] distances = calculateDistancesFromEpipolar(fm,
            leftPoints, rightPoints);
        
        double[] combinedDist = combineDistances(distances);
        
        double[] meanAndStDev = MiscMath.getAvgAndStDev(combinedDist);
        
        double t = chiSquaredStatFactor * meanAndStDev[1];
        
        List<Double> errors = new ArrayList<Double>();

        List<Integer> inlierIndexes = new ArrayList<Integer>();

        for (int i = 0; i < combinedDist.length; ++i) {

            if (combinedDist[i] > t) {
                continue;
            }

            inlierIndexes.add(Integer.valueOf(i));

            errors.add(combinedDist[i]);
        }

        EpipolarTransformationFit fit = null;

        if (errors.size() > 0) {
            fit = new EpipolarTransformationFit(fm, inlierIndexes,
                ErrorType.DIST_TO_EPIPOLAR_LINE, errors, t);
        } else {
            fit = new EpipolarTransformationFit(fm, new ArrayList<Integer>(),
                ErrorType.DIST_TO_EPIPOLAR_LINE, new ArrayList<Double>(), t);
        }

        fit.setNMaxMatchable(leftPoints.numColumns());

        fit.calculateErrorStatistics();

        return fit;
    }
    
    /**
     * given rectified correspondence pairs in x1, x2, calculate the
     * square root of the sum of squares of the differences in the
     * y coordinates, or both x and y if useXToo is true.
     * Also uses MAD or Tukey fences to determine inliers.
     * 
     * @param x1 the image 1 set of correspondence points. format is 3 x N
     * where N is the number of points. NOTE: since intrinsic parameters are not
     * known, users of this method should presumably center the coordinates in
     * some manner (e.g. subtract the image center or centroid of points) since
     * internally an identity matrix is used for K.
     * @param x2 the image 2 set of correspondence points. format is 3 x N where
     * N is the number of points. NOTE: since intrinsic parameters are not
     * known, users of this method should presumably center the coordinates in
     * some manner (e.g. subtract the image center or centroid of points).
     * @param outErrors array of length x1[0].length to be used to return each
     * error as the square root of the sum of differences between each x1 and x2.
     * @param outInlierIndexes indexes of points which are within the robust 
     * statistics range for inliers (defined as not being outliers).
     * @param useXToo if true, the errors are not only the differences between the
     * y coordinates of x1 and x2, but also the x-xoordinates of x1 and x2.
     * @return square root of the sum of squared errors
     */
    public static double calculateRectificationErrors(
    double[][] x1, double[][] x2,
    double[] outErrors, TIntList outInlierIndexes, boolean useXToo) {
        
        int n = x1[0].length;
        if (x1.length !=3 && x1.length != 2) {
            throw new IllegalArgumentException("x1.length must be 3 or 2");
        }
        if (x2.length !=3 && x2.length != 2) {
            throw new IllegalArgumentException("x2.length  must be 3 or 2");
        }
        if (x2[0].length != n) {
            throw new IllegalArgumentException("x2 must be size 3Xn and the same as size of x1");
        }
        if (outErrors.length !=n) {
            throw new IllegalArgumentException("outErrors.length  must be n");
        }
        if (outInlierIndexes == null) {
            throw new IllegalArgumentException("outInlierIndexes cannot be null");
        }
        double diff;
        double sum = 0;
        double sumI;
        int i;
        for (i = 0; i < n; ++i) {
            diff = x1[1][i] - x2[1][i];
            sumI = diff*diff;
            if (useXToo) {
                diff = x1[0][i] - x2[0][i];
                sumI += diff*diff;
            }
            sum += sumI;
            outErrors[i] = Math.sqrt(sumI);
        }
        sum = Math.sqrt(sum);
        
        // use MAD or Tukey fences
                
        //median of absolute deviation of x, median, min, and max.
        double[] mADMinMax = MiscMath0.calculateMedianOfAbsoluteDeviation(outErrors);
        double kMAD = 1.4826;
        double s = kMAD*mADMinMax[0];
        double r0 = mADMinMax[1] - 3*s;
        double r1 = mADMinMax[1] + 3*s;
        
        int[] inliers2 = MiscMath0.findInliersUsingTukeyFences(outErrors);
        
        for (i = 0; i < n; ++i) {
            diff = outErrors[i];
            if (diff >= r0 && diff <= r1) {
                outInlierIndexes.add(i);
            }
        }
        
        System.out.printf("inliers using MAD=\n%s\ninliers using Tukey fences=\n%s\n", 
            Arrays.toString(outInlierIndexes.toArray()),
            Arrays.toString(inliers2));
        
        return sum;
    }
    
    /**
     * For use upon data that have been unit normal standardized.
     * given the fundamental matrix solution and the correspondence pairs,
     * use the Sampson distance as an error for each correspondence pair, then use
     * the standard deviation of the errors,
     * and the chi-squared statistic factor to remove
     * outliers and return the results.
     *
     * @param fm
     * @param leftPoints points from left image in matrix of size 3 X nPoints.
     * row 0 is x, row 1 is y, row 2 is all 1's
     * @param rightPoints points from right image in matrix of size 3 X nPoints.
     * row 0 is x, row 1 is y, row 2 is all 1's
     * @param chiSquaredStatFactor
     * @return
     */
    public EpipolarTransformationFit calculateSampsonsError2(
        DenseMatrix fm, DenseMatrix leftPoints, DenseMatrix rightPoints,
        double chiSquaredStatFactor) {

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

        // dist^2:
        double[] dist = calculateEpipolarSampsonsDistanceSquared(fm, 
            leftPoints, rightPoints);
        for (int i = 0; i < dist.length; ++i) {
            dist[i] = Math.sqrt(dist[i]);
        }
        
        double[] meanAndStDev = MiscMath.getAvgAndStDev(dist);
        
        double t = chiSquaredStatFactor * meanAndStDev[1];
        
        //System.out.printf("t=%.3e\n", t);

        List<Double> errors = new ArrayList<Double>();

        List<Integer> inlierIndexes = new ArrayList<Integer>();

        double d;
        for (int i = 0; i < dist.length; ++i) {

            d = dist[i];
            
            if (d > t) {
                continue;
            }

            inlierIndexes.add(Integer.valueOf(i));

            errors.add(d);
        }

        EpipolarTransformationFit fit = null;

        if (errors.size() > 0) {
            fit = new EpipolarTransformationFit(fm, inlierIndexes,
                ErrorType.SAMPSONS, errors, t);
        } else {
            fit = new EpipolarTransformationFit(fm, new ArrayList<Integer>(),
                ErrorType.SAMPSONS, new ArrayList<Double>(), t);
        }

        fit.setNMaxMatchable(leftPoints.numColumns());

        fit.calculateErrorStatistics();

        return fit;
    }
    
    /**
     * find the distance of the given points from their respective projected
     * epipolar lines.
     Luong et al. 1993, "On determining the fundamental matrix?: analysis of 
     different methods and experimental results"
     
     * @param fm
     * @param matchedLeftPoints points from left image in matrix of size 3 X
     * nPoints. row 0 is x, row 1 is y, row 2 is all 1's
     * @param matchedRightPoints points from right image in matrix of size 3 X
     * nPoints. row 0 is x, row 1 is y, row 2 is all 1's
     * @return array of size[2][matchedLeftPoints.numColumns()]
     */
    public double[][] calculateDistancesFromEpipolar(
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

        double[][] distances = new double[2][n];
        distances[0] = new double[n];
        distances[1] = new double[n];

        // F * u_1 = 3 x n
        DenseMatrix rightEpipolarLines = MatrixUtil.multiply(fm, matchedLeftPoints);

        // F^T * u_2
        DenseMatrix leftEpipolarLines = MatrixUtil.multiply(
            algorithms.matrix.MatrixUtil.transpose(fm),
            matchedRightPoints);

        float[] output = new float[2];

        for (int i = 0; i < matchedLeftPoints.numColumns(); i++) {

            calculatePerpDistFromLines(matchedLeftPoints, matchedRightPoints,
                rightEpipolarLines, leftEpipolarLines,
                i, i, output);

            distances[0][i] = output[0];
            distances[1][i] = output[1];
        }

        return distances;
    }
    
    /**
     * For use upon data that have been unit normal standardized.
     Sampson's distance measures the distance between a correspondence pair
     and its Sampson correction (Torr and Zissermann, 1997).
     Torr & Murray 1997 describe the Sampson distance further:
     "it represents the sum of squares of the algebraic residuals divided by 
     their standard deviations" and provides a first order fit to the 
     Kendall & Stuart (1983) minimization of the orthogonal distance of each 
     point to a curve/surface of the maximum likelihood solution. 
     * 
     perpendicular geometric distances...
     <pre>
     implemented from eqn 11 from Fathy, Husseina, & Tolbaa, 2017
     "Fundamental Matrix Estimation: A Study of Error Criteria"
     https://arxiv.org/pdf/1706.07886.pdf
     which references Torr and Zisserman 1997, 
     Machine Vision and Applications 9 (5), 321–333,
     " Performance characterization of fundamental matrix estimation under image 
     degradation"
     
     </pre>
     * @param fm
     * @param x1 points from left image in matrix of size 3 X
     * nPoints. row 0 is x, row 1 is y, row 2 is all 1's
     * @param x2 points from right image in matrix of size 3 X
     * nPoints. row 0 is x, row 1 is y, row 2 is all 1's
     * @return double array of length [n]
     */
    public double[] calculateEpipolarSampsonsDistanceSquared(
        DenseMatrix fm, DenseMatrix x1, DenseMatrix x2) {

        int m = fm.numRows();
        int n = x1.numColumns();
        if (m != 3 || fm.numColumns() != 3) {
            throw new IllegalArgumentException("fm must be 3x3");
        }
        if (x1.numRows() != 3 || x2.numRows() != 3) {
            throw new IllegalArgumentException("x1 and x2 must have 3 rows");
        }
        if (x2.numColumns() != n) {
            throw new IllegalArgumentException("x1 and x2 must have the same number of columns");
        }
        
        /*
        R_i = u2_i^T * F * u1_i
        
        line2_i = F * u1_i = (a2_i, b2_i, c2_i)
        
        line1_i = F^T * u2_i = (a1_i, b1_i, c1_i)
        
        (dist_i)^2 = (R_i)^2 / ( (a1_i)^2 + (b1_i)^2 + (a2_i)^2 + (b2_i)^2 ) 
        */
        
        double[] distances = new double[n];

        // 3 x n.  left points projected onto right image
        DenseMatrix fX1 = MatrixUtil.multiply(fm, x1);

        // 3 X n.  right points projected onto left image
        DenseMatrix fTX2 = MatrixUtil.multiply(algorithms.matrix.MatrixUtil.transpose(fm),
            x2);
        
        /*R_i is found in (row i, col i) of the result of u2^T * F * u1                                          \\//
        >> x2_00  x2_10  x2_20   *  f00 f01 f02  *  x1_00  x1_01
           x2_01  x2_11  x2_21      f10 f11 f12     x1_10  x1_11
                                    f20 f21 f22     x1_20  x1_21
        
        e.g. R_0 = (x2_00*f00 + x2_10*f10 + x2_20*f20) * x1_00   
                    + (x2_00*f01 + x2_10*f11 + x2_20*f21) * x1_10
                    + (x2_00*f02 + x2_10*f12 + x2_20*f22) * x1_20
        */
        DenseMatrix r = MatrixUtil.multiply(MatrixUtil.transpose(x2), fm);
        r = MatrixUtil.multiply(r, x1);

        double a1, b1, c1, a2, b2, c2, denom, ri;
        for (int i = 0; i < x1.numColumns(); i++) {
            
            a1 = fTX2.get(0, i);
            b1 = fTX2.get(1, i);
            c1 = fTX2.get(2, i);
            
            a2 = fX1.get(0, i);
            b2 = fX1.get(1, i);
            c2 = fX1.get(2, i);

            denom = a1*a1 + b1*b1 + a2*a2 + b2*b2;

            ri = r.get(i, i);
            ri *= ri;
            ri /= denom;
           
            distances[i] = ri;
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
     * @param fTu2
     * @param leftIdx
     * @param rightIdx
     * @param output an output variable to hold as element 0, the distance of
     * the right image point at rightIdx from from it's left epipolar line
     * projected into the right image. element 1 holds the reverse for the left
     * image point at leftIdx.
     */
    public void calculatePerpDistFromLines(DenseMatrix u1, DenseMatrix u2,
            DenseMatrix fu1, DenseMatrix fTu2,
            int leftIdx, int rightIdx, float[] output) {

        // see references for eqn (1) of within Fathy et al. 2017,
        // "Fundamental Matrix Estimation: A Study of Error Criteria"
        double a = fu1.get(0, leftIdx);
        double b = fu1.get(1, leftIdx);
        double c = fu1.get(2, leftIdx);

        double aplusb = Math.sqrt((a * a) + (b * b));

        //dist = (a*x + b*y + c)/sqrt(a^2 + b^2)
        double x = u2.get(0, rightIdx);
        double y = u2.get(1, rightIdx);

        double d = (a * x + b * y + c) / aplusb;

        // find the reverse distance by projection:
        double aRev = fTu2.get(0, rightIdx);
        double bRev = fTu2.get(1, rightIdx);
        double cRev = fTu2.get(2, rightIdx);

        double xL = u1.get(0, leftIdx);
        double yL = u1.get(1, leftIdx);

        double dRev = (aRev * xL + bRev * yL + cRev)
                / Math.sqrt(aRev * aRev + bRev * bRev);

        output[0] = (float) dRev;
        output[1] = (float) d;
    }

    public EpipolarTransformationFit calculateError(DenseMatrix fm,
        DenseMatrix x1, DenseMatrix x2, ErrorType errorType, double threshhold) {

        if (errorType.equals(ErrorType.SAMPSONS)) {
            return calculateSampsonsError(fm, x1, x2, threshhold);
        } else {
            return calculateEpipolarDistanceError(fm, x1, x2, threshhold);
        }
    }
    
    /**
     * given the fundamental matrix solution and the correspondence pairs,
     * use the errorType to calculate errors for each point, then use
     * the standard deviation and the chi-squared statistic factor to remove
     * outliers and return the results.
     * 
     * @param fm
     * @param x1
     * @param x2
     * @param errorType
     * @param chiSqStatFactor
     * @return 
     */
    public EpipolarTransformationFit calculateError2(DenseMatrix fm,
            DenseMatrix x1, DenseMatrix x2, ErrorType errorType, 
            double chiSqStatFactor) {

        if (errorType.equals(ErrorType.SAMPSONS)) {
            return calculateSampsonsError2(fm, x1, x2, chiSqStatFactor);
        } else {
            return calculateEpipolarDistanceError2(fm, x1, x2, chiSqStatFactor);
        }
    }
    
    /**
     * performs repmat(X(rowToTransposeAndReplicate,:)',1,3)
     * @param x matrix of size 3 X N.
     * @param rowToTransposeAndReplicate
     * @return matrix of size N X 3
     */
    private double[][] repmat3(double[][] x, int rowToTransposeAndReplicate) {
        int n = x[0].length;
        double[][] out = new double[n][3];
        for (int i = 0; i < n; ++i) {
            out[i] = new double[3];
            out[i][0] = x[rowToTransposeAndReplicate][i];
            out[i][1] = x[rowToTransposeAndReplicate][i];
            out[i][2] = x[rowToTransposeAndReplicate][i];
        }
        return out;
    }
    
    /**
     * performs repmat(X(rowToReplicate,:),3,1)
     * @param x matrix of size 3XN
     * @param rowToReplicate
     * @return matrix of size 3 X N
     */
    private double[][] repmat_3(double[][] x, int rowToReplicate) {
        int n = x[0].length;
        double[][] out = new double[3][n];
        out[0] = Arrays.copyOf(x[rowToReplicate], x[rowToReplicate].length);
        out[1] = Arrays.copyOf(out[0], out[0].length);
        out[2] = Arrays.copyOf(out[0], out[0].length);
        return out;
    }
    
    /**
     * make an array composed of the concatenation of a, b, and c horizontally
     * as blocks of columns.
     * @param a
     * @param b
     * @param c
     * @return concatenated array of size [a.length][a[0].length + b[0].length + c[0].length]
     */
    private double[][] concatenateAsColumnBlocks(double[][] a, double[][] b, double[][] c) {
        int n0 = a[0].length;
        int n1 = b[0].length;
        int n2 = c[0].length;
        int m = a.length;
        if (b.length != m || c.length != m) {
            throw new IllegalArgumentException("a, b, and c must have same number of rows");
        }
        int n = n0 + n1 + n2;
        int j, nc;
        double[][] out = new double[m][n];
        for (int i = 0; i < m; ++i) {
            out[i] = new double[n];
            nc = 0;
            for (j = 0; j < n0; ++j) {
                out[i][nc] = a[i][j];
                nc++;
            }
            for (j = 0; j < n1; ++j) {
                out[i][nc] = b[i][j];
                nc++;
            }
            for (j = 0; j < n2; ++j) {
                out[i][nc] = c[i][j];
                nc++;
            }
        }
        return out;
    }
    
    /**
     * make an array composed of the concatenation of a, b, and c vertically
     * as blocks of rows.
     * @param a
     * @param b
     * @param c
     * @return concatenated array of size [a.length + b.length + c.length][a[0].length]
     */
    private double[][] concatenateAsRowBlocks(double[][] a, double[][] b, double[][] c, double[][] d) {
        int m0 = a.length;
        int m1 = b.length;
        int m2 = c.length;
        int m3 = d.length;
        int n = a[0].length;
        if (b[0].length != n || c[0].length != n || d[0].length != n) {
            throw new IllegalArgumentException("a, b, c and d must have same number of columns");
        }
        int m = m0 + m1 + m2 + m3;
        int i;
        double[][] out = new double[m][n];        
        int mc = 0;
        for (i = 0; i < m0; ++i) {
            out[mc] = Arrays.copyOf(a[i], n);
            mc++;
        }
        for (i = 0; i < m1; ++i) {
            out[mc] = Arrays.copyOf(b[i], n);
            mc++;
        }
        for (i = 0; i < m2; ++i) {
            out[mc] = Arrays.copyOf(c[i], n);
            mc++;
        }
        for (i = 0; i < m3; ++i) {
            out[mc] = Arrays.copyOf(d[i], n);
            mc++;
        }
        
        return out;
    }
    
    /**
     * make an array composed of the concatenation of a, b, and c vertically
     * as rows.
     * @param a
     * @param b
     * @param c
     * @return concatenated array of size [4][a.length + b.length + c.length]
     */
    private double[][] concatenateAsRows(double[] a, double[] b, double[] c, double[] d) {
        int m0 = a.length;
        int m1 = b.length;
        int m2 = c.length;
        int m3 = d.length;
        if (m1 != m0 || m2 != m0 || m3 != m0) {
            throw new IllegalArgumentException("a, b, c and d must have same lengths");
        }
        double[][] out = new double[4][m0];        
        out[0] = Arrays.copyOf(a, m0);
        out[1] = Arrays.copyOf(b, m0);
        out[2] = Arrays.copyOf(c, m0);
        out[3] = Arrays.copyOf(d, m0);
        
        return out;
    }

    private double[] combineDistances(double[][] distances) {
        double[] c = new double[distances[0].length];
        double d;
        for (int i = 0; i < distances[0].length; ++i) {
            d = distances[0][i]*distances[0][i] + distances[1][i]*distances[1][i];
            c[i] = Math.sqrt(d);
        }
        return c;
    }
}
