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
 * Sampson's distance used in error analysis of EpipolarTransformer.
 
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

    private EpipolarFeatureTransformationFit combineErrors(EpipolarTransformationFit
        distanceErrors, EpipolarFeatureTransformationFit featureErrors) {

        /*
        in order to have the distance errors count as much as the SSD errors,
        need to scale them up or SSD down by a factor.

        will use the descriptor size and the average of the maximum SSD that
        could be calculated and the minimum and make a factor for the
        distances of SSDFactor/tolerance.
        */

        // sum square diffs / size = (d0*d0) + (d3*d3).../n
        double maxSSD = 255. * 255.;

        double distScaleFactor = (maxSSD/2.)/distanceErrors.getTolerance();

        Map<Integer, Double> indexSSDErrorsMap = new HashMap<Integer, Double>();
        Map<Integer, FeatureComparisonStat> indexFeatureMap =
            new HashMap<Integer, FeatureComparisonStat>();
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
        double costTerm2 = 1./(double)outputDistances.size();
        for (int i = 0; i < outputDistances.size(); ++i) {
            double err = outputDistances.get(i).doubleValue() * costTerm2 * costTerm2
                * distScaleFactor;
            outputDistances.set(i, Double.valueOf(err));
        }

        EpipolarFeatureTransformationFit fit =
            new EpipolarFeatureTransformationFit(
            distanceErrors.getFundamentalMatrix(),
            outputInliers, fcs,
            distanceErrors.getErrorType(), outputDistances,
            distanceErrors.getTolerance());

        return fit;
    }

    /**
     * evaluate fit for already matched point lists
     * @param fm
     * @param leftPoints
     * @param rightPoints
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

            float dist = (float)Math.sqrt(leftPtD*leftPtD + rightPtD*rightPtD);

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

        return fit;
    }

    /**
     * find the distance of the given points from their respective projective
     * epipolar lines.
     * @param fm
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

        DenseMatrix rightEpipolarLines = MatrixUtil.multiply(fm, matchedLeftPoints);

        DenseMatrix leftEpipolarLines = MatrixUtil.multiply(fm.transpose(),
            matchedRightPoints);

        float[] output = new float[2];

        for (int i = 0; i < matchedLeftPoints.numColumns(); i++) {

            calculatePerpDistFromLines(matchedLeftPoints,
                matchedRightPoints, 
                rightEpipolarLines, leftEpipolarLines,
                i, i, output);

            distances.add(output[0], output[1]);
        }

        return distances;
    }

    public void calculatePerpDistFromLines(DenseMatrix leftPoints,
        DenseMatrix rightPoints, 
        DenseMatrix epipolarLinesFromLeft,
        DenseMatrix epipolarLinesFromRight, 
        int leftIdx, int rightIdx,
        float[] output) {
        
        double a = epipolarLinesFromLeft.get(0, leftIdx);
        double b = epipolarLinesFromLeft.get(1, leftIdx);
        double c = epipolarLinesFromLeft.get(2, leftIdx);

        double aplusb = Math.sqrt((a*a) + (b*b));

        //dist = (a*x + b*y + c)/sqrt(a^2 + b^2)

        double x = rightPoints.get(0, rightIdx);
        double y = rightPoints.get(1, rightIdx);

        double d = (a*x + b*y + c)/aplusb;

        // find the reverse distance by projection:
        double aRev = epipolarLinesFromRight.get(0, rightIdx);
        double bRev = epipolarLinesFromRight.get(1, rightIdx);
        double cRev = epipolarLinesFromRight.get(2, rightIdx);

        double xL = leftPoints.get(0, leftIdx);
        double yL = leftPoints.get(1, leftIdx);

        double dRev = (aRev*xL + bRev*yL + cRev)/
            Math.sqrt((aRev*aRev + bRev*bRev));

        output[0] = (float)dRev;
        output[1] = (float)d;
    }

    public EpipolarTransformationFit calculateError(DenseMatrix fm,
        DenseMatrix x1, DenseMatrix x2, ErrorType errorType, double tolerance) {

        if (errorType.equals(ErrorType.SAMPSONS)) {
            return calculateSampsonsError(fm, x1, x2, tolerance);
        } else {
            return calculateEpipolarDistanceError(fm, x1, x2, tolerance);
        }
    }

    /**
     Return the "algebraic distance" needed for use in calculating
     Sampson's distance.
     The algebraic distance has no geometrical significance,
     (it isn't the perpendicular distance).
     * The topic is discussed in
     "The Development and Comparison of Robust Methods
      for Estimating the Fundamental Matrix" by Torr and Murray, 1997
      as Equation (2) on page 274.
      The return here contains 2 matrices which contain parts of equation
      2, formatted for use in calculating Sampson's distance.
     <pre>
     The method is adapted from code from the book
     "Multiple View Geometry in Computer Vision" by 
     Hartley and Zisserman, 2004.
       
     obtained from 
     http://www.robots.ox.ac.uk/~vgg/hzbook/

        license is MIT license

        Acknowledgements: These functions are written by: David Capel,
        Andrew Fitzgibbon, Peter Kovesi, Tomas Werner, Yoni Wexler,
        and Andrew Zisserman
     </pre>
     * @param fm 3 X 3 matrix
     * @param x1 3 X nData matrix
     * left image points in format
     * matrix of size 3 X xy.getN();
        row 0 is x
        row 1 is y
        row 2 is all 1's
     * @param x2 3 X nData matrix
     * @return algebraic distances in format of 2 dimensional array of size
     * nData X 2, where [i][0] is first dimension distance and [i][1] is 2nd
        //  (i.e. [i][0] are x distances and [i][1] are y distances.
     */
    double[][] algebraicDistance(DenseMatrix fm, DenseMatrix x1,
        DenseMatrix x2) {

        // N = size(X1,2);
        int n = x1.numColumns();
        
        if (x2.numColumns() != n || x1.numRows() != x2.numRows()) {
            throw new IllegalArgumentException("x1 must be same size as x2");
        }
        
        // (x2_i * F * x1_i^T)^2
        DenseMatrix fm2 = MatrixUtil.transpose(fm);
        DenseMatrix temp1 = MatrixUtil.multiply(MatrixUtil.transpose(x2), fm2);
        DenseMatrix temp2 = MatrixUtil.multiply(temp1, x1);
        
        
        //% d = vgg_H_algebraic_distance(H,X1,X2)
        //    % For sets of homg points X1 and X2, returns the algebraic distances
        //    %  d = (p2'_x p2'_y) * p1_w - (p1_x p1_y) * p2'_w

        //Dx = [ X1' .* repmat(X2(3,:)',1,3) , zeros(N,3) , -X1' .* repmat(X2(1,:)',1,3) ];
        //Dy = [ zeros(N,3) , X1' .* repmat(X2(3,:)',1,3) , -X1' .* repmat(X2(2,:)',1,3) ];
        //h = reshape(H',9,1);
        //d = [Dx * h , Dy * h]';

        // NOTE: '.*' is matlab notation to operate on each field
        // NOTE: the ' is mathematica syntax for conjugate transpose, a.k.a.
        //       the Hermitian. it's a matrix with signs reversed for imaginary
        //        components of complex numbers and then the matrix transposed.
        
        //Dx:
        // 0) repmat(X2(3,:)',1,3).   X2 size is 3 X nData.
        //    extract row 2 for all columns.   this is 1 X nData of all "1's"
        //    transpose X2 extract.  result is nData X 1 matrix.
        //    replicate it to new matrix once along rows, and 3 times along
        //    columns.
        //    result is a nData X 3 matrix.
        //
        // 1) X1' is X1 transposed, so is nData X 3
        // 2) .* is pointwise multiplication of the matrices, which must be same size
        
        // nData X 3
        DenseMatrix x2ConjExtr = MatlabFunctions.exRowConjRepl(x2, 2, 1, 3);
        
        DenseMatrix x1Conj = MatrixUtil.transpose(x1);
        
        // nData X 3
        DenseMatrix dXCol0 = x2ConjExtr.copy();
        for (int r = 0; r < dXCol0.numRows(); ++r) {
            for (int c = 0; c < dXCol0.numColumns(); ++c) {
                double value = x1Conj.get(r, c) * dXCol0.get(r, c);
                dXCol0.set(r, c, value);
            }
        }
        
        //nData X 3
        DenseMatrix dXCol1 = new DenseMatrix(n, 3);
        for (int r = 0; r < dXCol1.numRows(); ++r) {
            for (int c = 0; c < dXCol1.numColumns(); ++c) {
                dXCol1.set(r, c, 0.);
            }
        }
        
        // -X1' .* repmat(X2(1,:)',1,3)
        // nData X 3
        DenseMatrix x2ConjExtr0 = MatlabFunctions.exRowConjRepl(x2, 0, 1, 3);
        
        //nData X 3
        DenseMatrix dXCol2 = x2ConjExtr0.copy();
        for (int r = 0; r < dXCol2.numRows(); ++r) {
            for (int c = 0; c < dXCol2.numColumns(); ++c) {
                double value = -x1Conj.get(r, c) * dXCol2.get(r, c);
                dXCol2.set(r, c, value);
            }
        }
        
        //Dy = [ zeros(N,3) , X1' .* repmat(X2(3,:)',1,3) , -X1' .* repmat(X2(2,:)',1,3) ];
        DenseMatrix dYCol0 = dXCol1.copy();
        
        DenseMatrix dYCol1 = dXCol0.copy();
        
        DenseMatrix x2ConjExtr1 = MatlabFunctions.exRowConjRepl(x2, 1, 1, 3);
        
        DenseMatrix dYCol2 = x2ConjExtr1.copy();
        for (int r = 0; r < dYCol2.numRows(); ++r) {
            for (int c = 0; c < dYCol2.numColumns(); ++c) {
                double value = -x1Conj.get(r, c) * dYCol2.get(r, c);
                dYCol2.set(r, c, value);
            }
        }
        
        
        // h = reshape(H',9,1);
        //reshape draws all from column 0 to copy to column 0 of output, etc
        DenseMatrix fmT = MatrixUtil.transpose(fm);
        double[] h2 = new double[9];
        int count = 0;
        for (int col = 0; col < fmT.numColumns(); ++col) {
            for (int row = 0; row < fmT.numRows(); ++row) {
                h2[count] = fmT.get(row, col);
                count++;
            }
        }
        assert(count == 9);

        // d = [Dx * h , Dy * h]';
       
        assert(dXCol0.numRows() == n);
        assert(dXCol0.numColumns() == 3);
        assert(dYCol0.numRows() == n);
        assert(dYCol0.numColumns() == 3);
        
        assert(dXCol1.numRows() == n);
        assert(dXCol1.numColumns() == 3);
        assert(dYCol1.numRows() == n);
        assert(dYCol1.numColumns() == 3);
        
        assert(dXCol2.numRows() == n);
        assert(dXCol2.numColumns() == 3);
        assert(dYCol2.numRows() == n);
        assert(dYCol2.numColumns() == 3);
        
        int nc2 = dXCol0.numColumns() + dXCol1.numColumns() + dXCol2.numColumns();
        assert(nc2 == 9);
        
        // nData X 9
        double[][] dXAll = new double[dXCol0.numRows()][nc2];
        for (int r = 0; r < dXCol0.numRows(); r++) {
            for (int z = 0; z < 3; z++) {
                for (int c = 0; c < dXCol0.numColumns(); c++) {
                    int c2 = z * dXCol0.numColumns() + c;
                    switch (z) {
                        case 0:
                            dXAll[r][c2] = dXCol0.get(r, c);
                            break;
                        case 1:
                            dXAll[r][c2] = dXCol1.get(r, c);
                            break;
                        default:
                            dXAll[r][c2] = dXCol2.get(r, c);
                            break;
                    }
                }
            }
        }
        
        //nData X 9
        double[][] dYAll = new double[dYCol0.numRows()][nc2];
        for (int r = 0; r < dYCol0.numRows(); r++) {
            for (int z = 0; z < 3; z++) {
                for (int c = 0; c < dYCol0.numColumns(); c++) {
                    int c2 = z * dYCol0.numColumns() + c;
                    switch (z) {
                        case 0:
                            dYAll[r][c2] = dYCol0.get(r, c);
                            break;
                        case 1:
                            dYAll[r][c2] = dYCol1.get(r, c);
                            break;
                        default:
                            dYAll[r][c2] = dYCol2.get(r, c);
                            break;
                    }
                }
            }
        }
        
        // d = [Dx * h , Dy * h]';
        
        // nData in length each
        double[] dCol0 = MatrixUtil.multiply(dXAll, h2);
        double[] dCol1 = MatrixUtil.multiply(dYAll, h2);
        
        assert(n == dCol0.length);
        assert(n == dCol1.length);
        
        double[][] d = new double[dCol0.length][2];
        for (int r = 0; r < dCol0.length; ++r) {
            d[r][0] = dCol0[r];
            d[r][1] = dCol1[r];
        }
        
        // nData X 2
        return d;
    }

    /**
     * calculate the Sampson's error for the correspondence and given
     * fundamental matrix.
     * Sampson error is a first approximation for pairs of correspondence to their
     * closest locations such that x'Fx=0.
     * 
     <pre>
     The method is adapted from code from the book
     "Multiple View Geometry in Computer Vision" by 
     Hartley and Zisserman, 2004.
       
     obtained from 
     http://www.robots.ox.ac.uk/~vgg/hzbook/

        license is MIT license

        Acknowledgements: These functions are written by: David Capel,
        Andrew Fitzgibbon, Peter Kovesi, Tomas Werner, Yoni Wexler,
        and Andrew Zisserman
     </pre>
     
     * @param fm
     * @param x1 left image points in format
     * matrix of size 3 X xy.getN();
        row 0 is x
        row 1 is y
        row 2 is all 1's
     * @param x2
     * @param tolerance .001
     * @return 
     */
    public EpipolarTransformationFit calculateSampsonsError(DenseMatrix fm,
        DenseMatrix x1, DenseMatrix x2, double tolerance) {

        if (fm == null) {
            throw new IllegalArgumentException("fm cannot be null");
        }
        if (x1 == null) {
            throw new IllegalArgumentException("x1 cannot be null");
        }
        if (x2 == null) {
            throw new IllegalArgumentException("x2 cannot be null");
        }
        if (fm.numRows() != 3 || fm.numColumns() != 3) {
            throw new IllegalArgumentException("fm should have 3 rows and 3 columns");
        }
        if (x1.numRows() != 3 || x2.numRows() != 3) {
            throw new IllegalArgumentException("x1 and x2 must "
                + "have 3 rows");
        }
        if (x1.numColumns() != x2.numColumns()) {
            throw new IllegalArgumentException("x1 and x2 must be same sizes");
        }

        int n = x1.numColumns();
        
        /*        
        geometric error of the final solution or the 7-point sample trial,
        can be approximated by Sampson's error:
             (x2_i * F * x1_i^T)^2                 (x2_i * F * x1_i^T)^2
           ---------------------------------  +  ---------------------------
             (F*x1_i^T)_x^2 + (F*x1_i^T)_y^2     (x2_i*F)_x^2 + (x2_i*F)_y^2
        */
        
        // x = A./B divides each element of A by the corresponding element of B. 
        //      The sizes of A and B must be the same or be compatible.
        // repmat result is 3 X nData.  for this project, the 3rd row is probably always all 1's.
        //p1 = X1 ./ repmat(X1(3,:),3,1);
        //p2 = X2 ./ repmat(X2(3,:),3,1);
        DenseMatrix p1 = MatlabFunctions.exRowRepl(x1, 2, 3, 1); 
        DenseMatrix p2 = MatlabFunctions.exRowRepl(x2, 2, 3, 1);
        assert(p1.numRows() == 3);
        assert(p2.numRows() == 3);
        assert(p1.numColumns() == n);
        assert(p2.numColumns() == n);
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < n; ++col) {
                double v = x1.get(row, col)/p1.get(row, col);
                p1.set(row, col, v);
                v = x2.get(row, col)/p2.get(row, col);
                p2.set(row, col, v);
            }
        }
        
        //alg = vgg_H_algebraic_distance(H,p1,p2);
        //algebraic distances in format of 2 dimensional array of size
        //nData X 2, where [i][0] is first dimension distance and [i][1] is 2nd
        //  (i.e. [i][0] are x distances and [i][1] are y distances.
        double[][] alg = algebraicDistance(fm, p1, p2);
        
        // NOTE: matlab notation '...' makes a comment out of anything else on same
        //       line and a continuation out of next line
        //G1 = [ H(1,1) - p2(1,:) * H(3,1) ; ...
        //  H(1,2) - p2(1,:) * H(3,2) ; ...
        //  -p1(1,:) * H(3,1) - p1(2,:) * H(3,2) - H(3,3) ; ...
        //  zeros(1,N) ];
        
        //G2 = [ H(2,1) - p2(2,:) * H(3,1) ; ...
        //  H(2,2) - p2(2,:) * H(3,2) ; ...
        //  zeros(1,N) ; ...
        //  -p1(1,:) * H(3,1) - p1(2,:) * H(3,2) - H(3,3) ];
                
        DenseMatrix g1 = new DenseMatrix(4, n);
        DenseMatrix g2 = new DenseMatrix(4, n);
        for (int i = 0; i < n; ++i) {
            double v1 = fm.get(0, 0) - p2.get(0, i) * fm.get(2, 0);
            double v2 = fm.get(0, 1) - p2.get(0, i) * fm.get(2, 1);
            double v3 = -p1.get(0, i) * fm.get(2, 0) - p1.get(1, i) * fm.get(2, 1) - fm.get(2, 2);
            g1.set(0, i, v1);
            g1.set(1, i, v2);
            g1.set(2, i, v3);
            g1.set(3, i, 0);
            v1 = fm.get(1, 0) - p2.get(1, i) * fm.get(2, 0);
            v2 = fm.get(1, 1) - p2.get(1, i) * fm.get(2, 1);
            g2.set(0, i, v1);
            g2.set(1, i, v2);
            g2.set(2, i, 0);
            g2.set(3, i, v3);
        }
        
        //magG1 = sqrt(sum(G1 .* G1));
        //magG2 = sqrt(sum(G2 .* G2));
        //magG1G2 = sum(G1 .*  G2);

        // 1 X nData       
        double[] magG1 = MatlabFunctions.sumEachColumn(MatlabFunctions.piecewiseMult(g1, g1));
        double[] magG2 = MatlabFunctions.sumEachColumn(MatlabFunctions.piecewiseMult(g2, g2));
        assert(magG1.length == n);
        assert(magG2.length == n);
        for (int i = 0; i < n; ++i) {
            magG1[i] = Math.sqrt(magG1[i]);
            magG2[i] = Math.sqrt(magG2[i]);
        }
        // 1 X nData
        double[] magG1G2 = MatlabFunctions.sumEachColumn(MatlabFunctions.piecewiseMult(g1, g2));
        
        // NOTE: matlab './' :
        //     x = A./B divides each element of A by the corresponding element of B. 
        //     The sizes of A and B must be the same or be compatible.
        
        //alpha = acos( magG1G2 ./ (magG1 .* magG2) );
        double[] alpha = new double[n];
        for (int i = 0; i < n; ++i) {
            double v1 = magG1G2[i];
            double v2 = magG1[i] * magG2[i];
            double v = v1/(v2 + eps);
            if (!(v < 1.)) {
                System.out.format("v1=%.3f, v2=%.3f, v=%.3f\n", 
                        (float)v1, (float)v2, (float)v);
            }
            assert(v < 1.);
            alpha[i] = Math.acos(v);
        }
        
        //alg is ndata X 2
        //magG1 is 1 X nData 
        //D1 = alg(1,:) ./ magG1;  <==== result is 1 X nData
        //D2 = alg(2,:) ./ magG2;
        
        double[] d1 = new double[n];
        double[] d2 = new double[n];
        for (int c = 0; c < n; c++) {
            d1[c] = alg[c][0] / magG1[c];
            d2[c] = alg[c][1] / magG1[c];
        }
       
        //d = (D1.*D1 + D2.*D2 - 2 * D1 .* D2 .* cos(alpha)) ./ sin(alpha);
        double[] d = new double[n];
        
        List<Integer> outputInliers = new ArrayList<Integer>();
        List<Double> outputDistances = new ArrayList<Double>();
        
        for (int i = 0; i < n; ++i) {
            d[i] = (
                (d1[i] * d1[i]) + (d2[i] * d2[i]) 
                - (2. * d1[i] * d2[i] * Math.cos(alpha[i])/Math.sin(alpha[i])));
                        
            if (d[i] < tolerance) {
                outputInliers.add(Integer.valueOf(i));
                outputDistances.add(Double.valueOf(d[i]));
            }
        }
      
        EpipolarTransformationFit fit = new EpipolarTransformationFit(fm,
            outputInliers, ErrorType.SAMPSONS, outputDistances, tolerance);

        fit.setNMaxMatchable(x1.numColumns());

        return fit;
    }

}
