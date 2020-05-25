package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.transform.Distances;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.imageProcessing.transform.EpipolarTransformer;
import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import no.uib.cipr.matrix.DenseMatrix;

/**
 *
 * @author nichole
 */
public class DegenerateFilter {

    /**
     * if a point maps to multiple indexes, that is, was matched to more than
     * one point in the other image, the smallest error pair of points among the conflicting
     * is kept.   NOTE, this could be improved with a bipartite matching of
     * min cost between the degenerately mapped point pairs.
     * @param xy1
     * @param xy2
     * @param outputInliers
     * @param outputDistances
     */
    public void filterForDegenerate(DenseMatrix xy1, DenseMatrix xy2, List<Integer> outputInliers, List<Double> outputDistances) {
        filterForDegenerate(xy1, outputInliers, outputDistances);
        filterForDegenerate(xy2, outputInliers, outputDistances);
    }

    public void filterForDegenerate(DenseMatrix xy1, List<Integer> outputInliers, List<Double> outputDistances) {
        Map<PairInt, List<Integer>> pointIndexes = new HashMap<PairInt, List<Integer>>();
        for (int i = 0; i < outputInliers.size(); ++i) {
            int idx = outputInliers.get(i);
            int x1 = (int) Math.round(xy1.get(0, idx));
            int y1 = (int) Math.round(xy1.get(1, idx));
            PairInt p1 = new PairInt(x1, y1);
            List<Integer> oIndexes = pointIndexes.get(p1);
            if (oIndexes == null) {
                oIndexes = new ArrayList<Integer>();
                pointIndexes.put(p1, oIndexes);
            }
            oIndexes.add(Integer.valueOf(i));
        }
        List<Integer> removeListIndexes = new ArrayList<Integer>();
        for (Map.Entry<PairInt, List<Integer>> entry : pointIndexes.entrySet()) {
            List<Integer> oIndexes = entry.getValue();
            if (oIndexes.size() < 2) {
                continue;
            }
            double minError = Double.MAX_VALUE;
            int minIdx = -1;
            for (Integer index : oIndexes) {
                int idx = index.intValue();
                double error = outputDistances.get(idx);
                if (error < minError) {
                    minError = error;
                    minIdx = idx;
                }
            }
            assert (minIdx > -1);
            for (int ii = 0; ii < oIndexes.size(); ++ii) {
                int idx = oIndexes.get(ii);
                if (idx == minIdx) {
                    continue;
                }
                removeListIndexes.add(Integer.valueOf(ii));
            }
        }
        Collections.sort(removeListIndexes);
        for (int i = removeListIndexes.size() - 1; i > -1; --i) {
            int ii = removeListIndexes.get(i);
            outputDistances.remove(ii);
            outputInliers.remove(ii);
        }
    }

    public EpipolarTransformationFit calculateErrorThenFilter(DenseMatrix fm, 
        DenseMatrix x1, DenseMatrix x2, ErrorType errorType, double tolerance) {
        //if (errorType.equals(ErrorType.SAMPSONS)) {
        //    return calculateSampsonsErrorThenFilter(fm, x1, x2, tolerance);
        //} else {
        return calculateEpipolarDistanceErrorThenFilter(fm, 
            x1, x2, tolerance);
        //}
    }
    
     /**
     * evaluate fit for already matched point lists
     * @param fm
     * @param leftPoints
     * @param rightPoints
     * @param tolerance
     * @return
     */
    public EpipolarTransformationFit calculateEpipolarDistanceErrorThenFilter(
        final DenseMatrix fm, final DenseMatrix leftPoints, final DenseMatrix rightPoints,
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

        Distances distancesObj = new Distances();
        
        PairFloatArray distances = distancesObj.calculateDistancesFromEpipolar(fm,
            leftPoints, rightPoints);
        
        List<Double> errors = new ArrayList<Double>();

        List<Integer> inlierIndexes = new ArrayList<Integer>();

        for (int i = 0; i < distances.getN(); ++i) {

            float rightPtD = distances.getX(i);

            float leftPtD = distances.getY(i);

            float dist = (float)Math.sqrt(leftPtD*leftPtD + rightPtD*rightPtD);

            if (dist > tolerance) {
                continue;
            }

            inlierIndexes.add(Integer.valueOf(i));

            errors.add(Double.valueOf(dist));
        }

        filterForDegenerate(leftPoints, rightPoints, inlierIndexes, errors);

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

}
