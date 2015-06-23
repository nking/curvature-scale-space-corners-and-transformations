package algorithms.misc;

import algorithms.util.PairIntArray;

/**
 * class to calculate a running median of the k previous points of curveY.
 *
 * @author nichole
 */
public class AverageUtil {

    /**
     * calculate a running average of size k points of curveY.
     * runtime complexity is O(N).
     * @param curveY
     * @param kPoints
     * @return
     */
    public int[] calculateBoxCarAverage(int[] curveY, final int kPoints) {

        if (curveY == null) {
            throw new IllegalArgumentException("curveY cannot be null");
        }
        if (curveY.length < kPoints) {
            throw new IllegalArgumentException(
            "curveY.length must be equal to or greater than kPoints");
        }

        int[] output = new int[curveY.length];

        long sum = 0;

        for (int i = 0; i < kPoints; ++i) {

            sum += curveY[i];

            output[i] = Math.round((float)sum/(i + 1.f));
        }

        for (int i = kPoints; i < curveY.length; ++i) {

            int idx = i - kPoints;

            sum += (curveY[i] - curveY[idx]);

            output[i] = Math.round((float)sum/(float)kPoints);
        }
        
        return output;
    }

    /**
     * calculate a running average of size k points of curveY.
     * runtime complexity is O(N).
     * @param curve
     * @param kPoints
     * @return
     */
    public PairIntArray calculateBoxCarAverage(PairIntArray curve, final int kPoints) {

        if (curve == null) {
            throw new IllegalArgumentException("curve cannot be null");
        }
        if (curve.getN() < kPoints) {
            throw new IllegalArgumentException(
            "curve length must be equal to or greater than kPoints");
        }

        PairIntArray output = new PairIntArray(curve.getN());

        long sum = 0;

        for (int i = 0; i < kPoints; ++i) {

            sum += curve.getY(i);

            int y = Math.round((float)sum/(i + 1.f));

            int x = curve.getX(i);

            output.add(x, y);
        }

        for (int i = kPoints; i < curve.getN(); ++i) {

            int idx = i - kPoints;

            sum += (curve.getY(i) - curve.getY(idx));

            int y = Math.round((float)sum/(float)kPoints);

            int x = curve.getX(i);

            output.add(x, y);

        }

        return output;
    }

    /**
     * bin the curve by making bins of size binSize and averaging the points
     * within the bins.
     * runtime complexity is O(N).
     * @param curveY
     * @param binSize
     * @return
     */
    public int[] bin(int[] curveY, final int binSize) {

        if (curveY == null) {
            throw new IllegalArgumentException("curveY cannot be null");
        }
        if (curveY.length < binSize) {
            throw new IllegalArgumentException(
            "curveY.length must be equal to or greater than binSize");
        }

        int n = curveY.length;

        int nBins = (int)Math.ceil((float)n/(float)binSize);

        int[] output = new int[nBins];

        long sumY[] = new long[nBins];

        for (int i = 0; i < n; ++i) {

            int binNumber = i/binSize;

            sumY[binNumber] += curveY[i];
        }

        for (int i = 0; i < nBins; ++i) {
            output[i] = Math.round(sumY[i]/binSize);
        }

        return output;
    }

    public PairIntArray bin(PairIntArray curve, final int binSize) {

        if (curve == null) {
            throw new IllegalArgumentException("curve cannot be null");
        }
        if (curve.getN() < binSize) {
            throw new IllegalArgumentException(
            "curve length must be equal to or greater than binSize");
        }

        int n = curve.getN();

        int nBins = (int)Math.ceil((float)n/(float)binSize);

        PairIntArray output = new PairIntArray(nBins);

        long sumY[] = new long[nBins];

        long sumX[] = new long[nBins];

        for (int i = 0; i < n; ++i) {

            int binNumber = i/binSize;

            sumX[binNumber] += curve.getX(i);
            
            sumY[binNumber] += curve.getY(i);
        }
        
        for (int i = 0; i < nBins; ++i) {
            
            int x = Math.round(sumX[i]/binSize);
            
            int y = Math.round(sumY[i]/binSize);
            
            output.add(x, y);
        }

        return output;
    }
}
