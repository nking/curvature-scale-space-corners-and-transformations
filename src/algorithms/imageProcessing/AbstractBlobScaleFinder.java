package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public abstract class AbstractBlobScaleFinder implements IBlobScaleFinder {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    protected boolean debug = false;

    @Override
    public void setToDebug() {
        debug = true;
    }

    @Override
    public abstract TransformationParameters solveForScale(ISegmentedImageHelper 
        img1Helper, SegmentationType type1, boolean useBinned1, 
        ISegmentedImageHelper img2Helper, SegmentationType type2, 
        boolean useBinned2, float[] outputScaleRotTransXYStDev);

    /**
     * extract the local points surrounding (x, y) on the
     * perimeter and return an object when creating descriptors.
     * Note that the perimeter is expected to be a closed curve.
     * @param theEdgeIndex
     * @param perimeter
     * @param blob
     * @return
     */
    protected BlobPerimeterRegion extractBlobPerimeterRegion(int theEdgeIndex, 
        CurvatureScaleSpaceImagePoint peakDetail, PairIntArray perimeter, 
        Set<PairInt> blob) {
        
        if (perimeter == null || (perimeter.getN() < 4)) {
            throw new IllegalArgumentException(
            "perimeter cannot be null and must have at least 4 points");
        }
        
        if (blob == null) {
            throw new IllegalArgumentException("blob cannot be null");
        }
        
        //TODO: consider asserting that this is a closed curve
        // because of averaging for some peaks, sometimes
        // (x[detailIdx, y[detailIdx]) != (peakDetail.getXCoord(), peakDetail.getYCoord())
        // so for those, have to use the preceding index and then
        // the next + 1 (instead of next)
        
        int detailIdx = peakDetail.getCoordIdx();

        // inflection points are found in between extremes of curvature,
        // so detailIdx must not be one of the curve endpoints.

        int xm = peakDetail.getXCoord();
        int ym = peakDetail.getYCoord();

        int xPrev, yPrev;
        int xNext, yNext;
        
        if (detailIdx == 0) {
            
            // wrap around closed curve
            xPrev = perimeter.getX(perimeter.getN() - 1);
            yPrev = perimeter.getY(perimeter.getN() - 1);
            
            if ((xm != perimeter.getX(detailIdx)) || (ym != perimeter.getY(detailIdx))) {
                xNext = perimeter.getX(detailIdx + 2);
                yNext = perimeter.getY(detailIdx + 2);
            } else {
                xNext = perimeter.getX(detailIdx + 1);
                yNext = perimeter.getY(detailIdx + 1);
            }
            
        } else if ((xm != perimeter.getX(detailIdx)) || (ym != perimeter.getY(detailIdx))) {
            
            xPrev = perimeter.getX(detailIdx - 1);
            yPrev = perimeter.getY(detailIdx - 1);
            
            if ((detailIdx + 2) < perimeter.getN()) {
                xNext = perimeter.getX(detailIdx + 2);
                yNext = perimeter.getY(detailIdx + 2);
            } else {
                // this is a closed curve, so wrap around
                xNext = perimeter.getX(0);
                yNext = perimeter.getY(0);
            }
            
        } else {
            
            xPrev = perimeter.getX(detailIdx - 1);
            yPrev = perimeter.getY(detailIdx - 1);
            
            if ((detailIdx + 1) < perimeter.getN()) {
                xNext = perimeter.getX(detailIdx + 1);
                yNext = perimeter.getY(detailIdx + 1);
            } else {
                // this is a closed curve, so wrap around
                xNext = perimeter.getX(0);
                yNext = perimeter.getY(0);
            }
            
        }
        
        BlobPerimeterRegion region = new BlobPerimeterRegion(theEdgeIndex, 
            xPrev, yPrev, xm, ym, xNext, yNext, blob);
        
        return region;
    }

    protected String printToString(List<FeatureComparisonStat> compStats) {
        
        StringBuilder sb = new StringBuilder();
        
        for (FeatureComparisonStat compStat : compStats) {
            
            sb.append(String.format(
            " (%d,%d) (%d,%d) theta1=%.0f theta2=%.0f intSqDiff=%.1f(%.1f)", 
            compStat.getImg1Point().getX(), compStat.getImg1Point().getY(), 
            compStat.getImg2Point().getX(), compStat.getImg2Point().getY(), 
            compStat.getImg1PointRotInDegrees(), 
            compStat.getImg2PointRotInDegrees(), 
            compStat.getSumIntensitySqDiff(), 
            compStat.getImg2PointIntensityErr()));
        }
        
        return sb.toString();
    }

    protected double calculateCombinedIntensityStat(List<FeatureComparisonStat> 
        compStats) {
        
        double sum = 0;
        
        for (FeatureComparisonStat compStat : compStats) {
            sum += compStat.getSumIntensitySqDiff();
        }
        
        sum /= (double) compStats.size();
        
        return sum;
    }

    protected TransformationParameters calculateTransformation(int binFactor1, 
        SegmentationType type1, boolean useBinned1, int binFactor2, 
        SegmentationType type2, boolean useBinned2, 
        List<FeatureComparisonStat> compStats, 
        float[] outputScaleRotTransXYStDev) {
        
        assert (compStats.isEmpty() == false);
        
        removeIntensityOutliers(compStats);
        
        MatchedPointsTransformationCalculator tc = 
            new MatchedPointsTransformationCalculator();
        
        int centroidX1 = 0;
        int centroidY1 = 0;
        
        PairIntArray matchedXY1 = new PairIntArray();
        PairIntArray matchedXY2 = new PairIntArray();
        
        float[] weights = new float[compStats.size()];
        
        double sum = 0;
        
        for (int i = 0; i < compStats.size(); ++i) {
            
            FeatureComparisonStat compStat = compStats.get(i);
            
            int x1 = compStat.getImg1Point().getX() * binFactor1;
            int y1 = compStat.getImg1Point().getY() * binFactor1;
            
            matchedXY1.add(x1, y1);
            
            int x2 = compStat.getImg2Point().getX() * binFactor2;
            int y2 = compStat.getImg2Point().getY() * binFactor2;
            
            matchedXY2.add(x2, y2);
            
            weights[i] = compStat.getSumIntensitySqDiff();
            
            sum += weights[i];
        }
        
        double tot = 0;
        
        for (int i = 0; i < compStats.size(); ++i) {
            
            double div = (sum - weights[i]) / ((compStats.size() - 1) * sum);
            
            weights[i] = (float) div;
            
            tot += div;
        }
        
        assert (Math.abs(tot - 1.) < 0.03);
        
        TransformationParameters params = tc.calulateEuclidean(matchedXY1, 
            matchedXY2, weights, centroidX1, centroidY1, 
            outputScaleRotTransXYStDev);
        
        return params;
    }

    protected void removeOutliers(List<FeatureComparisonStat> compStats) {
        
        if (compStats.size() < 2) {
            return;
        }
        
        //TODO: improve w/ a more robust outlier removal
        double[] errDivInt = new double[compStats.size()];
        
        float[] weights = new float[compStats.size()];
        
        double sum = 0;
        
        for (int i = 0; i < compStats.size(); ++i) {
            
            FeatureComparisonStat compStat = compStats.get(i);
            
            weights[i] = compStat.getSumIntensitySqDiff();
            
            sum += weights[i];
            
            errDivInt[i] = compStat.getImg2PointIntensityErr() 
                / compStat.getSumIntensitySqDiff();
        }
        
        for (int i = 0; i < compStats.size(); ++i) {
            
            double div = (sum - weights[i]) / ((compStats.size() - 1) * sum);
            
            weights[i] = (float) div;
        }
        
        float[] wghtsMeanAndStDev = MiscMath.getAvgAndStDev(weights);
        
        float maxWeight = MiscMath.findMax(weights);
        
        /*
        if all stats have intensities < 5 times their errors and
        if the stdev is approx 0.15 times the mean or less, should filter here
         */
        boolean doNotFilter = true;
        if ((wghtsMeanAndStDev[1] / wghtsMeanAndStDev[0]) > 0.15) {
            doNotFilter = false;
        }
        if (doNotFilter) {
            for (int i = 0; i < errDivInt.length; ++i) {
                if (errDivInt[i] < 5) {
                    doNotFilter = false;
                    break;
                }
            }
        }
        //TODO: may need revision
        if (doNotFilter && (weights.length <= 4)) {
            return;
        }
        
        List<FeatureComparisonStat> filteredCompStats = 
            new ArrayList<FeatureComparisonStat>();
        
        for (int i = 0; i < compStats.size(); ++i) {
            float w = weights[i];
            float diffW = Math.abs(maxWeight - w);
            if (diffW < (1.5 * wghtsMeanAndStDev[1])) {
                filteredCompStats.add(compStats.get(i));
            }
        }
        
        compStats.clear();
        
        compStats.addAll(filteredCompStats);
    }

    protected float[] calculateThetaDiff(List<FeatureComparisonStat> compStats) {
        
        if (compStats == null || compStats.isEmpty()) {
            return null;
        }
        
        float[] values = new float[compStats.size()];
        
        for (int i = 0; i < compStats.size(); ++i) {
            
            FeatureComparisonStat stat = compStats.get(i);
            
            float diff = AngleUtil.getAngleDifference(
                stat.getImg1PointRotInDegrees(), stat.getImg2PointRotInDegrees());
            
            values[i] = diff;
        }
        
        return values;
    }

    protected void removeDiscrepantThetaDiff(List<FeatureComparisonStat> compStats) {
        
        if (compStats == null || compStats.isEmpty()) {
            return;
        }
        
        float[] values = calculateThetaDiff(compStats);
        
        // 20 degree wide bins
        HistogramHolder hist = Histogram.createSimpleHistogram(20.f, values, 
            Errors.populateYErrorsBySqrt(values));
        
        int yMaxIdx = MiscMath.findYMaxIndex(hist.getYHist());
        
        float thetaDiff;
        if (yMaxIdx == -1) {
            float[] thetaDiffMeanStDev = MiscMath.getAvgAndStDev(values);
            thetaDiff = thetaDiffMeanStDev[0];
        } else {
            thetaDiff = hist.getXHist()[yMaxIdx];
        }
        
        //TODO: consider a bin larger than 20 degrees... 25
        List<Integer> remove = new ArrayList<Integer>();
        
        for (int i = 0; i < values.length; ++i) {
            float diffRot = Math.abs(values[i] - thetaDiff);
            if (diffRot > 20) {
                remove.add(Integer.valueOf(i));
            }
        }
        
        for (int i = remove.size() - 1; i > -1; --i) {
            int idx = remove.get(i);
            compStats.remove(idx);
        }
    }

    protected float[] calcIntensitySSDMeanAndStDev(List<FeatureComparisonStat> compStats) {
        
        int n = compStats.size();
        
        float[] ssds = new float[n];
        
        for (int i = 0; i < n; ++i) {
            
            FeatureComparisonStat stat = compStats.get(i);
            
            ssds[i] = stat.getSumIntensitySqDiff();
        }
        
        float[] meanStDv = MiscMath.getAvgAndStDev(ssds);
        
        return meanStDv;
    }

    protected void removeIntensityOutliers(List<FeatureComparisonStat> compStats) {
        
        if (compStats.size() < 3) {
            return;
        }
        
        int n = compStats.size();
        
        float[] meanStDv = calcIntensitySSDMeanAndStDev(compStats);
        
        List<Integer> rm = new ArrayList<Integer>();
        
        for (int i = 0; i < n; ++i) {
            
            FeatureComparisonStat stat = compStats.get(i);
            
            float diff = Math.abs(stat.getSumIntensitySqDiff() - meanStDv[0]);
            
            if (diff > (1.25 * meanStDv[1])) {
                rm.add(Integer.valueOf(i));
            }
        }
        
        for (int i = rm.size() - 1; i > -1; --i) {
            
            int idx = rm.get(i).intValue();
            
            compStats.remove(idx);
        }
    }

    protected float calculateDiffThetaMean(List<FeatureComparisonStat> 
        comparisonStats) {
        
        float[] values = calculateThetaDiff(comparisonStats);
        
        if (values == null || values.length == 0) {
            return Float.POSITIVE_INFINITY;
        }
        
        double sum = 0;
        
        for (int i = 0; i < values.length; ++i) {
            sum += values[i];
        }
        
        return (float) (sum / ((float) values.length));
    }
    
}
