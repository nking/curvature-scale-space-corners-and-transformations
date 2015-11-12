package algorithms.imageProcessing;

import algorithms.compGeometry.NearestPoints;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public abstract class AbstractBlobScaleFinder {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    protected boolean debug = false;

    public void setToDebug() {
        debug = true;
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
        int binFactor2, List<FeatureComparisonStat> compStats, 
        float[] outputScaleRotTransXYStDev) {
        
        assert (compStats.isEmpty() == false);
        
        FeatureMatcher.removeIntensityOutliers(compStats);
        
        if (compStats.size() < 2) {
            return null;
        }
        
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

        if (sum > 0) {
            
            double tot = 0;

            for (int i = 0; i < compStats.size(); ++i) {

                double div = (sum - weights[i]) / ((compStats.size() - 1) * sum);

                weights[i] = (float) div;

                tot += div;
            }
 
            assert(Math.abs(tot - 1.) < 0.03);
            
        } else {
            float a = 1.f/weights.length;
            Arrays.fill(weights, a);
        }
        
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
        
        if (sum > 0) {
            for (int i = 0; i < compStats.size(); ++i) {

                double div = (sum - weights[i]) / ((compStats.size() - 1) * sum);

                weights[i] = (float) div;
            }
        } else {
            float a = 1.f/weights.length;
            Arrays.fill(weights, a);
            Arrays.fill(errDivInt, 0);
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
        
            //NOTE: it's possible that the diff should be (360-d1)+d2 when d1>d2
            float diff = AngleUtil.getAngleDifference(
                stat.getImg1PointRotInDegrees(), stat.getImg2PointRotInDegrees());
            
            values[i] = diff;
        }
        
        return values;
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
    
    protected boolean rotationIsConsistent(TransformationParameters params, 
        List<FeatureComparisonStat> comparisonStats, float tolerance) {
                
        if (comparisonStats == null || comparisonStats.isEmpty()) {
            return false;
        }
        
        float rotationInDegrees = params.getRotationInDegrees();
        
        for (FeatureComparisonStat stat : comparisonStats) {
            
            float a1 = stat.getImg1PointRotInDegrees();
            
            float a2 = stat.getImg2PointRotInDegrees();
            
            float theta = AngleUtil.getAngleDifference(a1, a2);
            if (theta < 0) {
                theta += 360;
            }
            
            float diff = AngleUtil.getAngleDifference(theta, rotationInDegrees);
            
            if (Math.abs(diff) > tolerance) {
                
                if (a1 > a2) {
                    
                    theta = (360.f - a1) + a2;
                    
                    diff = AngleUtil.getAngleDifference(theta, rotationInDegrees);
                    
                    if (Math.abs(diff) > tolerance) {
                        return false;
                    }
                    
                } else if (a1 < a2) {
                    
                    theta = a2 - a1;
                    
                    diff = AngleUtil.getAngleDifference(theta, rotationInDegrees);
                    
                    if (Math.abs(diff) > tolerance) {
                        return false;
                    } 
                }
            }
        }
            
        return true;
    }
    
    public boolean areSimilar(TransformationParameters params0, 
        TransformationParameters params1, float toleranceTranslation) {
        
        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();
        
        if (tc.areSimilarByScaleAndRotation(params0, params1)) {
            
            assert(params0.getOriginX() == params1.getOriginX());
            assert(params0.getOriginY() == params1.getOriginY());
            
            float diffX = params0.getTranslationX() - params1.getTranslationX();
            if (Math.abs(diffX) > toleranceTranslation) {
                return false;
            }
            float diffY = params0.getTranslationY() - params1.getTranslationY();
            if (Math.abs(diffY) > toleranceTranslation) {
                return false;
            }

            return true;
        }
        
        return false;
    }

    protected <T extends CornerRegion> int countMaxMatchable(
        List<List<T>> corners1List, List<List<T>> corners2List) {
        
        int n1 = 0;
        int n2 = 0;
        
        for (List<T> list : corners1List) {
            n1 += list.size();
        }
        
        for (List<T> list : corners2List) {
            n2 += list.size();
        }
        
        return Math.max(n1, n2);
    }
    
    protected <T extends CornerRegion> int[] convertToXPoints(
        List<List<T>> cornersList) {
        
        int n = 0;
        for (List<T> list : cornersList) {
            n += list.size();
        }
        
        int[] x = new int[n];
        n = 0;
        for (List<T> list : cornersList) {
            for (T cr : list) {
                x[n] = cr.getX()[cr.getKMaxIdx()];
                n++;
            }
        }
        
        return x;
    }
    
    protected <T extends CornerRegion> int[] convertToYPoints(
        List<List<T>> cornersList) {
        
        int n = 0;
        for (List<T> list : cornersList) {
            n += list.size();
        }
        
        int[] y = new int[n];
        n = 0;
        for (List<T> list : cornersList) {
            for (T cr : list) {
                y[n] = cr.getY()[cr.getKMaxIdx()];
                n++;
            }
        }
        
        return y;
    }

    protected <T extends CornerRegion> int evaluate(TransformationParameters params, 
        List<List<T>> corners1List, List<List<T>> corners2List, 
        int tolTransXY) {
        
        //TODO: move this method to a class for utility methods
        
        int nMatched = 0;
        
        int[] xPoints = convertToXPoints(corners2List);
        int[] yPoints = convertToYPoints(corners2List);
        
        Transformer transformer = new Transformer();
        
        NearestPoints np = new NearestPoints(xPoints, yPoints);
        
        for (int i = 0; i < corners1List.size(); ++i) {
            
            List<T> corners1 = corners1List.get(i);
            
            for (int ii = 0; ii < corners1.size(); ++ii) {
                
                T cr = corners1.get(ii);
                
                double[] xyTr = transformer.applyTransformation(params, 
                    cr.getX()[cr.getKMaxIdx()], cr.getY()[cr.getKMaxIdx()]);
                
                Set<Integer> indexes = np.findNeighborIndexes(
                    (int)Math.round(xyTr[0]), (int)Math.round(xyTr[1]), 
                    tolTransXY);
                
                if (indexes != null && indexes.size() > 0) {
                    nMatched++;
                }
            }
        }
        
        return nMatched;
    }
    
    protected int[] convertToXPoints2(List<List<CurvatureScaleSpaceContour>> cList) {
        
        int n = 0;
        for (List<CurvatureScaleSpaceContour> list : cList) {
            n += list.size();
        }
        
        int[] x = new int[n];
        n = 0;
        for (List<CurvatureScaleSpaceContour> list : cList) {
            for (CurvatureScaleSpaceContour cr : list) {
                x[n] = cr.getPeakDetails()[0].getXCoord();
                n++;
            }
        }
        
        return x;
    }
    
    protected int[] convertToYPoints2(List<List<CurvatureScaleSpaceContour>> cList) {
        
        int n = 0;
        for (List<CurvatureScaleSpaceContour> list : cList) {
            n += list.size();
        }
        
        int[] y = new int[n];
        n = 0;
        for (List<CurvatureScaleSpaceContour> list : cList) {
            for (CurvatureScaleSpaceContour cr : list) {
                y[n] = cr.getPeakDetails()[0].getYCoord();
                n++;
            }
        }
        
        return y;
    }

}
