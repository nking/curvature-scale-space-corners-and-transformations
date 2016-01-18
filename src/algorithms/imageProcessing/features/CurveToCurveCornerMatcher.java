package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Logger;

/**
 * class to invoke methods needed to solve for euclidean scale between
 * image1 and image2 using methods specific to corners on closed curves.
 *
 * @author nichole
 */
public class CurveToCurveCornerMatcher <T extends CornerRegion> {
    
    private final Logger log = Logger.getLogger(this.getClass().getName());

    private List<FeatureComparisonStat> solutionStats = null;
    
    private List<PairInt> matched1 = null;
    
    private List<PairInt> matched2 = null;

    /**
     * @return the solutionStats
     */
    public List<FeatureComparisonStat> getSolutionStats() {
        return solutionStats;
    }

    /**
     * @return the matched1
     */
    public List<PairInt> getMatched1() {
        return matched1;
    }

    /**
     * @return the matched2
     */
    public List<PairInt> getMatched2() {
        return matched2;
    }

    private enum State {
        INITIALIZED, FAILED, SOLVED
    }

    private State state = null;
    
    private final int dither;

    public CurveToCurveCornerMatcher(int dither) {
        this.dither = dither;
    }
    
    private void resetDefaults() {
        state = null;
        solutionStats = null;
        matched1 = null;
        matched2 = null;
    }
    
    public boolean matchCorners(
        final IntensityFeatures features1, final IntensityFeatures features2,
        final List<List<T>> cornerLists1, final List<List<T>> cornerLists2, 
        GreyscaleImage img1, GreyscaleImage img2, int binFactor1, int binFactor2) {
        
        int n1 = 0;
        int n2 = 0;
        for (List<T> list : cornerLists1) {
            n1 += list.size();
        }
        for (List<T> list : cornerLists2) {
            n2 += list.size();
        }

        if (n1 == 0 || n2 == 0) {
            state = State.FAILED;
            return false;
        }
        
        if (state != null) {
            resetDefaults();
        }
/*        
MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
float[] xPoints1 = new float[cornerLists1.size()];
float[] yPoints1 = new float[cornerLists1.size()];
double[][] xy1 = new double[cornerLists1.size()][2];
for (int i = 0; i < cornerLists1.size(); ++i) {
xy1[i] = curveHelper.calculateXYCentroids0(cornerLists1.get(i));
xPoints1[i] = (float)xy1[i][0];
yPoints1[i] = (float)xy1[i][1];
}
float[] xPoints2 = new float[cornerLists2.size()];
float[] yPoints2 = new float[cornerLists2.size()];
double[][] xy2 = new double[cornerLists2.size()][2];
for (int i = 0; i < cornerLists2.size(); ++i) {
xy2[i] = curveHelper.calculateXYCentroids0(cornerLists2.get(i));
xPoints2[i] = (float)xy2[i][0];
yPoints2[i] = (float)xy2[i][1];
}
StringBuilder sb = new StringBuilder("xy1:\n");
for (int i = 0; i < xy1.length; ++i) {
    sb.append(String.format("[%2d] (%3d, %3d)\n", i,
        (int)Math.round(xy1[i][0]), (int)Math.round(xy1[i][1])));
}
sb.append("xy2:\n");
for (int i = 0; i < xy2.length; ++i) {
    sb.append(String.format("[%2d] (%3d, %3d)\n", i,
        (int)Math.round(xy2[i][0]), (int)Math.round(xy2[i][1])));
}
System.out.println(sb.toString());
*/
        Map<PairInt, List<FeatureComparisonStat>> statsMap = 
            new HashMap<PairInt, List<FeatureComparisonStat>>();
        
        for (int idx1 = 0; idx1 < cornerLists1.size(); ++idx1) {

            List<T> corners1 = cornerLists1.get(idx1);

            if (corners1.size() < 1) {
                continue;
            }

            int maxNEval = Integer.MIN_VALUE;
            Integer maxNEvalIndex2 = null;
            double minCost = Double.MAX_VALUE;
            List<FeatureComparisonStat> minCostStats = null;

            for (int idx2 = 0; idx2 < cornerLists2.size(); ++idx2) {

                List<T> corners2 = cornerLists2.get(idx2);

                if (corners2.size() < 1) {
                    continue;
                }

                Integer index2 = Integer.valueOf(idx2);
        
                CornerMatcher<T> matcher = new CornerMatcher<T>(dither);
        
                boolean matched = matcher.matchCorners(features1, features2,
                    corners1, corners2, img1, img2, binFactor1, binFactor2);

                if (!matched) {
                    continue;
                }
                
                List<FeatureComparisonStat> stats2 = matcher.getSolutionStats();
                
                int nEval = stats2.size();

                if (nEval < 1) {
                    continue;
                }
                
                double cost = MiscStats.calculateCombinedIntensityStat(stats2);

                if (cost < minCost) {
                    maxNEval = nEval;
                    maxNEvalIndex2 = index2;
                    minCost = cost;
                    minCostStats = stats2;
                } else if ((maxNEval == 2) && (nEval > 3)) {
                    //TODO: may need to revise this
                    double avgCost = (cost + minCost)/2.;
                    if ((Math.abs(cost - avgCost)/(0.1*avgCost)) < 2) {
                        maxNEval = nEval;
                        maxNEvalIndex2 = index2;
                        minCost = cost;
                        minCostStats = stats2;
                    }
                }
            }

            if (minCostStats != null) {
                statsMap.put(new PairInt(idx1, maxNEvalIndex2.intValue()),
                    minCostStats);
            }
        }
        
        // might need to consider homology except that want to be able
        // to identify an object whose position has changed relative to
        // other objects due to motion or perspective change or different
        // camera perspective, etc.
        
        List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();
        for (Entry<PairInt, List<FeatureComparisonStat>> entry : statsMap.entrySet()) {
            stats.addAll(entry.getValue());
        }
        
        MiscStats.filterForDegeneracy(stats);
        
        if (stats.isEmpty()) {
            state = State.FAILED;
            return false;
        }
        
        solutionStats = stats;
    
        matched1 = new ArrayList<PairInt>();
    
        matched2 = new ArrayList<PairInt>();
        
        for (FeatureComparisonStat stat : stats) {
            PairInt p1 = stat.getImg1Point();
            PairInt p2 = stat.getImg2Point();
            matched1.add(p1);
            matched2.add(p2);
        }

        return true;
    }
    
}
