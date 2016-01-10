package algorithms.imageProcessing;

import algorithms.imageProcessing.util.MiscStats;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class CornerMatcher<T extends CornerRegion> {
    
    private final Logger log = Logger.getLogger(this.getClass().getName());

    private List<FeatureComparisonStat> solutionStats = null;
    
    private List<PairInt> matched1 = null;
    
    private List<PairInt> matched2 = null;

    private enum State {
        INITIALIZED, FAILED, SOLVED
    }

    private State state = null;
    
    private final int dither;

    public CornerMatcher(int dither) {
        this.dither = dither;
    }
    
    private void resetDefaults() {
        state = null;
        solutionStats = null;
        matched1 = null;
        matched2 = null;
    }
    
    /**
     *
     * @param features1
     * @param features2
     * @param corners1
     * @param corners2
     * @param img1 image from which to extract descriptors for features1
     * @param img2 image from which to extract descriptors for features2
     * @return
     */
    @SuppressWarnings({"unchecked"})
    public boolean matchCorners(
        final IntensityFeatures features1, final IntensityFeatures features2,
        final List<T> corners1,final List<T> corners2, GreyscaleImage img1, 
        GreyscaleImage img2, int binFactor1, int binFactor2) {

        if (state != null) {
            resetDefaults();
        }
        
        List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();
        
        FeatureMatcher featureMatcher = new FeatureMatcher();

        for (int i = 0; i < corners1.size(); ++i) {

            T region1 = corners1.get(i);

            FeatureComparisonStat best = null;
            int bestIdx2 = -1;
            FeatureComparisonStat best2nd = null;
            int bestIdx2_2nd = -1;

            for (int j = 0; j < corners2.size(); ++j) {

                T region2 = corners2.get(j);

                FeatureComparisonStat compStat = 
                    featureMatcher.ditherAndRotateForBestLocation2(
                    features1, features2, region1, region2, dither,
                    img1, img2);
                
                if ((compStat == null) ||
                    (compStat.getSumIntensitySqDiff() > compStat.getImg2PointIntensityErr())
                    ) {
                    continue;
                }
                
                if ((best2nd != null) && (compStat.getSumIntensitySqDiff() 
                    >= best2nd.getSumIntensitySqDiff())) {
                    continue;
                }
                
                if (best == null) {
                    best = compStat;
                    bestIdx2 = j;
                } else if (best2nd == null) {
                    if (compStat.getSumIntensitySqDiff() < best.getSumIntensitySqDiff()) {
                        // first becomes second and this becomes first
                        best2nd = best;
                        bestIdx2_2nd = bestIdx2;
                        best = compStat;
                        bestIdx2 = j;
                    } else {
                        best2nd = compStat;
                        bestIdx2_2nd = j;
                    }
                } else {
                    // we know it's better than 2nd best
                    if (compStat.getSumIntensitySqDiff() < best.getSumIntensitySqDiff()) {
                        // first becomes second and this becomes first
                        best2nd = best;
                        bestIdx2_2nd = bestIdx2;
                        best = compStat;
                        bestIdx2 = j;
                    } else {
                        // replaces 2nd best
                        best2nd = compStat;
                        bestIdx2_2nd = j;
                    }
                }
            }
            
            if (best == null) {
                continue;
            }

            if (best2nd == null) {
                stats.add(best);
            } else {
                
                //TODO: the ratio threshold may need to be revised.
                // see Mikolajczyk and Schmid 2005 and the Brown & Lowe paper
                
                float ratio = best.getSumIntensitySqDiff()/best2nd.getSumIntensitySqDiff();
                
                if (ratio < 0.8) {
                    stats.add(best);
                }
            }
        }
        
        MiscStats.filterForDegeneracy(stats);
        
        assignInstanceResults(stats);
        
        return !stats.isEmpty();
    }
    
    private void assignInstanceResults(List<FeatureComparisonStat> stats) {
        
        this.solutionStats = stats;
        
        matched1 = new ArrayList<PairInt>();
        
        matched2 = new ArrayList<PairInt>();
        
        for (FeatureComparisonStat stat : stats) {
            
            matched1.add(stat.getImg1Point().copy());
            
            matched2.add(stat.getImg2Point().copy());
        }
    }
    
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

}
