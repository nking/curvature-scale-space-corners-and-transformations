package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

/**
 * create lists of possibly degenerately matched points between 2 images.
 * It uses the cosine similarity of the descriptors along with a threshold
 * of 0.95 to keep matches above the threshold even if there are more than
 * one.
 * @author nichole
 */
public class CosineSimilarityCornerMatcher<T extends CornerRegion> {
    
    private final Logger log = Logger.getLogger(this.getClass().getName());

    private Map<PairInt, List<PairInt>> matched1Matched2 = null;
    
    private Map<PairInt, List<FeatureComparisonStat>> solutionStats = null;

    private enum State {
        INITIALIZED, FAILED, SOLVED
    }

    private State state = null;
    
    private final int dither;

    public CosineSimilarityCornerMatcher(int dither) {
        this.dither = dither;
    }
    
    private void resetDefaults() {
        state = null;
        solutionStats = null;
        matched1Matched2 = null;
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
        GreyscaleImage img2) {

        if (state != null) {
            resetDefaults();
        }
        
        matched1Matched2 = new HashMap<PairInt, List<PairInt>>();
        solutionStats = new HashMap<PairInt, List<FeatureComparisonStat>>();
        
        FeatureMatcher featureMatcher = new FeatureMatcher();

        for (int i = 0; i < corners1.size(); ++i) {

            T region1 = corners1.get(i);

            List<FeatureComparisonStat> index1Stats = new ArrayList<FeatureComparisonStat>();
            List<PairInt> index1Matches = new ArrayList<PairInt>();

            for (int j = 0; j < corners2.size(); ++j) {

                T region2 = corners2.get(j);

                FeatureComparisonStat compStat = 
                    featureMatcher.ditherAndRotateForBestLocation3(
                    features1, features2, region1, region2, dither,
                    img1, img2);
                
                if ((compStat == null) ||
                    (compStat.getSumIntensitySqDiff() > compStat.getImg2PointIntensityErr())
                    ) {
                    continue;
                }
                
                index1Stats.add(compStat);
                index1Matches.add(new PairInt(region2.getX()[region2.getKMaxIdx()],
                    region2.getY()[region2.getKMaxIdx()]));                
            }
            
            if (index1Stats.isEmpty()) {
                continue;
            }
            
            PairInt key = new PairInt(region1.getX()[region1.getKMaxIdx()],
                region1.getY()[region1.getKMaxIdx()]);   

            matched1Matched2.put(key, index1Matches);
            solutionStats.put(key, index1Stats);
        }
        
        return !matched1Matched2.isEmpty();
    }
    
    /**
     * @return the solutionStats
     */
    public Map<PairInt, List<FeatureComparisonStat>> getSolutionStats() {
        return solutionStats;
    }

    /**
     * @return the matched1
     */
    public Map<PairInt, List<PairInt>> getMatched1Matched2() {
        return matched1Matched2;
    }

}
