package algorithms.imageProcessing.util;

import algorithms.compGeometry.PointInPolygon;
import algorithms.imageProcessing.features.FeatureComparisonStat;
import algorithms.imageProcessing.transform.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import com.climbwithyourfeet.clustering.DTClusterFinder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class MiscStats {

    public static double calculateCombinedIntensityStat(List<FeatureComparisonStat> 
        compStats) {
        
        if (compStats.isEmpty()) {
            return Double.POSITIVE_INFINITY;
        }
        
        double sum = 0;
        
        for (FeatureComparisonStat compStat : compStats) {
            sum += compStat.getSumIntensitySqDiff();
        }
        
        sum /= (double) compStats.size();
        
        return sum;
    }
   
}
