package algorithms.imageProcessing.util;

import algorithms.imageProcessing.FeatureComparisonStat;
import algorithms.imageProcessing.FeatureMatcher;
import algorithms.imageProcessing.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.TransformationParameters;
import algorithms.imageProcessing.Transformer;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import com.climbwithyourfeet.clustering.DTClusterFinder;
import java.util.ArrayList;
import java.util.Arrays;
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
    
    public static int[] filterForScaleAndRotation(
        List<TransformationParameters> paramsList, int shiftRotation) {
        
        if (paramsList.size() < 4) {
            int[] indexes = new int[paramsList.size()];
            for (int i = 0; i < paramsList.size(); ++i) {
                indexes[i] = i;
            }
            return indexes;
        }
        
        float minRD = Float.MAX_VALUE;
        float maxRD = Float.MIN_VALUE;
        float[] rotations = new float[paramsList.size()];
        for (int i = 0; i < paramsList.size(); ++i) {
            float r = shiftRotation + paramsList.get(i).getRotationInDegrees();
            if (r >= 360) {
                r = 360 - r;
            }
            rotations[i] = r;
            if (r < minRD) {
                minRD = r;
            }
            if (r > maxRD) {
                maxRD = r;
            }
        }

        float minS = Float.MAX_VALUE;
        float maxS = Float.MIN_VALUE;
        float[] scales = new float[paramsList.size()];
        for (int i = 0; i < paramsList.size(); ++i) {
            scales[i] = paramsList.get(i).getScale();
            if (scales[i] < minS) {
                minS = scales[i];
            }
            if (scales[i] > maxS) {
                maxS = scales[i];
            }
        }
        int diffS = (int)(maxS - minS);
        if (diffS == 0) {
            diffS = 1;
        }
        // cluster finder scales by O(dimension * lg2(dimension)), so start at 0
        float factor = (int)(maxRD - minRD)/diffS;
        float offset = -1 * (minS * factor);
        maxS = Float.MIN_VALUE;
        for (int i = 0; i < scales.length; ++i) {
            float s = scales[i];
            scales[i] = (s * factor) + offset;
            if (scales[i] > maxS) {
                maxS = scales[i];
            }
        }

        Set<PairIntWithIndex> points = new HashSet<PairIntWithIndex>();

        for (int i = 0; i < scales.length; ++i) {
            float x = rotations[i];
            float y = scales[i];
            PairIntWithIndex p = new PairIntWithIndex(Math.round(x),
                Math.round(y), i);
            points.add(p);
        }

        int dimen1 = (int)(Math.ceil(maxRD) + 1);
        int dimen2 = (int)(Math.ceil(maxS) + 1);

        DTClusterFinder<PairIntWithIndex> cFinder 
            = new DTClusterFinder<PairIntWithIndex>(points, dimen1, dimen2);
        //cFinder.setToDebug();
        cFinder.calculateCriticalDensity();

        cFinder.findClusters();

        int n = cFinder.getNumberOfClusters();
        
        if (n == 0) {
            return new int[0];
        }
        
        int maxN = Integer.MIN_VALUE;
        int maxNIdx = -1;
        for (int i = 0; i < n; ++i) {
            int ns = cFinder.getCluster(i).size();
            if (ns > maxN) {
                maxN = ns;
                maxNIdx = i;
            }
        }
        int[] indexes = new int[maxN];
        int count = 0;
        for (PairIntWithIndex p : cFinder.getCluster(maxNIdx)) {
            indexes[count] = p.getPixIndex();
            count++;
        }
        
        return indexes;
    }

    public static int[] filterForScaleAndRotationUsingHist(
        List<TransformationParameters> paramsList, int shiftRotation) {
        
        if (paramsList.size() < 4) {
            int[] indexes = new int[paramsList.size()];
            for (int i = 0; i < paramsList.size(); ++i) {
                indexes[i] = i;
            }
            return indexes;
        }
        
        float minRD = Float.MAX_VALUE;
        float maxRD = Float.MIN_VALUE;
        float[] rotations = new float[paramsList.size()];
        for (int i = 0; i < paramsList.size(); ++i) {
            float r = shiftRotation + paramsList.get(i).getRotationInDegrees();
            if (r >= 360) {
                r = 360 - r;
            }
            rotations[i] = r;
            if (r < minRD) {
                minRD = r;
            }
            if (r > maxRD) {
                maxRD = r;
            }
        }

        float minS = Float.MAX_VALUE;
        float maxS = Float.MIN_VALUE;
        float[] scales = new float[paramsList.size()];
        for (int i = 0; i < paramsList.size(); ++i) {
            scales[i] = paramsList.get(i).getScale();
            if (scales[i] < minS) {
                minS = scales[i];
            }
            if (scales[i] > maxS) {
                maxS = scales[i];
            }
        }
        int diffS = (int)(maxS - minS);
        if (diffS == 0) {
            diffS = 1;
        }
        // cluster finder scales by O(dimension * lg2(dimension)), so start at 0
        float factor = (int)(maxRD - minRD)/diffS;
        float offset = -1 * (minS * factor);
        maxS = Float.MIN_VALUE;
        for (int i = 0; i < scales.length; ++i) {
            float s = scales[i];
            scales[i] = (s * factor) + offset;
            if (scales[i] > maxS) {
                maxS = scales[i];
            }
        }

        HistogramHolder rHist = Histogram.createSimpleHistogram(25.f, rotations, 
            Errors.populateYErrorsBySqrt(rotations));
        
        HistogramHolder sHist = Histogram.createSimpleHistogram(25.f, scales, 
            Errors.populateYErrorsBySqrt(scales));
        
        int[] rYMaxIdx = MiscMath.findYMaxIndexes(rHist.getYHist());
        
        int[] sYMaxIdx = MiscMath.findYMaxIndexes(sHist.getYHist());
        
        Set<Integer> keep = new HashSet<Integer>();
        
        if (rYMaxIdx != null) {
            for (int vIdx : rYMaxIdx) {
                float v = rHist.getXHist()[vIdx];
                for (int i = 0; i < rotations.length; ++i) {
                    float diff = Math.abs(rotations[i] - v);
                    if (diff <= 20) {
                        keep.add(Integer.valueOf(i));
                    }
                }
            }
        }
        
        if (sYMaxIdx != null) {
            for (int vIdx : sYMaxIdx) {
                float v = sHist.getXHist()[vIdx];
                for (int i = 0; i < scales.length; ++i) {
                    float diff = Math.abs(scales[i] - v);
                    if (diff <= 20) {
                        keep.add(Integer.valueOf(i));
                    }
                }
            }
        }
        
        int[] indexes = new int[keep.size()];
        int count = 0;
        for (Integer index : keep) {
            indexes[count] = index.intValue();
            count++;
        }
        
        return indexes;
    }
    
    public static int[] filterForRotationUsingHist(
        List<TransformationParameters> paramsList, int shiftRotation) {
        
        if (paramsList.size() < 4) {
            int[] indexes = new int[paramsList.size()];
            for (int i = 0; i < paramsList.size(); ++i) {
                indexes[i] = i;
            }
            return indexes;
        }
        
        float minRD = Float.MAX_VALUE;
        float maxRD = Float.MIN_VALUE;
        float[] rotations = new float[paramsList.size()];
        for (int i = 0; i < paramsList.size(); ++i) {
            float r = shiftRotation + paramsList.get(i).getRotationInDegrees();
            if (r >= 360) {
                r = 360 - r;
            }
            rotations[i] = r;
            if (r < minRD) {
                minRD = r;
            }
            if (r > maxRD) {
                maxRD = r;
            }
        }

        HistogramHolder rHist = Histogram.createSimpleHistogram(25.f, rotations, 
            Errors.populateYErrorsBySqrt(rotations));
        
        int[] rYMaxIdx = MiscMath.findYMaxIndexes(rHist.getYHist());
        
        Set<Integer> keep = new HashSet<Integer>();
        
        if (rYMaxIdx != null) {
            for (int vIdx : rYMaxIdx) {
                float v = rHist.getXHist()[vIdx];
                for (int i = 0; i < rotations.length; ++i) {
                    float diff = Math.abs(rotations[i] - v);
                    if (diff <= 20) {
                        keep.add(Integer.valueOf(i));
                    }
                }
            }
        }
        
        int[] indexes = new int[keep.size()];
        int count = 0;
        for (Integer index : keep) {
            indexes[count] = index.intValue();
            count++;
        }
        
        return indexes;
    }

    public static int[] filterForScaleUsingHist(
        List<TransformationParameters> paramsList) {
        
        if (paramsList.size() < 3) {
            int[] indexes = new int[paramsList.size()];
            for (int i = 0; i < paramsList.size(); ++i) {
                indexes[i] = i;
            }
            return indexes;
        }
        
        float minS = Float.MAX_VALUE;
        float maxS = Float.MIN_VALUE;
        float[] scales = new float[paramsList.size()];
        for (int i = 0; i < paramsList.size(); ++i) {
            scales[i] = paramsList.get(i).getScale();
            if (scales[i] < minS) {
                minS = scales[i];
            }
            if (scales[i] > maxS) {
                maxS = scales[i];
            }
        }
        int diffS = (int)(maxS - minS);
        if (diffS == 0) {
            diffS = 1;
        }
        // cluster finder scales by O(dimension * lg2(dimension)), so start at 0
        float factor = 360/diffS;
        float offset = -1 * (minS * factor);
        maxS = Float.MIN_VALUE;
        for (int i = 0; i < scales.length; ++i) {
            float s = scales[i];
            scales[i] = (s * factor) + offset;
            if (scales[i] > maxS) {
                maxS = scales[i];
            }
        }

        HistogramHolder sHist = Histogram.createSimpleHistogram(25.f, scales, 
            Errors.populateYErrorsBySqrt(scales));
                
        int[] sYMaxIdx = MiscMath.findYMaxIndexes(sHist.getYHist());
        
        Set<Integer> keep = new HashSet<Integer>();
        
        if (sYMaxIdx != null) {
            for (int vIdx : sYMaxIdx) {
                float v = sHist.getXHist()[vIdx];
                for (int i = 0; i < scales.length; ++i) {
                    float diff = Math.abs(scales[i] - v);
                    if (diff <= 20) {
                        keep.add(Integer.valueOf(i));
                    }
                }
            }
        }
        
        int[] indexes = new int[keep.size()];
        int count = 0;
        for (Integer index : keep) {
            indexes[count] = index.intValue();
            count++;
        }
        
        return indexes;
    }
    
    public static int[] filterForTranslationXUsingHist(
        List<TransformationParameters> paramsList) {        
        
        if (paramsList.size() < 3) {
            int[] indexes = new int[paramsList.size()];
            for (int i = 0; i < paramsList.size(); ++i) {
                indexes[i] = i;
            }
            return indexes;
        }
        
        float minX = Float.MAX_VALUE;
        float maxX = Float.MIN_VALUE;
        float[] translationXs = new float[paramsList.size()];
        for (int i = 0; i < paramsList.size(); ++i) {
            translationXs[i] = paramsList.get(i).getTranslationX();
            if (translationXs[i] < minX) {
                minX = translationXs[i];
            }
            if (translationXs[i] > maxX) {
                maxX = translationXs[i];
            }
        }
        
        HistogramHolder tXHist = Histogram.createSimpleHistogram(translationXs, 
            Errors.populateYErrorsBySqrt(translationXs));
                
        int[] txMaxIdx = MiscMath.findYMaxIndexes(tXHist.getYHist());
        
        float binSize = tXHist.getXHist()[1] - tXHist.getXHist()[0];
        
        Set<Integer> keep = new HashSet<Integer>();
        
        if (txMaxIdx != null) {
            for (int vIdx : txMaxIdx) {
                float v = tXHist.getXHist()[vIdx];
                for (int i = 0; i < translationXs.length; ++i) {
                    float diff = Math.abs(translationXs[i] - v);
                    if (diff <= 1.1*binSize) {
                        keep.add(Integer.valueOf(i));
                    }
                }
            }
        }
        
        int[] indexes = new int[keep.size()];
        int count = 0;
        for (Integer index : keep) {
            indexes[count] = index.intValue();
            count++;
        }
        
        return indexes;
    }
    
    public static int[] filterForTranslationYUsingHist(
        List<TransformationParameters> paramsList) {
        
        if (paramsList.size() < 3) {
            int[] indexes = new int[paramsList.size()];
            for (int i = 0; i < paramsList.size(); ++i) {
                indexes[i] = i;
            }
            return indexes;
        }
        
        float minY = Float.MAX_VALUE;
        float maxY = Float.MIN_VALUE;
        float[] translationYs = new float[paramsList.size()];
        for (int i = 0; i < paramsList.size(); ++i) {
            translationYs[i] = paramsList.get(i).getTranslationY();
            if (translationYs[i] < minY) {
                minY = translationYs[i];
            }
            if (translationYs[i] > maxY) {
                maxY = translationYs[i];
            }
        }
        
        HistogramHolder tYHist = Histogram.createSimpleHistogram(translationYs, 
            Errors.populateYErrorsBySqrt(translationYs));
                
        int[] tyMaxIdx = MiscMath.findYMaxIndexes(tYHist.getYHist());
        
        float binSize = tYHist.getXHist()[1] - tYHist.getXHist()[0];
        
        Set<Integer> keep = new HashSet<Integer>();
        
        if (tyMaxIdx != null) {
            for (int vIdx : tyMaxIdx) {
                float v = tYHist.getXHist()[vIdx];
                for (int i = 0; i < translationYs.length; ++i) {
                    float diff = Math.abs(translationYs[i] - v);
                    if (diff <= 1.1*binSize) {
                        keep.add(Integer.valueOf(i));
                    }
                }
            }
        }
        
        int[] indexes = new int[keep.size()];
        int count = 0;
        for (Integer index : keep) {
            indexes[count] = index.intValue();
            count++;
        }
        
        return indexes;
    }
    
    public static int[] filterForTranslation(List<TransformationParameters> paramsList) {
        
        if (paramsList.size() < 4) {
            int[] indexes = new int[paramsList.size()];
            for (int i = 0; i < paramsList.size(); ++i) {
                indexes[i] = i;
            }
            return indexes;
        }
        
        float minTX = Float.MAX_VALUE;
        float maxTX = Float.MIN_VALUE;
        float[] transX = new float[paramsList.size()];
        for (int i = 0; i < paramsList.size(); ++i) {
            float v = paramsList.get(i).getTranslationX();
            transX[i] = v;
            if (v < minTX) {
                minTX = v;
            }
            if (v > maxTX) {
                maxTX = v;
            }
        }
        if (minTX != 0) {
            // shift the values so that '1' is the lowest value... 
            // dt in clustering runtime is O(dimension * lg_2(dimension))
            float offset = -1*minTX;
            for (int i = 0; i < transX.length; ++i) {
                transX[i] += (offset + 1);
            }
            minTX += (offset + 1);
            maxTX += (offset + 1);
        }
        
        float minTY = Float.MAX_VALUE;
        float maxTY = Float.MIN_VALUE;
        float[] transY = new float[paramsList.size()];
        for (int i = 0; i < paramsList.size(); ++i) {
            float v = paramsList.get(i).getTranslationY();
            transY[i] = v;
            if (v < minTY) {
                minTY = v;
            }
            if (v > maxTY) {
                maxTY = v;
            }
        }
        if (minTY != 0) {
            // shift the values so that '1' is the lowest value... 
            // dt in clustering runtime is O(dimension * lg_2(dimension))
            float offset = -1*minTY;
            for (int i = 0; i < transY.length; ++i) {
                transY[i] += (offset + 1);
            }
            minTY += (offset + 1);
            maxTY += (offset + 1);
        }
        
        if ((maxTX > 1000) || (maxTY > 1000)) {
            //TODO: implement a histogram version
        }
        
        Set<PairIntWithIndex> points = new HashSet<PairIntWithIndex>();

        for (int i = 0; i < transX.length; ++i) {
            float x = transX[i];
            float y = transY[i];
            PairIntWithIndex p = new PairIntWithIndex(Math.round(x),
                Math.round(y), i);
            points.add(p);
        }

        int dimen1 = (int)(Math.ceil(maxTX) + 1);
        int dimen2 = (int)(Math.ceil(maxTY) + 1);

        DTClusterFinder<PairIntWithIndex> cFinder 
            = new DTClusterFinder<PairIntWithIndex>(points, dimen1, dimen2);
        //cFinder.setToDebug();
        cFinder.setThreshholdFactor(0.75f);
        cFinder.calculateCriticalDensity();

        cFinder.findClusters();

        int n = cFinder.getNumberOfClusters();
        
        if (n == 0) {
            return new int[0];
        }
        
        int maxN = Integer.MIN_VALUE;
        int maxNIdx = -1;
        for (int i = 0; i < n; ++i) {
            int ns = cFinder.getCluster(i).size();
            if (ns > maxN) {
                maxN = ns;
                maxNIdx = i;
            }
        }
        int[] indexes = new int[maxN];
        int count = 0;
        for (PairIntWithIndex p : cFinder.getCluster(maxNIdx)) {
            indexes[count] = p.getPixIndex();
            count++;
        }
        
        return indexes;
    }

    public static boolean standardDeviationsAreSmall(TransformationParameters params) {

        if (params.getStandardDeviations() == null) {
            return false;
        }
        
        // calculation
        float tS = (params.getStandardDeviations()[0] / params.getScale());
        float tR = (float) (2. * Math.PI / params.getStandardDeviations()[1]);

        // consider comparing stdev in translations to a fraction of the image
        int tTx = Math.round(params.getStandardDeviations()[2]);
        int tTy = Math.round(params.getStandardDeviations()[3]);

        float tXConstraint = 20;
        float tYConstraint = 20;
        if (params.getNumberOfPointsUsed() < 3) {
            tXConstraint = 10;
            tYConstraint = 10;
        }

        //TODO: review these limits
        if ((tS < 0.2) && (tR >= 18.) && (tTx < tXConstraint)
            && (tTy < tYConstraint)) {

            return true;
        }

        return false;
    }

    /**
     * comparing each parameter to others and keeping only if there are other
     * similar params and among the similar, only returning average of similar.
     * 
     * @param paramsMap
     * @return 
     */
    public static List<TransformationParameters> filterToSimilarParamSets(
        Map<PairInt, TransformationParameters> paramsMap, int binFactor1,
        int binFactor2) {
        
        int transTol = 30;
        if (binFactor1 != 1 || binFactor2 != 1) {
            transTol /= ((binFactor1 + binFactor2)/2.f);
        }
        
        List<Set<TransformationParameters>> similarSets = 
            new ArrayList<Set<TransformationParameters>>();
        
        Set<PairInt> alreadyCombined = new HashSet<PairInt>();
        
        for (Map.Entry<PairInt, TransformationParameters> entry : paramsMap.entrySet()) {
            
            if (alreadyCombined.contains(entry.getKey())) {
                continue;
            }
                        
            Set<TransformationParameters> set = new HashSet<TransformationParameters>();
            
            TransformationParameters params0 = entry.getValue();
            
            for (Map.Entry<PairInt, TransformationParameters> entry2 : paramsMap.entrySet()) {
                if (entry2.getKey().equals(entry.getKey()) || 
                    alreadyCombined.contains(entry2.getKey())) {
                    continue;
                }
                TransformationParameters compare = entry2.getValue();
                float diff = AngleUtil.getAngleDifference(
                    params0.getRotationInDegrees(), compare.getRotationInDegrees());
                if (Math.abs(diff) > 20) {
                    continue;
                }
                float avg = (params0.getScale() + compare.getScale())/2.f;
                if (Math.abs(params0.getScale() - compare.getScale()) > 0.2*avg) {
                    continue;
                }
                                
                //TODO: may need to revise translation tolerance...
                double diffX = Math.abs(params0.getTranslationX()- compare.getTranslationX());
                double diffY = Math.abs(params0.getTranslationY()- compare.getTranslationY());
                if ((diffX > transTol) || (diffY > transTol)){
                    continue;
                }
                if (!alreadyCombined.contains(entry.getKey())) {
                    alreadyCombined.add(entry.getKey());
                }
                alreadyCombined.add(entry2.getKey());
                
                set.add(params0);
                set.add(compare);
            }
            if (!set.isEmpty()) {
                similarSets.add(set);
            }
        }
                  
        // use histograms to remove translation outliers in similar sets
        
        List<TransformationParameters> combinedParams = 
            new ArrayList<TransformationParameters>();

        for (Set<TransformationParameters> similarSet : similarSets) {
            
            List<TransformationParameters> similar = 
                new ArrayList<TransformationParameters>(similarSet);
            int[] keep = filterForTranslationXUsingHist(similar);
            filter(similar, keep);
            keep = filterForTranslationYUsingHist(similar);
            filter(similar, keep);
            
            if (similar.isEmpty()) {
                continue;
            } else if (similar.size() == 1) {
                combinedParams.add(similar.get(0));
                continue;
            }
            
            //average them
            MatchedPointsTransformationCalculator
                tc = new MatchedPointsTransformationCalculator();
            TransformationParameters params =  tc.averageWithoutRemoval(similar);
            combinedParams.add(params);
        }
        
        // sometimes, the averaged params in similarSets are similar to one
        // another too, so one more round of checking for similarity
        List<TransformationParameters> combinedParams2 = new ArrayList<TransformationParameters>();
        
        Set<TransformationParameters> alreadyCombined2 = new HashSet<TransformationParameters>();
        
        for (TransformationParameters params0 : combinedParams) {
            if (alreadyCombined2.contains(params0)) {
                continue;
            }
            Set<TransformationParameters> set = new HashSet<TransformationParameters>();
            set.add(params0);
            for (TransformationParameters compare : combinedParams) {
                if (compare.equals(params0)) {
                    continue;
                }
                if (alreadyCombined2.contains(compare)) {
                    continue;
                }
                float diff = AngleUtil.getAngleDifference(
                    params0.getRotationInDegrees(), compare.getRotationInDegrees());
                if (Math.abs(diff) > 20) {
                    continue;
                }
                float avg = (params0.getScale() + compare.getScale())/2.f;
                if (Math.abs(params0.getScale() - compare.getScale()) > 0.2*avg) {
                    continue;
                }
                //TODO: may need to revise translation tolerance...
                double diffX = Math.abs(params0.getTranslationX() - compare.getTranslationX());
                double diffY = Math.abs(params0.getTranslationY() - compare.getTranslationY());
                if ((diffX > transTol) || (diffY > transTol)) {
                    continue;
                }
                set.add(compare);
            }
            alreadyCombined2.addAll(set);
            if (set.size() > 1) {
                MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();
                TransformationParameters params =  tc.averageWithoutRemoval(
                    new ArrayList<TransformationParameters>(set));
                combinedParams2.add(params);
            } else {
                combinedParams2.addAll(set);
            }
        }
        
        return combinedParams2;
    }
    
    /**
     * comparing each parameter to others and if there are other
     * similar params, averaging those.
     * 
     * @param paramsMap
     * @param binFactor1
     * @param binFactor2
     * @return 
     */
    public static List<TransformationParameters> filterToSimilarParamSets2(
        Map<PairInt, TransformationParameters> paramsMap, int binFactor1,
        int binFactor2) {
        
        int transTol = 30;
        if (binFactor1 != 1 || binFactor2 != 1) {
            transTol /= ((binFactor1 + binFactor2)/2.f);
        }
        
        List<TransformationParameters> combinedParams = 
            new ArrayList<TransformationParameters>();
        
        List<Set<TransformationParameters>> similarSets = 
            new ArrayList<Set<TransformationParameters>>();
        
        Set<PairInt> alreadyCombined = new HashSet<PairInt>();
        
        for (Map.Entry<PairInt, TransformationParameters> entry : paramsMap.entrySet()) {
            
            if (alreadyCombined.contains(entry.getKey())) {
                continue;
            }
                      
            Set<TransformationParameters> set = new HashSet<TransformationParameters>();
            
            TransformationParameters params0 = entry.getValue();
            
            for (Map.Entry<PairInt, TransformationParameters> entry2 : paramsMap.entrySet()) {
                if (entry2.getKey().equals(entry.getKey()) || 
                    alreadyCombined.contains(entry2.getKey())) {
                    continue;
                }
                TransformationParameters compare = entry2.getValue();
                float diff = AngleUtil.getAngleDifference(
                    params0.getRotationInDegrees(), compare.getRotationInDegrees());
                if (Math.abs(diff) > 20) {
                    continue;
                }
                float avg = (params0.getScale() + compare.getScale())/2.f;
                if (Math.abs(params0.getScale() - compare.getScale()) > 0.2*avg) {
                    continue;
                }
                                
                //TODO: may need to revise translation tolerance...
                double diffX = Math.abs(params0.getTranslationX()- compare.getTranslationX());
                double diffY = Math.abs(params0.getTranslationY()- compare.getTranslationY());
                if ((diffX > transTol) || (diffY > transTol)){
                    continue;
                }
                if (!alreadyCombined.contains(entry.getKey())) {
                    alreadyCombined.add(entry.getKey());
                }
                alreadyCombined.add(entry2.getKey());
                
                set.add(params0);
                set.add(compare);
            }
            if (set.isEmpty()) {
                combinedParams.add(params0);
            } else {
                similarSets.add(set);
            }
        }
        
        // use histograms to remove translation outliers in similar sets

        for (Set<TransformationParameters> similarSet : similarSets) {
            
            List<TransformationParameters> similar = 
                new ArrayList<TransformationParameters>(similarSet);
            int[] keep = filterForTranslationXUsingHist(similar);
            Set<Integer> unique = MiscStats.indexesNotPresent(keep, similar.size());
            for (Integer ind : unique) {
                combinedParams.add(similar.get(ind.intValue()));
            }
            
            filter(similar, keep);
            keep = filterForTranslationYUsingHist(similar);
            unique = MiscStats.indexesNotPresent(keep, similar.size());
            for (Integer ind : unique) {
                combinedParams.add(similar.get(ind.intValue()));
            }
            filter(similar, keep);
            
            if (similar.isEmpty()) {
                continue;
            } else if (similar.size() == 1) {
                combinedParams.add(similar.get(0));
                continue;
            }
            
            //average them
            MatchedPointsTransformationCalculator
                tc = new MatchedPointsTransformationCalculator();
            TransformationParameters params =  tc.averageWithoutRemoval(similar);
            combinedParams.add(params);
        }
        
        return combinedParams;
    }
    
    public static Set<Integer> indexesNotPresent(int[] indexes, int totalNumber) {
        int n = totalNumber = indexes.length;
        Set<Integer> set = new HashSet<Integer>();
        for (int index : indexes) {
            set.add(Integer.valueOf(index));
        }
        Set<Integer> diff = new HashSet<Integer>();
        for (int i = 0; i < totalNumber; ++i) {
            if (!set.contains(Integer.valueOf(i))) {
                diff.add(Integer.valueOf(i));
            }
        }
        return diff;
    }

    private static void filter(List<TransformationParameters> params, 
        int[] indexesToKeep) {
                
        if (indexesToKeep.length < 2) {
            return;
        }
        
        List<TransformationParameters> paramsList2 = 
            new ArrayList<TransformationParameters>();
        
        for (int i = 0; i < indexesToKeep.length;++i) {
            int idx = indexesToKeep[i];
            paramsList2.add(params.get(idx));
        }
                
        params.clear();
        params.addAll(paramsList2);
    }
    
    public static TransformationParameters calculateTransformation(int binFactor1, 
        int binFactor2, List<FeatureComparisonStat> compStats, 
        float[] outputScaleRotTransXYStDev, boolean removeIntensityOutliers) {
        
        assert (compStats.isEmpty() == false);
        
        Logger log = Logger.getLogger(MiscStats.class.getName());
                
        if (removeIntensityOutliers) {
            log.info("filter for intensity outliers");
            FeatureMatcher.removeIntensityOutliers(compStats);
        }
        
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
        
        log.info("calc euclidean transform");
        
        TransformationParameters params = tc.calulateEuclidean(matchedXY1, 
            matchedXY2, weights, centroidX1, centroidY1, 
            outputScaleRotTransXYStDev);
        
        return params;
    }

    /**
     * for a map of points whose values lie between 0 and maxValueForWrapAround,
     * find the boundaries of the range, knowing that each point in the 
     * contiguous range is separated by <= toleranceInValue.  If the contiguous
     * portion wraps around from maxValueForWrapAround to 0 or beyond,
     * the returned int[]{start, end} will have a start > end.
     * @param pointValueMap
     * @param maxValueForWrapAround
     * @param toleranceInValue
     * @return 
     */
    public static int[] determineStartEndValues(Map<PairInt, Float> pointValueMap, 
        int maxValueForWrapAround, int toleranceInValue) {
        
        int minTheta = Integer.MAX_VALUE;
        int maxTheta = Integer.MIN_VALUE;
                
        for (Map.Entry<PairInt, Float> entry : pointValueMap.entrySet()) {
            int v = Math.round(entry.getValue().floatValue());
            if (v < minTheta) {
                minTheta = v;
            }
            if (v > maxTheta) {
                maxTheta = v;
            }
        }
                
        if ((minTheta == 0) && (maxTheta == maxValueForWrapAround)) {
            
            // sort values and find where delta > toleranceInValues.
            // that's the end of wrap around and beginning of wrap around
            
            int[] values = new int[pointValueMap.size()];
            int count = 0;
            for (Map.Entry<PairInt, Float> entry : pointValueMap.entrySet()) {
                values[count] = Math.round(entry.getValue().floatValue());
                count++;
            }
            Arrays.sort(values);
            
            return determineStartEndValues(values, maxValueForWrapAround, 
                toleranceInValue);
        }
        
        return new int[]{minTheta, maxTheta};
    }
    
    /**
     * for a map of points whose values lie between 0 and maxValueForWrapAround,
     * find the boundaries of the range, knowing that each point in the 
     * contiguous range is separated by <= toleranceInValue.  If the contiguous
     * portion wraps around from maxValueForWrapAround to 0 or beyond,
     * the returned int[]{start, end} will have a start > end.
     * @param sortedValues
     * @param maxValueForWrapAround
     * @param toleranceInValue
     * @return 
     */
    public static int[] determineStartEndValues(int[] sortedValues, 
        int maxValueForWrapAround, int toleranceInValue) {
        
        int[] startEndIndexes = determineStartEndIndexes(sortedValues, 
            maxValueForWrapAround, toleranceInValue);
        
        return new int[]{sortedValues[startEndIndexes[0]], 
            sortedValues[startEndIndexes[1]]};
    }
    
    /**
     * for a map of points whose values lie between 0 and maxValueForWrapAround,
     * find the boundaries of the range, knowing that each point in the 
     * contiguous range is separated by <= toleranceInValue.  If the contiguous
     * portion wraps around from maxValueForWrapAround to 0 or beyond,
     * the returned int[]{start, end} will have a start > end.
     * @param sortedValues
     * @param maxValueForWrapAround
     * @param toleranceInValue
     * @return 
     */
    public static int[] determineStartEndIndexes(int[] sortedValues, 
        int maxValueForWrapAround, int toleranceInValue) {
        
        int minTheta = sortedValues[0];
        int maxTheta = sortedValues[sortedValues.length - 1];
         
        if ((minTheta == 0) && (maxTheta == maxValueForWrapAround)) {
            
            // find where delta > toleranceInValues.
            // that's the end of wrap around and beginning of wrap around
            
            int endIdx = -1;
            int startIdx = -1;
            for (int i = 0; i < (sortedValues.length - 1); ++i) {
                int diff = sortedValues[i + 1] - sortedValues[i];
                if (diff > (toleranceInValue + 1)) {
                    endIdx = i;
                    startIdx = i;
                    break;
                }
            }
            if (startIdx > -1 && endIdx > -1) {
                return new int[]{startIdx, endIdx};
            }
        }
        
        return new int[]{0, sortedValues.length - 1};
    }
    
    /**
     * for a map of points whose values lie between 0 and maxValueForWrapAround,
     * find the boundaries of the range, knowing that each point in the 
     * contiguous range is separated by <= toleranceInValue.  If the contiguous
     * portion wraps around from maxValueForWrapAround to 0 or beyond,
     * the returned int[]{start, end} will have a start > end.
     * @param sortedValues
     * @param maxValueForWrapAround
     * @param toleranceInValue
     * @return 
     */
    public static int[] determineStartEndIndexes(float[] sortedValues, 
        int maxValueForWrapAround, int toleranceInValue) {
        
        float minTheta = sortedValues[0];
        float maxTheta = sortedValues[sortedValues.length - 1];
         
        if ((minTheta == 0) && (maxTheta == maxValueForWrapAround)) {
            
            // find where delta > toleranceInValues.
            // that's the end of wrap around and beginning of wrap around
            
            int endIdx = -1;
            int startIdx = -1;
            for (int i = 0; i < (sortedValues.length - 1); ++i) {
                float diff = sortedValues[i + 1] - sortedValues[i];
                if (diff > (toleranceInValue + 1)) {
                    endIdx = i;
                    startIdx = i;
                    break;
                }
            }
            if (startIdx > -1 && endIdx > -1) {
                return new int[]{startIdx, endIdx};
            }
        }
        
        return new int[]{0, sortedValues.length - 1};
    }
    
    public static double[] filterStatsForTranslation(
        TransformationParameters params, List<FeatureComparisonStat> compStats, 
        float sigmaFactor) {
        
        assert(!compStats.isEmpty());
                
        if (compStats.size() < 2) {
            return null;
        }
        
        Transformer transformer = new Transformer();
        
        double sumTx = 0;
        double sumTy = 0;
        double[] diffX = new double[compStats.size()];
        double[] diffY = new double[compStats.size()];
        
        for (int i = 0; i < compStats.size(); ++i) {
            
            FeatureComparisonStat compStat = compStats.get(i);
            
            int x1 = compStat.getImg1Point().getX();
            int y1 = compStat.getImg1Point().getY();
                        
            int x2 = compStat.getImg2Point().getX();
            int y2 = compStat.getImg2Point().getY();
            
            double[] xy1Tr = transformer.applyTransformation(params, x1, y1);
            
            diffX[i] = xy1Tr[0] - x1;
            diffY[i] = xy1Tr[1] - y1;
            
            sumTx += diffX[i];
            sumTy += diffY[i];
        }
        
        double length = (double)compStats.size();
        
        double avgTx = sumTx/length;
        double avgTy = sumTy/length;
        
        double stDevX = 0;
        double stDevY = 0;
        for (int i = 0; i < diffX.length; ++i) {
            double d = diffX[i] - avgTx;
            stDevX += (d*d);
            d = diffY[i] - avgTy;
            stDevY += (d*d);
        }
        stDevX = (Math.sqrt(stDevX/(length - 1.0f)));
        stDevY = (Math.sqrt(stDevY/(length - 1.0f)));

        double sumDist = 0;
        double sumSSD = 0;
        for (int i = (compStats.size() - 1); i > -1; --i) {
                        
            double dx = Math.abs(diffX[i] - avgTx);
            double dy = Math.abs(diffY[i] - avgTy);
            
            if ((dx > (sigmaFactor * stDevX)) || (dy > (sigmaFactor * stDevY))) {
                sumDist += Math.sqrt(diffX[i]*diffX[i] + diffY[i]*diffY[i]); 
                sumSSD += compStats.get(i).getSumIntensitySqDiff();
                compStats.remove(i);
            }
        }
        length = compStats.size();
        
        double[] sumDistSSD = new double[2];
        sumDistSSD[0] = (sumDist/params.getScale())/length;
        sumDistSSD[1] = sumSSD/length;
        
        return sumDistSSD;
    }

    public static List<Integer> filterForDegeneracy(List<FeatureComparisonStat> 
        stats) {
           
        Map<PairInt, List<Integer>> pointIndexes = new HashMap<PairInt, List<Integer>>();
        
        // filter for same pt1 first
        FeatureComparisonStat[] statsCp = new FeatureComparisonStat[stats.size()];
        
        for (int i = 0; i < stats.size(); ++i) {
            FeatureComparisonStat stat = stats.get(i);
            statsCp[i] = stat;
            List<Integer> indexes = pointIndexes.get(stat.getImg1Point());
            if (indexes == null) {
                indexes = new ArrayList<Integer>();
            }
            indexes.add(Integer.valueOf(i));
            pointIndexes.put(stat.getImg1Point(), indexes);
        }
        
        for (Entry<PairInt, List<Integer>> entry : pointIndexes.entrySet()) {
            List<Integer> indexes = entry.getValue();
            if (indexes.size() < 2) {
                continue;
            }
            float minCost = Float.MAX_VALUE;
            Integer minCostIndex = null;
            for (Integer index : indexes) {
                FeatureComparisonStat stat = statsCp[index.intValue()];
                if (stat == null) {
                    continue;
                }
                float cost = stat.getSumIntensitySqDiff();
                if (cost < minCost) {
                    minCost = cost;
                    minCostIndex = index;
                }
            }
            if (minCostIndex != null) {
                for (Integer index : indexes) {
                    if (index.equals(minCostIndex)) {
                        continue;
                    }
                    statsCp[index.intValue()] = null;
                } 
            }
        }
        
        pointIndexes.clear();
        
        // filter for same point2
        
        for (int i = 0; i < statsCp.length; ++i) {
            FeatureComparisonStat stat = statsCp[i];
            if (stat == null) {
                continue;
            }
            List<Integer> indexes = pointIndexes.get(stat.getImg2Point());
            if (indexes == null) {
                indexes = new ArrayList<Integer>();
            }
            indexes.add(Integer.valueOf(i));
            pointIndexes.put(stat.getImg2Point(), indexes);
        }
        
        for (Entry<PairInt, List<Integer>> entry : pointIndexes.entrySet()) {
            List<Integer> indexes = entry.getValue();
            if (indexes.size() < 2) {
                continue;
            }
            float minCost = Float.MAX_VALUE;
            Integer minCostIndex = null;
            for (Integer index : indexes) {
                FeatureComparisonStat stat = statsCp[index.intValue()];
                if (stat == null) {
                    continue;
                }
                float cost = stat.getSumIntensitySqDiff();
                if (cost < minCost) {
                    minCost = cost;
                    minCostIndex = index;
                }
            }
            if (minCostIndex != null) {
                for (Integer index : indexes) {
                    if (index.equals(minCostIndex)) {
                        continue;
                    }
                    statsCp[index.intValue()] = null;
                } 
            }
        }
        
        stats.clear();
        // store removed indexes
        List<Integer> removed = new ArrayList<Integer>();
        
        for (int i = 0; i < statsCp.length; ++i) {
            FeatureComparisonStat stat = statsCp[i];
            if (stat == null) {
                removed.add(Integer.valueOf(i));
                continue;
            }
            stats.add(stat);
        }
        
        return removed;
    }
}
