package algorithms.imageProcessing.util;

import algorithms.imageProcessing.FeatureComparisonStat;
import algorithms.imageProcessing.FeatureMatcher;
import algorithms.imageProcessing.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.TransformationParameters;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import com.climbwithyourfeet.clustering.DTClusterFinder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
        
        return combinedParams;
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

}
