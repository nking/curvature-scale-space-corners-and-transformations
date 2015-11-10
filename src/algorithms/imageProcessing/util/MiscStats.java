package algorithms.imageProcessing.util;

import algorithms.imageProcessing.TransformationParameters;
import com.climbwithyourfeet.clustering.DTClusterFinder;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class MiscStats {
    
    public static int[] filterForScaleAndRotation(
        List<TransformationParameters> paramsList, int shiftRotation) {
        
        if (paramsList.size() < 4) {
            return new int[0];
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
        float offset = -1 * (minRD * factor);
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

    public static int[] filterForTranslation(List<TransformationParameters> paramsList) {
        
        if (paramsList.size() < 4) {
            return new int[0];
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

}
