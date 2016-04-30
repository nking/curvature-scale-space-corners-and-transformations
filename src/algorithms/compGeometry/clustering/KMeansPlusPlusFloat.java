package algorithms.compGeometry.clustering;

import algorithms.misc.MiscMath;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;

/**
 * k-means clustering is a method of cluster analysis to partition n
 * observations into k clusters in which each observation belongs to the cluster
 * with the nearest mean.
 * This results in a partitioning of the data space into Voronoi cells, which
 * for this single parameter analysis, is 1-D.
 * 
 * Kmeans++ calculates the initial seed centers first and then proceeds with
 * the standard Kmeans algorithm.
 * 
 * The characteristic clustered in this implementation is the intensity of the
 * pixel rather than the location so distance is the difference between 
 * intensities.  It's tailored for image segmentation.
 * 
 * Useful reading:
 * http://en.wikipedia.org/wiki/K-means_clustering
 * 
 * @author nichole
 * @deprecated 
 */
public class KMeansPlusPlusFloat {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * final solution for centers of groups (== seed centers)
     */
    protected float[] center = null;
    protected int[] numberOfPointsPerSeedCell = null;
    
    protected int[] lastImgSeedIndexes = null;
    
    /**
     * this is k and is chosen by the user
     */
    protected int nSeeds = 0;
    
    // After final iteration, standard deviations are stored in seedVariances 
    // instead of variances
    protected float[] seedVariances = null;
    
    protected final static int nMaxIter = 20;
    protected int nIter = 0;
    
    private float minValue = Float.POSITIVE_INFINITY;
    private float maxValue = Float.POSITIVE_INFINITY;
    
    private final ThreadLocalRandom sr;
    
    public KMeansPlusPlusFloat() {
         sr = ThreadLocalRandom.current();
    }
    
    protected void init(final int k, final float[] values) {
        this.nSeeds = k;
        this.nIter = 0;
        this.seedVariances = new float[nSeeds];
        
        minValue = MiscMath.findMin(values);
        maxValue = MiscMath.findMax(values);
    }
    
    /**
     * note that values may have needed some pre-processing
       to group close numbers (similarity of floating point numbers depends
       upon context so a decision about numerical resolution would be
       difficult to make here).
       * Also, note that the minimum and maximum in values are used in binning,
       * so the algorithm may need to one day offer ability to pass those in
       * as argument for some use cases.
     * @param k
     * @param values
     */
    public void computeMeans(final int k, final float[] values) {
         
        init(k, values);
        
        // starter seeds, sorted by increasing value
        float[] seeds = createStartSeeds(values);
        
        int[] imgSeedIndexes = null;

        boolean hasConverged = false;

        while (!hasConverged && (nIter < nMaxIter) ) {

            imgSeedIndexes = binPoints(values, seeds);

            seeds = calculateMeanOfSeedPoints(values, imgSeedIndexes);

            if (seeds == null) {
                nIter = 0;
                lastImgSeedIndexes = null;
                seeds = createStartSeeds(values);
                continue;
            }

            hasConverged = calculateVarianceFromSeedCenters(values, seeds, 
                imgSeedIndexes);

            nIter++;
        }

        // store final numbers
        center = seeds;

        // calculate final stats
        calculateFinalStats(values, imgSeedIndexes);
    }

    /**
     * choose seeds sequentially by distance weighted probabilities
     * 
     * @param img
     * @return 
     */
    private float[] createStartSeeds(final float[] values) {
                
        float[] seed = new float[nSeeds];
        
        int[] indexes = new int[nSeeds];

        int index = getModeIdx(values);

        int nSeedsChosen = 0;
        seed[nSeedsChosen] = values[index];
        indexes[nSeedsChosen] = index;
        
        log.fine(String.format("choose seed %d) %f", nSeedsChosen, 
            seed[nSeedsChosen]));
        
        nSeedsChosen++;
        
        for (int n = 1; n < nSeeds; n++) {

            float[] distOfSeeds = new float[values.length];
            int[] indexOfDistOfSeeds = new int[values.length];

            float minAllDist = Float.MAX_VALUE;

            for (int xyIndex = 0; xyIndex < values.length; xyIndex++) {
                
                if (contains(indexes, nSeedsChosen, xyIndex)) {
                    continue;
                }
                
                float pt = values[xyIndex];
                
                float minDist = Float.MAX_VALUE;

                for (int seedIndex = 0; seedIndex < nSeedsChosen; seedIndex++) {
                    
                    float dist = Math.abs(pt - seed[seedIndex]);
                    
                    if (dist < minDist) {
                        minDist = dist;
                    }
                }
                
                distOfSeeds[xyIndex] = minDist;
                indexOfDistOfSeeds[xyIndex] = xyIndex;

                if (minDist < minAllDist) {
                    minAllDist = minDist;
                }
            }
            
            index = chooseRandomlyFromNumbersPresentByProbability(distOfSeeds, 
                indexOfDistOfSeeds, indexes, nSeedsChosen);

            seed[nSeedsChosen] = values[index];
            indexes[nSeedsChosen] = index;

            log.fine(String.format("choose seed %d) %f", nSeedsChosen, 
                seed[nSeedsChosen]));
            
            nSeedsChosen++;
        }
        
        Arrays.sort(seed);

        return seed;
    }
    
    protected boolean contains(int[] array, int nArray, int value) {
        for (int i = 0; i < nArray; i++) {
            if (array[i] == value) {
                return true;
            }
        }
        return false;
    }
    
    /**
     * calculate the mean value of all points within a seed bin and return them
     *   as new seed bin centers.  note that if there is a bin without points
     *   in it, null is returned.
     *
     * @param values
     * @param imgSeedIndexes
     * @return
     */
    protected float[] calculateMeanOfSeedPoints(final float[] values, 
        final int[] imgSeedIndexes) {

        float[] sum = new float[nSeeds];
        int[] nSum = new int[nSeeds];

        for (int xyIndex = 0; xyIndex < values.length; xyIndex++) {

            int seedIndex = imgSeedIndexes[xyIndex];

            sum[seedIndex] += values[xyIndex];
            
            nSum[seedIndex]++;
        }

        for (int i = 0; i < nSeeds; i++) {
            
            if (nSum[i] == 0) {
                return null;
            } else {
                sum[i] /= nSum[i];
            }

            log.fine(String.format("seed mean = %d) %f number of points=%d", 
                i, sum[i], nSum[i]));
            
        }

        return sum;
    }

    /**
     * Calculate the variance of the points from their seed centers and compare 
     * results with the last iteration and return true when solution has 
     * converged.  The solution has converged if each seed's variation differs 
     * from the last iteration by less than 2 sigma.
     *
     * @param values
     * @param seed
     * @param imgSeedIndexes
     * @return
     */
    protected boolean calculateVarianceFromSeedCenters(final float[] values,
        float[] seed, int[] imgSeedIndexes) {

        if (lastImgSeedIndexes == null) {
            lastImgSeedIndexes = imgSeedIndexes;
            return false;
        }
        
        boolean hasChanged = false;
        for (int i = 0; i < values.length; ++i) {
            if (lastImgSeedIndexes[i] != imgSeedIndexes[i]) {
                hasChanged = true;
                break;
            }
        }
        if (!hasChanged) {
            return true;
        }
         
        float[] sumVariance = new float[nSeeds];
        int[] nSumVariance = new int[nSeeds];

        for (int xyIndex = 0; xyIndex < values.length; xyIndex++) {

            int seedIndex = imgSeedIndexes[xyIndex];

            float d = values[xyIndex] - seed[seedIndex];
            
            sumVariance[seedIndex] += (d * d);
            
            nSumVariance[seedIndex]++;
        }

        for (int i = 0; i < sumVariance.length; i++) {
            if ((float)nSumVariance[i] > 0) {
                sumVariance[i] /= (float)nSumVariance[i];
            } else {
                sumVariance[i] = 0;
            }
        }

        // store in the instance fields

        boolean allAreBelowCriticalLimit = true;

        for (int i = 0; i < nSeeds; i++) {
            float diff = seedVariances[i] - sumVariance[i];
            if (diff < 0 ) {
                allAreBelowCriticalLimit = false;
            }
            
            seedVariances[i] = sumVariance[i];
        }

        lastImgSeedIndexes = imgSeedIndexes;
        
        return allAreBelowCriticalLimit;
    }
    
    protected void calculateFinalStats(final float[] values, final int[] 
        imgSeedIndexes) {

        float[] sumStDev = new float[nSeeds];
        int[] nSumStDev = new int[nSeeds];

        for (int xyIndex = 0; xyIndex < values.length; xyIndex++) {

            int seedIndex = imgSeedIndexes[xyIndex];

            float d = values[xyIndex] - center[seedIndex];
            
            sumStDev[seedIndex] += (d * d);
            
            nSumStDev[seedIndex]++;
        }

        for (int i = 0; i < sumStDev.length; i++) {

            if ((float)(nSumStDev[i] - 1) > 0) {
                // divide by N-1 rather because mean was calc'ed from the data
                sumStDev[i] = 
                    (float)Math.sqrt(sumStDev[i]/(float)(nSumStDev[i] - 1.));
            } else {
                sumStDev[i] = 0;
            }
            seedVariances[i] = sumStDev[i];

            log.fine(String.format("seed %d) %f stDev=%.2f number of points=%d", 
                i, center[i], seedVariances[i], nSumStDev[i]));
            
        }

        lastImgSeedIndexes = imgSeedIndexes;
        
        numberOfPointsPerSeedCell = nSumStDev;
    }
    
     /**
     *
     * @param values
     * @param seed array of pixel intensities of voronoi-like seeds
     */
    protected int[] binPoints(final float[] values, final float[] seed) {

        if (values == null) {
            throw new IllegalArgumentException("values cannot be null");
        }
        if (seed == null) {
            throw new IllegalArgumentException("seed cannot be null");
        }
        
        //TODO: review to improve this:

        int[] seedNumber = new int[values.length];

        for (int seedIndex = 0; seedIndex < seed.length; seedIndex++) {

            float bisectorBelow = ((seedIndex - 1) > -1) ?
                ((seed[seedIndex - 1] + seed[seedIndex])/2) : minValue;
                
            float bisectorAbove = ((seedIndex + 1) > (seed.length - 1)) ?
                maxValue : ((seed[seedIndex + 1] + seed[seedIndex])/2);
                       
            for (int xyIndex = 0; xyIndex < values.length; xyIndex++) {

                float pt = values[xyIndex];

                boolean isInCell = (pt >= bisectorBelow) &&  (pt <= bisectorAbove);

                if (isInCell) {
                    seedNumber[xyIndex] = seedIndex;
                    //break;
                }
            }
        }

        return seedNumber;
    }
    
    int chooseRandomlyFromNumbersPresentByProbability(float[] distOfSeeds, 
        int[] indexOfDistOfSeeds, int[] indexesAlreadyChosen, 
        int nIndexesAlreadyChosen) {
                
        // we want to choose randomly from the indexes based upon probabilities 
        // that scale by distance
        // so create an array that represents by number, the probability of a 
        //  value.  for example, distOfSeeds={2,3,4}
        //  we'd have 
        //  distIndexDistr={0,0,1,1,1,2,2,2,2}
        // and then randomly choose from that.
        // here, we skip storing every value in a large array and instead,
        // find the value for the position once the position has been 
        // drawn randomly

        int chosenIndex = -1;

        /*
        for distOfSeeds being floats and possibly having a maximum value smaller
        than the number of bins,
        need to rescale the numbers to integers.
        If one knew the desired numerical resolution already, could pick the
        factor from that.  
        Without such knowledge the factor is Integer.MAX_VALUE/maxValue.
        Need to make sure the sum is never larger than Long.MAX_VALUE.
        
        TODO: should adjust for minValue too... <===========
        */

        float factor = (float)Integer.MAX_VALUE/maxValue;

        long nDistDistr = 0;
        for (int i = 0; i < distOfSeeds.length; i++) {            
            int nValues = Math.round(factor * distOfSeeds[i]);
            // value should be present nValues number of times
            nDistDistr += nValues;
        }
        
        if (nDistDistr < 1) {
            throw new IllegalStateException("distOfSeeds is in error: " + 
                Arrays.toString(distOfSeeds));
        }
                
        while ((chosenIndex == -1) || 
            contains(indexesAlreadyChosen, nIndexesAlreadyChosen, chosenIndex)) {
            
            long chosen = sr.nextLong(nDistDistr);

            // walk thru same iteration to obtain the chosen index
            long n = 0;
            for (int i = 0; i < distOfSeeds.length; i++) {  
                
                int nValues = Math.round(factor * distOfSeeds[i]);
                // value should be present nValues number of times
                
                if ((chosen >= n) && (chosen < (n + nValues))) {
                    chosenIndex = indexOfDistOfSeeds[i];
                    break;
                }
                n += nValues;
            }
        }

        return chosenIndex;
    }
    
    public float[] getStandardDeviationsFromCenters() {
        return this.seedVariances;
    }

    public float[] getCenters() {
        return this.center;
    }
    
    public int[] getImgPixelSeedIndexes() {
        return lastImgSeedIndexes;
    }

    public int[] getNumberOfPointsPerSeedCell() {
        return numberOfPointsPerSeedCell;
    }

    private int getModeIdx(float[] values) {
        
        //TODO: numbers given to main invocation may have needed some pre-processing
        // to group close numbers (similarity of floating point numbers depends
        // upon context so a decision about numerical resolution would be
        // difficult to make here).
        
        Map<Float, Integer> counts = new HashMap<Float, Integer>();
        int maxCounts = 0;
        int maxCountsIdx = -1;
        
        for (int idx = 0; idx < values.length; idx++) {
            float v = values[idx];
            Float key = Float.valueOf(v);
            Integer value = counts.get(key);
            int freq = (value == null) ? 1 : (value.intValue() + 1);
            counts.put(key, Integer.valueOf(freq));
            if (freq > maxCounts) {
                maxCounts = freq;
                maxCountsIdx = idx;
            }
        }
        
        return maxCountsIdx;
    }

    /**
     * @return the minValue
     */
    public float getMinValue() {
        return minValue;
    }

    /**
     * @return the maxValue
     */
    public float getMaxValue() {
        return maxValue;
    }

}
