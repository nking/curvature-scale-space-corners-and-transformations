package algorithms.compGeometry.clustering;

import algorithms.misc.MiscMath;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
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
 */
public class KMeansPlusPlusFloat {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * final solution for centers of groups (== seed centers)
     */
    protected float[] center = null;
    protected int[] numberOfPointsPerSeedCell = null;
    
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
    
    public KMeansPlusPlusFloat() {
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
     * @throws IOException
     * @throws NoSuchAlgorithmException 
     */
    public void computeMeans(final int k, final float[] values) throws IOException, 
        NoSuchAlgorithmException {
         
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
    private float[] createStartSeeds(final float[] values) throws 
        NoSuchAlgorithmException {
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(System.currentTimeMillis());
        
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
                indexOfDistOfSeeds, sr, indexes, nSeedsChosen);

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

        /*
        calculate stdev or variance of points within each seed
        calculate that solution has converged by comparing for each seed:
        that changes are very little to none compared to previous solution.
        can define this as something like change change in variation should be
        be very small, near zero.
        */
 
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
            if (diff < 0) {
                diff *= -1;
            }
            // TODO:  may want to change the critical factor
            if (diff > 0.0*seedVariances[i]) {
                allAreBelowCriticalLimit = false;
            }
            
            seedVariances[i] = sumVariance[i];
        }

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

        numberOfPointsPerSeedCell = nSumStDev;
    }
    
     /**
     *
     * @param values
     * @param seed array of pixel intensities of voronoi-like seeds
     * @throws IOException
     */
    protected int[] binPoints(final float[] values, final float[] seed) throws 
        IOException {

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
        int[] indexOfDistOfSeeds, SecureRandom sr, 
        int[] indexesAlreadyChosen, int nIndexesAlreadyChosen) {
        
        //TODO: update these notes.
        
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

        int nDistDistr = 0;
        for (int i = 0; i < distOfSeeds.length; i++) {            
            int nValues = (int)Math.ceil(distOfSeeds[i]);
            // value should be present nValues number of times
            nDistDistr += nValues;
        }
        
        if (nDistDistr < 1) {
            throw new IllegalStateException("distOfSeeds is in error: " + 
                Arrays.toString(distOfSeeds));
        }
                
        while ((chosenIndex == -1) || 
            contains(indexesAlreadyChosen, nIndexesAlreadyChosen, chosenIndex)){
            
            int chosen = sr.nextInt(nDistDistr);

            // walk thru same iteration to obtain the chosen index
            int n = 0;
            for (int i = 0; i < distOfSeeds.length; i++) {            
                int nValues = (int)Math.ceil(distOfSeeds[i]);
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
