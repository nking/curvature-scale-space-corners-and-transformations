package algorithms.compGeometry.clustering;

import algorithms.compGeometry.NearestPoints1D;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.search.KDTree;
import algorithms.search.KDTreeNode;
import algorithms.util.PairInt;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
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
 */
public class KMeansPlusPlus {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * final solution for centers of groups (== seed centers)
     */
    protected int[] center = null;
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
    
    protected int imgModeIdx = -1;
    
    private final ThreadLocalRandom sr;
    private int[] imgPixelSeedIndexes = null;
    
    public KMeansPlusPlus() {
         sr = ThreadLocalRandom.current();
    }
    
    protected void init(int k) {
        this.nSeeds = k;
        this.nIter = 0;
        this.seedVariances = new float[nSeeds];
    }
    
    /**
     * note that an internal method binPoints makes an assumption that the
     * minimum values in an img pixel is 0 and a maximum is 255.
     * @param k
     * @param img
     * @throws IOException
     * @throws NoSuchAlgorithmException 
     */
    public void computeMeans(int k, GreyscaleImage img) throws IOException, 
        NoSuchAlgorithmException {
         
        init(k);
        
        // starter seeds, sorted by increasing value
        int[] seeds = createStartSeeds(img);
        
        int[] imgSeedIndexes = null;

        boolean hasConverged = false;

        while (!hasConverged && (nIter < nMaxIter) ) {

            imgSeedIndexes = binPoints(img, seeds);

            seeds = calculateMeanOfSeedPoints(img, imgSeedIndexes);

            if (seeds == null) {
                nIter = 0;
                seeds = createStartSeeds(img);
                continue;
            }

            hasConverged = calculateVarianceFromSeedCenters(img, seeds, 
                imgSeedIndexes);

            nIter++;
        }

        // store final numbers
        center = seeds;

        // calculate final stats
        calculateFinalStats(img, imgSeedIndexes);
    }

    /**
     * choose seeds sequentially by distance weighted probabilities
     * 
     * @param img
     * @return 
     */
    private int[] createStartSeeds(GreyscaleImage img) throws 
        NoSuchAlgorithmException {
                
        int[] seed = new int[nSeeds];
        
        int[] indexes = new int[nSeeds];

        int index = getModeIdx(img);

        int nSeedsChosen = 0;
        seed[nSeedsChosen] = img.getValue(index);
        indexes[nSeedsChosen] = index;
        
        log.fine(String.format("choose seed %d) %d", nSeedsChosen, 
            seed[nSeedsChosen]));
        
        nSeedsChosen++;
        
        for (int n = 1; n < nSeeds; n++) {

            int[] distOfSeeds = new int[img.getNPixels()];
            int[] indexOfDistOfSeeds = new int[img.getNPixels()];

            populateDistanceArrays(img, seed, nSeedsChosen, distOfSeeds, indexOfDistOfSeeds);
 
            index = chooseRandomlyFromNumbersPresentByProbability(img,
                distOfSeeds, 
                indexOfDistOfSeeds, seed, indexes, nSeedsChosen);

            seed[nSeedsChosen] = img.getValue(index);
            indexes[nSeedsChosen] = index;

            log.fine(String.format("choose seed %d) %d", nSeedsChosen, 
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
     * @param img
     * @param imgSeedIndexes
     * @return
     */
    protected int[] calculateMeanOfSeedPoints(final GreyscaleImage img, 
        final int[] imgSeedIndexes) {

        int[] sum = new int[nSeeds];
        int[] nSum = new int[nSeeds];

        for (int xyIndex = 0; xyIndex < img.getNPixels(); xyIndex++) {

            int seedIndex = imgSeedIndexes[xyIndex];

            sum[seedIndex] += img.getValue(xyIndex);
            
            nSum[seedIndex]++;
        }

        for (int i = 0; i < nSeeds; i++) {
            
            if (nSum[i] == 0) {
                return null;
            } else {
                sum[i] /= nSum[i];
            }

            log.fine(String.format("seed mean = %d) %d number of points=%d", 
                i, sum[i], nSum[i]));
            
        }

        return sum;
    }

    /**
     * calculate the variance of the points from their seed centers and compare 
     * results with the last iteration and return true when solution has 
     * converged.  the solution has converged if each seed's variation differs 
     * from the last iteration by less than 2 sigma.
     *
     * @param img
     * @param seed
     * @param imgSeedIndexes
     * @return
     */
    protected boolean calculateVarianceFromSeedCenters(final GreyscaleImage img,
        int[] seed, int[] imgSeedIndexes) {

        /*
        calculate stdev or variance of points within each seed
        calculate that solution has converged by comparing for each seed:
        that changes are very little to none compared to previous solution.
        can define this as something like change change in variation should be
        be very small, near zero.
        */
 
        float[] sumVariance = new float[nSeeds];
        int[] nSumVariance = new int[nSeeds];

        for (int xyIndex = 0; xyIndex < img.getNPixels(); xyIndex++) {

            int seedIndex = imgSeedIndexes[xyIndex];

            int d = img.getValue(xyIndex) - seed[seedIndex];
            
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
    
    protected void calculateFinalStats(final GreyscaleImage img, 
        final int[] imgSeedIndexes) {

        float[] sumStDev = new float[nSeeds];
        int[] nSumStDev = new int[nSeeds];

        for (int xyIndex = 0; xyIndex < img.getNPixels(); xyIndex++) {

            int seedIndex = imgSeedIndexes[xyIndex];

            int d = img.getValue(xyIndex) - center[seedIndex];
            
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

            log.fine(String.format("seed %d) %d stDev=%.2f number of points=%d", 
                i, center[i], seedVariances[i], nSumStDev[i]));
            
        }
        
        imgPixelSeedIndexes = imgSeedIndexes;

        numberOfPointsPerSeedCell = nSumStDev;
    }
    
     /**
     *
     * @param img
     * @param seed array of pixel intensities of voronoi-like seeds
     * @throws IOException
     */
    protected int[] binPoints(final GreyscaleImage img,
        int[] seed) throws IOException {
        
        int nc = seed.length;
        int[] values = new int[nc];
        for (int i = 0; i < nc; ++i) {
            values[i] = seed[i];
        }
        
        NearestPoints1D np = new NearestPoints1D(values);
                    
        int[] seedNumber = new int[img.getNPixels()];
        
        for (int pixIdx = 0; pixIdx < img.getNPixels(); pixIdx++) {

            int v = img.getValue(pixIdx);
            
            int seedIdx = np.findClosestValueIndex(v);
                                    
            seedNumber[pixIdx] = seedIdx;
        }

        return seedNumber;
                
    }
    
    /**
     * an algorithm to choose randomly from values using a probability 
     * function based upon the value.  In other words, if the values were
     * [2, 3, 4] their presence in a bin to choose from using their indexes
     * as the bin values would be [0, 0, 1, 1, 1, 2, 2, 2, 2] and choosing
     * '2' from the bin would be more likely and '2' is the index so the result
     * would more often be values[2] = 4.
     * The additional arguments are present for uses specific to KMeansPlusPlus
     * where one doesn't want two similar seeds as a result.
     * @param distOfSeeds
     * @param indexOfDistOfSeeds
     * @param sr
     * @param indexesAlreadyChosen
     * @param nIndexesAlreadyChosen
     * @return 
     */
    int chooseRandomlyFromNumbersPresentByProbability(GreyscaleImage img,
        int[] distOfSeeds, 
        int[] indexOfDistOfSeeds, int[] valuesAlreadyChosen, 
        int[] indexesAlreadyChosen, int nIndexesAlreadyChosen) {
        
        Set<Integer> chosenValues = new HashSet<Integer>();
        for (int i = 0; i < nIndexesAlreadyChosen; ++i) {
            chosenValues.add(Integer.valueOf(valuesAlreadyChosen[i]));
        }
        
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
        
        long nDistDistr = 0;
        for (int i = 0; i < distOfSeeds.length; i++) {            
            int nValues = distOfSeeds[i];
            // value should be present nValues number of times
            nDistDistr += nValues;
        }
        
        if (nDistDistr < 1) {
            throw new IllegalStateException("distOfSeeds is in error: " + 
                Arrays.toString(distOfSeeds));
        }
                
        while ((chosenIndex == -1) || 
            contains(indexesAlreadyChosen, nIndexesAlreadyChosen, chosenIndex)
            || ((chosenIndex > -1) &&
                chosenValues.contains(
                    Integer.valueOf(img.getValue(chosenIndex))))
            ) {
            
            long chosen = sr.nextLong(nDistDistr);

            // walk thru same iteration to obtain the chosen index
            long n = 0;
            for (int i = 0; i < distOfSeeds.length; i++) {            
                int nValues = distOfSeeds[i];
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

    public int[] getCenters() {
        return this.center;
    }
    
    public int[] getImgPixelSeedIndexes() {
        return imgPixelSeedIndexes;
    }

    public int[] getNumberOfPointsPerSeedCell() {
        return numberOfPointsPerSeedCell;
    }

    private int getModeIdx(GreyscaleImage img) {
        
        if (imgModeIdx > -1) {
            return imgModeIdx;
        }
        
        Map<Integer, Integer> counts = new HashMap<Integer, Integer>();
        int maxCounts = 0;
        int maxCountsIdx = -1;
        
        for (int idx = 0; idx < img.getNPixels(); idx++) {
            int v = img.getValue(idx);
            Integer key = Integer.valueOf(v);
            Integer value = counts.get(key);
            int f = (value == null) ? 1 : (value.intValue() + 1);
            counts.put(key, Integer.valueOf(f));
            if (f > maxCounts) {
                maxCounts = f;
                maxCountsIdx = idx;
            }
        }
        
        imgModeIdx = maxCountsIdx;
        
        return maxCountsIdx;
    }

    private void populateDistanceArrays(GreyscaleImage img, int[] seeds, 
        int nSeedsChosen, int[] distOfSeeds, int[] indexOfDistOfSeeds) {
        
        NearestPoints1D np = new NearestPoints1D(seeds, nSeedsChosen);
         
        for (int pixIdx = 0; pixIdx < img.getNPixels(); pixIdx++) {
            
            int v = img.getValue(pixIdx);
            
            int seedIdx = np.findClosestValueIndex(v);
            
            int diff = v - seeds[seedIdx];

            distOfSeeds[pixIdx] = (diff * diff);

            indexOfDistOfSeeds[pixIdx] = pixIdx;
        }
    }

}
