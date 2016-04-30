package algorithms.compGeometry.clustering;

import algorithms.compGeometry.NearestPoints1D;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.util.TrioInt;
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
    
    protected int imgModeIdx = -1;
    
    private final ThreadLocalRandom sr;
    
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
                lastImgSeedIndexes = null;
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
        
        int index = getModeIdx(img);

        seed[0] = img.getValue(index);
        
        log.fine(String.format("choose seed %d) %d", 0, seed[0]));
        
        Set<Integer> alreadyChosenValues = new HashSet<Integer>();
        
        alreadyChosenValues.add(seed[0]);

        for (int n = 1; n < nSeeds; n++) {

            int[] distOfSeeds = new int[img.getNPixels()];
            int[] indexOfDistOfSeeds = new int[img.getNPixels()];

            populateDistanceArrays(img, alreadyChosenValues, distOfSeeds, indexOfDistOfSeeds);
 
            int selectedSeedValue = chooseRandomlyFromNumbersPresentByProbability(img,
                distOfSeeds, 
                indexOfDistOfSeeds, alreadyChosenValues);

            seed[n] = selectedSeedValue;
            
            alreadyChosenValues.add(Integer.valueOf(selectedSeedValue));

            log.fine(String.format("choose seed %d) %d", n, 
                seed[n]));
            
        }
        
        Arrays.sort(seed);

        return seed;
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

        if (lastImgSeedIndexes == null) {
            lastImgSeedIndexes = imgSeedIndexes;
            return false;
        }
        
        boolean hasChanged = false;
        for (int i = 0; i < img.getNPixels(); ++i) {
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
            if (diff < 0 ) {
                allAreBelowCriticalLimit = false;
            }
            
            seedVariances[i] = sumVariance[i];
        }

        lastImgSeedIndexes = imgSeedIndexes;
        
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
        
        lastImgSeedIndexes = imgSeedIndexes;

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
        int[] indexOfDistOfSeeds, Set<Integer> alreadyChosenValues) {
        
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
        
        // make the cumulative array.
        long nDistDistr = 0;
        long[] nDistrC = new long[distOfSeeds.length];
        for (int i = 0; i < distOfSeeds.length; i++) {            
            int nValues = distOfSeeds[i];
            // value should be present nValues number of times
            nDistDistr += nValues;
            nDistrC[i] = nDistDistr;
        }
        
        if (nDistDistr < 1) {
            throw new IllegalStateException("distOfSeeds is in error: " + 
                Arrays.toString(distOfSeeds));
        }
        
        int chosenValue = -1;
                
        while ((chosenValue == -1) || 
            alreadyChosenValues.contains(chosenValue)) {
            
            long chosen = sr.nextLong(nDistDistr);

            int chosenIdx = findChosen(nDistrC, chosen);
            int pixIdx = indexOfDistOfSeeds[chosenIdx];
            chosenValue = img.getValue(pixIdx);
        }
        
        return chosenValue;
    }
    
    public float[] getStandardDeviationsFromCenters() {
        return this.seedVariances;
    }

    public int[] getCenters() {
        return this.center;
    }
    
    public int[] getImgPixelSeedIndexes() {
        return lastImgSeedIndexes;
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

    private void populateDistanceArrays(GreyscaleImage img, 
        Set<Integer> seeds, int[] distOfSeeds, 
        int[] indexOfDistOfSeeds) {
                 
        for (int pixIdx = 0; pixIdx < img.getNPixels(); pixIdx++) {
            
            int v = img.getValue(pixIdx);
            
            int seedValue = findNearestNeighborBruteForce(seeds, v);
                        
            int diff = v - seedValue;

            distOfSeeds[pixIdx] = (diff * diff);

            indexOfDistOfSeeds[pixIdx] = pixIdx;
        }
    }

    private int findNearestNeighborBruteForce(Set<Integer> seeds, int v) {
        
        int minDistSq = Integer.MAX_VALUE;
        int minDistV = -1;
        
        for (Integer seed : seeds) {
            int diff = seed.intValue() - v;
            int distSq = diff * diff;
            if (distSq < minDistSq) {
                minDistSq = distSq;
                minDistV = seed;
            }
        }
        
        return minDistV;
    }

    private int findChosen(long[] nDistrC, long chosen) {
        
        // find bin where cosen is found.
        // the next bin is too high in value
        
        int idx = Arrays.binarySearch(nDistrC, chosen);
        
        //nDistrC is ordered by increasing value, but there may be more than
        // one sequential item with same value (if a point is a seed, for
        // example, it's distance is 0, so cumulative sum is same).
        // so need to search before found bin also to find earliest bin.
        
        // if it's negative, (-(insertion point) - 1)
        if (idx < 0) {
            // idx = -*idx2 - 1
            idx = -1*(idx + 1);
        }
        if (idx > (nDistrC.length - 1)) {
            idx = nDistrC.length - 1;
        }
        
        long v = nDistrC[idx];
        for (int i = (idx - 1); idx > -1; --idx) {
            if (nDistrC[i] == v) {
                idx = i;
            } else {
                break;
            }
        }
        
        return idx;
    }
    
}
