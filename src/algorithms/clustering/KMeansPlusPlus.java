package algorithms.clustering;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.search.NearestNeighbor1D;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
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
 * runtime complexity is O(N*log_2(N))
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
     * assignments of pixels to center bins
     */
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
    
    protected int minValue = -1;
    protected int maxValue = -1;
    
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
        
        this.maxValue = Integer.MIN_VALUE;
        this.minValue = Integer.MAX_VALUE;
        
        int n = img.getNPixels();
        int[] values = new int[n];
        for (int i = 0; i < n; ++i) {
            values[i] = img.getValue(i);
            if (values[i] > maxValue) {
                maxValue = values[i];
            }
            if (values[i] < minValue) {
                minValue = values[i];
            }
        }
        
        computeMeans(k, values);
    }
    
    /**
     * note that an internal method binPoints makes an assumption that the
     * minimum values in an img pixel is 0 and a maximum is 255.
     * @param k
     * @param values
     * @throws IOException
     * @throws NoSuchAlgorithmException 
     */
    public void computeMeans(int k, int[] values) throws IOException, 
        NoSuchAlgorithmException {
         
        init(k);
        
        // starter seeds, sorted by increasing value
        int[] seeds = createStartSeeds(values);
        
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
    private int[] createStartSeeds(int[] values) throws 
        NoSuchAlgorithmException {
                
        int[] seed = new int[nSeeds];
        
        int index = getModeIdx(values);

        seed[0] = values[index];
        
        log.fine(String.format("choose seed %d) %d", 0, seed[0]));
        
        TIntSet alreadyChosenValues = new TIntHashSet();
        
        alreadyChosenValues.add(seed[0]);

        int nV = values.length;
        
        for (int n = 1; n < nSeeds; n++) {

            int[] distOfSeeds = new int[nV];
            int[] indexOfDistOfSeeds = new int[nV];

            populateDistanceArrays(values, alreadyChosenValues, 
                distOfSeeds, indexOfDistOfSeeds);
 
            int selectedSeedValue = 
                chooseRandomlyFromNumbersPresentByProbability(values,
                distOfSeeds, 
                indexOfDistOfSeeds, alreadyChosenValues);

            seed[n] = selectedSeedValue;
            
            alreadyChosenValues.add(selectedSeedValue);

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
     * @param values image values
     * @param imgSeedIndexes
     * @return
     */
    protected int[] calculateMeanOfSeedPoints(final int[] values, 
        final int[] imgSeedIndexes) {

        int[] sum = new int[nSeeds];
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
     * @param values
     * @param seed
     * @param imgSeedIndexes
     * @return
     */
    protected boolean calculateVarianceFromSeedCenters(final int[] values,
        int[] seed, int[] imgSeedIndexes) {

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

            int d = values[xyIndex] - seed[seedIndex];
            
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
    
    protected void calculateFinalStats(final int[] values, 
        final int[] imgSeedIndexes) {

        float[] sumVar = new float[nSeeds];
        int[] nSumVar = new int[nSeeds];

        for (int xyIndex = 0; xyIndex < values.length; xyIndex++) {

            int seedIndex = imgSeedIndexes[xyIndex];

            int d = values[xyIndex] - center[seedIndex];
            
            sumVar[seedIndex] += (d * d);
            
            nSumVar[seedIndex]++;
        }

        for (int i = 0; i < sumVar.length; i++) {

            if ((nSumVar[i] - 1.f) > 0) {
                // divide by N-1 rather because mean was calc'ed from the data
                sumVar[i] = (sumVar[i]/((float)nSumVar[i] - 1.f));
            } else {
                sumVar[i] = 0;
            }
            seedVariances[i] = sumVar[i];

            log.fine(String.format("seed %d) %d stDev=%.2f number of points=%d", 
                i, center[i], seedVariances[i], nSumVar[i]));
            
        }
        
        lastImgSeedIndexes = imgSeedIndexes;

        numberOfPointsPerSeedCell = nSumVar;
    }
    
     /**
     *
     * @param values
     * @param seed array of pixel intensities of voronoi-like seeds
     * @throws IOException
     */
    protected int[] binPoints(final int[] values,
        int[] seed) throws IOException {
        
        int nc = seed.length;
        
        //reverse mapping of seed and index value.
        // key = seed value, value = index of seed array
        TIntIntMap seedValueIndexMap = new TIntIntHashMap(); 
        
        NearestNeighbor1D nn1d = new NearestNeighbor1D(maxValue + 1);
        
        for (int i = 0; i < seed.length; ++i) {
            int s = seed[i];
            nn1d.insert(s);
            seedValueIndexMap.put(s, i);
        }
        
        TIntIterator iter;
             
        int[] seedNumber = new int[values.length];
        
        for (int pixIdx = 0; pixIdx < values.length; pixIdx++) {

            int v = values[pixIdx];
            
            TIntSet nearest = nn1d.findClosest(v);
            assert(nearest != null && !nearest.isEmpty());
            
            iter = nearest.iterator();
            int seedValue = iter.next();
           
            assert(seedValueIndexMap.containsKey(seedValue));
            
            int seedIdx = seedValueIndexMap.get(seedValue);
                                    
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
    int chooseRandomlyFromNumbersPresentByProbability(
        int[] values,
        int[] distOfSeeds, 
        int[] indexOfDistOfSeeds, TIntSet alreadyChosenValues) {
        
        // we want to choose randomly from the indexes based upon probabilities 
        // that scale by distance
        // so create an array that represents by number, the probability of a 
        //  value.  
        // for example, distOfSeeds={2,3,4}
        //  we'd have 
        //  distIndexDistr={0,0,1,1,1,2,2,2,2}
        // and then randomly choose an index from that.
        //
        // after choosing a number from between 0 and
        //   the length of the suggested array distIndexDistr
        //   need to find the value in that index.
        // making a lookup array of the array distIndexDistr then is
        // a cumulative array of distOfSeeds
        //           distIndex 0 starts at index 2 in distIndexDistr
        //                     1                 5
        //                     2                 9
        // then a binary search of the later returns the index

        // value = end index of ramps in the replicated distances
        long[] distOfSeedsC = new long[distOfSeeds.length];
               
        long nDistDistr = 0;
        for (int i = 0; i < distOfSeeds.length; i++) {            
            int nValues = distOfSeeds[i];
            // value should be present nValues number of times
            nDistDistr += nValues;
            distOfSeedsC[i] = nDistDistr;
        }
        
        if (nDistDistr < 1) {
            throw new IllegalStateException("distOfSeeds is in error: " + 
                Arrays.toString(distOfSeeds));
        }
        
        int chosenValue = -1;
                
        while ((chosenValue == -1) || 
            alreadyChosenValues.contains(chosenValue)) {
            
            long chosen = sr.nextLong(nDistDistr);

            int chosenIdx = findChosen(distOfSeedsC, chosen);
            int distOfSeedsIdx = indexOfDistOfSeeds[chosenIdx];
            chosenValue = values[distOfSeedsIdx];
        }
        
        return chosenValue;
    }
    
    public float[] getSeedVariances() {
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

    private int getModeIdx(int[] values) {
        
        if (imgModeIdx > -1) {
            return imgModeIdx;
        }
        
        TIntIntMap countMap = new TIntIntHashMap();
        int maxCounts = 0;
        int maxCountsIdx = -1;
        
        for (int idx = 0; idx < values.length; idx++) {
            int v = values[idx];
            int f;
            if (countMap.containsKey(v)) {
                f = countMap.get(v) + 1;
            } else {
                f = 1;
            }
            countMap.put(v, f);
            if (f > maxCounts) {
                maxCounts = f;
                maxCountsIdx = idx;
            }
        }
        
        imgModeIdx = maxCountsIdx;
        
        return maxCountsIdx;
    }

    private void populateDistanceArrays(int[] values, 
        TIntSet seeds, int[] distOfSeeds, 
        int[] indexOfDistOfSeeds) {
       
        NearestNeighbor1D nn1d = new NearestNeighbor1D(maxValue + 1);
        
        TIntIterator iter = seeds.iterator();
        while (iter.hasNext()) {
            int s = iter.next();
            nn1d.insert(s);
        }
        
        for (int pixIdx = 0; pixIdx < values.length; pixIdx++) {
            
            int v = values[pixIdx];
            
            TIntSet nearest = nn1d.findClosest(v);
            
            assert(nearest != null && !nearest.isEmpty());
            
            int seedValue = nearest.iterator().next();
                        
            int diff = v - seedValue;

            distOfSeeds[pixIdx] = (diff * diff);

            indexOfDistOfSeeds[pixIdx] = pixIdx;
        }
    }

    private int findChosen(long[] nDistrC, long chosen) {
        
        // find bin where chosen is found.
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
        for (int i = (idx - 1); i > -1; --i) {
            if (nDistrC[i] == v) {
                idx = i;
            } else {
                break;
            }
        }
        
        return idx;
    }

}
