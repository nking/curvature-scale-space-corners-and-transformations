package algorithms.compGeometry.clustering;

import algorithms.QuickSort;
import algorithms.imageProcessing.Image;
import algorithms.search.KDTree;
import algorithms.search.KDTreeNode;
import algorithms.util.PairInt;
import algorithms.util.TrioInt;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
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
 * Note that this color kmeans internally uses chromaticity, that is 
 * r/(r+g+b) and g/(r+g+b) and multiplies both vectors by a factor 
 * integerFactor = 255 to
 * keep the values in integers.
 * 
 * Useful reading:
 * http://en.wikipedia.org/wiki/K-means_clustering
 * 
 * @author nichole
 */
public class KMeansPlusPlusColor {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    protected int integerFactor = 255;
    
    /**
     * final solution for centers of groups (== seed centers)
     * in [0][index where 0 is for rPrime and [1][index] where 1 is for gPrime
     */
    protected int[][] center = null;
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
    
    private final ThreadLocalRandom sr;
    
    public KMeansPlusPlusColor() {
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
    public void computeMeans(int k, Image img) throws IOException, 
        NoSuchAlgorithmException {
         
        init(k);
        
        // starter seeds, sorted by increasing value
        // [r' or g'][index]
        int[][] seeds = createStartSeeds(img);
        
// TODO: error in crete start seeds        
        
        int[] imgSeedIndexes = null;

        boolean hasConverged = false;

        int nIter2 = 0;
        
        while (!hasConverged && (nIter < nMaxIter) ) {

            imgSeedIndexes = binPoints(img, seeds);

            seeds = calculateMeanOfSeedPoints(img, imgSeedIndexes);

            nIter2++;
            
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
    private int[][] createStartSeeds(Image img) throws 
        NoSuchAlgorithmException {

        // [r/(r+g+b), g/(r+g+b)][index]
        int[][] seed = new int[2][nSeeds];
        for (int i = 0; i < 2; ++i) {
            seed[i] = new int[nSeeds];
        }
        
        int[] indexes = new int[nSeeds];

        int index = getModeIdx(img);

        int nSeedsChosen = 0;
        
        int rgbSum = img.getR(index) + img.getG(index) + img.getB(index);
        if (rgbSum == 0) {
            seed[0][nSeedsChosen] = 0;
            seed[1][nSeedsChosen] = 0;
        } else {
            seed[0][nSeedsChosen] = 
                Math.round(integerFactor * (float)img.getR(index)/(float)rgbSum);
            seed[1][nSeedsChosen] = 
                Math.round(integerFactor * (float)img.getG(index)/(float)rgbSum);
        }
   
        indexes[nSeedsChosen] = index;
        
        log.fine(String.format("choose seed %d) r'=%d  g'=%d", nSeedsChosen, 
            seed[0][nSeedsChosen], seed[1][nSeedsChosen]));
        
        nSeedsChosen++;
        
        for (int n = 1; n < nSeeds; n++) {

            //TODO: may need to use 2-dimensional distance, that is dist from r' and g'
            int[] distOfSeeds = new int[img.getNPixels()];
            int[] indexOfDistOfSeeds = new int[img.getNPixels()];

            int minAllDist = Integer.MAX_VALUE;

            for (int pixIdx = 0; pixIdx < img.getNPixels(); pixIdx++) {
                
                if (contains(indexes, nSeedsChosen, pixIdx)) {
                    continue;
                }
                
                int rPt = img.getR(pixIdx);
                int gPt = img.getG(pixIdx);
                int bPt = img.getB(pixIdx);
                int sumPt = (rPt + gPt + bPt);
                int rPrimePt, gPrimePt;
                if (sumPt == 0) {
                    rPrimePt = 0;
                    gPrimePt = 0;
                } else {
                    rPrimePt = Math.round(integerFactor * (float)rPt/(float)sumPt);
                    gPrimePt = Math.round(integerFactor * (float)gPt/(float)sumPt);
                }
                
                int minDist = Integer.MAX_VALUE;

                for (int seedIndex = 0; seedIndex < nSeedsChosen; seedIndex++) {
                    
                    int diffR = rPrimePt - seed[0][seedIndex];
                    int diffG = gPrimePt - seed[1][seedIndex];
                    
                    int dist = (int)Math.round(Math.sqrt(
                        diffR * diffR + diffG * diffG));
                    
                    if (dist < minDist) {
                        minDist = dist;
                    }
                }
                                
                distOfSeeds[pixIdx] = minDist;
                indexOfDistOfSeeds[pixIdx] = pixIdx;

                if (minDist < minAllDist) {
                    minAllDist = minDist;
                }
            }
            
            index = chooseRandomlyFromNumbersPresentByProbability(distOfSeeds, 
                indexOfDistOfSeeds, indexes, nSeedsChosen);

            int r2 = img.getR(index);
            int g2 = img.getG(index);
            int b2 = img.getB(index);
            int rgbSum2 = r2 + g2 + b2;            
            if (rgbSum2 == 0) {
                seed[0][nSeedsChosen] = 0;
                seed[1][nSeedsChosen] = 0;
            } else {
                seed[0][nSeedsChosen] = Math.round(integerFactor * (float)r2/(float)rgbSum2);
                seed[1][nSeedsChosen] = Math.round(integerFactor * (float)g2/(float)rgbSum2);
            }
            indexes[nSeedsChosen] = index;

            log.fine(String.format("choose seed %d) r'=%d  g'=%d", nSeedsChosen, 
                seed[0][nSeedsChosen], seed[1][nSeedsChosen]));
            
            nSeedsChosen++;
        }
        
        QuickSort.sortByDimension1Then2(seed);
    
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
    protected int[][] calculateMeanOfSeedPoints(final Image img, 
        final int[] imgSeedIndexes) {

        int[][] sum = new int[2][nSeeds];
        for (int i = 0; i < 2; ++i) {
            sum[i] = new int[nSeeds];
        }
        int[] nSum = new int[nSeeds];

        for (int pixIdx = 0; pixIdx < img.getNPixels(); pixIdx++) {

            int seedIndex = imgSeedIndexes[pixIdx];

            int r = img.getR(pixIdx);
            int g = img.getG(pixIdx);
            int b = img.getB(pixIdx);
            int rgbSum = r + g + b;
            if (rgbSum == 0) {
                sum[0][seedIndex] = 0;
                sum[1][seedIndex] = 0;
            } else {
                sum[0][seedIndex] = Math.round(integerFactor * (float)r/(float)rgbSum);
                sum[1][seedIndex] = Math.round(integerFactor * (float)g/(float)rgbSum);
            }
            
            nSum[seedIndex]++;
        }

        for (int i = 0; i < nSeeds; i++) {
            
            if (nSum[i] == 0) {
                return null;
            } else {
                sum[0][i] /= nSum[i];
                sum[1][i] /= nSum[i];
            }

            log.fine(String.format("seed mean = %d) r'=%d g'=%d number of points=%d", 
                i, sum[0][i], sum[1][i], nSum[i]));
            
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
    protected boolean calculateVarianceFromSeedCenters(final Image img,
        int[][] seed, int[] imgSeedIndexes) {

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
            
            int r = img.getR(xyIndex);
            int g = img.getG(xyIndex);
            int b = img.getB(xyIndex);
            int rgbSum = r + g + b;
            int rPrime, gPrime;
            if (rgbSum == 0) {
                rPrime = 0;
                gPrime = 0;
            } else {
                rPrime = Math.round(integerFactor * (float)r/(float)rgbSum);
                gPrime = Math.round(integerFactor * (float)g/(float)rgbSum);
            }
            
            int diffR = rPrime - seed[0][seedIndex];
            int diffG = gPrime - seed[1][seedIndex];
            
            int d = (int)Math.round(Math.sqrt(diffR * diffR + diffG * diffG));
            
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
    
    protected void calculateFinalStats(final Image img, 
        final int[] imgSeedIndexes) {

        float[] sumStDev = new float[nSeeds];
        int[] nSumStDev = new int[nSeeds];

        for (int xyIndex = 0; xyIndex < img.getNPixels(); xyIndex++) {

            int seedIndex = imgSeedIndexes[xyIndex];

            int r = img.getR(xyIndex);
            int g = img.getG(xyIndex);
            int b = img.getB(xyIndex);
            int rgbSum = r + g + b;
            int rPrime, gPrime;
            if (rgbSum == 0) {
                rPrime = 0;
                gPrime = 0;
            } else {
                rPrime = Math.round(integerFactor * (float)r/(float)rgbSum);
                gPrime = Math.round(integerFactor * (float)g/(float)rgbSum);
            }

            int diffR = rPrime - center[0][seedIndex];
            int diffG = gPrime - center[1][seedIndex];
            
            int d = (int)Math.round(Math.sqrt(diffR * diffR + diffG * diffG));
                        
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

            log.fine(String.format("seed %d) r'=%d g'=%d stDev=%.2f number of points=%d", 
                i, center[0][i], center[1][i], seedVariances[i], nSumStDev[i]));
            
        }

        numberOfPointsPerSeedCell = nSumStDev;
    }
    
     /**
     *
     * @param img
     * @param seed array of pixel intensities of voronoi-like seeds
     * @throws IOException
     */
    protected int[] binPoints(final Image img, int[][] seed) throws IOException {

        if (img == null) {
            throw new IllegalArgumentException("img cannot be null");
        }
        if (seed == null) {
            throw new IllegalArgumentException("seed cannot be null");
        }
        
        Map<PairInt, Integer> pointIndexMap = new HashMap<PairInt, Integer>();
        
        int nc = seed[0].length;
        int[] xc = new int[nc];
        int[] yc = new int[nc];
        for (int i = 0; i < nc; ++i) {
            xc[i] = seed[0][i];
            yc[i] = seed[1][i];
            pointIndexMap.put(new PairInt(xc[i], yc[i]), Integer.valueOf(i));
        }
        
        KDTree kdTree = new KDTree(xc, yc);
        
        int[] seedNumber = new int[img.getNPixels()];
        
        for (int pixIdx = 0; pixIdx < img.getNPixels(); pixIdx++) {

            int r = img.getR(pixIdx);
            int g = img.getG(pixIdx);
            int b = img.getB(pixIdx);
            int rgbSum = r + g + b;
            int rPrime, gPrime;
            if (rgbSum == 0) {
                rPrime = 0;
                gPrime = 0;
            } else {
                rPrime = Math.round(integerFactor * (float)r/(float)rgbSum);
                gPrime = Math.round(integerFactor * (float)g/(float)rgbSum);
            }
            
            KDTreeNode node = kdTree.findNearestNeighbor(rPrime, gPrime);
            
            Integer index = pointIndexMap.get(new PairInt(node.getX(), node.getY()));
            
            assert(index != null);
            
            seedNumber[pixIdx] = index.intValue();
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
    int chooseRandomlyFromNumbersPresentByProbability(int[] distOfSeeds, 
        int[] indexOfDistOfSeeds, int[] indexesAlreadyChosen, int nIndexesAlreadyChosen) {
        
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
            throw new IllegalStateException("distOfSeeds is in error:]n" + 
                Arrays.toString(distOfSeeds));
        }
                
        while ((chosenIndex == -1) || 
            contains(indexesAlreadyChosen, nIndexesAlreadyChosen, chosenIndex)) {
            
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
    
    public int getChromaFactor() {
        return integerFactor;
    }

    /**
     * dimension 0 is rPrime values where rPrime is r/(r + g + b)
     * and dimension 1 is gPrime values where gPrime is g/(r + g + b)
     * and both values have been multiplied by 
     * integerFactor = 255 to keep them in integer range.
     * @return 
     */
    public int[][] getCenters() {
        return this.center;
    }

    public int[] getNumberOfPointsPerSeedCell() {
        return numberOfPointsPerSeedCell;
    }

    private int getModeIdx(Image img) {
        
        Map<TrioInt, Integer> counts = new HashMap<TrioInt, Integer>();
        int maxCounts = 0;
        int maxCountsIdx = -1;
        
        for (int idx = 0; idx < img.getNPixels(); idx++) {
            int r = img.getR(idx);
            int g = img.getG(idx);
            int b = img.getB(idx);
            
            TrioInt trio = new TrioInt(r, g, b);
            
            Integer value = counts.get(trio);
            int f = (value == null) ? 1 : (value.intValue() + 1);
            counts.put(trio, Integer.valueOf(f));
            if (f > maxCounts) {
                maxCounts = f;
                maxCountsIdx = idx;
            }
        }
        
        return maxCountsIdx;
    }

}
