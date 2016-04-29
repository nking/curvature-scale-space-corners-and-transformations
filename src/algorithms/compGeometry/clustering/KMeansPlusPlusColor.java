package algorithms.compGeometry.clustering;

import algorithms.QuickSort;
import algorithms.imageProcessing.Image;
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
    
    protected int integerFactor = 51;//255;
    
    /**
     * final solution for centers of groups (== seed centers)
     * in [0][index where 0 is for rPrime and [1][index] where 1 is for gPrime
     */
    protected int[][] center = null;
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
    protected final static int nMaxIter2 = 200;
    protected int nIter = 0;
    
    protected int imgModeIdx = -1;
    
    // key = r,g chromaticity and value = number of pixels with key
    protected Map<PairInt, Integer> countsMap = null;
    
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
         
        int[] imgSeedIndexes = null;

        boolean hasConverged = false;
        
        int nIter2 = 0;
        
        while (!hasConverged && (nIter < nMaxIter) && (nIter2 < nMaxIter2)) {
            
            imgSeedIndexes = binPoints(img, seeds);
            
            Set<Integer> zeroSumIndexes = new HashSet<Integer>();
            
            int[][] tmpSeeds = calculateMeanOfSeedPoints(img, imgSeedIndexes, 
                zeroSumIndexes);
            
            if (tmpSeeds == null) {
                nIter2++;
                nIter = 0;
                lastImgSeedIndexes = null;
                if ((nIter2 % nSeeds) == 0) {
                    seeds = createStartSeeds(img);
                } else {
                    seeds = replaceZeroSumIndexes(img, seeds, zeroSumIndexes);
                }
                continue;
            }
            
            seeds = tmpSeeds;

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
    private int[][] createStartSeeds(Image img) throws NoSuchAlgorithmException {

        // [r/(r+g+b), g/(r+g+b)][index]
        int[][] seed = new int[2][nSeeds];
        for (int i = 0; i < 2; ++i) {
            seed[i] = new int[nSeeds];
        }
        
        int[] indexes = new int[nSeeds];
        
        PairInt[] indexColors = new PairInt[nSeeds];
        
        int selectedPixIdx = getModeIdx(img);

        int nSeedsChosen = 0;
        PairInt clr = calculateRGChromaticity(img, selectedPixIdx);
        seed[0][nSeedsChosen] = clr.getX();
        seed[1][nSeedsChosen] = clr.getY();
   
        indexes[nSeedsChosen] = selectedPixIdx;
        indexColors[nSeedsChosen] = clr;
        
        log.fine(String.format("choose seed %d) r'=%d  g'=%d", nSeedsChosen, 
            seed[0][nSeedsChosen], seed[1][nSeedsChosen]));
        
        nSeedsChosen++;
      
        for (int n = 1; n < nSeeds; n++) {
            
            int[] distOfSeeds = new int[img.getNPixels()];
            int[] indexOfDistOfSeeds = new int[img.getNPixels()];

            int[] xs = new int[nSeedsChosen];
            int[] ys = new int[nSeedsChosen];
            for (int seedIndex = 0; seedIndex < nSeedsChosen; seedIndex++) {
                xs[seedIndex] = seed[0][seedIndex];
                ys[seedIndex] = seed[1][seedIndex];
            }
            
            populateDistanceArrays(img, xs, ys, distOfSeeds, indexOfDistOfSeeds);
            
            selectedPixIdx = chooseRandomlyFromNumbersPresentByProbability(img,
                distOfSeeds, 
                indexOfDistOfSeeds, indexes, indexColors, nSeedsChosen);
            
            PairInt selectedPixIdxColors = calculateRGChromaticity(img, selectedPixIdx);
     
            while (contains(indexes, nSeedsChosen, selectedPixIdx, 
                indexColors, selectedPixIdxColors)) {
                
                selectedPixIdx = chooseRandomlyFromNumbersPresentByProbability(img,
                    distOfSeeds, 
                    indexOfDistOfSeeds, indexes, indexColors, nSeedsChosen);
                
                selectedPixIdxColors = calculateRGChromaticity(img, selectedPixIdx);
            }
            
            indexes[nSeedsChosen] = selectedPixIdx;
            indexColors[nSeedsChosen] = selectedPixIdxColors;

            log.fine(String.format("choose seed %d) r'=%d  g'=%d", nSeedsChosen, 
                seed[0][nSeedsChosen], seed[1][nSeedsChosen]));
            
            nSeedsChosen++;
        }
        
        QuickSort.sortByDimension1Then2(seed);
    
        return seed;
    }
    
    protected boolean contains(int[] array, int nArray, int value,
        PairInt[] arrayColors, PairInt valueColors) {
        
        for (int i = 0; i < nArray; i++) {
            if ((array[i] == value) || 
                ((valueColors != null) && arrayColors[i].equals(valueColors))) {
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
     * @param outputZeroSums output list holding seeds indexes where there
     * were no points in their bins.
     * @return
     */
    protected int[][] calculateMeanOfSeedPoints(final Image img, 
        final int[] imgSeedIndexes, Set<Integer> outputZeroSums) {

        int[][] sum = new int[2][nSeeds];
        for (int i = 0; i < 2; ++i) {
            sum[i] = new int[nSeeds];
        }
        int[] nSum = new int[nSeeds];

        for (int pixIdx = 0; pixIdx < img.getNPixels(); pixIdx++) {

            int seedIndex = imgSeedIndexes[pixIdx];

            PairInt clr = calculateRGChromaticity(img, pixIdx);
            sum[0][seedIndex] = clr.getX();
            sum[1][seedIndex] = clr.getY();
        
            nSum[seedIndex]++;
        }
        
        for (int i = 0; i < nSeeds; i++) {
            
            if (nSum[i] == 0) {
                outputZeroSums.add(Integer.valueOf(i));
            } else {
                sum[0][i] /= nSum[i];
                sum[1][i] /= nSum[i];
            }

            log.fine(String.format("seed mean = %d) r'=%d g'=%d number of points=%d", 
                i, sum[0][i], sum[1][i], nSum[i]));
        }
        
        if (outputZeroSums.size() > 0) {
            return null;
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
            
            PairInt clr = calculateRGChromaticity(img, xyIndex);
            int rPrime = clr.getX();
            int gPrime = clr.getY();
            
            int diffR = rPrime - seed[0][seedIndex];
            int diffG = gPrime - seed[1][seedIndex];
            
            int d = diffR * diffR + diffG * diffG;
            
            sumVariance[seedIndex] += d;
            
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

        return allAreBelowCriticalLimit;
    }
    
    protected void calculateFinalStats(final Image img, 
        final int[] imgSeedIndexes) {

        float[] sumStDev = new float[nSeeds];
        int[] nSumStDev = new int[nSeeds];

        for (int xyIndex = 0; xyIndex < img.getNPixels(); xyIndex++) {

            int seedIndex = imgSeedIndexes[xyIndex];

            PairInt clr = calculateRGChromaticity(img, xyIndex);
            int rPrime = clr.getX();
            int gPrime = clr.getY();

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
        
        KDTree kdTree = null;
        if (nc > 2) {
            kdTree = new KDTree(xc, yc);
        }
                    
        int[] seedNumber = new int[img.getNPixels()];
        
        for (int pixIdx = 0; pixIdx < img.getNPixels(); pixIdx++) {

            PairInt clr = calculateRGChromaticity(img, pixIdx);
            int rPrime = clr.getX();
            int gPrime = clr.getY();
            
            KDTreeNode node;
            if (kdTree != null) {
                node = kdTree.findNearestNeighbor(rPrime, gPrime);
            } else {
                node = findNearestNeighborBruteForce(xc, yc, rPrime, gPrime);
            }
            
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
     * @param indexColors the colors of the pixels at indexOfDistOfSeeds
     * @param sr
     * @param indexesAlreadyChosen
     * @param nIndexesAlreadyChosen
     * @return 
     */
    int chooseRandomlyFromNumbersPresentByProbability(Image img,
        int[] distOfSeeds, 
        int[] indexOfDistOfSeeds, int[] indexesAlreadyChosen, 
        PairInt[] indexColors, int nIndexesAlreadyChosen) {
        
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
        
        PairInt chosenIndexColors = (chosenIndex == -1) ? null :
            calculateRGChromaticity(img, chosenIndex);
                
        while ((chosenIndex == -1) || 
            contains(indexesAlreadyChosen, nIndexesAlreadyChosen, chosenIndex,
                indexColors, chosenIndexColors)) {
                        
            long chosen = sr.nextLong(nDistDistr);

            // walk thru same iteration to obtain the chosen index
            long n = 0;
            for (int i = 0; i < distOfSeeds.length; i++) {            
                int nValues = distOfSeeds[i];
                // value should be present nValues number of times
                
                if ((chosen >= n) && (chosen < (n + nValues))) {
                    chosenIndex = indexOfDistOfSeeds[i];
                    chosenIndexColors = calculateRGChromaticity(img, chosenIndex);
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
        
        if (imgModeIdx > -1) {
            return imgModeIdx;
        }
        
        Map<PairInt, Integer> counts = new HashMap<PairInt, Integer>();
        int maxCounts = 0;
        int maxCountsIdx = -1;
        
        for (int idx = 0; idx < img.getNPixels(); idx++) {
            
            PairInt clr = calculateRGChromaticity(img, idx);
            
            Integer value = counts.get(clr);
            int f = (value == null) ? 1 : (value.intValue() + 1);
            counts.put(clr, Integer.valueOf(f));
            if (f > maxCounts) {
                maxCounts = f;
                maxCountsIdx = idx;
            }
        }
        
        if (counts.size() < nSeeds) {
            throw new IllegalStateException("there are only " + counts.size() +
                " unique chromaticity values in the image, so cannot form " 
                + nSeeds + " clusters.");
        }
        
        imgModeIdx = maxCountsIdx;
        
        countsMap = counts;
        
        // can return any of the pixels with mode colors since he algorithm
        // uses distance from colors rather than location in x,y coordinates.
        
        return maxCountsIdx;
    }

    private KDTreeNode findNearestNeighborBruteForce(int[] xs, int[] ys, 
        int x, int y) {
        
        int minDistSq = Integer.MAX_VALUE;
        int minDistIdx = -1;
        
        for (int i = 0; i < xs.length; ++i) {
            int diffX = xs[i] - x;
            int diffY = ys[i] - y;
            int distSq = diffX * diffX + diffY * diffY;
            if (distSq < minDistSq) {
                minDistSq = distSq;
                minDistIdx = i;
            }
        }
        
        KDTreeNode node = new KDTreeNode();
        node.setX(xs[minDistIdx]);
        node.setY(ys[minDistIdx]);
        
        return node;
    }

    private int[][] replaceZeroSumIndexes(Image img, int[][] seeds, 
        Set<Integer> zeroSumIndexes) {
        
        if (zeroSumIndexes.size() == 0) {
            throw new IllegalStateException("zeroSumIndexes should not be"
                + " empty");
        }
        
        assert(seeds[0].length == nSeeds);
        
        Set<Integer> current = new HashSet<Integer>();
        for (int i = 0; i < nSeeds; ++i) {
            Integer index = Integer.valueOf(i);
            if (!zeroSumIndexes.contains(index)) {
                current.add(index);
            }
        }
        
  if (current.size() == 0) {
      int z = 1;
  }
  
        assert(current.size() == (nSeeds - zeroSumIndexes.size()));
       
        for (Integer seedIndex : zeroSumIndexes) {
            
            int[] distOfSeeds = new int[img.getNPixels()];
            int[] indexOfDistOfSeeds = new int[img.getNPixels()];

            PairInt[] currentSeedColors = new PairInt[current.size()];
            int[] currentSeedIndexes = new int[current.size()];
            
            int[] xs = new int[current.size()];
            int[] ys = new int[current.size()];
            int count = 0;
            for (int i = 0; i < seeds[0].length; ++i) {
                Integer index = Integer.valueOf(i);
                int idx = index.intValue();
                if (current.contains(index)) {
                    xs[count] = seeds[0][idx];
                    ys[count] = seeds[1][idx];
                    
                    currentSeedIndexes[count] = idx;
                    currentSeedColors[count] = calculateRGChromaticity(img, idx);
                        
                    count++;
                }
            }
            
            populateDistanceArrays(img, xs, ys, distOfSeeds, indexOfDistOfSeeds);
           
            int selectedPixIdx = chooseRandomlyFromNumbersPresentByProbability(
                img, distOfSeeds, indexOfDistOfSeeds, 
                currentSeedIndexes, currentSeedColors, current.size());
            
            PairInt selectedPixIdxColors = calculateRGChromaticity(img, selectedPixIdx);
            
            while (contains(currentSeedIndexes, current.size(), selectedPixIdx,
                currentSeedColors, selectedPixIdxColors)) {
                
                selectedPixIdx = chooseRandomlyFromNumbersPresentByProbability(
                    img, distOfSeeds, indexOfDistOfSeeds, 
                    currentSeedIndexes, currentSeedColors, current.size());
                
                selectedPixIdxColors = calculateRGChromaticity(img, selectedPixIdx);
            }
            
            seeds[0][seedIndex.intValue()] = selectedPixIdxColors.getX();
            seeds[1][seedIndex.intValue()] = selectedPixIdxColors.getY();
            
            current.add(seedIndex);
        }
        
        assert(current.size() == nSeeds);
        
        QuickSort.sortByDimension1Then2(seeds);
        
        return seeds;
    }

    private void populateDistanceArrays(Image img, int[] xs, int[] ys, 
        int[] distOfSeeds, int[] indexOfDistOfSeeds) {
        
        KDTree kdTree = null;
        if (xs.length > 2) {
            kdTree = new KDTree(xs, ys);
        }

        for (int pixIdx = 0; pixIdx < img.getNPixels(); pixIdx++) {
            
            PairInt clr = calculateRGChromaticity(img, pixIdx);
            int rPrime = clr.getX();
            int gPrime = clr.getY();
            
            KDTreeNode seedNode;
            if (kdTree != null) {
                seedNode = kdTree.findNearestNeighbor(rPrime, gPrime);
            } else {
                seedNode = findNearestNeighborBruteForce(xs, ys, rPrime, gPrime);
            }

            int diffR = rPrime - seedNode.getX();
            int diffG = gPrime - seedNode.getY();

            int distSq = diffR * diffR + diffG * diffG;

            distOfSeeds[pixIdx] = distSq;

            indexOfDistOfSeeds[pixIdx] = pixIdx;
        }
            
    }

    private PairInt calculateRGChromaticity(Image img, int pixIdx) {
        
        int r = img.getR(pixIdx);
        int g = img.getG(pixIdx);
        int b = img.getB(pixIdx);
        int rgbSum = r + g + b;
        int rPrime, gPrime;
        if (rgbSum == 0) {
            rPrime = 0;
            gPrime = 0;
        } else {
            rPrime = Math.round(integerFactor * (float) r / (float) rgbSum);
            gPrime = Math.round(integerFactor * (float) g / (float) rgbSum);
        }

        return new PairInt(rPrime, gPrime);
    }

}
