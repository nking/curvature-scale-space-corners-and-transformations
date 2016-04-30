package algorithms.compGeometry.clustering;

import algorithms.QuickSort;
import algorithms.imageProcessing.Image;
import algorithms.util.TrioInt;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
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
 * Note that this color kmeans internally uses r, g, b.
 * 
 * Useful reading:
 * http://en.wikipedia.org/wiki/K-means_clustering
 * 
 * @author nichole
 */
public class KMeansPlusPlusColor {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
            
    /**
     * final solution for centers of groups (== seed centers)
     * in [0][index where 0 is for r and [1][index] where 1 is for g
     * and [2][index] is for b
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
    
    protected final static int nMaxIter = 100;
    protected final static int nMaxIter2 = 200;
    protected int nIter = 0;
    
    protected int imgModeIdx = -1;
    
    protected TrioInt[] imageColorsAscending = null;
    protected int[] imageColorsAscendingCounts = null;
    
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
        // [r,g,b][index]
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
                //if ((nIter2 % nSeeds) == 0) {
                    seeds = createStartSeeds(img);
                /*} else if (nIter2 < 50) {
                    seeds = replaceZeroSumIndexes(img, seeds, zeroSumIndexes);
                } else {
                    seeds = replaceZeroSumIndexesUsingFrequencyMap(img, seeds, zeroSumIndexes);
                }*/
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
     * choose seeds sequentially by distance weighted probabilities.  Note that
     * the starter point is chosen to be the mode of the image's colors and
     * thereafter is a color randomly chosen from a distribution weighted
     * by distances from nearest seeds.
     * 
     * @param img
     * @return 
     */
    private int[][] createStartSeeds(Image img) throws NoSuchAlgorithmException {

        // [(r,g, or b)][index]
        int[][] seed = new int[3][nSeeds];
        for (int i = 0; i < 3; ++i) {
            seed[i] = new int[nSeeds];
        }
                        
        int selectedPixIdx = getModeIdx(img);

        TrioInt clr = calculateRGB(img, selectedPixIdx);
        seed[0][0] = clr.getX();
        seed[1][0] = clr.getY();
        seed[2][0] = clr.getZ();
           
        log.fine(String.format("choose seed %d) r=%d  g=%d  b=%d", 0, 
            seed[0][0], seed[1][0], seed[2][0]));
                
        int nPix = img.getNPixels();
      
        for (int n = 1; n < nSeeds; n++) {
            
            Set<TrioInt> colorsAlreadyChosen = new HashSet<TrioInt>();
            for (int seedIndex = 0; seedIndex < n; seedIndex++) {
                colorsAlreadyChosen.add(new TrioInt(seed[0][seedIndex], 
                    seed[1][seedIndex], seed[2][seedIndex]));
            }
            
            int[] distOfSeedsX = new int[nPix];
            int[] distOfSeedsY = new int[nPix];
            int[] distOfSeedsZ = new int[nPix];
            int[] indexOfDistOfSeedsX = new int[nPix];
            int[] indexOfDistOfSeedsY = new int[nPix];
            int[] indexOfDistOfSeedsZ = new int[nPix];
            
            populateDistanceArrays(img, colorsAlreadyChosen, 
                distOfSeedsX, indexOfDistOfSeedsX,
                distOfSeedsY, indexOfDistOfSeedsY,
                distOfSeedsZ, indexOfDistOfSeedsZ
            );
            
            // thi selects until one not already selected is found
            TrioInt selectedColor = 
                chooseRandomlyFromNumbersPresentByProbability(img,
                distOfSeedsX, indexOfDistOfSeedsX,
                distOfSeedsY, indexOfDistOfSeedsY,
                distOfSeedsZ, indexOfDistOfSeedsZ,
                colorsAlreadyChosen);
            
            seed[0][n] = selectedColor.getX();
            seed[1][n] = selectedColor.getY();
            seed[2][n] = selectedColor.getZ();
            
            log.fine(String.format("choose seed %d) r=%d g=%d b=%d", n, 
                seed[0][n], seed[1][n], seed[2][n]));            
        }
        
        QuickSort.sortByDimension1FirstSecondThird(seed);
    
        return seed;
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

        int[][] sum = new int[3][nSeeds];
        for (int i = 0; i < 3; ++i) {
            sum[i] = new int[nSeeds];
        }
        int[] nSum = new int[nSeeds];

        for (int pixIdx = 0; pixIdx < img.getNPixels(); pixIdx++) {

            int seedIndex = imgSeedIndexes[pixIdx];

            TrioInt clr = calculateRGB(img, pixIdx);
            sum[0][seedIndex] += clr.getX();
            sum[1][seedIndex] += clr.getY();
            sum[2][seedIndex] += clr.getZ();
        
            nSum[seedIndex]++;
        }
        
        for (int i = 0; i < nSeeds; i++) {
            
            if (nSum[i] == 0) {
                outputZeroSums.add(Integer.valueOf(i));
            } else {
                sum[0][i] /= nSum[i];
                sum[1][i] /= nSum[i];
                sum[2][i] /= nSum[i];
            }

            log.fine(String.format("seed mean = %d) r'=%d g'=%d b=%d number of points=%d", 
                i, sum[0][i], sum[1][i], sum[2][i], nSum[i]));
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
                        
            int diffR = img.getR(xyIndex) - seed[0][seedIndex];
            int diffG = img.getG(xyIndex) - seed[1][seedIndex];
            int diffB = img.getB(xyIndex) - seed[2][seedIndex];
            
            int d = diffR * diffR + diffG * diffG + diffB * diffB;
            
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
        
        lastImgSeedIndexes = imgSeedIndexes;

        return allAreBelowCriticalLimit;
    }
    
    protected void calculateFinalStats(final Image img, 
        final int[] imgSeedIndexes) {

        float[] sumStDev = new float[nSeeds];
        int[] nSumStDev = new int[nSeeds];

        for (int xyIndex = 0; xyIndex < img.getNPixels(); xyIndex++) {

            int seedIndex = imgSeedIndexes[xyIndex];

            int diffR = img.getR(xyIndex) - center[0][seedIndex];
            int diffG = img.getG(xyIndex) - center[1][seedIndex];
            int diffB = img.getB(xyIndex) - center[2][seedIndex];
            
            int d = diffR * diffR + diffG * diffG + diffB * diffB;
            
            sumStDev[seedIndex] += d;
            
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

            log.fine(String.format("seed %d) r=%d g=%d b=%d stDev=%.2f number of points=%d", 
                i, center[0][i], center[1][i], center[2][i], seedVariances[i], nSumStDev[i]));
            
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
    protected int[] binPoints(final Image img, int[][] seed) throws IOException {
        
        if (img == null) {
            throw new IllegalArgumentException("img cannot be null");
        }
        if (seed == null) {
            throw new IllegalArgumentException("seed cannot be null");
        }
        
        int[] seedNumber = new int[img.getNPixels()];
        
        for (int pixIdx = 0; pixIdx < img.getNPixels(); pixIdx++) {

            TrioInt clr = calculateRGB(img, pixIdx);

            int seedIdx = findNearestNeighbor(seed, clr.getX(), clr.getY(),
                clr.getZ());
     
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
     * @param indexColors the colors of the pixels at indexOfDistOfSeeds
     * @param sr
     * @param indexesAlreadyChosen
     * @param nIndexesAlreadyChosen
     * @return 
     */
    TrioInt chooseRandomlyFromNumbersPresentByProbability(Image img,
        int[] distOfSeedsX, int[] indexOfDistOfSeedsX,
        int[] distOfSeedsY, int[] indexOfDistOfSeedsY,
        int[] distOfSeedsZ, int[] indexOfDistOfSeedsZ,
        Set<TrioInt> colorsAlreadyChosen) {
                
        // we want to choose randomly from the indexes based upon probabilities 
        // that scale by distance
        // so create an array that represents by number, the probability of a 
        //  value.  for example, distOfSeeds={2,3,4}
        //  we'd have 
        //  distIndexDistr={0,0,1,1,1,2,2,2,2}
        // and then randomly choose from that.
        // here, we skip storing every value in a potentially very large array 
        // and instead,
        // find the value for the position once the position has been 
        // drawn randomly
        
        // ****TODO: this could be improved by storing the sums at large intervals *****
        long nDistDistrX = 0;
        for (int i = 0; i < distOfSeedsX.length; i++) {            
            int nValues = distOfSeedsX[i];
            // value should be present nValues number of times
            nDistDistrX += nValues;
        }
        long nDistDistrY = 0;
        for (int i = 0; i < distOfSeedsY.length; i++) {            
            int nValues = distOfSeedsY[i];
            nDistDistrY += nValues;
        }
        long nDistDistrZ = 0;
        for (int i = 0; i < distOfSeedsZ.length; i++) {            
            int nValues = distOfSeedsZ[i];
            nDistDistrZ += nValues;
        }
        
        if (nDistDistrX < 1) {
            throw new IllegalStateException("distOfSeeds is in error:]n" + 
                Arrays.toString(distOfSeedsX));
        }
        if (nDistDistrY < 1) {
            throw new IllegalStateException("distOfSeeds is in error:]n" + 
                Arrays.toString(distOfSeedsX));
        }
        if (nDistDistrZ < 1) {
            throw new IllegalStateException("distOfSeeds is in error:]n" + 
                Arrays.toString(distOfSeedsX));
        }
                
        TrioInt chosenIndexColors = null;
                
        while ((chosenIndexColors == null) || 
            colorsAlreadyChosen.contains(chosenIndexColors)) {
            
            long nDistDistr;
            int[] distOfSeeds;
            int[] indexOfDistOfSeeds;
            
            int xYZ = sr.nextInt(3);
            
            /*
            NOTE: could choose randomly from all three and chose the one
            among them that is more frequent in the frequency data stored
            from getModeIdx method.
            */
            
            switch (xYZ) {
                case 0:
                    nDistDistr = nDistDistrX;
                    distOfSeeds = distOfSeedsX;
                    indexOfDistOfSeeds = indexOfDistOfSeedsX;
                    break;
                case 1:
                    nDistDistr = nDistDistrY;
                    distOfSeeds = distOfSeedsY;
                    indexOfDistOfSeeds = indexOfDistOfSeedsY;
                    break;
                default:
                    nDistDistr = nDistDistrZ;
                    distOfSeeds = distOfSeedsZ;
                    indexOfDistOfSeeds = indexOfDistOfSeedsZ;
                    break;
            }
                        
            long chosen = sr.nextLong(nDistDistr);

            // walk thru same iteration to obtain the chosen index
            long n = 0;
            for (int i = 0; i < distOfSeeds.length; i++) {            
                int nValues = distOfSeeds[i];
                // value should be present nValues number of times
                if ((chosen >= n) && (chosen < (n + nValues))) {
                    int chosenIndex = indexOfDistOfSeeds[i];
                    chosenIndexColors = calculateRGB(img, chosenIndex);
                    break;
                }
                n += nValues;
            }
        }
        
        return chosenIndexColors;
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
    TrioInt chooseRandomlyFromNumbersPresentByProbability(Image img,
        int[] counts, TrioInt[] colors,
        Set<TrioInt> colorsAlreadyChosen, int nIndexesAlreadyChosen) {
               
        // we want to choose randomly from the indexes based upon probabilities 
        // that scale by distance
        // so create an array that represents by number, the probability of a 
        //  value.  for example, distOfSeeds={2,3,4}
        //  we'd have 
        //  distIndexDistr={0,0,1,1,1,2,2,2,2}
        // and then randomly choose from that.
        // here, we skip storing every value in a potentially very large array 
        // and instead,
        // find the value for the position once the position has been 
        // drawn randomly
        
        // ****TODO: this could be improved by storing the sums at large intervals *****
        long nCountDistr = 0;
        for (int i = 0; i < counts.length; i++) {            
            int nValues = counts[i];
            // value should be present nValues number of times
            nCountDistr += nValues;
        }
        
        if (nCountDistr < 1) {
            throw new IllegalStateException("distOfSeeds is in error:]n");
        }
       
        TrioInt chosenIndexColor = null;
        
        while ((chosenIndexColor == null) || 
            colorsAlreadyChosen.contains(chosenIndexColor)) {
            
            long chosen = sr.nextLong(nCountDistr);

            // walk thru same iteration to obtain the chosen index
            long n = 0;
            for (int i = 0; i < counts.length; i++) {            
                int nValues = counts[i];
                // value should be present nValues number of times
                if ((chosen >= n) && (chosen < (n + nValues))) {
                    chosenIndexColor = colors[i];
                    break;
                }
                n += nValues;
            }
        }
        
        return chosenIndexColor;
    }
    
    public float[] getStandardDeviationsFromCenters() {
        return this.seedVariances;
    }
    
    /**
     * accessed s [0,1,2][seedIndex] where the first dimension 0, 1 and 2 hold
     * r, g, and b, respectively.
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
        
        Map<TrioInt, Integer> counts = new HashMap<TrioInt, Integer>();
        int maxCounts = 0;
        int maxCountsIdx = -1;
        
        for (int idx = 0; idx < img.getNPixels(); idx++) {
            
            TrioInt clr = calculateRGB(img, idx);
            
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
        
        int[] colorCounts = new int[counts.size()];
        TrioInt[] colors = new TrioInt[counts.size()];
        int count = 0;
        for (Entry<TrioInt, Integer> entry : counts.entrySet()) {
            colorCounts[count] = entry.getValue();
            colors[count] = entry.getKey();
            count++;
        }
        QuickSort.sortBy1stArg(colorCounts, colors);
        
        this.imageColorsAscending = colors;
        this.imageColorsAscendingCounts = colorCounts;
        
        // can return any of the pixels with mode colors since he algorithm
        // uses distance from colors rather than location in x,y coordinates.
        
        return maxCountsIdx;
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
  
        assert(current.size() == (nSeeds - zeroSumIndexes.size()));
       
        int nPix = img.getNPixels();
        
        for (Integer seedIndex : zeroSumIndexes) {
            
            Set<TrioInt> colorsAlreadyChosen = new HashSet<TrioInt>();
            for (int i = 0; i < seeds[0].length; ++i) {
                Integer index = Integer.valueOf(i);
                int idx = index.intValue();
                if (current.contains(index)) {                    
                    colorsAlreadyChosen.add(new TrioInt(
                        seeds[0][idx], seeds[1][idx], seeds[2][idx]));
                }
            }
            
            int[] distOfSeedsX = new int[nPix];
            int[] distOfSeedsY = new int[nPix];
            int[] distOfSeedsZ = new int[nPix];
            int[] indexOfDistOfSeedsX = new int[nPix];
            int[] indexOfDistOfSeedsY = new int[nPix];
            int[] indexOfDistOfSeedsZ = new int[nPix];
            
            populateDistanceArrays(img, colorsAlreadyChosen, 
                distOfSeedsX, indexOfDistOfSeedsX,
                distOfSeedsY, indexOfDistOfSeedsY,
                distOfSeedsZ, indexOfDistOfSeedsZ
            );
            
            TrioInt selectedClr = chooseRandomlyFromNumbersPresentByProbability(img,
                distOfSeedsX, indexOfDistOfSeedsX,
                distOfSeedsY, indexOfDistOfSeedsY,
                distOfSeedsZ, indexOfDistOfSeedsZ,
                colorsAlreadyChosen);
            
            seeds[0][seedIndex.intValue()] = selectedClr.getX();
            seeds[1][seedIndex.intValue()] = selectedClr.getY();
            seeds[2][seedIndex.intValue()] = selectedClr.getZ();
            
            current.add(seedIndex);
        }
        
        assert(current.size() == nSeeds);
        
        QuickSort.sortByDimension1FirstSecondThird(seeds);
        
        return seeds;
    }

    private int[][] replaceZeroSumIndexesUsingFrequencyMap(Image img, int[][] seeds, 
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
  
        assert(current.size() == (nSeeds - zeroSumIndexes.size()));
       
        int nPix = img.getNPixels();
        
        for (Integer seedIndex : zeroSumIndexes) {
            
            Set<TrioInt> selectedColorsSet = new HashSet<TrioInt>();
            for (int i = 0; i < seeds[0].length; ++i) {
                Integer index = Integer.valueOf(i);
                int idx = index.intValue();
                if (current.contains(index)) {                    
                    selectedColorsSet.add(calculateRGB(img, idx));
                }
            }
        
            TrioInt selectedColor = chooseRandomlyFromNumbersPresentByProbability(img,
                imageColorsAscendingCounts, imageColorsAscending,
                selectedColorsSet, current.size());
            
            seeds[0][seedIndex.intValue()] = selectedColor.getX();
            seeds[1][seedIndex.intValue()] = selectedColor.getY();
            seeds[2][seedIndex.intValue()] = selectedColor.getZ();
            
            current.add(seedIndex);
        }
        
        assert(current.size() == nSeeds);
        
        QuickSort.sortByDimension1FirstSecondThird(seeds);
        
        return seeds;
    }

    private void populateDistanceArrays(Image img, Set<TrioInt> seeds,
        int[] distOfSeedsX, int[] indexOfDistOfSeedsX,
        int[] distOfSeedsY, int[] indexOfDistOfSeedsY,
        int[] distOfSeedsZ, int[] indexOfDistOfSeedsZ) {
                    
        /*        
        to randomly sample from cumulative probability distributions made
        separately for distance in x, y, and z,
        one way would be to first select randomly from 0, 1, 2 (for x-axis,
            y-axis, or z-axis),
            then choose a point randomly within that chosen distribution.
            All three cumulative distributions result in different lengths 
            and do not have the same values (pixel indexes) as a relationship.
        
                pt1       pt2
            X:  ---       --------------
            Y   -------   --------
            Z   ---       --------------
            
        */
        
        for (int pixIdx = 0; pixIdx < img.getNPixels(); pixIdx++) {
            
            TrioInt clr = calculateRGB(img, pixIdx);
            
            TrioInt seed = findNearestNeighbor(seeds, clr.getX(), 
                clr.getY(), clr.getZ());

            int diffR = clr.getX() - seed.getX();
            int diffG = clr.getY() - seed.getY();
            int diffB = clr.getZ() - seed.getZ();

            distOfSeedsX[pixIdx] = diffR * diffR;
            distOfSeedsY[pixIdx] = diffG * diffG;
            distOfSeedsZ[pixIdx] = diffB * diffB;
            
            indexOfDistOfSeedsX[pixIdx] = pixIdx;
            indexOfDistOfSeedsY[pixIdx] = pixIdx;
            indexOfDistOfSeedsZ[pixIdx] = pixIdx;
        }
            
    }

    private TrioInt calculateRGB(Image img, int pixIdx) {
        
        int r = img.getR(pixIdx);
        int g = img.getG(pixIdx);
        int b = img.getB(pixIdx);

        return new TrioInt(r, g, b);
    }

    private int findNearestNeighbor(int[][] seed, int xs, int ys, int zs) {
        
        //TODO: adapt the KD-Tree to three dimensions
        return findNearestNeighborBruteForce(seed, xs, ys, zs);
    }
    
    private int findNearestNeighborBruteForce(int[][] seed, int xs, int ys, 
        int zs) {
        
        long minDistSq = Long.MAX_VALUE;
        int minDistIdx = -1;
        
        for (int i = 0; i < seed[0].length; ++i) {
            int diffX = seed[0][i] - xs;
            int diffY = seed[1][i] - ys;
            int diffZ = seed[2][i] - zs;
            long distSq = diffX * diffX + diffY * diffY + diffZ * diffZ;
            if (distSq < minDistSq) {
                minDistSq = distSq;
                minDistIdx = i;
            }
        }
        
        return minDistIdx;
    }

    private TrioInt findNearestNeighbor(Set<TrioInt> seeds, int xs, int ys, int zs) {
        
        //TODO: adapt the KD-Tree to three dimensions
        return findNearestNeighborBruteForce(seeds, xs, ys, zs);
    }
    
    private TrioInt findNearestNeighborBruteForce(Set<TrioInt> seeds, int xs, int ys, 
        int zs) {
        
        long minDistSq = Long.MAX_VALUE;
        TrioInt minDistClr = null;
        
        for (TrioInt seed : seeds) {
            int diffX = seed.getX() - xs;
            int diffY = seed.getY() - ys;
            int diffZ = seed.getZ() - zs;
            long distSq = diffX * diffX + diffY * diffY + diffZ * diffZ;
            if (distSq < minDistSq) {
                minDistSq = distSq;
                minDistClr = seed;
            }
        }
        
        return minDistClr;
    }
    
    public int[] getImgPixelSeedIndexes() {
        return lastImgSeedIndexes;
    }
}
