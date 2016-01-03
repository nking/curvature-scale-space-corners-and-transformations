package algorithms.imageProcessing;

import algorithms.compGeometry.PerimeterFinder;
import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import thirdparty.HungarianAlgorithm;

/**
 *
 * @author nichole
 */
public class BlobCornerFinderForParameters {
    
    public List<FeatureComparisonStat> extractFeatures(
        TransformationParameters parameters, 
        GreyscaleImage img1, GreyscaleImage img2,
        GreyscaleImage segImg1, GreyscaleImage segImg2, 
        int binFactor1, int binFactor2,
        int smallestGroupLimit, int largestGroupLimit,
        RotatedOffsets rotatedOffsets, boolean filterOutImageBoundaryBlobs,
        boolean debug, String debugTag) {
        
        List<FeatureComparisonStat> allStats = new ArrayList<FeatureComparisonStat>();
        
        List<Set<PairInt>> blobs1 = new ArrayList<Set<PairInt>>();
        
        List<Set<PairInt>> blobs2 = new ArrayList<Set<PairInt>>();
        
        extractBlobs(blobs1, blobs2, segImg1, segImg2, binFactor1, binFactor2,
            smallestGroupLimit, largestGroupLimit, filterOutImageBoundaryBlobs,
            debug, debugTag);
        
        Transformer transformer = new Transformer();
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
                
        double[][] xyCen1 = new double[blobs1.size()][];
        double[][] xyCenTr1 = new double[blobs1.size()][];
        double[][] xyCen2 = new double[blobs2.size()][];
        for (int i = 0; i < blobs1.size(); ++i) {
            xyCen1[i] = curveHelper.calculateXYCentroids(blobs1.get(i));
            xyCenTr1[i] = transformer.applyTransformation(parameters, 
                xyCen1[i][0], xyCen1[i][1]);
        }
        for (int i = 0; i < blobs2.size(); ++i) {
            xyCen2[i] = curveHelper.calculateXYCentroids(blobs2.get(i));
        }
        
        float[] stdevs = parameters.getStandardDeviations();
        int transXYTol = Math.round(Math.max(stdevs[2], stdevs[3]));
        
        List<Integer> matchedIndexes1 = new ArrayList<Integer>();
        List<Integer> matchedIndexes2 = new ArrayList<Integer>();
        
        //bipartiteMatching(xyCenTr1, xyCen2, transXYTol, matchedIndexes1, matchedIndexes2);
        greedyDegenerateMatching(xyCenTr1, xyCen2, transXYTol, matchedIndexes1, 
            matchedIndexes2);
        
        if (matchedIndexes1.isEmpty()) {
            return allStats;
        }
        
        CornerRegion[][] blobs1Index1BPRs = 
            BlobsAndPerimeters.extractBlobPerimeterAsCornerRegions(blobs1, 
            matchedIndexes1, img1.getWidth(), img1.getHeight());
        
        CornerRegion[][] blobs1Index1BPRsTR = transform(blobs1Index1BPRs, 
            parameters);
        
        CornerRegion[][] blobs2Index2BPRs = 
            BlobsAndPerimeters.extractBlobPerimeterAsCornerRegions(blobs2, 
            matchedIndexes2, img2.getWidth(), img2.getHeight());
        
        Map<Integer, IntensityFeatureComparisonStats> blobs1Index1Matches = 
            new HashMap<Integer, IntensityFeatureComparisonStats>();
        
        Map<Integer, IntensityFeatureComparisonStats> blobs2Index2Matches = 
            new HashMap<Integer, IntensityFeatureComparisonStats>();
                
        FeatureMatcher featureMatcher = new FeatureMatcher();
        final int blockHalfWidth = 5;
        final boolean useNormalizedIntensities = true;
        
        IntensityFeatures features1 = new IntensityFeatures(blockHalfWidth, 
            useNormalizedIntensities, rotatedOffsets);
        features1.calculateGradientWithGreyscale(img1);
        
        IntensityFeatures features2 = new IntensityFeatures(blockHalfWidth,
            useNormalizedIntensities, rotatedOffsets);
        features2.calculateGradientWithGreyscale(img2);
        
        int dither = 1;
        double solutionScale = parameters.getScale();
        
        for (int i = 0; i < matchedIndexes1.size(); ++i) {
            
            Integer index1 = matchedIndexes1.get(i);
            Integer index2 = matchedIndexes2.get(i);
            
            CornerRegion[] cr1 = blobs1Index1BPRs[index1.intValue()];
            CornerRegion[] cr1Tr = blobs1Index1BPRsTR[index1.intValue()];
            CornerRegion[] cr2 = blobs2Index2BPRs[index2.intValue()];
            
            // compare pixels whose SSD error is relatively high, meaning,
            // they aren't featureless patches as seen by autocorrelation
            List<FeatureComparisonStat> stats = 
                featureMatcher.findSimilarFeaturesAsStats(
                img1, cr1, cr1Tr, img2, cr2, features1, features2, 
                parameters, dither, transXYTol, rotatedOffsets);
            
            if (stats.isEmpty()) {
                continue;
            }
            
            double solutionCost = calculateCombinedIntensityStat(stats);
            
            IntensityFeatureComparisonStats ifcs = 
                new IntensityFeatureComparisonStats(index1.intValue(), 
                    index2.intValue(), solutionCost, solutionScale);
            ifcs.addAll(stats);
                   
            IntensityFeatureComparisonStats ifcs1 = blobs1Index1Matches.get(index1);
            IntensityFeatureComparisonStats ifcs2 = blobs2Index2Matches.get(index2);
            
            if (ifcs1 != null) {
                if (ifcs2 != null && ifcs2.equals(ifcs1)) {
                    if (ifcs.getCost() < ifcs1.getCost()) {
                        blobs1Index1Matches.put(index1, ifcs);
                        blobs2Index2Matches.put(index2, ifcs);
                    }
                    continue;
                }
                if (ifcs.getCost() < ifcs1.getCost()) {
                    blobs1Index1Matches.put(index1, ifcs);
                    blobs2Index2Matches.put(index2, ifcs);
                }
                continue;
            } 
            if (ifcs2 != null) {
                if (ifcs.getCost() < ifcs2.getCost()) {
                    blobs1Index1Matches.put(index1, ifcs);
                    blobs2Index2Matches.put(index2, ifcs);
                }
                continue;
            }
            blobs1Index1Matches.put(index1, ifcs);
            blobs2Index2Matches.put(index2, ifcs);
        }
        
        for (Entry<Integer, IntensityFeatureComparisonStats> entry :
            blobs1Index1Matches.entrySet()) {
            
            allStats.addAll(entry.getValue().getComparisonStats());
        }
        
        return allStats;
    }
    
    private void extractBlobs(List<Set<PairInt>> outputBlobs1, 
        List<Set<PairInt>> outputBlobs2,
        GreyscaleImage img1, GreyscaleImage img2, int binFactor1, int binFactor2,
        int smallestGroupLimit, int largestGroupLimit, 
        boolean filterOutImageBoundaryBlobs, boolean debug, String debugTag) {
        
        List<Set<PairInt>> blobs1 = 
            BlobsAndPerimeters.extractBlobsKeepBounded(
                img1, smallestGroupLimit, largestGroupLimit, binFactor1,
                filterOutImageBoundaryBlobs);
        
        outputBlobs1.addAll(blobs1);
        
        List<Set<PairInt>> blobs2 = 
            BlobsAndPerimeters.extractBlobsKeepBounded(
                img2, smallestGroupLimit, largestGroupLimit, binFactor2,
                filterOutImageBoundaryBlobs);
        
        outputBlobs2.addAll(blobs2);
        
        if (debug) {
            Image img0 = ImageIOHelper.convertImage(img1);
            int c = 0;
            for (int i = 0; i < blobs1.size(); ++i) {
                Set<PairInt> blobSet = blobs1.get(i);
                int clr = ImageIOHelper.getNextColorRGB(c);
                for (PairInt p : blobSet) {
                    int x = p.getX();
                    int y = p.getY();
                    ImageIOHelper.addPointToImage(x, y, img0, 0, clr);
                }
                c++;
            }
            
            MiscDebug.writeImage(img0, "blobs1_ext_seg" + debugTag + "_" + 
                MiscDebug.getCurrentTimeFormatted());
            
            img0 = ImageIOHelper.convertImage(img2);
            c = 0;
            for (int i = 0; i < blobs2.size(); ++i) {
                Set<PairInt> blobSet = blobs2.get(i);
                int clr = ImageIOHelper.getNextColorRGB(c);
                for (PairInt p : blobSet) {
                    int x = p.getX();
                    int y = p.getY();
                    ImageIOHelper.addPointToImage(x, y, img0, 0, clr);
                }
                c++;
            }
            
            MiscDebug.writeImage(img0, "blobs2_ext_seg" + debugTag + "_" + 
                MiscDebug.getCurrentTimeFormatted());
        }
    }

    private void bipartiteMatching(double[][] xy1, double[][] xy2, 
        int transXYTol, List<Integer> outputMatchedIndexes1, 
        List<Integer> outputMatchedIndexes2) {
        
        int n1 = xy1.length;
        int n2 = xy2.length;
        
        float[][] cost = new float[n1][n2];
        float[][] costCopy = new float[n1][n2];
        for (int i = 0; i < n1; ++i) {
            cost[i] = new float[n2];
            costCopy[i] = new float[n2];
            for (int j = 0; j < n2; ++j) {
                double dist = distance(xy1[i], xy2[i]);
                if (dist > transXYTol) {
                    cost[i][j] = Float.MAX_VALUE;
                } else {
                    cost[i][j] = (float)dist;
                }
                costCopy[i][j] = cost[i][j];
            }
        }
        
        boolean transposed = false;
        if (cost.length > cost[0].length) {
            cost = MatrixUtil.transpose(cost);
            transposed = true;
        }

        // one pass thru to count for array sizes
        HungarianAlgorithm b = new HungarianAlgorithm();
        int[][] match = b.computeAssignments(cost);

        for (int i = 0; i < match.length; i++) {
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }
            if (costCopy[idx1][idx2] == Float.MAX_VALUE) {
                continue;
            }
            outputMatchedIndexes1.add(Integer.valueOf(idx1));
            outputMatchedIndexes2.add(Integer.valueOf(idx2));
        }
        
    }
    
    private void greedyDegenerateMatching(double[][] xy1, double[][] xy2, 
        int transXYTol, List<Integer> outputMatchedIndexes1, 
        List<Integer> outputMatchedIndexes2) {
        
        int n1 = xy1.length;
        int n2 = xy2.length;
        
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                double dist = distance(xy1[i], xy2[j]);
                if (dist > transXYTol) {
                    continue;
                }
                outputMatchedIndexes1.add(Integer.valueOf(i));
                outputMatchedIndexes2.add(Integer.valueOf(j));
            }
        }
    }

    private double distance(double[] xy1, double[] xy2) {

        double diffX = xy1[0] - xy2[0];
        double diffY = xy1[1] - xy2[1];

        double dist = Math.sqrt(diffX * diffX + diffY * diffY);
        
        return dist;
    }

    private CornerRegion[][] transform(CornerRegion[][] crs, 
        TransformationParameters parameters) {
        
        CornerRegion[][] tr = new CornerRegion[crs.length][];
        
        Transformer transformer = new Transformer();
        
        for (int i = 0; i < crs.length; ++i) {
            
            CornerRegion[] crA = crs[i];
            if (crA == null) {
                continue;
            }
            
            tr[i] = new CornerRegion[crA.length];
            
            for (int j = 0; j < crA.length; ++j) {
                
                CornerRegion cr = crA[j];
                int kIdx = cr.getKMaxIdx();
                
                double[] xy = transformer.applyTransformation(parameters, 
                    cr.getX()[kIdx], cr.getY()[kIdx]);
                
                CornerRegion crTr = new CornerRegion(i, 1, 0);
                crTr.setFlagThatNeighborsHoldDummyValues();
                crTr.set(0, Float.MIN_VALUE, 
                    (int)Math.round(xy[0]), (int)Math.round(xy[1]));
                crTr.setIndexWithinCurve(j);
                
                tr[i][j] = crTr;
            }
        }
        
        return tr;
    }
    
    protected double calculateCombinedIntensityStat(List<FeatureComparisonStat> 
        compStats) {
        
        double sum = 0;
        
        for (FeatureComparisonStat compStat : compStats) {
            sum += compStat.getSumIntensitySqDiff();
        }
        
        sum /= (double) compStats.size();
        
        return sum;
    }
}
