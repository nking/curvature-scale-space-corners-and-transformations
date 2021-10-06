package algorithms.imageProcessing.matching;

import algorithms.QuickSort;
import algorithms.imageProcessing.features.HOGs;
import algorithms.imageProcessing.features.RANSACSolver;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.QuadInt;
import algorithms.VeryLongBitString;
import algorithms.bipartite.Graph;
import algorithms.bipartite.MinCostUnbalancedAssignment;
import algorithms.imageProcessing.transform.EpipolarTransformer;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import no.uib.cipr.matrix.DenseMatrix;

/**
 * a class to hold various methods related to matching
 * the descriptors of ORB.
 * See also ObjectMatcher.
 *
 * ORB features can be used to match 2 images as long as the number of points
 * that are possible true matches are larger than the number of points
 * which are not matches in the images.  In other words, for the sparse
 * feature matching approach of keypoints, need the number of possible true
 * matches to be larger than the number of possible false matches.
 * If that is not the case, such as in finding an object which has changed
 * location, then the more dense approach of using blob detecter MSER is 
 * recommended.
 * 
 * NOTE that methods are being added specifically for the sparse matching.
 * 
 * @see ORB
 * @see ObjectMatcher
 *
 * @author nichole
 */
public class ORBMatcher {

    /**
     * greedy matching of d1 to d2 by min cost, with unique mappings for
     * all indexes.
     *
     * @param d1
     * @param d2
     * @return matches - two dimensional int array of indexes in d1 and
     * d2 which are matched.
     */
    public static int[][] matchDescriptors(VeryLongBitString[] d1, 
        VeryLongBitString[] d2, List<PairInt> keypoints1, List<PairInt> keypoints2) {
    
        int n1 = d1.length;
        int n2 = d2.length;
        //[n1][n2]
        int[][] cost = ORB.calcDescriptorCostMatrix(d1, d2);
        int[][] matches = greedyMatch(keypoints1, keypoints2, cost);
        // greedy or optimal match can be performed here.
        // NOTE: some matching problems might benefit from using the spatial
        //   information at the same time.  for those, will consider adding
        //   an evaluation term for these descriptors to a specialization of
        //   PartialShapeMatcher.java
        return matches;
    }

    /**
     * greedy matching of d1 to d2 by min difference, with unique mappings for
     * all indexes. NOTE that if 2 descriptors match equally well, either one
     * might get the assignment. Consider using instead, matchDescriptors2 which
     * matches by descriptor and relative spatial location.
     *
     * @param d1
     * @param d2
     * @param keypoints2
     * @param keypoints1
     * @param hogs1
     * @param hogs2
     * @param orientations1
     * @param orientations2
     * @return matches array of objects encapsulating a pair of matched points
     */
    public static QuadInt[] matchDescriptors(ORB.Descriptors d1,
        ORB.Descriptors d2, List<PairInt> keypoints1,
        List<PairInt> keypoints2, HOGs hogs1, HOGs hogs2,
        TDoubleList orientations1, TDoubleList orientations2) {

        int n1 = d1.descriptors.length;
        int n2 = d2.descriptors.length;
        if (n1 == 0 || n2 == 0) {
            return null;
        }

        if (d1.descriptors[0].getCapacity() != d2.descriptors[0].getCapacity()) {
            throw new IllegalArgumentException("d1 and d2 must have same bitstring"
                + " capacities (== 256) "
                + d1.descriptors[0].getCapacity() + " "
                + d2.descriptors[0].getCapacity()
            );
        }

        if (n1 != keypoints1.size()) {
            throw new IllegalArgumentException("number of descriptors in " + " d1 bitstrings must be same as keypoints1 length");
        }
        if (n2 != keypoints2.size()) {
            throw new IllegalArgumentException("number of descriptors in " + " d2 bitstrings must be same as keypoints2 length");
        }
        //[n1][n2]
        int[][] cost = ORB.calcDescriptorCostMatrix(
            d1.descriptors, d2.descriptors);

        // pairs of indexes of matches
        int[][] matches = greedyMatch(keypoints1, keypoints2, cost, 
            hogs1, hogs2, orientations1, orientations2);

        if (matches.length < 7) {
            // 6!/(1!5!) + 6!/(2!4!) + 6!/(3!3!) + 6!/(2!4!) + 6!/(1!5!) = 6 + 15 + 20 + 15 + 6=62
            // 5!/(1!4!) + 5!/(2!3!) + 5!/(3!2!) + 5!/(1!4!) = 5 + 10 + 10 + 5 = 20
            // 4!/(1!3!) + 4!/(2!2!) + 4!/(1!3!) = 4+6+4=14
            // 3!/(1!2!) + 3!/(2!1!) = 3 + 3
            // 2!/(1!1!) = 2
            /*
            considering how to filter for outliers in these few number of points.
            
            can use affine projection.
            can iterate over subsamples of the input point to remove points from it,
                fit and evaluate the affine projection
            
            apply the best affine projection to keypoints1 to fins close matches
               in keypoints2 where close match is (x,y) and descriptor cost.
            
            if there are a large number of matches, proceed to RANSAC below,
            else return either the matches that are the best fitting subset
            or return the subset and additional points found through the projection.
            ------
            the projection algorithm solves for rotation and real world scene coordinates.
            It does not solve for translation,
            but one could make a rough translation if the image scales are the
            same by subtracting the 2nd image correspondence points from the
            1st image correspondence points rotated.
            
            OrthographicProjectionResults re = Reconstruction.calculateAffineReconstruction(
                double[][] x, int mImages).
            
            where OrthographicProjectionResults results = new OrthographicProjectionResults();
            results.XW = s;
            results.rotationMatrices = rotStack;
            
            Considering the paper 
            "Outlier Correction in Image Sequences for the Affine Camera"
               by Huartley, and Heydeon 2003
               Proceedings of the Ninth IEEE International Conference on Computer Vision (ICCV’03)
               
               excerpt from the abstract:
                  In this paper, we present an outlier correction scheme that 
                  iteratively updates the elements of the image measurement matrix
            */

            QuadInt[] qs = new QuadInt[matches.length];
            for (int i = 0; i < matches.length; ++i) {
                int idx1 = matches[i][0];
                int idx2 = matches[i][1];
                QuadInt q = new QuadInt(
                    keypoints1.get(idx1).getX(), keypoints1.get(idx1).getY(),
                    keypoints2.get(idx2).getX(), keypoints2.get(idx2).getY()
                );
                qs[i] = q;
            }

            return qs;
        }
        
        // ransac to remove outliers.  
        //     note, the method normalizes a copy of the matched points for the algorithm
        EpipolarTransformationFit fit = fitWithRANSAC(matches, 
            keypoints1, keypoints2);
        
        if (fit == null) {
            return null;
        }
         
        List<Integer> inliers = fit.getInlierIndexes();
        QuadInt[] qs = new QuadInt[inliers.size()];
        for (int i = 0; i < inliers.size(); ++i) {
            int idx = inliers.get(i);
            QuadInt q = new QuadInt(
                keypoints1.get(idx).getX(), keypoints1.get(idx).getY(),
                keypoints2.get(idx).getX(), keypoints2.get(idx).getY()
            );
            qs[i] = q;
        }

        System.out.println("fit=" + fit.toString());
        
        return qs;
    }
    
    /**
     * greedy matching of d1 to d2 by min difference, with unique mappings for
     * all indexes.
     * NOTE that if 2 descriptors match equally well, either one
     * might get the assignment.
     * Consider using instead, matchDescriptors2 which matches
     * by descriptor and relative spatial location.
     *
     * @param d1
     * @param d2
     * @param keypoints2
     * @param keypoints1
     * @return matches array of objects encapsulating a pair of
     * matched points
     */
    public static QuadInt[] matchDescriptors(ORB.Descriptors d1, 
        ORB.Descriptors d2, List<PairInt> keypoints1, 
        List<PairInt> keypoints2) {
        
        int n1 = d1.descriptors.length;
        int n2 = d2.descriptors.length;
        if (n1 == 0 || n2 == 0) {
            return null;
        }
        
        if (d1.descriptors[0].getCapacity() != d2.descriptors[0].getCapacity()) {
            throw new IllegalArgumentException("d1 and d2 must have same bitstring" 
                + " capacities (== 256) " + 
                d1.descriptors[0].getCapacity() + " " +
                d2.descriptors[0].getCapacity()
            );
        }
        
        if (n1 != keypoints1.size()) {
            throw new IllegalArgumentException("number of descriptors in " + " d1 bitstrings must be same as keypoints1 length");
        }
        if (n2 != keypoints2.size()) {
            throw new IllegalArgumentException("number of descriptors in " + " d2 bitstrings must be same as keypoints2 length");
        }
        //[n1][n2]
        int[][] cost = ORB.calcDescriptorCostMatrix(
            d1.descriptors, d2.descriptors);
       
        // pairs of indexes of matches
        int[][] matches = greedyMatch(keypoints1, keypoints2, cost);
        
        if (matches.length < 7) {
            
            QuadInt[] qs = new QuadInt[matches.length];
            for (int i = 0; i < matches.length; ++i) {
                int idx1 = matches[i][0];
                int idx2 = matches[i][1];
                QuadInt q = new QuadInt(
                    keypoints1.get(idx1).getX(), keypoints1.get(idx1).getY(),
                    keypoints2.get(idx2).getX(), keypoints2.get(idx2).getY()
                );
                qs[i] = q;
            }

            return qs;
        }
         
        // ransac to remove outliers.
        //NOTE: the fundamental matrix in the fit has not been de-normalized.
        EpipolarTransformationFit fit = fitWithRANSAC(matches, 
            keypoints1, keypoints2);
        
        if (fit == null) {
            return null;
        }
        
        System.out.println("fit=" + fit.toString());

        int i, idx;
        List<Integer> inliers = fit.getInlierIndexes();
        QuadInt[] qs = new QuadInt[inliers.size()];
        for (i = 0; i < inliers.size(); ++i) {
            idx = inliers.get(i);
            QuadInt q = new QuadInt(
                keypoints1.get(idx).getX(), keypoints1.get(idx).getY(),
                keypoints2.get(idx).getX(), keypoints2.get(idx).getY()
            );
            qs[i] = q;
        }
        
        return qs;        
    }

    /**
     * finds best match for each point if a close second best does not exist,
     * then sorts by lowest cost to keep the unique best starter points.
     * returns matching indexes (no ransac performed in this method)
     * @param keypoints1
     * @param keypoints2
     * @param cost
     * @return 
     */
    private static int[][] greedyMatch(List<PairInt> keypoints1,
        List<PairInt> keypoints2, int[][] cost) {
        
        int n1 = keypoints1.size();
        int n2 = keypoints2.size();
        
        /*
        -- for each keypoint, finding best match, but only keeping it if there is
           no close 2nd best.
        -- sorting the results by lowest cost and keepint the unique of those.
        -- return correspondence
        */
        
        //nearest neighbor distance ratio (Mikolajczyk and Schmid 2005):
        // using a ratio of 0.8 or 0.9.
        int[] bestMatch = findGreedyBestIsolated(cost, 0.8f);
        
        //int[] bestMatch = minCostBipartiteUnbalanced(cost);
        
        assert(bestMatch.length == n1);
        
        int nBest = 0;
        for (int idx : bestMatch) {
            if (idx > -1) {
                nBest++;
            }
        }
        
        PairInt[] indexes = new PairInt[nBest];
        int[] costs = new int[nBest];
        int count = 0;
        int idx1, idx2;
        for (idx1 = 0; idx1 < bestMatch.length; ++idx1) {
            idx2 = bestMatch[idx1];
            if (idx2 > -1) {
                indexes[count] = new PairInt(idx1, idx2);
                costs[count] = cost[idx1][idx2];
                count++;
            }
        }
        
        assert(count == nBest);
        QuickSort.sortBy1stArg(costs, indexes);
        Set<PairInt> set1 = new HashSet<PairInt>();
        Set<PairInt> set2 = new HashSet<PairInt>();
        List<PairInt> matches = new ArrayList<PairInt>();
        PairInt index12, p1, p2;
        // visit lowest costs (== differences) first
        for (int i = 0; i < nBest; ++i) {
            index12 = indexes[i];
            idx1 = index12.getX();
            idx2 = index12.getY();
            p1 = keypoints1.get(idx1);
            p2 = keypoints2.get(idx2);
            if (set1.contains(p1) || set2.contains(p2)) {
                continue;
            }
            //System.out.println("p1=" + p1 + " " + " p2=" + p2 + " cost=" + costs[i]);
            matches.add(index12);
            set1.add(p1);
            set2.add(p2);
        }
        int[][] results = new int[matches.size()][2];
        for (int i = 0; i < matches.size(); ++i) {
            results[i][0] = matches.get(i).getX();
            results[i][1] = matches.get(i).getY();
        }
        return results;
    }
    
    /**
     * finds best match for each point if a close second best does not exist,
     * then sorts by lowest cost to keep the unique best starter points.
     * returns matching indexes (no ransac performed in this method)
     * @param keypoints1
     * @param keypoints2
     * @param cost
     * @return 
     */
    private static int[][] greedyMatch(List<PairInt> keypoints1,
        List<PairInt> keypoints2, int[][] cost, HOGs hogs1, HOGs hogs2,
        TDoubleList orientations1, TDoubleList orientations2) {
        
        int n1 = keypoints1.size();
        int n2 = keypoints2.size();
        
        /*
        -- for each keypoint, finding best match, but only keeping it if there is
           no close 2nd best.
        -- sorting the results by lowest cost and keepint the unique of those.
        -- return correspondence
        */
        
        //nearest neighbor distance ratio (Mikolajczyk and Schmid 2005):
        // using a ratio of 0.8 or 0.9.
        TIntList outputI1 = new TIntArrayList();
        TIntList outputI2 = new TIntArrayList();
        TFloatList outputCost = new TFloatArrayList();
        
        findGreedyBest(cost, 0.8f,
            hogs1, hogs2, keypoints1, keypoints2,
            orientations1, orientations2,
            outputI1, outputI2, outputCost);
                
        int nBest = outputI1.size();
        
        PairInt[] indexes = new PairInt[nBest];
        float[] costs = new float[nBest];
        int count = 0;
        for (int i = 0; i < nBest; ++i) {
            int idx1 = outputI1.get(i);
            int idx2 = outputI2.get(i);
            if (idx2 > -1) {
                indexes[count] = new PairInt(idx1, idx2);
                costs[count] = outputCost.get(i);
                count++;
            }
        }
        
        assert (count == nBest);
        QuickSort.sortBy1stArg(costs, indexes);
        Set<PairInt> set1 = new HashSet<PairInt>();
        Set<PairInt> set2 = new HashSet<PairInt>();
        List<PairInt> matches = new ArrayList<PairInt>();
        // visit lowest costs (== differences) first
        for (int i = 0; i < nBest; ++i) {
            PairInt index12 = indexes[i];
            int idx1 = index12.getX();
            int idx2 = index12.getY();
            PairInt p1 = keypoints1.get(idx1);
            PairInt p2 = keypoints2.get(idx2);
            if (set1.contains(p1) || set2.contains(p2)) {
                continue;
            }
            //System.out.println("p1=" + p1 + " " + " p2=" + p2 + " cost=" + costs[i]);
            matches.add(index12);
            set1.add(p1);
            set2.add(p2);
        }
        int[][] results = new int[matches.size()][2];
        for (int i = 0; i < matches.size(); ++i) {
            results[i][0] = matches.get(i).getX();
            results[i][1] = matches.get(i).getY();
        }
        return results;
    }
    
    private static int[] minCostBipartiteUnbalanced(int[][] cost) {
        
        TObjectIntMap<PairInt> weights = new TObjectIntHashMap<PairInt>();
            
        int i, j;
        for (i = 0; i < cost.length; ++i) {
            for (j = 0; j < cost[i].length; ++j) {
                weights.put(new PairInt(i, j), cost[i][j]);
            }
        }
        boolean createSourceAndSinkEdges = true;
        Graph g = new Graph(cost.length, cost[0].length, weights, createSourceAndSinkEdges);
        
        MinCostUnbalancedAssignment bipartite = 
            new MinCostUnbalancedAssignment();
        
        TIntIntMap map = bipartite.flowAssign(g);
        
        int[] bestMatch = new int[cost.length];
        Arrays.fill(bestMatch, -1);
        
        TIntIntIterator iter = map.iterator();
        for (i = 0; i < map.size(); ++i) {
            iter.advance();
            bestMatch[iter.key()] = iter.value();
        }
        return bestMatch;
    }
    
    //@param ratioLimit Mikolajczyk and Schmid 2005) 0.8 or 0.9.
    private static int[] findGreedyBestIsolated(int[][] cost, float ratioLimit) {
        int n1 = cost.length;
        int n2 = cost[0].length;
                
        // best match cost
        int bc;
        // best match index
        int bcIdx;
        // 2nd best match cost
        int bc2;
        // 2nd best match index
        int bc2Idx;
        int c;
        int[] bestMatch = new int[n1];
        int i, j;
        for (i = 0; i < n1; ++i) {
            bc = Integer.MAX_VALUE;
            bc2 = Integer.MAX_VALUE;
            bcIdx = -1;
            bc2Idx = -1;
            for (j = 0; j < n2; ++j) {
                c = cost[i][j];
                if (c >= bc2) {
                    continue;
                }
                if (c < bc) {
                    bc2 = bc;
                    bc2Idx = bcIdx;
                    bc = c;
                    bcIdx = j;
                } else if (c == bc) {
                    if (c < bc2) {
                        bc2 = bc;
                        bc2Idx = bcIdx;
                        bc = c;
                        bcIdx = j;
                    } else {
                        assert(c == bc2 && bc == bc2);
                    }
                } else {
                    // c > bc
                    if (c < bc2) {
                        bc2 = c;
                        bc2Idx = j;
                    }
                }
            }
            if (bc2Idx == -1) {
                bestMatch[i] = bcIdx;
            } else {
                float ratio = (float)bc/(float)bc2;
                if (ratio < ratioLimit) {
                    bestMatch[i] = bcIdx;
                } else {
                    bestMatch[i] = -1;
                }
            }
        }
        
        return bestMatch;
    }  
    
    private static void findGreedyBest(int[][] cost, float ratioLimit,
        HOGs hogs1, HOGs hogs2, 
        List<PairInt> keypoints1, List<PairInt> keypoints2,
        TDoubleList orientations1, TDoubleList orientations2,
        TIntList outputI1, TIntList outputI2, TFloatList outputCost) {
        
        int n1 = cost.length;
        int n2 = cost[0].length;
        
        int[] h1 = new int[hogs1.getNumberOfBins()];
        int[] h2 = new int[h1.length];
                
        for (int i = 0; i < n1; ++i) {
            float bc = Integer.MAX_VALUE;
            float bc2 = Integer.MAX_VALUE;
            int bcIdx = -1;
            int bc2Idx = -1;
            
            int or1 = (int)Math.round(orientations1.get(i) * 180./Math.PI);
            if (or1 < 0) {
                or1 += 180;
            }
            hogs1.extractBlock(keypoints1.get(i).getX(), 
                keypoints1.get(i).getY(), h1);
            
            for (int j = 0; j < n2; ++j) {
                
                float c = cost[i][j]/255.f;
                
                int or2 = (int)Math.round(orientations2.get(j) * 180./Math.PI);
                if (or2 < 0) {
                    or2 += 180;
                }
                hogs2.extractBlock(keypoints2.get(j).getX(), 
                    keypoints2.get(j).getY(), h2);
                
                float intersection = hogs1.intersection(h1, or1, h2, or2);
                                
                c += (1.f - intersection);
                
                if (c >= bc2) {
                    continue;
                }
                if (c < bc) {
                    bc2 = bc;
                    bc2Idx = bcIdx;
                    bc = c;
                    bcIdx = j;
                } else if (c == bc) {
                    if (c < bc2) {
                        bc2 = bc;
                        bc2Idx = bcIdx;
                        bc = c;
                        bcIdx = j;
                    } else {
                        assert(c == bc2 && bc == bc2);
                    }
                } else {
                    // c > bc
                    if (c < bc2) {
                        bc2 = c;
                        bc2Idx = j;
                    }
                }
            }
            if (bc2Idx == -1) {
                outputI1.add(i);
                outputI2.add(bcIdx);
                outputCost.add(bc);
            } else {
                float ratio = bc/bc2;
                if (ratio < 0.8) {
                    outputI1.add(i);
                    outputI2.add(bcIdx);
                    outputCost.add(bc);
                }
            }
        }        
    }  
    
    public static double distance(int x, int y, PairInt b) {
        int diffX = x - b.getX();
        int diffY = y - b.getY();
        double dist = Math.sqrt(diffX * diffX + diffY * diffY);
        return dist;
    }

    public static int distance(PairInt p1, PairInt p2) {
        int diffX = p1.getX() - p2.getX();
        int diffY = p1.getY() - p2.getY();
        return (int) Math.sqrt(diffX * diffX + diffY * diffY);
    }

    /**
     * calculate the fundamental matrix given the correspondence matches.
     * the correspondence is normalized and the fundamental matrix is calculated,
     * then the errors are estimated using the Sampson's distance as errors
     * and a 3.8*sigma as inlier threshold.
     * @param matches
     * @param keypoints1
     * @param keypoints2
     * @return epipolar fit to the matches.  note that the correspondence
     * is unit standard normalized and the fundamental matrix returned
     * in the fit is not de-normalized.
     */
    private static EpipolarTransformationFit fitWithRANSAC(int[][] matches, 
        List<PairInt> keypoints1, List<PairInt> keypoints2) {
        
        int n0 = matches.length;

        double[][] left = new double[3][n0];
        double[][] right = new double[3][n0];
        for (int i = 0; i < 3; ++i) {
            left[i] = new double[n0];
            right[i] = new double[n0];
        }
        Arrays.fill(left[2], 1.0);
        Arrays.fill(right[2], 1.0);
        int i, idx1, idx2;
        for (i = 0; i < n0; ++i) {
            idx1 = matches[i][0];
            idx2 = matches[i][1];
            left[0][i] = keypoints1.get(idx1).getX();
            left[1][i] = keypoints1.get(idx1).getY();
            right[0][i] = keypoints2.get(idx2).getX();
            right[1][i] = keypoints2.get(idx2).getY();
        }
    
        // normalize left and right
        boolean useToleranceAsStatFactor = true;
        final double tolerance = 3.8;
        ErrorType errorType = ErrorType.SAMPSONS;
        
        EpipolarTransformer.NormalizedXY normXY1 = EpipolarTransformer.normalize(new DenseMatrix(left));
        EpipolarTransformer.NormalizedXY normXY2 = EpipolarTransformer.normalize(new DenseMatrix(right));
        DenseMatrix leftM = normXY1.getXy();
        DenseMatrix rightM = normXY2.getXy();
        
        boolean reCalcIterations = false;
        RANSACSolver solver = new RANSACSolver();
        
        EpipolarTransformationFit fit = solver.calculateEpipolarProjection(
            leftM, rightM, errorType, useToleranceAsStatFactor, tolerance,
            reCalcIterations, false);                
        
        return fit;        
    }

}
