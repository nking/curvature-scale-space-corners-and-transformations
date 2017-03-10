package algorithms.imageProcessing.matching;

import algorithms.QuickSort;
import algorithms.imageProcessing.VanishingPoints;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.util.PairInt;
import algorithms.util.VeryLongBitString;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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

    // vanishing points for dataset2
    private VanishingPoints vp2 = null;

    public void setVanishingPointsForSet2(VanishingPoints vp) {
        vp2 = vp;
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
     * greedy matching of d1 to d2 by min cost, with unique mappings for
     * all indexes.
     *
     * @param d1
     * @param d2
     * @return matches - two dimensional int array of indexes in d1 and
     * d2 which are matched.
     */
    public static int[][] matchDescriptors(VeryLongBitString[] d1, VeryLongBitString[] d2, List<PairInt> keypoints1, List<PairInt> keypoints2) {
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
     * @return matches - two dimensional int array of indexes in d1 and
     * d2 which are matched.
     */
    public static int[][] matchDescriptors(ORB.Descriptors[] d1, ORB.Descriptors[] d2, List<PairInt> keypoints1, List<PairInt> keypoints2) {
        if (d1.length != d2.length) {
            throw new IllegalArgumentException("d1 and d2 must" + " be same length");
        }
        int n1 = d1[0].descriptors.length;
        int n2 = d2[0].descriptors.length;
        if (n1 != keypoints1.size()) {
            throw new IllegalArgumentException("number of descriptors in " + " d1 bitstrings must be same as keypoints1 length");
        }
        if (n2 != keypoints2.size()) {
            throw new IllegalArgumentException("number of descriptors in " + " d2 bitstrings must be same as keypoints2 length");
        }
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

    private static int[][] greedyMatch(List<PairInt> keypoints1,
        List<PairInt> keypoints2, int[][] cost) {
        int n1 = keypoints1.size();
        int n2 = keypoints2.size();
        // for the greedy match, separating the index information from the cost
        // and then sorting by cost
        int nTot = n1 * n2;
        PairInt[] indexes = new PairInt[nTot];
        int[] costs = new int[nTot];
        int count = 0;
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                indexes[count] = new PairInt(i, j);
                costs[count] = cost[i][j];
                count++;
            }
        }
        assert (count == nTot);
        QuickSort.sortBy1stArg(costs, indexes);
        Set<PairInt> set1 = new HashSet<PairInt>();
        Set<PairInt> set2 = new HashSet<PairInt>();
        List<PairInt> matches = new ArrayList<PairInt>();
        // visit lowest costs (== differences) first
        for (int i = 0; i < nTot; ++i) {
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

}
