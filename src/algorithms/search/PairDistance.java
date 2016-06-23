package algorithms.search;

import algorithms.util.PairInt;

/**
 *
 * @author nichole
 */
public class PairDistance implements Comparable<PairDistance> {

    final PairInt p2;
    final int dist;

    public PairDistance(PairInt p2, int dist) {
        this.p2 = p2;
        this.dist = dist;
    }

    @Override
    public int compareTo(PairDistance other) {

        if (dist < other.dist) {
            return -1;
        } else if (dist == other.dist) {
            return 0;
        }
        return 1;
    }

}
