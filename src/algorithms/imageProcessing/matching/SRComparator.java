package algorithms.imageProcessing.matching;

import java.util.Comparator;

/**
 *
 * @author nichole
 */
class SRComparator implements Comparator<SR> {

    @Override
    public int compare(SR o1, SR o2) {
        double d1 = o1.calcSalukDist();
        double d2 = o2.calcSalukDist();
        if (d1 < d2) {
            return -1;
        } else if (d1 > d2) {
            return 1;
        }
        return 0;
    }


}
