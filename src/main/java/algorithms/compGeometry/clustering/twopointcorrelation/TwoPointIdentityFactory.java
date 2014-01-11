package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.util.Util;

/**
 * creates an instance of ITwoPointIdentity
 *
 * @author nichole
 */
public class TwoPointIdentityFactory {
    
    public static ITwoPointIdentity create(int indexerNXY) {

        if (indexerNXY <= TwoPointHashMap.nMax) {

            long requiredMemory = TwoPointHashMap.checkRequiredMemory((long)indexerNXY);

            long totAvail = Util.getAvailableHeapMemory();

            java.util.logging.Logger.getLogger(TwoPointIdentityFactory.class.getName()).fine(
                " memory required = "  + (requiredMemory/1000000.) + " [MB]\n" +
                "       memory available = " + (totAvail/1000000.) + " [MB] ");

            if (requiredMemory < totAvail) {

                return new TwoPointHashMap(indexerNXY);

            } else {

                return new TwoPointBinarySearchTree();
            } 

        } else {

            return new TwoPointBinarySearchTree();
        }
    }
}
