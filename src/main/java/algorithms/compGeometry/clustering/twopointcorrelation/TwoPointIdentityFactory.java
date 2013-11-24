package algorithms.compGeometry.clustering.twopointcorrelation;

/**
 * creates an instance of ITwoPointIdentity
 *
 * @author nichole
 */
public class TwoPointIdentityFactory {

    public static ITwoPointIdentity create(int indexerNXY) {

        if (indexerNXY <= TwoPointHashMap.nMax) {

            return new TwoPointHashMap(indexerNXY);

        } else {

            return new TwoPointBinarySearchTree();
        }
    }
}
