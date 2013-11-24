package algorithms.compGeometry.clustering.twopointcorrelation;

/**
 *
 * @author nichole
 */
public interface ITwoPointIdentity {

    public long approximateMemoryUsed();

    public boolean storeIfDoesNotContain(int index0, int index1);

}
