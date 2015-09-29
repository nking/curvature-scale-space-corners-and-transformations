package algorithms.compGeometry.clustering.twopointcorrelation;

/**
 * interface for classes used to hold a point identity pair uniquely.
 * 
 * @author nichole
 */
public interface ITwoPointIdentity {

    /**
     *
     * @return
     */
    public long approximateMemoryUsed();

    /**
     *
     * @param index0
     * @param index1
     * @return
     */
    public boolean storeIfDoesNotContain(int index0, int index1);

}
