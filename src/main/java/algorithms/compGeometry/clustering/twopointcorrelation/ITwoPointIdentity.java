package algorithms.compGeometry.clustering.twopointcorrelation;

/**
 * interface for classes used to hold a point identity pair uniquely.
 * 
 * @author nichole
 */
public interface ITwoPointIdentity {

    public long approximateMemoryUsed();

    public boolean storeIfDoesNotContain(int index0, int index1);

}
