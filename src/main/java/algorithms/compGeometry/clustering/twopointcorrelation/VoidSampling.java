package algorithms.compGeometry.clustering.twopointcorrelation;

/**
 * enumeration of types of sampling for the background void points
 * 
 * @author nichole
 */
public enum VoidSampling {

    /**
     *
     */
    COMPLETE, 

    /**
     *
     */
    LEAST_COMPLETE, 

    /**
     *
     */
    COMPLETE_ON_SUBSET,

    /**
     *
     */
    FOR_SPARSE_BACKGROUND,

    /**
     *
     */
    N_A;
    
    /**
     *
     * @param name
     * @return
     */
    public static VoidSampling resolve(String name) {
        for (VoidSampling vs : VoidSampling.values()) {
            if (vs.name().equals(name)) {
                return vs;
            }
        }
        return null;
    }
    
}