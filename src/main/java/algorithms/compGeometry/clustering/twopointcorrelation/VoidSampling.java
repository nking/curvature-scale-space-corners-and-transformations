package algorithms.compGeometry.clustering.twopointcorrelation;

public enum VoidSampling {

    COMPLETE, 
    LEAST_COMPLETE, 
    COMPLETE_ON_SUBSET,
    FOR_SPARSE_BACKGROUND,
    N_A;
    
    public static VoidSampling resolve(String name) {
        for (VoidSampling vs : VoidSampling.values()) {
            if (vs.name().equals(name)) {
                return vs;
            }
        }
        return null;
    }
    
}
