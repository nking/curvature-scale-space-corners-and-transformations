package algorithms.compGeometry.clustering.twopointcorrelation;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.logging.Logger;

public class VoidReader extends AbstractVoidFinder {

    public VoidReader(ObjectInputStream ois) throws IOException {
        
        super();
        
        deserializeTwoPointBackground(ois);
        
        this.sampling = VoidSampling.N_A;
    }
    
    public void setSampling(VoidSampling sampling) {
        // do nothing
    }
    
    public void findVoids(AxisIndexer indexer) throws TwoPointVoidStatsException {
        
        if (indexer == null) {
            throw new IllegalArgumentException("indexer cannot be null");
        }
        
        this.indexer = indexer;
        
        constructLogger();
        
        //initializeVariables();
        
        findVoidsImpl();
        
        condenseArrays();
    }
    
    @Override
    protected void findVoidsImpl() {
        // do nothing.  the voids were already 'found' and read in from ObjectInputStream
    }

    @Override
    public void constructLogger() {
        this.log = Logger.getLogger(this.getClass().getName());
    }

    protected void deserializeTwoPointBackground(ObjectInputStream ois) throws IOException {
        
        this.nTwoPointSurfaceDensities = ois.readInt();

        this.allTwoPointSurfaceDensities = new float[nTwoPointSurfaceDensities];
        this.allTwoPointSurfaceDensitiesErrors = new float[nTwoPointSurfaceDensities];
        this.point1 = new int[nTwoPointSurfaceDensities];
        this.point2 = new int[nTwoPointSurfaceDensities];

        for (int i = 0; i < nTwoPointSurfaceDensities; i++) {
            this.allTwoPointSurfaceDensities[i] = ois.readFloat();
        }
        for (int i = 0; i < nTwoPointSurfaceDensities; i++) {
            this.point1[i] = ois.readInt();
            this.point2[i] = ois.readInt();
        }
        for (int i = 0; i < nTwoPointSurfaceDensities; i++) {
            this.allTwoPointSurfaceDensitiesErrors[i] = ois.readFloat();
        }
    }

}
