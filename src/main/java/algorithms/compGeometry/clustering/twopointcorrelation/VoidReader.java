package algorithms.compGeometry.clustering.twopointcorrelation;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.logging.Logger;

/**
 * a class to read a persisted set of background void points
 * 
 * @author nichole
 */
public class VoidReader extends AbstractVoidFinder {

    /**
     *
     * @param ois
     * @throws IOException
     */
    public VoidReader(ObjectInputStream ois) throws IOException {
        
        super();
        
        this.sampling = VoidSampling.N_A;
        
        deserializeTwoPointBackground(ois);        
    }
    
    /**
     *
     * @param sampling
     */
    public void setSampling(VoidSampling sampling) {
        // do nothing
    }
    
    /**
     *
     * @param indexer
     * @throws TwoPointVoidStatsException
     */
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
    
    /**
     *
     */
    @Override
    protected void findVoidsImpl() {
        // do nothing.  the voids were already 'found' and read in from ObjectInputStream
    }

    /**
     *
     */
    @Override
    public void constructLogger() {
        this.log = Logger.getLogger(this.getClass().getName());
    }

    /**
     *
     * @param ois
     * @throws IOException
     */
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
        try {
            String s = ois.readUTF();
            if (s != null && s.length() > 0) {
                sampling = VoidSampling.resolve(s);
            }
        } catch(IOException e) {
        }
    }

}
