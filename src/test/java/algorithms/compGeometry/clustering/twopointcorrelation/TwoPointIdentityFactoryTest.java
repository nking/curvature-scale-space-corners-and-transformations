package algorithms.compGeometry.clustering.twopointcorrelation;

import junit.framework.TestCase;

public class TwoPointIdentityFactoryTest extends TestCase {

    public void testCreate() throws Exception {
                        
        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        int nPoints = 5000;
        float[] x = new float[nPoints];
        float[] y = new float[nPoints];
        indexer.sortAndIndexXThenY(x, y, x.length);
        
        TwoPointBinarySearchTree bt = (TwoPointBinarySearchTree) TwoPointIdentityFactory.create(indexer.getNXY());
        
        assertNotNull(bt);
        
    }
}
