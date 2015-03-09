package algorithms.compGeometry.clustering.twopointcorrelation;

import junit.framework.TestCase;

public class TwoPointIdentityFactoryTest extends TestCase {

    public void testCreate() throws Exception {
                        
        AxisIndexer indexer = new AxisIndexer();
        int nPoints = 46350;
        float[] x = new float[nPoints];
        float[] y = new float[nPoints];
        indexer.sortAndIndexX(x, y, x.length);
        
        ITwoPointIdentity bt = TwoPointIdentityFactory.create(indexer.getNXY());
        
        assertTrue(bt instanceof TwoPointHash);
        
       
        nPoints = 1000;
        x = new float[nPoints];
        y = new float[nPoints];
        indexer.sortAndIndexX(x, y, x.length);
        
        ITwoPointIdentity hm = TwoPointIdentityFactory.create(indexer.getNXY());
        assertTrue(hm instanceof TwoPointHashMap);
        
        // make coverage reports happy: work around for class built purely for static methods. constructor not directly used...
        TwoPointIdentityFactory tf = new TwoPointIdentityFactory();
        assertNotNull(tf);
    }
}
