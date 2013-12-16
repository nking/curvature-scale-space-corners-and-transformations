package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class CompleteSamplingVoidFinderTest extends BaseTwoPointTest {

    boolean debug = true;

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public void testFindVoids_0() throws Exception {

        log.info("testFindVoids_0()");

        float xmin = 0;
        float xmax = 3;
        float ymin = 0;
        float ymax = 3;

        int numberOfBackgroundPoints = 9;

        float[] xb = new float[numberOfBackgroundPoints];
        float[] yb = new float[numberOfBackgroundPoints];
        // make a uniform grid of background points:
        int nDiv = (int) Math.ceil(Math.sqrt(numberOfBackgroundPoints));
        double divXSz = (xmax - xmin)/nDiv;
        double divYSz = (ymax - ymin)/nDiv;
        int c = 0;
        for (int j = 0; j < nDiv; j++) {
            float yStart = (float) (ymin + j*divYSz);
            if (yStart > ymax) {
                yStart = ymax;
            }
            for (int jj = 0; jj < nDiv; jj++) {
                float xStart = (float)(xmin + jj*divXSz);
                if (xStart > xmax) {
                    xStart = xmax;
                }
                if (c > (numberOfBackgroundPoints - 1)) {
                    break;
                }
                xb[c] = xStart;
                yb[c] = yStart;
                c++;
            }
        }
        
        float[] xbe = new float[numberOfBackgroundPoints];
        float[] ybe = new float[numberOfBackgroundPoints];
        for (int i = 0; i < numberOfBackgroundPoints; i++) {
            // simulate x error as a percent error of 0.03 for each bin
            xbe[i] = xb[i] * 0.03f;
            ybe[i] = (float) (Math.sqrt(yb[i]));
        }
        
        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(xb, yb, xbe, ybe, xb.length);
        
        CompleteSamplingVoidFinder voidFinder = new CompleteSamplingVoidFinder();
        voidFinder.setSampling(VoidSampling.COMPLETE);
        
        voidFinder.findVoids(indexer);
        
        float[] linearDensities = voidFinder.getTwoPointDensities();
        
        assertNotNull(linearDensities);
        
        assertTrue(linearDensities.length == 20);
        
        // count values of 2 and 1.6817929
        int count0 = 0;
        int count1 = 0;
        
        for (int i = 0; i < linearDensities.length; i++) {
            float d = linearDensities[i];
            if (d == 2.) {
                count0++;
            } else if (Math.abs(d - 1.6817929) < 0.0001) {
                count1++;
            }
        }
        
        assertTrue(count0 == 12);
      
        assertTrue(count1 == 8);
      
    }
}
