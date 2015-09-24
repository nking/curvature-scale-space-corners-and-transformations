package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.clustering.distanceTransform.DistanceTransform;
import algorithms.curves.GEVYFit;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import java.security.SecureRandom;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;

/**
 * class to test the TwoPointCorrelation class on larger datasets with and without clusters.
 * 
 * @author nichole
 */
public class FindClusters5Test extends BaseTwoPointTest {

    boolean debug = true;

    boolean writeToTmpData = false;

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public void testFindClustersStats() throws Exception {

        log.info("testFindClustersStats()");
        
        /* Goal of this test is to examine the substructure created by increasing numbers of randomly 
           placed points.
           
           1000 x 1000 unit^2 space to place
           
               N=100,450,900,4500,9000,14500,19000,24500,29000
               random points                
               
           And track, N, linear density, and the number of groups found.
           
                  N=(100)  calcLinDens=(0.1236)  expectedLinDens=(0.0100)  nGroups=(0)
                  N=(450)  calcLinDens=(0.1489)  expectedLinDens=(0.0212)  nGroups=(0)
                  N=(900)  calcLinDens=(0.1564)  expectedLinDens=(0.0300)  nGroups=(0)
                  N=(4500)  calcLinDens=(0.1580)  expectedLinDens=(0.0671)  nGroups=(125)
                  N=(9000)  calcLinDens=(0.1536)  expectedLinDens=(0.0949)  nGroups=(665)
                  N=(14500)  calcLinDens=(0.1716)  expectedLinDens=(0.1204)  nGroups=(1391)
                  N=(19000)  calcLinDens=(0.1747)  expectedLinDens=(0.1378)  nGroups=(2160)
                  N=(24500)  calcLinDens=(0.1553)  expectedLinDens=(0.1565)  nGroups=(2922)
                  N=(29000)  calcLinDens=(0.1544)  expectedLinDens=(0.1703)  nGroups=(2980)
               SEED=1386750505246
               
            Can see that after N=1000, begin to see groups for from randomly close points.
                roughly nGroups is less than or equal to (0.2*N)
         */
        
        int[] numberOfBackgroundPoints = new int[]{100,450,900,4500,9000,14500,19000,24500,29000};
        
        int[] nGroupsFound = new int[numberOfBackgroundPoints.length];
        float[] expectedLinearDensities = new float[nGroupsFound.length];
        float[] calcLinearDensities = new float[nGroupsFound.length];
            
        float xmin = 0;
        float xmax = 1000;
        float ymin = 0;
        float ymax = 1000;

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);

        long seed = System.currentTimeMillis();
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
     
        //seed = 1443059507225L;

        sr.setSeed(seed);
        log.info("SEED=" + seed);

        AxisIndexer indexer = null;
        
        for (int ii = 0; ii < numberOfBackgroundPoints.length; ii++) { 
            
            int xyStartOffset = 0;
            
            float[] xb = new float[numberOfBackgroundPoints[ii]];
            float[] yb = new float[numberOfBackgroundPoints[ii]];
            
            double expectedDensity = Math.sqrt(numberOfBackgroundPoints[ii])/(xmax-xmin);

            createRandomPointsInRectangle(sr, numberOfBackgroundPoints[ii],
                xmin, xmax, ymin, ymax, xb, yb, xyStartOffset);
                    
            float[] xbe = new float[numberOfBackgroundPoints[ii]];
            float[] ybe = new float[numberOfBackgroundPoints[ii]];
            for (int i = 0; i < numberOfBackgroundPoints[ii]; i++) {
                // simulate x error as a percent error of 0.03 for each bin
                xbe[i] = xb[i] * 0.03f;
                ybe[i] = (float) (Math.sqrt(yb[i]));
            }
            
            indexer = new AxisIndexer();
            
            indexer.sortAndIndexX(xb, yb, xbe, ybe, xbe.length);
                        
            log.info(" " + ii + " (" + indexer.nXY + " points) ... ");

            Set<PairInt> points = new HashSet<PairInt>();
            for (int k = 0; k < indexer.getNXY(); ++k) {
                PairInt p = new PairInt(Math.round(indexer.getX()[k]),
                    Math.round(indexer.getY()[k]));
                points.add(p);
            }
            
            TwoPointCorrelation twoPtC = new TwoPointCorrelation(indexer);

            //twoPtC.setDebug(true);
            
            //twoPtC.setAllowRefinement();
                
//twoPtC.setUseDownhillSimplexHistogramFitting();
              
            //twoPtC.logPerformanceMetrics();
            
            //twoPtC.setBackground(0.5f, 0.003f);
            twoPtC.calculateBackground();
                
            twoPtC.findClusters();

            String plotLabel = "";
                
            if (twoPtC.backgroundStats != null) {
                TwoPointVoidStats stats = (TwoPointVoidStats)twoPtC.backgroundStats;
                HistogramHolder histogram = stats.statsHistogram;
    
                GEVYFit bestFit = stats.bestFit;
                if (bestFit != null) {
                    
                    // label needs:  x10, peak,  mean/peak, median/mean and x80/median
                    plotLabel = String.format(
                        "(%d) k=%.4f s=%.4f m=%.4f chiSq=%.6f chst=%.1f",
                        numberOfBackgroundPoints[ii], bestFit.getK(), bestFit.getSigma(), bestFit.getMu(), bestFit.getChiSqSum(), bestFit.getChiSqStatistic()
                    );
                    if (debug) {
                        log.info(plotLabel + " findVoid sampling=" + stats.getSampling().name());
                    }
                }
            }
            
            plotter.addPlot(twoPtC, plotLabel);
            //plotter.addPlotWithoutHull(twoPtC, plotLabel);
            plotter.writeFile();
            
            // ----- a look at the distance transform ----
            int[] minMaxXY = MiscMath.findMinMaxXY(points);
            int w = minMaxXY[1] + 1;
            int h = minMaxXY[3] + 1;
            DistanceTransform dtr = new DistanceTransform();
            int[][] dt = dtr.applyMeijsterEtAl(points, w, h);
            float[] values = new float[dt.length*dt[0].length];
            int count2 = 0;
            for (int i0 = 0; i0 < dt.length; ++i0) {
                for (int j0 = 0; j0 < dt[0].length; ++j0) {
                    int v = dt[i0][j0];
                    values[count2] = (float)(1./Math.sqrt(v));
                    count2++;
                }
            }
            float xl = 0.1f;
            /*if (i==8) {
                xl = 0.05f;
            }*/
            float[] vErrors = Errors.populateYErrorsBySqrt(values);
            HistogramHolder hist = Histogram.createSimpleHistogram(
                0, xl, 40, values, vErrors);
            hist.plotHistogram("clstr", "_cluster_" + ii);
            
            
            nGroupsFound[ii] = twoPtC.getNumberOfGroups();
            expectedLinearDensities[ii] = (float)expectedDensity;
            calcLinearDensities[ii] = twoPtC.getBackgroundDensity();
                        
        }
        
        log.info("\n start computing stats for all sets");
        
        for (int i = 0; i < numberOfBackgroundPoints.length; i++) {
            
            String str = String.format("N=(%d)  calcLinDens=(%.4f)  expectedLinDens=(%.4f)  nGroups=(%d)",
                numberOfBackgroundPoints[i], calcLinearDensities[i], expectedLinearDensities[i], nGroupsFound[i]);
            
            log.info(str.toString());
            
            // assertion that shouldn't change too much w/ other component improvements
            assertTrue(nGroupsFound[i] <= 0.2 * numberOfBackgroundPoints[i]);
        }

        log.info("SEED=" + seed);
    }
}
