package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION;
import algorithms.curves.GEVYFit;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.misc.Statistic;
import algorithms.util.ArrayPair;
import algorithms.util.ResourceFinder;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.logging.Logger;

/**
 * class to test the TwoPointCorrelation class on larger datasets with and without clusters.
 * 
 * @author nichole
 */
public class FindClusters2Test extends BaseTwoPointTest {

    boolean debug = true;

    boolean writeToTmpData = false;

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public void test_Find_Clusters_Stats() throws Exception {

        log.info("test_Find_Clusters_Stats()");

        /* case 0:
               sampling the grid, there are 12 void separations of size 1
                                             8 void separations of size 1.414 
               
               the linear densities are then 2./1 = 2.0 is the most frequent
                                             2./1.414 = 1.414 is the next
                                             
               the peak of the GEV is found to be 1.71 which is what we'd expect
                              
           case 1:
               add 4 points closer than 1 to create a small group:
            
           3  |
              |
           2  *   *   *
              |   
           1  *   *   *           add them around (1,1) to lower left: delta of 0.01  (1-delta, 1), (1-delta, 1-delta), (1, 1-delta), (1-delta*2, 1-delta*2)
              |     
              *---*---*-----
           0  0   1   2   3  
           
              The background density is found as 1.77.  critical density = 2.5*1.77 = 4.425.
                            
         */

        long seed = System.currentTimeMillis();
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
     
        seed = 1386750505246l;

        sr.setSeed(seed);
        log.info("SEED=" + seed);

        
        float xmin = 0;
        float xmax = 3;
        float ymin = 0;
        float ymax = 3;
        
        int numberOfBackgroundPoints = 9;

        int nSwitches = 2;

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);
        
        DoubleAxisIndexer indexer = null;

        int count = 0;
                                
        for (int ii = 0; ii < 1; ii++) { 
            
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
            double tpdiv = 1./divYSz; // divYSz is the separation of 2 points
            double expectedDensity = Math.sqrt(numberOfBackgroundPoints)/(xmax-xmin);
            log.info("grid division=" + divYSz + "   and  density from grid size=" + tpdiv);
            
            float[] xbe = new float[numberOfBackgroundPoints];
            float[] ybe = new float[numberOfBackgroundPoints];
            for (int i = 0; i < numberOfBackgroundPoints; i++) {
                // simulate x error as a percent error of 0.03 for each bin
                xbe[i] = xb[i] * 0.03f;
                ybe[i] = (float) (Math.sqrt(yb[i]));
            }
            
            for (int i = 0; i < nSwitches; i++) {              
                                
                switch(i) {
                    
                    case 0: {
                        indexer = new DoubleAxisIndexer();
                        indexer.sortAndIndexXThenY(xb, yb, xbe, ybe, xbe.length);
                        break;
                    }

                    case 1: {
                        xb = Arrays.copyOf(xb, xb.length + 7);
                        yb = Arrays.copyOf(yb, yb.length + 7);
                        
                        float delta = 0.15f;
                        xb[numberOfBackgroundPoints] = 1.0f - delta;
                        yb[numberOfBackgroundPoints] = 1.0f;
                        
                        xb[numberOfBackgroundPoints + 1] = 1.0f - delta;
                        yb[numberOfBackgroundPoints + 1] = 1.0f - delta;
                        
                        xb[numberOfBackgroundPoints + 2] = 1.0f;
                        yb[numberOfBackgroundPoints + 2] = 1.0f - delta;
                        
                        xb[numberOfBackgroundPoints + 3] = 1.0f - 2.0f*delta;
                        yb[numberOfBackgroundPoints + 3] = 1.0f - 2.0f*delta;
                        
                        xb[numberOfBackgroundPoints + 4] = 1.0f - 2.0f*delta;
                        yb[numberOfBackgroundPoints + 4] = 1.0f - 1.0f*delta;
                        
                        xb[numberOfBackgroundPoints + 5] = 1.0f - 1.0f*delta;
                        yb[numberOfBackgroundPoints + 5] = 1.0f - 2.0f*delta;
                        
                        xb[numberOfBackgroundPoints + 6] = 1.0f - 2.0f*delta;
                        yb[numberOfBackgroundPoints + 6] = 1.0f;
                        
                        xbe = new float[xb.length];
                        ybe = new float[xb.length];
                        for (int j = 0; j < xb.length; j++) {
                            // simulate x error as a percent error of 0.03 for each bin
                            xbe[j] = xb[j] * 0.03f;
                            ybe[j] = (float) (Math.sqrt(yb[j]));
                        }
                        
                        indexer = new DoubleAxisIndexer();
                        indexer.sortAndIndexXThenY(xb, yb, xbe, ybe, xbe.length);
                        break;
                    }
                                        
                    default:
                        break;
                }

                log.info(" " + count + " (" + indexer.nXY + " points) ... ");

                TwoPointCorrelation twoPtC = new TwoPointCorrelation(indexer);

                twoPtC.setDebug(true);
                twoPtC.useFindMethodForDataWithBackgroundPoints();
                
                twoPtC.logPerformanceMetrics();
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
                            "(%d %d) k=%.4f s=%.4f m=%.4f chSq=%.6f chst=%.1f",
                            i, ii, bestFit.getK(), bestFit.getSigma(), bestFit.getMu(), bestFit.getChiSqSum(), bestFit.getChiSqStatistic()
                        );
                        if (debug) {
                            log.info(plotLabel + " findVoid sampling=" + stats.getSampling().name());
                        }
                    }
                }
                
                plotter.addPlot(twoPtC, plotLabel);
                //plotter.addPlotWithoutHull(twoPtC, plotLabel);
                plotter.writeFile();

                log.info(
                    "expected density = " + expectedDensity + "  calc density = " + twoPtC.getBackgroundSurfaceDensity()
                    + " npoints=" + numberOfBackgroundPoints + " xmax=" + xmax + "  (2/griddiv) = " + tpdiv
                    + "  r=exp/calc=" + (expectedDensity/twoPtC.getBackgroundSurfaceDensity()));
                // is expectedDensity always <= found surface density?
                // assertTrue( Math.abs(expectedDensity - twoPtC.getBackgroundSurfaceDensity()) < 0.25*expectedDensity);
                
                count++;
            }
        }
        
        log.info("\n start computing stats for all sets");

        count = 0;

        log.info("SEED=" + seed);
    }
    
}
