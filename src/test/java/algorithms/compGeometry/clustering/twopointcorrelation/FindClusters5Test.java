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
public class FindClusters5Test extends BaseTwoPointTest {

    boolean debug = true;

    boolean writeToTmpData = false;

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public void test_Find_Clusters_Stats() throws Exception {

        log.info("test_Find_Clusters_Stats()");
        
        /* Goal of this test is to examine the substructure created by increasing numbers of randomly 
           placed points.
           
           1000 x 1000 unit^2 space to place
           
               N=100,450,900,4500,9000,14500,19000,24500,29000
               random points                
               
           And track, N, linear density, and the number of groups found.
           
              N    calculateLinearDensity    expectedLinearDensity   nGroups
              100  0.1241                    0.0100                    0
              450  0.0907                    0.0212                    2
              900  0.0945                    0.0300                   10
             4500  0.0918                    0.0671                  475
             9000  0.0892                    0.0949                  975
            14500  0.0911                    0.1204                  706
            19000  0.0918                    0.1378                  235
            24500  0.0904                    0.1565                   48
            29000  0.0911                    0.1703                   15
               SEED=1386750505246
               
            Can see that after N=9000, the number of clusters introduced is large and the groups
            implied by the clusters quickly saturates to essentially one large group by N=29000
            
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
     
        seed = 1386750505246l;

        sr.setSeed(seed);
        log.info("SEED=" + seed);

        DoubleAxisIndexer indexer = null;
        
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
            
            indexer = new DoubleAxisIndexer();
            
            indexer.sortAndIndexXThenY(xb, yb, xbe, ybe, xbe.length);
                        
            log.info(" " + ii + " (" + indexer.nXY + " points) ... ");

            TwoPointCorrelation twoPtC = new TwoPointCorrelation(indexer);

            twoPtC.setDebug(true);
            
            //twoPtC.setAllowRefinement();
                
//twoPtC.setUseDownhillSimplexHistogramFitting();
              
            twoPtC.logPerformanceMetrics();
                
            //twoPtC.setBackground(0.125f, 0.015f);
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
            
            nGroupsFound[ii] = twoPtC.getNumberOfGroups();
            expectedLinearDensities[ii] = (float)expectedDensity;
            calcLinearDensities[ii] = twoPtC.getBackgroundSurfaceDensity();
            
        }
        
        log.info("\n start computing stats for all sets");
        
        for (int i = 0; i < numberOfBackgroundPoints.length; i++) {
            
            String str = String.format("N=(%d)  calcLinDens=(%.4f)  expectedLinDens=(%.4f)  nGroups=(%d)",
                numberOfBackgroundPoints[i], calcLinearDensities[i], expectedLinearDensities[i], nGroupsFound[i]);
            
            log.info(str.toString());
            
        }

        log.info("SEED=" + seed);
    }
}
