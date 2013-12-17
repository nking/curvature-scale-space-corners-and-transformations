package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION;
import algorithms.curves.GEVYFit;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
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

    public void est_Find_Clusters_Stats() throws Exception {

        log.info("test_Find_Clusters_Stats()");

        float xmin = 0;
        float xmax = 3;
        float ymin = 0;
        float ymax = 3;
        
        int numberOfBackgroundPoints = 9;
        
        /*
         * case 0:
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

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);

        long seed = System.currentTimeMillis();
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
     
        seed = 1386750505246l;

        sr.setSeed(seed);
        log.info("SEED=" + seed);

        int nSwitches = 2;

        int nIterPerBackground = 1;

        int m = nIterPerBackground * nSwitches;

        DoubleAxisIndexer indexer = null;

        int count = 0;
                                
        for (int ii = 0; ii < nIterPerBackground; ii++) { 
            
            float[] xb = new float[numberOfBackgroundPoints];
            float[] yb = new float[numberOfBackgroundPoints];

            int xyStartOffset = 0;
            
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
            float t0 = MiscMath.findMax(xb);
            float t1 = MiscMath.findMax(yb);
            double diag = Math.sqrt(2)*(divYSz*divYSz);
            double dens = 2./diag;
            System.out.println("grid diag=" + diag + " 2./diag=" + dens);
            
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
                            "  (%d %d) best k=%.4f sigma=%.4f mu=%.4f chiSqSum=%.6f chst=%.1f",
                            i, ii, bestFit.getK(), bestFit.getSigma(), bestFit.getMu(), bestFit.getChiSqSum(), bestFit.getChiSqStatistic()
                        );
                        if (debug) {
                            log.info(plotLabel + " findVoid sampling=" + stats.getSampling().name());
                        }
                    }
                    
                }
                
                /*
                System.out.println("nGroups = " + twoPtC.getNumberOfGroups());
                for (int j = 0; j < twoPtC.getNumberOfGroups(); j++) {
                    System.out.println("group " + j);
                    ArrayPair group = twoPtC.getGroup(j);
                    for (int jj = 0; jj < group.getX().length; jj++) {
                        System.out.println("   " + group.getX()[jj] + "," + group.getY()[jj]);
                    }
                }*/
                
                plotter.addPlot(twoPtC, plotLabel);
                //plotter.addPlotWithoutHull(twoPtC, plotLabel);
                plotter.writeFile();

                /*
                if (i == 1) {
                    
                    TwoPointVoidStats stats = (TwoPointVoidStats)twoPtC.backgroundStats;
                    
                    String voidDensFileName = "void_densities_001.txt";
                    
                    writeVoidDensitiesToTestResources(voidDensFileName, 
                        stats.voidFinder.getTwoPointDensities(), stats.voidFinder.getTwoPointDensityErrors());
                }
                */
                
                count++;
            }
        }
        
        log.info("\n start computing stats for all sets");

        count = 0;

        log.info("SEED=" + seed);
    }
    
    public void test_Find_Clusters_Stats_2() throws Exception {

        log.info("test_Find_Clusters_Stats_2()");

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

        int nSwitches = 1;

        int nIterPerBackground = 1;

        int m = nIterPerBackground * nSwitches;

        DoubleAxisIndexer indexer = null;

        int count = 0;
        
        int nClusters = 3;
        
        /* 
         *    case 0:    9000 background points in a 1000 x 1000 area
         *               grid divisions are 10.5 units each so the most frequent expected linear density is 2./10.5 = 0.19
         *               The peak of the GEV and hence the density is found to be 0.098.<=====?? histogram needs improving!!!
         *    case 1:    a single cluster of 100 points in a 5x5 area of the 1000 x 1000 area without background points.
         *               Expect that the most frequent separation is sqrt(100)/5 = 
         *    case 2:  
         *
         */        
        int numberOfBackgroundPoints = 9000;
        
        CLUSTER_SEPARATION clusterSeparation = CLUSTER_SEPARATION.LARGE;
        
        for (int ii = 0; ii < nIterPerBackground; ii++) { 
            
            float[] xb = new float[numberOfBackgroundPoints];
            float[] yb = new float[numberOfBackgroundPoints];

            int xyStartOffset = 0;
            
            /*
            createRandomPointsInRectangle(sr, numberOfBackgroundPoints,
                xmin, xmax, ymin, ymax, xb, yb, xyStartOffset);
            */
            
            /* 
             *           |  10 x 10 and 9000 points.  94.9 points x dim in 10 cells = 1 pt/cell
             *           |
             *           |
             *    --------
             */
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
            float t0 = MiscMath.findMax(xb);
            float t1 = MiscMath.findMax(yb);
            double dens = 2./divYSz;
            System.out.println("grid division=" + divYSz + " 2./div=" + dens);
            
            float[] xbe = new float[numberOfBackgroundPoints];
            float[] ybe = new float[numberOfBackgroundPoints];
            for (int i = 0; i < numberOfBackgroundPoints; i++) {
                // simulate x error as a percent error of 0.03 for each bin
                xbe[i] = xb[i] * 0.03f;
                ybe[i] = (float) (Math.sqrt(yb[i]));
            }
            
            float[] xclust = null;
            float[] yclust = null;
            float[] xcluste = null;
            float[] ycluste = null;
            
            for (int i = 0; i < nSwitches; i++) {              
                                
                switch(i) {
                    
                    case 0: {
                        indexer = new DoubleAxisIndexer();
                        indexer.sortAndIndexXThenY(xb, yb, xbe, ybe, xbe.length);
                        break;
                    }
                    
                    case 1: {
                        int[] clusterNumbers = new int[]{100};
                        
                        int tot = clusterNumbers[0];
                        
                        xclust = new float[tot];
                        yclust = new float[tot];
                        
                        float[] xbc = new float[clusterNumbers.length];
                        float[] ybc = new float[clusterNumbers.length];
                        
                        generator.createRandomClusters(sr, 100, 105, 100, 105,
                            clusterNumbers, CLUSTER_SEPARATION.SMALL, xclust, yclust, xbc, ybc, 0);
                        
                        xcluste = new float[tot];
                        ycluste = new float[tot];
                        for (int j = xcluste.length; j < tot; j++) {
                            // simulate x error as a percent error of 0.03 for each bin
                            xcluste[j] = xclust[j] * 0.03f;
                            ycluste[j] = (float) (Math.sqrt(yclust[j]));
                        }
                        
                        indexer = new DoubleAxisIndexer();
                        indexer.sortAndIndexXThenY(xclust, yclust, xcluste, ycluste, xcluste.length);
                        
                        break;
                    }
                    
                    case 2: {
                        int[] clusterNumbers = new int[]{100};
                        
                        int tot = numberOfBackgroundPoints + 100;
                        
                        xb = Arrays.copyOf(xb, tot);
                        yb = Arrays.copyOf(yb, tot);
                        
                        float[] xbc = new float[clusterNumbers.length];
                        float[] ybc = new float[clusterNumbers.length];
                        
                        generator.createRandomClusters(sr, 100, 105, 100, 105,
                            clusterNumbers, CLUSTER_SEPARATION.SMALL, xb, yb,xbc, ybc, numberOfBackgroundPoints);
                        
                        xbe = Arrays.copyOf(xbe, tot);
                        ybe = Arrays.copyOf(ybe, tot);
                        for (int j = numberOfBackgroundPoints; j < tot; j++) {
                            // simulate x error as a percent error of 0.03 for each bin
                            xbe[j] = xb[j] * 0.03f;
                            ybe[j] = (float) (Math.sqrt(yb[j]));
                        }
                        
                        indexer = new DoubleAxisIndexer();
                        indexer.sortAndIndexXThenY(xb, yb, xbe, ybe, xbe.length);
                        
                        break;
                    }
                    
                    case 3: {
                        int[] clusterNumbers = new int[]{1000, 300, 100};
                        
                        int tot = numberOfBackgroundPoints + 1000 + 300 + 100;
                        
                        xb = Arrays.copyOf(xb, tot);
                        yb = Arrays.copyOf(yb, tot);
                        
                        float[] xbc = new float[clusterNumbers.length];
                        float[] ybc = new float[clusterNumbers.length];
                        
                        generator.createRandomClusters(sr, xmin, xmax, ymin, ymax,
                            clusterNumbers, clusterSeparation, xb, yb,xbc, ybc, numberOfBackgroundPoints);
                        
                        xbe = Arrays.copyOf(xbe, tot);
                        ybe = Arrays.copyOf(ybe, tot);
                        for (int j = numberOfBackgroundPoints; j < tot; j++) {
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
               
                /* to zoom in to confirm density estimate visually:
                DoubleAxisIndexer tmp = new DoubleAxisIndexer();
                tmp.sortAndIndexXThenY(
                    Arrays.copyOf(indexer.getX(), 100), Arrays.copyOf(indexer.getY(), 100),
                    Arrays.copyOf(indexer.getXErrors(), 100), Arrays.copyOf(indexer.getYErrors(), 100),
                    100);
                indexer = tmp;
                */
                
                log.info(" " + count + " (" + indexer.nXY + " points) ... ");


                if (writeToTmpData) {
                    // write to tmpdata if need to use in tests improve fits, histogram etc
                    String str = String.valueOf(count);
                    while (str.length() < 3) {
                        str = "0" + str;
                    }
                    String fileNamePostfix = "_clusters_" + str + ".dat";
                    String fileName = CreateClusterDataTest.indexerFileNamePrefix + fileNamePostfix;
                    String filePath = ResourceFinder.getAFilePathInTmpData(fileName);
                    CreateClusterDataTest.writeIndexer(filePath, indexer);
                }

                log.info(" " + count + " (" + indexer.nXY + " points) ... ");

                TwoPointCorrelation twoPtC = new TwoPointCorrelation(indexer);

                twoPtC.setDebug(true);
                
                twoPtC.useFindMethodForDataWithBackgroundPoints();
                
                //twoPtC.setAllowRefinement();
                
//twoPtC.setUseDownhillSimplexHistogramFitting();
              
                twoPtC.logPerformanceMetrics();
                //twoPtC.setBackground(0.35f, 0.02f);
                //twoPtC.setBackground(0.5f, 0.02f);
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
                            "  (%d %d) best k=%.4f sigma=%.4f mu=%.4f chiSqSum=%.6f chst=%.1f",
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

                /*
                if (i == 2) {
                    
                    TwoPointVoidStats stats = (TwoPointVoidStats)twoPtC.backgroundStats;
                    
                    String voidDensFileName = "void_densities_002.txt";
                    
                    writeVoidDensitiesToTestResources(voidDensFileName, 
                        stats.voidFinder.getTwoPointDensities(), stats.voidFinder.getTwoPointDensityErrors());
                }*/

                count++;
            }
        }
        
        log.info("\n start computing stats for all sets");

        count = 0;

        log.info("SEED=" + seed);
    }
}
