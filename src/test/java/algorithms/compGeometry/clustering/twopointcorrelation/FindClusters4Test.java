package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION;
import algorithms.curves.GEVYFit;
import algorithms.misc.HistogramHolder;
import algorithms.util.ArrayPair;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.logging.Logger;

/**
 * class to test the TwoPointCorrelation class on larger datasets with and without clusters.
 * 
 * @author nichole
 */
public class FindClusters4Test extends BaseTwoPointTest {

    boolean debug = true;

    boolean writeToTmpData = false;

    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public void testFindClustersStats() throws Exception {

        log.info("testFindClustersStats()");
        
        /* 
         *    case 0:    9000 background points in a 1000 x 1000 area
         *               grid divisions are 10.5 units each so the most frequent expected linear density is 2./10.5 = 0.19
         *               The peak of the GEV and hence the density is found to be near 0.1. 
         *               
         *    case 1:    a single cluster of 100 points in a 5x5 area of the 1000 x 1000 area without background points.
         *               Expect that the most frequent separation is sqrt(100)/5 = 2.0.  ld = 2./2 = 1.
         *               The density is found to be 1.6
         *               
         *    case 2:  
         *               The 9000 background points in a 1000 x 1000 area with a single cluster of points in a 5 x 5 area.
         *               Expect the most frequent separation to be near 0.1 still.
         *               The peak of the GEV and hence the density is found to be 0.12.
         *               
         *    case 3:    The 9000 background points in a 1000 x 1000 area with 3 large clusters of varying density.
         *    
         *    case 4:    9000 random points in a 1000 x 1000 area.
         *    
         *               The rest here is reasoning about whether the peak found and the clusters implied by 
         *               2.5* that density are the criteria to apply to find clusters.
         *    
         *               Can see from the counting in DoubleAxisInderStats that the avg is 111.111 in a cell of size 111.1
         *               and that the st.dev. is 10.5 and that all cell counts are < avg + 2.5*st.dev.
         *               so the density from those stats is 0.1.  the peak of the GEV is near 0.1 too, it's at 0.089.
         *               
         *               The random placement makes shorter separations more frequent than in the grid data
         *               which can be seen as the peak at 0.1 while the grid data peak is at 1.6.
         *                              
         *               background stats:
         *               
         *                   In TwoPointVoidStats, some basic stats on cell counts are performed
         *                         Statistic statistic = stats.calculateCellDensities(nCellsPerDimensionForStats, indexer)
         *                   The average counts in a cell of size 111.11 is 111.11 and stdev is 10.5 which is sqrt(111.11).
         *                   that says that this random background is evenly distributed  within an error of sqrt(counts).
         *                   
         *                   those counts convert to linear densities of Math.sqrt(111.11)/(111.11) = 0.095 points/unit dimension
         *                   
         *                   structure smaller than 111.11 units of space by a factor of 2.5 might be visible
         *                   with simple stats using cell size of 44.
         *                   Using the same stats.calculateCellDensities(25, indexer) for smaller cells,
         *                   One can see that for cells of size 40 on one side, the avg counts are 14.4 with a stdev of 3.9 which is sqrt(14.4).
         *                   
         *                   If one uses the program to place points into a group when their separation implies a linear
         *                   density that is smaller than a critical limit:
         *                   
         *                   One can see that if we use the background density as an error estimate and look for points above
         *                   that threshold by a factor of 2.5 we get a number which = 0.25 for a critical linear density.
         *                   Using that critical density finds many clusters in the data, maybe more so than would be expected.
         *                   
         *                   If however, one uses a critical linear density of 0.475 or so, there are no clusters found more dense than that limit.
         *                   
         *                   If one applies the density found from the grid dataset which has same features, but is in a grid of
         *                   evenly spaced intervals rather than randomly distributed, that density is 0.125.
         *                   For that, the number of clusters matches what one might expect by eye.
         *                   
         *                   So far, one can say that randomly placing points within a space eventually introduces substructure
         *                   that is clustered more than if the same points were in grid formation as a contrast (grid spacing is even
         *                   and so the voids are always the edges or diagonals of a square.  there's no 'clumping' on scales finer
         *                   than those voids and those are the peak of the linear density distribution).
         *                   
         *                   Below, in test test_Find_Clusters_Stats_3(), one can see the effects of substructure introduced
         *                   by random placement of points within an area.
         *                   
         *                   
         */
        
        int numberOfBackgroundPoints = 9000;
            
        float xmin = 0;
        float xmax = 1000;
        float ymin = 0;
        float ymax = 1000;

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);

        long seed = System.currentTimeMillis();
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
     
        //seed = 1387770646833l;

        sr.setSeed(seed);
        log.info("SEED=" + seed);

        int nSwitches = 4;//6;

        int nIterPerBackground = 1;

        AxisIndexer indexer = null;

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
            
            float[] xclust = null;
            float[] yclust = null;
            float[] xcluste = null;
            float[] ycluste = null;
            
            float[] xg5 = null;
            float[] yg5 = null;
            
            for (int i = 0; i < nSwitches; i++) {
                                
                switch(i) {
                    
                    case 0: {
                        indexer = new AxisIndexer();
                        indexer.sortAndIndexX(xb, yb, xbe, ybe, xbe.length);
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
                        for (int j = 0; j < tot; j++) {
                            // simulate x error as a percent error of 0.03 for each bin
                            xcluste[j] = xclust[j] * 0.03f;
                            ycluste[j] = (float) (Math.sqrt(yclust[j]));
                        }
                        
                        indexer = new AxisIndexer();
                        indexer.sortAndIndexX(xclust, yclust, xcluste, ycluste, xcluste.length);
                        
                        break;
                    }
                    
                    case 2: {
                        int[] clusterNumbers = new int[]{100};
                        
                        int tot = numberOfBackgroundPoints + clusterNumbers[0];
                        
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
                        
                        indexer = new AxisIndexer();
                        indexer.sortAndIndexX(xb, yb, xbe, ybe, xbe.length);
                        
                        break;
                    }
                    
                    case 3: {
                        int[] clusterNumbers = new int[]{1000};
                        
                        int tot = numberOfBackgroundPoints + 1000;
                        
                        xb = Arrays.copyOf(xb, tot);
                        yb = Arrays.copyOf(yb, tot);
                                                
                        generator.createRandomPointsAroundCenter(sr, 50, clusterNumbers[0], 400, 400, xb, yb, numberOfBackgroundPoints);
                        
                        xbe = Arrays.copyOf(xbe, tot);
                        ybe = Arrays.copyOf(ybe, tot);
                        for (int j = numberOfBackgroundPoints; j < tot; j++) {
                            // simulate x error as a percent error of 0.03 for each bin
                            xbe[j] = xb[j] * 0.03f;
                            ybe[j] = (float) (Math.sqrt(yb[j]));
                        }
                        
                        indexer = new AxisIndexer();
                        indexer.sortAndIndexX(xb, yb, xbe, ybe, xbe.length);
                        
                        break;
                    }
                    
                    //   ======================= make the background random instead of a grid =======
                    
                    case 4: {
                        
                        xb = new float[numberOfBackgroundPoints];
                        yb = new float[numberOfBackgroundPoints];

                        xyStartOffset = 0;
                        
                        createRandomPointsInRectangle(sr, numberOfBackgroundPoints,
                            xmin, xmax, ymin, ymax, xb, yb, xyStartOffset);
                    
                        xbe = new float[numberOfBackgroundPoints];
                        ybe = new float[numberOfBackgroundPoints];
                        for (int j = 0; j < numberOfBackgroundPoints; j++) {
                            // simulate x error as a percent error of 0.03 for each bin
                            xbe[j] = xb[j] * 0.03f;
                            ybe[j] = (float) (Math.sqrt(yb[j]));
                        }
                        
                        indexer = new AxisIndexer();
                        indexer.sortAndIndexX(xb, yb, xbe, ybe, xbe.length);
                        
                        break;
                    }
                    
                    case 5: {
                        // add to existing background
                        
                        int[] clusterNumbers = new int[]{1000};
                        
                        int tot = numberOfBackgroundPoints + clusterNumbers[0];
                        
                        xb = Arrays.copyOf(xb, tot);
                        yb = Arrays.copyOf(yb, tot);
                        
                        generator.createRandomPointsAroundCenter(sr, 50, clusterNumbers[0], 400, 400, xb, yb, numberOfBackgroundPoints);
                        
                        xbe = Arrays.copyOf(xbe, tot);
                        ybe = Arrays.copyOf(ybe, tot);
                        int ci = 0;
                        xg5 = new float[clusterNumbers[0]];
                        yg5 = new float[clusterNumbers[0]];
                        for (int j = numberOfBackgroundPoints; j < tot; j++) {
                            // simulate x error as a percent error of 0.03 for each bin
                            xbe[j] = xb[j] * 0.03f;
                            ybe[j] = (float) (Math.sqrt(yb[j]));
                            
                            // for test assertion, store group
                            xg5[ci] = xb[j];
                            yg5[ci] = yb[j];
                            ci++;
                        }
                        
                        indexer = new AxisIndexer();
                        indexer.sortAndIndexX(xb, yb, xbe, ybe, xbe.length);
                        break;
                    }
                    
                    default:
                        break;
                }
               
                log.info(" " + count + " (" + indexer.nXY + " points) ... ");

                TwoPointCorrelation twoPtC = new TwoPointCorrelation(indexer);

                //twoPtC.setDebug(true);
                
                if (i == 1) {
                    twoPtC.useFindMethodForDataWithBackgroundPoints();
                }
                
                //twoPtC.setAllowRefinement();
                
//twoPtC.setUseDownhillSimplexHistogramFitting();
              
                twoPtC.logPerformanceMetrics();
                
                //twoPtC.setBackground(0.13f, 0.015f);
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
            
                if (i == 2) {  
                    assertTrue(twoPtC.getNumberOfGroups() == 1);
                    /*
                    log.info("nGroups = " + twoPtC.getNumberOfGroups());
                    for (int j = 0; j < twoPtC.getNumberOfGroups(); j++) {
                        log.info("group " + j);
                        ArrayPair group = twoPtC.getGroup(j);
                        for (int jj = 0; jj < group.getX().length; jj++) {
                            log.info("   " + group.getX()[jj] + "," + group.getY()[jj]);
                        }
                    }*/
                    int n0 = twoPtC.getGroup(0).getX().length;
                    assertTrue(n0 >= 100 && (n0 <= 105));
                } else if (i == 3) {
                    //assertTrue(twoPtC.getNumberOfGroups() == 2);
                    ArrayPair hull;
                    if (twoPtC.getGroup(0).getX().length > twoPtC.getGroup(1).getX().length) {
                        hull = twoPtC.getGroupHull(0);
                    } else {
                        hull = twoPtC.getGroupHull(1);
                    }
                    float[] areaAndCenter = twoPtC.calculateAreaAndCentroidOfHull(hull.getX(), hull.getY());
                    assertNotNull(areaAndCenter);
                    float radius = (float) Math.sqrt(areaAndCenter[0]/(2.*Math.PI));
                    //assertTrue(radius <= 50.f);
                } else if (i == 4 || ii == 5) {
                    // assertion that shouldn't change too much w/ other component improvements
                    assertTrue(twoPtC.getNumberOfGroups() <= 0.2 * indexer.getNumberOfPoints());    
                }
                
                if (i == 5) {
                    log.info("assert case 5 results");
                    // assert that the points added as a group were found as a group:
                    int largestGroupId = -1;
                    int maxN = Integer.MIN_VALUE;
                    for (int j = 0; j < twoPtC.getNumberOfGroups(); j++) {
                        if (twoPtC.getGroup(j).getX().length > maxN) {
                            maxN = twoPtC.getGroup(j).getX().length;
                            largestGroupId = j;
                        }
                    }
                    ArrayPair group = twoPtC.getGroup(largestGroupId);
                    boolean[] found = new boolean[xg5.length];
                    for (int j = 0; j < xg5.length; j++) {
                        float xx = xg5[j];
                        float yy = yg5[j];
                        for (int jj = 0; jj < group.getX().length; jj++) {
                            if ((Math.abs(group.getX()[jj] - xx) < 0.0001f) && (Math.abs(group.getY()[jj] - yy)) < 0.0001f) {
                                found[j] = true;
                                break;
                            }
                        }
                    }
                    int cf = 0;
                    for (int j = 0; j < found.length; j++) {
                        if (found[j]) {
                            cf++;
                        } else {
                            log.info("did not find point (" + xg5[j] + "," + yg5[j] + ")");
                        }
                    }
                    log.info("====> expected " + found.length + " and found=" + cf + " and group size=" + group.getX().length);
                    // Missing some of the points if use the calculated background density of 0.16655745 f.
                    //     expected 1000 and found=976 and group size=1039.  so recovered 98% of the large group. <====
                    // if use a manual density of 0.13f we recover each point of the group plus those randomly placed near:
                    //     expected 1000 and found=1000 and group size=1082
                    // Note that the smaller 100 point group was not recovered.
                }
               
                if (false && ((i == 0) || (i == 1))) {
                    
                    TwoPointVoidStats stats = (TwoPointVoidStats)twoPtC.backgroundStats;
                    
                    String voidDensFileName = (i == 0) ? "void_densities_002.txt" : "void_densities_003.txt";
                    
                    writeVoidDensitiesToTestResources(voidDensFileName, 
                        stats.voidFinder.getTwoPointDensities(), stats.voidFinder.getTwoPointDensityErrors());
                }

                log.info(
                    "expected density = " + expectedDensity + "  calc density = " + twoPtC.getBackgroundDensity()
                    + " npoints=" + numberOfBackgroundPoints + " xmax=" + xmax + "  (2/griddiv) = " + tpdiv
                    + "  r=exp/calc=" + (expectedDensity/twoPtC.getBackgroundDensity()));
                
                
                if (true) {
                    if (twoPtC.backgroundStats != null && twoPtC.backgroundStats instanceof TwoPointVoidStats) {
                        HistogramHolder histogram = ((TwoPointVoidStats)twoPtC.backgroundStats).statsHistogram;
                        if (histogram != null) {
                            StringBuilder xsb = new StringBuilder();
                            StringBuilder ysb = new StringBuilder();
                            StringBuilder xesb = new StringBuilder();
                            StringBuilder yesb = new StringBuilder();
    
                            for (int z = 0; z < histogram.getYHist().length; z++) {
                                if (z > 0) {
                                    xsb.append("f, ");
                                    ysb.append("f, ");
                                    xesb.append("f, ");
                                    yesb.append("f, ");
                                }
                                xsb.append(histogram.getXHist()[z]);
                                ysb.append(histogram.getYHist()[z]);
                                xesb.append(histogram.getXErrors()[z]);
                                yesb.append(histogram.getYErrors()[z]);
                            }
                            log.info("float[] x = new float[]{" + xsb.append("f").toString() + "};");
                            log.info("float[] y = new float[]{" + ysb.append("f").toString() + "};");
                            log.info("float[] xe = new float[]{" + xesb.append("f").toString() + "};");
                            log.info("float[] ye = new float[]{" + yesb.append("f").toString() + "};");
                        }
                    }
                }
                
                count++;
            }
        }
        
        log.info("\n start computing stats for all sets");

        count = 0;

        log.info("SEED=" + seed);
    }
    
}
