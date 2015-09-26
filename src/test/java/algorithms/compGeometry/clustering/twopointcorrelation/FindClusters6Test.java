package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.misc.HistogramHolder;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 * @author nichole
 */
public class FindClusters6Test extends TestCase {

    public FindClusters6Test(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    boolean debug = true;

    boolean persistHistogram = false;

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public void testFindClusters() throws Exception {

        log.info("testFindClusters()");
        
        String[] fileNames = {
            "Aggregation.txt", "Compound.txt", "Pathbased.txt" , "Spiral.txt",
            "D31.txt", "R15.txt" , "Jain.txt", "Flame.txt",
            "a1.txt", "a2.txt", "a3.txt"/*,
            "s1.txt", "s2.txt", "s3.txt", "s4.txt",
            "birch1.txt", "birch2.txt", "birch3.txt" */
        };
        
        // The 'a' datasets are strangely constructed.  
        // the critical threshhold is only held by single isolated pairs of points, 
        // never by 3 or more points 
        // (that is 2 or more pairs within critical distance of one another).
        // TO HANDLE these datasets properly with a histogram, one would need
        // to plot triplet or quadruplet point average separations instead of
        // pairs. OR, could use the useFindMethodForDataWithoutBackgroundPoints
        // which uses the smallest measured bin in the histogram as the 
        // representative bckground density, use findGroups to find the first
        // round of clusters, then take extra steps: make a
        // histogram of the centers of those clusters ==> the peak of that 
        // histogram is then the the critical separation for "associated" points.
        // one doesn't usually want to take these extra steps, but can see
        // that it's necessary for this unusually constructed dataset.
        // One would need to determine whether this step was necessary
        // or not from the start... it's similar to an aglomerative approach
        // to bottom up clustering, that is one more iteration to add
        // clusters together...
        
        
        int[] expectedNGroups = {
            7, 6, 3, 3, 
            31, 15 , 2, 2,
            20, 35, 50/*,
            15, 15, 15, 15,
            100, 100, 100*/
        };
        
        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter();
        
        for (int i = 0; i < fileNames.length; i++) {

            String fileName = fileNames[i];
            
            //NOTE:  for i=8, distance transform needs alot of memory for array size, so have divided numbers there by 10
            AxisIndexer indexer = CreateClusterDataTest.getUEFClusteringDataset(
                fileName);
            
            TwoPointCorrelation twoPtC = new TwoPointCorrelation(indexer);
            
            // s1.txt 0.0008 finds centers of each "cluster"
            // dataset is more like saturated ccd wells surrounding by wells 
            //  with too low exposure resulting in poor S/N
            //  --> use useFindMethodForDataWithoutBackgroundPoints
            //     could benefit from voronoi afterwards
            
            // birch1.txt is not currently fit by this code.
            //     how would it be better solved?
            //     binning (smoothing) of points to reduce number of points?
            //     a kernel technique?  can see from the regularity of the
            //     data that that should be a yes, but is the technique
            //     useful for other datasets too? 
            //        -- outline the answer if yes or no
            
            // birch2.txt is 0.001
            // birch3.txt is 0.0007, but not clearly delineated as clusters...
            
            //twoPtC.setBackground(0.32f, 0.0001f);

            //twoPtC.setAllowRefinement();

            if (fileName.startsWith("a") || fileName.equals("Flame.txt")) {
                twoPtC.useFindMethodForDataWithoutBackgroundPoints();
            }

            twoPtC.calculateBackground();
            
            twoPtC.findClusters();

            TwoPointVoidStats stats = (TwoPointVoidStats)twoPtC.backgroundStats;

            String plotLabel = "";

            if (twoPtC.backgroundStats != null && twoPtC.backgroundStats 
                instanceof TwoPointVoidStats 
                && ((TwoPointVoidStats)twoPtC.backgroundStats).bestFit != null) {

                // centroid of area defined by the top portion of the fit where y >= ypeak/2
                float[] areaAndXYTopCentroid = 
                    TwoPointVoidStats.calculateCentroidOfTop(
                    stats.bestFit.getOriginalScaleX(), 
                    stats.bestFit.getOriginalScaleYFit(), 0.5f);

                if (areaAndXYTopCentroid != null) {

                    plotLabel = String.format(
                        "  (%4d)  k=%.4f  sigma=%.4f  mu=%.4f chst=%.1f",
                        twoPtC.indexer.nXY, 
                        ((TwoPointVoidStats)twoPtC.backgroundStats).bestFit.getK(),
                        ((TwoPointVoidStats)twoPtC.backgroundStats).bestFit.getSigma(),
                        ((TwoPointVoidStats)twoPtC.backgroundStats).bestFit.getMu(),
                        areaAndXYTopCentroid[1], 
                        ((TwoPointVoidStats)twoPtC.backgroundStats).bestFit.getChiSqStatistic());

                } else {

                    plotLabel = String.format("  (%4d)", twoPtC.indexer.nXY);
                }
                
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
                    log.info("fileName=" + fileNames[i]);
                    log.info("float[] x = new float[]{" + xsb.append("f").toString() + "};");
                    log.info("float[] y = new float[]{" + ysb.append("f").toString() + "};");
                    log.info("float[] xe = new float[]{" + xesb.append("f").toString() + "};");
                    log.info("float[] ye = new float[]{" + yesb.append("f").toString() + "};");
                }
            }

            //xmin, xmax, ymin, ymax
            float[] minMaxes = indexer.findXYMinMax();

            Float xmin = minMaxes[0];
            Float xmax = minMaxes[1];
            Float ymin = minMaxes[2];
            Float ymax = minMaxes[3];
            
            System.out.println("fileName=" + fileNames[i] + " " 
                + xmin + " " + xmax + " " + ymin + " " + ymax);
            
            plotter.addPlotWithoutHull(twoPtC, plotLabel, xmin, xmax, ymin, ymax);
            plotter.writeFile();
            
            int nExpected = expectedNGroups[i];
            
            int nGroups = twoPtC.getNumberOfGroups();

            /*assertTrue(twoPtC.getNumberOfGroups() == nExpected);

            ArrayPair centroids = twoPtC.getHullCentroids();

            assertTrue(centroids.getX().length == nExpected);
            assertTrue(centroids.getY().length == nExpected);
            */
            
            log.info(twoPtC.indexer.nXY + " points ... ");
            
        }

        log.info("  END ");
    }
    
}
