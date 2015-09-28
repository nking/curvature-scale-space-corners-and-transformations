package com.climbwithyourfeet.clustering;

import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION;
import algorithms.compGeometry.clustering.twopointcorrelation.AxisIndexer;
import algorithms.compGeometry.clustering.twopointcorrelation.BaseTwoPointTest;
import algorithms.compGeometry.clustering.twopointcorrelation.CreateClusterDataTest;
import algorithms.misc.MiscMath;
import algorithms.util.ResourceFinder;
import com.climbwithyourfeet.clustering.util.PairInt;
import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import javax.imageio.ImageIO;

/**
 *
 * @author nichole
 */
public class DTClusterFinderTest extends BaseTwoPointTest {
    
    private Logger log = Logger.getLogger(this.getClass().getName());

    public void testFindRanGenClusters() throws Exception {
        
        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        long seed = System.currentTimeMillis();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");

        //seed = 1443413022277L;

        log.info("SEED=" + seed);
        
        sr.setSeed(seed);

        int nSwitches = 3;

        int nIterPerBackground = 3;

        AxisIndexer indexer = null;

        int count = 0;
        
        // these were generated w/ a specifc seed so may need more tolerance
        float[] r0s = new float[]{
            0.01f, 0.18f, 0.4f,
            0.01f, 0.18f, 0.4f,
            0.02f, 0.18f, 0.4f
        };
        float[] r1s = new float[]{
            0.05f, 0.285f, 0.5f,
            0.05f, 0.28f, 0.5f,
            0.04f, 0.28f, 0.5f
        };
        
        ClusterPlotter plotter = new ClusterPlotter();

        //TODO: improve these simulated clusters and assert the numbers
        
        for (int ii = 0; ii < nIterPerBackground; ii++) {
            for (int i = 0; i < nSwitches; i++) {
                switch(i) {
                    case 0:
                        //~100
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            //10, 100, 110, 100.0f);
                            3, 33, 33, 0.1f);
                        break;
                    case 1:
                        int[] clusterNumbers = new int[]{2000, 300, 1000};
                        int nBackgroundPoints = 10000;
                        CLUSTER_SEPARATION clusterSeparation = CLUSTER_SEPARATION.LARGE;

                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            clusterNumbers, nBackgroundPoints, clusterSeparation);
                        break;
                    default: {
                         // 100*100
                         indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 30, 60, 100.0f);
                         break;
                    }
                }
                
                log.info(" " + count + " (" + indexer.getNXY() + " points) ... ");

                Set<PairInt> points = new HashSet<PairInt>();
                for (int k = 0; k < indexer.getNXY(); ++k) {
                    PairInt p = new PairInt(Math.round(indexer.getX()[k]),
                        Math.round(indexer.getY()[k]));
                    points.add(p);
                }

                int[] minMaxXY = MiscMath.findMinMaxXY(points);
                int width = minMaxXY[1] + 1;
                int height = minMaxXY[3] + 1;
                    
                DTClusterFinder<PairInt> clusterFinder = 
                    new DTClusterFinder<PairInt>(points, width, height);
                
                clusterFinder.setToDebug();

                clusterFinder.calculateCriticalDensity();
                //clusterFinder.setCriticalDensity(0.139f);
                
                clusterFinder.findClusters();

                int nGroups = clusterFinder.getNumberOfClusters();

                List<Set<PairInt>> groupList = new ArrayList<Set<PairInt>>();
                for (int k = 0; k < nGroups; ++k) {
                    Set<PairInt> set = clusterFinder.getCluster(k);
                    groupList.add(set);
                }
                
                plotter.addPlotWithoutHull(
                    (int)Math.floor(minMaxXY[0] - 1), 
                    (int)Math.ceil(minMaxXY[1] + 1), 
                    (int)Math.floor(minMaxXY[2] - 1), 
                    (int)Math.ceil(minMaxXY[3] + 1), 
                    points, groupList, clusterFinder.getCriticalDensity(), 
                    "ran" + ii + "_" + i);
                   
                /*
                BufferedImage img = createImage(points, width, height);

                addAlternatingColorPointSets(img, groupList, 1);

                writeImage(img, "dt_ran_" + ii + "_" + i + ".png");
                */
                
                plotter.writeFile();
                
                float r0 = r0s[count];
                float r1 = r1s[count];
                
                float critDens = clusterFinder.getCriticalDensity();
            
                log.info("i=" + i + " ro-" + r0 + " r1=" + r1 + " critDens=" + critDens);
            
                /*
                if (i == 0) {
                    assertTrue(critDens >= r0 && (critDens <= 0.17));
                } else {
                    assertTrue(critDens >= r0 && (critDens <= (r1 + 0.1f*r1)));
                }
                */
                
                count++;
            }
        }
        
        plotter.writeFile();

        log.info("SEED=" + seed);
    }
    
    public void testFindClustersOtherData() throws Exception {
        
        String[] fileNames = {
            "Aggregation.txt", "Compound.txt", "Pathbased.txt" , "Spiral.txt",
            "D31.txt", "R15.txt" , "Jain.txt", "Flame.txt",
            //"a1.txt", "a2.txt", "a3.txt"
            /*,
            "s1.txt", "s2.txt", "s3.txt", "s4.txt",
            "birch1.txt", "birch2.txt", "birch3.txt" */
        };
        
        ClusterPlotter plotter = new ClusterPlotter();
        
        for (int i = 0; i < fileNames.length; i++) {

            String fileName = fileNames[i];
            
            //NOTE:  for i=8, distance transform needs alot of memory for array size, so have divided numbers there by 10
            AxisIndexer indexer = CreateClusterDataTest.getUEFClusteringDataset(
                fileName);
            
            Set<PairInt> points = new HashSet<PairInt>();
            for (int k = 0; k < indexer.getNXY(); ++k) {
                PairInt p = new PairInt(Math.round(indexer.getX()[k]),
                    Math.round(indexer.getY()[k]));
                points.add(p);
            }

            int[] minMaxXY = MiscMath.findMinMaxXY(points);
            int width = minMaxXY[1] + 1;
            int height = minMaxXY[3] + 1;

            DTClusterFinder<PairInt> clusterFinder = 
                new DTClusterFinder<PairInt>(points, width, height);

            clusterFinder.calculateCriticalDensity();
            clusterFinder.findClusters();
            //clusterFinder.setCriticalDensity(dens);

            int nGroups = clusterFinder.getNumberOfClusters();

            List<Set<PairInt>> groupList = new ArrayList<Set<PairInt>>();
            for (int k = 0; k < nGroups; ++k) {
                Set<PairInt> set = clusterFinder.getCluster(k);
                groupList.add(set);
            }
            
            plotter.addPlotWithoutHull(
                (int)Math.floor(minMaxXY[0] - 1), 
                (int)Math.ceil(minMaxXY[1] + 1), 
                (int)Math.floor(minMaxXY[2] - 1), 
                (int)Math.ceil(minMaxXY[3] + 1), 
                points, groupList, clusterFinder.getCriticalDensity(), 
                "other_" + i);

            /*
            BufferedImage img = createImage(points, width, height);

            addAlternatingColorPointSets(img, groupList, 0);

            writeImage(img, "dt_other_" + i + ".png");
            */
        }
        
        plotter.writeFile2();
    }
    
    public void testNoClusters() throws Exception {
        
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
               
            Can see that after N=1000, begin to see groups form from randomly close points.
                roughly nGroups is less than or equal to (0.2*N)
         */
        
        int[] numberOfBackgroundPoints = new int[]{100,450,900,4500,9000,14500,19000,24500,29000};
        
        int[] nGroupsFound = new int[numberOfBackgroundPoints.length];
        float[] expectedLinearDensities = new float[nGroupsFound.length];
        float[] calcLinearDensities = new float[nGroupsFound.length];
        
        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        long seed = System.currentTimeMillis();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");

        //seed = 1387775326745l;

        log.info("SEED=" + seed);
        
        sr.setSeed(seed);

        AxisIndexer indexer = null;

        int count = 0;
        
        ClusterPlotter plotter = new ClusterPlotter();
        
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
                        
            log.info(" " + ii + " (" + indexer.getNumberOfPoints() + " points) ... ");

            Set<PairInt> points = new HashSet<PairInt>();
            for (int k = 0; k < indexer.getNXY(); ++k) {
                PairInt p = new PairInt(Math.round(indexer.getX()[k]),
                    Math.round(indexer.getY()[k]));
                points.add(p);
            }

            int[] minMaxXY = MiscMath.findMinMaxXY(points);
            int width = minMaxXY[1] + 1;
            int height = minMaxXY[3] + 1;

            DTClusterFinder<PairInt> clusterFinder 
                = new DTClusterFinder<PairInt>(points, width, height);

            clusterFinder.calculateCriticalDensity();
            clusterFinder.findClusters();
            //clusterFinder.setCriticalDensity(dens);

            int nGroups = clusterFinder.getNumberOfClusters();
            
            float critDensity = clusterFinder.getCriticalDensity();

            List<Set<PairInt>> groupList = new ArrayList<Set<PairInt>>();
            for (int k = 0; k < nGroups; ++k) {
                Set<PairInt> set = clusterFinder.getCluster(k);
                groupList.add(set);
            }
            
            plotter.addPlotWithoutHull(
                (int)Math.floor(minMaxXY[0] - 1), 
                (int)Math.ceil(minMaxXY[1] + 1), 
                (int)Math.floor(minMaxXY[2] - 1), 
                (int)Math.ceil(minMaxXY[3] + 1), 
                points, groupList, critDensity, "no_clusters_" + ii);
                        
            float frac = ((float)nGroups/(float)points.size());
            log.info("nPoints=" + points.size() + " nGroups=" + nGroups
                + " frac=" + frac);
            
            assertTrue(frac < 0.2);
        }
        
        plotter.writeFile3();
        
        log.info("SEED=" + seed);
    }
    
    public static String writeDataset(float[] values, String fileName) throws 
        IOException {
        
        if (values == null) {
            throw new IllegalArgumentException("values cannot be null");
        }

        String outFilePath = ResourceFinder.findDirectory("bin") + "/" +
            fileName;
        
        FileOutputStream fs = null;
        ObjectOutputStream os = null;

        try {
            File file = new File(outFilePath);
            file.delete();
            file.createNewFile();

            fs = new FileOutputStream(file);
            os = new ObjectOutputStream(fs);
            
            os.writeInt(values.length);

            int count = 0;

            for (float v : values) {

                os.writeFloat(v);

                if ((count % 10) == 0) {
                    os.flush();
                }

                count++;
            }

            os.flush();

        } finally {

            if (os != null) {
                os.close();
            }
            if (fs != null) {
                fs.close();
            }
        }
        
        return outFilePath;        
    }
    
    private BufferedImage createImage(Set<PairInt> points, int width, int height) {
        
        BufferedImage outputImage = new BufferedImage(width, height, 
            BufferedImage.TYPE_INT_RGB);
        
        for (PairInt p : points) {
            int rgb = Color.WHITE.getRGB();
            outputImage.setRGB(p.getX(), p.getY(), rgb);
        }
        
        return outputImage;
    }

    private void addAlternatingColorPointSets(BufferedImage img, 
        List<Set<PairInt>> groups, int nExtraForDot) {
        
        int width = img.getWidth();
        int height = img.getHeight();
        
        int clrCount = -1;
        int clr = -1;
        
        for (int i = 0; i < groups.size(); ++i) {
                                    
            clr = getNextColorRGB(clrCount);
            
            Set<PairInt> group = groups.get(i);
            
            for (PairInt p : group) {
                
                int x = p.getX();
                int y = p.getY();

                for (int dx = (-1 * nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                    float xx = x + dx;
                    if ((xx > -1) && (xx < (width - 1))) {
                        for (int dy = (-1 * nExtraForDot); dy < (nExtraForDot + 1); dy++) {
                            float yy = y + dy;
                            if ((yy > -1) && (yy < (height - 1))) {
                                img.setRGB((int) xx, (int) yy, clr);
                            }
                        }
                    }
                }
                img.setRGB(p.getX(), p.getY(), clr);
                clrCount++;
            }
        }
    }
    
    public int getNextColorRGB(int clrCount) {
        
        if (clrCount == -1) {
            clrCount = 0;
        }
        
        clrCount = clrCount % 6;
        
        int c = Color.BLUE.getRGB();
        switch (clrCount) {
            case 0:
                c = Color.BLUE.getRGB();
                break;
            case 1:
                c = Color.PINK.getRGB();
                break;
            case 2:
                c = Color.GREEN.getRGB();
                break;
            case 3:
                c = Color.RED.getRGB();
                break;
            case 4:
                c = Color.CYAN.getRGB();
                break;
            case 5:
                c = Color.MAGENTA.getRGB();
                break;
            default:
                break;
        }
        
        return c;
    }

    private void writeImage(BufferedImage img, String fileName) throws IOException {
        
        String outFilePath = ResourceFinder.findDirectory("bin") + "/" +
            fileName;
        
        ImageIO.write(img, "PNG", new File(outFilePath));
    }
}
