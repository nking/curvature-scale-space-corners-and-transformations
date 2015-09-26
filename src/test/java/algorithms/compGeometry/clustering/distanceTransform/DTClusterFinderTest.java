package algorithms.compGeometry.clustering.distanceTransform;

import algorithms.compGeometry.clustering.twopointcorrelation.AxisIndexer;
import algorithms.compGeometry.clustering.twopointcorrelation.BaseTwoPointTest;
import algorithms.compGeometry.clustering.twopointcorrelation.CreateClusterDataTest;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.security.NoSuchAlgorithmException;
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

    public void estFindRanGenClusters() throws Exception {
        
        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        long seed = System.currentTimeMillis();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");

        seed = 1387775326745l;

        log.info("SEED=" + seed);
        
        //sr.setSeed(seed);

        int nSwitches = 3;

        int nIterPerBackground = 3;

        AxisIndexer indexer = null;

        int count = 0;

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
                        //~1000
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 33, 33, 10f);
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
                    
                DTClusterFinder clusterFinder = new DTClusterFinder(points,
                    width, height);

                clusterFinder.calculateCriticalDensity();
                clusterFinder.findClusters();
                //clusterFinder.setCriticalDensity(dens);

                int nGroups = clusterFinder.getNumberOfClusters();

                List<Set<PairInt>> groupList = new ArrayList<Set<PairInt>>();
                for (int k = 0; k < nGroups; ++k) {
                    Set<PairInt> set = clusterFinder.getCluster(k);
                    groupList.add(set);
                }
                   
                BufferedImage img = createImage(points, width, height);

                addAlternatingColorPointSets(img, groupList, 1);

                writeImage(img, "dt_ran_" + ii + "_" + i + ".png");
             
                count++;
            }
        }

        count = 0;

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

            DTClusterFinder clusterFinder = new DTClusterFinder(points,
                width, height);

            clusterFinder.calculateCriticalDensity();
            clusterFinder.findClusters();
            //clusterFinder.setCriticalDensity(dens);

            int nGroups = clusterFinder.getNumberOfClusters();

            List<Set<PairInt>> groupList = new ArrayList<Set<PairInt>>();
            for (int k = 0; k < nGroups; ++k) {
                Set<PairInt> set = clusterFinder.getCluster(k);
                groupList.add(set);
            }

            BufferedImage img = createImage(points, width, height);

            addAlternatingColorPointSets(img, groupList, 0);

            writeImage(img, "dt_other_" + i + ".png");
        }
    }
    
    public void testNoClusters() {
        
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
