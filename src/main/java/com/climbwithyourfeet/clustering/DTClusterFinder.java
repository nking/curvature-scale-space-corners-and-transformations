package com.climbwithyourfeet.clustering;

import com.climbwithyourfeet.clustering.util.MiscMath;
import com.climbwithyourfeet.clustering.util.PairInt;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.Arrays;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;

/**
 * main class to cluster finder whose logic is based upon distance transform,
 * density threshold, and a signal to noise argument of ~ 3 as a factor.
 * 
 * runtime complexity ~ O(N_pixels) + ~O(N_points * lg2(N_points))
 * 
 * @author nichole
 */
public class DTClusterFinder {
    
    private final Set<PairInt> points;
    private final int width;
    private final int height;
    
    private float critDens = Float.POSITIVE_INFINITY;
    
    private DTGroupFinder groupFinder = null;

    private enum STATE {
        INIT, HAVE_CLUSTER_DENSITY, HAVE_GROUPS
    }
    
    private STATE state = null;
    
    private float threshholdFactor = 2.5f;
    
    private boolean debug = false;

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public DTClusterFinder(Set<PairInt> thePoints, int width, int height) {
        
        this.points = thePoints;
        this.width = width;
        this.height = height;
        
        state = STATE.INIT;
    }
    
    public void setToDebug() {
        debug = true;
    }
    
    public void setThreshholdFactor(float factor) {
        this.threshholdFactor = factor;
    }
    
    public void calculateCriticalDensity() {
        
        if (state.compareTo(STATE.HAVE_CLUSTER_DENSITY) > -1) {
            return;
        }
        
        DistanceTransform dtr = new DistanceTransform();
        int[][] dt = dtr.applyMeijsterEtAl(points, width, height);
        
        CriticalDensitySolver densSolver = new CriticalDensitySolver();
        
        if (debug) {
            
            densSolver.setToDebug();

            log.info("print dist trans for " + points.size() + " points " +
                "within width=" + width + " height=" + height);
            
            int[] minMax = MiscMath.findMinMaxValues(dt);
            
            log.info("min and max =" + Arrays.toString(minMax));
            
            try {
                writeDebugImage(dt, Long.toString(System.currentTimeMillis()));
            } catch (IOException ex) {
                Logger.getLogger(DTClusterFinder.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        densSolver.setThreshholdFactor(threshholdFactor);
        
        this.critDens = densSolver.findCriticalDensity(dt, points.size(), 
            width, height);  
        
        this.state = STATE.HAVE_CLUSTER_DENSITY;
    }
    
    public void setCriticalDensity(float dens) {
        
        if (state.compareTo(STATE.HAVE_CLUSTER_DENSITY) > -1) {
            throw new IllegalStateException("cluster density is already set");
        }
        
        this.critDens = dens;
        
        this.state = STATE.HAVE_CLUSTER_DENSITY;
    }
    
    public void findClusters() {
        
        if (state.compareTo(STATE.HAVE_CLUSTER_DENSITY) < 0) {
            calculateCriticalDensity();
        } else if (state.compareTo(STATE.HAVE_GROUPS) >= 0) {
            return;
        }
        
        groupFinder = new DTGroupFinder();
        
        groupFinder.setThreshholdFactor(threshholdFactor);
        
        groupFinder.calculateGroups(critDens, points);
        
    }
    
    public int getNumberOfClusters() {
        
        if (groupFinder == null) {
            return 0;
        }
        
        return groupFinder.getNumberOfGroups();
    }
    
    public Set<PairInt> getCluster(int idx) {
        
        if (groupFinder == null) {
            throw new IllegalArgumentException(
                "findClusters was not successfully invoked");
        }
        
        if ((idx < 0) || (idx > (groupFinder.getNumberOfGroups() - 1))) {
            throw new IllegalArgumentException("idx is out of bounds");
        }
        
        return groupFinder.getGroup(idx);
    }
    
    public float getCriticalDensity() {
        return critDens;
    }
    
    private void writeDebugImage(int[][] dt, String fileSuffix) throws IOException {
        
        BufferedImage outputImage = new BufferedImage(width, height, 
            BufferedImage.TYPE_BYTE_GRAY);

        WritableRaster raster = outputImage.getRaster();
        
        for (int i = 0; i < dt.length; ++i) {
            for (int j = 0; j < dt[0].length; ++j) {
                int v = dt[i][j];
                raster.setSample(i, j, 0, v);
            }
        }
        
        // write to an output directory.  we have user.dir from system properties
        // but no other knowledge of users's directory structure
        URL baseDirURL = this.getClass().getClassLoader().getResource(".");
        String baseDir = null;
        if (baseDirURL != null) {
            baseDir = baseDirURL.getPath();
        } else {
            baseDir = System.getProperty("user.dir");
        }
        if (baseDir == null) {
            return;
        }
        File t = new File(baseDir + "/bin");
        if (t.exists()) {
            baseDir = t.getPath();
        } else if ((new File(baseDir + "/target")).exists()) {
            baseDir = baseDir + "/target";
        }
        
        // no longer need to use file.separator
        String outFilePath = baseDir + "/distance_transform_" + fileSuffix + ".png";
        
        ImageIO.write(outputImage, "PNG", new File(outFilePath));
        
        Logger.getLogger(this.getClass().getName()).info("wrote " + outFilePath);
    }
    
}
