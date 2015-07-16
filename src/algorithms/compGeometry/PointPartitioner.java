package algorithms.compGeometry;

import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author nichole
 */
public class PointPartitioner {
 
    /**
     * make subsets of size subsetSize out of points by randomly distributing
     * the points into the subsets.  Note that the last subset will be smaller
     * than subsetSize if (points.getN() % subsetSize) != 0.
     * 
     * @param points
     * @param subsetSize
     * @return 
     */
    public List<PairIntArray> randomSubsets(PairIntArray points, int subsetSize) {
        
        if (points == null) {
            throw new IllegalArgumentException("points cannot be null");
        }
        
        if (points.getN() == 0) {
            return new ArrayList<PairIntArray>();
        }
        
        try {
            SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
            long seed = System.currentTimeMillis();
            sr.setSeed(seed);
            
            int n = (int)Math.ceil((float)points.getN()/(float)subsetSize);
            
            List<PairIntArray> output = new ArrayList<PairIntArray>();
            for (int i = 0; i < n; ++i) {
                output.add(new PairIntArray());
            }
            
            for (int i = 0; i < points.getN(); ++i) {
                
                int index = sr.nextInt(n);
                
                while (output.get(index).getN() == subsetSize) {
                    index++;
                    if (index == n) {
                        index = 0;
                    }
                }
                
                output.get(index).add(points.getX(i), points.getY(i));
            }
            
            return output;
        
        } catch (NoSuchAlgorithmException ex) {
            
            throw new RuntimeException(ex);
        }
    }
    
    public PairIntArray[] partition(PairIntArray points, int nDimensions) {
        
        PairIntArray[] partitions = new PairIntArray[nDimensions*nDimensions];
        for (int i = 0; i < partitions.length; i++) {
            partitions[i] = new PairIntArray();
        }
        
        int minX = MiscMath.findMin(points.getX());
        int maxX = MiscMath.findMax(points.getX());
        int minY = MiscMath.findMin(points.getY());
        int maxY = MiscMath.findMax(points.getY());
        
        float binX = (float)(maxX - minX)/(float)nDimensions;
        float binY = (float)(maxY - minY)/(float)nDimensions;
    
        for (int i = 0; i < points.getN(); i++) {
            
            int x = points.getX(i);
            int y = points.getY(i);
            
            int col = (int)((x - minX)/binX);
            
            int row = (int)((y - minY)/binY);
            
            if (col > (nDimensions - 1)) {
                col = nDimensions - 1;
            }
            
            if (row > (nDimensions - 1)) {
                row = nDimensions - 1;
            }
            
            int idx = (row * nDimensions) + col;
            
            partitions[idx].add(x, y);
            
        }
        
        return partitions;
    }
    
    public PairIntArray[] partitionVerticalOnly(PairIntArray points, int n) {
        
        PairIntArray[] partitions = new PairIntArray[n];
        for (int i = 0; i < partitions.length; i++) {
            partitions[i] = new PairIntArray();
        }
        
        int minX = MiscMath.findMin(points.getX());
        int maxX = MiscMath.findMax(points.getX());
        
        float binX = (float)(maxX - minX)/(float)n;
    
        for (int i = 0; i < points.getN(); i++) {
            
            int x = points.getX(i);
            int y = points.getY(i);
            
            int col = (int)((x - minX)/binX);
                        
            if (col > (n - 1)) {
                col = n - 1;
            }
            
            int idx = col;
            
            partitions[idx].add(x, y);
            
        }
        
        return partitions;
    }
    
    public PairIntArray[] partitionHorizontalOnly(PairIntArray points, int n) {
        
        PairIntArray[] partitions = new PairIntArray[n];
        for (int i = 0; i < partitions.length; i++) {
            partitions[i] = new PairIntArray();
        }
        
        int minY = MiscMath.findMin(points.getY());
        int maxY = MiscMath.findMax(points.getY());
        
        float bin = (float)(maxY - minY)/(float)n;
    
        for (int i = 0; i < points.getN(); i++) {
            
            int x = points.getX(i);
            int y = points.getY(i);
            
            int row = (int)((x - minY)/bin);
                        
            if (row > (n - 1)) {
                row = n - 1;
            }
            
            int idx = row;
            
            partitions[idx].add(x, y);
            
        }
        
        return partitions;
    }
}
