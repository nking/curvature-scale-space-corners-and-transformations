package algorithms.compGeometry;

import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;

/**
 *
 * @author nichole
 */
public class PointPartitioner {
 
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
