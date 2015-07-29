package algorithms.compGeometry;

import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
            
            // fill up the first bin if it's smaller than subsetSize
            int lastIdx = n - 1;
            while ((output.get(0).getN() < subsetSize) && (lastIdx > 0)) {
                int nToAdd = subsetSize - output.get(0).getN();
                PairIntArray p = output.get(lastIdx);
                int canAdd = (p.getN() < nToAdd) ? p.getN() : nToAdd;
                for (int i = 0; i < canAdd; ++i) {
                    int idx = p.getN() - 1 - i;
                    output.get(0).add(p.getX(idx), p.getY(idx));
                }
                for (int i = 0; i < nToAdd; ++i) {
                    p.removeRange(p.getN() - 1, p.getN() - 1);
                }
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
    
    public PairIntArray[] reduceByBinSampling(PairIntArray points, 
        int nBinsPerDimension, int nPerBin) {
        
        PairIntArray[] binnedPoints = partition(points, nBinsPerDimension);
        
        PairIntArray reduced = new PairIntArray();
        PairIntArray remaining = new PairIntArray();
     
        try {
            SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
            long seed = System.currentTimeMillis();
            sr.setSeed(seed);
            
            for (PairIntArray binOfPoints : binnedPoints) {
                if (binOfPoints.getN() <= nPerBin) {
                    reduced.addAll(binOfPoints);
                    continue;  
                } 
                Set<Integer> chosen = new HashSet<Integer>();
                while (chosen.size() < nPerBin) {
                    int idx = sr.nextInt(binOfPoints.getN());
                    while (chosen.contains(Integer.valueOf(idx))) {
                        idx = sr.nextInt(binOfPoints.getN());
                    }
                    reduced.add(binOfPoints.getX(idx), binOfPoints.getY(idx));
                    chosen.add(Integer.valueOf(idx));
                }
                for (int idx = 0; idx < binOfPoints.getN(); ++idx) {
                    if (chosen.contains(Integer.valueOf(idx))) {
                        continue;
                    }
                    remaining.add(binOfPoints.getX(idx), binOfPoints.getY(idx));
                }
            }
                   
        } catch (NoSuchAlgorithmException ex) {
            
            throw new RuntimeException(ex);
        }
        
        return new PairIntArray[]{reduced, remaining};
    }
}
