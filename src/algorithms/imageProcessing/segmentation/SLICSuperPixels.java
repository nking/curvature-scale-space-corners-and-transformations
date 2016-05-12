package algorithms.imageProcessing.segmentation;

import algorithms.QuickSort;
import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.Gaussian1DFirstDeriv;
import algorithms.imageProcessing.ImageExt;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

/**
 * a variant of kmeans whose goal is to make k super-pixels in the image
 * based upon CIE Lab color similarity and x,y proximity.
 * 
 * The code follows the algorithm of:
 * "SLIC Superpixels Compared to State-of-the-Art Superpixel Methods"
   by Achanta, Appu Shaji,Smith,  Lucchi, Fua, and Su ̈sstrunk,
 *  
 * @author nichole
 */
public class SLICSuperPixels {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    protected int[] labels = null;
    
    protected double[] distances = null;
    
    protected final int k;
    
    protected final int s;
    
    protected final float[][] seedDescriptors;
    
    protected final ImageExt img;
    
    protected final double threshold;
    
    // range of CIE LAB values is (0,0,0) to (28.51 3.28 2.15)
    protected static final double maxClr = 28.78;
    
    public SLICSuperPixels(ImageExt img, int nClusters) {
        
        //TOOD: after have an implementation as authors suggest,
        //  change to use deltaE instead of sqrt sum diffs of CIE lab
        //  and compare differences in results and runtime (many more flops...)
        
        this.k = nClusters;
        
        double sampling = Math.sqrt(( (float)img.getNPixels()/(float)k));
        
        this.s = (sampling < 1) ? 1 : (int)Math.round(sampling);
        
        this.img = img;
        
        // l, a, b, x, y
        seedDescriptors = new float[k][];
        for (int i = 0; i < k; ++i) {
            seedDescriptors[i] = new float[5];
        }
        
        // max error would be ( maxClr * maxClr * k) + 2*( s/2 * s/2 * k)
        double maxError = (maxClr * maxClr * k) + (s*s/2) * k;
        maxError = Math.sqrt(maxError);
        threshold = 0.1 * maxError;
    }
    
    protected void init() {
        
        if (labels != null) {
            throw new IllegalStateException("variables have been initialized");
        }
        
        int nPix = img.getNPixels();
                
        labels = new int[nPix];
        Arrays.fill(labels, -1);
        
        distances = new double[nPix];
        Arrays.fill(distances, Double.MAX_VALUE);
        
        // init cluster centers, seedDescriptors for grid with cell size s
        populateSeedDescriptors();
    }
    
    private void populateSeedDescriptors() {
        
        double[][] gradient = calcGradient();
        
        /*
        sampled on a regular grid spaced S pixels apart. 
        To produce roughly equally sized superpixels, 
        the grid interval is S. 
        
        ** The centers are moved to seed locations corresponding to the lowest 
        gradient position in a 3 X 3 neighborhood. This is done to avoid 
        centering a superpixel on an edge and to reduce the chance of 
        seeding a superpixel with a noisy pixel.
        */
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        // determine the centers of each s x s cell within search range of 3x3 in gradient
        int nX = Math.round((float)w/(float)s);
        int nY = Math.round((float)h/(float)s);
        
        int sHalf = s/2;
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        // kCurrent = (iNy * nX) + iNx;
        for (int iNx = 0; iNx < nX; ++iNx) {
            int x1 = sHalf * iNx*s;
            for (int iNy = 0; iNy < nY; ++iNy) {
                
                int kCurrent = (iNy * nX) + iNx;
                
                int y1 = sHalf * iNy*s;
                
                // find smallest gradient within range (x1, y1) +-1 to set center for seed
                double minG = Double.MAX_VALUE;
                int minX2 = -1;
                int minY2 = -1;
                for (int m = 0; m < dxs.length; ++m) {
                    int x2 = x1 + dxs[m];
                    int y2 = y1 + dys[m];
                    if (x2 < 0 || y2 < 0 || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }
                    double g = gradient[x2][y2];
                    if (g < minG) {
                        minG = g;
                        minX2 = x2;
                        minY2 = y2;
                    }
                }
                // [l, a, b, x, y]
                seedDescriptors[kCurrent][3] = minX2;
                seedDescriptors[kCurrent][4] = minY2;
                
                float[] lab = img.getCIELAB(minX2, minY2);
                System.arraycopy(lab, 0, seedDescriptors[kCurrent], 0, 3);
            }
        }        
    }
    
   /**
    * O(N)
    * 
    * @return 
    */
   private double[][] calcGradient() {
       
       int nPix = img.getNPixels();
       int width = img.getWidth();
       int height = img.getHeight();
              
       double[][] gradient = new double[width][];              
       for (int i = 0; i < width; ++i) {
           gradient[i] = new double[height];
       }
       
       CIEChromaticity cieC = new CIEChromaticity();
       
       float[] kernel = Gaussian1DFirstDeriv.getBinomialKernelSigmaZeroPointFive();
              
       int h = (kernel.length - 1) >> 1;
       
       for (int i = 0; i < nPix; ++i) {
           final int x1 = img.getCol(i);
           final int y1 = img.getRow(i);
           
           float rXSum = 0;
           float gXSum = 0;
           float bXSum = 0;
           float rYSum = 0;
           float gYSum = 0;
           float bYSum = 0;
           
           for (int g = 0; g < kernel.length; ++g) {
               float gg = kernel[g];
               if (gg == 0) {
                   continue;
               }
                
               int x2, y2;
               // calc for X gradient first
               int delta = g - h;
                x2 = x1 + delta;
                y2 = y1;
                // edge corrections.  use replication
                if (x2 < 0) {
                    x2 = -1 * x2 - 1;
                } else if (x2 >= width) {
                    int diff = x2 - width;
                    x2 = width - diff - 1;
                }
                rXSum += gg * img.getR(x2, y2);
                gXSum += gg * img.getG(x2, y2);
                bXSum += gg * img.getB(x2, y2);
               
                // calc for y
                y2 = y1 + delta;
                x2 = x1;
                // edge corrections.  use replication
                if (y2 < 0) {
                    y2 = -1 * y2 - 1;
                } else if (y2 >= height) {
                    int diff = y2 - height;
                    y2 = height - diff - 1;
                }
                rYSum += gg * img.getR(x2, y2);
                gYSum += gg * img.getG(x2, y2);
                bYSum += gg * img.getB(x2, y2);
           }
           
           double rC = Math.sqrt(rXSum * rXSum + rYSum * rYSum);
           double gC = Math.sqrt(gXSum * gXSum + gYSum * gYSum);
           double bC = Math.sqrt(bXSum * bXSum + bYSum * bYSum);
           
           // presumably, greyscale gradient is fine
           double gr = (rC + gC + bC)/3;
           
           gradient[x1][y1] = gr;
       }
       
       return gradient;
   }
    
    public void calculate() {
        
        init();
        
        int nIterMax = 20;
        
        int nIter = 0;
        
        while (nIter < nIterMax) {
         
            assignPixelsNearSeeds();
            
            double l2Norm = adjustClusterCenters();
            
            System.out.println("l2Norm=" + l2Norm + " nIter=" + nIter);
            
            if (l2Norm < threshold) {
                break;
            }
            
            nIter++;
        }
        
        assignTheUnassigned();
    }

    private void assignPixelsNearSeeds() {

        for (int kCurrent = 0; kCurrent < k; ++kCurrent) {

            int x1 = (int) seedDescriptors[kCurrent][3];
            int y1 = (int) seedDescriptors[kCurrent][4];

            for (int x2 = (x1 - s); x2 <= (x1 + s); ++x2) {
                for (int y2 = (y1 - s); y2 <= (y1 + s); ++y2) {
                    if (x1 == x2 && y1 == y2) {
                        continue;
                    }
                    int pixIdx2 = img.getInternalIndex(x2, y2);

                    double dist = calcDist(img, seedDescriptors[kCurrent],
                        x2, y2);

                    if (dist < distances[pixIdx2]) {
                        distances[pixIdx2] = dist;
                        labels[pixIdx2] = kCurrent;
                    }
                }
            }
        }

    }

    private double adjustClusterCenters() {
        
        // L2 norm is the residuals added in quadratur
        
        double[][] meanDescriptors = new double[k][];        
        // l, a, b, x, y
        for (int i = 0; i < k; ++i) {
            meanDescriptors[i] = new double[5];
        } 
        
        // calculate new colors
        int[] count = new int[k];
        
        for (int i = 0; i < img.getNPixels(); ++i) {
            int label = labels[i];         
            int x = img.getCol(i);
            int y = img.getRow(i);
            float[] lab = img.getCIELAB(i);

            for (int m = 0; m < 3; ++m) {
                meanDescriptors[label][m] += lab[m];
            }
            meanDescriptors[label][3] += x;
            meanDescriptors[label][4] += y;
            count[label]++;
        }
        
        // l, a, b, x, y
        double[] sqDiffSum = new double[5];
        
        double diff;
        for (int kCurrent = 0; kCurrent < k; ++kCurrent) {            
            for (int m = 0; m < 5; ++m) {
                meanDescriptors[kCurrent][m] /= (float)count[kCurrent];
                diff = meanDescriptors[kCurrent][m] - seedDescriptors[kCurrent][m];
                sqDiffSum[m] += (diff * diff);
            }
        }
        
        double l2Norm = 0;
        for (double sd : sqDiffSum) {
            l2Norm += sd;
        }
        l2Norm = Math.sqrt(l2Norm);
        
        // calc sqrt of sum of sq diffs with old centers and reset old centers
        for (int kCurrent = 0; kCurrent < k; ++kCurrent) {            
            for (int m = 0; m < 5; ++m) {
                seedDescriptors[kCurrent][m] = (float)meanDescriptors[kCurrent][m];
            }
        }
        
        return l2Norm;
    }

    private double calcDist(ImageExt img, float[] desc1, int x2, int y2) {
        
        float[] lab2 = img.getCIELAB(x2, y2);
        
        double dClr = 0;
        for (int i = 0; i < 3; ++i) {
            float diff = desc1[i] - lab2[i];
            dClr += (diff * diff);
        }
        dClr = Math.sqrt(dClr)/maxClr;
        
        float diffX = desc1[3] - x2;
        float diffY = desc1[4] - y2;
        double dXY = Math.sqrt(diffX * diffX + diffY * diffY)/(double)s;
        
        double dComb = Math.sqrt(dClr*dClr + dXY * dXY);
        
        return dComb;
    }

    private void assignTheUnassigned() {
        
        Map<PairInt, Set<Integer>> unassignedMap = new HashMap<PairInt, Set<Integer>>();
        
        int w = img.getWidth();
        int h = img.getHeight();
        int n = img.getNPixels();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int pixIdx = img.getInternalIndex(i, j);
                if (labels[pixIdx] == -1) {
                    PairInt p = new PairInt(i, j);
                    
                    Set<Integer> adjLabels = unassignedMap.get(p);
                    if (adjLabels == null) {
                        adjLabels = new HashSet<Integer>();
                        unassignedMap.put(p, adjLabels);
                    }
                    
                    for (int m = 0; m < dxs.length; ++m) {
                        int x2 = i + dxs[m];
                        int y2 = j + dys[m];
                        if (x2 < 0 || y2 < 0 || (x2 > (w - 1)) || (y2 > (h - 1))) {
                            continue;
                        }
                        int pixIdx2 = img.getInternalIndex(x2, y2);
                        if (labels[pixIdx2] > -1) {
                            adjLabels.add(Integer.valueOf(labels[pixIdx2]));
                        }
                    }
                }
            }
        }
        
        ArrayDeque<PairInt> queue0 = populateByNumberOfNeighbors(unassignedMap);
        
        ArrayDeque<PairInt> queue1 = new ArrayDeque<PairInt>();
        
        while (true) {
        
            while (!queue0.isEmpty()) {

                PairInt p = queue0.poll();

                Set<Integer> adjLabels = unassignedMap.get(p);
                assert(adjLabels != null);

                double minD = Double.MAX_VALUE;
                Integer minLabel2 = null;

                for (Integer label2 : adjLabels) {

                    double dist = calcDist(img, seedDescriptors[label2.intValue()], 
                        p.getX(), p.getY());

                    if (dist < minD) {
                        minD = dist;
                        minLabel2 = label2;
                    }
                }

                int pixIdx1 = img.getInternalIndex(p.getX(), p.getY());
                labels[pixIdx1] = minLabel2.intValue();
                distances[pixIdx1] = minD;

                unassignedMap.remove(p);

                for (int m = 0; m < dxs.length; ++m) {
                    int x2 = p.getX() + dxs[m];
                    int y2 = p.getY() + dys[m];
                    if (x2 < 0 || y2 < 0 || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }
                    int pixIdx2 = img.getInternalIndex(x2, y2);
                    if (labels[pixIdx2] == -1) {
                        queue1.add(new PairInt(x2, y2));
                    }
                }
            }
            if (queue1.isEmpty()) {
                break;
            }
            queue0.addAll(queue1);
            queue1.clear();
        }
       
        assert(unassignedMap.isEmpty());       
    }

    private ArrayDeque<PairInt> populateByNumberOfNeighbors(
        Map<PairInt, Set<Integer>> unassignedMap) {
        
        int n = unassignedMap.size();
        
        PairInt[] points = new PairInt[n];
        int[] nN = new int[n];
        
        int count = 0;
        for (Entry<PairInt, Set<Integer>> entry : unassignedMap.entrySet()) {
            points[count] = entry.getKey();
            nN[count] = entry.getValue().size();
            count++;
        }
        
        QuickSort.sortBy1stArg(nN, points);
        
        ArrayDeque<PairInt> queue = new ArrayDeque<PairInt>();
        
        for (int i = (n - 1); i > -1; --i) {
            
            int nP = nN[i];
            if (nP == 0) {
                break;
            }
            
            queue.add(points[i]);
        }
        
        return queue;
    }
}
