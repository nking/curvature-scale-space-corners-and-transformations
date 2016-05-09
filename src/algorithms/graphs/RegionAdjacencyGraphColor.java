package algorithms.graphs;

import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.segmentation.ColorSpace;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.awt.Color;
import java.util.Set;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * class to hold a list of region nodes and calculate an adjacency map
 * for them and pairwise color differences or similarity.
 * 
 * @author nichole
 */
public class RegionAdjacencyGraphColor extends RegionAdjacencyGraph {
    
    /*
    the edge weights between 2 points need 1/2 of npixel X npixel storage,
    so need to use compressed data structures.
    
    The ia4j api has CRS matrix and array implementations.
        Matrix a = Matrices.CRS_FACTORY.createMatrix(100000, 100000);
    
    There is also LAML, but I could not find a benchmark results including it.
    http://web.engr.illinois.edu/~mqian2/upload/projects/java/LAML/doc/index.html
    
    The convention adopted here is to always access the values for a pair of
    points by the smallest index first.
    
       4
       3
       2
       1  
       0  
         0  1  2  3
    
    For each pixel in img, the color difference or similarity to it's adjacent 
    pixels are calculated.  
    For each pixel, using a pattern of only calculating the neighbors for
    offsets (+1,0), (+1,+1), (0,+1), (-1, +1) calculates the edge without repeating
    a pair.
    
    The conversion from col, row to a single index will use the convention in 
    Image.java which is index = (row * width) + col.
    
    Also note that the convention used with sparse matrix diffOrSim will be to
    to place the smallest index first: [smallestIndex][largerIndex] for
    any 2 pixel indexes.
    
    */
    protected SimpleMatrix diffOrSim = null;
    
    protected static final int[] dxNbrs = new int[]{1, 1, 0, -1};
    protected static final int[] dyNbrs = new int[]{0, 1, 1,  1};
        
    protected ColorSpace colorSpace = null;
    
    protected final ImageExt img;
    
     /**
     * constructor.  Note that the similarity matrix created is nPix X nPix
       in size, and currently the data structure is limited in size to max value of
       an integer, so that means images with fewer than 46340 pixels can be processed
       with this.
       A 250 X 180 image will need at least 2.75 GB of memory.
       In the future, this method might create an overlapping tiling handling to handle larger
       images and use a compact data structure (espec since normalized cuts does not need 
       a high significance in digits for edge weights).  
       
       Currently, best use would be pre-shrinking the image using pyramidal decimation or other binning
       to reduce the image to size closer to 100 x 100 pixels or smaller.
       
       (Note that if the statement in the normalized cuts paper regarding bit significance
       is followed, I should be able to factor the doubles to integers.
        in that case, can use compression as I do in Image.java and
        find a linear algebra eigen solver which will use the get and set methods
        of my image class to keep the data small.
        Also note, that instead of using an eigen solver, could alternatively,
        use my two-point clustering code for the remaining logic.  the distance transform
        and histograms make the remaining logic very fast, but it is dependent upon
        the bounds of the data (that is maximum values of data points which would be
        differences here converted to integers).
       
     * @param img 
     * @param labels double array of labels for each pixel using the convention
     * labels[col][row].  Note that the largest label must be less than 
     * the number of pixels in the image.
     */
    public RegionAdjacencyGraphColor(ImageExt img, int[][] labels) {

        super(img, labels);  
        
        this.img = img;

        if (img.getNPixels() > 46340) {
            throw new IllegalStateException("an internal data structure is limited by " 
            + "max integer size currently.  Please reduce image size to less than " 
            + "46340 pixels and try again.");
        }
        
    }
   
    public void populateEdgesWithColorSimilarity(ColorSpace clrSpace) {
        
        populatePairDifferences(clrSpace);

        double sigma;
        if (clrSpace.equals(clrSpace.HSV)) {
            // max possible = sqrt(1*1 + 1*1 + 1*1) = sqrt3
            sigma = 1.7320508;

        } else if (clrSpace.equals(clrSpace.CIELAB)) {
            // max possible deltaE2000 = 19.22
            sigma = 19.22;

        } else {
             // rgb
            // max possible = math.sqrt(3*255*255) = sqrt(3)*255
            sigma = 441.6729559;
        }
        
        // offsets for neighbors of each pixel to calculate: (+1,0), (+1,+1), (0,+1)
        
        for (int row = 0; row < imageHeight; ++row) {
            for (int col = 0; col < imageWidth; ++col) {
                int idx1 = calculatePixelIndex(row, col);                
                for (int k = 0; k < dxNbrs.length; ++k) {
                    int col2 = col + dxNbrs[k];
                    int row2 = row + dyNbrs[k];
                    if ((col2 > (imageWidth - 1)) || (row2 > (imageHeight - 1))) {
                        continue;
                    }
                    int idx2 = calculatePixelIndex(row2, col2); 
                    
                    if (idx1 < idx2) {
                        double d = diffOrSim.get(idx1, idx2);
                        double similarity = Math.exp(-1*d*d/sigma);
                        diffOrSim.set(idx1, idx2, similarity);
                    } else {
                        double d = diffOrSim.get(idx2, idx1);
                        double similarity = Math.exp(-1*d*d/sigma);
                        diffOrSim.set(idx2, idx1, similarity);
                    }
                }
            }
        }
        
    }
    
    public boolean isWithinImageBounds(int col, int row) {
        if ((col > (imageWidth - 1)) || (row > (imageHeight - 1)) ||
            (col < 0) || (row < 0)) {
            return false;
        }
        return true;
    }
    
    public int calculatePixelIndex(int row, int col) {
        return (row * imageWidth) + col;
    }
    public int getRowFromPixelIndex(int pixIdx) {
        return pixIdx/imageWidth;
    }
    public int getColFromPixelIndex(int pixIdx) {
        int row = pixIdx/imageWidth;
        return pixIdx - (row * imageWidth);
    }
        
    private void populatePairDifferences(ColorSpace clrSpace) {
        
        if (diffOrSim != null) {
            throw new IllegalStateException("this method is expected to be inoked only once");
        }
        
        this.colorSpace = clrSpace;
        
        int nPix = img.getNPixels();
        int nCols = img.getWidth();
        int nRows = img.getHeight();

        // this project is built for single threaded use, so using UJMP
        //    may be slower due to use of mutexes and synchronization needed
        //    for multi-threaded use,
        //    but need an API which handles memory well,
        //    so trying it 
        
        diffOrSim = new SimpleMatrix(nPix, nPix);
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        // offsets for neighbors of each pixel to calculate: (+1,0), (+1,+1), (0,+1)

        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {
                int r1 = img.getR(col, row);
                int g1 = img.getG(col, row);
                int b1 = img.getB(col, row);
                int idx1 = calculatePixelIndex(row, col);                

                for (int k = 0; k < dxNbrs.length; ++k) {
                    final int col2 = col + dxNbrs[k];
                    final int row2 = row + dyNbrs[k];
                    if ((col2 > (nCols - 1)) || (row2 > (nRows - 1))) {
                        continue;
                    }
                    int r2 = img.getR(col2, row2);
                    int g2 = img.getG(col2, row2);
                    int b2 = img.getB(col2, row2);
                    
                    double dist;
                    
                    if (clrSpace.equals(clrSpace.HSV)) {

                        float[] hsb1 = new float[3];
                        Color.RGBtoHSB(r1, g1, b1, hsb1);

                        float[] hsb2 = new float[3];
                        Color.RGBtoHSB(r2, g2, b2, hsb2);

                        float sumDiff = 0;
                        for (int m = 0; m < hsb1.length; ++m) {
                            float diff = hsb1[m] - hsb2[m];
                            sumDiff += (diff * diff);
                        }
                        dist = Math.sqrt(sumDiff);

                    } else if (clrSpace.equals(clrSpace.CIELAB)) {

                        float[] cieLAB1 = cieC.rgbToCIELAB(r1, g1, b1);

                        float[] cieLAB2 = cieC.rgbToCIELAB(r2, g2, b2);

                        dist = Math.abs(cieC.calcDeltaECIE2000(cieLAB1, cieLAB2));

                    } else {

                        //RGB
                        int diffR = r1 - r2;
                        int diffG = g1 - g2;
                        int diffB = b1 - b2;
                        dist = Math.sqrt(diffR * diffR + diffG * diffG + diffB * diffB);
                    }
                    
                    int idx2 = calculatePixelIndex(row2, col2); 

                    if (idx1 < idx2) {
                        diffOrSim.set(idx1, idx2, dist);
                    } else {
                        diffOrSim.set(idx2, idx1, dist);
                    }
                }
            }
        }
    }
    
    /**
     * calculate the normalized cute between regions regionIndex1 and
     * regionIndex2 using
     * normalized_cut(A, B) = 2 - normalized_assoc(A, B);
     * 
     * @param regionIndex1
     * @param regionIndex2
     * @return 
     */
    public double calculateNormalizedCut(int regionIndex1, int regionIndex2) {
        
        double sum = calculateNormalizedAssoc(regionIndex1) + 
            calculateNormalizedAssoc(regionIndex2);
        
        return 2 - sum;
    }
    
    /**
     <pre>
                               assoc(A, A) 
         normalized_assoc(A) = ----------- 
                               assoc(A, V)
     </pre>
     * @param regionIndex1
     * @param regionIndex2
     * @return 
     */
    double calculateNormalizedAssoc(int regionIndex) {
        double sum0 = calculateSelfAssoc(regionIndex);
        double sum1 = calculateAssocOfRegion(regionIndex);
        if (sum1 == 0) {
            return 0;
        }
        return sum0/sum1;
    }
    
    /**
     * calculate the sum of all of regionIndex's edge's for nodes within 
     * regionIndex.
     * @param regionIndex
     * @return 
     */
    double calculateSelfAssoc(int regionIndex) {
        
        Region region = regions.get(regionIndex);
        
        Set<PairInt> points = region.getPoints();
        
        // uniquely adding relationships by only
        // adding the neighbors (+1,0), (+1,+1), (0,+1)
        // after testing that it exists within points
        double sum = 0;
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            int idx1 = calculatePixelIndex(y, x);
            for (int k = 0; k < dxNbrs.length; ++k) {
                int x2 = x + dxNbrs[k];
                int y2 = y + dyNbrs[k];
                PairInt p2 = new PairInt(x2, y2);
                if (points.contains(p2)) {
                    int idx2 = calculatePixelIndex(y2, x2);
                    if (idx1 < idx2) {
                        sum += diffOrSim.get(idx1, idx2);
                    } else {
                        sum += diffOrSim.get(idx2, idx1);
                    }
                }
            }
        }
        return sum;
    }
    
    /**
     * calculate the sum of all of regionIndex's edge's for all nodes it is
     * adjacent to.
     * @param regionIndex
     * @return 
     */
    double calculateAssocOfRegion(int regionIndex) {
        
        double sum = calculateSelfAssoc(regionIndex);
        
        Region region = regions.get(regionIndex);
        
        Set<PairInt> points = region.getPoints();
        
        Set<PairInt> perimeter = region.getPerimeter();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        for (PairInt p : perimeter) {
            int x = p.getX();
            int y = p.getY();
            int idx1 = calculatePixelIndex(y, x);
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if (x2 < 0 || y2 < 0 || (x2 > (imageWidth - 1)) || (y2 > (imageHeight - 1))) {
                    continue;
                }
                PairInt p2 = new PairInt(x2, y2);
                if (points.contains(p2)) {
                    continue;
                }
                int idx2 = calculatePixelIndex(y2, x2);
                if (idx1 < idx2) {
                    sum += diffOrSim.get(idx1, idx2);
                } else {
                    sum += diffOrSim.get(idx2, idx1);
                }
            }
        }
        
        return sum;
    }
    
    /**
     * get the edge weights matrix of differences or similarity between 
     * image pixels.  Note that it has to be be accessed using convention
     * (smallestIndex, largerIndex).  Also note that the matrix is not copied,
     * so any modifications will be present in this instance's too.
     * 
     * @return 
     */
    public SimpleMatrix getEdgeMatrix() {
        return diffOrSim;
    }
}
