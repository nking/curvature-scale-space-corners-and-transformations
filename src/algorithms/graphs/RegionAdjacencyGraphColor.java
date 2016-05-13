package algorithms.graphs;

import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.segmentation.ColorSpace;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;

/**
 *
 * class to hold a list of region nodes and calculate an adjacency map
 * for them and pairwise color differences or similarity.
 * 
 * @author nichole
 */
public class RegionAdjacencyGraphColor extends RegionAdjacencyGraph {
    
    /*
       4     
       3
       2       
       1      \|/_
       0 |/_  \|/_\|
         0  1  2   3
    
    For each pixel in img, the color difference or similarity to it's adjacent 
    pixels are calculated.  
    For each pixel, using a pattern of only calculating the neighbors for
    offsets (+1,0), (+1,+1), (0,+1), (-1,+1) calculates the edge without repeating
    a pair.
    
    The conversion from col, row to a single index will use the convention in 
    Image.java which is index = (row * width) + col.
    
    Also note that the normalized cuts class using the sparse matrix diffOrSim
    needs a symmetric matrix, so the pairs are stored in both [i][j] and [j][i].
    
    */
    protected FlexCompRowMatrix diffOrSim = null;
    
    protected ColorSpace colorSpace = null;
    
    protected final ImageExt img;
    
    protected final List<NormalizedCutsNode> nodes;
    
     /**
     * constructor.  
     *        
     * @param img 
     * @param labels double array of labels for each pixel using the convention
     * labels[col][row]. Note that the labeled regions must be contiguous.
     */
    public RegionAdjacencyGraphColor(ImageExt img, int[][] labels) {

        super(img, labels);  
        
        this.img = img;

        /*
        may make some improvements:
        (if the statement in the normalized cuts paper regarding bit significance
       is followed, I should be able to factor the doubles to integers.
        in that case, can use compression as I do in Image.java and
        find a linear algebra eigen solver which will use the get and set methods
        of my image class to keep the data small.
        Also note, that instead of using an eigen solver, could alternatively,
        use my two-point clustering code for the remaining logic.  the distance transform
        and histograms make the remaining logic very fast, but it is dependent upon
        the bounds of the data (that is maximum values of data points which would be
        differences here converted to integers).
        */
        
        int n = regions.size();
        this.nodes = new ArrayList<NormalizedCutsNode>(n);
        for (int i = 0; i < n; ++i) {
            NormalizedCutsNode node = new NormalizedCutsNode(i);
            nodes.add(node);
        }
    }
    
    public RAGCSubGraph createANodesGraph() {
        RAGCSubGraph g = new RAGCSubGraph(nodes, diffOrSim,
        adjacencyMap);
        return g;
    }
    
    public int[][] relabelUsingNodes() {
        
        int n0 = labels.length;
        
        int[][] labels2 = new int[n0][];
        for (int i = 0; i < n0; ++i) {
            labels2[i] = new int[labels[i].length];
            for (int j = 0; j < labels[i].length; ++j) {
                int v = labels[i][j];
                labels2[i][j] = nodes.get(v).getnCutsLabel();
            }
        }
        
        return labels2;
    }
    
   
    public void populateEdgesWithColorSimilarity(ColorSpace clrSpace) {
        
        populatePairDifferences(clrSpace);

        double sigma;
        if (clrSpace.equals(ColorSpace.HSV)) {
            // max possible = sqrt(1*1 + 1*1 + 1*1) = sqrt3
            sigma = 1.7320508;

        } else if (clrSpace.equals(ColorSpace.CIELAB)) {
            // max possible deltaE2000 = 19.22
            sigma = 19.22;

        } else {
             // rgb
            // max possible = math.sqrt(3*255*255) = sqrt(3)*255
            sigma = 441.6729559;
        }
        
        Set<PairInt> added = new HashSet<PairInt>();
        
        for (Entry<Integer, Set<Integer>> entry : adjacencyMap.entrySet()) {
            Integer index1 = entry.getKey();
            int idx1 = index1.intValue();
            for (Integer index2 : entry.getValue()) {
                int idx2 = index2.intValue();
                assert(idx1 != idx2);
                PairInt p;
                if (idx1 < idx2) {
                    p = new PairInt(idx1, idx2);
                } else {
                    p = new PairInt(idx2, idx1);
                }
                if (added.contains(p)) {
                    continue;
                }
                added.add(p);
                
                // set both [i][j] and [j][i] to make symmetric matrix
                double d = diffOrSim.get(idx1, idx2);
                assert(d == diffOrSim.get(idx2, idx1));
                double similarity = Math.exp(-1*d*d/sigma);
                diffOrSim.set(idx1, idx2, similarity);
                diffOrSim.set(idx2, idx1, similarity);
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
    
    public int getImageWith() {
        return imageWidth;
    }
    public int getImageHeight() {
        return imageHeight;
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
        
        int nNodes = regions.size();
        int nCols = img.getWidth();
        int nRows = img.getHeight();

        diffOrSim = new FlexCompRowMatrix(nNodes, nNodes);
                
        // calculate the average colors for each node
        float[][] nodeColors = new float[3][nNodes];
        for (int i = 0; i < 3; ++i) {
            nodeColors[i] = new float[nNodes];
        }
        int[] count = new int[nNodes];
        
        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {
                int label = labels[col][row];
                nodeColors[0][label] += img.getR(col, row);
                nodeColors[1][label] += img.getG(col, row);
                nodeColors[2][label] += img.getB(col, row);
                count[label]++;
            }
        }
        for (int i = 0; i < nNodes; ++i) {
            nodeColors[0][i] /= (float)count[i];
            nodeColors[1][i] /= (float)count[i];
            nodeColors[2][i] /= (float)count[i];
        }
        
        // calculate pairwise differences for adjacent nodes (similarity or distance edgeS)
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        Set<PairInt> added = new HashSet<PairInt>();
        
        for (Entry<Integer, Set<Integer>> entry : adjacencyMap.entrySet()) {
            Integer index1 = entry.getKey();
            int idx1 = index1.intValue();
            for (Integer index2 : entry.getValue()) {
                int idx2 = index2.intValue();
                assert(idx1 != idx2);
                PairInt p;
                if (idx1 < idx2) {
                    p = new PairInt(idx1, idx2);
                } else {
                    p = new PairInt(idx2, idx1);
                }
                if (added.contains(p)) {
                    continue;
                }
                added.add(p);
                
                double diff;
                if (clrSpace.equals(ColorSpace.HSV)) {
                    float sumDiff = 0;
                    for (int m = 0; m < 3; ++m) {
                        float d = nodeColors[m][idx1] - nodeColors[m][idx2];
                        sumDiff += (d * d);
                    }
                    diff = Math.sqrt(sumDiff);
                } else if (clrSpace.equals(ColorSpace.CIELAB)) {
                    diff = Math.abs(cieC.calcDeltaECIE2000(
                        nodeColors[0][idx1], nodeColors[1][idx1], nodeColors[2][idx1],
                        nodeColors[0][idx2], nodeColors[1][idx2], nodeColors[2][idx2]
                    ));
                } else {
                    //RGB
                    float diffR = nodeColors[0][idx1] - nodeColors[0][idx2];
                    float diffG = nodeColors[1][idx1] - nodeColors[1][idx2];
                    float diffB = nodeColors[2][idx1] - nodeColors[2][idx2];
                    diff = Math.sqrt(diffR * diffR + diffG * diffG + diffB * diffB);
                }
                // set both [i][j] and [j][i] to make matrix symmetric
                diffOrSim.set(idx1, idx2, diff);
                diffOrSim.set(idx2, idx1, diff);
            }
        }                                       
    }
    
    /*
       calculate the normalized cute between regions regionIndex1 and
       regionIndex2 using
       normalized_cut(A, B) = 2 - normalized_assoc(A, B);
     
                               assoc(A, A) 
         normalized_assoc(A) = ----------- 
                               assoc(A, V)
                               
         where assoc(A, A) is the sum of all edges within region A
         and assoc(A, V) is the sum of all edges from A to another region.
    */
    
    public double[] calculateCentroidOfRegion(int regionIndex) {
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        return curveHelper.calculateXYCentroids(regions.get(regionIndex).getPoints());
    }
    
    /**
     * get the edge weights matrix of differences or similarity between 
     * image pixels as a sparse symmetric matrix. 
     * @return 
     */
    public FlexCompRowMatrix getEdgeMatrix() {
        return diffOrSim;
    }
    
    public Set<PairInt> getRegionPoints(int regionIndex) {
        return regions.get(regionIndex).getPoints();
    }
    
}
