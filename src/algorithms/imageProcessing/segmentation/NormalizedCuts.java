package algorithms.imageProcessing.segmentation;

import algorithms.graphs.NormalizedCutsNode;
import algorithms.graphs.RAGCSubGraph;
import algorithms.graphs.RegionAdjacencyGraphColor;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.util.MatrixUtil;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Logger;
import no.uib.cipr.matrix.DenseVectorSub;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.VectorEntry;
import no.uib.cipr.matrix.sparse.ArpackSym;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;

/**
 * A hierarchical segmentation method that determines graph cuts
 * which maximize the difference between regions and minimize the differences
 * within a region using
 * cuts normalized by associations between adjacent regions.
 * 
 * following the algorithm of
 * "Normalized Cuts and Image Segmentation"
        by Jianbo Shi and Jitendra Malik
         
    Most of the structure additionally follows the scipy skimage implementation
    of normalized cuts.
    the skimage API license is at
    https://github.com/scikit-image/scikit-image/blob/master/LICENSE.txt  
    (copy that here)
    * 
 * @author nichole
 */
public class NormalizedCuts {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * using the recursive 2-way Ncut pattern in normalized cuts 
     * to segment the image into regions.
     * @param img
     * @param labels of contiguous pixels
     * @return 
     */
    public int[] normalizedCut(ImageExt img, int[] labels) {
       
        //TODO: note, as authors mention, the edge weights with values > 0.01
        // are significant and the remaining are zeroes of different precision
        
        log.fine("input labels=" + Arrays.toString(labels));
        
        RegionAdjacencyGraphColor rag = new RegionAdjacencyGraphColor(
            img, labels);

        log.fine("rag.nNodes=" + rag.getNumberOfRegions() + " at start");

        ColorSpace colorSpace = ColorSpace.RGB;
                        
        rag.populateEdgesWithColorSimilarity(colorSpace);
        
        RAGCSubGraph nodesGraph = rag.createANodesGraph();
        
        // add an edge to self for every node, weight = edge_max (which is 1.0 if not specified)
        // only doing this in the similarity matrix, not the graph adjacency map
        FlexCompRowMatrix w = nodesGraph.getEdgeMatrix();
        int n = nodesGraph.getNumberOfNodes();
        for (int i = 0; i < n; ++i) {
            w.set(i, i, 1.0);
        }
        
        performNormalizedCuts(nodesGraph);
        
        return rag.relabelUsingNodes();
    }
    
     /**
      * Perform a 2-way normalized cut recursively using the graph's similarity
      * edges with a lancos eigen value decomposition to find the smallest
      * normalized cut each time and labeling the nCuts field in nodes as a result.
      * 
      * @param rag 
      */
    private void performNormalizedCuts(RAGCSubGraph graph) {
        
        // number of normalizd cuts to perform before determinging optimal among them
        int numCuts = 10;
        
        nCutRelabel(graph, numCuts);        
    }

    /**
     * perform 2-way recursive splitting of graph
     * @param regionsList
     * @param rag
     * @param numCuts 
     */
    private void nCutRelabel(RAGCSubGraph graph, int numCuts) {

        double thresh = 0.06;
        
        FlexCompRowMatrix w = graph.getEdgeMatrix();
        FlexCompRowMatrix d = createD(w, graph);
        
        log.fine("w.nnz=" + MatrixUtil.countNodes(w));
        
        if (w.numRows()> 2) {
        
            FlexCompRowMatrix d2 = (FlexCompRowMatrix) d.copy();
            Iterator<MatrixEntry> iter = d.iterator();
            while (iter.hasNext()) {
                MatrixEntry entry = iter.next();
                int col = entry.column();
                int row = entry.row();
                double v = entry.get();
                v = (v == 0) ? Double.MAX_VALUE : 1./Math.sqrt(v);
                d2.set(col, row, v);
            }
        
            //D^(-1/2)(D - W)D^(-1/2)z = lamdba z;
            // tmp = diag - w
            // tmp = d2 mult tmp mult d2
            FlexCompRowMatrix tmp = MatrixUtil.sparseMatrixSubtract(d, w);
            tmp = MatrixUtil.sparseMatrixMultiply(d2, tmp);        
            tmp = MatrixUtil.sparseMatrixMultiply(tmp, d2);

            // PRECISION CORRECTIONS needed for perfectly symmetric matrix
            FlexCompRowMatrix tmp2 = (FlexCompRowMatrix) tmp.copy();
            iter = tmp2.iterator();
            while (iter.hasNext()) {
                MatrixEntry entry = iter.next();
                int col = entry.column();
                int row = entry.row();
                double v = entry.get();
                tmp.set(col, row, v);
                tmp.set(row, col, v);
            }
                            
            // 2nd smallest eigenvector
            int m = w.numRows();
            //int nEig = 2;
            int nEig = Math.min(100, m - 2);

            log.fine("nEigen =" + nEig + " eig input=" + tmp.toString());

            /*
            info status:
                0 = Normal Exit
                1 = Maximum number of iterations taken
                2 = 
            */
//TODO: when have a return tpe,
// catch errors here and return a default
            ArpackSym arpackSym = new ArpackSym(tmp);
            Map<Double, DenseVectorSub> rMap = arpackSym.solve(nEig, ArpackSym.Ritz.SM);

            double secondSmallestEigenValue = find2ndSmallestEigenValue(rMap);
            DenseVectorSub eigenVector = rMap.get(Double.valueOf(secondSmallestEigenValue));
            
            if (eigenVector != null) {
           
                log.finest("secondSmallestEigenValue=" + secondSmallestEigenValue);

                log.finest("eigen vector size=" + eigenVector.size());

                log.finest("eigenVector=" + eigenVector);

                MinCut minCut = getMinNCut(eigenVector, d, w, numCuts);

                if (minCut != null) {

                    log.fine("mCut=" + minCut.mCut);
                    log.fine("cut=" + Arrays.toString(minCut.minMask));

                    if (minCut.mCut < thresh) {

                        RAGCSubGraph[] subGraphs = graph.partition(minCut.minMask);

                        log.fine("len(sub1)=" + subGraphs[0].getNumberOfNodes());

                        log.fine("len(sub2)=" + subGraphs[1].getNumberOfNodes());

                        nCutRelabel(subGraphs[1], numCuts);
                        nCutRelabel(subGraphs[0], numCuts);

                        return;
                    }
                }
            }
        }
        
        /* from scipy:
        # The N-cut wasn't small enough, or could not be computed.
        # The remaining graph is a region.
        # Assign `ncut label` by picking any label from the existing nodes, since
        # `labels` are unique, `new_label` is also unique.
        */
        int label = graph.getNodes().get(0).getLabel();
        for (NormalizedCutsNode node : graph.getNodes()) {
            node.setNCutsLabel(label);
        } 
    }

    /**
     * Use the eigenvector with the second smallest eigenvalue to bipartition 
     * the graph by finding the splitting point such that Ncut is minimized.
     * @param graph
     * @param eigenVector
     * @param numCuts 
     */
    private MinCut getMinNCut(
        DenseVectorSub eigenVector, FlexCompRowMatrix d, FlexCompRowMatrix w,
        int numCuts) {
        
        double minEV = Double.MAX_VALUE;
        double maxEV = Double.MIN_VALUE;
        for (VectorEntry entry : eigenVector) {
            double v = entry.get();
            if (minEV > v) {
                minEV = v;
            }
            if (maxEV < v) {
                maxEV = v;
            }
        }
        
        /*
        "In the ideal case, the eigenvector should only take on two discrete 
        values and the signs of the values can tell us exactly how to partition 
        the graph. However, our eigenvectors can take on continuous values and 
        we need to choose a splitting point to partition it into two parts. 
        There are many different ways of choosing such a splitting point. 
           One can take 0 or the median value as the splitting point 
           or one can search for the splitting point such that the resulting 
        partition has the best N cut(A, B) value. 
        We take the latter approach in our work. 
        *Currently, the search is done by checking l evenly spaced 
        possible splitting points, and computing the best Ncut among them. 
        In our experiments, the values in the eigenvectors are usually well 
        separated and this method of choosing a splitting point is very reliable 
        even with a small l.  
        * After the graph is broken into two pieces, we can recursively run our 
        algorithm on the two partitioned parts. Or, equivalently, we could take 
        advantage of the special properties of the other top eigenvectors as 
        explained in the previous section to subdivide the graph based on those 
        eigenvectors. The recursion stops once the Ncut value exceeds certain 
        limit.
        We also impose a stability criterion on the partition. 
        As we saw earlier, and as we see in the eigenvectors with the seventh 
        to ninth smallest eigenvalues (Fig. 3g-h), sometimes an eigenvector 
        can take on the shape of a continuous function, rather that the 
        discrete indicator function that we seek. From the view of segmentation, 
        such an eigenvector is attempting to subdivide an image region where 
        there is no sure way of breaking it. In fact, if we are forced to 
        partition the image based on this eigenvector, we will see there are 
        many different splitting points which have similar Ncut values. Hence, 
        
        the partition will be highly uncertain and unstable. In our current 
        segmentation scheme, we simply choose to ignore all those eigenvectors 
        which have smoothly varying eigenvector values. We achieve this by 
        imposing a stability criterion which measures the degree of smoothness 
        in the eigenvec- tor values. The simplest measure is based on first 
        computing the histogram of the eigenvector values and then computing 
        the ratio between the minimum and maximum values in the bins. When the 
        eigenvector values are continuously varying, the values in the 
        histogram bins will stay relatively the same and the ratio will be 
        relatively high. In our experiments, we find that simple thresholding 
        on the ratio described above can be used to exclude unstable 
        eigenvectors. We have set that value to be 0.06 in all our experiments.
        Fig. 4 shows the final segmentation for the image shown in Fig. 2.
        */
        
        // 1e-5 * abs(b) + 1e-8
        if ((maxEV - minEV) < ((1e-5)*Math.abs(minEV) + 1e-8)) {
            return null;
        }
        
        //min_mask = np.zeros_like(ev, dtype=np.bool)
        boolean[] minMask = null;
        
        double mCut = Double.MAX_VALUE;
        
        double delta = (maxEV - minEV)/(double)numCuts;
        double t = minEV;
        while (t < maxEV) {
            
            boolean[] mask = new boolean[eigenVector.size()];
            for (int i = 0; i < mask.length; ++i) {
                mask[i] = (eigenVector.get(i) > t) ? true : false;
            }
            
            double cost = nCutCost(mask, d, w);
            
            if (Double.isNaN(cost)) {
                // happens when mask has all same values
                break;
            }
            
            log.fine("cut.length=" + mask.length + " cost=" + cost + 
                "  cut_mask=" + Arrays.toString(mask));
            
            if (cost < mCut) {
                mCut = cost;
                minMask = mask;
            }
            
            t += delta;
        }
       
        MinCut m = new MinCut();
        m.minMask = minMask;
        m.mCut = mCut;
    
        return m;
    }
    
    private class MinCut {
        double mCut;
        boolean[] minMask;
    }

    /**
     <pre>
     normalized cut is the measure of the disassociation
                                   cut(A, B)      cut(A,B)
           normalized_cut(A, B) = -----------  + -----------
                                  assoc(A, V)    assoc(B, V)
                                  
      where
           cut(A, B) is the sum of all edges connecting set A with set B
      and
           assoc(A, V) is the sum of all of A's edge's to any node
      and 
           assoc(A, A) is the sum of all edges within A alone and not to any
              other regions
            
      the normalized_cut can be rewritten using similarity:
                                    assoc(A, A)   assoc(B, B)
           normalized_assoc(A, B) = ----------- + ------------
                                    assoc(A, V)   assoc(B, V)
                                    
           normalized_cut(A, B) = 2 - normalized_assoc(A, B)
           
      </pre>
     * @param mask
     * @param d
     * @param w
     * @return 
     */
    private double nCutCost(boolean[] cut, FlexCompRowMatrix d, FlexCompRowMatrix w) {

        double cutCost = 0;
        for (MatrixEntry entry : w) {
            int row = entry.row();
            int col = entry.column();
            // it's a summetric matrix so only count edges once
            if (row <= col) {
                if (cut[col] != cut[row]) {
                    cutCost += entry.get();
                }
            }
        } 
        
        double sumAssocA = 0;
        double sumAssocB = 0;
        for (MatrixEntry entry : d) {
            int i = entry.row();
            double v = entry.get();
            if (cut[i]) {
                sumAssocA += v;
            } else {
                sumAssocB += v;
            }
        }
        
        return (cutCost / sumAssocA) + (cutCost / sumAssocB);
    }
    
    private double find2ndSmallestEigenValue(Map<Double, DenseVectorSub> eig) {
        
        int n = eig.size();
        if (n < 2) {
            return -1;
        }
        
        double min1 = Double.MAX_VALUE;
        double min2 = Double.MAX_VALUE;
        
        for (Entry<Double, DenseVectorSub> entry : eig.entrySet()) {
            double v = entry.getKey().doubleValue();
            if (min1 == Double.MAX_VALUE) {
                min1 = v;
            } else if (min2 == Double.MAX_VALUE) {
                if (v <= min1) {
                    min2 = min1;
                    min1 = v;
                } else {
                    min2 = v;
                }
            } else {
                if (v <= min1) {
                    min2 = min1;
                    min1 = v;
                } else if (v < min2) {
                    min2 = v;
                }
            }
        }
        
        return min2;
    }        
    
    private FlexCompRowMatrix createD(FlexCompRowMatrix w, RAGCSubGraph graph) {

        //D is an N X N diagonal matrix with d on the diagonal
        //    d(i) = summation over j of w(i, j) where w is "weight" of the edge
        //    and j is over all nodes
        FlexCompRowMatrix d = new FlexCompRowMatrix(w.numRows(), w.numColumns());
        
        for (MatrixEntry entry : w) {
            int col = entry.column();
            double v0 = entry.get();
            double v1 = d.get(col, col);
            d.set(col, col, (v0 + v1));
        }
        
        return d;
    }
    
}
