package algorithms.imageProcessing.segmentation;

import algorithms.graphs.NormalizedCutsNode;
import algorithms.graphs.RAGCSubGraph;
import algorithms.graphs.RegionAdjacencyGraphColor;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.MiscMath;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
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
TODO: one day, revise this code to follow the paper implementation more closely
      (e.g. using the Lanczos method based on power method) and improvements since the publication.

 * A hierarchical segmentation method that determines graph cuts
 * which maximize the difference between regions and minimize the differences
 * within a region using
 * cuts normalized by associations between adjacent regions.
 * 
 * implementation of the algorithm of
 * "Normalized Cuts and Image Segmentation"
        by Jianbo Shi and Jitendra Malik

    grouping algorithm:
      1. Given an image or image sequence, set up a
        weighted graph G(V; E) and set the weight on
        the edge connecting two nodes to be a measure of
        the similarity between the two nodes.
      2. Solve (D - W)*x = lambda*D*x for eigenvectors with the
        smallest eigenvalues.
        (NOTE from mmds.org book chap 10, is that since the Laplacian is symmetric, the
         2nd smallest value is the minimum of x^T * L * x where x is [x1, x2, . . . , xn] 
         is a column vector with n components, and the minimum is taken under the constraints:
         the length of x = 1 (ie. sum of (x_i)^2 = 1), and x is orthogonal to the
         eigenvector associated with the smallest eigenvalue (which is eigenvalue=0 and eigenvector
         is all 1's).
      3. Use the eigenvector with the second smallest
        eigenvalue to bipartition the graph.
        NOTE: The second smallest eigenvalue of D-W is sometimes known as the Fiedler value.
      4. Decide if the current partition should be subdivided
        and recursively repartition the segmented parts if
        necessary

         
    Much of the logic here is adopted from the scipy skimage 
    implementation of normalized cuts.
    * https://github.com/scikit-image/scikit-image/blob/master/skimage/future/graph/
    The skimage API license is at
    https://github.com/scikit-image/scikit-image/blob/master/LICENSE.txt  
    and is:
    ---begin scipy license ----
    Copyright (C) 2011, the scikit-image team
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    1. Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
    
    2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in
    the documentation and/or other materials provided with the
    distribution.
    
    3. Neither the name of skimage nor the names of its contributors may be
    used to endorse or promote products derived from this software without
    specific prior written permission.
    
    THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
    IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
    INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
    IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.

    --- end scipy license ----

    <pre>
    useful documentation for netlib-java and matrix-toolkit-java (mtj) is
    https://github.com/fommil/netlib-java
    https://github.com/fommil/matrix-toolkits-java/wiki
    https://github.com/fommil/matrix-toolkits-java/blob/master/src/main/java/no/uib/cipr/matrix/sparse/ArpackSym.java
    http://www.javadoc.io/doc/com.googlecode.matrix-toolkits-java/mtj/1.0.4
    
    their license information is in file lib/NETLIB-JAVA.txt
    </pre>

 * @author nichole
 */
public class NormalizedCuts {
    
    //TODO: needs to be edited to use spatial location and adjacency
    //   map.  look at DNCuts
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    // number of normalized cuts to perform before 
    // determining optimal among them.
    // cannot be smaller than 2
    private int numCuts = 3;//5;

    private ColorSpace colorSpace = ColorSpace.RGB;
    
    private ColorOption colorOption = ColorOption.RGB;
    
    private enum ColorOption {
        RGB, HSV, LOW_THRESHOLD_RGB, LOW_THRESHOLD_HSV, HSV_COLOR_HISTOGRAMS,
        CIELAB, POLAR_CIELAB
    }
    
    private double thresh = 0.06;
    double sigma = 22;

    public void setThreshold(double thresh) {
        this.thresh = thresh;
    }
    public void setSigma(double sigma) {
        this.sigma = sigma;
    }
    public void setToLowThresholdRGB() {
        colorOption = ColorOption.LOW_THRESHOLD_RGB;
        colorSpace = ColorSpace.RGB;
        // for superpixels w/ n=200
        thresh = 1.e-16;
        sigma = 1.7320508;
    }
    
    /**
     * NOTE: pre-processing to filter out largest discrepant
     * color histograms (smallest intersection values)
     * may improve results.
     */
    public void setToColorHistogramsOfHSV() {
        colorOption = ColorOption.HSV_COLOR_HISTOGRAMS;
        colorSpace = ColorSpace.HSV;
        thresh = 0.5;
        sigma = 1.732;
    }
    
    public void setToLowThresholdHSV() {
        colorOption = ColorOption.LOW_THRESHOLD_HSV;
        colorSpace = ColorSpace.HSV;
        // for superpixels w/ n=200  thresh=0.05, sigma=0.1
        thresh = 5e-6; //1e-13 and 0.001 or 5e-6 and 0.0025
        sigma = 0.0025; //change smaller
    }
   
    public void setColorSpaceToRGB() {
        colorOption = ColorOption.RGB;
        colorSpace = ColorSpace.RGB;
        // for superpixels w/ n=200 thresh=0.06 sigma=22
        thresh = 0.06;
        sigma = 22;
    }
    
    public void setColorSpaceToHSV() {
        colorOption = ColorOption.HSV;
        colorSpace = ColorSpace.HSV;
        // for superpixels w/ n=200  thresh=0.05, sigma=0.1
        thresh = 0.05;
        sigma = 0.1;
    }
    
    public void setColorSpaceToCIELAB() {
        colorOption = ColorOption.CIELAB;
        colorSpace = ColorSpace.CIELAB;
        // for superpixels w/ n=200 
        thresh = 0.06; 
        sigma = 10;
    }
    
    /**
     * NOTE: this one especially, should improve when location
     * and adjacency maps are used in the similarity score
     */
    public void setColorSpaceToPolarCIELAB() {
        colorOption = ColorOption.POLAR_CIELAB;
        colorSpace = ColorSpace.POLAR_CIELAB;
        // for superpixels w/ n=200 
        thresh = 0.01; 
        sigma = 4.;//4: 5
    }
    
    /**
     * using the recursive 2-way Ncut pattern in normalized cuts 
     * to segment the image into regions.
     * works well when the labels are large segmentation
     * cells such as the results of a "super pixels" algorithm.
     * 
     * @param img
     * @param labels of contiguous pixels.  note that label value range is 
     * compressed from minimum value to minimum value plus number of values
     * internally, but given array is not modified.
     * @return 
     */
    public int[] normalizedCut(ImageExt img, int[] labels) {
       
        // compress labels
        int[] labels2 = Arrays.copyOf(labels, labels.length);
        LabelToColorHelper.condenseLabels(labels2);
        labels = labels2;
        
        //TODO: note, as authors mention, the edge weights with values > 0.01
        // are significant and the remaining are zeroes of different precision
        
        log.fine("input labels=" + Arrays.toString(labels));
        
        RegionAdjacencyGraphColor rag = new RegionAdjacencyGraphColor(
            img, labels);

        log.info("rag.nNodes=" + rag.getNumberOfRegions() + " at start");
         
        if (colorOption.equals(ColorOption.LOW_THRESHOLD_HSV)) {
            rag.populateEdgesWithLowThreshHSVSimilarity(sigma);
        } else if (colorOption.equals(ColorOption.LOW_THRESHOLD_RGB)) {
            rag.populateEdgesWithLowThreshRGBSimilarity(sigma);
        } else if (colorOption.equals(ColorOption.HSV_COLOR_HISTOGRAMS)) {
            rag.populateEdgesWithHSVColorHistogramSimilarity(sigma);
        } else {
            rag.populateEdgesWithColorSimilarity(colorSpace,
                sigma);
        }
        
        RAGCSubGraph nodesGraph = rag.createANodesGraph();

        // TODO: should normalize matrix data before use to zero mean and unit variance.
        
        // add an edge to self for every node, weight = edge_max (which is 1.0 if not specified)
        // only doing this in the similarity matrix, not the graph adjacency map
        FlexCompRowMatrix w = nodesGraph.getEdgeMatrix();
        int n = nodesGraph.getNumberOfNodes();
        for (int i = 0; i < n; ++i) {
            w.set(i, i, 1.0);
        }
        
        performNormalizedCuts(nodesGraph);
        
        int[] labeled = rag.relabelUsingNodes();
        
        labelTheUnlabeled(labeled);
        
        return labeled;
    }
    
     /**
      * Perform a 2-way normalized cut recursively using the graph's similarity
      * edges with a lancos eigen value decomposition to find the smallest
      * normalized cut each time and labeling the nCuts field in nodes as a result.
      * 
      * @param rag 
      */
    private void performNormalizedCuts(RAGCSubGraph graph) {
        
        nCutRelabel(graph);        
    }

    /**
     * perform 2-way recursive splitting of graph
     * @param regionsList
     * @param rag
     */
    private void nCutRelabel(RAGCSubGraph graph) {
        
        FlexCompRowMatrix w = graph.getEdgeMatrix();
        if (w.numRows() < 5) {
            return;
        }
        
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
            //NOTE: improve upon this one day, following "splitting point" of page 892:3 of paper
            //  and K-way cut on page 893
            // NOTE: The second smallest eigenvalue of D-W is sometimes known as the Fiedler value.
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
            Map<Double, DenseVectorSub> rMap = null;
            try {
                rMap = arpackSym.solve(nEig, ArpackSym.Ritz.SM);
            } catch (Throwable t) {
                log.severe("ERROR: " +  "array size=" +
                    " rows=" + tmp.numRows() +
                    " cols=" + tmp.numColumns()
                    + " => " + t.getMessage()
                );
                return;
            }

            double secondSmallestEigenValue = find2ndSmallestEigenValue(rMap);
            DenseVectorSub eigenVector = rMap.get(Double.valueOf(secondSmallestEigenValue));
            
            if (eigenVector != null) {
           
                log.finest("secondSmallestEigenValue=" + secondSmallestEigenValue);

                log.finest("eigen vector size=" + eigenVector.size());

                log.finest("eigenVector=" + eigenVector);

                MinCut minCut = getMinNCut(eigenVector, d, w);

                if (minCut != null) {

                    log.fine("mCut=" + minCut.mCut);
                    log.fine("cut=" + Arrays.toString(minCut.minMask));
                    log.info("mCut=" + minCut.mCut + " thresh=" + thresh);
                    
                    if (minCut.mCut < thresh) {

                        RAGCSubGraph[] subGraphs = graph.partition(minCut.minMask);

                        log.fine("len(sub1)=" + subGraphs[0].getNumberOfNodes());

                        log.fine("len(sub2)=" + subGraphs[1].getNumberOfNodes());

                        nCutRelabel(subGraphs[1]);
                        nCutRelabel(subGraphs[0]);

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
     */
    private MinCut getMinNCut(
        DenseVectorSub eigenVector, FlexCompRowMatrix d, 
        FlexCompRowMatrix w) {
        
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
    
    /**
     * this is a fudge to assign unlabeled pixels with a
     * value higher than the largest labeled.  it may
     * need to be modified or the reason for unlabeled to
     * be corrected upstream.
     * @param labels 
     */
    private void labelTheUnlabeled(int[] labels) {
        
        int maxLabel = Integer.MIN_VALUE;
        TIntSet unlabeled = new TIntHashSet();
        for (int i = 0; i < labels.length; ++i) {
            if (labels[i] == -1) {
                unlabeled.add(i);
            }
            if (labels[i] > maxLabel) {
                maxLabel = labels[i];
            }
        }
        if (unlabeled.size() > 0) {
            maxLabel++;
            TIntIterator iter = unlabeled.iterator();
            while (iter.hasNext()) {
                int idx = iter.next();
                labels[idx] = maxLabel; 
            }
        }
    }
}
