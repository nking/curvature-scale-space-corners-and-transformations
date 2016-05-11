package algorithms.graphs;

import algorithms.compGeometry.clustering.KMeansPlusPlus;
import algorithms.imageProcessing.DFSConnectedGroupsFinder;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.segmentation.ColorSpace;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseVectorSub;
import no.uib.cipr.matrix.MatrixEntry;
import no.uib.cipr.matrix.sparse.ArpackSym;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author nichole
 */
public class RegionAdjacencyGraphColorTest extends TestCase {
    
    public RegionAdjacencyGraphColorTest() {
    }
    
    public void estEig() throws Exception {
        
        /*
        no.uib.cipr.matrix
        */
        
        /*
        looking at LA4J tests for eidgen decomposition and
        comparing results to EJML to see how to extract the vectors and
        values
        */
        double[][] input = new double[][]{
                {16.0, -11.0, 99.0},
                {7.0, -2.0, -42.0},
                {8.0, -1.0, -7.0}
        };
        double[][] inputCp = new double[][]{
                {16.0, -11.0, 99.0},
                {7.0, -2.0, -42.0},
                {8.0, -1.0, -7.0}
        };
        double[][][] output = new double[][][]{
                {
                        {0.982, 0.644, -0.353},
                        {-0.026, -0.896, -2.320},
                        {0.187, -0.359, -0.186}
                },
                {
                        {35.149, 0.0, 0.0},
                        {0.0, -23.856, 0.0},
                        {0.0, 0.0, -4.293}
                }
        };
        
        SimpleMatrix ejmlM = new SimpleMatrix(input);
        org.ejml.interfaces.decomposition.EigenDecomposition ejmlEVD = ejmlM.eig().getEVD();
        for (int i = 0; i < ejmlEVD.getNumberOfEigenvalues(); ++i) {
            System.out.println("EJML eigenvector i=" + i + " " + 
                ejmlEVD.getEigenVector(i));
        }
        for (int i = 0; i < ejmlEVD.getNumberOfEigenvalues(); ++i) {
            System.out.println("EJML eigenvalue i=" + i + " " + 
                ejmlEVD.getEigenvalue(i));
        }
        
         /*
        double[][][] output = new double[][][]{
                {
                        {0.982, 0.644, -0.353},
                        {-0.026, -0.896, -2.320},
                        {0.187, -0.359, -0.186}
                },
                {
                        {35.149, 0.0, 0.0},
                        {0.0, -23.856, 0.0},
                        {0.0, 0.0, -4.293}
                }
        };
        */       
    }
    
    public void testLinAlg() throws IOException, Exception {
        
        String fileName = "tmp3.png";
        //String fileName = "android_statues_02.jpg";
                 
        System.out.println("fileName=" + fileName);
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        String fileNameRoot = fileName.substring(0, fileName.lastIndexOf("."));
            
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        // use a faster cluster code here for inital seeds
        // implement the SLIC algorithm
        
        int k = 8;
        KMeansPlusPlus kmpp = new KMeansPlusPlus();
        kmpp.computeMeans(k, img.copyToGreyscale());
        int[] pixAssignments = kmpp.getImgPixelSeedIndexes();
        
        System.out.println("have initial labels (super pixels)");
        System.out.flush();
        
        int w = img.getWidth();
        int h = img.getHeight();
        //int n = img.getNPixels();
        
        int nMaxLabel = Integer.MIN_VALUE;
        for (int i = 0; i < pixAssignments.length; ++i) {
            if (pixAssignments[i] > nMaxLabel) {
                nMaxLabel = pixAssignments[i];
            }
        }
        List<Set<PairInt>> pixelSets = new ArrayList<Set<PairInt>>();
        for (int i = 0; i <= nMaxLabel; ++i) {
            pixelSets.add(new HashSet<PairInt>());
        }
        for (int i = 0; i < pixAssignments.length; ++i) {
            int label = pixAssignments[i];
            int col = img.getCol(i);
            int row = img.getRow(i);
            pixelSets.get(label).add(new PairInt(col, row));
        }
        
        List<Set<PairInt>> contigPixelSets = new ArrayList<Set<PairInt>>();
        for (int i = 0; i < pixelSets.size(); ++i) {
            Set<PairInt> set = pixelSets.get(i);
            DFSConnectedGroupsFinder finder = new DFSConnectedGroupsFinder();
            finder.setMinimumNumberInCluster(1);
            finder.findConnectedPointGroups(set);
            for (int j = 0; j < finder.getNumberOfGroups(); ++j) {
                Set<PairInt> g = finder.getXY(j);
                contigPixelSets.add(g);
            }
        }
        
        int[][] labels = new int[w][h];
        for (int i = 0; i < w; ++i) {
            labels[i] = new int[h];
        }
        for (int i = 0; i < contigPixelSets.size(); ++i) {
            Set<PairInt> set = contigPixelSets.get(i);
            for (PairInt p : set) {
                labels[p.getX()][p.getY()] = i;
            }
        }
        
        RegionAdjacencyGraphColor rag = new RegionAdjacencyGraphColor(img, labels);
        rag.populateEdgesWithColorSimilarity(ColorSpace.HSV);
        
        System.out.println("created region agency graph and similarity matrix");
        System.out.flush();
        FlexCompRowMatrix weights = rag.getEdgeMatrix();
        FlexCompRowMatrix diag = createD(weights, rag);

        FlexCompRowMatrix d2 = (FlexCompRowMatrix) diag.copy();
        Iterator<MatrixEntry> iter = diag.iterator();
        while (iter.hasNext()) {
            MatrixEntry entry = iter.next();
            int col = entry.column();
            int row = entry.row();
            double v = entry.get();
            v = (v == 0) ? Double.MAX_VALUE : 1./Math.sqrt(v);
            d2.set(col, row, v);
        }
        
        System.out.println("weight matrix " + 
            " nRowsXCols=" + (weights.numRows() * weights.numColumns())
            + " number of nodes=" + MatrixUtil.countNodes(weights)
        );
        
        //D^(-1/2)(D - W)D^(-1/2)z = lamdba z;
        // tmp = diag - w
        // tmp = d2 mult tmp mult d2
        FlexCompRowMatrix tmp = MatrixUtil.sparseMatrixSubtract(diag, weights);
        
        System.out.println("tmp matrix after subtract " + 
            " nRowsXCols=" + (tmp.numRows() * tmp.numColumns())
            + " number of nodes=" + MatrixUtil.countNodes(tmp)
        );
        
        tmp = MatrixUtil.sparseMatrixMultiply(d2, tmp);
        
        System.out.println("tmp matrix after multiply " + 
            " nRowsXCols=" + (tmp.numRows() * tmp.numColumns())
            + " number of nodes=" + MatrixUtil.countNodes(tmp)
        );
        
        tmp = MatrixUtil.sparseMatrixMultiply(tmp, d2);
 
        // PRECISION CORRECTIONS needed for perfectly symmetric matrix
        FlexCompRowMatrix tmp2 = new FlexCompRowMatrix(tmp.numRows(), tmp.numRows());
        //matrix is not symetric
        iter = tmp.iterator();
        while (iter.hasNext()) {
            MatrixEntry entry = iter.next();
            int col = entry.column();
            int row = entry.row();
            double v = entry.get();
            tmp2.set(col, row, v);
            tmp2.set(row, col, v);
        }
       
        int m = weights.numRows();
        int nEig = Math.min(100, m - 2);
        ArpackSym arpackSym = new ArpackSym(tmp2);
        Map<Double, DenseVectorSub> rMap = arpackSym.solve(nEig, ArpackSym.Ritz.SM);
        
        assertNotNull(rMap);
        assertTrue(rMap.size() > 1);
        /*for (Map.Entry<Double, DenseVectorSub> result : rMap.entrySet()) {           
            System.out.println("resulting eigenvalue=" + result.getKey().toString());
            System.out.println("resulting eigenvector=" + result.getValue().toString());
        }*/
        
        /*
        SimpleMatrix tmp = diag.minus(w);
        tmp = d2.mult(tmp).mult(d2);

        SimpleEVD eigVD = tmp.eig();
        
        System.out.println("n eigenvalues=" + eigVD.getNumberOfEigenvalues());
        */
    }
    
    private FlexCompRowMatrix createD(FlexCompRowMatrix w, RegionAdjacencyGraphColor rag) {

        //D is an N X N diagonal matrix with d on the diagonal
        //    d(i) = summation over j of w(i, j) where w is "weight" of the edge
        //    and j is over all nodes
        FlexCompRowMatrix d = new FlexCompRowMatrix(w.numRows(), w.numColumns());
        
        for (int i = 0; i < w.numRows(); ++i) {
            
            //d(i) = summation over j of w(i, j)
            
            Integer index1 = Integer.valueOf(i);
            
            Set<Integer> indexes2 = rag.getAdjacentIndexes(index1);
            
            assert(indexes2 != null);
        
            int idx1 = index1.intValue();
            
            double dSum = 0;
            
            for (Integer index2 : indexes2) {
                int idx2 = index2.intValue();
                assert(idx1 != idx2);
                dSum += w.get(idx1, idx2);
            }
            
            d.set(i, i, dSum);
        }
        
        return d;
    }
   
}
