package algorithms.graphs;

import algorithms.compGeometry.clustering.KMeansPlusPlus;
import algorithms.compGeometry.clustering.KMeansPlusPlusColor;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.segmentation.ColorSpace;
import algorithms.misc.Misc;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import junit.framework.TestCase;
import no.uib.cipr.matrix.sparse.FlexCompColMatrix;
import org.ejml.simple.SimpleEVD;
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
        /*
        DefaultSparseRowDoubleMatrix2D ujmpM = new DefaultSparseRowDoubleMatrix2D(inputCp.length, inputCp[0].length);
        for (int i = 0; i < inputCp.length; ++i) {
            for (int j = 0; j < inputCp[i].length; ++j) {
                ujmpM.setDouble(inputCp[i][j]);
            }
        }
       
        Matrix[] eigVD = ujmpM.eigSymm();
        
        //For LAML output, the 2nd EVD matrix is the diagonal matrix of eigenvalues.
        //the first matrix is the eigen vectors, accessed by first diemsnion (rows)        
        for (int row = 0; row < eigVD[0].getRowCount(); ++row) {
            System.out.println("UJMP eigenvector i=");
            for (int col = 0; col < eigVD[0].getColumnCount(); ++col) {
                System.out.println("UJMP eigenvector i=" + eigVD[0].
            }
        }
        for (int row = 0; row < eigVD[1].getRowCount(); ++row) {
            System.out.println("LAML eigenvalue i=" + row + " " + eigVD[1].getEntry(row, row));
        } 
        */
    }
    
    public void testLinAlg() throws IOException, Exception {
        
        //String fileName = "tmp3.png";
        String fileName = "android_statues_02.jpg";
                 
        System.out.println("fileName=" + fileName);
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        String fileNameRoot = fileName.substring(0, fileName.lastIndexOf("."));
            
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        // use a faster cluster code here for inital seeds
        
        int k = 150;
        KMeansPlusPlus kmpp = new KMeansPlusPlus();
        kmpp.computeMeans(k, img.copyToGreyscale());
        int[] pixAssignements = kmpp.getImgPixelSeedIndexes();
        
        System.out.println("have initial labels (super pixels)");
        
        int w = img.getWidth();
        int h = img.getHeight();
        int n = img.getNPixels();
        
        int[][] labels = new int[w][h];
        for (int i = 0; i < w; ++i) {
            labels[i] = new int[h];
        }
        
        for (int i = 0; i < n; ++i) {
            int label = pixAssignements[i];
            labels[img.getCol(i)][img.getRow(i)] = label;
        }

        RegionAdjacencyGraphColor rag = new RegionAdjacencyGraphColor(img, labels);
        rag.populateEdgesWithColorSimilarity(ColorSpace.HSV);
        
        FlexCompColMatrix weights = rag.getEdgeMatrix();
        FlexCompColMatrix diag = createD(weights, rag, n);
        
        FlexCompColMatrix d2 = diag.copy();
        for (int i = 0; i < n; ++i) {
            double v = d2.get(i, i);
            v = 1./Math.sqrt(v);
            d2.set(i, i, v);
        }

        //D^(-1/2)(D - W)D^(-1/2)z = lamdba z;

        /*
        SimpleMatrix tmp = diag.minus(w);
        tmp = d2.mult(tmp).mult(d2);

        SimpleEVD eigVD = tmp.eig();
        
        System.out.println("n eigenvalues=" + eigVD.getNumberOfEigenvalues());
        */
        int z = 1;
    }
    
    private FlexCompColMatrix createD(FlexCompColMatrix w, RegionAdjacencyGraphColor rag,
        int nPix) {

        //D is an N X N diagonal matrix with d on the diagonal
        //    d(i) = summation over j of w(i, j) where w is "weight" of the edge
        //    and j is over all nodes
        FlexCompColMatrix d = new FlexCompColMatrix(nPix, nPix);

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        for (int i = 0; i < nPix; ++i) {
            int imgCol = rag.getColFromPixelIndex(i);
            int imgRow = rag.getRowFromPixelIndex(i);

            double dISum = 0;

            for (int k = 0; k < dxs.length; ++k) {
                int col2 = imgCol + dxs[k];
                int row2 = imgRow + dys[k];
                if (!rag.isWithinImageBounds(col2, row2)) {
                    continue;
                }
                int idx2 = rag.calculatePixelIndex(row2, col2);
                
                dISum += w.get(i, idx2);
            }

            d.set(i, i, dISum);
        }
        return d;
    }
   
}
