package algorithms.graphs;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.segmentation.ColorSpace;
import algorithms.imageProcessing.segmentation.SLICSuperPixels;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.Iterator;
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
        
        SLICSuperPixels slic = new SLICSuperPixels(img, 200);
        slic.calculate();
        int[] labels = slic.getLabels();
        
        System.out.println("have initial labels (super pixels)");
        System.out.flush();
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        RegionAdjacencyGraphColor rag = new RegionAdjacencyGraphColor(img, labels);
        rag.populateEdgesWithColorSimilarity(ColorSpace.HSV);
        
        System.out.println("created region agency graph and similarity matrix");
        System.out.flush();
        FlexCompRowMatrix weights = rag.getEdgeMatrix();

    }
   
}
