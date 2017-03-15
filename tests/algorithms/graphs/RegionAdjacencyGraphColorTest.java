package algorithms.graphs;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.segmentation.ColorSpace;
import algorithms.imageProcessing.segmentation.SLICSuperPixels;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;

/**
 *
 * @author nichole
 */
public class RegionAdjacencyGraphColorTest extends TestCase {
    
    public RegionAdjacencyGraphColorTest() {
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
        rag.populateEdgesWithColorSimilarity(ColorSpace.HSV, 22);
        
        System.out.println("created region agency graph and similarity matrix");
        System.out.flush();
        FlexCompRowMatrix weights = rag.getEdgeMatrix();

    }
   
}
