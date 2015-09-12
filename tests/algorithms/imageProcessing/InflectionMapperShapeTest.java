package algorithms.imageProcessing;

import algorithms.Rotate;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import algorithms.util.PairIntArray;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class InflectionMapperShapeTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());

    public InflectionMapperShapeTest() {
    }
    
    public void testContours() throws Exception {
        
        String[] fileSuffixes = new String[]{
            //"closed_curve_circle",
            "closed_curve_triangle", 
            //"closed_curve_square"
        };

        for (String fileSuffix : fileSuffixes) {
        
            String fileName = fileSuffix + ".png";
        
            String filePath = ResourceFinder.findFileInTestResources(fileName);

            ImageExt img = ImageIOHelper.readImageExt(filePath);

            GreyscaleImage imgGrey = readAsBinary(img);

            BlobsAndContours blobsAndContours = new BlobsAndContours(imgGrey, 
                0, 1 << 30, fileSuffix);

            boolean setToExtractWeakCurvesTooIfNeeded = false;

            int edgeIndex = 0;

            List<PairIntArray> curves = blobsAndContours.getBlobOrderedPerimeters();

            assertNotNull(curves);

            assertTrue(curves.size() == 1);
            
            PairIntArray edge = curves.get(0);
            /*
            if (fileSuffix.contains("triangle")) {
                int n = edge.getN();
                Rotate rotate = new Rotate();
                rotate.rotate2(edge.getX(), (n-1)/2);
                rotate.rotate2(edge.getY(), (n-1)/2);
            }*/

            ScaleSpaceCurveImage sscImg =
                CurvatureScaleSpaceInflectionSingleEdgeMapper.createScaleSpaceImage(
                edge, edgeIndex);

            List<CurvatureScaleSpaceContour> c =
                CurvatureScaleSpaceInflectionSingleEdgeMapper.populateContours(
                sscImg, edgeIndex, setToExtractWeakCurvesTooIfNeeded);

            MiscDebug.printScaleSpaceCurve(sscImg, fileSuffix);

            MiscDebug.debugPlot(c, img, 0, 0, fileSuffix);
        }
        
    }
    
    public static void main(String[] args) {
        
        try {
            
            InflectionMapperShapeTest test = new InflectionMapperShapeTest();
                        
            test.testContours();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }

    private GreyscaleImage readAsBinary(ImageExt img) {
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        GreyscaleImage out = new GreyscaleImage(w, h);
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int r = img.getR(i, j);
                if (r > 125) {
                    out.setValue(i, j, 255);
                }
            }
        }
        
        return out;
    }
}
