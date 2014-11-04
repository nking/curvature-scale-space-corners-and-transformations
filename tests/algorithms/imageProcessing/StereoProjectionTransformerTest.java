package algorithms.imageProcessing;

import Jama.Matrix;
import algorithms.ResourceFinder;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class StereoProjectionTransformerTest {
    
    public StereoProjectionTransformerTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    @Test
    public void testC() throws Exception {
        
        String fileName1 = "books_illum3_v6_695x555.png";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        GreyscaleImage img1 = ImageIOHelper.readImageAsGrayScaleG(filePath1);

/*        String fileName2 = "books_illum3_v6_695x555.png";
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        GreyscaleImage img2 = ImageIOHelper.readImageAsGrayScaleG(filePath2);

        CurvatureScaleSpaceInflectionMapper mapper = new 
            CurvatureScaleSpaceInflectionMapper(img1, img2);
        
        mapper.useLineDrawingLineMode();

        TransformationParameters transformationParams
            = mapper.createEuclideanTransformation();
  */
        
        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img1);
               
        detector.findCorners();

        int z = 1;
    }
    
    public void testCalculateEpipolarProjection() {
        
        PairFloatArray leftXY = new PairFloatArray();
        
        PairFloatArray rightXY = new PairFloatArray();
        
        populateWithTestData0(leftXY, rightXY);
        
        StereoProjectionTransformer transformer = 
            new StereoProjectionTransformer();
        
        transformer.calculateEpipolarProjection(leftXY, rightXY);
    }
    
    private void populateWithTestData0(PairFloatArray leftXY, 
        PairFloatArray rightXY) {
        
        /*
        extracted from the books V0 and V6 images of the middlebury vision
        database.
        http://vision.middlebury.edu/stereo/data/
        
        The References listed with the data are:
        
        D. Scharstein and C. Pal. "Learning conditional random fields for 
        stereo."  In IEEE Computer Society Conference on Computer Vision and 
        Pattern Recognition (CVPR 2007), Minneapolis, MN, June 2007.

        H. Hirschmüller and D. Scharstein. "Evaluation of cost functions for 
        stereo matching."  In IEEE Computer Society Conference on Computer 
        Vision and Pattern Recognition (CVPR 2007), Minneapolis, MN, June 2007.
        
        */
        
    }
   
    public static void main(String[] args) {
        
        try {
            StereoProjectionTransformerTest test = 
                new StereoProjectionTransformerTest();
            
            test.testC();
            
        } catch(Exception e) {
            System.out.println(e.getMessage());
        }
    }
   
}
