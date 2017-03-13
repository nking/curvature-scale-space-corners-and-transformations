package algorithms.imageProcessing.scaleSpace;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;

/**
 * @author nichole
 */
public class CSSContourInflectionMakerTest extends TestCase {
    
    public CSSContourInflectionMakerTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
     
    public void testScaleSpaceImagesFigure2() throws Exception {
        
        System.out.println("testScaleSpaceImagesFigure2");
                                
        // IEEE 'TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE, 
        // VOL. PAMI-8, NO. 1. JANUARY 1986
        // Scale-Based Description and Recognition of Planar Curves and 
        // Two-Dimensional Shapes by FARZIN MOKHTARIAN AND ALAN MACKWORTH
        String filePath = ResourceFinder.findFileInTestResources("africa2.png");
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        int imageWidth = img.getWidth();
        int imageHeight = img.getHeight();
        
        CSSContourInflectionMaker maker = new
            CSSContourInflectionMaker(img);
                        
        maker.setToUseLineDrawingLineMode();
        maker.setToDebug();
        
        maker.findContours();
        
    }
    
}
