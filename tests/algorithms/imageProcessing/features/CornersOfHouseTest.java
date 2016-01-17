package algorithms.imageProcessing.features;

import algorithms.imageProcessing.*;
import algorithms.imageProcessing.scaleSpace.CurvatureScaleSpaceCornerDetector;
import algorithms.util.ResourceFinder;
import algorithms.util.PairIntArray;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class CornersOfHouseTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public CornersOfHouseTest(String testName) {
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
   
    public void testProcess() throws Exception {
        
        /*
        objectives:
            (1) all true corners should be detected
            (2) no false corners should be detected.
            (3) corner points should be well localized
            (4) corner detectors should be robust with respect to noise
            (5) corner detectors should be efficient        
        */
        
        String fileName = "house.gif";
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img);
       
        detector.findCorners();
        
        PairIntArray corners = detector.getCorners();
        
        PairIntArray expectedCorners = getExpectedHouseCorners();
        
        int foundExpectedCount = 0;
        for (int i = 0; i < expectedCorners.getN(); i++) {
            int x = expectedCorners.getX(i);
            int y = expectedCorners.getY(i);
            for (int j = 0; j < corners.getN(); j++) {
                int xx = corners.getX(j);
                int yy = corners.getY(j);
                int diffX = xx - x;
                int diffY = yy - y;
                if (diffX < 0) {
                    diffX *= -1;
                }
                if (diffY < 0) {
                    diffY *= -1;
                }
                if ((diffX <= 4) && (diffY <= 4)) {
                    foundExpectedCount++;
                    break;
                }
            }
        }
        
        Image image = ImageIOHelper.readImageAsGrayScale(filePath);
        List<PairIntArray> edges = detector.getEdgesInOriginalReferenceFrame();
        corners = detector.getCornersInOriginalReferenceFrame();
        ImageIOHelper.addAlternatingColorCurvesToImage(edges, image);
        ImageIOHelper.addCurveToImage(corners, image, 2, 255, 0, 0);
        String dirPath = ResourceFinder.findDirectory("bin");
        String sep = System.getProperty("file.separator");
        ImageIOHelper.writeOutputImage(dirPath + sep + "corners_house.png", image);
        
        log.info(foundExpectedCount + " out of " + expectedCorners.getN() 
            + " expected");
        
        log.info((corners.getN() - foundExpectedCount) 
            + " beyond expected found");
       
    }
    
    protected PairIntArray getExpectedHouseCorners() {
        
        // tolerance should be one or 2 pix in each direction because of my 
        // error in locating these coordinates manually with gimp
        
        PairIntArray a = new PairIntArray();
        a.add(66, 34);
        a.add(94, 23);
        a.add(114, 29);
        a.add(65, 86);
        a.add(66, 95);
        a.add(63, 109);
        a.add(29, 135);
        a.add(17, 138);
        a.add(40, 135);
        a.add(38, 134);
        a.add(37, 156);
        a.add(57, 142);
        a.add(54, 175);       
        a.add(33, 213);
        a.add(30, 247);
        a.add(22, 230);
        a.add(150, 222);
        a.add(119, 209);
        a.add(239, 206);
        a.add(253, 192);
        a.add(239, 198);
        a.add(100, 52);
        a.add(113, 50);
        a.add(110, 72);
        a.add(174, 98);
        a.add(174, 110);
        a.add(154, 120);
        a.add(172, 114);
        a.add(154, 120);
        a.add(114, 129);
        a.add(111, 166);
        a.add(135, 162);
        a.add(118, 158);
        a.add(129, 129);
        a.add(152, 199);
        a.add(172, 192);
        a.add(196, 188);
        a.add(194, 121);
        a.add(205, 107);
        a.add(190, 101);
        a.add(190, 107);
        
        return a;
    }
 
    public static void main(String[] args) {
        try {
            CornersOfHouseTest test = new CornersOfHouseTest("CornersOfHouseTest");
            test.testProcess();
            int z = 1;
        } catch (Exception e) {
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
