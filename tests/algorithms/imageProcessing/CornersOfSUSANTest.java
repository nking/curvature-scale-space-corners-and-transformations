package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import algorithms.util.PairIntArray;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class CornersOfSUSANTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public CornersOfSUSANTest(String testName) {
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
        
        //String fileName = "susan-in.gif";
        String fileName = "susan-in_plus.png";
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img);
               
        //detector.useLineDrawingMode();
               
        detector.findCorners();
        
        PairIntArray corners = detector.getCorners();
        
        PairIntArray expectedCorners = getExpectedImageCorners();
        
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
        
        log.info(foundExpectedCount + " out of " + expectedCorners.getN() 
            + " expected");
        
        log.info((corners.getN() - foundExpectedCount) 
            + " beyond expected found");
       
        Image image = ImageIOHelper.readImageAsGrayScale(filePath);
        List<PairIntArray> edges = detector.getEdgesInOriginalReferenceFrame();
        corners = detector.getCornersInOriginalReferenceFrame();
        ImageIOHelper.addAlternatingColorCurvesToImage(edges, image);
        ImageIOHelper.addCurveToImage(corners, image, 2, 2550, 0, 0);
        String dirPath = ResourceFinder.findDirectory("bin");
        String sep = System.getProperty("file.separator");
        ImageIOHelper.writeOutputImage(dirPath + sep + "corners_susan-in.png", 
            image);
    }
    
    protected PairIntArray getExpectedImageCorners() {
        
        // tolerance should be one or 2 pix in each direction because of my 
        // error in locating these coordinates manually with gimp
        
        PairIntArray a = new PairIntArray();
        a.add(0,23);
        a.add(0,46);
        a.add(0,69);
        a.add(0,92);
        a.add(0,115);
        a.add(0,138);
        a.add(0,161);
        a.add(0,184);
        a.add(0,207);
        a.add(0,230);
        a.add(0,253);
        a.add(0,276);
        a.add(0,300);
        a.add(0,323);
        a.add(0,346);
        a.add(93,24);
        a.add(93,46);
        a.add(93,69);
        a.add(93,92);
        a.add(93,115);
        a.add(93,138);
        a.add(93,161);
        a.add(93,184);
        a.add(93,207);
        a.add(93,230);
        a.add(93,253);
        a.add(93,276);
        a.add(93,300);
        a.add(93,323);
        a.add(93,346);
        a.add(143,14);
        a.add(143,157);
        a.add(214,57);
        a.add(215,23);
        a.add(289,21);
        a.add(288,14);
        a.add(301,0);
        a.add(301,94);
        a.add(407,94);
        a.add(407,0);
        a.add(324,94);
        a.add(417,0);
        a.add(417,23);
        a.add(417,45);
        a.add(417,68);
        a.add(418,104);
        a.add(511,104);
        a.add(152,100);
        a.add(287,159);
        a.add(287,142);
        a.add(173,224);
        a.add(174,166);
        a.add(112,116);
        a.add(144,196);
        a.add(202,209);
        a.add(253,209);
        a.add(252,166);
        a.add(202,166);
        a.add(227,180);
        a.add(274,208);
        a.add(274,167);
        a.add(344,208);
        a.add(343,165);
        a.add(330,208);
        a.add(288,166);
        a.add(496,275);
        a.add(392,275);
        a.add(496,367);
        a.add(392,367);
        a.add(328,345);
        a.add(328,288);
        a.add(237,288);
        a.add(237,345);
        a.add(248,302);
        a.add(248,330);
        a.add(316,330);
        a.add(316,302);
        a.add(129,345);
        a.add(158,345);
        a.add(187,346);
        a.add(187,318);
        a.add(187,288);
        a.add(158,288);
        a.add(129,288);
        a.add(129,318);

        return a;
    }
 
    public static void main(String[] args) {
        try {
            CornersOfSUSANTest test = new CornersOfSUSANTest("CornersOfSUSANTest");
            test.testProcess();
        } catch (Exception e) {
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
