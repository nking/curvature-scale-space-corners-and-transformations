package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class CornersOfBloxTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public CornersOfBloxTest(String testName) {
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
        
        String fileName = "blox.gif";
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        /*boolean useSturges = false;
        ImageStatistics stats = ImageStatisticsHelper.examineImageBorders(img, 
            useSturges);
        
        ImageStatisticsHelper.plotHistogram(stats.getHistogram(), "original");
        log.info(stats.toString());*/
        
        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img);
               
        detector.findCorners();
        
        PairIntArray corners = detector.getCorners();
        
        PairIntArray expectedCorners = getExpectedLabCorners();
        
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
                int tolerance = 10;
                if ((diffX <= tolerance) && (diffY <= tolerance)) {
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
        ImageIOHelper.addCurveToImage(corners, image, 2, 255, 0, 0);
        String dirPath = ResourceFinder.findDirectory("bin");
        String sep = System.getProperty("file.separator");
        ImageIOHelper.writeOutputImage(dirPath + sep + "corners_blox.png", image);
        
        Image img2 = detector.getImage().copyToColorGreyscale();
        
        try {
            for (PairIntArray edge : detector.getEdges()) {
                debugAddCurveToImage(edge, img2, 0, 255, 255, 0);
            }
            debugAddCurveToImage(expectedCorners, img2, 2, 255, 0, 255);
            ImageIOHelper.writeOutputImage(
                dirPath + "/image_with_expected_corners.png", img2);
 
        } catch (IOException ex) {
            throw new RuntimeException("ERROR: " + ex.getMessage());
        }
    }
    
    protected PairIntArray getExpectedLabCorners() {
        
        // tolerance should be one or 2 pix in each direction because of my 
        // error in locating these coordinates manually with gimp
        
        PairIntArray a = new PairIntArray();
        
        return a;
    }
 
    private void debugAddCurveToImage(PairIntArray edge, Image input, 
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        for (int i = 0; i < edge.getN(); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                float xx = x + dx;
                if ((xx > -1) && (xx < (input.getWidth() - 1))) {
                    for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); dy++) {
                        float yy = y + dy;
                        if ((yy > -1) && (yy < (input.getHeight() - 1))) {
                            input.setRGB((int)xx, (int)yy, rClr, gClr, bClr);
                        }
                    }
                }
            }
        }
    }
    
    public static void main(String[] args) {
        try {
            CornersOfBloxTest test = new CornersOfBloxTest("CornersOfBloxTest");
            test.testProcess();
        } catch (Exception e) {
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
