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
public class CornersOfLabTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public CornersOfLabTest(String testName) {
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
        
        String fileName = "lab.gif";
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleG(filePath);
        
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
        ImageIOHelper.writeOutputImage(dirPath + sep + "corners_lab.png", image);
     
        Image img2 = detector.getImage().copyImageToGreen();
        
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
        a.add(8,350);
        a.add(43,319);
        a.add(17,293);
        a.add(34,277);
        a.add(60,302);
        a.add(69,249);
        a.add(93,272);
        a.add(87,239);
        a.add(98,250);
        a.add(119,211);
        a.add(131,222);
        a.add(108,260);
        a.add(138,232);
        a.add(181,274);
        a.add(152,301);
        a.add(125,258);
        a.add(134,250);
        a.add(151,284);
        a.add(160,276);
        a.add(103,294);
        a.add(82,314);
        a.add(118,350);
        a.add(140,330);
        a.add(117,337);
        a.add(96,319);
        a.add(110,306);
        a.add(130,325);
        a.add(70,325);
        a.add(32,362);
        a.add(53,361);
        a.add(69,346);
        a.add(102,378);
        a.add(86,394);
        a.add(120,373);
        a.add(82,410);
        a.add(16,456);
        a.add(65,448);
        a.add(76,460);
        a.add(27,422);
        a.add(193,481);
        a.add(198,474);
        a.add(170,425);
        a.add(151,406);
        a.add(193,365);
        a.add(212,384);
        a.add(252,413);
        a.add(257,406);
        a.add(224,456);
        a.add(248,459);
        a.add(268,406);
        a.add(294,404);
        a.add(260,461);
        a.add(287,463);
        a.add(320,421);
        a.add(292,421);
        a.add(357,456);
        a.add(424,363);
        a.add(389,372);
        a.add(350,391);
        a.add(411,313);
        a.add(379,324);
        a.add(319,395);
        a.add(313,383);
        a.add(354,335);
        a.add(325,342);
        a.add(285,389);
        a.add(280,382);
        a.add(257,386);
        a.add(293,346);
        a.add(313,341);
        a.add(485,375);
        a.add(452,414);
        a.add(209,357);
        a.add(198,348);
        a.add(221,322);
        a.add(189,292);
        a.add(165,314);
        a.add(179,313);
        a.add(193,325);
        a.add(201,317);
        a.add(73,197);
        a.add(75,162);
        a.add(113,132);
        a.add(107,170);
        a.add(85,112);
        a.add(130,79);
        a.add(145,10);
        a.add(92,47);
        a.add(195,32);
        a.add(213,48);
        a.add(231,12);
        a.add(257,7);
        a.add(240,37);
        a.add(226,55);
        a.add(209,85);
        a.add(244,118);
        a.add(264,92);
        a.add(293,88);
        a.add(315,61);
        a.add(282,63);
        a.add(267,30);
        a.add(259,42);
        a.add(291,53);
        a.add(278,98);
        a.add(298,119);
        a.add(275,143);
        a.add(257,126);
        a.add(216,143);
        a.add(218,155);
        a.add(167,186);
        a.add(186,187);
        a.add(210,209);
        a.add(243,178);
        a.add(250,188);
        a.add(273,210);
        a.add(217,220);
        a.add(241,244);
        a.add(252,272);
        a.add(268,260);
        a.add(256,179);
        a.add(270,169);
        a.add(380,21);
        a.add(384,36);
        a.add(394,63);
        a.add(402,54);
        a.add(406,104);
        a.add(424,33);
        a.add(432,135);
        a.add(361,204);
        a.add(502,203);
        a.add(425,272);
        a.add(467,237);
        a.add(426,240);
        a.add(413,193);
        a.add(472,193);
        a.add(416,174);
        a.add(466,203);
        a.add(463,191);
        a.add(429,192);
        a.add(414,205);
        a.add(420,231);
        a.add(437,233);
        a.add(356,90);
        a.add(386,120);
        a.add(351,155);
        a.add(322,125);
        a.add(496,155);
        a.add(502,161);
        a.add(294,273);
        a.add(240,289);
        a.add(230,302);
                
        a.add(129,145);
        a.add(121,175);
        a.add(157,173);
        a.add(306,325);
        a.add(381,298);
        a.add(448,438);
        a.add(142,472);
        a.add(234,499);
        a.add(239,391);
        
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
            CornersOfLabTest test = new CornersOfLabTest("CornersOfLabTest");
            test.testProcess();
        } catch (Exception e) {
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
