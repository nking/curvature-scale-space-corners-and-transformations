package algorithms.imageProcessing;

import algorithms.imageProcessing.features.BlobPerimeterCornerHelper;
import algorithms.imageProcessing.features.CornerRegion;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class HoughTransformTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public HoughTransformTest(String testName) {
        super(testName);
    }

    public void testLines1() throws Exception {

        String fileName1, fileName2;

        for (int i = 5; i < 6; ++i) {
            //fileName1 = "valve_gaussian.png";
            //fileName2 = "valve_gaussian.png";
            switch(i) {
                case 0: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    fileName2 = "brown_lowe_2003_image2.jpg";
                    break;
                }
                case 1: {
                    fileName1 = "venturi_mountain_j6_0001.png";
                    fileName2 = "venturi_mountain_j6_0010.png";
                    break;
                }
                case 2: {
                    fileName1 = "books_illum3_v0_695x555.png";
                    fileName2 = "books_illum3_v6_695x555.png";
                    break;
                }
                case 3: {
                    fileName1 = "campus_010.jpg";
                    fileName2 = "campus_011.jpg";
                    break;
                }
                case 4: {
                    fileName1 = "merton_college_I_001.jpg";
                    fileName2 = "merton_college_I_002.jpg";
                    break;
                }
                default: {
                    fileName1 = "checkerboard_01.jpg";
                    fileName2 = "checkerboard_02.jpg";
                    break;
                }
            }

            System.out.println("fileName1=" + fileName1);

            String bin = ResourceFinder.findDirectory("bin");
            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
            int idx = fileName1.lastIndexOf(".");
            String fileNameRoot = fileName1.substring(0, idx);

            ImageProcessor imageProcessor = new ImageProcessor();

            boolean useBinned = true;
            int binnedImageMaxDimension = 512;
            int binFactor1 = 1;
            int binFactor2 = 1;
            SegmentationType type = SegmentationType.GREYSCALE_HIST;

            GreyscaleImage img1 = ImageIOHelper.readImage(filePath1).copyToGreyscale();
            GreyscaleImage img2 = ImageIOHelper.readImage(filePath2).copyToGreyscale();

            if (useBinned) {
                binFactor1 = (int) Math.ceil(Math.max((float) img1.getWidth() / binnedImageMaxDimension,
                (float) img1.getHeight() / binnedImageMaxDimension));

                binFactor2 = (int) Math.ceil(Math.max((float) img2.getWidth() / binnedImageMaxDimension,
                (float) img2.getHeight() / binnedImageMaxDimension));

                img1 = imageProcessor.binImage(img1, binFactor1);
                img2 = imageProcessor.binImage(img2, binFactor2);
            }

            BlobPerimeterCornerHelper bph = new BlobPerimeterCornerHelper(
                ImageIOHelper.readImageExt(filePath1), fileNameRoot);
            bph.createBinnedGreyscaleImage(binnedImageMaxDimension);
            bph.applySegmentation(type, useBinned);
            
            List<List<CornerRegion>> cornerRegionLists =
                bph.generatePerimeterCorners(type, useBinned);

            GreyscaleImage segImg1 = useBinned ?
                bph.getBinnedSegmentationImage(type) :
                bph.getSegmentationImage(type);
            
            Image tmp1SegImg1 = segImg1.copyToColorGreyscale();
            
            List<PairIntArray> edgeLists = bph.getBlobPerimeters(type, useBinned);

            for (int ii = 0; ii < edgeLists.size(); ++ii) {

                //NOTE: in testable method for this, should allow ability to
                // pass in junctions and not delete corners that are in
                // junctions.
                // For these blob perimeters, there are not junctions.
                
                PairIntArray edge = edgeLists.get(ii);

                List<CornerRegion> cornerRegions = cornerRegionLists.get(ii);

                if (cornerRegions.size() < 2) {
                    continue;
                }
                
                ImageIOHelper.addCurveToImage(edge, tmp1SegImg1, 1, 255, 255, 255);

                for (int j = 0; j < cornerRegions.size(); ++j) {
                    CornerRegion cr = cornerRegions.get(j);
                    int x = cr.getX()[cr.getKMaxIdx()];
                    int y = cr.getY()[cr.getKMaxIdx()];
                    System.out.println("plotting: (" + x + "," + y + ")");
                    ImageIOHelper.addPointToImage(x, y, tmp1SegImg1, 0, 0, 2,
                        255, 0, 0);
                }
            }

            ImageIOHelper.writeOutputImage(
                bin + "/seg_1_hough1_" + fileNameRoot + ".png", tmp1SegImg1);

            int z0 = 1;
        }
    }
    
    private double distance(CornerRegion cr, PairInt p) {
        int x1 = cr.getX()[cr.getKMaxIdx()];
        int y1 = cr.getY()[cr.getKMaxIdx()];

        int diffX = x1 - p.getX();
        int diffY = y1 - p.getY();
        
        return Math.sqrt(diffX * diffX + diffY * diffY);
    }
    
    private double distance(PairInt p1, PairInt p2) {

        int diffX = p1.getX() - p2.getX();
        int diffY = p1.getY() - p2.getY();
        
        return Math.sqrt(diffX * diffX + diffY * diffY);
    }

    // assuming cr was tested as further from endPoints than 2 pixels
    private boolean isInBetween(PairInt[] endPoints, double distBetweenEndPoints,
        CornerRegion cr) {
        
        double d0 = distance(cr, endPoints[0]);
        
        double d1 = distance(cr, endPoints[1]);
        
        return (d0 + d1) < distBetweenEndPoints;
    }

    public void testFindContiguousLines()  throws Exception {
        
        /*
        public List<Set<PairInt>> findContiguousLines(Set<PairInt> points, 
        GreyscaleImage theta360) {
        */
        /*
         *           Y
         *          90
         *     135   |    +45
         *           |
         *   180------------ 0   X
         *           |
         *    225    |   315
         *          270
         *
        Most of the theta image values are w.r.t. lines at origin of above diagram, for
        example.  A line at 225 degrees would look like this in the image:
                 /
               /  
             /
        
        The results of atan2 on gY and gX resulted in values from -pi to pi radians
        then those were transformed to 0 to 360 degrees
       
        for blox.gif,
        
        math.atan2(gY, gX)
        
        
        There is a vertical line with expected theta being 90 or 270,
           but theta image values are *orthogonal* to that:  
           (198, 163) math.atan2(-1,-19)*180./math.pi=183
           (198, 157) math.atan2(0,-17)*180./math.pi=180
           might be a deceptive feature... consider whether the 2-1D binomial 
           smoothing on a vertical feature produces a highly local
           appearance of angle 180 preferentially along an axis as an artifact
           while globally, the whole vertical line has that behavior...
           comparing other points now...
        
        A line that is /  expected 215
            has theta image value of 239
            so is *not orthogonal* to expected  
            (180, 84) math.atan2(-10,-6)*180./math.pi=239
        
        A line that is /  expected 45
            has theta image value of about 45
            so is *not orthogonal* to expected 
            (60, 72) math.atan2(18,16)*180./math.pi=48
        
        a line that is expected to have angle 135
            has theta image values of about 128
            so is *not orthogonal*
            (46, 118) math.atan2(19,-15)*180./math.pi=128
        
        a line that is expected to have angle 315
            has theta image values of about 305
            so is *not orthogonal*
            (83, 235) math.atan2(-32,22)*180./math.pi=305
        
        a horizontal line -- expected 180 or 0
            has theta image values of about 90
            so is *orthogonal* to expected
            (38, 228) math.atan2(46,0)*180./math.pi=90
            The gX and gY values 
        
        For what angles subtended from horizontal or vertical can a significant
        gX or gY be measured.
            The first blur is one sigma, which has a FWHM of about 2.35 
            and a FWZI ~ 4.5 pix.
            The gradient uses sobel which is a sigma of sqrt(1)/2 so total sigma is sqrt(1.25),
            and then FWHM is 2.35*sqrt(1.25) which is about 2.6.
            So the minimum measurable angle would be a slope of 1 pixel vs 3 pixels,
            math.atan2(1,3)*180./math.pi = 18.4 degrees.
        
        So, for any angle which has gX or gY < 3 pixels, the theta angle
        in the image is present at values orthogonal to expected at
        angles that are horizontal or vertical (and their actual lines in
        the image are vertical or horizontal, respectively).
        
        */
        /*
        7
        6       @ @ @ @
        5     @         @
        4   @             @
        3 @                 @
        2   @             @
        1     @         @
        0       @ @ @ @
          0 1 2 3 4 5 6 7 8 9 
        */
    }
}
