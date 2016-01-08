package algorithms.imageProcessing;

import algorithms.compGeometry.RotatedOffsets;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class IntensityFeaturesTest extends TestCase {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public void test0() throws Exception {
        
        int cellDim = 2;
        int nCellsAcross = 6;
        int nColsHalf = nCellsAcross / 2;
        int range0 = cellDim * nColsHalf;
        float[] output = new float[nCellsAcross * nCellsAcross];
        float[] xT = new float[cellDim * cellDim];
        float[] yT = new float[xT.length];
        
        int count = 0;
        
        int xCenter = 10;
        int yCenter = 10;
        
        int w = 20;
        int h = 20;
        
        int rotation = 45;
        
        int tc = 0;
                
        for (int dx = -range0; dx < range0; dx += cellDim) {
            for (int dy = -range0; dy < range0; dy += cellDim) {
                
                // --- calculate values for the cell ---
                boolean withinBounds = IntensityFeatures.transformCellCoordinates(
                    rotation, xCenter, yCenter, dx, dy, cellDim, w, h, xT, yT);
                
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < xT.length; ++i) {
                    sb.append(String.format("(%.1f,%.1f) ", xT[i] - xCenter, 
                        yT[i] - yCenter));
                }
                                
                //log.fine(String.format("prev %d %d (%d, %d): %s", count, tc, 
                //    dx, dy, sb.toString()));
                
                tc += xT.length;
                
                count++;
            }
        }
    }
    
    public void testCalculateOrientation_gradient_0() throws Exception {
        
        RotatedOffsets rOffsets = RotatedOffsets.getInstance();
        
        int blockHalfWidths = 5;
        boolean useNormalizedIntensities = true;
        
        String fileName = "valve_gaussian_subimage.png";
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
        
        // values are cached so create a new instance for each "img"
        IntensityFeatures featuresGS = new IntensityFeatures(blockHalfWidths,
            useNormalizedIntensities, rOffsets);
        IntensityFeatures featuresGradient = new IntensityFeatures(blockHalfWidths,
            useNormalizedIntensities, rOffsets);
        featuresGradient.calculateGradientWithGreyscale(img);
        //MiscDebug.writeImage(featuresGradient.theta, "_theta");
        
        int x, y, expected, orientationGradient;
        
        /*
        x=52  y=53
        x=47  y=11
        x=104 y=46
        x=40  y=43
        
        98, 72
        96, 76
        32, 32
        83, 84
        */
        //makeCannyCorners(filePath);
        //make2ndDerivPoints(filePath);
        
        Logger log = Logger.getLogger(this.getClass().getName());
        
        x = 52; y = 53; expected=0;
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.fine(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
        
        x = 47; y = 11; expected=270;
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.fine(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
        
        x = 40; y = 43; expected=135;
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.fine(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
        
        x = 104; y = 46; expected=200;
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.fine(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
        
        x = 98; y = 72; expected=150;
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.fine(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
        
        x = 83; y = 84; expected=220;
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.fine(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
        
    }
    
    public void testCalculateOrientation_gradient_1() throws Exception {
        
        RotatedOffsets rOffsets = RotatedOffsets.getInstance();
        
        int blockHalfWidths = 5;
        boolean useNormalizedIntensities = true;
        
        String fileName = "checkerboard_subimage.jpg";
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
        
        // values are cached so create a new instance for each "img"
        IntensityFeatures featuresGS = new IntensityFeatures(blockHalfWidths,
            useNormalizedIntensities, rOffsets);
        IntensityFeatures featuresGradient = new IntensityFeatures(blockHalfWidths,
            useNormalizedIntensities, rOffsets);
        featuresGradient.calculateGradientWithGreyscale(img);
        //MiscDebug.writeImage(featuresGradient.gXY, "_gradient");
        
        int x, y, expected, orientationGradient;
        img = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
        
        /*
        x=26 y=30
        
        (25,29)
        (28,27)
        (24,29)
        (28,28)
        */
        //makeCannyCorners(filePath);
        make2ndDerivPoints(filePath);
        
        Logger log = Logger.getLogger(this.getClass().getName());
        
        // this corner is near intersection of 2 white and 2 black squares,
        // so might expect closer to 270 if centered, but since it is 
        // below the intersection, the dominating angle is the adjacent
        // square edge, so orientation is 180
        x = 26; y = 30; expected=180;
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.info(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
        
        x = 25; y = 29; expected=270;
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.info(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
        
        x = 28; y = 27; expected=90;
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.info(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
        
        x = 28; y = 28; expected=90;
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.info(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
        
        x = 24; y = 29; expected=270;
        orientationGradient = featuresGradient.calculateOrientation(x, y);
        log.info(String.format("(%d,%d) gr exp=%d angle=%d\n", x, y, 
            expected, orientationGradient));
        assertTrue(Math.abs(orientationGradient - expected) <= 20);
       /*
                                       90
                             *        (28,27)
                             *        (28,28) 90
         (24,29) (25,29)          *    *
          270      270   (26,30)       *  
                           180
        */
    }

    private void printDebug(GreyscaleImage img) {
        
        StringBuilder sb = new StringBuilder();
        for (int row = 0; row < img.getHeight(); ++row) {
            sb.append(String.format("row %4d: ", row));
            for (int col = 0; col < img.getWidth(); ++col) {
                sb.append(String.format("%4d ", img.getValue(col, row)));
            }
            sb.append("\n");
        }
        
        log.fine(sb.toString());
    }

    private void makeCannyCorners(String filePath) throws Exception {
        boolean useBinned = false;
        boolean filterOutImageBoundaryBlobs = true;
        boolean filterOutZeroPixels = false;
        boolean doNotAddPoints = true;
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);

        ImageExt imgCp = img.copyToImageExt();
        
        /*        
        SegmentationType type = SegmentationType.GREYSCALE_CANNY;
        BlobPerimeterCornerHelper img1Helper = new BlobPerimeterCornerHelper(img);
        img1Helper.applySegmentation(type, useBinned);

        img1Helper.getBlobs(type, useBinned, filterOutImageBoundaryBlobs, filterOutZeroPixels);

        List<List<CornerRegion>> cornersList = img1Helper.getPerimeterCorners(type, useBinned);
        
        for (List<CornerRegion> list : cornersList) {
            for (CornerRegion cr : list) {
                int x = cr.getX()[cr.getKMaxIdx()];
                int y = cr.getY()[cr.getKMaxIdx()];
                System.out.println("x=" + x + " y=" + y);
            }
            ImageIOHelper.addCornerRegionsToImage(list, imgCp, 0, 0, 1);
        }
        
        MiscDebug.writeImage(imgCp, "canny_corners");
        */
        
        CurvatureScaleSpaceCornerDetector detector = new CurvatureScaleSpaceCornerDetector(img);
        detector.findCorners();
        Set<CornerRegion> corners = detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);
        for (CornerRegion cr : corners) {
            int x = cr.getX()[cr.getKMaxIdx()];
            int y = cr.getY()[cr.getKMaxIdx()];
            System.out.println("x=" + x + " y=" + y);
        }
        ImageIOHelper.addCornerRegionsToImage(corners, imgCp, 0, 0, 1);
        
        MiscDebug.writeImage(imgCp, "canny_corners");
    }
    
    private void make2ndDerivPoints(String filePath) throws Exception {
       
        ImageExt img = ImageIOHelper.readImageExt(filePath);

        ImageExt imgCp = img.copyToImageExt();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Set<PairInt> pixels = imageProcessor.extract2ndDerivPoints(img.copyToGreyscale(),
            10, true);
        
        for (PairInt p : pixels) {
            System.out.println(p.toString());
        }
        ImageIOHelper.addCurveToImage(pixels, imgCp, 1, 255, 0, 0);
        
        MiscDebug.writeImage(imgCp, "_pts_2nd_deriv");
    }
}
