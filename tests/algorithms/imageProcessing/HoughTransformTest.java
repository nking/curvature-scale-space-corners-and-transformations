package algorithms.imageProcessing;

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

            BlobPerimeterHelper bph = new BlobPerimeterHelper(
                ImageIOHelper.readImageExt(filePath1), fileNameRoot);
            bph.createBinnedGreyscaleImage(binnedImageMaxDimension);
            bph.applySegmentation(type, useBinned);
            BlobCornerHelper bch = new BlobCornerHelper(bph, fileNameRoot);
            
            List<List<CornerRegion>> cornerRegionLists =
                bch.generatePerimeterCorners(type, useBinned);

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

}
