package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PointMatcher1Test extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public PointMatcher1Test() {
    }
    
    public void testPatchesAreSimilar() throws Exception {
        
        //TODO: add to these tests
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        String fileName1 = "brown_lowe_2003_image1.jpg";
        String fileName2 = "brown_lowe_2003_image2.jpg";
        
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
        HistogramEqualizationForColor hEq = new HistogramEqualizationForColor(img1);
        hEq.applyFilter();
        hEq = new HistogramEqualizationForColor(img2);
        hEq.applyFilter();
        img1 = imageProcessor.binImage(img1, 3);
        img2 = imageProcessor.binImage(img2, 3);
        
        Transformer transformer = new Transformer();
        
        PointMatcher pointMatcher = new PointMatcher();
        
        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();

        // transformation from point pairs:
        /*   img1          img2
         (127, 88)      (31, 82)
         (162, 55)      (69, 56)
        */
        TransformationParameters params = tc.calulateEuclidean(
            127, 88, 162, 55, 31, 82, 69, 56, 0, 0);
        
        pointMatcher.debug = true;
        
        log.info("SIMILAR:");
        boolean similar = pointMatcher.patchesAreSimilar(params, 
            img1, img2, transformer, 162, 55);
        //assertTrue(similar);
        
        log.info("NOT SIMILAR:");
        // something should like different
        TransformationParameters paramsWrong = params.copy();
        paramsWrong.setTranslationY(-30);
        similar = pointMatcher.patchesAreSimilar(paramsWrong, 
            img1, img2, transformer, 162, 55);
        //assertFalse(similar);
        
        log.info("SIMILAR:");
        similar = pointMatcher.patchesAreSimilar(params, 
            img1, img2, transformer, 127, 88);
        //assertTrue(similar);
        
        log.info("NOT SIMILAR:");
        similar = pointMatcher.patchesAreSimilar(paramsWrong, 
            img1, img2, transformer, 127, 88);
        //assertFalse(similar);


        log.info("SIMILAR:");
        similar = pointMatcher.patchesAreSimilar(params, 
            img1, img2, transformer, 157, 82);
        //assertTrue(similar);
        
        log.info("SIMILAR:");
        similar = pointMatcher.patchesAreSimilar(params, 
            img1, img2, transformer, 134, 89);
        //assertTrue(similar);
        
        log.info("SIMILAR:");
        similar = pointMatcher.patchesAreSimilar(params, 
            img1, img2, transformer, 146, 67);
        //assertTrue(similar);
        
        log.info("SIMILAR:");
        similar = pointMatcher.patchesAreSimilar(params, 
            img1, img2, transformer, 140, 60);
        //assertTrue(similar);
        
        log.info("SIMILAR:");
        similar = pointMatcher.patchesAreSimilar(params, 
            img1, img2, transformer, 117, 82);
        //assertTrue(similar);
        
        log.info("SIMILAR:");
        similar = pointMatcher.patchesAreSimilar(params, 
            img1, img2, transformer, 119, 111);
        //assertTrue(similar);
    }

    public void testSortByDescendingMatches3() throws Exception {
        
        PointMatcher matcher = new PointMatcher();
        
        TransformationPointFit[] fits = new TransformationPointFit[4];
        
        fits[0] = new TransformationPointFit(
            new TransformationParameters(),
            10, 10.0, 5.0, 
            10.0f, 10.0f);
        
        fits[1] = new TransformationPointFit(
            new TransformationParameters(),
            20, 5.0, 2.0, 
            10.0f, 10.0f);
        
        fits[2] = new TransformationPointFit(
            new TransformationParameters(),
            2, 25.0, 20.0, 
            10.0f, 10.0f);
        
        fits[3] = new TransformationPointFit(
            new TransformationParameters(),
            20, 7.0, 6.0, 
            10.0f, 10.0f);
      
        matcher.sortByDescendingMatches(fits, 0, fits.length - 1);
        
        assertTrue(fits[0].getNumberOfMatchedPoints() == 20);
        assertTrue(Math.abs(fits[0].getMeanDistFromModel() - 5) < 0.1);
        
        assertTrue(fits[1].getNumberOfMatchedPoints() == 20);
        assertTrue(Math.abs(fits[1].getMeanDistFromModel() - 7) < 0.1);
        
        assertTrue(fits[2].getNumberOfMatchedPoints() == 10);
        
        assertTrue(fits[3].getNumberOfMatchedPoints() == 2);
    }
    
    public static void main(String[] args) {

        try {
            PointMatcher1Test test = new PointMatcher1Test();

            test.testSortByDescendingMatches3();

        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            fail(e.getMessage());
        }
    }
}
