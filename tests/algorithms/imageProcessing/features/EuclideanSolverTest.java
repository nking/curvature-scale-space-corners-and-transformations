package algorithms.imageProcessing.features;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.transform.EuclideanTransformationFit;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.util.ResourceFinder;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class EuclideanSolverTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public EuclideanSolverTest() {
    }
    
    public void test0() throws Exception {
                
        String fileName1, fileName2;
        
        FeatureMatcherSettings settings = new FeatureMatcherSettings();
        settings.setDebug(true);
        settings.setStartWithBinnedImages(true);
        settings.setToUse2ndDerivCorners();
        
        for (int i = 0; i < 5; ++i) {
            switch(i) {
                case 0: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    fileName2 = "brown_lowe_2003_image2.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 1: {
                    fileName1 = "venturi_mountain_j6_0001.png";
                    fileName2 = "venturi_mountain_j6_0010.png";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 2: {
                    fileName1 = "books_illum3_v0_695x555.png";
                    fileName2 = "books_illum3_v6_695x555.png";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 3: {
                    fileName1 = "campus_010.jpg";
                    fileName2 = "campus_011.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 4: {
                    fileName1 = "merton_college_I_001.jpg";
                    fileName2 = "merton_college_I_002.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                default: {
                    fileName1 = "checkerboard_01.jpg";
                    fileName2 = "checkerboard_02.jpg";
                    settings.setUseNormalizedFeatures(false);
                    break;
                }
            }
            runEpipolarSolver(fileName1, fileName2, settings, false);
        }
    }
    
    public void testRot90() throws Exception {
                
        String fileName1, fileName2;
        
        FeatureMatcherSettings settings = new FeatureMatcherSettings();
        settings.setDebug(true);
        settings.setStartWithBinnedImages(true);
        settings.setToUse2ndDerivCorners();
        
        for (int i = 0; i < 5; ++i) {
            switch(i) {
                case 0: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    fileName2 = "brown_lowe_2003_image2.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 1: {
                    fileName1 = "venturi_mountain_j6_0001.png";
                    fileName2 = "venturi_mountain_j6_0010.png";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 2: {
                    fileName1 = "books_illum3_v0_695x555.png";
                    fileName2 = "books_illum3_v6_695x555.png";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 3: {
                    fileName1 = "campus_010.jpg";
                    fileName2 = "campus_011.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 4: {
                    fileName1 = "merton_college_I_001.jpg";
                    fileName2 = "merton_college_I_002.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                default: {
                    fileName1 = "checkerboard_01.jpg";
                    fileName2 = "checkerboard_02.jpg";
                    settings.setUseNormalizedFeatures(false);
                    break;
                }
            }
            runEpipolarSolver(fileName1, fileName2, settings, true);
        }
    }
    
    private void runEpipolarSolver(String fileName1, String fileName2, 
        FeatureMatcherSettings settings, boolean rotateBy90) throws Exception {
        
        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);
        
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
        
        if (rotateBy90) {
            TransformationParameters params90 = new TransformationParameters();
            params90.setRotationInDegrees(90);
            params90.setOriginX(0);
            params90.setOriginY(0);
            params90.setTranslationX(0);
            params90.setTranslationY(img1.getWidth() - 1);
            Transformer transformer = new Transformer();
            img1 = (ImageExt) transformer.applyTransformation(img1,
                params90, img1.getHeight(), img1.getWidth());
            int z = 1;
        }
        
        log.info("fileName1=" + fileName1);
        
        settings.setDebugTag(fileName1Root);
        
        EuclideanSolver solver = new EuclideanSolver(img1, img2, settings);
        
        EuclideanTransformationFit fit = solver.solve();
        
        assertNotNull(fit);

        TransformationParameters params = fit.getTransformationParameters();
        
        float scale = params.getScale();
        
        float rotationInDegrees = params.getRotationInDegrees();
        
        log.info("scale for " + fileName1 + " =" + scale + " rotationDeg=" + 
            rotationInDegrees);

        assertTrue(Math.abs(scale - 1) < 0.20);
        
        if (fileName1.contains("checkerboard")) {
            if (rotateBy90) {
                assertTrue(Math.abs(params.getTranslationX() - 464) < 10);
                assertTrue(Math.abs(params.getTranslationY() - -11) < 10);
                float diffRot = AngleUtil.getAngleDifference(270, params.getRotationInDegrees());
                assertTrue(Math.abs(diffRot) < 20);
            } else {
                assertTrue(Math.abs(params.getTranslationX() - -85) < 10);
                assertTrue(Math.abs(params.getTranslationY() - -11) < 10);
                float diffRot = AngleUtil.getAngleDifference(360, params.getRotationInDegrees());
                assertTrue(Math.abs(diffRot) < 20);
            }
        } else if (fileName1.contains("brown")) {
            if (rotateBy90) {
                assertTrue(Math.abs(params.getTranslationX() - 276) < 35);
                assertTrue(Math.abs(params.getTranslationY() - 5) < 35);
                float diffRot = AngleUtil.getAngleDifference(260, params.getRotationInDegrees());
                assertTrue(Math.abs(diffRot) < 20);
            } else {
                assertTrue(Math.abs(params.getTranslationX() - -250) < 35);
                assertTrue(Math.abs(params.getTranslationY() - -82) < 35);
                float diffRot = AngleUtil.getAngleDifference(350, params.getRotationInDegrees());
                assertTrue(Math.abs(diffRot) < 20);
            }
        } else if (fileName1.contains("venturi")) {
            if (rotateBy90) {
                assertTrue(Math.abs(params.getTranslationX() - 610) < 15);
                assertTrue(Math.abs(params.getTranslationY() - 0) < 15);
                float diffRot = AngleUtil.getAngleDifference(270, params.getRotationInDegrees());
                assertTrue(Math.abs(diffRot) < 20);
            } else {
                assertTrue(Math.abs(params.getTranslationX() - -40) < 25);
                assertTrue(Math.abs(params.getTranslationY() - 2) < 25);
                float diffRot = AngleUtil.getAngleDifference(345, params.getRotationInDegrees());
                assertTrue(Math.abs(diffRot) < 20);
            }
        } else if (fileName1.contains("books")) {
            if (rotateBy90) {
                assertTrue(Math.abs(params.getTranslationX() - 610) < 40);
                assertTrue(Math.abs(params.getTranslationY() - 5) < 30);
                float diffRot = AngleUtil.getAngleDifference(270, params.getRotationInDegrees());
                assertTrue(Math.abs(diffRot) < 20);
            } else {
                assertTrue(Math.abs(params.getTranslationX() - -52) < 40);
                assertTrue(Math.abs(params.getTranslationY() - -26) < 40);
                float diffRot = AngleUtil.getAngleDifference(354, params.getRotationInDegrees());
                assertTrue(Math.abs(diffRot) < 15);
            }
        } else if (fileName1.contains("campus")) {
            if (rotateBy90) {
                assertTrue(Math.abs(params.getTranslationX() - 737) < 25);
                assertTrue(Math.abs(params.getTranslationY() - 2) < 25);
                float diffRot = AngleUtil.getAngleDifference(270, params.getRotationInDegrees());
                assertTrue(Math.abs(diffRot) < 20);
            } else {
                assertTrue(Math.abs(params.getTranslationX() - 258) < 40);
                assertTrue(Math.abs(params.getTranslationY() - 7) < 40);
                float diffRot = AngleUtil.getAngleDifference(360, params.getRotationInDegrees());
                assertTrue(Math.abs(diffRot) < 20);
            }
        } else if (fileName1.contains("merton")) {
            if (rotateBy90) {
                assertTrue(Math.abs(params.getTranslationX() - 1050) < 45);
                assertTrue(Math.abs(params.getTranslationY() - -18) < 35);
                float diffRot = AngleUtil.getAngleDifference(269, params.getRotationInDegrees());
                assertTrue(Math.abs(diffRot) < 20);
            } else {
                assertTrue(Math.abs(params.getTranslationX() - -42) < 40);
                assertTrue(Math.abs(params.getTranslationY() - -7) < 40);
                float diffRot = AngleUtil.getAngleDifference(357, params.getRotationInDegrees());
                assertTrue(Math.abs(diffRot) < 15);
            }
        }
    }
}
