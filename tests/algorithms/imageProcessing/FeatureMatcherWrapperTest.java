package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.util.Collection;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class FeatureMatcherWrapperTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public FeatureMatcherWrapperTest() {
    }
    
    public void test0() throws Exception {
                
        String fileName1, fileName2;
        
        // TODO: follow up on changes needed for repeated patterns
        //    with small projection effects.
        //    need to use matching of top k solutions...
        
        for (int i = 0; i < 4; ++i) {
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
                default: {
                    fileName1 = "merton_college_I_001.jpg";
                    fileName2 = "merton_college_I_002.jpg";
                    break;
                }
            }
            runCorrespondenceList(fileName1, fileName2, false, false);
        }
    }
    
    public void estRot90() throws Exception {
                
        String fileName1, fileName2;
        
        for (int i = 0; i < 4; ++i) {
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
                default: {
                    fileName1 = "campus_010.jpg";
                    fileName2 = "campus_011.jpg";
                    break;
                }
            }
            runCorrespondenceList(fileName1, fileName2, true, false);
        }
    }
    
    public void estParameters() throws Exception {
                
        String fileName1, fileName2;
        
        for (int i = 0; i < 4; ++i) {
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
                default: {
                    fileName1 = "campus_010.jpg";
                    fileName2 = "campus_011.jpg";
                    break;
                }
            }
            runCorrespondenceList(fileName1, fileName2, false, true);
        }
    }

    private void runCorrespondenceList(String fileName1, String fileName2, 
        boolean rotateBy90, boolean useScale) throws Exception {
        
        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);
        idx = fileName2.lastIndexOf(".");
        String fileName2Root = fileName2.substring(0, idx);
        
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
        
        FeatureMatcherWrapper wrapper = null;
        if (!useScale) {
            wrapper = new FeatureMatcherWrapper(img1, img2, fileName1Root);
        } else {
            TransformationParameters parameters = new TransformationParameters();
            parameters.setOriginX(0);
            parameters.setOriginY(0);
            
            if (fileName1Root.contains("brown")) {
                parameters.setRotationInDegrees(350);
                parameters.setScale(0.97f);
                parameters.setTranslationX(-222.85f);
                parameters.setTranslationY(-72.91f);
                float[] stdevs = new float[]{0.019f, 0.028f, 6.55f, 13.1f};
                parameters.setStandardDeviations(stdevs);
            } else if (fileName1Root.contains("venturi")) {
                parameters.setRotationInDegrees(0);
                parameters.setScale(1);
                parameters.setTranslationX(-30.18f);
                parameters.setTranslationY(2.28f);
                float[] stdevs = new float[]{0f, 0f, 0.22f, 0.08f};
                parameters.setStandardDeviations(stdevs);
            } else if (fileName1Root.contains("books")) {
                parameters.setRotationInDegrees(0);
                parameters.setScale(1);
                parameters.setTranslationX(-64f);
                parameters.setTranslationY(0.f);
                float[] stdevs = new float[]{0, 0, 0, 0};
                parameters.setStandardDeviations(stdevs);
            } else if (fileName1Root.contains("campus")) {
                parameters.setRotationInDegrees(0);
                parameters.setScale(1);
                parameters.setTranslationX(265.27f);
                parameters.setTranslationY(9.55f);
                float[] stdevs = new float[]{0.02f, 0.013f, 4.67f, 5.29f};
                parameters.setStandardDeviations(stdevs);
            }
            wrapper = new FeatureMatcherWrapper(img1, img2, parameters, fileName1Root);
        }
        
        CorrespondenceList cl = wrapper.matchFeatures();
        
        assertNotNull(cl);
        
        Collection<PairInt> m1 = cl.getPoints1();
        Collection<PairInt> m2 = cl.getPoints2();
        GreyscaleImage gsImg1 = img1.copyToGreyscale();
        GreyscaleImage gsImg2 = img2.copyToGreyscale();
        String name1 = "1_" + fileName1Root;
        String name2 = "2_" + fileName2Root;
        if (rotateBy90) {
            name1 = name1 + "_r90";
            name2 = name2 + "_r90";
        }
        name1 = name1 + "_matched";
        name2 = name2 + "_matched";
        MiscDebug.plotCorners(gsImg1, m1, name1, 2);
        MiscDebug.plotCorners(gsImg2, m2, name2, 2);
        
        float scale = cl.getScale();
        int rotationInDegrees = cl.getRotationInDegrees();
        
        log.info("scale for " + fileName1 + " =" + scale + " rotationDeg=" + 
            rotationInDegrees);

        if (fileName1.contains("brown_lowe_2003_image1")) {
            // one portion of image scale is ~ 0.9
            assertTrue(Math.abs(scale - 1) < 0.15);
        } else {
            assertTrue(Math.abs(scale - 1) < 0.12);
        }

    }
}
