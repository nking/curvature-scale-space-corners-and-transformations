package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class EpipolarSolverTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public EpipolarSolverTest() {
    }
    
    public void test0() throws Exception {
                
        String fileName1, fileName2;
        
        // TODO: follow up on changes needed for repeated patterns
        //    with small projection effects.
        //    need to use matching of top k solutions...
        
        for (int i = 0; i < 5; ++i) {
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
            runEpipolarSolver(fileName1, fileName2, false, false);
        }
    }
    
    public void testRot90() throws Exception {
                
        String fileName1, fileName2;
        
        for (int i = 0; i < 5; ++i) {
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
            runEpipolarSolver(fileName1, fileName2, true, false);
        }
    }
    
    public void testParameters() throws Exception {
                
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
            runEpipolarSolver(fileName1, fileName2, false, true);
        }
    }

    private void runEpipolarSolver(String fileName1, String fileName2, 
        boolean rotateBy90, boolean useParameters) throws Exception {
        
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
        
        EpipolarSolver solver = null;
        if (!useParameters) {
            solver = new EpipolarSolver(img1, img2, fileName1Root);
        } else {
            TransformationParameters parameters = new TransformationParameters();
            parameters.setOriginX(0);
            parameters.setOriginY(0);
            
            if (fileName1Root.contains("brown")) {
                parameters.setRotationInDegrees(352);
                parameters.setScale(0.98f);
                parameters.setTranslationX(-241);
                parameters.setTranslationY(-64);
                float[] stdevs = new float[]{0.035f, 0.036f, 16.1f, 18.3f};
                parameters.setStandardDeviations(stdevs);
            } else if (fileName1Root.contains("venturi")) {
                parameters.setRotationInDegrees(1);
                parameters.setScale(1);
                parameters.setTranslationX(-50);
                parameters.setTranslationY(-1);
                float[] stdevs = new float[]{0.01f, 0.01f, 4.6f, 4.8f};
                parameters.setStandardDeviations(stdevs);
            } else if (fileName1Root.contains("books")) {
                parameters.setRotationInDegrees(1);
                parameters.setScale(1);
                parameters.setTranslationX(-64f);
                parameters.setTranslationY(0.f);
                float[] stdevs = new float[]{0, 0, 0, 0};
                parameters.setStandardDeviations(stdevs);
            } else if (fileName1Root.contains("campus")) {
                parameters.setRotationInDegrees(0.8f);
                parameters.setScale(1);
                parameters.setTranslationX(258f);
                parameters.setTranslationY(3);
                float[] stdevs = new float[]{0.003f, 0.005f, 1.4f, 0.93f};
                parameters.setStandardDeviations(stdevs);
            } else if (fileName1Root.contains("merton")) {
                parameters.setRotationInDegrees(359.f);
                parameters.setScale(1.06f);
                parameters.setTranslationX(21);
                parameters.setTranslationY(-2);
                float[] stdevs = new float[]{0.01f, 0.014f, 8.95f, 11.43f};
                parameters.setStandardDeviations(stdevs);
            }
            solver = new EpipolarSolver(img1, img2, parameters, fileName1Root);
        }
        
        StereoProjectionTransformerFit fit = solver.solve();
        
        assertNotNull(fit);
        
    }
}
