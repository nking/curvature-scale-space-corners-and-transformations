package algorithms.imageProcessing.features.mser;

import algorithms.imageProcessing.CannyEdgeFilterAdaptiveDeltaE2000;
import algorithms.imageProcessing.EdgeFilterProducts;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class MSEREdgesTest extends TestCase {
    
    public MSEREdgesTest() {
    }

    public void test0() throws Exception {

        int maxDimension = 256;//512;

        String fileName1 = "";

        //for (int i = 6; i < 7; ++i) {
        for (int i = 0; i < 37; ++i) {

            switch(i) {
                case 0: {
                    fileName1 = "android_statues_01.jpg";
                    break;
                }
                case 1: {
                    fileName1 = "android_statues_02.jpg";
                    break;
                }
                case 2: {
                    fileName1 = "android_statues_03.jpg";
                    break;
                }
                case 3: {
                    fileName1 = "android_statues_04.jpg";
                    break;
                }
                case 4: {
                    fileName1 = "seattle.jpg";
                    break;
                }
                case 5: {
                    fileName1 = "stonehenge.jpg";
                    break;
                }
                case 6: {
                    fileName1 = "costa_rica.jpg";
                    break;
                }
                case 7: {
                    fileName1 = "cloudy_san_jose.jpg";
                    break;
                }
                case 8: {
                    fileName1 = "patagonia_snowy_foreground.jpg";
                    break;
                }
                case 9: {
                    fileName1 = "mt_rainier_snowy_field.jpg";
                    break;
                }
                case 10: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    break;
                }
                case 11: {
                    fileName1 = "brown_lowe_2003_image2.jpg";
                    break;
                }
                case 12: {
                    fileName1 = "venturi_mountain_j6_0001.png";
                    break;
                }
                case 13: {
                    fileName1 = "venturi_mountain_j6_0010.png";
                    break;
                }
                case 14: {
                    fileName1 = "campus_010.jpg";
                    break;
                }
                case 15: {
                    fileName1 = "campus_011.jpg";
                    break;
                }
                case 16: {
                    fileName1 = "merton_college_I_001.jpg";
                    break;
                }
                case 17: {
                    fileName1 = "merton_college_I_002.jpg";
                    break;
                }
                case 18: {
                    fileName1 = "arches.jpg";
                    break;
                }
                case 19: {
                    fileName1 = "stinson_beach.jpg";
                    break;
                }
                case 20: {
                    fileName1 = "norwegian_mtn_range.jpg";
                    break;
                }
                case 21: {
                    fileName1 = "halfdome.jpg";
                    break;
                }
                case 22: {
                    fileName1 = "halfdome2.jpg";
                    break;
                }
                case 23: {
                    fileName1 = "halfdome3.jpg";
                    break;
                }
                case 24: {
                    fileName1 = "new-mexico-sunrise_w725_h490.jpg";
                    break;
                }
                case 25: {
                    fileName1 = "arizona-sunrise-1342919937GHz.jpg";
                    break;
                }
                case 26: {
                    fileName1 = "sky_with_rainbow.jpg";
                    break;
                }
                case 27: {
                    fileName1 = "sky_with_rainbow2.jpg";
                    break;
                }
                case 28: {
                    fileName1 = "books_illum3_v0_695x555.png";
                    break;
                }
                case 29: {
                    fileName1 = "books_illum3_v6_695x555.png";
                    break;
                }
                case 30: {
                    fileName1 = "klein_matterhorn_snowy_foreground.jpg";
                    break;
                }
                case 31: {
                    fileName1 = "30.jpg";
                    break;
                }
                case 32: {
                    fileName1 = "arches_sun_01.jpg";
                    break;
                }
                case 33: {
                    fileName1 = "stlouis_arch.jpg";
                    break;
                }
                case 34: {
                    fileName1 = "contrail.jpg";
                    break;
                }
                case 35: {
                    fileName1 = "checkerboard_01.jpg";
                    break;
                }
                default: {
                    fileName1 = "checkerboard_02.jpg";
                    break;
                }
            }

            int idx = fileName1.lastIndexOf(".");
            String fileName1Root = fileName1.substring(0, idx);

            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            ImageExt img = ImageIOHelper.readImageExt(filePath1);

            ImageProcessor imageProcessor = new ImageProcessor();
            ImageSegmentation imageSegmentation = new ImageSegmentation();

            int w1 = img.getWidth();
            int h1 = img.getHeight();

            int binFactor1 = (int) Math.ceil(Math.max(
                (float) w1 / maxDimension,
                (float) h1 / maxDimension));

            img = imageProcessor.binImage(img, binFactor1);
            //MiscDebug.writeImage(img, "_"  + fileName1Root);
                
            MSEREdges mserE = new MSEREdges(img);
            mserE.setToDebug();
            //mserE.extractEdges();
            mserE.setToLowerContrast();
            mserE.mergeAndExtractEdges();
        }
    }
    
}
