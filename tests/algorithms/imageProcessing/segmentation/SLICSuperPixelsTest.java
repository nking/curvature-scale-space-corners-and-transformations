package algorithms.imageProcessing.segmentation;

import algorithms.imageProcessing.CannyEdgeColorAdaptive;
import algorithms.imageProcessing.CannyEdgeFilterAdaptive;
import algorithms.imageProcessing.EdgeFilterProducts;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.SIGMA;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SLICSuperPixelsTest extends TestCase {
    
    public SLICSuperPixelsTest() {
    }
    
    public void test0() throws Exception {

        String fileName1 = "";

        //for (int i = 0; i < 4; ++i) {
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
                    fileName1 = "cloudy_san_jose.jpg";
                    break;
                }
                case 7: {
                    fileName1 = "patagonia_snowy_foreground.jpg";
                    break;
                }
                case 8: {
                    fileName1 = "mt_rainier_snowy_field.jpg";
                    break;
                }
                case 9: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    break;
                }
                case 10: {
                    fileName1 = "brown_lowe_2003_image2.jpg";
                    break;
                }
                case 11: {
                    fileName1 = "venturi_mountain_j6_0001.png";
                    break;
                }
                case 12: {
                    fileName1 = "venturi_mountain_j6_0010.png";
                    break;
                }
                case 13: {
                    fileName1 = "campus_010.jpg";
                    break;
                }
                case 14: {
                    fileName1 = "campus_011.jpg";
                    break;
                }
                case 15: {
                    fileName1 = "merton_college_I_001.jpg";
                    break;
                }
                case 16: {
                    fileName1 = "merton_college_I_002.jpg";
                    break;
                }
                case 17: {
                    fileName1 = "arches.jpg";
                    break;
                }
                case 18: {
                    fileName1 = "stinson_beach.jpg";
                    break;
                }
                case 19: {
                    fileName1 = "norwegian_mtn_range.jpg";
                    break;
                }
                case 20: {
                    fileName1 = "halfdome.jpg";
                    break;
                }
                case 21: {
                    fileName1 = "halfdome2.jpg";
                    break;
                }
                case 22: {
                    fileName1 = "halfdome3.jpg";
                    break;
                }
                case 23: {
                    fileName1 = "costa_rica.jpg";
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

            int maxDimension = 256;//512;

            int binFactor1 = (int) Math.ceil(Math.max(
                (float) w1 / maxDimension,
                (float) h1 / maxDimension));
            img = imageProcessor.binImage(img, binFactor1);

            int w = img.getWidth();
            int h = img.getHeight();

            CannyEdgeColorAdaptive canny =
                new CannyEdgeColorAdaptive();
            canny.overrideToNotUseLineThinner();
            canny.applyFilter(img);
            EdgeFilterProducts edgeProducts = canny.getFilterProducts();
            
            float nPix = img.getNPixels();
            int x1 = 17;//11; 17
            float f10 = (float)w/(float)x1;
            f10 *= f10;
            float f11 = (float)h/(float)x1;
            f11 *= f11;
            int n10 = Math.round(nPix/f10);
            int n11 = Math.round(nPix/f11);
            //==> nClusters = nPix/((w/x0)^2)
            //==> nClusters = nPix/((h/x0)^2)
            int nc = (n10 + n11)/2;
            //nc = 40;
            //nc *= 2;
      nc = 100;      
            System.out.println("nc=" + nc);
            SLICSuperPixels slic = new SLICSuperPixels(img, nc);
            slic.setGradient(edgeProducts.getGradientXY());
            slic.calculate();
            int[] labels = slic.getLabels();

            Image img3 = img.createWithDimensions();
            ImageIOHelper.addAlternatingColorLabelsToRegion(img3, labels);
            String str = Integer.toString(i);
            while (str.length() < 3) {
                str = "0" + str;
            }
            MiscDebug.writeImage(img3, "_slic_" + "_" + str);
        
            NormalizedCuts norm = new NormalizedCuts();
            //norm.setColorSpaceToPolarCIELAB();
            norm.setColorSpaceToCIELAB();
            labels = norm.normalizedCut(img, labels);
            
            img3 = img.createWithDimensions();
            ImageIOHelper.addAlternatingColorLabelsToRegion(img3, labels);
            MiscDebug.writeImage(img3, "_norm_" + "_" + str);
            
        }
    }

}
