package algorithms.imageProcessing;

import algorithms.imageProcessing.features.mser.MSEREdges;
import algorithms.imageProcessing.segmentation.LabelToColorHelper;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.QuadInt;
import algorithms.util.ResourceFinder;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.set.TIntSet;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class VanishingPointsTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public VanishingPointsTest() {
    }

    public void test0() throws Exception {

        /*
        NOTE that the vanishing lines for the Merton college
        images could be improved w/ improved segmentation.
        Also, for the Seattle test images, can see that
        finding rectangles first and then small lines associated
        with them and then using the frequency of lines not
        orthogonal to the 2-d image plane would find the
        projected dimension.
        */
        
        int maxDimension = 256;//512;

        String fileName1 = "";

        //for (int i = 0; i < 5; ++i) {
        for (int i = 0; i < 38; ++i) {
System.out.println("index i=" + i);
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
                case 36: {
                    fileName1 = "checkerboard_02.jpg";
                    break;
                }
                default: {
                    fileName1 = "house_color.png";
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

            MSEREdges mserEdges = new MSEREdges(img);
            mserEdges.setToLowerContrast();
            mserEdges.mergeAndExtractEdges();
            List<TIntSet> pointSets = mserEdges.getLabeledSets();
            List<Set<PairInt>> contigSets = new ArrayList<Set<PairInt>>(pointSets.size());
            for (TIntSet set : pointSets) {
                Set<PairInt> set2 = imageProcessor.convertIndexesToPoints(set, 
                    img.getWidth());
                contigSets.add(set2);
            }
            
            EdgeFilterProducts products = mserEdges.getEdgeFilterProducts();       
            
            ImageExt img11 = img.createWithDimensions();

            VanishingPoints vp2 = new VanishingPoints();
            vp2.setToDebug();
            //vp2.dbgImg = img.copyImage();
            vp2.find(contigSets, img.getWidth(), img.getHeight());
            //vp2.correctLinesWithGradient(products.getGradientXY());            
            //MiscDebug.writeImage(img, "_lines_" + fileName1Root);
        
            TIntObjectMap<QuadInt> linesMap = vp2.getVanishingLines();
            TIntObjectIterator<QuadInt> iter = linesMap.iterator();
            for (int ii = 0; ii < linesMap.size(); ++ii) {
                iter.advance();
                QuadInt eps = iter.value();
                int segIdx = iter.key();
                int clr = ImageIOHelper.getNextColorRGB(ii);
                
                ImageIOHelper.drawLineInImage(
                    eps.getA(), eps.getB(), eps.getC(), eps.getD(), img, 1, clr);
            }
            MiscDebug.writeImage(img, fileName1 + "_lines_");
        }
    }

     public static void main(String[] args) {

        try {
            VanishingPointsTest test = new VanishingPointsTest();
            //test.test0();
            //test.testRot90();

        } catch(Exception e) {
            e.printStackTrace();
            System.err.println(e.getMessage());
            fail(e.getMessage());
        }
    }

}
