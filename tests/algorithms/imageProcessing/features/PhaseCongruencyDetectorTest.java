package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.features.UnsupervisedTextureFinder.TexturePatchesAndResponse;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PhaseCongruencyDetectorTest extends TestCase {

    public PhaseCongruencyDetectorTest() {
    }

    public void est0() throws Exception {

        String[] fileNames = new String[]{
           // "blox.gif", "lab.gif", "house.gif", "seattle.jpg", "merton_college_I_001.jpg",
           // "susan-in_plus.png", "lena.jpg",
           // "campus_010.jpg",
            "android_statues_01.jpg", "android_statues_02.jpg", 
            "android_statues_03.jpg", "android_statues_04.jpg"
        };

        ImageProcessor imageProcessor = new ImageProcessor();

        for (String fileName : fileNames) {

            System.out.println("fileName=" + fileName);

            String filePath = ResourceFinder.findFileInTestResources(fileName);

            GreyscaleImage img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();

            PhaseCongruencyDetector phaseCDetector = new PhaseCongruencyDetector();
            PhaseCongruencyDetector.PhaseCongruencyProducts products =
                phaseCDetector.phaseCongMono(img);

            assertNotNull(products);
            int[][] thinned = products.getThinned();

            GreyscaleImage pcImg = img.createWithDimensions();
            GreyscaleImage out2 = img.createWithDimensions();
            GreyscaleImage out = img.createWithDimensions();
            for (int i = 0; i < out.getWidth(); ++i) {
                for (int j = 0; j < out.getHeight(); ++j) {
                    int vPC = (int)Math.round(255. * products.getPhaseCongruency()[j][i]);
                    if (thinned[j][i] > 0) {
                        out.setValue(i, j, thinned[j][i]);
                        out2.setValue(i, j, vPC);
                    }
                    pcImg.setValue(i, j, vPC);
                }
            }
            MiscDebug.writeImage(out, "_thinned_" + fileName + "_");
            MiscDebug.writeImage(out2, "_pc_thinned_" + fileName + "_");
            MiscDebug.writeImage(pcImg, "_pc_" + fileName + "_");    
        }
    }

    public void est1() throws Exception {

        String[] fileNames = new String[]{
            //"blox.gif",
            //"lab.gif",
            "house.gif",
            //"seattle.jpg",
            //"merton_college_I_001.jpg",
            // "susan-in_plus.png",
            //"lena.jpg",
            //"campus_010.jpg",
            //"android_statues_01.jpg",
            //"android_statues_02.jpg",
            //"android_statues_03.jpg",
            //"android_statues_04.jpg"
        };

        ImageSegmentation imageSegmentation = new ImageSegmentation();

        for (String fileName : fileNames) {

            System.out.println("fileName=" + fileName);

            String filePath = ResourceFinder.findFileInTestResources(fileName);

            ImageExt img = ImageIOHelper.readImageExt(filePath);

        }
    }

    public void testUnsupervisedTextureExtraction() throws Exception {

        String[] fileNames = new String[]{
            //"seattle.jpg",
            //"merton_college_I_001.jpg",
            //"house.gif",
            //"lena.jpg",
            //"campus_010.jpg",           //*
            //"android_statues_01.jpg",
            "android_statues_02.jpg",
            //"android_statues_03.jpg",    //*
            //"android_statues_04.jpg"
        };

        ImageProcessor imageProcessor = new ImageProcessor();

        ImageSegmentation imageSegmentation = new ImageSegmentation();
        /*
        a look at O(N) patterns to make a single combined image for input to
        phase conguency that would result in closed curves for the main objects.
        */

        int maxDimension = 256;//512;

        for (String fileName : fileNames) {

            System.out.println("fileName=" + fileName);

            String filePath = ResourceFinder.findFileInTestResources(fileName);

            ImageExt img = ImageIOHelper.readImageExt(filePath);

            int w1 = img.getWidth();
            int h1 = img.getHeight();
            int binFactor1 = (int) Math.ceil(Math.max(
                (float) w1 / maxDimension,
                (float) h1 / maxDimension));
            img = imageProcessor.binImage(img, binFactor1);


            GreyscaleImage img2 = img.copyToGreyscale2();

            PhaseCongruencyDetector phaseCDetector
                = new PhaseCongruencyDetector();
            phaseCDetector.setToExtractNoise();
            phaseCDetector.setToDebug();
            //phaseCDetector.setToCreateCorners();
            PhaseCongruencyDetector.PhaseCongruencyProducts products =
                phaseCDetector.phaseCongMono(img2);

            UnsupervisedTextureFinder finder = new
                UnsupervisedTextureFinder();

            TexturePatchesAndResponse[] tpar =
                finder.createTextureImages(img, products, fileName);

            /*
            -- look into the texture stats of Malik et al 2001.
               -- whether a sure edge lies along path in between two
                  texture xlusters.
                     along the line between the pixels:
                        W_IC_i_j = 1 - argmax(local maxima of pc perpendicular
                                   to curve)
                     to search for curve points, they use a radii of 30
                     around the  textons of interest.

               -- different goal: adding other terms to the normalized cuts
                  weighting function is in this paper and the DNCuts paper.

                  The term they use for chi squared as the difference between
                  textons and intensities, could be used to make
                  a weighting function for the difference between color
                  histograms:
                   chi squared =
                      1/2 times sum over all bins of : (h1 - h2)^2/(h1 + h2)
                   W_TX_i_j = exp(- chi squared / sigma_TX)
                    (note, for the intensity weight component, sigma is 0.02,
                     and sigma_TX = 0.025,
                     so might expect similar for other weight sigmas.
                -- They also include suggestions for using position, that is
                   adjacency.

            -- explore making a frequency domain filter for spatial domain
               from a representative set of subset noise or centered on it.
               -- might be able to reduce computations due to sparse
                  data
               -- Malik et al. 2001 normalize their texton responses using:
                    F(x) = F(x) X log(1+(|F(x)|/0.03))/|F(x)|
            */


            /*
            assertNotNull(products);
            int[][] thinned = products.getThinned();

            GreyscaleImage pcImg = img2.createWithDimensions();
            GreyscaleImage out2 = img2.createWithDimensions();
            GreyscaleImage out = img2.createWithDimensions();
            for (int i = 0; i < out.getWidth(); ++i) {
                for (int j = 0; j < out.getHeight(); ++j) {
                    int vPC = (int)Math.round(255. * products.getPhaseCongruency()[j][i]);
                    if (thinned[j][i] > 0) {
                        out.setValue(i, j, thinned[j][i]);
                        out2.setValue(i, j, vPC);
                    }
                    pcImg.setValue(i, j, vPC);
                }
            }
            MiscDebug.writeImage(out, "_thinned_" + "_" + fileName + "_");
            MiscDebug.writeImage(out2, "_pc_thinned_"  + "_" + fileName + "_");
            MiscDebug.writeImage(pcImg, "_pc_" + "_" + fileName + "_");


            PhaseCongruencyDetectorPyramidal phaseCDetector0
                = new PhaseCongruencyDetectorPyramidal();
            phaseCDetector0.setToCreateCorners();
            PhaseCongruencyDetectorPyramidal.PhaseCongruencyProducts products0 =
                phaseCDetector0.phaseCongMono(img2);

            assertNotNull(products0);
            thinned = products0.getThinned();
            pcImg = img2.createWithDimensions();
            out2 = img2.createWithDimensions();
            out = img2.createWithDimensions();
            for (int i = 0; i < out.getWidth(); ++i) {
                for (int j = 0; j < out.getHeight(); ++j) {
                    int vPC = (int)Math.round(255.
                        * products0.getPhaseCongruency()[j][i]);
                    if (thinned[j][i] > 0) {
                        out.setValue(i, j, thinned[j][i]);
                        out2.setValue(i, j, vPC);
                    }
                    pcImg.setValue(i, j, vPC);
                }
            }
            MiscDebug.writeImage(out, "_thinned_" + "_" + fileName + "_0");
            MiscDebug.writeImage(out2, "_pc_thinned_"  + "_" + fileName + "_0");
            MiscDebug.writeImage(pcImg, "_pc_" + "_" + fileName + "_0");
            */
        }
    }
}
