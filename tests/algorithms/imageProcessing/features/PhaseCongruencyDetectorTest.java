package algorithms.imageProcessing.features;

import algorithms.compGeometry.clustering.KMeansHSV;
import algorithms.imageProcessing.ColorHistogram;
import algorithms.imageProcessing.FixedSizeSortedVector;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.GroupPixelHSV;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.MedianTransform;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.PeriodicFFT;
import algorithms.imageProcessing.SegmentationMergeThreshold;
import algorithms.imageProcessing.SummedAreaTable;
import algorithms.imageProcessing.features.UnsupervisedTextureFinder.TexturePatchesAndResponse;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.Complex;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import com.climbwithyourfeet.clustering.DTClusterFinder;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
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
            "android_statues_01.jpg",
            "android_statues_02.jpg", "android_statues_03.jpg", "android_statues_04.jpg"
        };

        ImageProcessor imageProcessor = new ImageProcessor();

        for (String fileName : fileNames) {

            System.out.println("fileName=" + fileName);

            String filePath = ResourceFinder.findFileInTestResources(fileName);

            GreyscaleImage img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();

            PhaseCongruencyDetector phaseCDetector = new PhaseCongruencyDetector();
            phaseCDetector.setToCreateCorners();
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

            // ----- make O1 edges
            ImageExt imgClr = ImageIOHelper.readImageExt(filePath);

            GreyscaleImage o1 = imageProcessor.createO1(imgClr);

            products =
                phaseCDetector.phaseCongMono(o1);
            thinned = products.getThinned();
            out = img.createWithDimensions();
            pcImg = img.createWithDimensions();
            out2 = img.createWithDimensions();
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

            MiscDebug.writeImage(out, "_thinned_o1_" + fileName);

            MiscDebug.writeImage(out2, "_pc_thinned_o1_" + fileName);

            MiscDebug.writeImage(pcImg, "_pc_o1_" + fileName);
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

            GreyscaleImage edgeImage = imageSegmentation.createColorEdges(img);
        }
    }

    public void est2() throws Exception {

        String[] fileNames = new String[]{
            //"seattle.jpg",
            "merton_college_I_001.jpg",
            //"merton_college_I_002.jpg",
            //"lena.jpg",
            //"campus_010.jpg",
            //"android_statues_01.jpg",
            //"android_statues_02.jpg",
            //"android_statues_03.jpg",
            //"android_statues_04.jpg"
        };

        ImageSegmentation imageSegmentation = new ImageSegmentation();

        boolean doDecimate = true;
        int minDimension = 300;//512;//300;

        for (String fileName : fileNames) {
            try {
            System.out.println("fileName=" + fileName);

            String filePath = ResourceFinder.findFileInTestResources(fileName);

            ImageExt img = ImageIOHelper.readImageExt(filePath);
            List<ImageExt> transformed = null;

            int selectIdx = -1;

            if (doDecimate) {
                MedianTransform mt = new MedianTransform();
                transformed = new ArrayList<ImageExt>();
                //List<ImageExt> coeffs = new ArrayList<ImageExt>();
                //mt.multiscalePyramidalMedianTransform(img, transformed, coeffs);
                mt.<ImageExt>multiscalePyramidalMedianTransform2(img, transformed);
                //GreyscaleImage r = mt.reconstructPyramidalMultiscaleMedianTransform(
                //    transformed.get(transformed.size() - 1), coeffs);

                // choose the first image which is smaller than 300 x 300
                for (int j = 0; j < transformed.size(); ++j) {
                    ImageExt tr = transformed.get(j);
                    //MiscDebug.writeImage(tr, "_tr_" + j);
                    if (selectIdx == -1) {
                        if (tr.getWidth() <= minDimension && tr.getHeight() <= minDimension) {
                            selectIdx = j;
                            img = transformed.get(selectIdx);
                        }
                    }
                }
            }

            GreyscaleImage edgeImage = imageSegmentation.createColorEdges(img);

            edgeImage = imageSegmentation.fillInGapsOf1(edgeImage,
                new HashSet<PairInt>(), 255);

// TODO: edit performSegmentationWithColorEdges
            List<Set<PairInt>> segmentedPoints =
                imageSegmentation.performSegmentationWithColorEdges(img,
                edgeImage, SegmentationMergeThreshold.DEFAULT, fileName);

            List<PairIntArray> perimeters = BlobsAndPerimeters.extractBoundsOfBlobs(
                segmentedPoints, false, false, 1, img.getWidth(), img.getHeight());

            if (selectIdx > -1) {
                ImageProcessor imageProcessor = new ImageProcessor();
                for (int jj = (selectIdx - 1); jj > -1; --jj) {
                    perimeters = imageProcessor.unbinZeroPointLists(perimeters, 2);
                    //Image outImg = transformed.get(jj);
                    //ImageIOHelper.addAlternatingColorCurvesToImage(perimeters, outImg);
                    //MiscDebug.writeImage(outImg, "_final_edges_" + (jj + 1) + "_" + fileName);
                }
            }

            Image outImg = ImageIOHelper.readImage(filePath);
            ImageIOHelper.addAlternatingColorCurvesToImage(perimeters, outImg, 2);
            MiscDebug.writeImage(outImg, "_final_edges_" + fileName);

            } catch (Throwable t) {
                int z = 1;
            }
        }
    }

    // use of phase congruency on grey, r-g, g-b, and r-b then combining all results
    // the phase ongruency is performed on small overlapping regions.
    public void est3() throws Exception {

        String[] fileNames = new String[]{
            //"seattle.jpg",
            "merton_college_I_001.jpg",
            // "lena.jpg",
            // "campus_010.jpg",
            //"android_statues_01.jpg",
            //"android_statues_02.jpg",
            //"android_statues_03.jpg",
            //"android_statues_04.jpg"
        };

        ImageSegmentation imageSegmentation = new ImageSegmentation();

        ImageProcessor imageProcessor = new ImageProcessor();

        for (String fileName : fileNames) {

            System.out.println("fileName=" + fileName);

            String filePath = ResourceFinder.findFileInTestResources(fileName);

            ImageExt img = ImageIOHelper.readImageExt(filePath);

            GreyscaleImage combined = imageSegmentation.createColorEdges_1(img, 100);

            MiscDebug.writeImage(combined, "_combined_" + fileName);
        }
    }

    public void test4() throws Exception {

        String[] fileNames = new String[]{
            //"seattle.jpg",
            //"merton_college_I_001.jpg",
            //"house.gif",
            //"lena.jpg",
            //"campus_010.jpg",
            //"android_statues_01.jpg",
            "android_statues_02.jpg",
            //"android_statues_03.jpg",
            //"android_statues_04.jpg"
        };

        ImageSegmentation imageSegmentation = new ImageSegmentation();
        /*
        a look at O(N) patterns to make a single combined image for input to
        phase conguency that would result in closed curves for the main objects.
        */

        for (String fileName : fileNames) {

            System.out.println("fileName=" + fileName);

            String filePath = ResourceFinder.findFileInTestResources(fileName);

            ImageExt img = ImageIOHelper.readImageExt(filePath);

            //GreyscaleImage combined = imageSegmentation.createColorEdges_2(img);

            //MiscDebug.writeImage(combined, "_MAX_SOBEL_EDGES_");

            //GreyscaleImage combinedCopy = combined.copyImage();
            //ImageProcessor imageProcessor = new ImageProcessor();
            //imageProcessor.applyAdaptiveMeanThresholding(combinedCopy);
            //MiscDebug.writeImage(combinedCopy, "_MAX_SOBEL_EDGES__AMT_");

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
            -- color histograms to look for color classes in the subset noise.
            -- look into the texture stats of Malik et al 2001.
               -- whether a sure edge lies along path in between two
                  texture xlusters.
                    along the line between the pixels:
                       W_IC_i_j = 1 - argmax(local maxima of pc perpendicular
                                             to curve)
                    to search for curve points, they use a radii of 30
                    around the  textons of interest.
                  -- this might be useful during some merging steps
                     after super pixels stage.
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
                    F(x) = F(x) X log(1-(|F(x)|/0.03))/|F(x)|
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
