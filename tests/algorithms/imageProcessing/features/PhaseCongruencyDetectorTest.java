package algorithms.imageProcessing.features;

import algorithms.imageProcessing.CannyEdgeFilterLite;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.MedianTransform;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.SegmentationMergeThreshold;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
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

        float cutOff = 0.5f;//0.3f;//0.5f;
        int nScale = 5;
        int minWavelength = 3;//nScale;// 3;
        float mult = 2.1f;
        float sigmaOnf = 0.55f;
        int k = 5;//2;
        float g = 10; 
        float deviationGain = 1.5f;
        int noiseMethod = -1;
        double tLow = 0.0001;
        double tHigh = 0.1;
        boolean increaseKIfNeeded = false;

        ImageProcessor imageProcessor = new ImageProcessor();
        
        String label = "label=n" + nScale + "_mw" + minWavelength + "_k" + k 
            + "_t0_" + tLow + "_t1_" + tHigh;
            
        for (String fileName : fileNames) {
            
            System.out.println("fileName=" + fileName);
        
            String filePath = ResourceFinder.findFileInTestResources(fileName);
        
            GreyscaleImage img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();
        
            PhaseCongruencyDetector phaseCDetector = new PhaseCongruencyDetector();
            phaseCDetector.setToCreateCorners();                
            PhaseCongruencyDetector.PhaseCongruencyProducts products =
                phaseCDetector.phaseCongMono(img, nScale, minWavelength, mult, 
                sigmaOnf, k, increaseKIfNeeded,
                cutOff, g, deviationGain, noiseMethod, tLow, tHigh);

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
            MiscDebug.writeImage(out, "_thinned_" + label + "_" + fileName + "_"); 
            MiscDebug.writeImage(out2, "_pc_thinned_" + label + "_" + fileName + "_");
            MiscDebug.writeImage(pcImg, "_pc_" + label + "_" + fileName + "_");
            
            // ----- make O1 edges
            ImageExt imgClr = ImageIOHelper.readImageExt(filePath);
            
            GreyscaleImage o1 = imageProcessor.createO1(imgClr);
            
            products =
                phaseCDetector.phaseCongMono(o1, nScale, minWavelength, mult, 
                sigmaOnf, k, increaseKIfNeeded,
                cutOff, g, deviationGain, noiseMethod, tLow, tHigh);
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
        
            MiscDebug.writeImage(out, "_thinned_o1_" + label + "_" + fileName);
            
            MiscDebug.writeImage(out2, "_pc_thinned_o1_" + label + "_" + fileName);
            
            MiscDebug.writeImage(pcImg, "_pc_o1_" + label + "_" + fileName);
        }
    }
    
    public void est1() throws Exception {

        String[] fileNames = new String[]{
            //"blox.gif", "lab.gif", "house.gif", "seattle.jpg", "merton_college_I_001.jpg",
            // "susan-in_plus.png", "lena.jpg",
            // "campus_010.jpg", 
            //"android_statues_01.jpg", 
            //"android_statues_02.jpg", 
            //"android_statues_03.jpg", 
            "android_statues_04.jpg"
        };
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
                     
        for (String fileName : fileNames) {
            
            System.out.println("fileName=" + fileName);
        
            String filePath = ResourceFinder.findFileInTestResources(fileName);
        
            ImageExt img = ImageIOHelper.readImageExt(filePath);            
            
            GreyscaleImage edgeImage = imageSegmentation.createColorEdges(img);
        }
    }
    
    public void test2() throws Exception {

        String[] fileNames = new String[]{
            //"blox.gif", "lab.gif", "house.gif", "seattle.jpg", 
            // "merton_college_I_001.jpg", "susan-in_plus.png", "lena.jpg",
            // "campus_010.jpg", 
            //"android_statues_01.jpg", 
            //"android_statues_02.jpg", 
            //"android_statues_03.jpg", 
            "android_statues_04.jpg"
        };
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        boolean doDecimate = true;
        int minDimension = 512;//300;
        
        for (String fileName : fileNames) {
            
            System.out.println("fileName=" + fileName);
        
            String filePath = ResourceFinder.findFileInTestResources(fileName);
        
            ImageExt img = ImageIOHelper.readImageExt(filePath);     
            
            if (doDecimate) {
                MedianTransform mt = new MedianTransform();
                List<ImageExt> transformed = new ArrayList<ImageExt>();
                //List<ImageExt> coeffs = new ArrayList<ImageExt>();
                //mt.multiscalePyramidalMedianTransform(img, transformed, coeffs);    
                mt.<ImageExt>multiscalePyramidalMedianTransform2(img, transformed);
                //GreyscaleImage r = mt.reconstructPyramidalMultiscaleMedianTransform(
                //    transformed.get(transformed.size() - 1), coeffs);
                
                // choose the first image which is smaller than 300 x 300
                int selectIdx = -1;
                for (int j = 0; j < transformed.size(); ++j) {
                    ImageExt tr = transformed.get(j);
                    MiscDebug.writeImage(tr, "_tr_" + j);
                    if (selectIdx == -1) {
                        if (tr.getWidth() <= minDimension && tr.getHeight() <= minDimension) {
                            selectIdx = j;
                        }
                    }
                }
                
                img = transformed.get(selectIdx);
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
            
            ImageIOHelper.addAlternatingColorCurvesToImage(perimeters, img);
            
            MiscDebug.writeImage(img, "_final_edges_");
        }
    }
    
    // use of phase congruency on grey, r-g, g-b, and r-b then combining all results
    // the phase ongruency is performed on small overlapping regions.
    public void est3() throws Exception {

        String[] fileNames = new String[]{
            //"blox.gif", "lab.gif", "house.gif", "seattle.jpg", "merton_college_I_001.jpg",
            // "susan-in_plus.png", "lena.jpg",
            // "campus_010.jpg", 
            //"android_statues_01.jpg", 
            //"android_statues_02.jpg", 
            //"android_statues_03.jpg", 
            "android_statues_04.jpg"
        };
        
        /*
        float cutOff = 0.5f;//0.3f;//0.5f;
        int nScale = 5;
        int minWavelength = 3;//nScale;// 3;
        float mult = 2.1f;
        float sigmaOnf = 0.55f;
        int k = 10;//2;
        float g = 10;
        float deviationGain = 1.5f;
        int noiseMethod = -1;
        double tLow = 0.0001;
        double tHigh = 0.1;
        boolean increaseKIfNeeded = true;
        */
        float cutOff = 0.5f;//0.3f;//0.5f;
        int nScale = 5;//4
        int minWavelength = 3;//nScale;// 3;
        float mult = 2.1f;
        float sigmaOnf = 0.55f;
        int k = 10;//2;
        float g = 10;
        float deviationGain = 1.5f;
        int noiseMethod = -1;
        double tLow = 0.001;
        double tHigh = 0.1;
        boolean increaseKIfNeeded = false;
                     
        for (String fileName : fileNames) {
            
            System.out.println("fileName=" + fileName);
        
            String filePath = ResourceFinder.findFileInTestResources(fileName);
        
            ImageExt img = ImageIOHelper.readImageExt(filePath);            
            
            /*
            0 grey
            1 r-g
            2 b-g
            3 r-b
            */
            int x0 = 600;
            int x1 = 1200;
            int y0 = 150;
            int y1 = 750;
            int width = x1 - x0;
            int height = y1 - y0;
            final int sz = 100;
            
            GreyscaleImage combined = new GreyscaleImage(width, height);
            for (int xOff = 0; xOff < width; xOff += sz) {
                if (xOff > 0) {
                    xOff -= 10;
                }
                for (int yOff = 0; yOff < height; yOff += sz) {
                    if (yOff > 0) {
                        yOff -= 10;
                    }
                    for (int clrIdx = 0; clrIdx < 4; ++clrIdx) {

                        float[] values = null;
                        String lbl = "";
                        int count = 0;
                        if (clrIdx == 0) {
                            lbl = "_grey_";
                            values = read_R_G_B(img, x0 + xOff, x0 + xOff + sz, y0 + yOff, y0 + yOff + sz);
                        } else if (clrIdx == 1) {
                            lbl = "_r-g_";
                            values = read_R_minus_G(img, x0 + xOff, x0 + xOff + sz, y0 + yOff, y0 + yOff + sz);
                        } else if (clrIdx == 2) {
                            lbl = "_b-g_";
                            values = read_B_minus_G(img, x0 + xOff, x0 + xOff + sz, y0 + yOff, y0 + yOff + sz);
                        } else if (clrIdx == 3) {
                            lbl = "_r-b_";
                            values = read_R_minus_B(img, x0 + xOff, x0 + xOff + sz, y0 + yOff, y0 + yOff + sz);
                        }
                        values = MiscMath.rescale(values, 0, 255);

                        GreyscaleImage img2 = new GreyscaleImage(sz, sz);
                        count = 0;
                        int e0 = sz;
                        if ((x0 + xOff + sz) > x1) {
                            e0 = x1 - x0 - xOff;//sz - (x0 + xOff + sz - x1);
                        }
                        int e1 = sz;
                        if ((y0 + yOff + sz) > y1) {
                            e1 = sz - (y0 + yOff + sz - y1);
                        }
                        for (int i = 0; i < e0; ++i) {
                            for (int j = 0; j < e1; ++j) {
                                int v = Math.round(values[count]);
                                img2.setValue(i, j, v);
                                count++;
                            }
                        }

                        PhaseCongruencyDetector phaseDetector = new PhaseCongruencyDetector();
                        PhaseCongruencyDetector.PhaseCongruencyProducts products
                            = phaseDetector.phaseCongMono(img2, nScale, minWavelength, mult,
                                sigmaOnf, k, increaseKIfNeeded,
                                cutOff, g, deviationGain, noiseMethod, tLow, tHigh);
                        int[][] thinned = products.getThinned();
                        assert(thinned.length == sz);
                        assert(thinned[0].length == sz);
                        for (int j = 0; j < thinned.length; ++j) {
                            for (int i = 0; i < thinned[j].length; ++i) {
                                if (((i + xOff) > (combined.getWidth() - 1)) ||
                                    ((j + yOff) > (combined.getHeight() - 1))) {
                                    continue;
                                }
                                if (thinned[j][i] > 0) {
                                    combined.setValue(i + xOff, j + yOff, 255);
                                }
                            }
                        }
                    }
                }
            }
            MiscDebug.writeImage(combined, "_COMBINED_EDGES_");
            
        }
    }
       
    public void est4() throws Exception {

        String[] fileNames = new String[]{
            //"blox.gif", "lab.gif", "house.gif", "seattle.jpg", "merton_college_I_001.jpg",
            // "susan-in_plus.png", "lena.jpg",
            // "campus_010.jpg", 
            //"android_statues_01.jpg", 
            //"android_statues_02.jpg", 
            //"android_statues_03.jpg", 
            "android_statues_04.jpg"
        };
        
        /*
        a look at O(N) patterns to make a single combined image for input to
        phase conguency that would result in closed curves for the main objects.
        */
                     
        for (String fileName : fileNames) {
            
            System.out.println("fileName=" + fileName);
        
            String filePath = ResourceFinder.findFileInTestResources(fileName);
        
            ImageExt img = ImageIOHelper.readImageExt(filePath);            
            
            GreyscaleImage[] gradients = new GreyscaleImage[4];
            GreyscaleImage[] masks = new GreyscaleImage[4];
            String[] labels = new String[4];
            /*
            0 grey    min:    0    max: 255
            1 r-g     min: -255    max: 255
            2 b-g       "
            3 r-b       "
            */
            for (int clrIdx = 0; clrIdx < 4; ++clrIdx) {
                masks[clrIdx] = new GreyscaleImage(img.getWidth(), 
                    img.getHeight());
                gradients[clrIdx] = new GreyscaleImage(img.getWidth(), 
                    img.getHeight(), GreyscaleImage.Type.Bits32FullRangeInt);
                if (clrIdx == 0) {
                    labels[clrIdx] = "_grey_";
                    for (int i = 0; i < img.getNPixels(); ++i) {
                        int v = (img.getR(i) + img.getG(i) + img.getB(i))/3;
                        gradients[clrIdx].setValue(i, v);
                    }
                } else if (clrIdx == 1) {
                    labels[clrIdx] = "_r-g_";
                    for (int i = 0; i < img.getNPixels(); ++i) {
                        int v = img.getR(i) - img.getG(i);
                        v = (v + 255)/2;
                        gradients[clrIdx].setValue(i, v);
                    }
                } else if (clrIdx == 2) {
                    labels[clrIdx] = "_b-g_";
                    for (int i = 0; i < img.getNPixels(); ++i) {
                        int v = img.getB(i) - img.getG(i);
                        v = (v + 255)/2;
                        gradients[clrIdx].setValue(i, v);
                    }
                } else if (clrIdx == 3) {
                    labels[clrIdx] = "_r-b_";
                    for (int i = 0; i < img.getNPixels(); ++i) {
                        int v = img.getR(i) - img.getB(i);
                        v = (v + 255)/2;
                        gradients[clrIdx].setValue(i, v);
                    }
                }                 
                CannyEdgeFilterLite filter = new CannyEdgeFilterLite();
                filter.setToUseSobel();
                filter.applyFilter(gradients[clrIdx]);
            }
            
            GreyscaleImage combined =  new GreyscaleImage(img.getWidth(), 
                img.getHeight(), GreyscaleImage.Type.Bits32FullRangeInt);
            for (int i = 0; i < img.getNPixels(); ++i) {
                int gradientMax = Integer.MIN_VALUE;
                int maxClrIdx = -1;
                for (int clrIdx = 0; clrIdx < 4; ++clrIdx) {
                    int v = gradients[clrIdx].getValue(i);
                    if (v > gradientMax) {
                        gradientMax = v;
                        maxClrIdx = clrIdx;
                    }
                }
                combined.setValue(i, gradientMax);
                masks[maxClrIdx].setValue(i, 255);
            }
            
            MiscDebug.writeImage(combined, "_MAX_SOBEL_EDGES_");
            
           // for (int clrIdx = 0; clrIdx < 4; ++clrIdx) {
           //     MiscDebug.writeImage(masks[clrIdx], "_mask_" + labels[clrIdx]);
           // }
            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.blur(combined, SIGMA.ONE);
            
            float[] values = new float[img.getNPixels()];
            GreyscaleImage gsImg = img.copyToGreyscale();
            for (int i = 0; i < img.getNPixels(); ++i) {
                values[i] = gsImg.getValue(i) + 10*combined.getValue(i);
            }
            values = MiscMath.rescale(values, 0, 255);
            
            for (int i = 0; i < img.getNPixels(); ++i) {
                gsImg.setValue(i, Math.round(values[i]));
            }
            MiscDebug.writeImage(gsImg, "_greyscale_plus_edges_");
            
            float cutOff = 0.5f;//0.3f;//0.5f;
            int nScale = 5;//4
            int minWavelength = 3;//nScale;// 3;
            float mult = 2.1f;
            float sigmaOnf = 0.55f;
            int k = 10;//2;
            float g = 10;
            float deviationGain = 1.5f;
            int noiseMethod = -1;
            double tLow = 0.00001;
            double tHigh = 0.1;
            boolean increaseKIfNeeded = false;
            
            PhaseCongruencyDetector phaseDetector = new PhaseCongruencyDetector();
            PhaseCongruencyDetector.PhaseCongruencyProducts products
                = phaseDetector.phaseCongMono(gsImg, nScale, minWavelength, mult,
                    sigmaOnf, k, increaseKIfNeeded,
                    cutOff, g, deviationGain, noiseMethod, tLow, tHigh);
            int[][] thinned = products.getThinned();
            for (int j = 0; j < thinned.length; ++j) {
                for (int i = 0; i < thinned[j].length; ++i) {
                    int v = 0;
                    if (thinned[j][i] > 0) {
                        v = 255;
                    }
                    gsImg.setValue(i, j, v);
                }
            }
            MiscDebug.writeImage(gsImg, "_greyscale_plus_edges_PC_");
        }
    }
    
    private float[] read_R_G_B(ImageExt img, int x0, int x1, int y0, int y1) {
        
        int nx = x1 - x0;
        int ny = y1 - y0;
        
        float[] values = new float[nx * ny];
        
        int count = 0;
        
        for (int x = x0; x < x1; ++x) {
            if (x > (img.getWidth() - 1)) {
                break;
            }
            for (int y = y0; y < y1; ++y) {
                if (y > (img.getHeight() - 1)) {
                    break;
                }
                int r = img.getR(x, y);
                int g = img.getG(x, y);
                int b = img.getB(x, y);
                values[count] = r + g + b;
                count++;
            }
        }
        if (count < nx * ny) {
            values = Arrays.copyOf(values, count);
        }
        return values;
    }
    
    private float[] read_R_minus_B(ImageExt img, int x0, int x1, int y0, int y1) {
        
        int nx = x1 - x0;
        int ny = y1 - y0;
        
        float[] values = new float[nx * ny];
        
        int count = 0;
        
        for (int x = x0; x < x1; ++x) {
            if (x > (img.getWidth() - 1)) {
                break;
            }
            for (int y = y0; y < y1; ++y) {
                if (y > (img.getHeight() - 1)) {
                    break;
                }
                int r = img.getR(x, y);
                int b = img.getB(x, y);
                values[count] = r - b;
                count++;
            }
        }
        if (count < nx * ny) {
            values = Arrays.copyOf(values, count);
        }
        return values;
    }
    
    private float[] read_B_minus_G(ImageExt img, int x0, int x1, int y0, int y1) {
        
        int nx = x1 - x0;
        int ny = y1 - y0;
        
        float[] values = new float[nx * ny];
        
        int count = 0;
        
        for (int x = x0; x < x1; ++x) {
            if (x > (img.getWidth() - 1)) {
                break;
            }
            for (int y = y0; y < y1; ++y) {
                if (y > (img.getHeight() - 1)) {
                    break;
                }
                int g = img.getG(x, y);
                int b = img.getB(x, y);
                values[count] = b - g;
                count++;
            }
        }
        if (count < nx * ny) {
            values = Arrays.copyOf(values, count);
        }
        return values;
    }
    
    private float[] read_R_minus_G(ImageExt img, int x0, int x1, int y0, int y1) {
        
        int nx = x1 - x0;
        int ny = y1 - y0;
        
        float[] values = new float[nx * ny];
        
        int count = 0;
        
        for (int x = x0; x < x1; ++x) {
            if (x > (img.getWidth() - 1)) {
                break;
            }
            for (int y = y0; y < y1; ++y) {
                if (y > (img.getHeight() - 1)) {
                    break;
                }
                int r = img.getR(x, y);
                int g = img.getG(x, y);
                values[count] = r - g;
                count++;
            }
        }
        if (count < nx * ny) {
            values = Arrays.copyOf(values, count);
        }
        return values;
    }
}
