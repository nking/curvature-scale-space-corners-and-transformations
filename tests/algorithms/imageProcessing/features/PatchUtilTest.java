package algorithms.imageProcessing.features;

import algorithms.Rotate;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.MiscDebug;
import algorithms.util.PixelHelper;
import algorithms.util.ResourceFinder;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.io.IOException;
import java.util.Arrays;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertTrue;

/**
 * 
 * @author nichole
 */
public class PatchUtilTest extends TestCase {
    
    public PatchUtilTest() {
    }

    public void test1() throws IOException {
        
        int w = 60;
        int h = 60;
        int xc1 = 226;
        int yc1 = 160;
        
        int dw = w-1;
        int dh = h-1;
        
        //int N_PIX_PER_CELL_DIM = 4;
        //int N_CELLS_PER_BLOCK_DIM = 2;
        int N_PIX_PER_CELL_DIM = 3;    //1;  3;
        int N_CELLS_PER_BLOCK_DIM = 3; //10; 3;
        //int angle = 20;
        int nBins = 8;//180/angle;
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        // make a test with spiral.png
        String filePath = ResourceFinder.findFileInTestResources(
            "android_statues_objects.png");
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        GreyscaleImage gsImg = img.copyToGreyscale2();
        
        GreyscaleImage img1 = gsImg.subImage(xc1, yc1, w, h);    

        //MiscDebug.writeImage(img1, "img1");
        //MiscDebug.writeImage(img1_rot, "img1_rot");
        
        HOGs hogs = new HOGs(img1, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, 
            nBins);
        
        PixelHelper ph = new PixelHelper();
        
        /*
               (....)
             (       )
            (  1   2  )
            -----------
            |          |
            | 3    4   |
            |          |
            ------------
        */
        int w2 = img1.getWidth();
        int h2 = img1.getHeight();
        TIntSet s1 = new TIntHashSet();
        TIntSet s2 = new TIntHashSet();
        TIntSet s3 = new TIntHashSet();
        TIntSet s4 = new TIntHashSet();
        
        for (int x = 20; x < 30; ++x) {
            for (int y = 17; y < 25; ++y) {
                long pixIdx = ph.toPixelIndex(x, y, w2);
                s1.add((int)pixIdx);
            }
        }
        for (int x = 30; x < 40; ++x) {
            for (int y = 17; y < 25; ++y) {
                long pixIdx = ph.toPixelIndex(x, y, w2);
                s2.add((int)pixIdx);
            }
        }
        for (int x = 20; x < 30; ++x) {
            for (int y = 30; y < 40; ++y) {
                long pixIdx = ph.toPixelIndex(x, y, w2);
                s3.add((int)pixIdx);
            }
        }
        for (int x = 30; x < 40; ++x) {
            for (int y = 30; y < 40; ++y) {
                long pixIdx = ph.toPixelIndex(x, y, w2);
                s4.add((int)pixIdx);
            }
        }
        
        PatchUtil p1 = new PatchUtil(w2, h2, nBins);
        PatchUtil p2 = new PatchUtil(w2, h2, nBins);
        PatchUtil p3 = new PatchUtil(w2, h2, nBins);
        PatchUtil p4 = new PatchUtil(w2, h2, nBins);
        
        p1.add(s1, hogs);
        p2.add(s2, hogs);
        p3.add(s3, hogs);
        p4.add(s4, hogs);
        
        System.out.println("HOGs:");
        System.out.println("p1=" + Arrays.toString(p1.getHistogram()));
        System.out.println("p2=" + Arrays.toString(p2.getHistogram()));
        System.out.println("p3=" + Arrays.toString(p3.getHistogram()));
        System.out.println("p4=" + Arrays.toString(p4.getHistogram()));
        
        String[] labels = new String[]{"p1", "p2", "p3", "p4"};
        PatchUtil[] ps = new PatchUtil[]{p1, p2, p3, p4};
        
        for (int i = 0; i < labels.length; ++i) {
            for (int j = (i + 1); j < labels.length; ++j) {
                double inter = ps[i].intersection(ps[j]);
                System.out.format("%s intersection %s = %.3f\n", 
                    labels[i], labels[j], (float)inter);
            }
        }
        for (int i = 0; i < labels.length; ++i) {
            for (int j = (i + 1); j < labels.length; ++j) {
                double[] diffAndErr = ps[i].diff(ps[j]);
                System.out.format("%s difference %s = %.3f(+-%.3f)\n", 
                    labels[i], labels[j], (float)diffAndErr[0], (float)diffAndErr[1]);
            }
        }
        
        // ==== HCPT =====
        
        GreyscaleImage luvTheta = 
            imageProcessor.createCIELUVTheta_WideRangeLightness(img, 255);
        img1 = luvTheta.subImage(xc1, yc1, w, h);    

        HCPT hcpt = new HCPT(img1, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, 
            nBins);
        
        p1 = new PatchUtil(w2, h2, nBins);
        p2 = new PatchUtil(w2, h2, nBins);
        p3 = new PatchUtil(w2, h2, nBins);
        p4 = new PatchUtil(w2, h2, nBins);
        p1.add(s1, hcpt);
        p2.add(s2, hcpt);
        p3.add(s3, hcpt);
        p4.add(s4, hcpt);
        ps = new PatchUtil[]{p1, p2, p3, p4};
        
        System.out.println("\nHCPT:");
        System.out.println("p1=" + Arrays.toString(p1.getHistogram()));
        System.out.println("p2=" + Arrays.toString(p2.getHistogram()));
        System.out.println("p3=" + Arrays.toString(p3.getHistogram()));
        System.out.println("p4=" + Arrays.toString(p4.getHistogram()));
        
        for (int i = 0; i < labels.length; ++i) {
            for (int j = (i + 1); j < labels.length; ++j) {
                double inter = ps[i].intersection(ps[j]);
                System.out.format("%s intersection %s = %.3f\n", 
                    labels[i], labels[j], (float)inter);
            }
        }
        for (int i = 0; i < labels.length; ++i) {
            for (int j = (i + 1); j < labels.length; ++j) {
                double[] diffAndErr = ps[i].diff(ps[j]);
                System.out.format("%s difference %s = %.3f(+-%.3f)\n", 
                    labels[i], labels[j], (float)diffAndErr[0], (float)diffAndErr[1]);
            }
        }
        
        // ==== HGS =====
        img1 = gsImg.subImage(xc1, yc1, w, h);    

        HGS hgs = new HGS(img1, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, 
            nBins);
        
        p1 = new PatchUtil(w2, h2, nBins);
        p2 = new PatchUtil(w2, h2, nBins);
        p3 = new PatchUtil(w2, h2, nBins);
        p4 = new PatchUtil(w2, h2, nBins);
        p1.add(s1, hgs);
        p2.add(s2, hgs);
        p3.add(s3, hgs);
        p4.add(s4, hgs);
        ps = new PatchUtil[]{p1, p2, p3, p4};
        
        System.out.println("\nHGS:");
        System.out.println("p1=" + Arrays.toString(p1.getHistogram()));
        System.out.println("p2=" + Arrays.toString(p2.getHistogram()));
        System.out.println("p3=" + Arrays.toString(p3.getHistogram()));
        System.out.println("p4=" + Arrays.toString(p4.getHistogram()));
        
        for (int i = 0; i < labels.length; ++i) {
            for (int j = (i + 1); j < labels.length; ++j) {
                double inter = ps[i].intersection(ps[j]);
                System.out.format("%s intersection %s = %.3f\n", 
                    labels[i], labels[j], (float)inter);
            }
        }
        for (int i = 0; i < labels.length; ++i) {
            for (int j = (i + 1); j < labels.length; ++j) {
                double[] diffAndErr = ps[i].diff(ps[j]);
                System.out.format("%s difference %s = %.3f(+-%.3f)\n", 
                    labels[i], labels[j], (float)diffAndErr[0], (float)diffAndErr[1]);
            }
        }
    }
}
