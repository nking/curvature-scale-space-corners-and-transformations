package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class DistanceTransformTest extends TestCase {

    public DistanceTransformTest() {
    }

    public void testApplyTransforms1() throws Exception {

        String filePath = ResourceFinder.findFileInTestResources("two_circles.png");
        Image img = ImageIOHelper.readImageAsGrayScale(filePath);
        GreyscaleImage img0 = img.copyToGreyscale();

        int w = img0.getWidth();
        int h = img0.getHeight();
        
        Set<PairInt> pointsM = new HashSet<PairInt>();
        double[][] dataFH = new double[w][h];
        for (int x = 0; x < w; ++x) {
            dataFH[x] = new double[h];
        }
        
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                if (img0.getValue(x, y) > 0) { 
                    pointsM.add(new PairInt(x, y));
                    dataFH[x][y] = 1;
                }
            }
        }

        String bin = ResourceFinder.findDirectory("bin");

        DistanceTransform dtr = new DistanceTransform();
        int[][] dtM = dtr.applyMeijsterEtAl(pointsM, w, h);

        GreyscaleImage outputImgM = ImageIOHelper.scaleToImgRange(dtM);
        
        ImageIOHelper.writeOutputImage(bin + "/m_dt_circ.png", outputImgM);
        
        // ----- inverse binary of that  ----------
        GreyscaleImage imgInv0 = new GreyscaleImage(w, h);
        Set<PairInt> pointsInvM = new HashSet<PairInt>();
        
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                PairInt p = new PairInt(x, y);
                if (!pointsM.contains(p)) {
                    pointsInvM.add(new PairInt(x, y));
                    imgInv0.setValue(x, y, 255);
                }
            }
        }
        ImageIOHelper.writeOutputImage(bin + "/dt_inv_input_circ.png", imgInv0);

        System.out.println("INV:");
        int[][] dtInvM = dtr.applyMeijsterEtAl(pointsInvM, w, h);
        
        GreyscaleImage outputImgInvM = ImageIOHelper.scaleToImgRange(dtInvM);
       
        ImageIOHelper.writeOutputImage(bin + "/m_inv_dt_circ.png", outputImgInvM);
    }

    public void testApplyTransforms2() throws Exception {

         String filePath = ResourceFinder.findFileInTestResources(
            "blox_binary.png");

        Image img = ImageIOHelper.readImageAsGrayScale(filePath);
        GreyscaleImage img0 = img.copyToGreyscale();

        int w = img0.getWidth();
        int h = img0.getHeight();

        double[][] dataM = new double[w][h];

        Set<PairInt> nonZeros = new HashSet<PairInt>();

        for (int i = 0; i < w; ++i) {
            dataM[i] = new double[h];
            for (int j = 0; j < h; ++j) {
                int v = img0.getValue(i, j);
                if (v > 0) {
                    dataM[i][j] = v;
                    nonZeros.add(new PairInt(i, j));
                }
            }
        }

        DistanceTransform dtr = new DistanceTransform();

        int[][] outDataM = dtr.applyMeijsterEtAl(nonZeros, w, h);

        GreyscaleImage outputImgM = ImageIOHelper.scaleToImgRange(outDataM);

        String bin = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(bin + "/m_dt_blox_binary.png", outputImgM);

         // ----- inverse binary of that  ----------
        GreyscaleImage imgInv0 = new GreyscaleImage(w, h);
        Set<PairInt> pointsInvM = new HashSet<PairInt>();
        
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                PairInt p = new PairInt(x, y);
                if (!nonZeros.contains(p)) {
                    pointsInvM.add(new PairInt(x, y));
                    imgInv0.setValue(x, y, 255);
                }
            }
        }
        ImageIOHelper.writeOutputImage(bin + "/dt_inv_input_blox.png", imgInv0);

        int[][] dtInvM = dtr.applyMeijsterEtAl(pointsInvM, w, h);
        
        GreyscaleImage outputImgInvM = ImageIOHelper.scaleToImgRange(dtInvM);
       
        ImageIOHelper.writeOutputImage(bin + "/m_inv_dt_blox.png", outputImgInvM);
    }

}
