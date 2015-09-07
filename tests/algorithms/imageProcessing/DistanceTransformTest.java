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

    public void testApplyTransforms3() throws Exception {

        // snapshot of spiral at http://www.agencia.fapesp.br/arquivos/survey-final-fabbri-ACMCSurvFeb2008.pdf
        
        String filePath = ResourceFinder.findFileInTestResources("spiral.png");
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
        
        ImageIOHelper.writeOutputImage(bin + "/m_dt_spiral.png", outputImgM);
        
        GreyscaleImage outputImgM2 = new GreyscaleImage(w, h);
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                outputImgM2.setValue(i, j, dtM[i][j]);
            }
        }
        
        ImageIOHelper.writeOutputImage(bin + "/m_dt_spiral_2.png", outputImgM2);
        
        
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
        ImageIOHelper.writeOutputImage(bin + "/dt_inv_input_spiral.png", imgInv0);

        System.out.println("INV:");
        int[][] dtInvM = dtr.applyMeijsterEtAl(pointsInvM, w, h);
        
        GreyscaleImage outputImgInvM = ImageIOHelper.scaleToImgRange(dtInvM);
       
        ImageIOHelper.writeOutputImage(bin + "/m_inv_dt_spiral.png", outputImgInvM);
    }
    
    public void testApplyTransforms4() throws Exception {

        int w = 9;
        int h = 9;
        
        Set<PairInt> pointsM = new HashSet<PairInt>();
        double[][] data = new double[w][h];
        for (int x = 0; x < w; ++x) {
            data[x] = new double[h];
        }
        
        /*
        making the pattern from Fig 1.a of Fabbri et al.
           0  1  2  3  4  5  6  7  8 
        0  -  -  -  -  -  -  -  -  - 
        1  -  -  -  -  -  -  -  -  - 
        2  -  -  -  1  1  1  -  -  -
        3  -  -  1  1  1  1  1  -  -
        4  -  -  1  1  1  1  1  -  -
        5  -  -  1  1  1  1  1  -  -
        6  -  -  -  1  1  1  -  -  -
        7  -  -  -  -  -  -  -  1  -
        8  -  -  -  -  -  -  1  -  -
           0  1  2  3  4  5  6  7  8  9
        */
        for (int x = 3; x <= 5; ++x) {
            for (int y = 2; y <= 6; ++y) {
                pointsM.add(new PairInt(x, y));
                data[x][y] = 1;
            }
        }
        for (int x = 2; x <= 2; ++x) {
            for (int y = 3; y <= 5; ++y) {
                pointsM.add(new PairInt(x, y));
                data[x][y] = 1;
            }
        }
        for (int x = 6; x <= 6; ++x) {
            for (int y = 3; y <= 5; ++y) {
                pointsM.add(new PairInt(x, y));
                data[x][y] = 1;
            }
        }
        pointsM.add(new PairInt(6, 8));
        data[6][8] = 1;
        pointsM.add(new PairInt(7, 7));
        data[7][7] = 1;
            
        DistanceTransform dtr = new DistanceTransform();
        
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
        
        /*
        StringBuilder sb = new StringBuilder("original data:\n");
        for (int j = 0; j < h; ++j) {
            sb.append("row ").append(j).append(": ");
            for (int i = 0; i < w; ++i) {
                sb.append(" ").append(imgInv0.getValue(i, j)/255);
            }
            sb.append("\n");
        }
        System.out.println(sb.toString());
        */
        
        int[][] dtInvM = dtr.applyMeijsterEtAl(pointsInvM, w, h);
        
        /*StringBuilder sb2 = new StringBuilder();
        for (int j = 0; j < h; ++j) {
            sb2.append("row ").append(j).append(": ");
            for (int i = 0; i < w; ++i) {
                sb2.append(" ").append(dtInvM[i][j]);
            }
            sb2.append("\n");
        }
        System.out.println(sb2.toString());
        */
            
        int[][] expected = new int[9][9];
        expected[2] = new int[]{0, 0, 0, 1, 1, 1, 0, 0, 0};
        expected[3] = new int[]{0, 0, 1, 2, 4, 2, 1, 0, 0};
        expected[4] = new int[]{0, 0, 1, 4, 8, 4, 1, 0, 0};
        expected[5] = new int[]{0, 0, 1, 2, 4, 2, 1, 0, 0};
        expected[6] = new int[]{0, 0, 0, 1, 1, 1, 0, 0, 1};
        expected[7] = new int[]{0, 0, 0, 0, 0, 0, 0, 1, 0};
        expected[8] = new int[]{0, 0, 0, 0, 0, 0, 0, 0, 0};
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int a = dtInvM[i][j];
                int b = expected[i][j];
                assertTrue(a == b);
            }
        }
    }
}
