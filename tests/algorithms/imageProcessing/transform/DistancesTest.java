package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.transform.EpipolarTransformer.NormalizedXY;
import algorithms.matrix.MatrixUtil;
import algorithms.imageProcessing.transform.Util;
import algorithms.util.FormatArray;
import algorithms.util.ResourceFinder;
import algorithms.util.PairIntArray;
import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;

/**
 *
 * @author nichole
 */
public class DistancesTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public DistancesTest() {
    }
    
    public void test0() throws Exception {
        System.out.println("test0");
        
        PairIntArray leftTrueMatches = new PairIntArray();
        PairIntArray rightTrueMatches = new PairIntArray();
        getMertonCollegeMoreThan7TrueMatches(leftTrueMatches, rightTrueMatches);

        DenseMatrix xy1 = Util.rewriteInto3ColumnMatrix(leftTrueMatches);
        DenseMatrix xy2 = Util.rewriteInto3ColumnMatrix(rightTrueMatches);

        boolean useToleranceAsStatFactor = true;
        final double tolerance = 3.8;
        ErrorType errorType = ErrorType.SAMPSONS;
        
        NormalizedXY normXY1 = EpipolarTransformer.normalize(xy1);
        NormalizedXY normXY2 = EpipolarTransformer.normalize(xy2);
        DenseMatrix leftM = normXY1.getXy();
        DenseMatrix rightM = normXY2.getXy();
        
        EpipolarTransformer spTransformer = new EpipolarTransformer();

        DenseMatrix normalizedFM 
            = spTransformer.calculateEpipolarProjection(leftM, rightM);

        assertNotNull(normalizedFM);
                
        DenseMatrix fm = EpipolarTransformer.denormalizeTheFundamentalMatrix(
            normalizedFM, normXY1.getNormalizationMatrix(),
            normXY2.getNormalizationMatrix());
        
        Distances distances = new Distances();
        EpipolarTransformationFit fit = null;
        if (useToleranceAsStatFactor) {
            fit = distances.calculateError2(normalizedFM,
                leftM, rightM, errorType, tolerance);
        } else {
            fit = distances.calculateError(normalizedFM,
                leftM, rightM, errorType, tolerance);
        }
        
        assertEquals(leftTrueMatches.getN(), fit.getInlierIndexes().size());
        
        System.out.println("fm normalized (>7pts) =" + 
            FormatArray.toString(normalizedFM, "%.3e"));
        
        System.out.println("fm de-normalized (>7pts) =" + 
            FormatArray.toString(fm, "%.3e"));

        
        String fileName1 = "merton_college_I_001.jpg";
        String fileName2 = "merton_college_I_002.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);

        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        int image1Width = img1.getWidth();
        int image1Height = img1.getHeight();
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        int image2Width = img2.getWidth();
        int image2Height = img2.getHeight();

        overplotEpipolarLines(fm, xy1, xy2,
            img1, img2,
            image1Width, image1Height, image2Width, image2Height,
            Integer.valueOf(0).toString());
    }
    
    public void test7() throws Exception {
        System.out.println("test7");
        
        PairIntArray leftTrueMatches = new PairIntArray();
        PairIntArray rightTrueMatches = new PairIntArray();
        
        getMertonCollege7TrueMatches(leftTrueMatches, rightTrueMatches);
        
        DenseMatrix xy1 = Util.rewriteInto3ColumnMatrix(leftTrueMatches);
        DenseMatrix xy2 = Util.rewriteInto3ColumnMatrix(rightTrueMatches);
        
        NormalizedXY normXY1 = EpipolarTransformer.normalize(xy1);
        NormalizedXY normXY2 = EpipolarTransformer.normalize(xy2);
        DenseMatrix leftM = normXY1.getXy();
        DenseMatrix rightM = normXY2.getXy();
        
        EpipolarTransformer spTransformer = new EpipolarTransformer();

        List<DenseMatrix> normalizedFMs 
            = spTransformer.calculateEpipolarProjectionFor7Points(leftM, rightM);

        assertNotNull(normalizedFMs);
        assertFalse(normalizedFMs.isEmpty());
         
        
        boolean useToleranceAsStatFactor = true;
        final double tolerance = 3.8;
        ErrorType errorType = ErrorType.DIST_TO_EPIPOLAR_LINE;
        Distances distances = new Distances();
        
        List<DenseMatrix> denormalizedFMs = new ArrayList<DenseMatrix>();
        List<EpipolarTransformationFit> fits = new ArrayList<EpipolarTransformationFit>();
        
        EpipolarTransformationFit bestFit = null;
        int bestFitIdx = -1;
        
        for (DenseMatrix nfm : normalizedFMs) {
            DenseMatrix fm = EpipolarTransformer.denormalizeTheFundamentalMatrix(
                nfm, normXY1.getNormalizationMatrix(),
                normXY2.getNormalizationMatrix());
            denormalizedFMs.add(fm);
            
            EpipolarTransformationFit fit = null;
            if (useToleranceAsStatFactor) {
                fit = distances.calculateError2(nfm,
                        leftM, rightM, errorType, tolerance);
            } else {
                fit = distances.calculateError(nfm,
                        leftM, rightM, errorType, tolerance);
            }
            
            fits.add(fit);
            
            if (fit.isBetter(bestFit)) {
                bestFit = fit;
                bestFitIdx = fits.size() - 1;
            }
        }
        
        assertEquals(leftTrueMatches.getN(), bestFit.getInlierIndexes().size());
                
        System.out.println("fm normalized (>7pts) =" + 
            FormatArray.toString(normalizedFMs.get(bestFitIdx), "%.3e"));
        
        System.out.println("fm de-normalized (>7pts) =" + 
            FormatArray.toString(denormalizedFMs.get(bestFitIdx), "%.3e"));
        
        //fm 7pt= -0.02  0.40  0.61
        //        -0.30 -0.03 -3.30
        //        -0.63  3.27  0.07
                        
        assertEquals(leftTrueMatches.getN(), bestFit.getInlierIndexes().size());
        
        String fileName1 = "merton_college_I_001.jpg";
        String fileName2 = "merton_college_I_002.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        int image1Width = img1.getWidth();
        int image1Height = img1.getHeight();
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        int image2Width = img2.getWidth();
        int image2Height = img2.getHeight();
        
        overplotEpipolarLines(normalizedFMs.get(bestFitIdx), xy1, xy2,
            img1, img2, image1Width, image1Height, image2Width, image2Height, 
            "_merton_7pt_"); 
    }
    
    /*
    for more datasets:
    http://www.robots.ox.ac.uk/~vgg/data/data-mview.html
    */
    
    private void overplotEpipolarLines(DenseMatrix fm, DenseMatrix xy1,
        DenseMatrix xy2, Image img1, Image img2, int image1Width,
        int image1Height, int image2Width, int image2Height, String outfileNumber) 
        throws IOException {
        
        for (int ii = 0; ii < xy1.numColumns(); ii++) {
            double x = xy1.get(0, ii);
            double y = xy1.get(1, ii);
            ImageIOHelper.addPointToImage((float) x, (float) y, img1, 3,
                255, 0, 0);
        }
        for (int ii = 0; ii < xy2.numColumns(); ii++) {
            double x2 = xy2.get(0, ii);
            double y2 = xy2.get(1, ii);
            ImageIOHelper.addPointToImage((float) x2, (float) y2, img2, 3,
                255, 0, 0);
        }

        EpipolarTransformer spTransformer = new EpipolarTransformer();

        Color clr = null;
        for (int ii = 0; ii < xy2.numColumns(); ii++) {
            clr = getColor(clr);
            DenseMatrix epipolarLinesInLeft = 
                MatrixUtil.multiply(
                algorithms.matrix.MatrixUtil.transpose(fm), xy2);
            PairIntArray leftLine = spTransformer.getEpipolarLine(
                epipolarLinesInLeft, image1Width, image1Height, ii);
            ImageIOHelper.addCurveToImage(leftLine, img1, 0,
                clr.getRed(), clr.getGreen(), clr.getBlue());
        }

        clr = null;
        for (int ii = 0; ii < xy1.numColumns(); ii++) {
            clr = getColor(clr);
            DenseMatrix epipolarLinesInRight = MatrixUtil.multiply(fm, xy1);
            PairIntArray rightLine = spTransformer.getEpipolarLine(
                epipolarLinesInRight, img2.getWidth(), img2.getHeight(), ii);
            ImageIOHelper.addCurveToImage(rightLine, img2, 0,
                clr.getRed(), clr.getGreen(), clr.getBlue());
        }

        String dirPath = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(
            dirPath + "/tmp_m_1_" + outfileNumber + ".png", img1);
        ImageIOHelper.writeOutputImage(
            dirPath + "/tmp_m_2_" + outfileNumber + ".png", img2);
    }
    
    private Color getColor(Color clr) {
        if ((clr == null) || clr.equals(Color.MAGENTA)) {
            return Color.BLUE;
        }
        if (clr.equals(Color.BLUE)) {
            return Color.PINK;
        } else if (clr.equals(Color.PINK)) {
            return Color.GREEN;
        } else if (clr.equals(Color.GREEN)) {
            return Color.RED;
        } else if (clr.equals(Color.RED)) {
            return Color.CYAN;
        } else if (clr.equals(Color.CYAN)) {
            return Color.MAGENTA;
        } else if (clr.equals(Color.MAGENTA)) {
            return Color.LIGHT_GRAY;
        } else {
            return Color.ORANGE;
        }
    }
    
    protected void getMertonCollegeMoreThan7TrueMatches(PairIntArray left, 
        PairIntArray right) {
       
        left.add(58, 103);  right.add(32, 100);
        left.add(486, 46);   right.add(474, 49);
        left.add(845, 127);   right.add(878, 151);
        left.add(949, 430);   right.add(998, 471);
        left.add(541, 428);   right.add(533, 460);
        left.add(225, 453);   right.add(213, 498);
        left.add(49, 509);   right.add(21, 571);
        left.add(373, 239);   right.add(365, 258);
        left.add(737, 305);   right.add(762, 335);
        left.add(84, 273);   right.add(60, 298);
     
        
        left.add(134, 24);   right.add(116, 12);
        left.add(134, 66);   right.add(116, 60);
        left.add(106, 204);   right.add(82, 220);
        left.add(104, 276);   right.add(84, 299);
        left.add(110, 371);   right.add(90, 410);
        left.add(100, 466);   right.add(76, 518);
        left.add(715, 286);   right.add(734, 316);
        left.add(736, 488);   right.add(761, 530);
        left.add(80, 443);   right.add(51, 493);
        left.add(204, 415);   right.add(189, 456);
        left.add(507, 481);   right.add(501, 517);
        left.add(817, 32);   right.add(850, 49);
        
    }
    
    protected void getMertonCollege7TrueMatches(PairIntArray left, 
        PairIntArray right) {
        
        /*
        58, 103   32, 100
        486, 46   474, 49
        845, 127  878, 151
        949, 430  998, 471
        541, 428  533, 460
        225, 453  213, 498
        //49, 509   21, 571        
        */
        
        left.add(58, 103);  right.add(32, 100);
        left.add(486, 46);   right.add(474, 49);
        left.add(845, 127);   right.add(878, 151);
        left.add(949, 430);   right.add(998, 471);
        left.add(541, 428);   right.add(533, 460);
        left.add(225, 453);   right.add(213, 498);
        left.add(49, 509);   right.add(21, 571);
        //left.add(26, 554);   right.add(35, 650);
        //left.add(312, 580);   right.add(414, 699);
    }
    
    public static void main(String[] args) {
        
        try {
            DistancesTest test = new DistancesTest();
            
            //test.testRANSAC();
                        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }
    }

}
