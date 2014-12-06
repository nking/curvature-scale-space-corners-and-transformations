package algorithms.imageProcessing;

import Jama.Matrix;
import Jama.SingularValueDecomposition;
import algorithms.KSelect;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import algorithms.misc.MiscMath;
import algorithms.util.LinearRegression;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class StereoProjectionTransformerTest {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public StereoProjectionTransformerTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    /*
    for more datasets:
    http://www.robots.ox.ac.uk/~vgg/data/data-mview.html
    */
    
    @Test
    public void testB() throws Exception {
        
        //String cwd = System.getProperty("user.dir") + "/";
        
        PairIntArray matched1 = new PairIntArray();
        PairIntArray matched2 = new PairIntArray();
        
        readBrownAndLoweMatches(matched1, matched2);
        
        PairIntArray xy1 = new PairIntArray();
        PairIntArray xy2 = new PairIntArray();
        
        readBrownAndLoweCorners(xy1, xy2);
        
       
        // diff in matched points and stdev
        LinearRegression linReg = new LinearRegression();
        
        linReg.plotTheLinearRegression(xy1, xy2);
    }
   
    public void testC() throws Exception {
        
        String cwd = System.getProperty("user.dir") + "/";
        
        //String fileName1 = "venturi_mountain_j6_0001.png";
        //String fileName1 = "lab.gif";
        String fileName1 = "brown_lowe_2003_image1.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        GreyscaleImage img1 = ImageIOHelper.readImageAsGrayScaleB(filePath1);
        String fileName2 = "brown_lowe_2003_image2.jpg";
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        //GreyscaleImage img2 = ImageIOHelper.readImageAsGrayScaleB(filePath2);
        
        /*
        CurvatureScaleSpaceInflectionMapper mapper = new 
            CurvatureScaleSpaceInflectionMapper(img1, img2);
        
        mapper.useOutdoorMode();

        TransformationParameters transformationParams
            = mapper.createEuclideanTransformation();
        */
        
        /*
         * wanting to edit contours to make contours more likely:
        CurvatureScaleSpaceImageMaker:
           -- cannyedge:  highThresh = 1 or 2 * low
           --             then blur w gaussian: 5, 10, 15, 20
        */        
        
        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img1);
        
        detector.useOutdoorMode();
        
        //detector.useSegmentationForSky();
        //detector.useLowestHighIntensityCutoff();
                       
        detector.findCorners();

        List<PairIntArray> edges = detector.getEdgesInOriginalReferenceFrame();
        
        Image image = ImageIOHelper.readImageAsGrayScale(filePath1);
        
        Image image2 = new Image(image.getWidth(), image.getHeight());
                                  
        ImageIOHelper.addAlternatingColorCurvesToImage(edges, image2);        
        
        String outFilePath = cwd + "tmp_edges.png";
                 
        ImageIOHelper.writeOutputImage(outFilePath, image2);
                 
                                                          
        outFilePath = cwd + "tmp.png";
                 
        ImageIOHelper.writeOutputImage(outFilePath, image);
        
        StringBuilder sb = new StringBuilder();
        PairIntArray corners = detector.getCornersInOriginalReferenceFrame();
        for (int i = 0; i < corners.getN(); i++) {
            int x = corners.getX(i);
            int y = corners.getY(i);
            int xe = (int)Math.sqrt(x);
            int ye = (int)Math.sqrt(y);
            sb.append(String.format("%d\t%d\t%d\t%d\n", x, y, xe, ye));
        }
        ResourceFinder.writeToCWD(sb.toString(), "tmp2.tsv");
        
    }
    
    public void estD() throws Exception {
        
        String fileName1 = "books_illum3_v0_695x555.png";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        GreyscaleImage img1 = ImageIOHelper.readImageAsGrayScaleG(filePath1);

/*        String fileName2 = "books_illum3_v6_695x555.png";
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        GreyscaleImage img2 = ImageIOHelper.readImageAsGrayScaleG(filePath2);

        CurvatureScaleSpaceInflectionMapper mapper = new 
            CurvatureScaleSpaceInflectionMapper(img1, img2);
       
        mapper.useLineDrawingLineMode();

        TransformationParameters transformationParams
            = mapper.createEuclideanTransformation();
  */
        
        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img1);
                       
        detector.findCorners();

        
        List<PairIntArray> edges = detector.getEdgesInOriginalReferenceFrame();
        
        Image image = ImageIOHelper.readImageAsGrayScale(filePath1);
        
        
        Image image2 = new Image(image.getWidth(), image.getHeight());
                                  
        ImageIOHelper.addAlternatingColorCurvesToImage(edges, image2);
        
        String cwd = System.getProperty("user.dir") + "/";
        
        String outFilePath = cwd + "tmp_edges.png";
                 
        ImageIOHelper.writeOutputImage(outFilePath, image2);
                 
                                                          
        outFilePath = cwd + "tmp.png";
                 
        ImageIOHelper.writeOutputImage(outFilePath, image);
        
        int z = 1;
    }
    
    public void testE() throws Exception {
        
        PairFloatArray matched1 = new PairFloatArray();
        PairFloatArray matched2 = new PairFloatArray();
        
        readBrownAndLoweMatches(matched1, matched2);
        
        StereoProjectionTransformer spTransformer = 
            new StereoProjectionTransformer();
        
        spTransformer.calculateEpipolarProjection(matched1, matched2);
        
        double[] leftEpipole = spTransformer.getLeftEpipole();
        
        double[] rightEpipole = spTransformer.getRightEpipole();
      
        String fileName1 = "brown_lowe_2003_image1.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        
        String fileName2 = "brown_lowe_2003_image2.jpg";
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        
        double rightEpipoleX = rightEpipole[0] / rightEpipole[2];
        double rightEpipoleY = rightEpipole[1] / rightEpipole[2];
        double leftEpipoleX = leftEpipole[0] / leftEpipole[2];
        double leftEpipoleY = leftEpipole[1] / leftEpipole[2];
        log.info("leftEpipole=" + Arrays.toString(leftEpipole));
        log.info("rightEpipole=" + Arrays.toString(rightEpipole));
        log.info("leftEpipole(X,Y) = ("  + leftEpipoleX + ", " + leftEpipoleY + ")");
        log.info("rightEpipole(X,Y) = ("  + rightEpipoleX + ", " + rightEpipoleY + ")");
        PairIntArray leftMatches = new PairIntArray();
        PairIntArray rightMatches = new PairIntArray();
        
        int nLimitTo = matched1.getN();
        
        for (int i = 0; i < nLimitTo; i++) {
            leftMatches.add((int)Math.round(matched1.getX(i)),
                (int)Math.round(matched1.getY(i)));
            rightMatches.add((int)Math.round(matched2.getX(i)),
                (int)Math.round(matched2.getY(i)));
        }
        
        Color clr = null;
        PairIntArray subsetLeft = new PairIntArray();
        for (int i = 0; i < nLimitTo; i++) {
            
            clr = getColor(clr);
            
            PairIntArray leftLine = spTransformer.getEpipolarLineInLeft(
                img1.getWidth(), i);
            
            ImageIOHelper.addCurveToImage(leftLine, img1, 0, 
                clr.getRed(), clr.getGreen(), clr.getBlue());
            
            subsetLeft.add(Math.round(matched1.getX(i)), Math.round(matched1.getY(i)));
        }
        ImageIOHelper.addCurveToImage(subsetLeft, img1, 2, 255, 0, 0);
        
        clr = null;
        PairIntArray subsetRight = new PairIntArray();
        for (int i = 0; i < nLimitTo; i++) {
            clr = getColor(clr);
            
            PairIntArray rightLine = spTransformer.getEpipolarLineInRight(
                img2.getWidth(), i);
            
            ImageIOHelper.addCurveToImage(rightLine, img2, 0, 
                clr.getRed(), clr.getGreen(), clr.getBlue()); 
            
            subsetRight.add(Math.round(matched2.getX(i)), Math.round(matched2.getY(i)));
        }
        ImageIOHelper.addCurveToImage(subsetRight, img2, 2, 255, 0, 0);
        
        
        String dirPath = ResourceFinder.findDirectory("bin");
        
        ImageIOHelper.writeOutputImage(
            dirPath + "/image1_matched_corners.png", img1);
       
        ImageIOHelper.writeOutputImage(
            dirPath + "/image2_epipolar_and_matches.png", img2);
        
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
    
    public void testCalculateEpipolarProjection() {
        
        PairFloatArray leftXY = new PairFloatArray();
        
        PairFloatArray rightXY = new PairFloatArray();
        
        populateWithTestData0(leftXY, rightXY);
        
        StereoProjectionTransformer transformer = 
            new StereoProjectionTransformer();
        
        transformer.calculateEpipolarProjection(leftXY, rightXY);
    }
    
    private void populateWithTestData0(PairFloatArray leftXY, 
        PairFloatArray rightXY) {
        
        /*
        extracted from the books V0 and V6 images of the middlebury vision
        database.
        http://vision.middlebury.edu/stereo/data/
        
        The References listed with the data are:
        
        D. Scharstein and C. Pal. "Learning conditional random fields for 
        stereo."  In IEEE Computer Society Conference on Computer Vision and 
        Pattern Recognition (CVPR 2007), Minneapolis, MN, June 2007.

        H. Hirschmüller and D. Scharstein. "Evaluation of cost functions for 
        stereo matching."  In IEEE Computer Society Conference on Computer 
        Vision and Pattern Recognition (CVPR 2007), Minneapolis, MN, June 2007.
        
        */
        
    }
    
    private void readBrownAndLoweMatches(PairIntArray imageCorners1XY,
        PairIntArray imageCorners2XY) throws IOException {
        
        String fileName1 = "brown_lowe_2003_matching.tsv";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        
        BufferedReader br = null;
        FileReader reader = null;

        try {
            reader = new FileReader(new File(filePath1));
            br = new BufferedReader(reader);
            String line = br.readLine();
            line = br.readLine();
            while (line != null) {
                String[] items = line.split("\\s+");
                if ((items != null) && (items.length == 4)) {
                    
                    Integer x1 = Integer.valueOf(items[0]);
                    Integer y1 = Integer.valueOf(items[1]);
                    Integer x2 = Integer.valueOf(items[2]);
                    Integer y2 = Integer.valueOf(items[3]);
                    imageCorners1XY.add(x1.intValue(), y1.intValue());
                    imageCorners2XY.add(x2.intValue(), y2.intValue());
//System.out.println(x1 + ", " + y1 + "  " + x2 + ", " + y2);                    
                }
                line = br.readLine();
            }
        } finally {
            if (reader != null) {
                reader.close();
            }
            if (br != null) {
                br.close();
            }
        }
    }
    
    private void readBrownAndLoweMatches(PairFloatArray imageCorners1XY,
        PairFloatArray imageCorners2XY) throws IOException {
        
        String fileName1 = "brown_lowe_2003_matching.tsv";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        
        BufferedReader br = null;
        FileReader reader = null;
        try {
            reader = new FileReader(new File(filePath1));
            br = new BufferedReader(reader);
            String line = br.readLine();
            line = br.readLine();
            while (line != null) {
                String[] items = line.split("\\s+");
                if ((items != null) && (items.length == 4)) {
                    Integer x1 = Integer.valueOf(items[0]);
                    Integer y1 = Integer.valueOf(items[1]);
                    Integer x2 = Integer.valueOf(items[2]);
                    Integer y2 = Integer.valueOf(items[3]);
                    imageCorners1XY.add(x1.intValue(), y1.intValue());
                    imageCorners2XY.add(x2.intValue(), y2.intValue());
//System.out.println(x1 + ", " + y1 + "  " + x2 + ", " + y2);
                }
                line = br.readLine();
            }
        } finally {
            if (reader != null) {
                reader.close();
            }
            if (br != null) {
                br.close();
            }
        }
    }

    private void readBrownAndLoweCorners(PairIntArray imageCorners1XY,
        PairIntArray imageCorners2XY) throws IOException {
                
        BufferedReader br = null;
        FileReader reader = null;
        
         String[] fileNames = new String[]{"brown_lowe_2003_image1.tsv", 
            "brown_lowe_2003_image2.tsv"};
        
        for (String fileName : fileNames) {
            
            String filePath1 = ResourceFinder.findFileInTestResources(fileName);
            PairIntArray xy = fileName.equals("brown_lowe_2003_image1.tsv") 
                ? imageCorners1XY : imageCorners2XY;
            
            try {
                reader = new FileReader(new File(filePath1));
                br = new BufferedReader(reader);
                String line = br.readLine();
                while (line != null) {
                    String[] items = line.split("\\s+");
                    if ((items != null) && (items.length == 2)) {
                        Integer x1 = Integer.valueOf(items[0]);
                        Integer y1 = Integer.valueOf(items[1]);
                        xy.add(x1.intValue(), y1.intValue());
                    }
                    line = br.readLine();
                }
            } finally {
                if (reader != null) {
                    reader.close();
                }
                if (br != null) {
                    br.close();
                }
            }
        }
        
    }
    
    public static void main(String[] args) {
        
        try {
            StereoProjectionTransformerTest test = 
                new StereoProjectionTransformerTest();
            
            //test.testB();
            
            //test.testC();
            
            test.testE();
                        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }
    }

}
