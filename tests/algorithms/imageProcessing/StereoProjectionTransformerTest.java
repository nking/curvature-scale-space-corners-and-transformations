package algorithms.imageProcessing;

import Jama.Matrix;
import Jama.SingularValueDecomposition;
import algorithms.KSelect;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import algorithms.misc.MiscMath;
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
        int n = matched1.getN();
        long sumX = 0;
        long sumY = 0;
        int[] dx = new int[n];
        int[] dy = new int[n];
        int dxMin = Integer.MAX_VALUE;
        int dxMax = Integer.MIN_VALUE;
        int dyMin = Integer.MAX_VALUE;
        int dyMax = Integer.MIN_VALUE;
        for (int i = 0; i < n; i++) {
            
            int diffX = matched1.getX(i) - matched2.getX(i);
            int diffY = matched1.getY(i) - matched2.getY(i);
            sumX += diffX;
            sumY += diffY;
            
            dx[i] = diffX;
            dy[i] = diffY;
            
            if (dxMin > diffX) {
                dxMin = diffX;
            }
            if (dyMin > diffY) {
                dyMin = diffY;
            }
            if (dxMax < diffX) {
                dxMax = diffX;
            }
            if (dyMax < diffY) {
                dyMax = diffY;
            }
            log.info("dx=" + dx[i] + "  dy=" + dy[i]);
        }
        double transX = (double)sumX/(double)matched1.getN();
        double transY = (double)sumY/(double)matched1.getN();
        
        sumX = 0;
        sumY = 0;
        for (int i = 0; i < matched1.getN(); i++) {
            double diffX = dx[i] - transX;
            double diffY = dy[i] - transY;
            sumX += (diffX * diffX);
            sumY += (diffY * diffY);
        }
        double stdDevX = (Math.sqrt(sumX/(n - 1.0f)));
        double stdDevY = (Math.sqrt(sumY/(n - 1.0f)));
        log.info("transX=" + transX + " (stdDev=" + stdDevX + ")");
        log.info("transY=" + transY + " (stdDev=" + stdDevY + ")");
        
        /*
        transX = 293.0882352941176 (stdDev=10.263203002789375)
        transY = 14.323529411764707 (stdDev=5.898125121074294)
      
        for 1000 points, for each possible pair w/ image 2 points,
        the real solution would be looking for a match within 
        2.5*stdev        
        */
        for (int i = 0; i < n; i++) {
            double diffX = Math.abs((dx[i] - transX));
            double diffY = Math.abs((dy[i] - transY));
            log.info("diffX=" + diffX + "  diffY=" + diffY);
            assertTrue(diffX <= 3.0*stdDevX);
            assertTrue(diffY <= 3.0*stdDevY);
        }
        
        /* linear regression w/ theil sen estimator:
        http://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator
        computing all O(n^2) lines and and then applying a linear time 
        median finding algorithm
        */
        int count = 0;
        
        float[] s = new float[n*n];
        for (int i = 0; i < n; i++) {
            if (dx[i] == 0) {
                continue;
            }
            s[count] = (float)(dy[i])/(float)(dx[i]);
            count++;
        }
        
        KSelect kSelect = new KSelect();
        float median = kSelect.findMedianOfMedians(s, 0, count - 1);
        log.info("thiel sen beta=" + median 
            + " transY/transX=" + transY/transX);
       
        int[] x1T = new int[xy1.getN()];
        int[] y1T = new int[xy1.getN()];
        for (int i = 0; i < xy1.getN(); i++) {
            x1T[i] = (int) (xy1.getX(i) - transX);
            y1T[i] = (int) (xy1.getY(i) - transY);
        }
        int[] x2 = new int[xy2.getN()];
        int[] y2 = new int[xy2.getN()];
        for (int i = 0; i < xy2.getN(); i++) {
            x2[i] = xy2.getX(i);
            y2[i] = xy2.getY(i);
        }
        
        /*
        y = y0 + b*(x-x0)
        y = y1[0] + deltay*(xMin - x1[0]);
        
        y = y1[0] + median*(0 - x1[0]);
        */
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(dxMin, dxMax, dyMin, dyMax, dx, dy, 
                    new int[]{0, dxMax}, 
                    new int[]{(int)((dy[0] + median*(0 - dx[0]))),
                    (int)((dy[0] + median*(dxMax - dx[0])))},
                    "diff x (img1 - img2) vs diff y");
        
        int xMax = MiscMath.findMax(x1T);
        int yMax1 = MiscMath.findMax(y1T);
        int yMax2 = MiscMath.findMax(y2);
        float yMin1 = MiscMath.findMin(y1T);
        float yMin2 = MiscMath.findMin(y2);
        int yMax = (yMax1 > yMax2) ? yMax1 : yMax2;
        float yMin = (yMin1 < yMin2) ? yMin1 : yMin2;
        plotter.addPlot(0, xMax, yMin, yMax,
            x1T, y1T, 
            new int[0], new int[0], "img1 corners Transformed");
        
        plotter.addPlot(0, xMax, yMin, yMax,
            x2, y2, 
            new int[0], new int[0], "img2 corners");
        
        plotter.writeFile();
        
        // === now, given transX and transY, find points in image2 within
        // projected +- 3*stdDev
        int[] xt = new int[xy1.getN()];
        int[] yt = new int[xy1.getN()];
        int n2 = 0;
        int[] xmatched1 = new int[xy1.getN()];
        int[] ymatched1 = new int[xy1.getN()];
        int[] xfound2 = new int[xy1.getN()];
        int[] yfound2 = new int[xy1.getN()];
        int xPlot1Min = Integer.MAX_VALUE;
        int xPlot1Max = Integer.MIN_VALUE;
        int yPlot1Min = Integer.MAX_VALUE;
        int yPlot1Max = Integer.MIN_VALUE;
        
        int xPlot2Min = Integer.MAX_VALUE;
        int xPlot2Max = Integer.MIN_VALUE;
        int yPlot2Min = Integer.MAX_VALUE;
        int yPlot2Max = Integer.MIN_VALUE;
        
        for (int i = 0; i < xy1.getN(); i++) {
            xt[i] = (int)(xy1.getX(i) - transX);
            yt[i] = (int)(xy1.getY(i) - transY);
            
            double minDiff = Double.MAX_VALUE;
            int minIdx2 = -1;
            for (int j = 0; j < xy2.getN(); j++) {
                double diffX = Math.abs(xy2.getX(j) - xt[i]);
                if (diffX > 1.0*stdDevX) {
                    continue;
                }
                double diffY = Math.abs(xy2.getY(j) - yt[i]);
                if (diffY > 1.0*stdDevY) {
                    continue;
                }
                double diff = Math.sqrt(diffX*diffX + diffY*diffY);
                if (diff < minDiff) {
                    minDiff = diff;
                    minIdx2 = j;
                }
            }
            if (minIdx2 > -1) {
                xfound2[n2] = xy2.getX(minIdx2);
                yfound2[n2] = xy2.getY(minIdx2);
                xmatched1[n2] = xy1.getX(i);
                ymatched1[n2] = xy1.getY(i);
                
                if (xPlot1Min > xmatched1[n2]) {
                    xPlot1Min = xmatched1[n2];
                }
                if (xPlot1Max < xmatched1[n2]) {
                    xPlot1Max = xmatched1[n2];
                }
                if (yPlot1Min > ymatched1[n2]) {
                    yPlot1Min = ymatched1[n2];
                }
                if (yPlot1Max < ymatched1[n2]) {
                    yPlot1Max = ymatched1[n2];
                }
                
                if (xPlot2Min > xfound2[n2]) {
                    xPlot2Min = xfound2[n2];
                }
                if (xPlot2Max < xfound2[n2]) {
                    xPlot2Max = xfound2[n2];
                }
                if (yPlot2Min > yfound2[n2]) {
                    yPlot2Min = yfound2[n2];
                }
                if (yPlot2Max < yfound2[n2]) {
                    yPlot2Max = yfound2[n2];
                }
                
                n2++;
            }
        }
        
        plotter.addPlot(xPlot1Min, xPlot1Max, yPlot1Min, yPlot1Max,
            xmatched1, ymatched1, 
            new int[0], new int[0], "the matched image1 points");
    
        plotter.addPlot(xPlot2Min, xPlot2Max, yPlot2Min, yPlot2Max,
            xfound2, yfound2, 
            new int[0], new int[0], "the matched image2 points");
    
        plotter.writeFile();
        
        log.info("original list matches=" + matched1.getN() 
            + " final matches=" + n2);
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
        Image img1 = ImageIOHelper.readImage(filePath1);
        
        String fileName2 = "brown_lowe_2003_image2.jpg";
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        Image img2 = ImageIOHelper.readImage(filePath2);
        
        PairIntArray leftMatches = new PairIntArray();
        PairIntArray rightMatches = new PairIntArray();
        
        int nLimitTo = 5;//matched1.getN();
        
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
