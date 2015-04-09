package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import algorithms.util.LinearRegression;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import java.awt.Color;
import java.util.List;
import java.util.logging.Logger;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.ejml.simple.*;

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
    
    //@Test
    public void estB() throws Exception {
        
        //String cwd = System.getProperty("user.dir") + "/";
        
        PairIntArray matched1 = new PairIntArray();
        PairIntArray matched2 = new PairIntArray();
        
        DataForTests.readBrownAndLoweMatches(matched1, matched2);
        
        PairIntArray xy1 = new PairIntArray();
        PairIntArray xy2 = new PairIntArray();
        
        DataForTests.readBrownAndLoweCorners(xy1, xy2);
        
       
        // diff in matched points and stdev
        LinearRegression linReg = new LinearRegression();
        
        linReg.plotTheLinearRegression(xy1, xy2);
    }
   
    //@Test
    public void estC() throws Exception {
        
        String cwd = System.getProperty("user.dir") + "/";
        
        //String fileName1 = "venturi_mountain_j6_0001.png";
        //String fileName1 = "lab.gif";
        String fileName1 = "brown_lowe_2003_image1.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
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
    
    //@Test
    public void estD() throws Exception {
        
        String fileName1 = "books_illum3_v0_695x555.png";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

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
    
    /* {        
        PairIntArray corners1Int = new PairIntArray();
        PairIntArray corners2Int = new PairIntArray();        
        DataForTests.readBrownAndLoweCorners(corners1Int, corners2Int);
        
        String fileName1 = "brown_lowe_2003_image1.jpg";
        String fileName2 = "brown_lowe_2003_image2.jpg";
       
        // temporary look at the edges that produced the corners:
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        GreyscaleImage img1 = ImageIOHelper.readImageAsGrayScaleB(filePath1);
        int image1Width = img1.getWidth();
        int image1Height = img1.getHeight();
       
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        GreyscaleImage img2 = ImageIOHelper.readImageAsGrayScaleB(filePath2);
        int image2Width = img2.getWidth();
        int image2Height = img2.getHeight();
        
        CurvatureScaleSpaceInflectionMapperForOpenCurves inflMapper = new
            CurvatureScaleSpaceInflectionMapperForOpenCurves(img1, img2);
        PairIntArray[] xyPeaks = inflMapper.createUnmatchedXYFromContourPeaks();        
        PairIntArray points1 = xyPeaks[0];
        PairIntArray points2 = xyPeaks[1];
        // there may be a couple redundant points per set
        DataForTests.writePointsToTestResources(points1, 
            "brown_lowe_2003_image1_infl_pts.tsv");        
        DataForTests.writePointsToTestResources(points2, 
            "brown_lowe_2003_image2_infl_pts.tsv");
        TransformationParameters params = inflMapper.createEuclideanTransformation();
    }
    */
    
    @Test
    public void testH() throws Exception {
        
        /*
        compare these calculated with the merton college I dataset:
            (0) epipolar projection solution with all points already matched.
                -- log the solution and plot the final 2 figures.
            (1) epipolar projection solution with 7 points already matched.
                -- log the solution and plot the final 2 figures.
            (2) epipolar projection solution with RANSAC on the already matched
                lists of points.
                -- log the solution and plot the final 2 figures.
            (3) epipolar projection solution for the stereo projection point
                matcher and solver.
                (not yet implemented).
                -- log the solution and plot the final 2 figures.
        
        Merton I College set is from
        http://www.robots.ox.ac.uk/~vgg/data/data-mview.html
        
        Note that this set was chosen to compare initial results with those from
        the "Programming Computer Vision with Python" by Jan Solem.
        */
                    
        String fileName1 = "merton_college_I_001.jpg";
        String fileName2 = "merton_college_I_002.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        
        // mRows = 3;  nCols = 484
        SimpleMatrix unnormXY1 = DataForTests.readMerton1UnnormalizedX1Data();
        // mRows = 3;  nCols = 484
        SimpleMatrix unnormXY2 = DataForTests.readMerton1UnnormalizedX2Data();
        
        for (int methodType = 0; methodType < 4; methodType++) {
        
            SimpleMatrix fm = null;
            SimpleMatrix input1 = null;
            SimpleMatrix input2 = null;
        
            StereoProjectionTransformer spTransformer = 
                new StereoProjectionTransformer();
            
            switch(methodType) {
                
                case 0:
                    // solution with all already matched points:
                    
                    fm = spTransformer
                        .calculateEpipolarProjectionForPerfectlyMatched(
                            unnormXY1, unnormXY2);
            
                    input1 = unnormXY1;
                    input2 = unnormXY2;
                    
                    break;
                    
                case 1:
                    // solution with 7 of the already matched points:
                    PairFloatArray xy1 = new PairFloatArray();
                    PairFloatArray xy2 = new PairFloatArray();
                    xy1.add(97.263000f, 466.206000f);
                    xy2.add(75.142000f, 519.539000f);
                    xy1.add(47.343000f, 33.856000f);
                    xy2.add(18.089000f, 19.859000f);
                    xy1.add(421.285000f, 44.954000f);
                    xy2.add(395.439000f, 45.559000f);
                    xy1.add(521.201000f, 493.853000f);
                    xy2.add(534.457000f, 538.037000f);
                    xy1.add(604.038000f, 239.050000f);
                    xy2.add(609.635000f, 260.938000f);
                    xy1.add(816.107000f, 33.034000f);
                    xy2.add(848.550000f, 50.388000f);
                    xy1.add(948.143000f, 344.338000f);
                    xy2.add(1000.549000f, 383.014000f);
                    
                    input1 = StereoProjectionTransformer
                        .rewriteInto3ColumnMatrix(xy1);
                    input2 = StereoProjectionTransformer
                        .rewriteInto3ColumnMatrix(xy2);
            
                    SimpleMatrix[] fms = spTransformer
                        .calculateEpipolarProjectionFor7Points(xy1, xy2);
                    if (fms != null) {
                        fm = fms[0];
                    }
                    
                    break;
                    
                case 2:
                    // use RANSAC on all points
                    
                    RANSACSolver solver = new RANSACSolver();
                    
                    PairFloatArray outputLeftXY = new PairFloatArray();
                    PairFloatArray outputRightXY = new PairFloatArray();
                    
                    input1 = unnormXY1;
                    input2 = unnormXY2;
                    
                    StereoProjectionTransformerFit fit = 
                        solver.calculateEpipolarProjection(
                        input1, input2, outputLeftXY, outputRightXY);
                    
                    if (fit == null) {
                        continue;
                    }
                    
                    fm = fit.getFundamentalMatrix();
                    
                    break;
                    
                case 3:
                    // scramble the points and use the matching
                    // method that solves for the epipolar projection
                    
                    if (true) {
                        continue;
                    }
                    
                    break;
            }
            
            Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
            int image1Width = img1.getWidth();
            int image1Height = img1.getHeight();
            Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
            
            for (int ii = 0; ii < input1.numCols(); ii++) {
                double x = input1.get(0, ii);
                double y = input1.get(1, ii);
                ImageIOHelper.addPointToImage((float) x, (float) y, img1, 3, 
                    255, 0, 0);
                double x2 = input2.get(0, ii);
                double y2 = input2.get(1, ii);
                ImageIOHelper.addPointToImage((float) x2, (float) y2, img2, 3, 
                    255, 0, 0);
            }

            Color clr = null;
            for (int ii = 0; ii < input2.numCols(); ii++) {
                clr = getColor(clr);
                SimpleMatrix epipolarLinesInLeft = fm.transpose().mult(input2);
                PairIntArray leftLine = spTransformer.getEpipolarLine(
                    epipolarLinesInLeft, image1Width, image1Height, ii);
                ImageIOHelper.addCurveToImage(leftLine, img1, 0,
                    clr.getRed(), clr.getGreen(), clr.getBlue());
            }

            clr = null;
            for (int ii = 0; ii < input1.numCols(); ii++) {
                clr = getColor(clr);
                SimpleMatrix epipolarLinesInRight = fm.mult(input1);
                PairIntArray rightLine = spTransformer.getEpipolarLine(
                    epipolarLinesInRight, img2.getWidth(), img2.getHeight(), ii);
                ImageIOHelper.addCurveToImage(rightLine, img2, 0,
                    clr.getRed(), clr.getGreen(), clr.getBlue());
            }

            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(
                dirPath + "/tmp_m_1_" + methodType + ".png", img1);
            ImageIOHelper.writeOutputImage(
                dirPath + "/tmp_m_2_" + methodType + ".png", img2);
        }
    
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
    
    public static void main(String[] args) {
        
        try {
            StereoProjectionTransformerTest test = 
                new StereoProjectionTransformerTest();
            
            //test.testB();
            
            //test.testC();
            
            //test.testE();
            
            //test.testF();
                        
            test.testH();
                        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }
    }

}
