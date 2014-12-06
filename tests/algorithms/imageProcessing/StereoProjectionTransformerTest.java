package algorithms.imageProcessing;

import Jama.Matrix;
import algorithms.util.ResourceFinder;
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
    
    @Test
    public void testF() throws Exception {
        
        /*
        test the fundamental matrix using the Merton I College set from
        http://www.robots.ox.ac.uk/~vgg/data/data-mview.html
        
        Note that this set was chosen to compare initial results with those from
        the "Programming Computer Vision with Python" by Jan Solem.
        */
        
        // mRows = 484;  nCols = 9
        Matrix A = readMerton1FundamentalMatrixAData();
        
        // mRows = 3;  nCols = 484
        Matrix unnormXY1 = readMerton1UnnormalizedX1Data();
        
        // mRows = 3;  nCols = 484
        Matrix unnormXY2 = readMerton1UnnormalizedX2Data();
        
        StereoProjectionTransformer spTransformer = 
            new StereoProjectionTransformer();
        
        StereoProjectionTransformer.NormalizedXY normXY1 = 
            spTransformer.normalize(unnormXY1);
        
        StereoProjectionTransformer.NormalizedXY normXY2 = 
            spTransformer.normalize(unnormXY2);

        
        Matrix fundMatrix = spTransformer.calculateFundamentalMatrix(
            normXY1, normXY2);
        fundMatrix = fundMatrix.times(1./fundMatrix.get(2, 2));
        
        double[][] leftRightEpipoles = 
            spTransformer.calculateEpipoles(fundMatrix);
                
       
        PairFloatArray matched1 = readMerton1UnnormalizedXY1Data();
        
        PairFloatArray matched2 = readMerton1UnnormalizedXY2Data();
       
        spTransformer = 
            new StereoProjectionTransformer();
        
        //spTransformer.calculateEpipolarProjection(matched1, matched2);
        spTransformer.calculateEpipolarProjectionWithoutNormalization(
            matched1, matched2);
        
        double[] leftEpipole = spTransformer.getLeftEpipole();
        
        double[] rightEpipole = spTransformer.getRightEpipole();
    
        int nPoints = spTransformer.getEpipolarLinesInLeft().getColumnDimension();
        
        String fileName1 = "merton_college_I_001.jpg";
        String fileName2 = "merton_college_I_002.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        
        int nLimitTo = 10;//nPoints;
        
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
        
        ImageIOHelper.writeOutputImage(dirPath + "/tmp1.png", img1);
       
        ImageIOHelper.writeOutputImage(dirPath + "/tmp2.png", img2);
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
    
    /**
     * dataset from 'Merton College I" at
     * http://www.robots.ox.ac.uk/~vgg/data/data-mview.html
     * 
     * @return 
     */
    private Matrix readMerton1FundamentalMatrixAData() {
        
        double[][] A = new double[484][];

        A[0]= new double[]{4.66704300e+05, 3.11380899e+05, 6.76587000e+02, 2.92965009e+05,
         1.95463611e+05, 4.24715000e+02, 6.89792000e+02, 4.60223000e+02,
         1.00000000e+00};
        A[1]= new double[]{4.65038845e+05, 2.89287075e+05, 6.75319000e+02, 2.71768409e+05,
         1.69059185e+05, 3.94656000e+02, 6.88621000e+02, 4.28371000e+02,
         1.00000000e+00};
        A[2]= new double[]{8.64275451e+05, 3.98196404e+05, 9.08129000e+02, 3.79547658e+05,
         1.74868455e+05, 3.98806000e+02, 9.51710000e+02, 4.38480000e+02,
         1.00000000e+00};
        A[3]= new double[]{8.62940029e+05, 3.84691605e+05, 9.06945000e+02, 3.66281741e+05,
         1.63285404e+05, 3.84960000e+02, 9.51480000e+02, 4.24162000e+02,
         1.00000000e+00};
        A[4]= new double[]{8.63234455e+05, 3.88322088e+05, 9.07284000e+02, 3.69698829e+05,
         1.66307335e+05, 3.88564000e+02, 9.51449000e+02, 4.28005000e+02,
         1.00000000e+00};
        A[5]= new double[]{8.62608147e+05, 3.73745451e+05, 9.06742000e+02, 3.55654550e+05,
         1.54095774e+05, 3.73851000e+02, 9.51327000e+02, 4.12185000e+02,
         1.00000000e+00};
        A[6]= new double[]{8.62069168e+05, 3.70813871e+05, 9.06447000e+02, 3.52380082e+05,
         1.51574174e+05, 3.70520000e+02, 9.51042000e+02, 4.09085000e+02,
         1.00000000e+00};
        A[7]= new double[]{5.29122826e+05, 1.86534570e+05, 7.17231000e+02, 1.72751283e+05,
         6.09009566e+04, 2.34166000e+02, 7.37730000e+02, 2.60076000e+02,
         1.00000000e+00};
        A[8]= new double[]{2.87722028e+05, 2.34346592e+05, 5.40106000e+02, 2.14730088e+05,
         1.74895418e+05, 4.03087000e+02, 5.32714000e+02, 4.33890000e+02,
         1.00000000e+00};
        A[9]= new double[]{4.41832506e+05, 2.82336691e+05, 6.59060000e+02, 2.64728773e+05,
         1.69165113e+05, 3.94883000e+02, 6.70398000e+02, 4.28393000e+02,
         1.00000000e+00};
        A[10]= new double[]{9.45635771e+05, 4.48745254e+05, 9.48040000e+02, 4.30358835e+05,
         2.04223963e+05, 4.31453000e+02, 9.97464000e+02, 4.73340000e+02,
         1.00000000e+00};
        A[11]= new double[]{5.06967976e+03, 4.29001708e+04, 8.47220000e+01, 2.72092720e+04,
         2.30247762e+05, 4.54708000e+02, 5.98390000e+01, 5.06364000e+02,
         1.00000000e+00};
        A[12]= new double[]{4.74406718e+04, 1.11476279e+05, 2.24362000e+02, 9.55615686e+04,
         2.24550953e+05, 4.51941000e+02, 2.11447000e+02, 4.96859000e+02,
         1.00000000e+00};
        A[13]= new double[]{5.12615473e+05, 1.21985821e+05, 7.08462000e+02, 1.10650566e+05,
         2.63312382e+04, 1.52925000e+02, 7.23561000e+02, 1.72184000e+02,
         1.00000000e+00};
        A[14]= new double[]{2.88145110e+05, 2.47983995e+05, 5.40541000e+02, 2.27205842e+05,
         1.95538326e+05, 4.26223000e+02, 5.33068000e+02, 4.58770000e+02,
         1.00000000e+00};
        A[15]= new double[]{7.40721422e+05, 2.52842364e+05, 8.42513000e+02, 2.36347591e+05,
         8.06763268e+04, 2.68827000e+02, 8.79181000e+02, 3.00105000e+02,
         1.00000000e+00};
        A[16]= new double[]{1.85637757e+04, 1.00095727e+04, 1.45319000e+02, 9.34595195e+03,
         5.03932968e+03, 7.31610000e+01, 1.27745000e+02, 6.88800000e+01,
         1.00000000e+00};
        A[17]= new double[]{3.61052540e+05, 2.75167586e+05, 5.99525000e+02, 2.56315536e+05,
         1.95344775e+05, 4.25610000e+02, 6.02231000e+02, 4.58976000e+02,
         1.00000000e+00};
        A[18]= new double[]{3.78836479e+04, 9.15960990e+04, 2.02223000e+02, 7.71327880e+04,
         1.86493721e+05, 4.11735000e+02, 1.87336000e+02, 4.52946000e+02,
         1.00000000e+00};
        A[19]= new double[]{3.08387352e+03, 3.59133351e+04, 7.00960000e+01, 2.01806825e+04,
         2.35014701e+05, 4.58704000e+02, 4.39950000e+01, 5.12345000e+02,
         1.00000000e+00};
        A[20]= new double[]{3.36254318e+05, 2.52048465e+05, 5.80710000e+02, 2.33014961e+05,
         1.74662629e+05, 4.02416000e+02, 5.79040000e+02, 4.34035000e+02,
         1.00000000e+00};
        A[21]= new double[]{3.33094220e+05, 1.13508297e+05, 5.79254000e+02, 1.03240381e+05,
         3.51811564e+04, 1.79536000e+02, 5.75040000e+02, 1.95956000e+02,
         1.00000000e+00};
        A[22]= new double[]{4.77814725e+04, 9.92963298e+04, 2.25174000e+02, 8.53320305e+04,
         1.77331443e+05, 4.02134000e+02, 2.12198000e+02, 4.40976000e+02,
         1.00000000e+00};
        A[23]= new double[]{8.59976638e+05, 3.60436679e+05, 9.05638000e+02, 3.41870051e+05,
         1.43285876e+05, 3.60022000e+02, 9.49581000e+02, 3.97992000e+02,
         1.00000000e+00};
        A[24]= new double[]{3.43151295e+05, 2.70718953e+05, 5.85170000e+02, 2.51653861e+05,
         1.98534788e+05, 4.29141000e+02, 5.86413000e+02, 4.62633000e+02,
         1.00000000e+00};
        A[25]= new double[]{4.63869316e+05, 2.70243701e+05, 6.74360000e+02, 2.53432534e+05,
         1.47646209e+05, 3.68433000e+02, 6.87866000e+02, 4.00741000e+02,
         1.00000000e+00};
        A[26]= new double[]{4.75426340e+04, 1.08451393e+05, 2.24807000e+02, 9.28526525e+04,
         2.11809878e+05, 4.39057000e+02, 2.11482000e+02, 4.82420000e+02,
         1.00000000e+00};
        A[27]= new double[]{1.82649973e+05, 4.60147151e+04, 4.32705000e+02, 4.21385967e+04,
         1.06159092e+04, 9.98280000e+01, 4.22112000e+02, 1.06342000e+02,
         1.00000000e+00};
        A[28]= new double[]{3.64534597e+05, 2.58202345e+05, 6.02217000e+02, 2.40348756e+05,
         1.70240666e+05, 3.97060000e+02, 6.05321000e+02, 4.28753000e+02,
         1.00000000e+00};
        A[29]= new double[]{3.68499708e+05, 1.92786097e+05, 6.05064000e+02, 1.78745477e+05,
         9.35133518e+04, 2.93494000e+02, 6.09026000e+02, 3.18621000e+02,
         1.00000000e+00};
        A[30]= new double[]{5.46636194e+03, 4.14075741e+04, 8.71300000e+01, 2.68430179e+04,
         2.03335283e+05, 4.27859000e+02, 6.27380000e+01, 4.75239000e+02,
         1.00000000e+00};
        A[31]= new double[]{2.63885888e+05, 1.20053197e+05, 5.17801000e+02, 1.09635762e+05,
         4.98780889e+04, 2.15129000e+02, 5.09628000e+02, 2.31852000e+02,
         1.00000000e+00};
        A[32]= new double[]{2.84904323e+05, 2.18078323e+05, 5.37234000e+02, 1.99858447e+05,
         1.52980462e+05, 3.76866000e+02, 5.30317000e+02, 4.05928000e+02,
         1.00000000e+00};
        A[33]= new double[]{7.22948383e+04, 1.41513699e+04, 2.75083000e+02, 1.40782596e+04,
         2.75575219e+03, 5.35680000e+01, 2.62811000e+02, 5.14440000e+01,
         1.00000000e+00};
        A[34]= new double[]{3.65411359e+03, 7.39332227e+03, 7.42540000e+01, 4.99949312e+03,
         1.01154118e+04, 1.01593000e+02, 4.92110000e+01, 9.95680000e+01,
         1.00000000e+00};
        A[35]= new double[]{1.86935768e+04, 3.17386361e+03, 1.45851000e+02, 4.08987279e+03,
         6.94393510e+02, 3.19100000e+01, 1.28169000e+02, 2.17610000e+01,
         1.00000000e+00};
        A[36]= new double[]{4.79441161e+05, 5.77472123e+04, 6.85191000e+02, 4.93777704e+04,
         5.94740047e+03, 7.05680000e+01, 6.99719000e+02, 8.42790000e+01,
         1.00000000e+00};
        A[37]= new double[]{2.02330509e+05, 6.62811079e+04, 4.55365000e+02, 6.03834591e+04,
         1.97809148e+04, 1.35899000e+02, 4.44326000e+02, 1.45556000e+02,
         1.00000000e+00};
        A[38]= new double[]{2.13343464e+03, 2.49131779e+03, 6.16190000e+01, 1.76743490e+03,
         2.06392169e+03, 5.10480000e+01, 3.46230000e+01, 4.04310000e+01,
         1.00000000e+00};
        A[39]= new double[]{1.46683142e+05, 1.00323288e+05, 3.86990000e+02, 9.11092624e+04,
         6.23137780e+04, 2.40371000e+02, 3.79036000e+02, 2.59240000e+02,
         1.00000000e+00};
        A[40]= new double[]{7.66073160e+02, 4.14415976e+03, 4.54050000e+01, 1.60260379e+03,
         8.66946721e+03, 9.49860000e+01, 1.68720000e+01, 9.12710000e+01,
         1.00000000e+00};
        A[41]= new double[]{1.98677788e+03, 5.79660178e+03, 6.02620000e+01, 3.26330459e+03,
         9.52098239e+03, 9.89810000e+01, 3.29690000e+01, 9.61900000e+01,
         1.00000000e+00};
        A[42]= new double[]{9.49250815e+05, 1.26322339e+05, 9.51115000e+02, 1.05881066e+05,
         1.40902105e+04, 1.06089000e+02, 9.98040000e+02, 1.32815000e+02,
         1.00000000e+00};
        A[43]= new double[]{6.91534237e+05, 1.61177410e+05, 8.12927000e+02, 1.46797915e+05,
         3.42145140e+04, 1.72567000e+02, 8.50672000e+02, 1.98268000e+02,
         1.00000000e+00};
        A[44]= new double[]{1.24695718e+05, 1.89947356e+05, 3.56008000e+02, 1.71299695e+05,
         2.60938586e+05, 4.89063000e+02, 3.50261000e+02, 5.33548000e+02,
         1.00000000e+00};
        A[45]= new double[]{2.19029391e+03, 2.51858939e+04, 6.13820000e+01, 1.32087761e+04,
         1.51885933e+05, 3.70170000e+02, 3.56830000e+01, 4.10314000e+02,
         1.00000000e+00};
        A[46]= new double[]{1.41554059e+05, 4.91524047e+04, 3.85047000e+02, 4.43885076e+04,
         1.54132062e+04, 1.20743000e+02, 3.67628000e+02, 1.27653000e+02,
         1.00000000e+00};
        A[47]= new double[]{1.82852387e+03, 1.37691324e+03, 5.87270000e+01, 1.14290915e+03,
         8.60632322e+02, 3.67070000e+01, 3.11360000e+01, 2.34460000e+01,
         1.00000000e+00};
        A[48]= new double[]{9.48663531e+05, 3.63152043e+05, 9.48143000e+02, 3.44527042e+05,
         1.31886275e+05, 3.44338000e+02, 1.00054900e+03, 3.83014000e+02,
         1.00000000e+00};
        A[49]= new double[]{1.35248353e+05, 9.55156634e+04, 3.71827000e+02, 8.66261360e+04,
         6.11774758e+04, 2.38154000e+02, 3.63740000e+02, 2.56882000e+02,
         1.00000000e+00};
        A[50]= new double[]{ 856.387527, 940.184637, 47.343,   612.421184,  672.346304, 33.856,
         18.089, 19.859, 1.};
        A[51]= new double[]{4.67608973e+05, 1.38990514e+05, 6.76781000e+02, 1.27443605e+05,
         3.78809072e+04, 1.84452000e+02, 6.90931000e+02, 2.05370000e+02,
         1.00000000e+00};
        A[52]= new double[]{5.51121129e+03, 3.87053492e+04, 8.73160000e+01, 2.52264973e+04,
         1.77166204e+05, 3.99672000e+02, 6.31180000e+01, 4.43279000e+02,
         1.00000000e+00};
        A[53]= new double[]{1.42862337e+05, 1.10984415e+05, 3.81895000e+02, 1.00441132e+05,
         7.80289650e+04, 2.68496000e+02, 3.74088000e+02, 2.90615000e+02,
         1.00000000e+00};
        A[54]= new double[]{3.34953527e+05, 4.58966009e+04, 5.86808000e+02, 3.98142893e+04,
         5.45550471e+03, 6.97510000e+01, 5.70806000e+02, 7.82140000e+01,
         1.00000000e+00};
        A[55]= new double[]{4.89685484e+03, 2.46785760e+04, 8.31060000e+01, 1.60436723e+04,
         8.08549567e+04, 2.72282000e+02, 5.89230000e+01, 2.96953000e+02,
         1.00000000e+00};
        A[56]= new double[]{3.31772137e+05, 2.32179200e+05, 5.76635000e+02, 2.14590496e+05,
         1.50173700e+05, 3.72968000e+02, 5.75359000e+02, 4.02645000e+02,
         1.00000000e+00};
        A[57]= new double[]{8.04020963e+04, 4.08831755e+04, 2.90638000e+02, 3.68780485e+04,
         1.87518958e+04, 1.33307000e+02, 2.76640000e+02, 1.40667000e+02,
         1.00000000e+00};
        A[58]= new double[]{6.00723076e+05, 2.39408419e+05, 7.57883000e+02, 2.25629325e+05,
         8.99209003e+04, 2.84658000e+02, 7.92633000e+02, 3.15891000e+02,
         1.00000000e+00};
        A[59]= new double[]{2.61258710e+05, 1.09431501e+05, 5.17551000e+02, 9.90030030e+04,
         4.14686547e+04, 1.96124000e+02, 5.04798000e+02, 2.11441000e+02,
         1.00000000e+00};
        A[60]= new double[]{1.84423823e+04, 2.08966283e+04, 1.45050000e+02, 1.76141597e+04,
         1.99581888e+04, 1.38536000e+02, 1.27145000e+02, 1.44065000e+02,
         1.00000000e+00};
        A[61]= new double[]{1.30468772e+05, 1.04645998e+05, 3.65203000e+02, 9.46816102e+04,
         7.59419397e+04, 2.65029000e+02, 3.57250000e+02, 2.86542000e+02,
         1.00000000e+00};
        A[62]= new double[]{8.03559613e+05, 4.48725279e+04, 8.84362000e+02, 2.81185259e+04,
         1.57020004e+03, 3.09460000e+01, 9.08632000e+02, 5.07400000e+01,
         1.00000000e+00};
        A[63]= new double[]{2.57432794e+05, 1.06619934e+05, 5.14230000e+02, 9.64135200e+04,
         3.99312107e+04, 1.92589000e+02, 5.00618000e+02, 2.07339000e+02,
         1.00000000e+00};
        A[64]= new double[]{1.61228690e+05, 1.17955878e+05, 4.05142000e+02, 1.07115031e+05,
         7.83660000e+04, 2.69163000e+02, 3.97956000e+02, 2.91147000e+02,
         1.00000000e+00};
        A[65]= new double[]{9.64117542e+05, 1.07482769e+05, 9.59959000e+02, 8.65352538e+04,
         9.64721449e+03, 8.61620000e+01, 1.00433200e+03, 1.11966000e+02,
         1.00000000e+00};
        A[66]= new double[]{7.08605887e+05, 2.28478328e+04, 8.34807000e+02, 9.12657715e+03,
         2.94271488e+02, 1.07520000e+01, 8.48826000e+02, 2.73690000e+01,
         1.00000000e+00};
        A[67]= new double[]{5.61394005e+05, 2.89004712e+05, 7.36695000e+02, 2.72705828e+05,
         1.40388512e+05, 3.57861000e+02, 7.62044000e+02, 3.92299000e+02,
         1.00000000e+00};
        A[68]= new double[]{5.41655955e+04, 3.45161755e+04, 2.40613000e+02, 3.09001854e+04,
         1.96906581e+04, 1.37264000e+02, 2.25115000e+02, 1.43451000e+02,
         1.00000000e+00};
        A[69]= new double[]{3.78300167e+04, 8.60196902e+04, 2.01982000e+02, 7.26912362e+04,
         1.65288788e+05, 3.88113000e+02, 1.87294000e+02, 4.25878000e+02,
         1.00000000e+00};
        A[70]= new double[]{6.08708202e+05, 4.03493468e+05, 7.66748000e+02, 3.85257924e+05,
         2.55375326e+05, 4.85283000e+02, 7.93883000e+02, 5.26240000e+02,
         1.00000000e+00};
        A[71]= new double[]{5.59771537e+05, 2.46424081e+05, 7.35915000e+02, 2.31347742e+05,
         1.01844505e+05, 3.04146000e+02, 7.60647000e+02, 3.34854000e+02,
         1.00000000e+00};
        A[72]= new double[]{1.80434037e+05, 3.00410013e+04, 4.36262000e+02, 2.72461343e+04,
         4.53629022e+03, 6.58770000e+01, 4.13591000e+02, 6.88600000e+01,
         1.00000000e+00};
        A[73]= new double[]{7.53392233e+05, 1.90079854e+05, 8.47984000e+02, 1.74250118e+05,
         4.39630718e+04, 1.96128000e+02, 8.88451000e+02, 2.24155000e+02,
         1.00000000e+00};
        A[74]= new double[]{1.30753836e+05, 9.86966518e+04, 3.65374000e+02, 8.95369647e+04,
         6.75850049e+04, 2.50199000e+02, 3.57863000e+02, 2.70125000e+02,
         1.00000000e+00};
        A[75]= new double[]{9.52702198e+05, 1.44951890e+05, 9.51490000e+02, 1.24369245e+05,
         1.89225522e+04, 1.24211000e+02, 1.00127400e+03, 1.52342000e+02,
         1.00000000e+00};
        A[76]= new double[]{1.00217563e+05, 6.11062348e+04, 3.19712000e+02, 5.58673919e+04,
         3.40643483e+04, 1.78227000e+02, 3.13462000e+02, 1.91129000e+02,
         1.00000000e+00};
        A[77]= new double[]{6.48670996e+05, 2.51513256e+04, 8.00284000e+02, 1.32589933e+04,
         5.14099224e+02, 1.63580000e+01, 8.10551000e+02, 3.14280000e+01,
         1.00000000e+00};
        A[78]= new double[]{8.62672726e+04, 4.64952234e+04, 3.00316000e+02, 4.20840065e+04,
         2.26818958e+04, 1.46504000e+02, 2.87255000e+02, 1.54821000e+02,
         1.00000000e+00};
        A[79]= new double[]{4.88915012e+05, 2.80862767e+05, 6.90206000e+02, 2.64561500e+05,
         1.51980350e+05, 3.73484000e+02, 7.08361000e+02, 4.06926000e+02,
         1.00000000e+00};
        A[80]= new double[]{5.58349119e+05, 2.36394870e+05, 7.34901000e+02, 2.21609368e+05,
         9.38253789e+04, 2.91683000e+02, 7.59761000e+02, 3.21669000e+02,
         1.00000000e+00};
        A[81]= new double[]{2.03687814e+04, 2.01622178e+04, 1.52558000e+02, 1.72613533e+04,
         1.70863027e+04, 1.29284000e+02, 1.33515000e+02, 1.32161000e+02,
         1.00000000e+00};
        A[82]= new double[]{9.07461024e+05, 7.30909209e+04, 9.34547000e+02, 5.39681538e+04,
         4.34683359e+03, 5.55790000e+01, 9.71017000e+02, 7.82100000e+01,
         1.00000000e+00};
        A[83]= new double[]{5.54780933e+03, 6.80636512e+03, 8.70190000e+01, 5.25511471e+03,
         6.44727088e+03, 8.24280000e+01, 6.37540000e+01, 7.82170000e+01,
         1.00000000e+00};
        A[84]= new double[]{1.09966256e+05, 3.75445193e+04, 3.41699000e+02, 3.39229352e+04,
         1.15819193e+04, 1.05409000e+02, 3.21822000e+02, 1.09876000e+02,
         1.00000000e+00};
        A[85]= new double[]{5.98491683e+05, 1.47285437e+05, 7.60281000e+02, 1.33981887e+05,
         3.29721887e+04, 1.70201000e+02, 7.87198000e+02, 1.93725000e+02,
         1.00000000e+00};
        A[86]= new double[]{1.09057600e+03, 2.11885535e+04, 4.96190000e+01, 8.45518943e+03,
         1.64273955e+05, 3.84694000e+02, 2.19790000e+01, 4.27025000e+02,
         1.00000000e+00};
        A[87]= new double[]{5.55395763e+05, 1.36001555e+05, 7.34667000e+02, 1.23505699e+05,
         3.02432395e+04, 1.63371000e+02, 7.55983000e+02, 1.85120000e+02,
         1.00000000e+00};
        A[88]= new double[]{6.20736432e+04, 2.84262382e+04, 2.57137000e+02, 2.59599958e+04,
         1.18882184e+04, 1.07538000e+02, 2.41403000e+02, 1.10549000e+02,
         1.00000000e+00};
        A[89]= new double[]{2.91708064e+05, 7.85111437e+04, 5.47383000e+02, 7.03286606e+04,
         1.89284571e+04, 1.31970000e+02, 5.32914000e+02, 1.43430000e+02,
         1.00000000e+00};
        A[90]= new double[]{3.34475461e+05, 4.99017181e+04, 5.86727000e+02, 4.29855583e+04,
         6.41318560e+03, 7.54040000e+01, 5.70070000e+02, 8.50510000e+01,
         1.00000000e+00};
        A[91]= new double[]{9.22989873e+05, 3.61875425e+04, 9.45511000e+02, 1.61655574e+04,
         6.33800880e+02, 1.65600000e+01, 9.76181000e+02, 3.82730000e+01,
         1.00000000e+00};
        A[92]= new double[]{6.05869480e+04, 1.01005319e+05, 2.50669000e+02, 8.89882657e+04,
         1.48353539e+05, 3.68175000e+02, 2.41701000e+02, 4.02943000e+02,
         1.00000000e+00};
        A[93]= new double[]{1.62499228e+03, 2.35099529e+04, 5.54870000e+01, 1.11634425e+04,
         1.61509694e+05, 3.81187000e+02, 2.92860000e+01, 4.23702000e+02,
         1.00000000e+00};
        A[94]= new double[]{6.79250346e+05, 1.73147445e+05, 8.06008000e+02, 1.58514052e+05,
         4.04067560e+04, 1.88095000e+02, 8.42734000e+02, 2.14821000e+02,
         1.00000000e+00};
        A[95]= new double[]{8.10756764e+05, 3.64153655e+05, 8.78711000e+02, 3.47026677e+05,
         1.55867997e+05, 3.76113000e+02, 9.22666000e+02, 4.14418000e+02,
         1.00000000e+00};
        A[96]= new double[]{3.82579152e+03, 5.56934525e+03, 7.56160000e+01, 4.01010910e+03,
         5.83766313e+03, 7.92590000e+01, 5.05950000e+01, 7.36530000e+01,
         1.00000000e+00};
        A[97]= new double[]{3.14122213e+04, 4.16932132e+04, 1.84880000e+02, 3.57280036e+04,
         4.74215197e+04, 2.10281000e+02, 1.69906000e+02, 2.25515000e+02,
         1.00000000e+00};
        A[98]= new double[]{1.28401268e+05, 8.28468094e+04, 3.61185000e+02, 7.57513620e+04,
         4.88761425e+04, 2.13084000e+02, 3.55500000e+02, 2.29375000e+02,
         1.00000000e+00};
        A[99]= new double[]{8.80925528e+05, 1.18674619e+05, 9.17909000e+02, 9.97531132e+04,
         1.34383240e+04, 1.03941000e+02, 9.59709000e+02, 1.29288000e+02,
         1.00000000e+00};
        A[100]= new double[]{5.52313944e+05, 2.70511431e+05, 7.31193000e+02, 2.54819185e+05,
         1.24804929e+05, 3.37348000e+02, 7.55360000e+02, 3.69959000e+02,
         1.00000000e+00};
        A[101]= new double[]{7.55815086e+05, 1.60367470e+05, 8.49323000e+02, 1.44522917e+05,
         3.06646097e+04, 1.62403000e+02, 8.89903000e+02, 1.88818000e+02,
         1.00000000e+00};
        A[102]= new double[]{8.20012424e+02, 3.06629817e+03, 4.67830000e+01, 1.28555610e+03,
         4.80712025e+03, 7.33430000e+01, 1.75280000e+01, 6.55430000e+01,
         1.00000000e+00};
        A[103]= new double[]{1.55219613e+05, 6.16392726e+04, 4.00726000e+02, 5.57208841e+04,
         2.21273246e+04, 1.43853000e+02, 3.87346000e+02, 1.53819000e+02,
         1.00000000e+00};
        A[104]= new double[]{1.29135075e+03, 2.65607240e+04, 5.24790000e+01, 1.11347905e+04,
         2.29022283e+05, 4.52505000e+02, 2.46070000e+01, 5.06121000e+02,
         1.00000000e+00};
        A[105]= new double[]{1.70424198e+05, 1.04460220e+05, 4.16626000e+02, 9.51108937e+04,
         5.82975012e+04, 2.32512000e+02, 4.09058000e+02, 2.50729000e+02,
         1.00000000e+00};
        A[106]= new double[]{9.67066837e+05, 3.54883472e+05, 9.56683000e+02, 3.35834003e+05,
         1.23240641e+05, 3.32228000e+02, 1.01085400e+03, 3.70952000e+02,
         1.00000000e+00};
        A[107]= new double[]{2.17197240e+04, 9.79769758e+03, 1.57032000e+02, 9.22374572e+03,
         4.16080199e+03, 6.66870000e+01, 1.38314000e+02, 6.23930000e+01,
         1.00000000e+00};
        A[108]= new double[]{2.58847531e+05, 1.29728335e+05, 5.13174000e+02, 1.17903660e+05,
         5.90905594e+04, 2.33748000e+02, 5.04405000e+02, 2.52796000e+02,
         1.00000000e+00};
        A[109]= new double[]{3.29506501e+04, 7.85043504e+04, 1.88643000e+02, 6.60939634e+04,
         1.57467718e+05, 3.78389000e+02, 1.74672000e+02, 4.16153000e+02,
         1.00000000e+00};
        A[110]= new double[]{8.44537884e+05, 1.35652904e+05, 8.98435000e+02, 1.17583031e+05,
         1.88866360e+04, 1.25087000e+02, 9.40010000e+02, 1.50988000e+02,
         1.00000000e+00};
        A[111]= new double[]{7.66404284e+05, 3.91861889e+04, 8.64980000e+02, 2.37962957e+04,
         1.21670267e+03, 2.68570000e+01, 8.86037000e+02, 4.53030000e+01,
         1.00000000e+00};
        A[112]= new double[]{2.54642371e+05, 2.62544132e+05, 5.08839000e+02, 2.39709802e+05,
         2.47148193e+05, 4.79000000e+02, 5.00438000e+02, 5.15967000e+02,
         1.00000000e+00};
        A[113]= new double[]{5.66704284e+05, 3.74557111e+05, 7.40249000e+02, 3.56561401e+05,
         2.35665429e+05, 4.65753000e+02, 7.65559000e+02, 5.05988000e+02,
         1.00000000e+00};
        A[114]= new double[]{3.55770586e+04, 2.43303608e+04, 1.98505000e+02, 2.14071717e+04,
         1.46398896e+04, 1.19443000e+02, 1.79225000e+02, 1.22568000e+02,
         1.00000000e+00};
        A[115]= new double[]{4.11336527e+03, 3.47740884e+04, 7.84260000e+01, 2.09562077e+04,
         1.77162244e+05, 3.99554000e+02, 5.24490000e+01, 4.43400000e+02,
         1.00000000e+00};
        A[116]= new double[]{1.94695962e+04, 2.53376243e+04, 1.48033000e+02, 2.14347979e+04,
         2.78951270e+04, 1.62975000e+02, 1.31522000e+02, 1.71162000e+02,
         1.00000000e+00};
        A[117]= new double[]{6.83548243e+05, 2.78383911e+04, 8.20272000e+02, 1.49997420e+04,
         6.10884000e+02, 1.80000000e+01, 8.33319000e+02, 3.39380000e+01,
         1.00000000e+00};
        A[118]= new double[]{1.18050575e+05, 7.23991463e+04, 3.46482000e+02, 6.61114158e+04,
         4.05454192e+04, 1.94039000e+02, 3.40712000e+02, 2.08955000e+02,
         1.00000000e+00};
        A[119]= new double[]{3.28461224e+04, 7.96112320e+04, 1.88709000e+02, 6.67992473e+04,
         1.61905576e+05, 3.83778000e+02, 1.74057000e+02, 4.21873000e+02,
         1.00000000e+00};
        A[120]= new double[]{7.56597234e+04, 8.96095085e+04, 2.79510000e+02, 7.98391307e+04,
         9.45594953e+04, 2.94950000e+02, 2.70687000e+02, 3.20595000e+02,
         1.00000000e+00};
        A[121]= new double[]{3.97219749e+04, 2.96160748e+04, 2.07889000e+02, 2.61450918e+04,
         1.94933660e+04, 1.36833000e+02, 1.91073000e+02, 1.42461000e+02,
         1.00000000e+00};
        A[122]= new double[]{9.12744326e+05, 1.69408733e+05, 9.30209000e+02, 1.49684893e+05,
         2.77820713e+04, 1.52549000e+02, 9.81225000e+02, 1.82119000e+02,
         1.00000000e+00};
        A[123]= new double[]{3.07951949e+05, 2.13988639e+05, 5.56144000e+02, 1.97269120e+05,
         1.37077718e+05, 3.56257000e+02, 5.53727000e+02, 3.84772000e+02,
         1.00000000e+00};
        A[124]= new double[]{1.30257444e+05, 1.80960276e+05, 3.64166000e+02, 1.63242624e+05,
         2.26784968e+05, 4.56384000e+02, 3.57687000e+02, 4.96917000e+02,
         1.00000000e+00};
        A[125]= new double[]{3.07677802e+05, 2.31227866e+05, 5.56133000e+02, 2.13257690e+05,
         1.60268698e+05, 3.85467000e+02, 5.53245000e+02, 4.15778000e+02,
         1.00000000e+00};
        A[126]= new double[]{1.27546269e+05, 2.48946650e+04, 3.70815000e+02, 2.27551501e+04,
         4.44138306e+03, 6.61560000e+01, 3.43962000e+02, 6.71350000e+01,
         1.00000000e+00};
        A[127]= new double[]{9.98331311e+04, 1.97410630e+04, 3.29578000e+02, 1.83525293e+04,
         3.62904013e+03, 6.05870000e+01, 3.02912000e+02, 5.98980000e+01,
         1.00000000e+00};
        A[128]= new double[]{1.92889121e+05, 7.76157504e+04, 4.44973000e+02, 7.00953915e+04,
         2.82053565e+04, 1.61702000e+02, 4.33485000e+02, 1.74428000e+02,
         1.00000000e+00};
        A[129]= new double[]{1.04747667e+03, 2.70149255e+04, 4.91750000e+01, 1.04363823e+04,
         2.69159303e+05, 4.89948000e+02, 2.13010000e+01, 5.49363000e+02,
         1.00000000e+00};
        A[130]= new double[]{5.48802086e+03, 1.76751826e+04, 8.70120000e+01, 1.21001109e+04,
         3.89706372e+04, 1.91846000e+02, 6.30720000e+01, 2.03135000e+02,
         1.00000000e+00};
        A[131]= new double[]{4.42282186e+04, 1.92682427e+04, 2.21296000e+02, 1.75557024e+04,
         7.64822880e+03, 8.78400000e+01, 1.99860000e+02, 8.70700000e+01,
         1.00000000e+00};
        A[132]= new double[]{4.29531927e+05, 6.50892809e+04, 6.57421000e+02, 5.65462614e+04,
         8.56875883e+03, 8.65470000e+01, 6.53359000e+02, 9.90070000e+01,
         1.00000000e+00};
        A[133]= new double[]{2.58782048e+04, 4.27280976e+04, 1.69247000e+02, 3.57385490e+04,
         5.90087381e+04, 2.33735000e+02, 1.52902000e+02, 2.52460000e+02,
         1.00000000e+00};
        A[134]= new double[]{7.02378397e+03, 1.86465004e+04, 9.58970000e+01, 1.34368678e+04,
         3.56717350e+04, 1.83456000e+02, 7.32430000e+01, 1.94443000e+02,
         1.00000000e+00};
        A[135]= new double[]{4.47582402e+05, 3.36536925e+05, 6.62606000e+02, 3.16820759e+05,
         2.38217328e+05, 4.69025000e+02, 6.75488000e+02, 5.07899000e+02,
         1.00000000e+00};
        A[136]= new double[]{7.47028217e+03, 3.89602038e+04, 9.76980000e+01, 2.76254702e+04,
         1.44076746e+05, 3.61292000e+02, 7.64630000e+01, 3.98782000e+02,
         1.00000000e+00};
        A[137]= new double[]{2.39246729e+05, 2.23568343e+05, 4.92969000e+02, 2.04225212e+05,
         1.90841866e+05, 4.20807000e+02, 4.85318000e+02, 4.53514000e+02,
         1.00000000e+00};
        A[138]= new double[]{3.74936210e+05, 1.46407670e+05, 6.08028000e+02, 1.35491883e+05,
         5.29078025e+04, 2.19725000e+02, 6.16643000e+02, 2.40791000e+02,
         1.00000000e+00};
        A[139]= new double[]{2.75908174e+05, 2.09930742e+05, 5.29102000e+02, 1.92138472e+05,
         1.46192741e+05, 3.68459000e+02, 5.21465000e+02, 3.96768000e+02,
         1.00000000e+00};
        A[140]= new double[]{2.38190700e+05, 1.98871785e+05, 4.92157000e+02, 1.81299190e+05,
         1.51371542e+05, 3.74606000e+02, 4.83973000e+02, 4.04082000e+02,
         1.00000000e+00};
        A[141]= new double[]{3.97655472e+03, 4.34636441e+04, 7.66890000e+01, 2.62430107e+04,
         2.86835454e+05, 5.06104000e+02, 5.18530000e+01, 5.66752000e+02,
         1.00000000e+00};
        A[142]= new double[]{1.26918966e+05, 1.56173440e+05, 3.59268000e+02, 1.41321471e+05,
         1.73895684e+05, 4.00037000e+02, 3.53271000e+02, 4.34699000e+02,
         1.00000000e+00};
        A[143]= new double[]{6.20552286e+04, 6.06574556e+04, 2.54002000e+02, 5.40648258e+04,
         5.28470339e+04, 2.21296000e+02, 2.44310000e+02, 2.38807000e+02,
         1.00000000e+00};
        A[144]= new double[]{1.06781208e+05, 1.55358502e+05, 3.30505000e+02, 1.39090677e+05,
         2.02366313e+05, 4.30508000e+02, 3.23085000e+02, 4.70064000e+02,
         1.00000000e+00};
        A[145]= new double[]{8.04466266e+03, 2.73277085e+04, 1.01533000e+02, 1.97057115e+04,
         6.69402761e+04, 2.48709000e+02, 7.92320000e+01, 2.69151000e+02,
         1.00000000e+00};
        A[146]= new double[]{1.81185490e+03, 3.08326436e+04, 5.82590000e+01, 1.47107665e+04,
         2.50335621e+05, 4.73015000e+02, 3.11000000e+01, 5.29234000e+02,
         1.00000000e+00};
        A[147]= new double[]{8.60520129e+05, 9.04326604e+04, 9.09584000e+02, 7.23281027e+04,
         7.60101074e+03, 7.64520000e+01, 9.46059000e+02, 9.94220000e+01,
         1.00000000e+00};
        A[148]= new double[]{1.14848035e+04, 2.53748267e+04, 1.17494000e+02, 1.97282833e+04,
         4.35881877e+04, 2.01828000e+02, 9.77480000e+01, 2.15967000e+02,
         1.00000000e+00};
        A[149]= new double[]{5.20038270e+02, 1.25143128e+04, 4.11260000e+01, 3.52253030e+03,
         8.47669267e+04, 2.78571000e+02, 1.26450000e+01, 3.04292000e+02,
         1.00000000e+00};
        A[150]= new double[]{1.04579558e+05, 1.61758227e+05, 3.27038000e+02, 1.44822340e+05,
         2.24003673e+05, 4.52884000e+02, 3.19778000e+02, 4.94616000e+02,
         1.00000000e+00};
        A[151]= new double[]{5.54294190e+05, 4.77492667e+04, 7.41886000e+02, 3.77007853e+04,
         3.24770652e+03, 5.04600000e+01, 7.47142000e+02, 6.43620000e+01,
         1.00000000e+00};
        A[152]= new double[]{4.30450763e+05, 1.52907424e+05, 6.49636000e+02, 1.41222565e+05,
         5.01659667e+04, 2.13133000e+02, 6.62603000e+02, 2.35374000e+02,
         1.00000000e+00};
        A[153]= new double[]{4.33388703e+05, 3.30934993e+05, 6.51878000e+02, 3.11777816e+05,
         2.38073094e+05, 4.68958000e+02, 6.64831000e+02, 5.07664000e+02,
         1.00000000e+00};
        A[154]= new double[]{3.01119020e+05, 6.08080975e+04, 5.57213000e+02, 5.42098862e+04,
         1.09471665e+04, 1.00314000e+02, 5.40402000e+02, 1.09129000e+02,
         1.00000000e+00};
        A[155]= new double[]{1.17949230e+05, 2.67967444e+04, 3.56425000e+02, 2.45101429e+04,
         5.56843001e+03, 7.40660000e+01, 3.30923000e+02, 7.51820000e+01,
         1.00000000e+00};
        A[156]= new double[]{9.66254925e+02, 1.87010259e+04, 4.79650000e+01, 7.09831234e+03,
         1.37381678e+05, 3.52361000e+02, 2.01450000e+01, 3.89889000e+02,
         1.00000000e+00};
        A[157]= new double[]{1.54041200e+05, 8.11689399e+04, 3.96004000e+02, 7.40891789e+04,
         3.90398160e+04, 1.90466000e+02, 3.88989000e+02, 2.04970000e+02,
         1.00000000e+00};
        A[158]= new double[]{1.71087771e+05, 1.68053549e+05, 4.16846000e+02, 1.52691709e+05,
         1.49983739e+05, 3.72025000e+02, 4.10434000e+02, 4.03155000e+02,
         1.00000000e+00};
        A[159]= new double[]{3.63073860e+03, 2.51990838e+04, 7.38000000e+01, 1.53010542e+04,
         1.06196724e+05, 3.11016000e+02, 4.91970000e+01, 3.41451000e+02,
         1.00000000e+00};
        A[160]= new double[]{6.37600388e+05, 1.17273230e+05, 7.85982000e+02, 1.03963692e+05,
         1.91219425e+04, 1.28158000e+02, 8.11215000e+02, 1.49206000e+02,
         1.00000000e+00};
        A[161]= new double[]{7.36046464e+05, 4.00588422e+05, 8.39136000e+02, 3.83349639e+05,
         2.08635507e+05, 4.37041000e+02, 8.77148000e+02, 4.77382000e+02,
         1.00000000e+00};
        A[162]= new double[]{8.86361007e+05, 5.01792753e+04, 9.26039000e+02, 3.11419300e+04,
         1.76302823e+03, 3.25360000e+01, 9.57153000e+02, 5.41870000e+01,
         1.00000000e+00};
        A[163]= new double[]{6.23971637e+03, 2.10622289e+04, 9.20190000e+01, 1.44859689e+04,
         4.88975418e+04, 2.13629000e+02, 6.78090000e+01, 2.28890000e+02,
         1.00000000e+00};
        A[164]= new double[]{6.83415681e+05, 3.31772098e+05, 8.09528000e+02, 3.14975772e+05,
         1.52908656e+05, 3.73099000e+02, 8.44215000e+02, 4.09834000e+02,
         1.00000000e+00};
        A[165]= new double[]{1.79257248e+05, 1.76579423e+05, 4.26677000e+02, 1.60646175e+05,
         1.58246371e+05, 3.82378000e+02, 4.20124000e+02, 4.13848000e+02,
         1.00000000e+00};
        A[166]= new double[]{1.28643261e+05, 1.53753529e+05, 3.61882000e+02, 1.38615298e+05,
         1.65672038e+05, 3.89934000e+02, 3.55484000e+02, 4.24872000e+02,
         1.00000000e+00};
        A[167]= new double[]{4.01514085e+05, 4.88963333e+04, 6.38325000e+02, 4.15324043e+04,
         5.05781083e+03, 6.60280000e+01, 6.29012000e+02, 7.66010000e+01,
         1.00000000e+00};
        A[168]= new double[]{2.66757652e+04, 6.56562413e+04, 1.70970000e+02, 5.46345322e+04,
         1.34470296e+05, 3.50163000e+02, 1.56026000e+02, 3.84022000e+02,
         1.00000000e+00};
        A[169]= new double[]{8.08299535e+05, 1.79302461e+05, 8.76309000e+02, 1.62453348e+05,
         3.60364985e+04, 1.76122000e+02, 9.22391000e+02, 2.04611000e+02,
         1.00000000e+00};
        A[170]= new double[]{1.45940471e+05, 1.73303587e+05, 3.85179000e+02, 1.56929418e+05,
         1.86352907e+05, 4.14182000e+02, 3.78890000e+02, 4.49930000e+02,
         1.00000000e+00};
        A[171]= new double[]{1.01636956e+05, 3.61850244e+04, 3.28886000e+02, 3.27708925e+04,
         1.16671690e+04, 1.06043000e+02, 3.09034000e+02, 1.10023000e+02,
         1.00000000e+00};
        A[172]= new double[]{9.59329526e+04, 2.83034585e+04, 3.21021000e+02, 2.57851505e+04,
         7.60748960e+03, 8.62850000e+01, 2.98837000e+02, 8.81670000e+01,
         1.00000000e+00};
        A[173]= new double[]{1.27744909e+05, 1.28885008e+05, 3.60791000e+02, 1.16701496e+05,
         1.17743036e+05, 3.29601000e+02, 3.54069000e+02, 3.57229000e+02,
         1.00000000e+00};
        A[174]= new double[]{5.24088906e+03, 2.05759790e+04, 8.51070000e+01, 1.38564853e+04,
         5.44012183e+04, 2.25016000e+02, 6.15800000e+01, 2.41766000e+02,
         1.00000000e+00};
        A[175]= new double[]{6.65566473e+05, 2.74942497e+05, 7.99126000e+02, 2.58961149e+05,
         1.06975678e+05, 3.10927000e+02, 8.32868000e+02, 3.44054000e+02,
         1.00000000e+00};
        A[176]= new double[]{1.26486912e+05, 3.75928104e+04, 3.66684000e+02, 3.39694442e+04,
         1.00959605e+04, 9.84770000e+01, 3.44948000e+02, 1.02521000e+02,
         1.00000000e+00};
        A[177]= new double[]{3.44071358e+05, 2.49089427e+04, 6.07313000e+02, 1.98654040e+04,
         1.43814996e+03, 3.50640000e+01, 5.66547000e+02, 4.10150000e+01,
         1.00000000e+00};
        A[178]= new double[]{3.17509633e+05, 1.27083605e+05, 5.64811000e+02, 1.16286201e+05,
         4.65436887e+04, 2.06859000e+02, 5.62152000e+02, 2.25002000e+02,
         1.00000000e+00};
        A[179]= new double[]{1.20284656e+05, 6.39043131e+04, 3.51834000e+02, 5.81833614e+04,
         3.09114052e+04, 1.70187000e+02, 3.41879000e+02, 1.81632000e+02,
         1.00000000e+00};
        A[180]= new double[]{3.02931022e+05, 6.89578984e+04, 5.57808000e+02, 6.13727927e+04,
         1.39706352e+04, 1.13010000e+02, 5.43074000e+02, 1.23623000e+02,
         1.00000000e+00};
        A[181]= new double[]{1.51726585e+05, 5.36689375e+04, 3.96303000e+02, 4.86581905e+04,
         1.72114424e+04, 1.27093000e+02, 3.82855000e+02, 1.35424000e+02,
         1.00000000e+00};
        A[182]= new double[]{2.24750662e+05, 9.61185140e+04, 4.80374000e+02, 8.69196776e+04,
         3.71727059e+04, 1.85779000e+02, 4.67866000e+02, 2.00091000e+02,
         1.00000000e+00};
        A[183]= new double[]{4.07217578e+05, 2.66247732e+05, 6.33430000e+02, 2.49391275e+05,
         1.63057453e+05, 3.87930000e+02, 6.42877000e+02, 4.20327000e+02,
         1.00000000e+00};
        A[184]= new double[]{1.72763825e+05, 2.18109861e+05, 4.19052000e+02, 1.97778077e+05,
         2.49689707e+05, 4.79726000e+02, 4.12273000e+02, 5.20484000e+02,
         1.00000000e+00};
        A[185]= new double[]{5.44277176e+05, 2.36854663e+05, 7.25820000e+02, 2.22566337e+05,
         9.68548326e+04, 2.96803000e+02, 7.49879000e+02, 3.26327000e+02,
         1.00000000e+00};
        A[186]= new double[]{8.43249584e+03, 5.78377073e+04, 1.02814000e+02, 4.13072059e+04,
         2.83322296e+05, 5.03642000e+02, 8.20170000e+01, 5.62547000e+02,
         1.00000000e+00};
        A[187]= new double[]{3.09080104e+04, 4.62962329e+04, 1.82984000e+02, 3.96040554e+04,
         5.93217923e+04, 2.34467000e+02, 1.68911000e+02, 2.53007000e+02,
         1.00000000e+00};
        A[188]= new double[]{6.80871520e+05, 3.18757445e+05, 8.08018000e+02, 3.02371845e+05,
         1.41558685e+05, 3.58837000e+02, 8.42644000e+02, 3.94493000e+02,
         1.00000000e+00};
        A[189]= new double[]{3.16851574e+03, 2.66801018e+04, 7.04600000e+01, 1.54688863e+04,
         1.30253877e+05, 3.43990000e+02, 4.49690000e+01, 3.78656000e+02,
         1.00000000e+00};
        A[190]= new double[]{6.29957749e+05, 2.87341623e+05, 7.78558000e+02, 2.71110865e+05,
         1.23661366e+05, 3.35063000e+02, 8.09134000e+02, 3.69069000e+02,
         1.00000000e+00};
        A[191]= new double[]{1.31068222e+05, 1.92577546e+05, 3.65060000e+02, 1.73776155e+05,
         2.55327990e+05, 4.84013000e+02, 3.59032000e+02, 5.27523000e+02,
         1.00000000e+00};
        A[192]= new double[]{1.46335940e+05, 3.02023829e+04, 3.95967000e+02, 2.73892754e+04,
         5.65289280e+03, 7.41120000e+01, 3.69566000e+02, 7.62750000e+01,
         1.00000000e+00};
        A[193]= new double[]{1.20736040e+05, 1.28704291e+05, 3.50777000e+02, 1.16102130e+05,
         1.23764554e+05, 3.37314000e+02, 3.44196000e+02, 3.66912000e+02,
         1.00000000e+00};
        A[194]= new double[]{4.46440213e+05, 1.06238459e+05, 6.64863000e+02, 9.59748791e+04,
         2.28389445e+04, 1.42931000e+02, 6.71477000e+02, 1.59790000e+02,
         1.00000000e+00};
        A[195]= new double[]{4.87471671e+04, 1.23484044e+05, 2.26269000e+02, 1.06561084e+05,
         2.69935556e+05, 4.94623000e+02, 2.15439000e+02, 5.45740000e+02,
         1.00000000e+00};
        A[196]= new double[]{1.81977271e+05, 1.35172392e+05, 4.29939000e+02, 1.23018428e+05,
         9.13778686e+04, 2.90643000e+02, 4.23263000e+02, 3.14399000e+02,
         1.00000000e+00};
        A[197]= new double[]{1.80042267e+03, 3.32176109e+04, 5.85770000e+01, 1.55438407e+04,
         2.86782242e+05, 5.05721000e+02, 3.07360000e+01, 5.67076000e+02,
         1.00000000e+00};
        A[198]= new double[]{8.85334881e+05, 5.85566283e+04, 9.25226000e+02, 3.95805911e+04,
         2.61788620e+03, 4.13640000e+01, 9.56885000e+02, 6.32890000e+01,
         1.00000000e+00};
        A[199]= new double[]{1.92290195e+03, 2.24677159e+04, 5.92300000e+01, 1.11669861e+04,
         1.30478140e+05, 3.43970000e+02, 3.24650000e+01, 3.79330000e+02,
         1.00000000e+00};
        A[200]= new double[]{4.29208793e+03, 2.47523440e+04, 7.85320000e+01, 1.57896499e+04,
         9.10584436e+04, 2.88902000e+02, 5.46540000e+01, 3.15188000e+02,
         1.00000000e+00};
        A[201]= new double[]{3.68849796e+05, 3.05579561e+05, 6.05011000e+02, 2.85139485e+05,
         2.36228404e+05, 4.67704000e+02, 6.09658000e+02, 5.05081000e+02,
         1.00000000e+00};
        A[202]= new double[]{6.42650809e+05, 3.20494369e+05, 7.85849000e+02, 3.04537628e+05,
         1.51875005e+05, 3.72396000e+02, 8.17779000e+02, 4.07832000e+02,
         1.00000000e+00};
        A[203]= new double[]{1.04486932e+05, 1.44094103e+04, 3.43261000e+02, 1.34399524e+04,
         1.85345463e+03, 4.41530000e+01, 3.04395000e+02, 4.19780000e+01,
         1.00000000e+00};
        A[204]= new double[]{4.41548686e+05, 1.06657188e+05, 6.61198000e+02, 9.66221233e+04,
         2.33393153e+04, 1.44687000e+02, 6.67801000e+02, 1.61309000e+02,
         1.00000000e+00};
        A[205]= new double[]{1.70440442e+04, 6.40197020e+04, 1.39568000e+02, 5.06469497e+04,
         1.90236695e+05, 4.14731000e+02, 1.22120000e+02, 4.58699000e+02,
         1.00000000e+00};
        A[206]= new double[]{1.38316422e+05, 2.79775147e+04, 3.85275000e+02, 2.52389101e+04,
         5.10512033e+03, 7.03020000e+01, 3.59007000e+02, 7.26170000e+01,
         1.00000000e+00};
        A[207]= new double[]{1.06074479e+05, 1.57253817e+05, 3.29386000e+02, 1.40896662e+05,
         2.08877179e+05, 4.37517000e+02, 3.22037000e+02, 4.77415000e+02,
         1.00000000e+00};
        A[208]= new double[]{1.73947751e+05, 1.65118727e+05, 4.20630000e+02, 1.49896620e+05,
         1.42288354e+05, 3.62471000e+02, 4.13541000e+02, 3.92551000e+02,
         1.00000000e+00};
        A[209]= new double[]{1.25790946e+04, 6.86898137e+04, 1.22228000e+02, 5.18884051e+04,
         2.83343514e+05, 5.04187000e+02, 1.02915000e+02, 5.61981000e+02,
         1.00000000e+00};
        A[210]= new double[]{5.16383443e+05, 2.33640814e+05, 7.08172000e+02, 2.18819755e+05,
         9.90063228e+04, 3.00091000e+02, 7.29178000e+02, 3.29921000e+02,
         1.00000000e+00};
        A[211]= new double[]{8.15003647e+05, 1.30090074e+05, 8.83536000e+02, 1.12893008e+05,
         1.80198699e+04, 1.22386000e+02, 9.22434000e+02, 1.47238000e+02,
         1.00000000e+00};
        A[212]= new double[]{1.18506226e+05, 1.57827605e+05, 3.47494000e+02, 1.42104889e+05,
         1.89256506e+05, 4.16692000e+02, 3.41031000e+02, 4.54188000e+02,
         1.00000000e+00};
        A[213]= new double[]{3.70307575e+05, 2.10427720e+04, 6.27955000e+02, 1.60428973e+04,
         9.11639550e+02, 2.72050000e+01, 5.89704000e+02, 3.35100000e+01,
         1.00000000e+00};
        A[214]= new double[]{1.83302547e+05, 9.18921505e+04, 4.32454000e+02, 8.33600308e+04,
         4.17895583e+04, 1.96666000e+02, 4.23866000e+02, 2.12490000e+02,
         1.00000000e+00};
        A[215]= new double[]{4.32614755e+05, 4.59145838e+04, 6.60651000e+02, 3.82283789e+04,
         4.05728212e+03, 5.83790000e+01, 6.54831000e+02, 6.94990000e+01,
         1.00000000e+00};
        A[216]= new double[]{5.27086293e+03, 2.27082169e+04, 8.47910000e+01, 1.53786289e+04,
         6.62550411e+04, 2.47392000e+02, 6.21630000e+01, 2.67814000e+02,
         1.00000000e+00};
        A[217]= new double[]{2.24154566e+05, 2.10415131e+05, 4.77694000e+02, 1.91577370e+05,
         1.79834737e+05, 4.08269000e+02, 4.69243000e+02, 4.40481000e+02,
         1.00000000e+00};
        A[218]= new double[]{4.80323364e+05, 3.98980160e+04, 6.94422000e+02, 3.18626077e+04,
         2.64666457e+03, 4.60650000e+01, 6.91688000e+02, 5.74550000e+01,
         1.00000000e+00};
        A[219]= new double[]{1.18711158e+03, 1.94258984e+04, 5.10410000e+01, 8.01928863e+03,
         1.31227669e+05, 3.44797000e+02, 2.32580000e+01, 3.80594000e+02,
         1.00000000e+00};
        A[220]= new double[]{4.58497333e+05, 1.98486667e+04, 6.93791000e+02, 1.29303476e+04,
         5.59763694e+02, 1.95660000e+01, 6.60858000e+02, 2.86090000e+01,
         1.00000000e+00};
        A[221]= new double[]{7.23679835e+05, 4.23475429e+05, 8.31923000e+02, 4.06064588e+05,
         2.37616647e+05, 4.66801000e+02, 8.69888000e+02, 5.09032000e+02,
         1.00000000e+00};
        A[222]= new double[]{5.66258400e+05, 3.20806005e+05, 7.39769000e+02, 3.04136675e+05,
         1.72304502e+05, 3.97329000e+02, 7.65453000e+02, 4.33657000e+02,
         1.00000000e+00};
        A[223]= new double[]{1.19568760e+04, 3.76644231e+04, 1.19842000e+02, 2.86757698e+04,
         9.03293073e+04, 2.87413000e+02, 9.97720000e+01, 3.14284000e+02,
         1.00000000e+00};
        A[224]= new double[]{9.31169940e+04, 5.73906060e+03, 3.27012000e+02, 6.36475435e+03,
         3.92277600e+02, 2.23520000e+01, 2.84751000e+02, 1.75500000e+01,
         1.00000000e+00};
        A[225]= new double[]{4.39963555e+05, 2.53538987e+05, 6.56517000e+02, 2.37780573e+05,
         1.37026454e+05, 3.54818000e+02, 6.70148000e+02, 3.86188000e+02,
         1.00000000e+00};
        A[226]= new double[]{1.49104398e+05, 2.80225164e+04, 4.00562000e+02, 2.56710214e+04,
         4.82458351e+03, 6.89640000e+01, 3.72238000e+02, 6.99580000e+01,
         1.00000000e+00};
        A[227]= new double[]{7.03989531e+05, 2.13205687e+05, 8.20868000e+02, 1.98331419e+05,
         6.00653626e+04, 2.31259000e+02, 8.57616000e+02, 2.59732000e+02,
         1.00000000e+00};
        A[228]= new double[]{7.84386128e+04, 9.39455958e+04, 2.84271000e+02, 8.39075255e+04,
         1.00495690e+05, 3.04091000e+02, 2.75929000e+02, 3.30479000e+02,
         1.00000000e+00};
        A[229]= new double[]{2.24039037e+03, 6.43807294e+03, 5.32640000e+01, 4.98144472e+03,
         1.43148734e+04, 1.18431000e+02, 4.20620000e+01, 1.20871000e+02,
         1.00000000e+00};
        A[230]= new double[]{3.16903255e+05, 1.63927960e+04, 5.83893000e+02, 1.25693620e+04,
         6.50188925e+02, 2.31590000e+01, 5.42742000e+02, 2.80750000e+01,
         1.00000000e+00};
        A[231]= new double[]{1.27101774e+05, 1.87723117e+04, 3.73452000e+02, 1.73813170e+04,
         2.56713569e+03, 5.10700000e+01, 3.40343000e+02, 5.02670000e+01,
         1.00000000e+00};
        A[232]= new double[]{8.91644002e+05, 1.94253515e+05, 9.18730000e+02, 1.75245465e+05,
         3.81789677e+04, 1.80569000e+02, 9.70518000e+02, 2.11437000e+02,
         1.00000000e+00};
        A[233]= new double[]{6.21801620e+05, 1.77936771e+05, 7.73312000e+02, 1.64723813e+05,
         4.71379015e+04, 2.04861000e+02, 8.04076000e+02, 2.30097000e+02,
         1.00000000e+00};
        A[234]= new double[]{3.81546101e+03, 2.03187452e+04, 7.53850000e+01, 1.25860865e+04,
         6.70255797e+04, 2.48673000e+02, 5.06130000e+01, 2.69533000e+02,
         1.00000000e+00};
        A[235]= new double[]{3.85533038e+03, 1.87860922e+04, 7.50780000e+01, 1.18732242e+04,
         5.78553490e+04, 2.31217000e+02, 5.13510000e+01, 2.50221000e+02,
         1.00000000e+00};
        A[236]= new double[]{9.20907772e+05, 3.85500695e+05, 9.35938000e+02, 3.66716779e+05,
         1.53511109e+05, 3.72702000e+02, 9.83941000e+02, 4.11887000e+02,
         1.00000000e+00};
        A[237]= new double[]{9.12944358e+05, 4.14712603e+05, 9.32456000e+02, 3.95644208e+05,
         1.79724687e+05, 4.04100000e+02, 9.79075000e+02, 4.44753000e+02,
         1.00000000e+00};
        A[238]= new double[]{7.37451093e+05, 1.25585467e+05, 8.42471000e+02, 1.09905441e+05,
         1.87165309e+04, 1.25557000e+02, 8.75343000e+02, 1.49068000e+02,
         1.00000000e+00};
        A[239]= new double[]{9.27899007e+05, 3.86853779e+05, 9.39636000e+02, 3.67553812e+05,
         1.53238208e+05, 3.72203000e+02, 9.87509000e+02, 4.11706000e+02,
         1.00000000e+00};
        A[240]= new double[]{3.79929231e+04, 1.00889821e+05, 2.02274000e+02, 8.50228630e+04,
         2.25777348e+05, 4.52661000e+02, 1.87829000e+02, 4.98778000e+02,
         1.00000000e+00};
        A[241]= new double[]{3.48090107e+05, 1.15275377e+05, 5.91271000e+02, 1.05058547e+05,
         3.47917487e+04, 1.78454000e+02, 5.88715000e+02, 1.94962000e+02,
         1.00000000e+00};
        A[242]= new double[]{3.79209036e+04, 9.76108892e+04, 2.02067000e+02, 8.23969456e+04,
         2.12095134e+05, 4.39064000e+02, 1.87665000e+02, 4.83062000e+02,
         1.00000000e+00};
        A[243]= new double[]{3.78303120e+04, 9.46050823e+04, 2.02085000e+02, 7.96384368e+04,
         1.99157778e+05, 4.25419000e+02, 1.87200000e+02, 4.68145000e+02,
         1.00000000e+00};
        A[244]= new double[]{2.52620020e+05, 5.72024386e+04, 5.09703000e+02, 5.14663797e+04,
         1.16538761e+04, 1.03842000e+02, 4.95622000e+02, 1.12227000e+02,
         1.00000000e+00};
        A[245]= new double[]{5.49615697e+05, 1.67090743e+05, 7.28569000e+02, 1.54337236e+05,
         4.69206458e+04, 2.04589000e+02, 7.54377000e+02, 2.29341000e+02,
         1.00000000e+00};
        A[246]= new double[]{5.44206640e+05, 1.61929035e+05, 7.25024000e+02, 1.49498748e+05,
         4.44834487e+04, 1.99171000e+02, 7.50605000e+02, 2.23343000e+02,
         1.00000000e+00};
        A[247]= new double[]{3.94698691e+05, 9.76947347e+04, 6.28266000e+02, 8.82563375e+04,
         2.18449660e+04, 1.40483000e+02, 6.28235000e+02, 1.55499000e+02,
         1.00000000e+00};
        A[248]= new double[]{1.28348932e+05, 7.86875604e+04, 3.60919000e+02, 7.19032681e+04,
         4.40821179e+04, 2.02193000e+02, 3.55617000e+02, 2.18020000e+02,
         1.00000000e+00};
        A[249]= new double[]{4.64284004e+05, 1.23712212e+05, 6.75735000e+02, 1.12633711e+05,
         3.00121596e+04, 1.63931000e+02, 6.87080000e+02, 1.83078000e+02,
         1.00000000e+00};
        A[250]= new double[]{1.29006708e+05, 7.65770059e+04, 3.61649000e+02, 7.01193438e+04,
         4.16220946e+04, 1.96568000e+02, 3.56718000e+02, 2.11744000e+02,
         1.00000000e+00};
        A[251]= new double[]{4.25339424e+04, 1.06635447e+05, 2.13034000e+02, 9.08481835e+04,
         2.27762491e+05, 4.55019000e+02, 1.99658000e+02, 5.00556000e+02,
         1.00000000e+00};
        A[252]= new double[]{3.27977924e+04, 9.48541149e+04, 1.88811000e+02, 7.90816751e+04,
         2.28711195e+05, 4.55259000e+02, 1.73707000e+02, 5.02376000e+02,
         1.00000000e+00};
        A[253]= new double[]{5.43985711e+03, 4.27530928e+04, 8.69890000e+01, 2.76228977e+04,
         2.17094729e+05, 4.41719000e+02, 6.25350000e+01, 4.91477000e+02,
         1.00000000e+00};
        A[254]= new double[]{1.89632186e+05, 8.45209316e+04, 4.41464000e+02, 7.65669631e+04,
         3.41266491e+04, 1.78248000e+02, 4.29553000e+02, 1.91456000e+02,
         1.00000000e+00};
        A[255]= new double[]{7.91636406e+04, 4.53220438e+04, 2.88200000e+02, 4.09263936e+04,
         2.34308047e+04, 1.48995000e+02, 2.74683000e+02, 1.57259000e+02,
         1.00000000e+00};
        A[256]= new double[]{8.05952851e+03, 3.02132547e+04, 1.01808000e+02, 2.16040931e+04,
         8.09886046e+04, 2.72903000e+02, 7.91640000e+01, 2.96767000e+02,
         1.00000000e+00};
        A[257]= new double[]{7.13108734e+05, 1.39260431e+05, 8.27788000e+02, 1.24384920e+05,
         2.42906820e+04, 1.44388000e+02, 8.61463000e+02, 1.68232000e+02,
         1.00000000e+00};
        A[258]= new double[]{1.14067100e+04, 6.84956177e+03, 1.17283000e+02, 6.29735824e+03,
         3.78147110e+03, 6.47490000e+01, 9.72580000e+01, 5.84020000e+01,
         1.00000000e+00};
        A[259]= new double[]{1.73941484e+04, 2.14370853e+04, 1.40619000e+02, 1.81015716e+04,
         2.23089354e+04, 1.46338000e+02, 1.23697000e+02, 1.52448000e+02,
         1.00000000e+00};
        A[260]= new double[]{6.55630231e+04, 4.70219661e+04, 2.60083000e+02, 4.26885781e+04,
         3.06163562e+04, 1.69342000e+02, 2.52085000e+02, 1.80796000e+02,
         1.00000000e+00};
        A[261]= new double[]{5.80424811e+05, 1.45065928e+05, 7.49292000e+02, 1.32195428e+05,
         3.30396842e+04, 1.70656000e+02, 7.74631000e+02, 1.93604000e+02,
         1.00000000e+00};
        A[262]= new double[]{6.92507595e+05, 4.11219995e+04, 8.16107000e+02, 2.80310007e+04,
         1.66451719e+03, 3.30340000e+01, 8.48550000e+02, 5.03880000e+01,
         1.00000000e+00};
        A[263]= new double[]{7.35299849e+05, 1.39842133e+05, 8.39928000e+02, 1.24548586e+05,
         2.36871256e+04, 1.42271000e+02, 8.75432000e+02, 1.66493000e+02,
         1.00000000e+00};
        A[264]= new double[]{3.06664188e+04, 5.56479601e+04, 1.83023000e+02, 4.68803810e+04,
         8.50701738e+04, 2.79791000e+02, 1.67555000e+02, 3.04049000e+02,
         1.00000000e+00};
        A[265]= new double[]{1.42825291e+04, 6.13995900e+03, 1.29535000e+02, 6.05007646e+03,
         2.60088540e+03, 5.48710000e+01, 1.10260000e+02, 4.74000000e+01,
         1.00000000e+00};
        A[266]= new double[]{3.03776069e+04, 7.33215108e+04, 1.81413000e+02, 6.15485918e+04,
         1.48557974e+05, 3.67564000e+02, 1.67450000e+02, 4.04169000e+02,
         1.00000000e+00};
        A[267]= new double[]{1.74167136e+05, 7.89901437e+04, 4.23551000e+02, 7.15113645e+04,
         3.24325995e+04, 1.73906000e+02, 4.11207000e+02, 1.86495000e+02,
         1.00000000e+00};
        A[268]= new double[]{3.23109682e+05, 1.15589465e+05, 5.70846000e+02, 1.05295383e+05,
         3.76684377e+04, 1.86028000e+02, 5.66019000e+02, 2.02488000e+02,
         1.00000000e+00};
        A[269]= new double[]{9.02079205e+04, 5.46296778e+04, 3.05571000e+02, 4.95954480e+04,
         3.00348720e+04, 1.68000000e+02, 2.95211000e+02, 1.78779000e+02,
         1.00000000e+00};
        A[270]= new double[]{6.48991883e+05, 1.58583908e+05, 7.89116000e+02, 1.45048513e+05,
         3.54432168e+04, 1.76366000e+02, 8.22429000e+02, 2.00964000e+02,
         1.00000000e+00};
        A[271]= new double[]{3.68242706e+05, 1.57616468e+05, 6.04038000e+02, 1.45733247e+05,
         6.23772289e+04, 2.39050000e+02, 6.09635000e+02, 2.60938000e+02,
         1.00000000e+00};
        A[272]= new double[]{2.52188314e+05, 7.41171166e+04, 5.09424000e+02, 6.66747755e+04,
         1.95954445e+04, 1.34684000e+02, 4.95046000e+02, 1.45492000e+02,
         1.00000000e+00};
        A[273]= new double[]{4.30038938e+05, 2.60941106e+05, 6.50614000e+02, 2.44157847e+05,
         1.48151279e+05, 3.69391000e+02, 6.60974000e+02, 4.01069000e+02,
         1.00000000e+00};
        A[274]= new double[]{1.85707644e+05, 1.97240527e+05, 4.34220000e+02, 1.79324505e+05,
         1.90460980e+05, 4.19295000e+02, 4.27681000e+02, 4.54241000e+02,
         1.00000000e+00};
        A[275]= new double[]{7.95587553e+05, 9.11553983e+04, 8.75972000e+02, 7.43607505e+04,
         8.51997219e+03, 8.18740000e+01, 9.08234000e+02, 1.04062000e+02,
         1.00000000e+00};
        A[276]= new double[]{7.54708477e+05, 1.50358944e+05, 8.49740000e+02, 1.34364114e+05,
         2.67690730e+04, 1.51283000e+02, 8.88164000e+02, 1.76947000e+02,
         1.00000000e+00};
        A[277]= new double[]{6.25041742e+05, 3.18117446e+05, 7.75310000e+02, 3.02027593e+05,
         1.53718128e+05, 3.74639000e+02, 8.06183000e+02, 4.10310000e+02,
         1.00000000e+00};
        A[278]= new double[]{8.24262533e+05, 3.41568361e+05, 8.86118000e+02, 3.23954362e+05,
         1.34244316e+05, 3.48265000e+02, 9.30195000e+02, 3.85466000e+02,
         1.00000000e+00};
        A[279]= new double[]{6.59745622e+05, 1.35566267e+05, 7.97641000e+02, 1.21536333e+05,
         2.49736055e+04, 1.46939000e+02, 8.27121000e+02, 1.69959000e+02,
         1.00000000e+00};
        A[280]= new double[]{8.50706010e+03, 4.05472815e+03, 1.03813000e+02, 3.92545924e+03,
         1.87099537e+03, 4.79030000e+01, 8.19460000e+01, 3.90580000e+01,
         1.00000000e+00};
        A[281]= new double[]{1.66592519e+05, 1.91933233e+04, 4.21285000e+02, 1.77765648e+04,
         2.04805929e+03, 4.49540000e+01, 3.95439000e+02, 4.55590000e+01,
         1.00000000e+00};
        A[282]= new double[]{1.61174908e+05, 5.41836239e+04, 4.06729000e+02, 4.94871150e+04,
         1.66365303e+04, 1.24882000e+02, 3.96271000e+02, 1.33218000e+02,
         1.00000000e+00};
        A[283]= new double[]{4.76589463e+05, 2.65130313e+05, 6.81781000e+02, 2.49350335e+05,
         1.38715473e+05, 3.56706000e+02, 6.99036000e+02, 3.88879000e+02,
         1.00000000e+00};
        A[284]= new double[]{1.07605678e+04, 4.83917752e+03, 1.14423000e+02, 4.77150300e+03,
         2.14581150e+03, 5.07380000e+01, 9.40420000e+01, 4.22920000e+01,
         1.00000000e+00};
        A[285]= new double[]{2.32124277e+05, 9.43261800e+04, 4.88800000e+02, 8.53598087e+04,
         3.46868703e+04, 1.79748000e+02, 4.74886000e+02, 1.92975000e+02,
         1.00000000e+00};
        A[286]= new double[]{3.03536860e+05, 1.54328448e+05, 5.53099000e+02, 1.41550727e+05,
         7.19691973e+04, 2.57931000e+02, 5.48793000e+02, 2.79025000e+02,
         1.00000000e+00};
        A[287]= new double[]{4.63573397e+05, 1.46719418e+05, 6.72559000e+02, 1.35235760e+05,
         4.28016625e+04, 1.96202000e+02, 6.89268000e+02, 2.18151000e+02,
         1.00000000e+00};
        A[288]= new double[]{2.91160148e+05, 1.46464283e+05, 5.42127000e+02, 1.33902829e+05,
         6.73580573e+04, 2.49321000e+02, 5.37070000e+02, 2.70166000e+02,
         1.00000000e+00};
        A[289]= new double[]{3.33921860e+05, 9.40219381e+04, 5.80214000e+02, 8.51428401e+04,
         2.39735573e+04, 1.47942000e+02, 5.75515000e+02, 1.62047000e+02,
         1.00000000e+00};
        A[290]= new double[]{9.28448298e+03, 4.37789645e+04, 1.06953000e+02, 3.21625609e+04,
         1.51655576e+05, 3.70498000e+02, 8.68090000e+01, 4.09329000e+02,
         1.00000000e+00};
        A[291]= new double[]{6.37367361e+04, 8.11461597e+04, 2.57414000e+02, 7.18546808e+04,
         9.14814872e+04, 2.90200000e+02, 2.47604000e+02, 3.15236000e+02,
         1.00000000e+00};
        A[292]= new double[]{4.41464930e+05, 1.30977297e+05, 6.59264000e+02, 1.19727032e+05,
         3.55215602e+04, 1.78795000e+02, 6.69633000e+02, 1.98672000e+02,
         1.00000000e+00};
        A[293]= new double[]{1.46415644e+05, 1.06340649e+05, 3.86578000e+02, 9.64773418e+04,
         7.00708126e+04, 2.54727000e+02, 3.78748000e+02, 2.75082000e+02,
         1.00000000e+00};
        A[294]= new double[]{3.66118559e+05, 2.35827939e+05, 6.02578000e+02, 2.19372932e+05,
         1.41304681e+05, 3.61056000e+02, 6.07587000e+02, 3.91365000e+02,
         1.00000000e+00};
        A[295]= new double[]{1.21837838e+05, 1.45702288e+05, 3.52274000e+02, 1.31500848e+05,
         1.57257998e+05, 3.80213000e+02, 3.45861000e+02, 4.13605000e+02,
         1.00000000e+00};
        A[296]= new double[]{7.07975126e+04, 8.48486854e+04, 2.70960000e+02, 7.52584144e+04,
         9.01949417e+04, 2.88033000e+02, 2.61284000e+02, 3.13141000e+02,
         1.00000000e+00};
        A[297]= new double[]{7.36878308e+05, 1.05283850e+05, 8.42163000e+02, 9.03524945e+04,
         1.29094022e+04, 1.03262000e+02, 8.74983000e+02, 1.25016000e+02,
         1.00000000e+00};
        A[298]= new double[]{1.47666739e+05, 1.18295109e+05, 3.88185000e+02, 1.07106269e+05,
         8.58023128e+04, 2.81560000e+02, 3.80403000e+02, 3.04739000e+02,
         1.00000000e+00};
        A[299]= new double[]{8.43741464e+05, 1.55106619e+05, 8.96311000e+02, 1.37572508e+05,
         2.52902192e+04, 1.46144000e+02, 9.41349000e+02, 1.73050000e+02,
         1.00000000e+00};
        A[300]= new double[]{3.26357195e+04, 9.15104715e+04, 1.88304000e+02, 7.63900520e+04,
         2.14197505e+05, 4.40761000e+02, 1.73314000e+02, 4.85972000e+02,
         1.00000000e+00};
        A[301]= new double[]{1.86994852e+05, 2.15457878e+05, 4.35921000e+02, 1.95843113e+05,
         2.25652958e+05, 4.56548000e+02, 4.28965000e+02, 4.94259000e+02,
         1.00000000e+00};
        A[302]= new double[]{5.58114533e+05, 1.95557960e+05, 7.34570000e+02, 1.81874055e+05,
         6.37269181e+04, 2.39376000e+02, 7.59784000e+02, 2.66221000e+02,
         1.00000000e+00};
        A[303]= new double[]{1.58724026e+05, 3.31543972e+04, 4.11728000e+02, 2.99635316e+04,
         6.25880562e+03, 7.77250000e+01, 3.85507000e+02, 8.05250000e+01,
         1.00000000e+00};
        A[304]= new double[]{6.42492122e+05, 8.29053077e+04, 7.87280000e+02, 7.01332284e+04,
         9.04978703e+03, 8.59380000e+01, 8.16091000e+02, 1.05306000e+02,
         1.00000000e+00};
        A[305]= new double[]{2.93287910e+04, 7.84160569e+04, 1.78361000e+02, 6.55380358e+04,
         1.75228305e+05, 3.98565000e+02, 1.64435000e+02, 4.39648000e+02,
         1.00000000e+00};
        A[306]= new double[]{7.19940119e+04, 4.71275465e+04, 2.73905000e+02, 4.25803032e+04,
         2.78732239e+04, 1.61999000e+02, 2.62843000e+02, 1.72058000e+02,
         1.00000000e+00};
        A[307]= new double[]{3.20878926e+05, 1.57409358e+05, 5.67648000e+02, 1.44483926e+05,
         7.08775810e+04, 2.55598000e+02, 5.65278000e+02, 2.77301000e+02,
         1.00000000e+00};
        A[308]= new double[]{9.31447363e+03, 4.61296236e+04, 1.07254000e+02, 3.37503987e+04,
         1.67147737e+05, 3.88628000e+02, 8.68450000e+01, 4.30097000e+02,
         1.00000000e+00};
        A[309]= new double[]{2.68435574e+05, 3.11902106e+04, 5.25415000e+02, 2.74763096e+04,
         3.19254214e+03, 5.37800000e+01, 5.10902000e+02, 5.93630000e+01,
         1.00000000e+00};
        A[310]= new double[]{5.67938654e+05, 3.38384051e+05, 7.41027000e+02, 3.21279851e+05,
         1.91422043e+05, 4.19195000e+02, 7.66421000e+02, 4.56642000e+02,
         1.00000000e+00};
        A[311]= new double[]{4.61854955e+04, 8.54042466e+04, 2.20886000e+02, 7.37009573e+04,
         1.36284664e+05, 3.52481000e+02, 2.09092000e+02, 3.86644000e+02,
         1.00000000e+00};
        A[312]= new double[]{3.71410640e+05, 2.31933978e+05, 6.06357000e+02, 2.15923470e+05,
         1.34837250e+05, 3.52512000e+02, 6.12528000e+02, 3.82504000e+02,
         1.00000000e+00};
        A[313]= new double[]{2.72862618e+04, 2.59081239e+04, 1.73569000e+02, 2.25276059e+04,
         2.13898118e+04, 1.43299000e+02, 1.57207000e+02, 1.49267000e+02,
         1.00000000e+00};
        A[314]= new double[]{1.32616326e+05, 8.95800366e+04, 3.68121000e+02, 8.11471233e+04,
         5.48134793e+04, 2.25251000e+02, 3.60252000e+02, 2.43344000e+02,
         1.00000000e+00};
        A[315]= new double[]{1.20271617e+05, 1.48162484e+05, 3.50203000e+02, 1.33513058e+05,
         1.64474603e+05, 3.88759000e+02, 3.43434000e+02, 4.23076000e+02,
         1.00000000e+00};
        A[316]= new double[]{1.30095139e+05, 1.10277142e+05, 3.64355000e+02, 9.99663965e+04,
         8.47380507e+04, 2.79974000e+02, 3.57056000e+02, 3.02664000e+02,
         1.00000000e+00};
        A[317]= new double[]{6.32042692e+05, 3.32589745e+05, 7.80185000e+02, 3.16217800e+05,
         1.66398249e+05, 3.90335000e+02, 8.10119000e+02, 4.26296000e+02,
         1.00000000e+00};
        A[318]= new double[]{4.29465345e+04, 2.18038266e+04, 2.18191000e+02, 1.95090023e+04,
         9.90466188e+03, 9.91160000e+01, 1.96830000e+02, 9.99300000e+01,
         1.00000000e+00};
        A[319]= new double[]{2.46224148e+05, 2.60397114e+04, 5.03212000e+02, 2.30785596e+04,
         2.44069900e+03, 4.71660000e+01, 4.89305000e+02, 5.17470000e+01,
         1.00000000e+00};
        A[320]= new double[]{1.06921810e+04, 9.73833546e+03, 1.13892000e+02, 8.30593912e+03,
         7.56496937e+03, 8.84740000e+01, 9.38800000e+01, 8.55050000e+01,
         1.00000000e+00};
        A[321]= new double[]{1.22474047e+05, 1.75114815e+05, 3.53116000e+02, 1.57816839e+05,
         2.25648350e+05, 4.55016000e+02, 3.46838000e+02, 4.95913000e+02,
         1.00000000e+00};
        A[322]= new double[]{3.81997823e+05, 1.62829342e+05, 6.14506000e+02, 1.50776705e+05,
         6.42696638e+04, 2.42549000e+02, 6.21634000e+02, 2.64976000e+02,
         1.00000000e+00};
        A[323]= new double[]{2.14784997e+05, 2.40871074e+04, 4.71455000e+02, 2.17169954e+04,
         2.43545688e+03, 4.76690000e+01, 4.55579000e+02, 5.10910000e+01,
         1.00000000e+00};
        A[324]= new double[]{2.78559523e+05, 2.80425422e+05, 5.21201000e+02, 2.63943193e+05,
         2.65711187e+05, 4.93853000e+02, 5.34457000e+02, 5.38037000e+02,
         1.00000000e+00};
        A[325]= new double[]{4.12876632e+05, 2.49022096e+05, 6.37539000e+02, 2.32677854e+05,
         1.40337143e+05, 3.59287000e+02, 6.47610000e+02, 3.90599000e+02,
         1.00000000e+00};
        A[326]= new double[]{1.22558640e+05, 1.15752103e+05, 3.53311000e+02, 1.04793220e+05,
         9.89733212e+04, 3.02097000e+02, 3.46886000e+02, 3.27621000e+02,
         1.00000000e+00};
        A[327]= new double[]{7.08053956e+04, 7.20754938e+04, 2.71041000e+02, 6.41457318e+04,
         6.52963697e+04, 2.45548000e+02, 2.61235000e+02, 2.65921000e+02,
         1.00000000e+00};
        A[328]= new double[]{7.12260623e+04, 6.44595973e+04, 2.71822000e+02, 5.77502806e+04,
         5.22640128e+04, 2.20394000e+02, 2.62032000e+02, 2.37139000e+02,
         1.00000000e+00};
        A[329]= new double[]{8.90011355e+05, 1.51364485e+05, 9.20188000e+02, 1.32761597e+05,
         2.25788027e+04, 1.37263000e+02, 9.67206000e+02, 1.64493000e+02,
         1.00000000e+00};
        A[330]= new double[]{1.39233564e+05, 8.14627440e+04, 3.75931000e+02, 7.43525182e+04,
         4.35021554e+04, 2.00752000e+02, 3.70370000e+02, 2.16696000e+02,
         1.00000000e+00};
        A[331]= new double[]{5.64193800e+05, 3.53036852e+05, 7.38577000e+02, 3.35749307e+05,
         2.10090714e+05, 4.39524000e+02, 7.63893000e+02, 4.77996000e+02,
         1.00000000e+00};
        A[332]= new double[]{3.54247175e+05, 1.37623389e+05, 5.93056000e+02, 1.26616175e+05,
         4.91897984e+04, 2.11972000e+02, 5.97325000e+02, 2.32058000e+02,
         1.00000000e+00};
        A[333]= new double[]{8.21947915e+03, 2.14206907e+04, 1.02509000e+02, 1.57454555e+04,
         4.10340517e+04, 1.96369000e+02, 8.01830000e+01, 2.08964000e+02,
         1.00000000e+00};
        A[334]= new double[]{8.14490232e+05, 1.11974102e+05, 8.84849000e+02, 9.44776599e+04,
         1.29885549e+04, 1.02639000e+02, 9.20485000e+02, 1.26546000e+02,
         1.00000000e+00};
        A[335]= new double[]{2.79187524e+04, 7.67352351e+04, 1.74405000e+02, 6.38699990e+04,
         1.75547937e+05, 3.98988000e+02, 1.60080000e+02, 4.39983000e+02,
         1.00000000e+00};
        A[336]= new double[]{5.41696167e+03, 3.75136361e+04, 8.66700000e+01, 2.44318284e+04,
         1.69195718e+05, 3.90903000e+02, 6.25010000e+01, 4.32833000e+02,
         1.00000000e+00};
        A[337]= new double[]{5.07105390e+05, 1.51479907e+05, 7.01442000e+02, 1.39566364e+05,
         4.16905447e+04, 1.93052000e+02, 7.22947000e+02, 2.15955000e+02,
         1.00000000e+00};
        A[338]= new double[]{4.14603404e+05, 3.03170297e+05, 6.38662000e+02, 2.84685309e+05,
         2.08170336e+05, 4.38534000e+02, 6.49175000e+02, 4.74696000e+02,
         1.00000000e+00};
        A[339]= new double[]{3.32807855e+05, 1.02194850e+05, 5.79175000e+02, 9.25340012e+04,
         2.84142883e+04, 1.61034000e+02, 5.74624000e+02, 1.76449000e+02,
         1.00000000e+00};
        A[340]= new double[]{8.00978978e+05, 3.27221952e+05, 8.73470000e+02, 3.10059662e+05,
         1.26667903e+05, 3.38121000e+02, 9.17008000e+02, 3.74623000e+02,
         1.00000000e+00};
        A[341]= new double[]{5.03791941e+05, 1.47213658e+05, 6.99851000e+02, 1.35361722e+05,
         3.95542140e+04, 1.88040000e+02, 7.19856000e+02, 2.10350000e+02,
         1.00000000e+00};
        A[342]= new double[]{2.49982793e+04, 8.97388803e+04, 1.65902000e+02, 7.35481495e+04,
         2.64023316e+05, 4.88105000e+02, 1.50681000e+02, 5.40915000e+02,
         1.00000000e+00};
        A[343]= new double[]{2.87778701e+05, 1.24178381e+05, 5.39960000e+02, 1.13108606e+05,
         4.88070988e+04, 2.12226000e+02, 5.32963000e+02, 2.29977000e+02,
         1.00000000e+00};
        A[344]= new double[]{1.02262571e+05, 1.66623117e+05, 3.23272000e+02, 1.49225814e+05,
         2.43143410e+05, 4.71732000e+02, 3.16336000e+02, 5.15427000e+02,
         1.00000000e+00};
        A[345]= new double[]{9.82203748e+04, 5.65895177e+04, 3.18493000e+02, 5.13563532e+04,
         2.95888839e+04, 1.66530000e+02, 3.08391000e+02, 1.77679000e+02,
         1.00000000e+00};
        A[346]= new double[]{5.92658524e+05, 1.29435931e+05, 7.58270000e+02, 1.16577722e+05,
         2.54604386e+04, 1.49154000e+02, 7.81593000e+02, 1.70699000e+02,
         1.00000000e+00};
        A[347]= new double[]{5.14620022e+03, 1.83057536e+04, 8.49250000e+01, 1.22351403e+04,
         4.35221043e+04, 2.01910000e+02, 6.05970000e+01, 2.15552000e+02,
         1.00000000e+00};
        A[348]= new double[]{6.14709138e+04, 1.10750943e+05, 2.52687000e+02, 9.71918040e+04,
         1.75108573e+05, 3.99524000e+02, 2.43269000e+02, 4.38293000e+02,
         1.00000000e+00};
        A[349]= new double[]{4.07062460e+05, 1.42550627e+05, 6.33170000e+02, 1.31373869e+05,
         4.60062749e+04, 2.04347000e+02, 6.42896000e+02, 2.25138000e+02,
         1.00000000e+00};
        A[350]= new double[]{7.94734937e+05, 3.56823819e+05, 8.70680000e+02, 3.39189928e+05,
         1.52291085e+05, 3.71603000e+02, 9.12775000e+02, 4.09822000e+02,
         1.00000000e+00};
        A[351]= new double[]{5.07965867e+05, 1.66575100e+05, 7.02401000e+02, 1.54073118e+05,
         5.05245462e+04, 2.13048000e+02, 7.23185000e+02, 2.37151000e+02,
         1.00000000e+00};
        A[352]= new double[]{7.28502472e+05, 3.35135368e+05, 8.34409000e+02, 3.18402086e+05,
         1.46475550e+05, 3.64690000e+02, 8.73076000e+02, 4.01644000e+02,
         1.00000000e+00};
        A[353]= new double[]{3.97870453e+05, 2.57114261e+05, 6.26396000e+02, 2.40550557e+05,
         1.55450042e+05, 3.78716000e+02, 6.35174000e+02, 4.10466000e+02,
         1.00000000e+00};
        A[354]= new double[]{6.95314607e+05, 3.75067268e+05, 8.16210000e+02, 3.57930149e+05,
         1.93075022e+05, 4.20164000e+02, 8.51882000e+02, 4.59523000e+02,
         1.00000000e+00};
        A[355]= new double[]{5.58610640e+05, 1.08230497e+05, 7.39032000e+02, 9.61123955e+04,
         1.86217226e+04, 1.27155000e+02, 7.55868000e+02, 1.46449000e+02,
         1.00000000e+00};
        A[356]= new double[]{1.63626981e+05, 1.42689059e+05, 4.07527000e+02, 1.29823686e+05,
         1.13211277e+05, 3.23337000e+02, 4.01512000e+02, 3.50134000e+02,
         1.00000000e+00};
        A[357]= new double[]{3.44466110e+04, 3.42794873e+04, 1.91436000e+02, 3.04001652e+04,
         3.02526736e+04, 1.68948000e+02, 1.79938000e+02, 1.79065000e+02,
         1.00000000e+00};
        A[358]= new double[]{8.26670837e+03, 8.37706640e+03, 1.02468000e+02, 6.87222371e+03,
         6.96396580e+03, 8.51830000e+01, 8.06760000e+01, 8.17530000e+01,
         1.00000000e+00};
        A[359]= new double[]{6.57174646e+05, 1.29709053e+05, 7.96362000e+02, 1.15917143e+05,
         2.28790064e+04, 1.40468000e+02, 8.25221000e+02, 1.62877000e+02,
         1.00000000e+00};
        A[360]= new double[]{1.62848390e+04, 7.43911354e+04, 1.36806000e+02, 5.82038426e+04,
         2.65882268e+05, 4.88960000e+02, 1.19036000e+02, 5.43771000e+02,
         1.00000000e+00};
        A[361]= new double[]{7.30853635e+03, 5.05319218e+04, 9.72630000e+01, 3.50316513e+04,
         2.42212199e+05, 4.66206000e+02, 7.51420000e+01, 5.19539000e+02,
         1.00000000e+00};
        A[362]= new double[]{4.13850481e+05, 1.40555593e+05, 6.37909000e+02, 1.29768419e+05,
         4.40731084e+04, 2.00025000e+02, 6.48761000e+02, 2.20338000e+02,
         1.00000000e+00};
        A[363]= new double[]{1.59429434e+05, 2.06989110e+05, 4.02625000e+02, 1.87549995e+05,
         2.43498364e+05, 4.73641000e+02, 3.95975000e+02, 5.14099000e+02,
         1.00000000e+00};
        A[364]= new double[]{2.27360204e+05, 2.44630113e+05, 4.80922000e+02, 2.23168723e+05,
         2.40120253e+05, 4.72056000e+02, 4.72759000e+02, 5.08669000e+02,
         1.00000000e+00};
        A[365]= new double[]{2.63663201e+04, 3.71469672e+04, 1.70815000e+02, 3.13143561e+04,
         4.41181535e+04, 2.02871000e+02, 1.54356000e+02, 2.17469000e+02,
         1.00000000e+00};
        A[366]= new double[]{9.09379757e+05, 1.48422106e+05, 9.30033000e+02, 1.28849650e+05,
         2.10298683e+04, 1.31776000e+02, 9.77793000e+02, 1.59588000e+02,
         1.00000000e+00};
        A[367]= new double[]{1.13569925e+05, 1.77110754e+05, 3.40466000e+02, 1.58961068e+05,
         2.47897625e+05, 4.76542000e+02, 3.33572000e+02, 5.20201000e+02,
         1.00000000e+00};
        A[368]= new double[]{7.34024854e+03, 4.07138418e+04, 9.72940000e+01, 2.85390318e+04,
         1.58296224e+05, 3.78281000e+02, 7.54440000e+01, 4.18462000e+02,
         1.00000000e+00};
        A[369]= new double[]{8.99329742e+05, 9.58724859e+04, 9.29123000e+02, 7.64106458e+04,
         8.14570921e+03, 7.89420000e+01, 9.67934000e+02, 1.03186000e+02,
         1.00000000e+00};
        A[370]= new double[]{6.92314443e+05, 3.33114584e+05, 8.14732000e+02, 3.16280187e+05,
         1.52181634e+05, 3.72206000e+02, 8.49745000e+02, 4.08864000e+02,
         1.00000000e+00};
        A[371]= new double[]{2.68601568e+05, 2.01513992e+05, 5.22375000e+02, 1.84129428e+05,
         1.38140132e+05, 3.58094000e+02, 5.14193000e+02, 3.85765000e+02,
         1.00000000e+00};
        A[372]= new double[]{4.37873274e+04, 3.04168903e+04, 2.17593000e+02, 2.70147926e+04,
         1.87658401e+04, 1.34245000e+02, 2.01235000e+02, 1.39788000e+02,
         1.00000000e+00};
        A[373]= new double[]{9.26983720e+04, 1.60928634e+05, 3.08182000e+02, 1.43501671e+05,
         2.49125496e+05, 4.77081000e+02, 3.00791000e+02, 5.22187000e+02,
         1.00000000e+00};
        A[374]= new double[]{1.32083973e+05, 5.41072476e+04, 3.71654000e+02, 4.88366039e+04,
         2.00055628e+04, 1.37415000e+02, 3.55395000e+02, 1.45585000e+02,
         1.00000000e+00};
        A[375]= new double[]{8.27285633e+05, 1.03500913e+05, 8.92057000e+02, 8.57382253e+04,
         1.07266273e+04, 9.24510000e+01, 9.27391000e+02, 1.16025000e+02,
         1.00000000e+00};
        A[376]= new double[]{5.55628310e+05, 3.27347883e+05, 7.32986000e+02, 3.10781811e+05,
         1.83096804e+05, 4.09984000e+02, 7.58034000e+02, 4.46595000e+02,
         1.00000000e+00};
        A[377]= new double[]{8.17170557e+03, 2.48933078e+04, 1.02131000e+02, 1.81133566e+04,
         5.51783660e+04, 2.26383000e+02, 8.00120000e+01, 2.43739000e+02,
         1.00000000e+00};
        A[378]= new double[]{7.89836971e+04, 6.94709006e+04, 2.85421000e+02, 6.24708435e+04,
         5.49468551e+04, 2.25749000e+02, 2.76727000e+02, 2.43398000e+02,
         1.00000000e+00};
        A[379]= new double[]{6.74318430e+05, 3.80493872e+05, 8.04602000e+02, 3.63010538e+05,
         2.04833917e+05, 4.33147000e+02, 8.38077000e+02, 4.72897000e+02,
         1.00000000e+00};
        A[380]= new double[]{2.31533945e+04, 2.23914423e+04, 1.62083000e+02, 1.91461943e+04,
         1.85161146e+04, 1.34031000e+02, 1.42849000e+02, 1.38148000e+02,
         1.00000000e+00};
        A[381]= new double[]{7.27909303e+05, 2.77058108e+05, 8.34020000e+02, 2.60806966e+05,
         9.92688019e+04, 2.98826000e+02, 8.72772000e+02, 3.32196000e+02,
         1.00000000e+00};
        A[382]= new double[]{7.01167073e+04, 1.32844051e+05, 2.69223000e+02, 1.17246632e+05,
         2.22137035e+05, 4.50185000e+02, 2.60441000e+02, 4.93435000e+02,
         1.00000000e+00};
        A[383]= new double[]{8.34228206e+05, 1.62441908e+05, 8.90903000e+02, 1.44923370e+05,
         2.82196508e+04, 1.54769000e+02, 9.36385000e+02, 1.82334000e+02,
         1.00000000e+00};
        A[384]= new double[]{6.30553028e+05, 1.10609762e+05, 7.82452000e+02, 9.73488544e+04,
         1.70766504e+04, 1.20800000e+02, 8.05868000e+02, 1.41363000e+02,
         1.00000000e+00};
        A[385]= new double[]{2.50276676e+05, 9.63185564e+04, 5.07605000e+02, 8.67701082e+04,
         3.33933297e+04, 1.75985000e+02, 4.93054000e+02, 1.89751000e+02,
         1.00000000e+00};
        A[386]= new double[]{3.17982718e+05, 2.22020658e+05, 5.64928000e+02, 2.04845245e+05,
         1.43026251e+05, 3.63928000e+02, 5.62873000e+02, 3.93007000e+02,
         1.00000000e+00};
        A[387]= new double[]{1.49263713e+05, 2.05433806e+05, 3.89575000e+02, 1.85857126e+05,
         2.55797848e+05, 4.85083000e+02, 3.83145000e+02, 5.27328000e+02,
         1.00000000e+00};
        A[388]= new double[]{2.81284439e+05, 1.12702315e+05, 5.35319000e+02, 1.02308657e+05,
         4.09920383e+04, 1.94706000e+02, 5.25452000e+02, 2.10533000e+02,
         1.00000000e+00};
        A[389]= new double[]{5.02608673e+05, 3.58221154e+05, 6.99413000e+02, 3.39464384e+05,
         2.41944339e+05, 4.72387000e+02, 7.18615000e+02, 5.12174000e+02,
         1.00000000e+00};
        A[390]= new double[]{7.79095823e+04, 1.34362598e+05, 2.83356000e+02, 1.18970513e+05,
         2.05176139e+05, 4.32694000e+02, 2.74953000e+02, 4.74183000e+02,
         1.00000000e+00};
        A[391]= new double[]{8.81288652e+05, 3.22224626e+05, 9.14943000e+02, 3.04198377e+05,
         1.11223727e+05, 3.15815000e+02, 9.63217000e+02, 3.52180000e+02,
         1.00000000e+00};
        A[392]= new double[]{3.10202802e+04, 6.89386464e+03, 1.91880000e+02, 6.99799285e+03,
         1.55521534e+03, 4.32870000e+01, 1.61665000e+02, 3.59280000e+01,
         1.00000000e+00};
        A[393]= new double[]{6.96313563e+04, 1.42660640e+05, 2.68208000e+02, 1.25722648e+05,
         2.57580411e+05, 4.84262000e+02, 2.59617000e+02, 5.31903000e+02,
         1.00000000e+00};
        A[394]= new double[]{9.86581814e+04, 1.51736881e+05, 3.18006000e+02, 1.35524621e+05,
         2.08437689e+05, 4.36838000e+02, 3.10240000e+02, 4.77151000e+02,
         1.00000000e+00};
        A[395]= new double[]{5.43148187e+04, 9.91688516e+04, 2.38758000e+02, 8.61644161e+04,
         1.57320348e+05, 3.78763000e+02, 2.27489000e+02, 4.15353000e+02,
         1.00000000e+00};
        A[396]= new double[]{2.25505821e+05, 1.97805766e+05, 4.78660000e+02, 1.80361313e+05,
         1.58206594e+05, 3.82836000e+02, 4.71119000e+02, 4.13249000e+02,
         1.00000000e+00};
        A[397]= new double[]{3.99957883e+05, 2.40209916e+05, 6.27834000e+02, 2.24401934e+05,
         1.34773115e+05, 3.52255000e+02, 6.37044000e+02, 3.82601000e+02,
         1.00000000e+00};
        A[398]= new double[]{6.42941376e+04, 6.84402422e+04, 2.58905000e+02, 6.06940830e+04,
         6.46080328e+04, 2.44408000e+02, 2.48331000e+02, 2.64345000e+02,
         1.00000000e+00};
        A[399]= new double[]{1.55121129e+05, 2.01437889e+05, 3.96782000e+02, 1.82675535e+05,
         2.37219613e+05, 4.67263000e+02, 3.90948000e+02, 5.07679000e+02,
         1.00000000e+00};
        A[400]= new double[]{1.09767134e+04, 6.46146093e+04, 1.15148000e+02, 4.79751240e+04,
         2.82406380e+05, 5.03269000e+02, 9.53270000e+01, 5.61144000e+02,
         1.00000000e+00};
        A[401]= new double[]{5.92548909e+04, 1.24629928e+05, 2.48161000e+02, 1.09010079e+05,
         2.29279273e+05, 4.56537000e+02, 2.38776000e+02, 5.02214000e+02,
         1.00000000e+00};
        A[402]= new double[]{2.11254465e+04, 7.09030248e+04, 1.53813000e+02, 5.72771227e+04,
         1.92238363e+05, 4.17031000e+02, 1.37345000e+02, 4.60969000e+02,
         1.00000000e+00};
        A[403]= new double[]{2.43665987e+04, 8.39392073e+04, 1.64139000e+02, 6.85211219e+04,
         2.36044789e+05, 4.61574000e+02, 1.48451000e+02, 5.11391000e+02,
         1.00000000e+00};
        A[404]= new double[]{3.04314368e+04, 3.65400301e+04, 1.81350000e+02, 3.16960152e+04,
         3.80584513e+04, 1.88886000e+02, 1.67805000e+02, 2.01489000e+02,
         1.00000000e+00};
        A[405]= new double[]{8.24253185e+04, 1.17008986e+05, 2.90707000e+02, 1.04541821e+05,
         1.48405038e+05, 3.68710000e+02, 2.83534000e+02, 4.02498000e+02,
         1.00000000e+00};
        A[406]= new double[]{5.73982097e+05, 7.09581332e+04, 7.51938000e+02, 5.95769262e+04,
         7.36515562e+03, 7.80480000e+01, 7.63337000e+02, 9.43670000e+01,
         1.00000000e+00};
        A[407]= new double[]{1.17773421e+05, 1.70530812e+05, 3.46632000e+02, 1.53379774e+05,
         2.22087268e+05, 4.51429000e+02, 3.39765000e+02, 4.91965000e+02,
         1.00000000e+00};
        A[408]= new double[]{4.13205968e+05, 1.30887034e+05, 6.38813000e+02, 1.19954072e+05,
         3.79966262e+04, 1.85448000e+02, 6.46834000e+02, 2.04891000e+02,
         1.00000000e+00};
        A[409]= new double[]{3.42588333e+05, 2.21525340e+05, 5.84515000e+02, 2.05113420e+05,
         1.32630961e+05, 3.49959000e+02, 5.86107000e+02, 3.78990000e+02,
         1.00000000e+00};
        A[410]= new double[]{4.76894611e+05, 3.48959971e+05, 6.82158000e+02, 3.30035305e+05,
         2.41498033e+05, 4.72088000e+02, 6.99097000e+02, 5.11553000e+02,
         1.00000000e+00};
        A[411]= new double[]{3.18009139e+05, 2.66855419e+05, 5.64991000e+02, 2.46502097e+05,
         2.06850723e+05, 4.37948000e+02, 5.62857000e+02, 4.72318000e+02,
         1.00000000e+00};
        A[412]= new double[]{1.63504221e+05, 2.03704679e+05, 4.07477000e+02, 1.84779427e+05,
         2.30210779e+05, 4.60498000e+02, 4.01260000e+02, 4.99917000e+02,
         1.00000000e+00};
        A[413]= new double[]{3.27706096e+05, 1.34884890e+05, 5.71866000e+02, 1.24048630e+05,
         5.10588177e+04, 2.16472000e+02, 5.73047000e+02, 2.35868000e+02,
         1.00000000e+00};
        A[414]= new double[]{7.68942177e+05, 1.17602668e+05, 8.60071000e+02, 1.01679738e+05,
         1.55509853e+04, 1.13730000e+02, 8.94045000e+02, 1.36736000e+02,
         1.00000000e+00};
        A[415]= new double[]{4.88561430e+04, 8.77790430e+04, 2.26657000e+02, 7.62690570e+04,
         1.37031383e+05, 3.53833000e+02, 2.15551000e+02, 3.87277000e+02,
         1.00000000e+00};
        A[416]= new double[]{1.92823346e+04, 8.21818041e+04, 1.47625000e+02, 6.53941848e+04,
         2.78711691e+05, 5.00656000e+02, 1.30617000e+02, 5.56693000e+02,
         1.00000000e+00};
        A[417]= new double[]{1.21546409e+05, 1.37969271e+05, 3.51592000e+02, 1.24815377e+05,
         1.41679929e+05, 3.61048000e+02, 3.45703000e+02, 3.92413000e+02,
         1.00000000e+00};
        A[418]= new double[]{5.93265671e+05, 3.16133936e+05, 7.56399000e+02, 2.99746230e+05,
         1.59726005e+05, 3.82169000e+02, 7.84329000e+02, 4.17946000e+02,
         1.00000000e+00};
        A[419]= new double[]{4.57383579e+04, 1.20474606e+05, 2.19671000e+02, 1.03411693e+05,
         2.72385882e+05, 4.96663000e+02, 2.08213000e+02, 5.48432000e+02,
         1.00000000e+00};
        A[420]= new double[]{2.37352854e+05, 2.14477725e+05, 4.91125000e+02, 1.95586485e+05,
         1.76736633e+05, 4.04703000e+02, 4.83284000e+02, 4.36707000e+02,
         1.00000000e+00};
        A[421]= new double[]{3.76713628e+04, 1.08792595e+05, 2.00546000e+02, 9.22002219e+04,
         2.66268610e+05, 4.90834000e+02, 1.87844000e+02, 5.42482000e+02,
         1.00000000e+00};
        A[422]= new double[]{9.48006917e+04, 6.34315471e+04, 3.10878000e+02, 5.79124099e+04,
         3.87494404e+04, 1.89911000e+02, 3.04945000e+02, 2.04040000e+02,
         1.00000000e+00};
        A[423]= new double[]{8.13572454e+05, 1.85404435e+05, 8.79728000e+02, 1.68404230e+05,
         3.83775177e+04, 1.82098000e+02, 9.24800000e+02, 2.10752000e+02,
         1.00000000e+00};
        A[424]= new double[]{7.08060519e+04, 7.81756099e+04, 2.71139000e+02, 6.95178335e+04,
         7.67533125e+04, 2.66206000e+02, 2.61143000e+02, 2.88323000e+02,
         1.00000000e+00};
        A[425]= new double[]{1.33999535e+05, 1.78904340e+05, 3.69447000e+02, 1.61449624e+05,
         2.15553273e+05, 4.45129000e+02, 3.62703000e+02, 4.84249000e+02,
         1.00000000e+00};
        A[426]= new double[]{8.59881572e+05, 3.30955684e+05, 9.04085000e+02, 3.12949394e+05,
         1.20449587e+05, 3.29037000e+02, 9.51107000e+02, 3.66067000e+02,
         1.00000000e+00};
        A[427]= new double[]{7.26013604e+05, 3.44922071e+05, 8.33129000e+02, 3.28211038e+05,
         1.55929903e+05, 3.76635000e+02, 8.71430000e+02, 4.14008000e+02,
         1.00000000e+00};
        A[428]= new double[]{3.81646147e+04, 1.07489265e+05, 2.01832000e+02, 9.11827057e+04,
         2.56812811e+05, 4.82216000e+02, 1.89091000e+02, 5.32568000e+02,
         1.00000000e+00};
        A[429]= new double[]{8.92173989e+04, 1.59658114e+05, 3.02674000e+02, 1.42062394e+05,
         2.54226352e+05, 4.81953000e+02, 2.94764000e+02, 5.27492000e+02,
         1.00000000e+00};
        A[430]= new double[]{2.17663023e+05, 1.90835605e+05, 4.70211000e+02, 1.74093479e+05,
         1.52636097e+05, 3.76089000e+02, 4.62905000e+02, 4.05851000e+02,
         1.00000000e+00};
        A[431]= new double[]{8.20588483e+05, 9.11241199e+04, 8.89441000e+02, 7.37748294e+04,
         8.19249422e+03, 7.99650000e+01, 9.22589000e+02, 1.02451000e+02,
         1.00000000e+00};
        A[432]= new double[]{1.71074058e+05, 4.71101964e+04, 4.18869000e+02, 4.32311511e+04,
         1.19049495e+04, 1.05850000e+02, 4.08419000e+02, 1.12470000e+02,
         1.00000000e+00};
        A[433]= new double[]{1.38489963e+05, 1.79834405e+05, 3.75156000e+02, 1.62844832e+05,
         2.11460115e+05, 4.41131000e+02, 3.69153000e+02, 4.79359000e+02,
         1.00000000e+00};
        A[434]= new double[]{2.79920101e+05, 1.46954634e+05, 5.32211000e+02, 1.34236849e+05,
         7.04727061e+04, 2.55224000e+02, 5.25957000e+02, 2.76121000e+02,
         1.00000000e+00};
        A[435]= new double[]{1.52249516e+05, 1.89852364e+05, 3.93434000e+02, 1.71996901e+05,
         2.14476992e+05, 4.44464000e+02, 3.86976000e+02, 4.82552000e+02,
         1.00000000e+00};
        A[436]= new double[]{5.24407500e+03, 3.40994500e+04, 8.50000000e+01, 2.23989250e+04,
         1.45648379e+05, 3.63059000e+02, 6.16950000e+01, 4.01170000e+02,
         1.00000000e+00};
        A[437]= new double[]{6.38810666e+03, 1.29583496e+04, 9.19880000e+01, 9.49556207e+03,
         1.92618595e+04, 1.36735000e+02, 6.94450000e+01, 1.40870000e+02,
         1.00000000e+00};
        A[438]= new double[]{7.30625144e+05, 2.91047119e+05, 8.35490000e+02, 2.74862632e+05,
         1.09492505e+05, 3.14313000e+02, 8.74487000e+02, 3.48355000e+02,
         1.00000000e+00};
        A[439]= new double[]{2.24693534e+05, 1.73704062e+05, 4.77868000e+02, 1.58386400e+05,
         1.22443938e+05, 3.36849000e+02, 4.70200000e+02, 3.63498000e+02,
         1.00000000e+00};
        A[440]= new double[]{3.80282417e+04, 1.96740172e+04, 2.06508000e+02, 1.75523461e+04,
         9.08075532e+03, 9.53160000e+01, 1.84149000e+02, 9.52700000e+01,
         1.00000000e+00};
        A[441]= new double[]{7.20862040e+05, 3.92667335e+05, 8.30469000e+02, 3.75673850e+05,
         2.04636729e+05, 4.32795000e+02, 8.68018000e+02, 4.72826000e+02,
         1.00000000e+00};
        A[442]= new double[]{6.41230118e+05, 2.55275262e+05, 7.84820000e+02, 2.40111192e+05,
         9.55888468e+04, 2.93879000e+02, 8.17041000e+02, 3.25266000e+02,
         1.00000000e+00};
        A[443]= new double[]{9.39756423e+04, 1.63601667e+05, 3.10123000e+02, 1.45990227e+05,
         2.54153565e+05, 4.81773000e+02, 3.03027000e+02, 5.27538000e+02,
         1.00000000e+00};
        A[444]= new double[]{7.30072370e+05, 5.92337747e+04, 8.37262000e+02, 4.50000654e+04,
         3.65104043e+03, 5.16070000e+01, 8.71976000e+02, 7.07470000e+01,
         1.00000000e+00};
        A[445]= new double[]{6.39326212e+04, 1.36774324e+05, 2.57771000e+02, 1.19696175e+05,
         2.56072143e+05, 4.82605000e+02, 2.48021000e+02, 5.30604000e+02,
         1.00000000e+00};
        A[446]= new double[]{2.17860492e+05, 1.09062443e+05, 4.69617000e+02, 9.97738027e+04,
         4.99474438e+04, 2.15071000e+02, 4.63911000e+02, 2.32237000e+02,
         1.00000000e+00};
        A[447]= new double[]{8.39690507e+04, 1.53435951e+05, 2.93893000e+02, 1.36231387e+05,
         2.48934486e+05, 4.76812000e+02, 2.85713000e+02, 5.22081000e+02,
         1.00000000e+00};
        A[448]= new double[]{1.61121493e+05, 1.92733377e+05, 4.04715000e+02, 1.74775108e+05,
         2.09065818e+05, 4.39011000e+02, 3.98111000e+02, 4.76220000e+02,
         1.00000000e+00};
        A[449]= new double[]{1.84789400e+05, 1.03497310e+05, 4.33112000e+02, 9.45045092e+04,
         5.29303220e+04, 2.21501000e+02, 4.26655000e+02, 2.38962000e+02,
         1.00000000e+00};
        A[450]= new double[]{1.81687747e+05, 1.85150197e+05, 4.29531000e+02, 1.68505233e+05,
         1.71716461e+05, 3.98366000e+02, 4.22991000e+02, 4.31052000e+02,
         1.00000000e+00};
        A[451]= new double[]{1.47994482e+04, 7.31254447e+04, 1.30979000e+02, 5.66314282e+04,
         2.79821134e+05, 5.01203000e+02, 1.12991000e+02, 5.58299000e+02,
         1.00000000e+00};
        A[452]= new double[]{1.62892710e+05, 1.63279888e+05, 4.06699000e+02, 1.48418173e+05,
         1.48770947e+05, 3.70560000e+02, 4.00524000e+02, 4.01476000e+02,
         1.00000000e+00};
        A[453]= new double[]{1.80181901e+05, 1.03117228e+05, 4.27954000e+02, 9.39981180e+04,
         5.37946672e+04, 2.23257000e+02, 4.21031000e+02, 2.40954000e+02,
         1.00000000e+00};
        A[454]= new double[]{1.83343602e+04, 7.52350731e+04, 1.44197000e+02, 5.97875326e+04,
         2.45338225e+05, 4.70220000e+02, 1.27148000e+02, 5.21752000e+02,
         1.00000000e+00};
        A[455]= new double[]{8.30410472e+05, 2.18079552e+05, 8.89119000e+02, 2.00480396e+05,
         5.26494745e+04, 2.14654000e+02, 9.33970000e+02, 2.45276000e+02,
         1.00000000e+00};
        A[456]= new double[]{7.22532520e+05, 3.80718246e+05, 8.31317000e+02, 3.63684648e+05,
         1.91633425e+05, 4.18441000e+02, 8.69142000e+02, 4.57970000e+02,
         1.00000000e+00};
        A[457]= new double[]{7.78969597e+05, 1.64208394e+05, 8.60992000e+02, 1.47941362e+05,
         3.11863437e+04, 1.63519000e+02, 9.04735000e+02, 1.90720000e+02,
         1.00000000e+00};
        A[458]= new double[]{7.27417884e+05, 3.75433310e+05, 8.34059000e+02, 3.58311691e+05,
         1.84931038e+05, 4.10841000e+02, 8.72142000e+02, 4.50128000e+02,
         1.00000000e+00};
        A[459]= new double[]{2.09897911e+04, 6.59818276e+04, 1.53368000e+02, 5.33447642e+04,
         1.67690332e+05, 3.89779000e+02, 1.36859000e+02, 4.30219000e+02,
         1.00000000e+00};
        A[460]= new double[]{2.32222878e+05, 1.76608445e+05, 4.85949000e+02, 1.60919628e+05,
         1.22381418e+05, 3.36740000e+02, 4.77875000e+02, 3.63430000e+02,
         1.00000000e+00};
        A[461]= new double[]{4.17080005e+05, 1.63415436e+05, 6.40310000e+02, 1.51449201e+05,
         5.93390642e+04, 2.32508000e+02, 6.51372000e+02, 2.55213000e+02,
         1.00000000e+00};
        A[462]= new double[]{1.89476464e+04, 6.72437406e+04, 1.46581000e+02, 5.35906569e+04,
         1.90189122e+05, 4.14583000e+02, 1.29264000e+02, 4.58748000e+02,
         1.00000000e+00};
        A[463]= new double[]{9.24193730e+05, 3.34448740e+05, 9.36173000e+02, 3.15930947e+05,
         1.14329609e+05, 3.20026000e+02, 9.87204000e+02, 3.57251000e+02,
         1.00000000e+00};
        A[464]= new double[]{2.04953919e+05, 1.84587307e+05, 4.56180000e+02, 1.68285238e+05,
         1.51562453e+05, 3.74564000e+02, 4.49283000e+02, 4.04637000e+02,
         1.00000000e+00};
        A[465]= new double[]{9.02779069e+05, 3.31508926e+05, 9.25744000e+02, 3.12885798e+05,
         1.14894595e+05, 3.20845000e+02, 9.75193000e+02, 3.58100000e+02,
         1.00000000e+00};
        A[466]= new double[]{2.11174649e+05, 2.96624178e+04, 4.67891000e+02, 2.68078262e+04,
         3.76553221e+03, 5.93970000e+01, 4.51333000e+02, 6.33960000e+01,
         1.00000000e+00};
        A[467]= new double[]{7.83923237e+05, 2.97634130e+05, 8.64660000e+02, 2.80390410e+05,
         1.06456540e+05, 3.09268000e+02, 9.06626000e+02, 3.44221000e+02,
         1.00000000e+00};
        A[468]= new double[]{8.98976036e+04, 1.22485033e+05, 3.03820000e+02, 1.09407768e+05,
         1.49067535e+05, 3.69757000e+02, 2.95891000e+02, 4.03150000e+02,
         1.00000000e+00};
        A[469]= new double[]{2.09536065e+05, 1.66120062e+05, 4.61460000e+02, 1.51240940e+05,
         1.19903723e+05, 3.33077000e+02, 4.54072000e+02, 3.59988000e+02,
         1.00000000e+00};
        A[470]= new double[]{6.58477144e+04, 1.11012500e+05, 2.61508000e+02, 9.75611690e+04,
         1.64478135e+05, 3.87455000e+02, 2.51800000e+02, 4.24509000e+02,
         1.00000000e+00};
        A[471]= new double[]{2.16577849e+05, 1.70588499e+05, 4.69035000e+02, 1.55537929e+05,
         1.22510136e+05, 3.36843000e+02, 4.61752000e+02, 3.63701000e+02,
         1.00000000e+00};
        A[472]= new double[]{8.01220334e+05, 1.98668991e+05, 8.73981000e+02, 1.81652699e+05,
         4.50422399e+04, 1.98149000e+02, 9.16748000e+02, 2.27315000e+02,
         1.00000000e+00};
        A[473]= new double[]{9.61772861e+04, 1.41047702e+05, 3.14234000e+02, 1.25871182e+05,
         1.84594946e+05, 4.11251000e+02, 3.06069000e+02, 4.48862000e+02,
         1.00000000e+00};
        A[474]= new double[]{1.07476302e+05, 1.75845897e+05, 3.31163000e+02, 1.57827695e+05,
         2.58227647e+05, 4.86309000e+02, 3.24542000e+02, 5.30995000e+02,
         1.00000000e+00};
        A[475]= new double[]{1.33552348e+05, 1.31080195e+05, 3.68868000e+02, 1.18551116e+05,
         1.16356647e+05, 3.27435000e+02, 3.62060000e+02, 3.55358000e+02,
         1.00000000e+00};
        A[476]= new double[]{5.75457770e+05, 3.22471840e+05, 7.45429000e+02, 3.06100127e+05,
         1.71530695e+05, 3.96512000e+02, 7.71982000e+02, 4.32599000e+02,
         1.00000000e+00};
        A[477]= new double[]{7.03297272e+05, 2.38918431e+05, 8.20611000e+02, 2.23355169e+05,
         7.58764020e+04, 2.60612000e+02, 8.57041000e+02, 2.91147000e+02,
         1.00000000e+00};
        A[478]= new double[]{7.89988304e+03, 2.56985293e+04, 1.01204000e+02, 1.83564325e+04,
         5.97139624e+04, 2.35161000e+02, 7.80590000e+01, 2.53928000e+02,
         1.00000000e+00};
        A[479]= new double[]{6.97121715e+05, 2.31497605e+05, 8.16936000e+02, 2.15912181e+05,
         7.16993198e+04, 2.53021000e+02, 8.53337000e+02, 2.83373000e+02,
         1.00000000e+00};
        A[480]= new double[]{5.60169690e+04, 9.61766981e+04, 2.41811000e+02, 8.40538314e+04,
         1.44313770e+05, 3.62839000e+02, 2.31656000e+02, 3.97735000e+02,
         1.00000000e+00};
        A[481]= new double[]{2.00724416e+05, 1.62572411e+05, 4.51204000e+02, 1.48301642e+05,
         1.20113716e+05, 3.33364000e+02, 4.44864000e+02, 3.60308000e+02,
         1.00000000e+00};
        A[482]= new double[]{2.25417820e+05, 2.27966360e+05, 4.79229000e+02, 2.07511076e+05,
         2.09857165e+05, 4.41160000e+02, 4.70376000e+02, 4.75694000e+02,
         1.00000000e+00};
        A[483]= new double[]{2.22055647e+05, 2.19953215e+05, 4.75771000e+02, 2.00250582e+05,
         1.98354601e+05, 4.29052000e+02, 4.66728000e+02, 4.62309000e+02,
         1.00000000e+00};
        
        return new Matrix(A);
    }

    private Matrix readMerton1UnnormalizedX1Data() {
        
        double[][] x1 = new double[3][];
        
        x1[0] = new double[]{676.587, 675.319, 908.129, 906.945, 907.284, 906.742, 906.447, 717.231,
        540.106, 659.06,  948.04,  84.722, 224.362, 708.462, 540.541, 842.513,
        145.319, 599.525, 202.223,  70.096, 580.71,  579.254, 225.174, 905.638,
        585.17,  674.36,  224.807, 432.705, 602.217, 605.064,  87.13,  517.801,
        537.234, 275.083,  74.254, 145.851, 685.191, 455.365,  61.619, 386.99,
         45.405,  60.262, 951.115, 812.927, 356.008,  61.382, 385.047,  58.727,
        948.143, 371.827,  47.343, 676.781,  87.316, 381.895, 586.808,  83.106,
        576.635, 290.638, 757.883, 517.551, 145.05,  365.203, 884.362, 514.23,
        405.142, 959.959, 834.807, 736.695, 240.613, 201.982, 766.748, 735.915,
        436.262, 847.984, 365.374, 951.49,  319.712, 800.284, 300.316, 690.206,
        734.901, 152.558, 934.547,  87.019, 341.699, 760.281,  49.619, 734.667,
        257.137, 547.383, 586.727, 945.511, 250.669,  55.487, 806.008, 878.711,
         75.616, 184.88,  361.185, 917.909, 731.193, 849.323,  46.783, 400.726,
         52.479, 416.626, 956.683, 157.032, 513.174, 188.643, 898.435, 864.98,
        508.839, 740.249, 198.505,  78.426, 148.033, 820.272, 346.482, 188.709,
        279.51,  207.889, 930.209, 556.144, 364.166, 556.133, 370.815, 329.578,
        444.973,  49.175,  87.012, 221.296, 657.421, 169.247,  95.897, 662.606,
         97.698, 492.969, 608.028, 529.102, 492.157,  76.689, 359.268, 254.002,
        330.505, 101.533,  58.259, 909.584, 117.494,  41.126, 327.038, 741.886,
        649.636, 651.878, 557.213, 356.425,  47.965, 396.004, 416.846,  73.8,
        785.982, 839.136, 926.039,  92.019, 809.528, 426.677, 361.882, 638.325,
        170.97,  876.309, 385.179, 328.886, 321.021, 360.791,  85.107, 799.126,
        366.684, 607.313, 564.811, 351.834, 557.808, 396.303, 480.374, 633.43,
        419.052, 725.82,  102.814, 182.984, 808.018,  70.46,  778.558, 365.06,
        395.967, 350.777, 664.863, 226.269, 429.939,  58.577, 925.226,  59.23,
         78.532, 605.011, 785.849, 343.261, 661.198, 139.568, 385.275, 329.386,
        420.63,  122.228, 708.172, 883.536, 347.494, 627.955, 432.454, 660.651,
         84.791, 477.694, 694.422,  51.041, 693.791, 831.923, 739.769, 119.842,
        327.012, 656.517, 400.562, 820.868, 284.271,  53.264, 583.893, 373.452,
        918.73,  773.312,  75.385,  75.078, 935.938, 932.456, 842.471, 939.636,
        202.274, 591.271, 202.067, 202.085, 509.703, 728.569, 725.024, 628.266,
        360.919, 675.735, 361.649, 213.034, 188.811,  86.989, 441.464, 288.2,
        101.808, 827.788, 117.283, 140.619, 260.083, 749.292, 816.107, 839.928,
        183.023, 129.535, 181.413, 423.551, 570.846, 305.571, 789.116, 604.038,
        509.424, 650.614, 434.22,  875.972, 849.74,  775.31,  886.118, 797.641,
        103.813, 421.285, 406.729, 681.781, 114.423, 488.8,  553.099, 672.559,
        542.127, 580.214, 106.953, 257.414, 659.264, 386.578, 602.578, 352.274,
        270.96,  842.163, 388.185, 896.311, 188.304, 435.921, 734.57,  411.728,
        787.28,  178.361, 273.905, 567.648, 107.254, 525.415, 741.027, 220.886,
        606.357, 173.569, 368.121, 350.203, 364.355, 780.185, 218.191, 503.212,
        113.892, 353.116, 614.506, 471.455, 521.201, 637.539, 353.311, 271.041,
        271.822, 920.188, 375.931, 738.577, 593.056, 102.509, 884.849, 174.405,
         86.67,  701.442, 638.662, 579.175, 873.47,  699.851, 165.902, 539.96,
        323.272, 318.493, 758.27,  84.925, 252.687, 633.17,  870.68,  702.401,
        834.409, 626.396, 816.21,  739.032, 407.527, 191.436, 102.468, 796.362,
        136.806,  97.263, 637.909, 402.625, 480.922, 170.815, 930.033, 340.466,
         97.294, 929.123, 814.732, 522.375, 217.593, 308.182, 371.654, 892.057,
        732.986, 102.131, 285.421, 804.602, 162.083, 834.02,  269.223, 890.903,
        782.452, 507.605, 564.928, 389.575, 535.319, 699.413, 283.356, 914.943,
        191.88,  268.208, 318.006, 238.758, 478.66,  627.834, 258.905, 396.782,
        115.148, 248.161, 153.813, 164.139, 181.35,  290.707, 751.938, 346.632,
        638.813, 584.515, 682.158, 564.991, 407.477, 571.866, 860.071, 226.657,
        147.625, 351.592, 756.399, 219.671, 491.125, 200.546, 310.878, 879.728,
        271.139, 369.447, 904.085, 833.129, 201.832, 302.674, 470.211, 889.441,
        418.869, 375.156, 532.211, 393.434,  85.,  91.988, 835.49,  477.868,
        206.508, 830.469, 784.82,  310.123, 837.262, 257.771, 469.617, 293.893,
        404.715, 433.112, 429.531, 130.979, 406.699, 427.954, 144.197, 889.119,
        831.317, 860.992, 834.059, 153.368, 485.949, 640.31,  146.581, 936.173,
        456.18,  925.744, 467.891, 864.66,  303.82,  461.46,  261.508, 469.035,
        873.981, 314.234, 331.163, 368.868, 745.429, 820.611, 101.204, 816.936,
        241.811, 451.204, 479.229, 475.771};
        x1[1] = new double[]{424.715, 394.656, 398.806, 384.96,  388.564, 373.851, 370.52,  234.166,
        403.087, 394.883, 431.453, 454.708, 451.941, 152.925, 426.223, 268.827,
         73.161, 425.61,  411.735, 458.704, 402.416, 179.536, 402.134, 360.022,
        429.141, 368.433, 439.057,  99.828, 397.06,  293.494, 427.859, 215.129,
        376.866,  53.568, 101.593,  31.91,  70.568, 135.899,  51.048, 240.371,
         94.986,  98.981, 106.089, 172.567, 489.063, 370.17,  120.743,  36.707,
        344.338, 238.154,  33.856, 184.452, 399.672, 268.496,  69.751, 272.282,
        372.968, 133.307, 284.658, 196.124, 138.536, 265.029,  30.946, 192.589,
        269.163,  86.162,  10.752, 357.861, 137.264, 388.113, 485.283, 304.146,
         65.877, 196.128, 250.199, 124.211, 178.227,  16.358, 146.504, 373.484,
        291.683, 129.284,  55.579,  82.428, 105.409, 170.201, 384.694, 163.371,
        107.538, 131.97,  75.404,  16.56,  368.175, 381.187, 188.095, 376.113,
         79.259, 210.281, 213.084, 103.941, 337.348, 162.403,  73.343, 143.853,
        452.505, 232.512, 332.228,  66.687, 233.748, 378.389, 125.087,  26.857,
        479.,   465.753, 119.443, 399.554, 162.975,  18.,   194.039, 383.778,
        294.95,  136.833, 152.549, 356.257, 456.384, 385.467,  66.156,  60.587,
        161.702, 489.948, 191.846,  87.84,  86.547, 233.735, 183.456, 469.025,
        361.292, 420.807, 219.725, 368.459, 374.606, 506.104, 400.037, 221.296,
        430.508, 248.709, 473.015,  76.452, 201.828, 278.571, 452.884,  50.46,
        213.133, 468.958, 100.314,  74.066, 352.361, 190.466, 372.025, 311.016,
        128.158, 437.041,  32.536, 213.629, 373.099, 382.378, 389.934,  66.028,
        350.163, 176.122, 414.182, 106.043,  86.285, 329.601, 225.016, 310.927,
         98.477,  35.064, 206.859, 170.187, 113.01,  127.093, 185.779, 387.93,
        479.726, 296.803, 503.642, 234.467, 358.837, 343.99,  335.063, 484.013,
         74.112, 337.314, 142.931, 494.623, 290.643, 505.721,  41.364, 343.97,
        288.902, 467.704, 372.396,  44.153, 144.687, 414.731,  70.302, 437.517,
        362.471, 504.187, 300.091, 122.386, 416.692,  27.205, 196.666,  58.379,
        247.392, 408.269,  46.065, 344.797,  19.566, 466.801, 397.329, 287.413,
         22.352, 354.818,  68.964, 231.259, 304.091, 118.431,  23.159,  51.07,
        180.569, 204.861, 248.673, 231.217, 372.702, 404.1,  125.557, 372.203,
        452.661, 178.454, 439.064, 425.419, 103.842, 204.589, 199.171, 140.483,
        202.193, 163.931, 196.568, 455.019, 455.259, 441.719, 178.248, 148.995,
        272.903, 144.388,  64.749, 146.338, 169.342, 170.656,  33.034, 142.271,
        279.791,  54.871, 367.564, 173.906, 186.028, 168.,   176.366, 239.05,
        134.684, 369.391, 419.295,  81.874, 151.283, 374.639, 348.265, 146.939,
         47.903,  44.954, 124.882, 356.706,  50.738, 179.748, 257.931, 196.202,
        249.321, 147.942, 370.498, 290.2,  178.795, 254.727, 361.056, 380.213,
        288.033, 103.262, 281.56,  146.144, 440.761, 456.548, 239.376,  77.725,
         85.938, 398.565, 161.999, 255.598, 388.628,  53.78,  419.195, 352.481,
        352.512, 143.299, 225.251, 388.759, 279.974, 390.335,  99.116,  47.166,
         88.474, 455.016, 242.549,  47.669, 493.853, 359.287, 302.097, 245.548,
        220.394, 137.263, 200.752, 439.524, 211.972, 196.369, 102.639, 398.988,
        390.903, 193.052, 438.534, 161.034, 338.121, 188.04,  488.105, 212.226,
        471.732, 166.53,  149.154, 201.91,  399.524, 204.347, 371.603, 213.048,
        364.69,  378.716, 420.164, 127.155, 323.337, 168.948,  85.183, 140.468,
        488.96,  466.206, 200.025, 473.641, 472.056, 202.871, 131.776, 476.542,
        378.281,  78.942, 372.206, 358.094, 134.245, 477.081, 137.415,  92.451,
        409.984, 226.383, 225.749, 433.147, 134.031, 298.826, 450.185, 154.769,
        120.8,  175.985, 363.928, 485.083, 194.706, 472.387, 432.694, 315.815,
         43.287, 484.262, 436.838, 378.763, 382.836, 352.255, 244.408, 467.263,
        503.269, 456.537, 417.031, 461.574, 188.886, 368.71,  78.048, 451.429,
        185.448, 349.959, 472.088, 437.948, 460.498, 216.472, 113.73,  353.833,
        500.656, 361.048, 382.169, 496.663, 404.703, 490.834, 189.911, 182.098,
        266.206, 445.129, 329.037, 376.635, 482.216, 481.953, 376.089,  79.965,
        105.85,  441.131, 255.224, 444.464, 363.059, 136.735, 314.313, 336.849,
         95.316, 432.795, 293.879, 481.773,  51.607, 482.605, 215.071, 476.812,
        439.011, 221.501, 398.366, 501.203, 370.56,  223.257, 470.22,  214.654,
        418.441, 163.519, 410.841, 389.779, 336.74,  232.508, 414.583, 320.026,
        374.564, 320.845,  59.397, 309.268, 369.757, 333.077, 387.455, 336.843,
        198.149, 411.251, 486.309, 327.435, 396.512, 260.612, 235.161, 253.021,
        362.839, 333.364, 441.16,  429.052};
        x1[2] = new double[]{1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
        
        return new Matrix(x1);
    }

    private Matrix readMerton1UnnormalizedX2Data() {
        
        double[][] x2 = new double[3][];
                
        x2[0]= new double[]{689.792,  688.621,  951.71,  951.48,  951.449,  951.327,  951.042,
        737.73,  532.714,  670.398,  997.464,  59.839,  211.447,  723.561,
        533.068,  879.181,  127.745,  602.231,  187.336,  43.995,  579.04,
        575.04,  212.198,  949.581,  586.413,  687.866,  211.482,  422.112,
        605.321,  609.026,  62.738,  509.628,  530.317,  262.811,  49.211,
        128.169,  699.719,  444.326,  34.623,  379.036,  16.872,  32.969,
        998.04,  850.672,  350.261,  35.683,  367.628,  31.136, 1000.549,
        363.74,   18.089,  690.931,  63.118,  374.088,  570.806,  58.923,
        575.359,  276.64,  792.633,  504.798,  127.145,  357.25,  908.632,
        500.618,  397.956, 1004.332,  848.826,  762.044,  225.115,  187.294,
        793.883,  760.647,  413.591,  888.451,  357.863, 1001.274,  313.462,
        810.551,  287.255,  708.361,  759.761,  133.515,  971.017,  63.754,
        321.822,  787.198,  21.979,  755.983,  241.403,  532.914,  570.07,
        976.181,  241.701,  29.286,  842.734,  922.666,  50.595,  169.906,
        355.5,   959.709,  755.36,  889.903,  17.528,  387.346,  24.607,
        409.058, 1010.854,  138.314,  504.405,  174.672,  940.01,  886.037,
        500.438,  765.559,  179.225,  52.449,  131.522,  833.319,  340.712,
        174.057,  270.687,  191.073,  981.225,  553.727,  357.687,  553.245,
        343.962,  302.912,  433.485,  21.301,  63.072,  199.86,  653.359,
        152.902,  73.243,  675.488,  76.463,  485.318,  616.643,  521.465,
        483.973,  51.853,  353.271,  244.31,  323.085,  79.232,  31.1,
        946.059,  97.748,  12.645,  319.778,  747.142,  662.603,  664.831,
        540.402,  330.923,  20.145,  388.989,  410.434,  49.197,  811.215,
        877.148,  957.153,  67.809,  844.215,  420.124,  355.484,  629.012,
        156.026,  922.391,  378.89,  309.034,  298.837,  354.069,  61.58,
        832.868,  344.948,  566.547,  562.152,  341.879,  543.074,  382.855,
        467.866,  642.877,  412.273,  749.879,  82.017,  168.911,  842.644,
        44.969,  809.134,  359.032,  369.566,  344.196,  671.477,  215.439,
        423.263,  30.736,  956.885,  32.465,  54.654,  609.658,  817.779,
        304.395,  667.801,  122.12,  359.007,  322.037,  413.541,  102.915,
        729.178,  922.434,  341.031,  589.704,  423.866,  654.831,  62.163,
        469.243,  691.688,  23.258,  660.858,  869.888,  765.453,  99.772,
        284.751,  670.148,  372.238,  857.616,  275.929,  42.062,  542.742,
        340.343,  970.518,  804.076,  50.613,  51.351,  983.941,  979.075,
        875.343,  987.509,  187.829,  588.715,  187.665,  187.2,   495.622,
        754.377,  750.605,  628.235,  355.617,  687.08,  356.718,  199.658,
        173.707,  62.535,  429.553,  274.683,  79.164,  861.463,  97.258,
        123.697,  252.085,  774.631,  848.55,  875.432,  167.555,  110.26,
        167.45,  411.207,  566.019,  295.211,  822.429,  609.635,  495.046,
        660.974,  427.681,  908.234,  888.164,  806.183,  930.195,  827.121,
        81.946,  395.439,  396.271,  699.036,  94.042,  474.886,  548.793,
        689.268,  537.07,  575.515,  86.809,  247.604,  669.633,  378.748,
        607.587,  345.861,  261.284,  874.983,  380.403,  941.349,  173.314,
        428.965,  759.784,  385.507,  816.091,  164.435,  262.843,  565.278,
        86.845,  510.902,  766.421,  209.092,  612.528,  157.207,  360.252,
        343.434,  357.056,  810.119,  196.83,  489.305,  93.88,  346.838,
        621.634,  455.579,  534.457,  647.61,  346.886,  261.235,  262.032,
        967.206,  370.37,  763.893,  597.325,  80.183,  920.485,  160.08,
        62.501,  722.947,  649.175,  574.624,  917.008,  719.856,  150.681,
        532.963,  316.336,  308.391,  781.593,  60.597,  243.269,  642.896,
        912.775,  723.185,  873.076,  635.174,  851.882,  755.868,  401.512,
        179.938,  80.676,  825.221,  119.036,  75.142,  648.761,  395.975,
        472.759,  154.356,  977.793,  333.572,  75.444,  967.934,  849.745,
        514.193,  201.235,  300.791,  355.395,  927.391,  758.034,  80.012,
        276.727,  838.077,  142.849,  872.772,  260.441,  936.385,  805.868,
        493.054,  562.873,  383.145,  525.452,  718.615,  274.953,  963.217,
        161.665,  259.617,  310.24,  227.489,  471.119,  637.044,  248.331,
        390.948,  95.327,  238.776,  137.345,  148.451,  167.805,  283.534,
        763.337,  339.765,  646.834,  586.107,  699.097,  562.857,  401.26,
        573.047,  894.045,  215.551,  130.617,  345.703,  784.329,  208.213,
        483.284,  187.844,  304.945,  924.8,   261.143,  362.703,  951.107,
        871.43,  189.091,  294.764,  462.905,  922.589,  408.419,  369.153,
        525.957,  386.976,  61.695,  69.445,  874.487,  470.2,   184.149,
        868.018,  817.041,  303.027,  871.976,  248.021,  463.911,  285.713,
        398.111,  426.655,  422.991,  112.991,  400.524,  421.031,  127.148,
        933.97,  869.142,  904.735,  872.142,  136.859,  477.875,  651.372,
        129.264,  987.204,  449.283,  975.193,  451.333,  906.626,  295.891,
        454.072,  251.8,   461.752,  916.748,  306.069,  324.542,  362.06,
        771.982,  857.041,  78.059,  853.337,  231.656,  444.864,  470.376,
        466.728};
        x2[1]= new double[]{460.223, 428.371, 438.48,  424.162, 428.005, 412.185, 409.085, 260.076,
        433.89,  428.393, 473.34,  506.364, 496.859, 172.184, 458.77,  300.105,
        68.88,  458.976, 452.946, 512.345, 434.035, 195.956, 440.976, 397.992,
        462.633, 400.741, 482.42,  106.342, 428.753, 318.621, 475.239, 231.852,
        405.928,  51.444,  99.568,  21.761,  84.279, 145.556,  40.431, 259.24,
        91.271,  96.19,  132.815, 198.268, 533.548, 410.314, 127.653,  23.446,
        383.014, 256.882,  19.859, 205.37,  443.279, 290.615,  78.214, 296.953,
        402.645, 140.667, 315.891, 211.441, 144.065, 286.542,  50.74,  207.339,
        291.147, 111.966,  27.369, 392.299, 143.451, 425.878, 526.24,  334.854,
        68.86,  224.155, 270.125, 152.342, 191.129,  31.428, 154.821, 406.926,
        321.669, 132.161,  78.21,  78.217, 109.876, 193.725, 427.025, 185.12,
        110.549, 143.43,  85.051,  38.273, 402.943, 423.702, 214.821, 414.418,
        73.653, 225.515, 229.375, 129.288, 369.959, 188.818,  65.543, 153.819,
        506.121, 250.729, 370.952,  62.393, 252.796, 416.153, 150.988,  45.303,
        515.967, 505.988, 122.568, 443.4,  171.162,  33.938, 208.955, 421.873,
        320.595, 142.461, 182.119, 384.772, 496.917, 415.778,  67.135,  59.898,
        174.428, 549.363, 203.135,  87.07,  99.007, 252.46,  194.443, 507.899,
        398.782, 453.514, 240.791, 396.768, 404.082, 566.752, 434.699, 238.807,
        470.064, 269.151, 529.234,  99.422, 215.967, 304.292, 494.616,  64.362,
        235.374, 507.664, 109.129,  75.182, 389.889, 204.97,  403.155, 341.451,
        149.206, 477.382,  54.187, 228.89,  409.834, 413.848, 424.872,  76.601,
        384.022, 204.611, 449.93,  110.023,  88.167, 357.229, 241.766, 344.054,
        102.521,  41.015, 225.002, 181.632, 123.623, 135.424, 200.091, 420.327,
        520.484, 326.327, 562.547, 253.007, 394.493, 378.656, 369.069, 527.523,
        76.275, 366.912, 159.79,  545.74,  314.399, 567.076,  63.289, 379.33,
        315.188, 505.081, 407.832,  41.978, 161.309, 458.699,  72.617, 477.415,
        392.551, 561.981, 329.921, 147.238, 454.188,  33.51,  212.49,  69.499,
        267.814, 440.481,  57.455, 380.594,  28.609, 509.032, 433.657, 314.284,
        17.55,  386.188,  69.958, 259.732, 330.479, 120.871,  28.075,  50.267,
        211.437, 230.097, 269.533, 250.221, 411.887, 444.753, 149.068, 411.706,
        498.778, 194.962, 483.062, 468.145, 112.227, 229.341, 223.343, 155.499,
        218.02,  183.078, 211.744, 500.556, 502.376, 491.477, 191.456, 157.259,
        296.767, 168.232,  58.402, 152.448, 180.796, 193.604,  50.388, 166.493,
        304.049,  47.4,  404.169, 186.495, 202.488, 178.779, 200.964, 260.938,
        145.492, 401.069, 454.241, 104.062, 176.947, 410.31,  385.466, 169.959,
        39.058,  45.559, 133.218, 388.879,  42.292, 192.975, 279.025, 218.151,
        270.166, 162.047, 409.329, 315.236, 198.672, 275.082, 391.365, 413.605,
        313.141, 125.016, 304.739, 173.05,  485.972, 494.259, 266.221,  80.525,
        105.306, 439.648, 172.058, 277.301, 430.097,  59.363, 456.642, 386.644,
        382.504, 149.267, 243.344, 423.076, 302.664, 426.296,  99.93,  51.747,
        85.505, 495.913, 264.976,  51.091, 538.037, 390.599, 327.621, 265.921,
        237.139, 164.493, 216.696, 477.996, 232.058, 208.964, 126.546, 439.983,
        432.833, 215.955, 474.696, 176.449, 374.623, 210.35,  540.915, 229.977,
        515.427, 177.679, 170.699, 215.552, 438.293, 225.138, 409.822, 237.151,
        401.644, 410.466, 459.523, 146.449, 350.134, 179.065,  81.753, 162.877,
        543.771, 519.539, 220.338, 514.099, 508.669, 217.469, 159.588, 520.201,
        418.462, 103.186, 408.864, 385.765, 139.788, 522.187, 145.585, 116.025,
        446.595, 243.739, 243.398, 472.897, 138.148, 332.196, 493.435, 182.334,
        141.363, 189.751, 393.007, 527.328, 210.533, 512.174, 474.183, 352.18,
        35.928, 531.903, 477.151, 415.353, 413.249, 382.601, 264.345, 507.679,
        561.144, 502.214, 460.969, 511.391, 201.489, 402.498,  94.367, 491.965,
        204.891, 378.99,  511.553, 472.318, 499.917, 235.868, 136.736, 387.277,
        556.693, 392.413, 417.946, 548.432, 436.707, 542.482, 204.04,  210.752,
        288.323, 484.249, 366.067, 414.008, 532.568, 527.492, 405.851, 102.451,
        112.47,  479.359, 276.121, 482.552, 401.17,  140.87,  348.355, 363.498,
        95.27,  472.826, 325.266, 527.538,  70.747, 530.604, 232.237, 522.081,
        476.22,  238.962, 431.052, 558.299, 401.476, 240.954, 521.752, 245.276,
        457.97,  190.72,  450.128, 430.219, 363.43,  255.213, 458.748, 357.251,
        404.637, 358.1,   63.396, 344.221, 403.15,  359.988, 424.509, 363.701,
        227.315, 448.862, 530.995, 355.358, 432.599, 291.147, 253.928, 283.373,
        397.735, 360.308, 475.694, 462.309};
        x2[2]= new double[]{1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
        
        return new Matrix(x2);
    }
   
    private PairFloatArray readMerton1UnnormalizedXY1Data() {
        
        Matrix m1 = readMerton1UnnormalizedX1Data();
        
        PairFloatArray xy = new PairFloatArray(m1.getColumnDimension());
        for (int i = 0; i < m1.getColumnDimension(); i++) {
            xy.add((float)m1.get(0, i), (float)m1.get(1, i));
        }
        
        return xy;
    }
    
    private PairFloatArray readMerton1UnnormalizedXY2Data() {
        
        Matrix m2 = readMerton1UnnormalizedX2Data();
        
        PairFloatArray xy = new PairFloatArray(m2.getColumnDimension());
        for (int i = 0; i < m2.getColumnDimension(); i++) {
            xy.add((float)m2.get(0, i), (float)m2.get(1, i));
        }
        
        return xy;
    }

    public static void main(String[] args) {
        
        try {
            StereoProjectionTransformerTest test = 
                new StereoProjectionTransformerTest();
            
            test.testB();
            
            //test.testC();
            
            //test.testE();
            
            //test.testF();
            
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }
    }

}
