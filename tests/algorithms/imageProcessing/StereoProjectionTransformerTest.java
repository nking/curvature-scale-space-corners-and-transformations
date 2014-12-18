package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import algorithms.util.LinearRegression;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import java.awt.Color;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import org.junit.After;
import static org.junit.Assert.assertTrue;
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
    
    //@Test
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
    
    public void estE() throws Exception {
        
        PairFloatArray matched1 = new PairFloatArray();
        PairFloatArray matched2 = new PairFloatArray();
        
        DataForTests.readBrownAndLoweMatches(matched1, matched2);
        
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
            leftMatches.add(Math.round(matched1.getX(i)),
                Math.round(matched1.getY(i)));
            rightMatches.add(Math.round(matched2.getX(i)),
                Math.round(matched2.getY(i)));
        }
        
        Color clr = null;
        PairIntArray subsetLeft = new PairIntArray();
        for (int i = 0; i < nLimitTo; i++) {
            
            clr = getColor(clr);
            
            PairIntArray leftLine = spTransformer.getEpipolarLineInLeft(
                img1.getWidth(), img1.getHeight(), i);
            
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
                img2.getWidth(), img2.getHeight(), i);
            
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
    
    //@Test
    public void estF() throws Exception {
        
        /*
        test the fundamental matrix using the Merton I College set from
        http://www.robots.ox.ac.uk/~vgg/data/data-mview.html
        
        Note that this set was chosen to compare initial results with those from
        the "Programming Computer Vision with Python" by Jan Solem.
        */
        
        // mRows = 3;  nCols = 484
        SimpleMatrix unnormXY1 = DataForTests.readMerton1UnnormalizedX1Data();
        
        // mRows = 3;  nCols = 484
        SimpleMatrix unnormXY2 = DataForTests.readMerton1UnnormalizedX2Data();
        
        StereoProjectionTransformer spTransformer = 
            new StereoProjectionTransformer();
        
        PairFloatArray matched1 = DataForTests.readMerton1UnnormalizedXY1Data();
        
        PairFloatArray matched2 = DataForTests.readMerton1UnnormalizedXY2Data();
       
        spTransformer = 
            new StereoProjectionTransformer();
        
        spTransformer.calculateEpipolarProjection(matched1, matched2);
        //spTransformer.calculateEpipolarProjectionWithoutNormalization(matched1, matched2);
        
        double[] leftEpipole = spTransformer.getLeftEpipole();
        
        double[] rightEpipole = spTransformer.getRightEpipole();
    
        int nPoints = spTransformer.getEpipolarLinesInLeft().numCols();
        
        String fileName1 = "merton_college_I_001.jpg";
        String fileName2 = "merton_college_I_002.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        
        int nLimitTo = nPoints;
        
        Color clr = null;
        PairIntArray subsetLeft = new PairIntArray();
        for (int i = 0; i < nLimitTo; i++) {
            
            clr = getColor(clr);
            
            PairIntArray leftLine = spTransformer.getEpipolarLineInLeft(
                img1.getWidth(), img1.getHeight(), i);
            
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
                img2.getWidth(), img2.getHeight(), i);
            
            ImageIOHelper.addCurveToImage(rightLine, img2, 0, 
                clr.getRed(), clr.getGreen(), clr.getBlue()); 
            
            subsetRight.add(Math.round(matched2.getX(i)), Math.round(matched2.getY(i)));
        }
        ImageIOHelper.addCurveToImage(subsetRight, img2, 2, 255, 0, 0);
    
        String dirPath = ResourceFinder.findDirectory("bin");
        
        ImageIOHelper.writeOutputImage(dirPath + "/tmp1.png", img1);
       
        ImageIOHelper.writeOutputImage(dirPath + "/tmp2.png", img2);
        
        StereoProjectionTransformerFit fit = spTransformer.evaluateFitForRight(
            3);
        
        assertTrue(fit.getNMatches() == unnormXY1.numCols());
    }
    
    //@Test
    public void estG() throws Exception {
        
        // ===== learning steps to matching points between images with 
        // projections such as epipolar projections.
        // the epipolar projection is sensitive to precisely matched input
        // points so this solves that stage and compares total to known
        // total result
        /*
        Best results so far from these steps:
        
        (1) extract corners in outdoor mode.
        (2) get a rough euclidean transformation from the PointMatcher.
        (3) use that to transform the image1 edges from outdoor mode 
        (4) mark edges from both sets that are not completely internal              
            to their union plus a buffer in skip lists
        (5)
        */
        
        /*
        PairIntArray corners1Int = new PairIntArray();
        PairIntArray corners2Int = new PairIntArray();        
        DataForTests.readBrownAndLoweCorners(corners1Int, corners2Int);
        */
        /*
        String fileName1 = "merton_college_I_001.jpg";
        String fileName2 = "merton_college_I_002.jpg";
        DataForTests.readUnmatchedMerton1CornersFromThisCode(corn1, corn2);
        */
        
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
        
        /*
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
        */
        
        PairIntArray points1 = new PairIntArray();
        PairIntArray points2 = new PairIntArray();
        DataForTests.readBrownAndLoweInflectionPointsImage1(points1);        
        DataForTests.readBrownAndLoweInflectionPointsImage2(points2);
        
        PointMatcher pointMatcher = new PointMatcher();
        
        TransformationPointFit fit = 
            pointMatcher.calculateRoughTransformationForUnmatched(
            points1, points2, image1Width, image1Height);
        
        log.info("FIT: " + fit.toString());
        
        if (true) {
            return;
        }
     
        // === get rough euclidean projection ====
        
        /*
        // ==== transform edges1 and corners1Int to image2 frame ====
        Transformer transformer = new Transformer();
        
        PairIntArray corners1Transformed = transformer.applyTransformation(
            fit2.getParameters(), corners1Int, 
            image1Width >> 1, image1Height >> 1);
           
        List<PairIntArray> edges1Transformed = transformer.applyTransformation(
            fit2.getParameters(), edges1,
            image1Width >> 1, image1Height >> 1);

        int buffer = 50;
        
        ExternalEdgeFinder extEF = new ExternalEdgeFinder();
        Set<Integer> edges2SkipList = extEF.findExteriorEdges(edges2, 
            edges1Transformed, buffer);
        
        Set<Integer> edges1SkipList = extEF.findExteriorEdges(edges1Transformed,
            image2Width, image2Height);
       */
        
        // edges are now roughly in same reference frame but there may be
        // position dependent increasing differences in 
        // translation and rotation due to projection effects
        // so registration:
        //     -- calc centroid of each edge not in a skip list
        //
        //
        
        /*
        double[] diffFromModel = new double[corners1Int.getN()];
        double[] avgDiffModel = new double[1];
        
        int centroidX1 = 517 >> 1;
        int centroidY1 = 374 >> 1;
        double transXTol = 517 * 0.02;
        double transYTol = 374 * 0.02;
        double tolerance = avgDiffModel[0] + 0.5*fit2.getStDevFromMean();
        
        int nIterMax = 15;
        int nIter = 0;
        
        log.info("Point matcher matched " + corners1Int.getN() + " from " 
            + outputMatched1.getN());
        
        PairFloatArray matched1 = new PairFloatArray();
        PairFloatArray matched2 = new PairFloatArray();
        for (int i = 0; i < outputMatched1.getN(); i++) {
            matched1.add(outputMatched1.getX(i), outputMatched1.getY(i));
            matched2.add(outputMatched2.getX(i), outputMatched2.getY(i));
        }
        
        // need the same number of corners, so truncating one without
        // using criteria for which to remove
                
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(System.currentTimeMillis());
                
        int nSubset = 16;
        int[] selected = new int[nSubset];
       
        int nKeep = 20;
        
        tolerance = Double.MAX_VALUE;
        
        nIterMax = 5;
        nIter = 0;
        
        double minRemoval = 5;
                                
        while (nIter < nIterMax) {
            
            StereoProjectionTransformer spTransformer = solveAndPlot(
                matched1, matched2, fileName1, fileName2, "ALL",
                matched1, matched2);

            StereoProjectionTransformerFit fitOfAllMatches = 
                spTransformer.evaluateFitInRightImageForMatchedPoints(
                matched1, matched2, tolerance);

            log.info("stereo projection FIT to all remaining corners) " 
                + fitOfAllMatches.toString());
            
 if (true) {
     break;
 }
 
 
            // keep the top fits and look at the image projections
            List<FitHolder> best = new ArrayList<FitHolder>();
            FitHolder ph = new FitHolder();
            ph.fit = fitOfAllMatches;
            ph.label = Long.MIN_VALUE; // not a subset bitstring, but used for label
            ph.setSelected(selected);

            best.add(ph);
            
            int nPoints = matched1.getN();
        
            Set<String> chosen = new HashSet<String>();
        
            for (int i = 0; i < 100; i++) {

                chooseRandomly(sr, selected, chosen, nPoints);

                PairFloatArray subset1 = new PairFloatArray();
                PairFloatArray subset2 = new PairFloatArray();
                for (int idx : selected) {
                    subset1.add(matched1.getX(idx), matched1.getY(idx));
                    subset2.add(matched2.getX(idx), matched2.getY(idx));
                }

                StereoProjectionTransformer spTransformer2 = new 
                    StereoProjectionTransformer();

                spTransformer2.calculateEpipolarProjection(subset1, subset2);

                StereoProjectionTransformerFit fit = 
                    spTransformer2.evaluateFitInRightImageForMatchedPoints(
                    matched1, matched2, tolerance);

                // look for idx=33 (471,156) and idx=36 (484,165)
                if ((Arrays.binarySearch(selected, 33) > -1) && 
                    (Arrays.binarySearch(selected, 36) > -1)) {

                    log.info("i=" + i + " " + Arrays.toString(selected));

                    log.info(i + ") " + fit.toString());
                }

                log.fine(i + ") " + fit.toString());

                if (!Double.isNaN(fit.getAvgDistance()) &&
                    fit.getAvgDistance() >= tolerance) {
                    continue;
                }

                ph = new FitHolder();
                ph.fit = fit;
                ph.label = (long)i; // not a subset bitstring, but used for label
                ph.setSelected(selected);

                best.add(ph);

                // trim best to size nKeep
                if ((i > 0) && ((i % 50) == 0) && (best.size() > nKeep)) {
                    Collections.sort(best, new FitComparator());
                    while (best.size() > nKeep) {
                        best.remove(best.size() - 1);
                    }
                }
            }

            Collections.sort(best, new FitComparator());
            while (best.size() > nKeep) {
                best.remove(best.size() - 1);
            }
         
            // ====== plot the best fitting projections =====
            
            for (FitHolder fh : best) {

                StereoProjectionTransformerFit fit = fh.fit;

                PairFloatArray subset1 = new PairFloatArray();
                PairFloatArray subset2 = new PairFloatArray();
                for (int idx : fh.getSelected()) {
                    subset1.add(matched1.getX(idx), matched1.getY(idx));
                    subset2.add(matched2.getX(idx), matched2.getY(idx));
                }

                spTransformer = solveAndPlot(
                   subset1, subset2, fileName1, fileName2, 
                    "_" + fh.label + "_", matched1, matched2);

                log.info(fh.label + ") " +  fit.toString());
            }

            // ==== remove outliers using best fit and iterate if any were removed =====
            
            FitHolder bestFitHolder = best.get(0);
        
            PairFloatArray subset1 = new PairFloatArray();
            PairFloatArray subset2 = new PairFloatArray();
            for (int idx : bestFitHolder.getSelected()) {
                subset1.add(matched1.getX(idx), matched1.getY(idx));
                subset2.add(matched2.getX(idx), matched2.getY(idx));
            }            
            
            LinkedHashSet<Integer> skipIndexes = new LinkedHashSet<Integer>();
            
            StereoProjectionTransformer spTransformer2 = new 
                StereoProjectionTransformer();

            spTransformer2.calculateEpipolarProjection(subset1, subset2);
                        
            float factor = 2.f;

            StereoProjectionTransformerFit fit = 
                spTransformer2.evaluateRightForMatchedAndStoreOutliers(
                matched1, matched2, factor, minRemoval, skipIndexes);
            
            if (skipIndexes.isEmpty()) {
                break;
            }
            
            log.info("removing " + skipIndexes + " points");
            
            List<Integer> rm = new ArrayList<Integer>(skipIndexes);
            
            for (int i = (rm.size() - 1); i > -1; i--) {
                int idx = rm.get(i).intValue();
                matched1.removeRange(idx, idx);
                matched2.removeRange(idx, idx);
            }
                
            nIter++;
        }
        */
    }
    
    @Test
    public void testH() throws Exception {
                    
        String fileName1 = "merton_college_I_001.jpg";
        String fileName2 = "merton_college_I_002.jpg";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        
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
                
        //8th point and 9th to check these w/ existing method
        //xy1.add(271.139000f, 266.206000f);
        //xy2.add(261.143000f, 288.323000f);
        //xy1.add(101.204000f, 235.161000f);  
        //xy2.add(78.059000f, 253.928000f);
        
        // mRows = 3;  nCols = 484
        SimpleMatrix unnormXY1 = DataForTests.readMerton1UnnormalizedX1Data();
        // mRows = 3;  nCols = 484
        SimpleMatrix unnormXY2 = DataForTests.readMerton1UnnormalizedX2Data();
        xy1 = new PairFloatArray();
        xy2 = new PairFloatArray();
        for (int j = 0; j < unnormXY1.numCols(); j++) {
            xy1.add((float) unnormXY1.get(0, j), (float) unnormXY1.get(1, j));
            xy2.add((float) unnormXY2.get(0, j), (float) unnormXY2.get(1, j));
        }
        /*
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(System.currentTimeMillis());
        Set<String> chosen = new HashSet<String>();
        int[] selected = new int[7];
        */
        for (int i = 0; i < 1; i++) {
            
            StereoProjectionTransformer spTransformer = new StereoProjectionTransformer();
/*
            chooseRandomly(sr, selected, chosen, 484);
            xy1 = new PairFloatArray();
            xy2 = new PairFloatArray();
            for (int j = 0; j < selected.length; j++) {
                int idx = selected[j];
                xy1.add((float)unnormXY1.get(0, idx), (float)unnormXY1.get(1, idx));
                xy2.add((float)unnormXY2.get(0, idx), (float)unnormXY2.get(1, idx));
            }
 */           
            //SimpleMatrix fm = spTransformer.calculateEpipolarProjectionFor7Points(xy1, xy2);
            SimpleMatrix fm = spTransformer.calculateEpipolarProjection(xy1, xy2);
            
            SimpleMatrix input1 = spTransformer.rewriteInto3ColumnMatrix(xy1);
            SimpleMatrix input2 = spTransformer.rewriteInto3ColumnMatrix(xy2);

            for (int row = 0; row < fm.numRows(); row++) {
                StringBuilder sb = new StringBuilder();
                for (int col = 0; col < fm.numCols(); col++) {
                    sb.append(fm.get(row, col)).append(", ");
                }
                System.out.println(row + ") " + sb.toString());
            }
            
            Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
            int image1Width = img1.getWidth();
            int image1Height = img1.getHeight();
            Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);

            for (int ii = 0; ii < input1.numCols(); ii++) {
                double x = input1.get(0, ii);
                double y = input1.get(1, ii);            
                ImageIOHelper.addPointToImage((float)x, (float)y, img1, 3, 255, 0, 0);
                double x2 = input2.get(0, ii);
                double y2 = input2.get(1, ii);            
                ImageIOHelper.addPointToImage((float)x2, (float)y2, img2, 3, 255, 0, 0);
                //System.out.println(String.format("%f, %f  %f, %f", x, y, x2, y2));
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
                dirPath + "/tmp_m_1_" + i + ".png", img1);
            ImageIOHelper.writeOutputImage(
                dirPath + "/tmp_m_2_" + i + ".png", img2);
        }
    }
    
    private void chooseRandomly(SecureRandom sr, int[] selected, 
        Set<String> chosen, int nPoints) {
        
        while (true) {
            for (int i = 0; i < selected.length; i++) {
                int sel = sr.nextInt(nPoints);
                while (contains(selected, i, sel)) {
                    sel = sr.nextInt(nPoints);
                }
                selected[i] = sel;
            }
            Arrays.sort(selected);
            String str = Arrays.toString(selected);
            if (!chosen.contains(str)) {
                chosen.add(str);
                break;
            }
        }
    }

    private StereoProjectionTransformer solveAndPlot(
        PairFloatArray matchedSubset1, PairFloatArray matchedSubset2, 
        String imageFileName1, 
        String imageFileName2, String fileNumber,
        PairFloatArray matched1, PairFloatArray matched2) throws Exception {
        
        StereoProjectionTransformer spTransformer 
            = new StereoProjectionTransformer();
        
        spTransformer.calculateEpipolarProjection(matchedSubset1, matchedSubset2);
                
        String fileName1 = imageFileName1;
        String fileName2 = imageFileName2;
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
                    
        SimpleMatrix theLeftPoints = spTransformer.rewriteInto3ColumnMatrix(matched1);
        SimpleMatrix theRightEpipolarLines = spTransformer.calculateEpipolarRightLines(
            theLeftPoints);
        
        SimpleMatrix theRightPoints = spTransformer.rewriteInto3ColumnMatrix(matched2);
        SimpleMatrix theLeftEpipolarLines = spTransformer.calculateEpipolarLeftLines(
            theRightPoints);
        
        Color clr = null;
        for (int i = 0; i < theLeftEpipolarLines.numCols(); i++) {
            
            clr = getColor(clr);
            
            PairIntArray leftLine = spTransformer.getEpipolarLine(
                theLeftEpipolarLines, img1.getWidth(), img1.getHeight(), i);
            
            ImageIOHelper.addCurveToImage(leftLine, img1, 0, 
                clr.getRed(), clr.getGreen(), clr.getBlue());
            
            ImageIOHelper.addPointToImage(matched2.getX(i), matched2.getY(i), 
                img2, 3, 255, 0, 0);
            
            ImageIOHelper.addPointToImage(matched2.getX(i), matched2.getY(i), 
                img2, 2, 
                clr.getRed(), clr.getGreen(), clr.getBlue());
        }
        
        clr = null;
        
        for (int i = 0; i < theRightEpipolarLines.numCols(); i++) {
            
            clr = getColor(clr);
            
            PairIntArray rightLine = spTransformer.getEpipolarLine(
                theRightEpipolarLines, img1.getWidth(), img1.getHeight(), i);
            
            ImageIOHelper.addCurveToImage(rightLine, img2, 0, 
                clr.getRed(), clr.getGreen(), clr.getBlue()); 
            
            ImageIOHelper.addPointToImage(matched1.getX(i), matched1.getY(i), 
                img1, 3, 255, 0, 0);
            
            ImageIOHelper.addPointToImage(matched1.getX(i), matched1.getY(i), 
                img1, 2, clr.getRed(), clr.getGreen(), clr.getBlue());
        }
    
        String dirPath = ResourceFinder.findDirectory("bin");
        
        ImageIOHelper.writeOutputImage(
            dirPath + "/tmp_" + fileNumber + "_1.png", img1);
       
        ImageIOHelper.writeOutputImage(
            dirPath + "/tmp_" + fileNumber + "_2.png", img2);
        
        return spTransformer;
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
    
    private boolean contains(int[] values, int lastIdx, int valueToCheck) {
        for (int i = 0; i < lastIdx; i++) {
            if (values[i] == valueToCheck) {
                return true;
            }
        }
        return false;
    }
    
    public static class FitHolder {
        public StereoProjectionTransformerFit fit;
        public Long label;
        private int[] selected = null;
        public void setSelected(int[] theSelectedIndexes) {
            selected = Arrays.copyOf(theSelectedIndexes, theSelectedIndexes.length);
        }
        public int[] getSelected() { return selected;}
    }
    public class FitComparator implements Comparator<FitHolder> {

        @Override
        public int compare(FitHolder o1, FitHolder o2) {
            
            if (o2 == null && o1 != null) {
                return -1;
            } else if (o1 == null && o2 != null) {
                return 1;
            } else if (o1 == null && o2 == null) {
                return 0;
            }
            
            if (o2.fit == null && o1.fit != null) {
                return -1;
            } else if (o1.fit == null && o2.fit != null) {
                return 1;
            } else if (o1.fit == null && o2.fit == null) {
                return 0;
            }
            
            if (o1.fit.getNMatches() > o2.fit.getNMatches()) {
                return -1;
            } else if (o1.fit.getNMatches() < o2.fit.getNMatches()) {
                return 1;
            }
            
            if (Double.isNaN(o2.fit.getAvgDistance()) && 
                !Double.isNaN(o1.fit.getAvgDistance())) {
                return -1;
            } else if (Double.isNaN(o1.fit.getAvgDistance()) && 
                !Double.isNaN(o2.fit.getAvgDistance())) {
                return 1;
            }
            
            if (o1.fit.getAvgDistance() < o2.fit.getAvgDistance()) {
                return -1;
            } else if (o1.fit.getAvgDistance() > o2.fit.getAvgDistance()) {
                return 1;
            } else {
                if (o1.fit.getStDevFromAvg() < o2.fit.getStDevFromAvg()) {
                    return -1;
                } else if (o1.fit.getStDevFromAvg() > o2.fit.getStDevFromAvg()) {
                    return 1;
                }
            }
            
            return Long.compare(o1.fit.getNMatches(), o2.fit.getNMatches());
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
            
            //test.testG();
            
            test.testH();
                        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }
    }

}
