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
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import org.junit.After;
import static org.junit.Assert.assertTrue;
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
        Matrix unnormXY1 = DataForTests.readMerton1UnnormalizedX1Data();
        
        // mRows = 3;  nCols = 484
        Matrix unnormXY2 = DataForTests.readMerton1UnnormalizedX2Data();
        
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
                
       
        PairFloatArray matched1 = DataForTests.readMerton1UnnormalizedXY1Data();
        
        PairFloatArray matched2 = DataForTests.readMerton1UnnormalizedXY2Data();
       
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
        
        assertTrue(fit.getNMatches() == unnormXY1.getColumnDimension());
    }
    
    @Test
    public void testG() throws Exception {
                
        PairFloatArray corners1 = new PairFloatArray();
        PairFloatArray corners2 = new PairFloatArray();
        DataForTests.readBrownAndLoweCorners(corners1, corners2);
        
        PairIntArray corn1 = new PairIntArray();
        PairIntArray corn2 = new PairIntArray();
        DataForTests.readBrownAndLoweCorners(corn1, corn2);
        
        String fileName1 = "brown_lowe_2003_image1.jpg";
        String fileName2 = "brown_lowe_2003_image2.jpg";
    
        PointMatcher pointMatcher = new PointMatcher();
        
        PairIntArray outputMatched1 = new PairIntArray();
        PairIntArray outputMatched2 = new PairIntArray();
        
        // side effect is corn1 and corn2 become matched
        TransformationPointFit fit2 = pointMatcher.calculateTransformationAndMatch(
            corn1, corn2, 517, 374, outputMatched1, outputMatched2);
        
        log.info(fit2.toString());
        
        double[] diffFromModel = new double[corn1.getN()];
        double[] avgDiffModel = new double[1];
        
        int centroidX1 = 517 >> 1;
        int centroidY1 = 374 >> 1;
        double transXTol = 517 * 0.02;
        double transYTol = 374 * 0.02;
        double tolerance = fit2.getMeanDistFromModel();
        //7.5 is mean dist from model;  10 is mean dist + stdev;
        
        pointMatcher.populateDiffFromModelWithMatchedPoints(
            outputMatched1, outputMatched2, 
            fit2.getParameters(), centroidX1, centroidY1, 
            diffFromModel, avgDiffModel);
        
        tolerance = avgDiffModel[0] + 0.5*fit2.getStDevFromMean();
        
        for (int i = (outputMatched1.getN() - 1); i > -1; i --) {
            double diff = diffFromModel[i];
            if (diff > tolerance) {
                outputMatched1.removeRange(i, i);
                outputMatched2.removeRange(i, i);
            }
        }
        int nIterMax = 15;
        int nIter = 0;
        
        log.info("Point matcher matched " + corn1.getN() + " from " 
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
        /*
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
         */
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
                    
        Matrix theLeftPoints = spTransformer.rewriteInto3ColumnMatrix(matched1);
        Matrix theRightEpipolarLines = spTransformer.calculateEpipolarRightLines(
            theLeftPoints);
        
        Matrix theRightPoints = spTransformer.rewriteInto3ColumnMatrix(matched2);
        Matrix theLeftEpipolarLines = spTransformer.calculateEpipolarLeftLines(
            theRightPoints);
        
        Color clr = null;
        for (int i = 0; i < theLeftEpipolarLines.getColumnDimension(); i++) {
            
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
        
        for (int i = 0; i < theRightEpipolarLines.getColumnDimension(); i++) {
            
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
            
            test.testG();
            
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }
    }

}
