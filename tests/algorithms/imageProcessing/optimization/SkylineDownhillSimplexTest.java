package algorithms.imageProcessing.optimization;

import algorithms.CountingSort;
import algorithms.MultiArrayMergeSort;
import algorithms.imageProcessing.CurvatureScaleSpaceCornerDetector;
import algorithms.imageProcessing.DFSContiguousValueFinder;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.SkylineExtractor;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SkylineDownhillSimplexTest extends TestCase {
    
    public SkylineDownhillSimplexTest() {
    }
    
    public void test0() throws Exception {
        
        boolean useBlueSkyImages = true;
        
        String[] fileNames;
        
        if (useBlueSkyImages) {
            fileNames = new String[]{
                "brown_lowe_2003_image1.jpg", 
                "venturi_mountain_j6_0001.png",
                /*"seattle.jpg",
                "arches.jpg",
                "stinson_beach.jpg",
                "cloudy_san_jose.jpg"*/
            };
        } else {
            fileNames = new String[]{
                 "stonehenge.jpg",
                 "norwegian_mtn_range.jpg",
                 "halfdome.jpg",
                 "costa_rica.jpg",
                 "new-mexico-sunrise_w725_h490.jpg",
                 "arizona-sunrise-1342919937GHz.jpg"
            };
        }
                
        List<ImageExt> images = new ArrayList<ImageExt>();
        List<GreyscaleImage> thetaImages = new ArrayList<GreyscaleImage>();
        List<Set<PairInt>> seedPoints = new ArrayList<Set<PairInt>>();
        List<Set<PairInt>> excludePoints = new ArrayList<Set<PairInt>>();
        List<Set<PairInt>> expectedSky = new ArrayList<Set<PairInt>>();
        List<String> fileNameRoots = new ArrayList<String>();
        
        for (String fileName : fileNames) {
            
            String filePath = ResourceFinder.findFileInTestResources(fileName);
            ImageExt img = ImageIOHelper.readImageExt(filePath);        
            int idx = fileName.lastIndexOf(".");
            String fileNameRoot = fileName.substring(0, idx);
            
            fileNameRoots.add(fileNameRoot);
            
            String filePath1 = ResourceFinder.findFileInTestResources(
            fileNameRoot + "_sky.png");
            GreyscaleImage skyMask = ImageIOHelper.readImageAsBinary(filePath1);
            // read the contiguous points above 0:
            Set<PairInt> skyPoints = readSkyPixels(skyMask);
            expectedSky.add(skyPoints);
            
            int width = img.getWidth();
            int height = img.getHeight();
            if (img.getWidth() != width || img.getHeight() != height) {
                throw new IllegalStateException(
                    "sky mask is not the same size a the test image");
            }
            CurvatureScaleSpaceCornerDetector detector = new
                CurvatureScaleSpaceCornerDetector(img);
            detector.useOutdoorModeAndExtractSkyline();
            detector.findCorners();
            
            images.add(img);
            thetaImages.add(detector.getTheta());
            
            SkylineExtractor skylineExtractor = new SkylineExtractor();
        
            SkylineExtractor.RemovedSets removedSets = 
                skylineExtractor.new RemovedSets();
            PairIntArray outputSkyCentroid = new PairIntArray();
        
            Set<PairInt> points = skylineExtractor.createBestSkyMaskPt1(
                detector.getTheta(), detector.getGradientXY(), 
                img, detector.getCannyEdgeFilterSettings(), outputSkyCentroid,
                removedSets);
        
            seedPoints.add(points);
            
            excludePoints.add(removedSets.getHighContrastRemoved());
            
        }
        
        SkylineANDedClauses skylineANDedClauses = new SkylineANDedClauses();
        ANDedClauses[] clauses;
         
        float[][] coeffLowerLimits;
        float[][] coeffUpperLimits;
        
        if (useBlueSkyImages) {
            clauses = skylineANDedClauses.getAllAndBlueClauses();
            coeffLowerLimits = skylineANDedClauses.getAllAndBlueCoeffLowerLimits();
            coeffUpperLimits = skylineANDedClauses.getAllAndBlueCoeffUpperLimits();
        } else {
            clauses = skylineANDedClauses.getAllAndRedClauses();
            coeffLowerLimits = skylineANDedClauses.getAllAndRedCoeffLowerLimits();
            coeffUpperLimits = skylineANDedClauses.getAllAndRedCoeffUpperLimits();
        }
        
        List<SetComparisonResults> resultsBeforeList = new ArrayList<SetComparisonResults>();
        
        // ---- get the comparison of points before refinement ----
        for (int i = 0; i < seedPoints.size(); i++) {
            Set<PairInt> points0 = new HashSet<PairInt>(seedPoints.get(i));
            Set<PairInt> exclude = excludePoints.get(i);
            ImageExt img = images.get(i);
            GreyscaleImage thetaImg = thetaImages.get(i);
         
            SkylineExtractor skylineExtractor = new SkylineExtractor();
            skylineExtractor.findClouds(points0, exclude, img, thetaImg,
                clauses);
              
            SetCompare setCompare = new SetCompare();
            SetComparisonResults results = setCompare.compare(
                expectedSky.get(i), points0);
                            
            resultsBeforeList.add(results);
            
            String fileNameRoot = fileNameRoots.get(i);
            
            try {
                String dirPath = ResourceFinder.findDirectory("bin");
                ImageExt clr = (ImageExt) img.copyImage();
                ImageIOHelper.addToImage(points0, 
                    thetaImg.getXRelativeOffset(),
                    thetaImg.getYRelativeOffset(), clr);
                ImageIOHelper.writeOutputImage(
                    dirPath + "/sky_before_optimization_" + fileNameRoot + ".png", clr);
            } catch (IOException e) {
                System.err.println("ERROR: " + e.getMessage());
            }
        }        
       
        // ---- get the comparison of points after refinement ----
        List<Set<PairInt>> finalSkyPoints = new ArrayList<Set<PairInt>>();
        finalSkyPoints.addAll(seedPoints);
        
        SkylineDownhillSimplex nelderMaed = new SkylineDownhillSimplex(images, 
            thetaImages, finalSkyPoints, excludePoints, expectedSky, clauses,
            coeffLowerLimits, coeffUpperLimits);
        
        SkylineFits fit = nelderMaed.fit();
        
        SetComparisonResults resultsBefore = new SetComparisonResults(
            resultsBeforeList);
        
        SetComparisonResults resultsAfter = fit.results;
        
        for (int i = 0; i < images.size(); i++) {
            
            Set<PairInt> points = finalSkyPoints.get(i);
            ImageExt img = images.get(i);
            GreyscaleImage thetaImg = thetaImages.get(i);
            
            String fileNameRoot = fileNameRoots.get(i);
            
            try {
                String dirPath = ResourceFinder.findDirectory("bin");
                ImageExt clr = (ImageExt) img.copyImage();
                ImageIOHelper.addToImage(points, 
                    thetaImg.getXRelativeOffset(),
                    thetaImg.getYRelativeOffset(), clr);
                ImageIOHelper.writeOutputImage(
                    dirPath + "/sky_after_optimization_" + fileNameRoot + ".png", clr);
            } catch (IOException e) {
                System.err.println("ERROR: " + e.getMessage());
            }
        }
        
        System.out.println("before: matched=" 
            + resultsBefore.numberMatchedDivExpected 
            + " overrun=" + resultsBefore.numberOverrunDivExpected);
        System.out.println("after : matched=" 
            + resultsAfter.numberMatchedDivExpected 
            + " overrun=" + resultsAfter.numberOverrunDivExpected);
        
        assertTrue(resultsAfter.numberMatchedDivExpected >= 
            resultsBefore.numberMatchedDivExpected);
        
    }
    
    private Set<PairInt> readSkyPixels(GreyscaleImage img) {
        
        DFSContiguousValueFinder zerosFinder = new DFSContiguousValueFinder(img);
        
        zerosFinder.findGroupsNotThisValue(0);
        
        int nGroups = zerosFinder.getNumberOfGroups();
        
        // ====== find the group(s) with the largest number of non-zero pixels
        int nMaxGroupN = Integer.MIN_VALUE;
        int[] groupIndexes = new int[nGroups];
        int[] groupN = new int[nGroups];
        for (int gId = 0; gId < nGroups; gId++) {
            int n = zerosFinder.getNumberofGroupMembers(gId);
            groupIndexes[gId] = gId;
            groupN[gId] = n;
            if (n > nMaxGroupN) {
                nMaxGroupN = n;
            }
        }
        
        int maxValue = MiscMath.findMax(groupN);
        if ((maxValue > groupN.length) || (nMaxGroupN > 10000000)) {
            MultiArrayMergeSort.sortByDecr(groupN, groupIndexes);
        } else {
            CountingSort.sortByDecr(groupN, groupIndexes, maxValue);
        }
        
        Set<PairInt> skyPoints = new HashSet<PairInt>();
        
        PairIntArray points = zerosFinder.getXY(groupIndexes[0]);
        
        for (int i = 0; i < points.getN(); i++) {
            skyPoints.add(new PairInt(points.getX(i), points.getY(i)));
        }
        
        return skyPoints;
    }

}
