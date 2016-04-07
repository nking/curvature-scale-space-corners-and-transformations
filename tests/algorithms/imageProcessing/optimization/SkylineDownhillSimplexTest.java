package algorithms.imageProcessing.optimization;

import algorithms.CountingSort;
import algorithms.MultiArrayMergeSort;
import algorithms.imageProcessing.CannyEdgeFilterSettings;
import algorithms.imageProcessing.scaleSpace.CurvatureScaleSpaceCornerDetector;
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
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
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
        
        if (true) {
            // re-enable upon need.  very time consuming!
            return;
        }
        
        boolean useAllImages = true;
        boolean useBlueSkyImages = false;
        
        String[] fileNames;
        
        if (useAllImages) {
            fileNames = new String[]{
                "brown_lowe_2003_image1.jpg",
                "venturi_mountain_j6_0001.png",
                "seattle.jpg",
                "arches.jpg",
                "stinson_beach.jpg",
                "cloudy_san_jose.jpg",
                "stonehenge.jpg",
                "norwegian_mtn_range.jpg",
                "halfdome.jpg",
                "costa_rica.jpg",
                "new-mexico-sunrise_w725_h490.jpg",
                "arizona-sunrise-1342919937GHz.jpg"
            };
        } else {
            if (useBlueSkyImages) {
                fileNames = new String[]{
                    "brown_lowe_2003_image1.jpg",
                    "venturi_mountain_j6_0001.png",
                    "seattle.jpg",
                    "arches.jpg",
                    "stinson_beach.jpg",
                    "cloudy_san_jose.jpg"};
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
        }
                
        List<ImageExt> images = new ArrayList<ImageExt>();
        List<GreyscaleImage> thetaImages = new ArrayList<GreyscaleImage>();
        List<Set<PairInt>> seedPoints = new ArrayList<Set<PairInt>>();
        List<Set<PairInt>> excludePoints = new ArrayList<Set<PairInt>>();
        List<Set<PairInt>> expectedSky = new ArrayList<Set<PairInt>>();
        List<Set<PairInt>> expectedBorderPoints = new ArrayList<Set<PairInt>>();
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
            //detector.useOutdoorModeAndExtractSkyline();
            detector.findCorners();
            
            images.add(img);
            thetaImages.add(detector.getEdgeFilterProducts().getTheta());
            
            SkylineExtractor skylineExtractor = new SkylineExtractor();
        
            SkylineExtractor.RemovedSets removedSets = 
                skylineExtractor.new RemovedSets();
            PairIntArray outputSkyCentroid = new PairIntArray();
            
            CannyEdgeFilterSettings settings = new CannyEdgeFilterSettings();
        
            Set<PairInt> points = skylineExtractor.extractSkyStarterPoints(
                detector.getEdgeFilterProducts().getTheta(), 
                detector.getEdgeFilterProducts().getGradientXY(), 
                img, settings, outputSkyCentroid, removedSets);
        
            seedPoints.add(points);
            
            Set<PairInt> embedded = new HashSet<PairInt>();
            Set<PairInt> borderPoints = new HashSet<PairInt>();
            SkylineExtractor.getEmbeddedAndBorderPoints(points,
                detector.getEdgeFilterProducts().getTheta().getWidth(), 
                detector.getEdgeFilterProducts().getTheta().getHeight(), 
                embedded, borderPoints);
            
            expectedBorderPoints.add(borderPoints);
            
            excludePoints.add(removedSets.getHighContrastRemoved());
            
        }
        
        SkylineANDedClauses skylineANDedClauses = new SkylineANDedClauses();
        ANDedClauses[] clauses;
        ANDedClauses[] coeffLowerLimits;
        ANDedClauses[] coeffUpperLimits;
         
        if (useAllImages) {
            clauses = skylineANDedClauses.getAllClauses();
            coeffLowerLimits = skylineANDedClauses.getAllClausesLowerLimits();
            coeffUpperLimits = skylineANDedClauses.getAllClausesUpperLimits();
        } else {
            if (useBlueSkyImages) {
                clauses = skylineANDedClauses.getGeneralAndBlueClauses();
                coeffLowerLimits = skylineANDedClauses.getGeneralAndBlueClausesLowerLimits();
                coeffUpperLimits = skylineANDedClauses.getGeneralAndBlueClausesUpperLimits();                
            } else {
                clauses = skylineANDedClauses.getGeneralAndRedClauses();
                coeffLowerLimits = skylineANDedClauses.getGeneralAndRedClausesLowerLimits();
                coeffUpperLimits = skylineANDedClauses.getGeneralAndRedClausesUpperLimits();
            }
        }
        
        List<SetComparisonResults> resultsBeforeList = new ArrayList<SetComparisonResults>();
        
        SkylineExtractor skylineExtractor = new SkylineExtractor();
        
        // ---- get the comparison of points before refinement ----
        for (int i = 0; i < seedPoints.size(); i++) {
            Set<PairInt> points0 = new HashSet<PairInt>(seedPoints.get(i));
            Set<PairInt> exclude = excludePoints.get(i);
            Set<PairInt> expectedBorder = expectedBorderPoints.get(i);
            ImageExt img = images.get(i);
            GreyscaleImage thetaImg = thetaImages.get(i);
         
            skylineExtractor.findClouds(points0, exclude, img, thetaImg,
                clauses);
            
            Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
            Set<PairInt> outputBorderPoints = new HashSet<PairInt>();
            SkylineExtractor.getEmbeddedAndBorderPoints(points0,
                thetaImg.getWidth(), thetaImg.getHeight(), 
                outputEmbeddedGapPoints, outputBorderPoints);
            
            SetCompare setCompare = new SetCompare();
            SetComparisonResults results = setCompare.compare(
                expectedSky.get(i), points0, expectedBorder, outputBorderPoints);
                            
            resultsBeforeList.add(results);
                        
            String fileNameRoot = fileNameRoots.get(i);
            
            System.out.println(fileNameRoot + " before optimization: " 
                + results.toString());
            
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
System.out.println("SIMPLEX");       
        // ---- get the comparison of points after refinement ----
        SkylineDownhillSimplex nelderMead = new SkylineDownhillSimplex(images, 
            thetaImages, seedPoints, excludePoints, expectedSky, 
            expectedBorderPoints,
            clauses, coeffLowerLimits, coeffUpperLimits
        );
        
        SkylineFits fit = nelderMead.fit();
        
        printFit(fit);
        
        SetComparisonResults resultsBefore = new SetComparisonResults(
            resultsBeforeList);
        
        SetComparisonResults resultsAfter = fit.results;
        
        // using the best fit, re-do the processing to see the results
        
        for (int i = 0; i < images.size(); i++) {
            
            skylineExtractor.findClouds(seedPoints.get(i), 
                excludePoints.get(i), images.get(i), thetaImages.get(i),
                fit.clauses);
            
            ImageExt img = images.get(i);
            GreyscaleImage thetaImg = thetaImages.get(i);
            
            String fileNameRoot = fileNameRoots.get(i);
            
            // to get the found border points, need to re-run the 
            // perimeter finder
            Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
            Set<PairInt> outputBorderPoints = new HashSet<PairInt>();
            SkylineExtractor.getEmbeddedAndBorderPoints(seedPoints.get(i),
                thetaImg.getWidth(), thetaImg.getHeight(), 
                outputEmbeddedGapPoints, outputBorderPoints);
            
            try {
                String dirPath = ResourceFinder.findDirectory("bin");
                ImageExt clr = (ImageExt) img.copyImage();
                ImageIOHelper.addToImage(seedPoints.get(i), 
                    thetaImg.getXRelativeOffset(),
                    thetaImg.getYRelativeOffset(), clr);
                ImageIOHelper.addToImage(outputBorderPoints, 
                    thetaImg.getXRelativeOffset(), 
                    thetaImg.getYRelativeOffset(), clr, 255, 0, 0);
                
                /*ImageIOHelper.addToImage(excludePoints.get(i), 
                    thetaImg.getXRelativeOffset(), 
                    thetaImg.getYRelativeOffset(), clr, 0, 255, 0);
                */
                ImageIOHelper.writeOutputImage(
                    dirPath + "/sky_after_optimization_" + fileNameRoot + ".png", clr);
            } catch (IOException e) {
                System.err.println("ERROR: " + e.getMessage());
            }
        }
        
        System.out.println("before : " + resultsBefore.toString());
        System.out.println("after : " + resultsAfter.toString());
        
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

    private void printFit(SkylineFits fit) {
        
        System.out.println("fit coefficients =");
        for (int i = 0; i < fit.clauses.length; i++) {
            ANDedClauses clause = fit.clauses[i];
            float[] coeff = clause.coefficients;
            System.out.println("clause " + i + ") " + Arrays.toString(coeff));
        }
        System.out.println("fit custom coefficients =");
        for (int i = 0; i < fit.clauses.length; i++) {
            ANDedClauses clause = fit.clauses[i];
            Map<Integer, Float> map = clause.customCoefficientVariables;
            if (!map.isEmpty()) {
                Iterator<Entry<Integer, Float>> iter = map.entrySet().iterator();
                while (iter.hasNext()) {
                    Entry<Integer, Float> entry = iter.next();
                    System.out.println("   clauseIndex=" + i + " coeffIndex=" 
                        + entry.getKey().toString() 
                        + " coeff=" + entry.getValue().toString());
                }
            }
        }
        System.out.println("fit results=" + fit.results.toString());
    }

}
