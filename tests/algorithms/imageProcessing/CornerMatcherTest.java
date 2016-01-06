package algorithms.imageProcessing;

import algorithms.compGeometry.RotatedOffsets;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.util.ArrayList;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class CornerMatcherTest extends TestCase {
    
    public CornerMatcherTest() {        
    }
    
    public void testMatch() throws Exception {
        
        int dither = 1;
        int binFactor1 = 1;
        int binFactor2 = 1;
        
        String fileName1 = "books_illum3_v0_695x555.png";
        String fileName2 = "books_illum3_v6_695x555.png";
                    
        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);
        idx = fileName2.lastIndexOf(".");
        String fileName2Root = fileName2.substring(0, idx);
        
        FeatureMatcherSettings settings = new FeatureMatcherSettings();
        settings.setDebug(true);
        settings.setStartWithBinnedImages(true);
        settings.setDebugTag(fileName1Root);
        
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
          
        final int blockHalfWidth = 5;
        final boolean useNormalizedIntensities = true;
        RotatedOffsets rotatedOffsets = RotatedOffsets.getInstance();
        
        IntensityFeatures features1 = new IntensityFeatures(blockHalfWidth, 
            useNormalizedIntensities, rotatedOffsets);
        features1.calculateGradientWithGreyscale(img1.copyToGreyscale());
        
        IntensityFeatures features2 = new IntensityFeatures(blockHalfWidth, 
            useNormalizedIntensities, rotatedOffsets);
        features2.calculateGradientWithGreyscale(img2.copyToGreyscale());
        
        final List<CornerRegion> corners1 = new ArrayList<CornerRegion>();
        final List<CornerRegion> corners2 = new ArrayList<CornerRegion>();
        populatePoints(corners1, corners2);
        
        CornerMatcher<CornerRegion> matcher = new CornerMatcher<CornerRegion>(dither);
        
        boolean matched = matcher.matchCorners(features1, features2, 
            corners1, corners2, img1.copyToGreyscale(), img2.copyToGreyscale(),
            binFactor1, binFactor2);
        
        assertTrue(matched);
        
        List<FeatureComparisonStat> stats = matcher.getSolutionStats();
        assertTrue((int)(0.8f*corners1.size()) <= stats.size());
        
        List<PairInt> matched1 = matcher.getMatched1();
        
        List<PairInt> matched2 = matcher.getMatched2();
        
        assertTrue((int)(0.8f*corners1.size()) <= matched1.size());
        
        assertEquals(matched1.size(), matched2.size());
        
    }
    
    private void populatePoints(List<CornerRegion> corners1, 
        List<CornerRegion> corners2) {
        
        /*
        img1       img2
        310, 57    246, 57
        397, 37    334, 37
        538, 43    478, 43
        530, 423   412, 424
        427, 438   320, 439
        164, 203    71, 203
        */
        
        CornerRegion cr = new CornerRegion(0, 1, 0);
        cr.setFlagThatNeighborsHoldDummyValues();
        cr.set(0, Float.MIN_VALUE, 310, 57);
        cr.setIndexWithinCurve(0);
        corners1.add(cr);
        
        cr = new CornerRegion(0, 1, 0);
        cr.setFlagThatNeighborsHoldDummyValues();
        cr.set(0, Float.MIN_VALUE, 246, 57);
        cr.setIndexWithinCurve(0);
        corners2.add(cr);
        
        cr = new CornerRegion(0, 1, 0);
        cr.setFlagThatNeighborsHoldDummyValues();
        cr.set(0, Float.MIN_VALUE, 397, 37);
        cr.setIndexWithinCurve(1);
        corners1.add(cr);
        
        cr = new CornerRegion(0, 1, 0);
        cr.setFlagThatNeighborsHoldDummyValues();
        cr.set(0, Float.MIN_VALUE, 334, 37);
        cr.setIndexWithinCurve(1);
        corners2.add(cr);
        
        
        cr = new CornerRegion(0, 1, 0);
        cr.setFlagThatNeighborsHoldDummyValues();
        cr.set(0, Float.MIN_VALUE, 538, 43);
        cr.setIndexWithinCurve(2);
        corners1.add(cr);
        
        cr = new CornerRegion(0, 1, 0);
        cr.setFlagThatNeighborsHoldDummyValues();
        cr.set(0, Float.MIN_VALUE, 478, 43);
        cr.setIndexWithinCurve(2);
        corners2.add(cr);
        
        
        cr = new CornerRegion(0, 1, 0);
        cr.setFlagThatNeighborsHoldDummyValues();
        cr.set(0, Float.MIN_VALUE, 530, 423);
        cr.setIndexWithinCurve(3);
        corners1.add(cr);
        
        cr = new CornerRegion(0, 1, 0);
        cr.setFlagThatNeighborsHoldDummyValues();
        cr.set(0, Float.MIN_VALUE, 412, 424);
        cr.setIndexWithinCurve(3);
        corners2.add(cr);
        
        
        cr = new CornerRegion(0, 1, 0);
        cr.setFlagThatNeighborsHoldDummyValues();
        cr.set(0, Float.MIN_VALUE, 427, 438);
        cr.setIndexWithinCurve(4);
        corners1.add(cr);
        
        cr = new CornerRegion(0, 1, 0);
        cr.setFlagThatNeighborsHoldDummyValues();
        cr.set(0, Float.MIN_VALUE, 320, 439);
        cr.setIndexWithinCurve(4);
        corners2.add(cr);
        
        
        cr = new CornerRegion(0, 1, 0);
        cr.setFlagThatNeighborsHoldDummyValues();
        cr.set(0, Float.MIN_VALUE, 164, 203);
        cr.setIndexWithinCurve(5);
        corners1.add(cr);
        
        cr = new CornerRegion(0, 1, 0);
        cr.setFlagThatNeighborsHoldDummyValues();
        cr.set(0, Float.MIN_VALUE, 71, 203);
        cr.setIndexWithinCurve(5);
        corners2.add(cr);
    }
}
