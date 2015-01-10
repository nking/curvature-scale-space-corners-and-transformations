package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import algorithms.util.PairIntArray;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Logger;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class PointMatcherTest {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public PointMatcherTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testNumberOfMatches() {
        
        PointMatcher matcher = new PointMatcher();
        
        PairIntArray set1 = new PairIntArray();
        set1.add(10, 10);
        set1.add(20, 20);
        
        PairIntArray set2 = new PairIntArray();
        set2.add(47, 240);
        set2.add(100, 259);
        
        int transX = 125;
        int transY = 14;
        double transXTol = 10.3;
        double transYTol = 5.9;
        float rotation = (float)(25.f * Math.PI/180.f);
        float scale = 4.f;
        int centroidX1 = 100;
        int centroidY1 = 100;
        
        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians(rotation);
        params.setScale(scale);
        params.setTranslationX(transX);
        params.setTranslationY(transY);
        
        TransformationPointFit fit = matcher.transform(set1, set2, params, 
            transXTol, transYTol, centroidX1, centroidY1);
        
        assertNotNull(fit);
        assertNotNull(fit.getParameters());
        
        assertTrue(fit.getNumberOfMatchedPoints() == 2);
        
        assertTrue(fit.getMeanDistFromModel() == 0);
        
        assertTrue(fit.getStDevFromMean() == 0);
        
        assertTrue(fit.getParameters().getScale() == 4);
        
        assertTrue((Math.abs(fit.getParameters().getRotationInDegrees()) 
            - 25) < 1);
        
        assertTrue(
            Math.abs(
                fit.getParameters().getRotationInRadians() - (25*Math.PI/180)) 
                < (5*Math.PI/180));
        
        assertTrue((Math.abs(fit.getParameters().getTranslationX()) 
            - 125) < 1);
        
        assertTrue((Math.abs(fit.getParameters().getTranslationY()) 
            - 14) < 1);
    }
    
    @Test
    public void testCalculateTranslation() {
        
        PointMatcher matcher = new PointMatcher();
        
        PairIntArray set1 = new PairIntArray();
        set1.add(10, 10);
        set1.add(20, 20);
        
        PairIntArray set2 = new PairIntArray();
        set2.add(47, 240);
        set2.add(100, 259);
        
        int transX = 125;
        int transY = 14;
        double transXTol = 10.3;
        double transYTol = 5.9;
        double rotation = 25*Math.PI/180.;
        double scale = 4;
        int centroidX1 = 100;
        int centroidY1 = 100;
        
        TransformationPointFit fit = 
            matcher.calculateTranslation(set1, set2,
            rotation, scale, 
            centroidX1, centroidY1, true, 1.0f);
        
        assertNotNull(fit);
        assertNotNull(fit.getParameters());
        
        assertTrue(fit.getNumberOfMatchedPoints() == 2);
        
        assertTrue(fit.getMeanDistFromModel() == 0);
        
        assertTrue(fit.getStDevFromMean() == 0);
    }
   
    @Test
    public void testSortByDescendingMatches() throws Exception {
        
        TransformationParameters params = new TransformationParameters();
        
        TransformationPointFit[] fits = new TransformationPointFit[5];
        fits[0] = new TransformationPointFit(params, 1, 10, 1, 100);
        fits[1] = new TransformationPointFit(params, 3, 10, 1, 100);
        fits[2] = new TransformationPointFit(params, 2, 10, 1, 100);
        fits[3] = new TransformationPointFit(params, 100, 10, 1, 100);
        fits[4] = new TransformationPointFit(params, 100, 0, 0, 100); //avg diff=0
                
        PointMatcher matcher = new PointMatcher();
        matcher.sortByDescendingMatches(fits, 0, fits.length - 1);
        
        assertTrue(fits[0].getNumberOfMatchedPoints() == 100);
        assertTrue(fits[0].getMeanDistFromModel() == 0);
        assertTrue(fits[1].getNumberOfMatchedPoints() == 100);
        assertTrue(fits[1].getMeanDistFromModel() == 10);
        assertTrue(fits[2].getNumberOfMatchedPoints() == 3);
        assertTrue(fits[3].getNumberOfMatchedPoints() == 2);
        assertTrue(fits[4].getNumberOfMatchedPoints() == 1);
    }
    
    @Test
    public void testCreateIntervals() throws Exception {
        
        PointMatcher matcher = new PointMatcher();
        double[] r = matcher.createIntervals(0, 360, 10);
        
        assertNotNull(r);
        
        assertTrue(r.length == 36);
        
        for (int i = 0; i < r.length; i++) {
            double expected = i * 10;
            assertTrue(expected == r[i]);
        }
   
        double[] s = matcher.createIntervals(1, 11, 1);
        
        assertNotNull(s);
        
        assertTrue(s.length == 10);
        
        for (int i = 0; i < s.length; i++) {
            double expected = i + 1;
            assertTrue(expected == s[i]);
        }
    }
}
