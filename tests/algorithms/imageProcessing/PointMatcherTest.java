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
        double rotation = 25 * Math.PI/180.;
        double scale = 4;
        int centroidX1 = 100;
        int centroidY1 = 100;
        
        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians((float)rotation);
        params.setScale((float)scale);
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
            transXTol, transYTol, rotation, scale, 
            centroidX1, centroidY1);
        
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
        fits[0] = new TransformationPointFit(params, 1, 10, 1);
        fits[1] = new TransformationPointFit(params, 3, 10, 1);
        fits[2] = new TransformationPointFit(params, 2, 10, 1);
        fits[3] = new TransformationPointFit(params, 100, 10, 1);
        fits[4] = new TransformationPointFit(params, 100, 0, 0); //avg diff=0
                
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
    public void testCalculateTransformation() throws Exception {
                
        PairIntArray set1 = new PairIntArray();
        set1.add(10, 10);
        set1.add(20, 20);
        set1.add(115, 120);
        
        PairIntArray set2 = new PairIntArray();
        set2.add(47, 240);
        set2.add(100, 259);
        set2.add(613, 461);
        
        PointMatcher matcher = new PointMatcher();
        
        TransformationPointFit fit = 
            matcher.calculateTransformation(set1, set2, 
            200, 200);
        
        assertNotNull(fit);
        assertNotNull(fit.getParameters());
                
        assertTrue(fit.getNumberOfMatchedPoints() == 3);
        
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
    public void testCalculateTransformation2() throws Exception {
        
        PairIntArray[] set1And2 = readCornersSets1And2();
        
        PointMatcher matcher = new PointMatcher();
        
        TransformationPointFit fit = 
            matcher.calculateTransformation(set1And2[0], set1And2[1], 
            517, 374);
        
        assertNotNull(fit);
        
        assertNotNull(fit.getParameters());
        
        assertTrue(Math.abs(fit.getParameters().getRotationInDegrees() - 0) < 5);
        
        assertTrue(Math.abs(fit.getParameters().getScale() - 1) < 1.0);
        
        assertTrue(Math.abs(fit.getParameters().getTranslationX() - -293.1) < 12);
        
        assertTrue(Math.abs(fit.getParameters().getTranslationY() - -14.3) < 6);
    }
    
    private PairIntArray[] readCornersSets1And2() throws IOException {
                
        /*
        corners were found in this project for the images
        below.  the images are from the Brown & Lowe 2003 paper, referenced
        in the file testresources/brown_lowe_2003.txt
        */
        String filePath1 = null;
        
        BufferedReader br = null;
        FileReader reader = null;
        
        PairIntArray xy1 = new PairIntArray();
        PairIntArray xy2 = new PairIntArray();
        
        String[] fileNames = new String[]{"brown_lowe_2003_image1.tsv", 
            "brown_lowe_2003_image2.tsv"};
        
        for (String fileName : fileNames) {
            
            filePath1 = ResourceFinder.findFileInTestResources(fileName);
            PairIntArray xy = fileName.equals("brown_lowe_2003_image1.tsv") 
                ? xy1 : xy2;
            
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
        
        return new PairIntArray[]{xy1, xy2};
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
