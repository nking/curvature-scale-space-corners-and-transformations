package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class CannyEdgeFilterAdaptiveTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public void test() throws Exception {
        
        Set<String> testFiles = new HashSet<String>();
        testFiles.add("blox.gif");
        //testFiles.add("two_circles_color.png");
        testFiles.add("house.gif");
        /*testFiles.add("lab.gif");
        testFiles.add("africa2.png");
        testFiles.add("susan-in.gif");
        testFiles.add("valve_gaussian.png");
        testFiles.add("lena.jpg");
        testFiles.add("android_statues_01.jpg");
        testFiles.add("android_statues_04.jpg");
        */
        testFiles.add("checkerboard_01.jpg");
        
        for (String fileName : testFiles) {
            
            String filePath = ResourceFinder.findFileInTestResources(fileName);
            
            int idx = fileName.lastIndexOf(".");
            String fileNameRoot = fileName.substring(0, idx);
            
            log.info("fileName=" + fileName);
          
            ImageExt img = ImageIOHelper.readImageExt(filePath);
            
            GreyscaleImage gsImg = img.copyToGreyscale();

            CannyEdgeFilterAdaptive filter = new CannyEdgeFilterAdaptive();
            //CannyEdgeFilterLite filter = new CannyEdgeFilterLite();
            //filter.setToNotUseNonMaximumSuppression();
            //filter.setToPerformHistogramEqualization();
            //filter.setOtsuScaleFactor(0.2f);
            //filter.override2LayerFactorBelowHighThreshold(10.f);
            filter.overrideDefaultNumberOfLevels(16);
            if (fileName.contains("africa") || fileName.contains("circle") 
                || fileName.contains("susan")) {
                filter.setToUseLineDrawingMode();
            }
            filter.setToDebug();
            filter.applyFilter(gsImg);
            
            for (int i = 0; i < gsImg.getNPixels(); ++i) {
                if (gsImg.getValue(i) > 0) {
                    gsImg.setValue(i, 255);
                }
            }
            
            MiscDebug.writeImage(gsImg, "_canny_adaptive_" + fileNameRoot);
        }
        
    }
    
    public void test1() throws Exception {
        
        //String fileName = "blox.gif";
        //String fileName = "lab.gif";
        //String fileName = "house.gif";
        //String fileName = "susan-in_plus.png";
        String fileName = "africa2.png";
        //String fileName = "lena.jpg";
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();
       
        CannyEdgeFilterAdaptive canny1 = new CannyEdgeFilterAdaptive();        
        canny1.setToUseLineDrawingMode();;
        //canny1.setToDebug();
        canny1.applyFilter(img);
        for (int i = 0; i < img.getNPixels(); ++i) {
            if (img.getValue(i) > 0) {
                img.setValue(i, 255);
            }
        }
        //MiscDebug.writeImage(img, "_gradient_canny_adaptive_");
        
        EdgeExtractorWithJunctions extractor = new EdgeExtractorWithJunctions(img);
        List<PairIntArray> edges = extractor.findEdges();
        int n = 0;
        for (PairIntArray edge : edges) {
            if (edge.getN() > 3) {
                n++;
            }
        }
        Map<Integer, Set<Integer>> junctionMap = extractor.getJunctionMap();
        assertEquals(1, n);
    }
    
}
