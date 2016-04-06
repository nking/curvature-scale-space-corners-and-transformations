package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import java.util.HashSet;
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
            filter.setToRestoreJunctions();
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
    
}
