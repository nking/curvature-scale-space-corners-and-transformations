package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ATrousWaveletTransformTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public void test() throws Exception {
        
        Set<String> testFiles = new HashSet<String>();
        testFiles.add("blox.gif");
        testFiles.add("house.gif");
        /*testFiles.add("lab.gif");
        testFiles.add("africa2.png");
        testFiles.add("susan-in.gif");
        testFiles.add("valve_gaussian.png");
        testFiles.add("lena.jpg");
        testFiles.add("android_statues_01.jpg");
        testFiles.add("android_statues_04.jpg");*/
        
        for (String fileName : testFiles) {
                        
            String filePath = ResourceFinder.findFileInTestResources(fileName);
            
            int idx = fileName.lastIndexOf(".");
            String fileNameRoot = fileName.substring(0, idx);
            
            log.info("fileName=" + fileName);
                        
            ImageExt img = ImageIOHelper.readImageExt(filePath);
            
            GreyscaleImage gsImg = img.copyToGreyscale();
            GreyscaleImage gsImgCp = img.copyToGreyscale();

            List<GreyscaleImage> outputTransformed = new ArrayList<GreyscaleImage>();
            List<GreyscaleImage> outputCoeff = new ArrayList<GreyscaleImage>();
            
            ATrousWaveletTransform at = new ATrousWaveletTransform();
            
            at.calculateWithB3SplineScalingFunction(gsImg, outputTransformed, outputCoeff);
            
            List<GreyscaleImage> outputTransformed2 = new ArrayList<GreyscaleImage>();
            List<GreyscaleImage> outputCoeff2 = new ArrayList<GreyscaleImage>();
            
            at.calculateWithTriangleScalingFunction(gsImg, 
                outputTransformed2, outputCoeff2);
                
            for (int i = 0; i < outputTransformed.size(); ++i) {
                MiscDebug.writeImage(outputTransformed.get(i), fileNameRoot 
                    + "_b3_tr_" + i + "_");
            }
            for (int i = 0; i < outputCoeff.size(); ++i) {
                MiscDebug.writeImage(outputCoeff.get(i), fileNameRoot + 
                    "_b3_coeff_" + i + "_");
            }
            
            for (int i = 0; i < outputTransformed2.size(); ++i) {
                MiscDebug.writeImage(outputTransformed2.get(i), fileNameRoot 
                    + "_b3_tr_" + i + "_");
            }
            for (int i = 0; i < outputCoeff2.size(); ++i) {
                MiscDebug.writeImage(outputCoeff2.get(i), fileNameRoot + 
                    "_b3_coeff_" + i + "_");
            }
            
            GreyscaleImage r = at.reconstruct(
                outputTransformed.get(outputTransformed.size() - 1), outputCoeff);
            
            MiscDebug.writeImage(r, fileNameRoot + "_reconstructed");
            
            r = r.subtract(gsImgCp);
            
            MiscDebug.writeImage(r, fileNameRoot + "_recon_resid");
            
            
            r = at.reconstruct(
                outputTransformed2.get(outputTransformed2.size() - 1), outputCoeff2);
            
            MiscDebug.writeImage(r, fileNameRoot + "_reconstructed_tri");
            
            r = r.subtract(gsImgCp);
            
            MiscDebug.writeImage(r, fileNameRoot + "_recon_resid_tri");
            
        }
        
    }
    
}
