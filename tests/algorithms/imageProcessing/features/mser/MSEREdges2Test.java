package algorithms.imageProcessing.features.mser;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MSEREdges2Test extends TestCase {
    
    public MSEREdges2Test() {
    }

    public void test0() throws Exception {

        int maxDimension = 256;//512;

        String fileName1 = "";
        
        //for (int i = 3; i < 4; ++i) {
        for (int i = 0; i < 11; ++i) {

            switch(i) {
                case 0: {
                    fileName1 = "101085";
                    break;
                }
                case 1: {
                    fileName1 = "101087";
                    break;
                }
                case 2: {
                    fileName1 = "126007";
                    break;
                }
                case 3: {
                    fileName1 = "167062";
                    break;
                }
                case 4: {
                    fileName1 = "216081";
                    break;
                }
                case 5: {
                    fileName1 = "227092";
                    break;
                }
                case 6: {
                    fileName1 = "229036";
                    break;
                }
                case 7: {
                    fileName1 = "241004";
                    break;
                }
                case 8: {
                    fileName1 = "37073";
                    break;
                }
                case 9: {
                    fileName1 = "42049";
                    break;
                }
                default: {
                    fileName1 = "62096";
                    break;
                }
            }

            String filePath1 = ResourceFinder.findTestResourcesDirectory();
            filePath1 = filePath1 + "/berkeleySegSubset/" + fileName1 + "/" + fileName1 + ".jpg";
            
            ImageExt img = ImageIOHelper.readImageExt(filePath1);

            ImageProcessor imageProcessor = new ImageProcessor();

            int w1 = img.getWidth();
            int h1 = img.getHeight();

            int binFactor1 = (int) Math.ceil(Math.max(
                (float) w1 / maxDimension,
                (float) h1 / maxDimension));

            img = imageProcessor.binImage(img, binFactor1);
       //     img = (ImageExt) img.copySubImage(
       //         img.getWidth()/3, 2*img.getWidth()/3, 
       //         img.getHeight()/3, 2*img.getHeight()/3);
            //MiscDebug.writeImage(img, "_"  + fileName1Root);
                
            MSEREdges mserE = new MSEREdges(img);
            mserE.setToDebug();
            //mserE.extractEdges();
            mserE.mergeAndExtractEdges();
            
            /*GreyscaleImage[] lma = imageProcessor.createLCHForLUV(img);
            for (int k = 0; k < lma.length; ++k) {
                MiscDebug.writeImage(lma[k], "_" + 
                    fileName1Root + "_lma_" + k + "_");
            }
            
            GreyscaleImage gradients = imageProcessor
                .createGradientWithColorAndGreyscale(img);
        
            MiscDebug.writeImage(gradients,  
                "_" + fileName1Root + "_SOBEL_");
            */
        }
    }
    
}
