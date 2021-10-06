package algorithms.imageProcessing.features.mser;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import gnu.trove.set.TIntSet;
import java.util.List;
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

            MSEREdgesWrapper msew = new MSEREdgesWrapper();
                
            MSEREdges mserE = msew.extractAndMergeEdges(img);
            
            List<TIntSet> edgeList = mserE.getEdges();
            Image im = mserE.getGsImg().copyToColorGreyscale();
            int[] clr = new int[]{255, 0, 0};
            for (int ii = 0; ii < edgeList.size(); ++ii) {
                ImageIOHelper.addCurveToImage(edgeList.get(ii), im, 0, clr[0],
                    clr[1], clr[2]);
            }
            MiscDebug.writeImage(im, "_" + fileName1 + "_mser_edges_");
        }
    }
    
}
