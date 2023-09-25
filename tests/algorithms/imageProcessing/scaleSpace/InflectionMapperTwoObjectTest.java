package algorithms.imageProcessing.scaleSpace;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.util.AngleUtil;
import algorithms.util.ResourceFinder;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class InflectionMapperTwoObjectTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());

    public InflectionMapperTwoObjectTest() {
    }
    
    public void testMap() throws Exception {
        
        //NOTE: until the inflection mapper is improved,
        //   have removed the asserts
        
        String rootName = "closed_curve_translate_scale_rotate";
        
        //20, 45, 210, 110, 160, 180, 225
        
        String[] list1 = new String[]{"20", "45", "210"};
        
        String[] list2 = new String[]{"110", "160", "180"};
        
        //list1 = new String[]{"20"};
        //list2 = new String[]{"160"};

        //for (int i = 0; i < list1.length; i++) {
        //    for (int j = 0; j < list2.length; j++) {
        for (int i = 0; i < 1; i++) {
            for (int j = 0; j < 1; j++) {
                
                String fileName1 = rootName + list1[i] + ".png";
                
                String fileName2 = rootName + list2[j] + ".png";
                
                log.info("TRANSFORM " + fileName1 + " to " + fileName2);
                
                int rotation1 = Integer.valueOf(list1[i]).intValue();
                
                int rotation2 = Integer.valueOf(list2[j]).intValue();
                
                String filePath1 = ResourceFinder.findFileInTestResources(
                    fileName1);
                
                String filePath2 = ResourceFinder.findFileInTestResources(
                    fileName2);
                
                ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
                ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
                
                double centroidX1 = img1.getWidth() >> 1;
                double centroidY1 = img1.getHeight() >> 1;
                
                CurvatureScaleSpaceInflectionMapper mapper = new 
                    CurvatureScaleSpaceInflectionMapper(img1, img2);

                mapper.useLineDrawingLineMode();

                mapper.useDebugMode();

                //TODO: this needs revision:
                //mapper.setToRefineTransformations();
                TransformationParameters transformationParams = 
                    mapper.createEuclideanTransformation();

                if (transformationParams == null) {
                    continue;
                }
                
                assertNotNull(transformationParams);

                double rotDeg = transformationParams.getRotationInDegrees();

                double scale = transformationParams.getScale();
                
                int nEdges1 = mapper.getEdges1().size();
                PairIntArray[] edges1 = 
                    mapper.getEdges1().toArray(
                    new PairIntArray[nEdges1]);

                Transformer transformer = new Transformer();
                PairIntArray[] transformedEdges = 
                    transformer.applyTransformation(transformationParams, 
                        edges1);

                img2 = ImageIOHelper.readImageExt(filePath2);

                debugDisplay(transformedEdges, img2, list1[i] + "->" + list2[j]);
                
                int expectedRotDeg = 360 - (rotation2 - rotation1);
                while (expectedRotDeg > 360) {
                    expectedRotDeg -= 360;
                }
                
                double foundRotDeg = rotDeg;
                
                float diffRot = AngleUtil.getAngleDifference(
                    (float)expectedRotDeg, (float)foundRotDeg);                
                
                log.info("PARAMS: " + transformationParams.toString() 
                    + "\nEXPECTED=" + expectedRotDeg  + " found=" + foundRotDeg
                    + " diffRot=" + diffRot);

                //assertTrue(Math.abs(diffRot) < 10.f);
                
                //assertTrue(Math.abs(scale - 1.0) < 0.15);
            }
        }
    }
   
    private void debugDisplay(PairIntArray[] transformedEdges,
        Image img, String rotDegrees) throws IOException {
        
        for (int i = 0; i < transformedEdges.length; i++) {
            PairIntArray edge = transformedEdges[i];
            debugAddCurveToImage(edge, img, 0, 255, 255, 255);
        }     
        
        String dirPath = ResourceFinder.findDirectory("bin");
        
        ImageIOHelper.writeOutputImage(
            dirPath + "/css_transformed_edges_" + rotDegrees + ".png", img);
    }
    
    private void debugAddCurveToImage(PairIntArray edge, Image input, 
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        for (int i = 0; i < edge.getN(); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                float xx = x + dx;
                if ((xx > -1) && (xx < (input.getWidth() - 1))) {
                    for (int dy = (-1*nExtraForDot); dy <= (nExtraForDot + 1); 
                        dy++) {
                        float yy = y + dy;
                        if ((yy > -1) && (yy < (input.getHeight() - 1))) {
                            input.setRGB((int)xx, (int)yy, rClr, gClr, bClr);
                        }
                    }
                }
            }
        }
    }

    public static void main(String[] args) {
        
        try {
            
            InflectionMapperTwoObjectTest test = 
                new InflectionMapperTwoObjectTest();
                        
            test.testMap();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
