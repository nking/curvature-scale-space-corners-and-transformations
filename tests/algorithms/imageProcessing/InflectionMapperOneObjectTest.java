package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class InflectionMapperOneObjectTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());

    public InflectionMapperOneObjectTest() {
    }
    
    public void estImproveLineDrawingMode() throws Exception {
        
        String fileName2 = "closed_curve_translate_scale_rotate60.png";
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        GreyscaleImage img2 = ImageIOHelper.readImageExt(filePath2).copyToGreyscale();
        
        /*
        ILineThinner lineThinner = new ZhangSuenLineThinner();
        lineThinner.useLineDrawingMode();
        lineThinner.applyFilter(img2);
        */
        
        CannyEdgeFilter filter = new CannyEdgeFilter();
        filter.useLineDrawingMode();
        CannyEdgeFilterSettings settings = new CannyEdgeFilterSettings();
        settings.setUseLineDrawingMode();        
        filter.setSetters(settings);
        filter.applyFilter(img2);
        
        img2.multiply(200);
        
        String dirPath = ResourceFinder.findDirectory("bin");
        
        ImageIOHelper.writeOutputImage(
            dirPath + "/line_drawing_thinned.png", img2);
    }
    
    public void testMap() throws Exception {
        
        String[] rotDegreesList = new String[]{"20", "45", "60", "110", "160",
            "135", "180", "210", "225", "255", "280", "315", "335"
        };
        
        for (boolean swapDueToScale : new boolean[]{true, false}) {
            //swapDueToScale = false;
            for (String rotDegrees : rotDegreesList) {

                /*
                if (!rotDegrees.equals("45")) {
                    continue;
                }
                */

                String fileName1 = "closed_curve.png";
                String filePath1 = ResourceFinder.findFileInTestResources(fileName1);

                //String fileName2 = "closed_curve_translate_scale.png";
                String fileName2 = "closed_curve_translate_scale_rotate" + rotDegrees 
                    + ".png";

                String filePath2 = ResourceFinder.findFileInTestResources(fileName2);

                if (swapDueToScale) {
                    String swap = filePath1;
                    filePath1 = filePath2;
                    filePath2 = swap;
                }
                
                ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
                ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
                
                CurvatureScaleSpaceInflectionMapper mapper = new 
                    CurvatureScaleSpaceInflectionMapper(img1, img2);
                
                mapper.useLineDrawingLineMode();

                mapper.useDebugMode();

  //              mapper.setToRefineTransformations();

                TransformationParameters transformationParams = 
                    mapper.createEuclideanTransformation();

                assertNotNull(transformationParams);

                double rotDeg = transformationParams.getRotationInDegrees();

                double scale = transformationParams.getScale();

                int nEdges2 = mapper.getEdges2InOriginalReferenceFrame().size();
                PairIntArray[] edges2 = 
                    mapper.getEdges2InOriginalReferenceFrame().toArray(
                    new PairIntArray[nEdges2]);
                
                int nEdges1 = mapper.getEdges1InOriginalReferenceFrame().size();
                PairIntArray[] edges1 = 
                    mapper.getEdges1InOriginalReferenceFrame().toArray(
                    new PairIntArray[nEdges1]);

                Transformer transformer = new Transformer();
                PairIntArray[] transformedEdges = 
                    transformer.applyTransformation(transformationParams, 
                        edges1);

                img2 = ImageIOHelper.readImageExt(filePath2);

                MiscDebug.writeImage(transformedEdges, img2, "transformed_edges_" + rotDegrees);
/*
MiscDebug.writeImage(edges1, img1, "check_1_" + rotDegrees + "_" + MiscDebug.getCurrentTimeFormatted());
MiscDebug.writeImage(edges2, img2, "check_2_" + rotDegrees + "_" + MiscDebug.getCurrentTimeFormatted());
*/
                double expectedRotDeg = Float.valueOf(rotDegrees).floatValue();

                if (!swapDueToScale) {
                    expectedRotDeg = 360 - expectedRotDeg;
                }

                double foundRotDeg = rotDeg;

                log.info("PARAMS: " + transformationParams.toString() 
                    + "\nEXPECTED=" + rotDegrees + " (" + expectedRotDeg + ")"
                    + " found=" + foundRotDeg);

                float diffRot = AngleUtil.getAngleDifference(
                    (float)expectedRotDeg, (float)foundRotDeg);
                
                assertTrue(Math.abs(diffRot) < 10.f);

                if (rotDegrees.equals("135")) {
                    assertTrue(Math.abs(scale - 1.0) < 0.15);
                } else {
                    if (swapDueToScale) {
                        assertTrue(Math.abs(scale - (1./1.3)) < 0.15);
                    } else {
                        assertTrue(Math.abs(scale - 1.3) < 0.15);
                    }
                }
            }
        }
    }
    
    public static void main(String[] args) {
        
        try {
            
            InflectionMapperOneObjectTest test = 
                new InflectionMapperOneObjectTest();
                        
            test.testMap();
            //test.testImproveLineDrawingMode();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
