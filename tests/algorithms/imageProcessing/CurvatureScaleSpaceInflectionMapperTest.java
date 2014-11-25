package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.util.List;
import java.util.logging.Logger;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class CurvatureScaleSpaceInflectionMapperTest {
    
    private Logger log = Logger.getLogger(this.getClass().getName());

    public CurvatureScaleSpaceInflectionMapperTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    @Test
    public void testMap() throws Exception {
        
        String[] rotDegreesList = new String[]{"20", "60", "135", "180", "225",
            "280", "335"
        };
        
        for (boolean swapDueToScale : new boolean[]{false, true}) {
        
            for (String rotDegrees : rotDegreesList) {

                /*
                if (!rotDegrees.equals("180")) {
                    continue;
                }
                */

                String fileName1 = "closed_curve.png";
                String filePath1 = ResourceFinder.findFileInTestResources(fileName1);

                //String fileName2 = "closed_curve_translate.png";
                //String fileName2 = "closed_curve_translate_scale.png";
                String fileName2 = "closed_curve_translate_scale_rotate" + rotDegrees 
                    + ".png";

                String filePath2 = ResourceFinder.findFileInTestResources(fileName2);

                if (swapDueToScale) {
                    String swap = filePath1;
                    filePath1 = filePath2;
                    filePath2 = swap;
                }
                
                GreyscaleImage img1 = ImageIOHelper.readImageAsGrayScaleG(filePath1);
                GreyscaleImage img2 = ImageIOHelper.readImageAsGrayScaleG(filePath2);
                
                double centroidX1 = img1.getWidth() >> 1;
                double centroidY1 = img1.getHeight() >> 1;

                CurvatureScaleSpaceInflectionMapper mapper = new 
                    CurvatureScaleSpaceInflectionMapper(img1, img2);

                mapper.useLineDrawingLineMode();

                mapper.useDebugMode();

                mapper.setToRefineTransformations();

                TransformationParameters transformationParams = 
                    mapper.createEuclideanTransformation();

                assertNotNull(transformationParams);

                float rotDeg = transformationParams.getRotationInDegrees();

                float scale = transformationParams.getScale();

                int nEdges1 = mapper.getEdges1InOriginalReferenceFrame().size();
                PairIntArray[] edges1 = 
                    mapper.getEdges1InOriginalReferenceFrame().toArray(
                    new PairIntArray[nEdges1]);

                Transformer transformer = new Transformer();
                PairIntArray[] transformedEdges = 
                    transformer.applyTransformation(transformationParams, 
                        edges1, centroidX1, centroidY1);

                img2 = ImageIOHelper.readImageAsGrayScaleG(filePath2);

                debugDisplay(transformedEdges, img2.copyImageToGreen(), rotDegrees);

                double expectedRotDeg = Float.valueOf(rotDegrees).floatValue();

                if (!swapDueToScale) {
                    if (rotDegrees.equals("20")) {
                        expectedRotDeg = 360 - expectedRotDeg;
                    } else if (rotDegrees.equals("60")) {
                        expectedRotDeg = 360 - expectedRotDeg;
                    } else if (rotDegrees.equals("135")) {
                        expectedRotDeg = 360 - expectedRotDeg;
                    } else if (rotDegrees.equals("180")) {
                        expectedRotDeg = 360 - expectedRotDeg;
                    } else if (rotDegrees.equals("225")) {
                        expectedRotDeg = 360 - expectedRotDeg;
                    } else if (rotDegrees.equals("280")) {
                        expectedRotDeg = 360 - expectedRotDeg;
                    } else if (rotDegrees.equals("335")) {
                        expectedRotDeg = 360 - expectedRotDeg;
                    }
                }

                double foundRotDeg = rotDeg;

                log.info("PARAMS: " + transformationParams.toString() 
                    + "\nEXPECTED=" + rotDegrees + " (" + expectedRotDeg + ")");

                assertTrue(Math.abs(expectedRotDeg - foundRotDeg) < 10.f);

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
    
    //@Test
    public void estMap2() throws Exception {
       
        String fileName1 = "closed_curves.png";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        GreyscaleImage img1 = ImageIOHelper.readImageAsGrayScaleG(filePath1);

        String fileName2 = "closed_curves_translate_scale_rotate60.png";

        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        GreyscaleImage img2 = ImageIOHelper.readImageAsGrayScaleG(filePath2);

        double centroidX1 = img1.getWidth() >> 1;
        double centroidY1 = img1.getHeight() >> 1;
            
        CurvatureScaleSpaceInflectionMapper mapper = new 
            CurvatureScaleSpaceInflectionMapper(img1, img2);

        mapper.useLineDrawingLineMode();

        mapper.useDebugMode();

        //mapper.doNotRefineTransformations();


        /*
        scale 2.1
        rotate 60
        translate 25, 0
        */

        TransformationParameters transformationParams = 
            mapper.createEuclideanTransformation();

        assertNotNull(transformationParams);

        float rotDeg = -1*transformationParams.getRotationInDegrees();

        float scale = transformationParams.getScale();

        int nEdges1 = mapper.getEdges1InOriginalReferenceFrame().size();
        PairIntArray[] edges1 = 
            mapper.getEdges1InOriginalReferenceFrame().toArray(
            new PairIntArray[nEdges1]);

        Transformer transformer = new Transformer();
        PairIntArray[] transformedEdges = 
            transformer.applyTransformation(transformationParams, edges1,
            centroidX1, centroidY1);

        img2 = ImageIOHelper.readImageAsGrayScaleG(filePath2);
            
        debugDisplay(transformedEdges, img2.copyImageToGreen(), "60");

        log.info("PARAMS: " + transformationParams.toString() + "\nEXPECTED=60");

        assertTrue(Math.abs(rotDeg - 60) < 10.f);

        assertTrue(Math.abs(scale - 2.1) < 0.15);
    }
    
    //@Test
    public void estTransformerGreyscaleImage() throws Exception {
        
        String[] rotDegrees = new String[]{"225"};
        
        String fileName1 = "closed_curve.png";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        GreyscaleImage img1 = ImageIOHelper.readImageAsGrayScaleG(filePath1);

        String fileName2 = "closed_curve_translate_scale_rotate225.png";

        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        GreyscaleImage img2 = ImageIOHelper.readImageAsGrayScaleG(filePath2);

        CurvatureScaleSpaceInflectionMapper mapper = new 
            CurvatureScaleSpaceInflectionMapper(img1, img2);

        mapper.useLineDrawingLineMode();
        
        mapper.useDebugMode();

        TransformationParameters transformationParams = 
            mapper.createEuclideanTransformation();

        assertNotNull(transformationParams);

        float rotDeg = -1*transformationParams.getRotationInDegrees();
        
        assertTrue(Math.abs(rotDeg - 225.f) < 10.f);

        float scale = transformationParams.getScale();

        assertTrue(Math.abs(scale - 1.3) < 0.15);

        List<PairIntArray> edges1 = 
            mapper.getEdges1InOriginalReferenceFrame();
        int nPointsInEdges1 = 0;
        for (int i = 0; i < edges1.size(); i++) {
            nPointsInEdges1 += edges1.get(i).getN();
        }

        Transformer transformer = new Transformer();
        
        img2 = ImageIOHelper.readImageAsGrayScaleG(filePath2);
        GreyscaleImage img2Masked = img2.copyImage();
        for (int i = 180; i < img2Masked.getWidth(); i++) {
            for (int j = 0; j < img2Masked.getHeight(); j++) {
                img2Masked.setValue(i, j, 0);
            }
        }
        
        img1 = ImageIOHelper.readImageAsGrayScaleG(filePath1);
        
        GreyscaleImage img1Transformed = transformer.applyTransformation(
            img1, transformationParams, img2.getWidth(), img2.getHeight());
        
        /*String dirPath = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(dirPath + "/transformed_to_225.png", 
            img1Transformed);
        */
        
        int residuals = 0;
        for (int i = 0; i < img2Masked.getWidth(); i++) {
            for (int j = 0; j < img2Masked.getHeight(); j++) {
                int pix2 = img2Masked.getValue(i, j);
                int pixT = img1Transformed.getValue(i, j);
                if (((pix2 != 0) && (pixT == 0)) || ((pix2 == 0) && (pixT != 0))) {
                    residuals++;
                }
            }
        }
        
        residuals /= 2;
        
        log.info("number of pixels not overlapping in transformed edge = " +
            residuals + " (n known in edges=" + nPointsInEdges1);
        
        int nBin = 0;
        
        int w = img2Masked.getWidth();
        int h = img2Masked.getHeight();
        while ((residuals > 0) && (img2Masked.getWidth() > 2)) {
            w >>= 1;
            h >>= 1;
            
            GreyscaleImage i2 = new GreyscaleImage(w, h);
            GreyscaleImage i1 = new GreyscaleImage(w, h);
            for (int i = 0; i < w; i++) {
                int x = i*2;
                for (int j = 0; j < h; j++) {
                    int y = j*2;
                    int sum = img2Masked.getValue(x, y) 
                        + img2Masked.getValue(x, y + 1)
                        + img2Masked.getValue(x + 1, y)
                        + img2Masked.getValue(x + 1, y + 1);
                    i2.setValue(i, j, sum);
                }
            }
            img2Masked = i2;
            
            for (int i = 0; i < w; i++) {
                int x = i*2;
                for (int j = 0; j < h; j++) {
                    int y = j*2;
                    int sum = img1Transformed.getValue(x, y) 
                        + img1Transformed.getValue(x, y + 1)
                        + img1Transformed.getValue(x + 1, y)
                        + img1Transformed.getValue(x + 1, y + 1);
                    i1.setValue(i, j, sum);
                }
            }
            img1Transformed = i1;
            
            residuals = 0;
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                    int pix2 = img2Masked.getValue(i, j);
                    int pixT = img1Transformed.getValue(i, j);
                    if (((pix2 != 0) && (pixT == 0)) || 
                        ((pix2 == 0) && (pixT != 0))) {
                        
                        residuals++;
                    }
                }
            }
            
            residuals /= 2;
            
            nBin++;
        
            log.info("number of pixels not overlapping in transformed edge = " +
                residuals + " for bin iteration " + nBin);
        }
       
        // scale is close, rotation is close, but translation is off by a couple
        // of pixels.
        
        // generous assert for now
        assertTrue(nBin < 7);
    }
    
    //@Test
    public void estTransformerImage() throws Exception {
        
        String[] rotDegrees = new String[]{"225"};
        
        String fileName1 = "closed_curve.png";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        GreyscaleImage img1 = ImageIOHelper.readImageAsGrayScaleG(filePath1);

        String fileName2 = "closed_curve_translate_scale_rotate225.png";

        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        GreyscaleImage img2 = ImageIOHelper.readImageAsGrayScaleG(filePath2);

        CurvatureScaleSpaceInflectionMapper mapper = new 
            CurvatureScaleSpaceInflectionMapper(img1, img2);

        mapper.useLineDrawingLineMode();

        TransformationParameters transformationParams = 
            mapper.createEuclideanTransformation();

        assertNotNull(transformationParams);

        float rotDeg = -1*transformationParams.getRotationInDegrees();
        
        assertTrue(Math.abs(rotDeg - 225.f) < 10.f);

        float scale = transformationParams.getScale();

        assertTrue(Math.abs(scale - 1.3) < 0.15);

        List<PairIntArray> edges1 = 
            mapper.getEdges1InOriginalReferenceFrame();
        int nPointsInEdges1 = 0;
        for (int i = 0; i < edges1.size(); i++) {
            nPointsInEdges1 += edges1.get(i).getN();
        }

        Transformer transformer = new Transformer();
        
        img2 = ImageIOHelper.readImageAsGrayScaleG(filePath2);
            
        Image img2Masked = img2.copyImageToGreen();
        for (int i = 180; i < img2Masked.getWidth(); i++) {
            for (int j = 0; j < img2Masked.getHeight(); j++) {
                img2Masked.setRGB(i, j, 0, 0, 0);
            }
        }
        
        img1 = ImageIOHelper.readImageAsGrayScaleG(filePath1);
        
        Image img1Transformed = transformer.applyTransformation(
            img1.copyImageToGreen(), transformationParams, img2Masked.getWidth(), 
            img2Masked.getHeight());
        
        /*String dirPath = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(dirPath + "/transformed_to_225.png", 
            img1Transformed);*/
        
        int residuals = 0;
        for (int i = 0; i < img2Masked.getWidth(); i++) {
            for (int j = 0; j < img2Masked.getHeight(); j++) {
                int pix2 = img2Masked.getG(i, j);
                int pixT = img1Transformed.getG(i, j);
                if (((pix2 != 0) && (pixT == 0)) 
                    || ((pix2 == 0) && (pixT != 0))) {
                    
                    residuals++;
                }
            }
        }
        
        residuals /= 2;
        
        log.info("number of pixels not overlapping in transformed edge = " +
            residuals + " (n known in edges=" + nPointsInEdges1);
        
        int nBin = 0;
        
        int w = img2Masked.getWidth();
        int h = img2Masked.getHeight();
        while ((residuals > 0) && (img2Masked.getWidth() > 2)) {
            w >>= 1;
            h >>= 1;
            
            Image i2 = new Image(w, h);
            Image i1 = new Image(w, h);
            for (int i = 0; i < w; i++) {
                int x = i*2;
                for (int j = 0; j < h; j++) {
                    int y = j*2;
                    int sum = img2Masked.getG(x, y) 
                        + img2Masked.getG(x, y + 1)
                        + img2Masked.getG(x + 1, y)
                        + img2Masked.getG(x + 1, y + 1);
                    i2.setRGB(i, j, 0, sum, 0);
                }
            }
            img2Masked = i2;
            
            for (int i = 0; i < w; i++) {
                int x = i*2;
                for (int j = 0; j < h; j++) {
                    int y = j*2;
                    int sum = img1Transformed.getG(x, y) 
                        + img1Transformed.getG(x, y + 1)
                        + img1Transformed.getG(x + 1, y)
                        + img1Transformed.getG(x + 1, y + 1);
                    i1.setRGB(i, j, 0, sum, 0);
                }
            }
            img1Transformed = i1;
            
            residuals = 0;
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                    int pix2 = img2Masked.getG(i, j);
                    int pixT = img1Transformed.getG(i, j);
                    if (((pix2 != 0) && (pixT == 0)) 
                        || ((pix2 == 0) && (pixT != 0))) {
                        
                        residuals++;
                    }
                }
            }
            
            residuals /= 2;
            
            nBin++;
        
            log.info("number of pixels not overlapping in transformed edge = " +
                residuals + " for bin iteration " + nBin);
        
        }
       
        // scale is close, rotation is close, but translation is off by a couple
        // of pixels.
        
        // generous assert for now
        assertTrue(nBin < 7);
    }
    
    private void debugDisplay(PairIntArray[] transformedEdges,
        Image img, String rotDegrees) throws IOException {
        
        for (int i = 0; i < transformedEdges.length; i++) {
            PairIntArray edge = transformedEdges[i];
            debugAddCurveToImage(edge, img, 0, 255, 255, 255);
        }     
        
        String dirPath = ResourceFinder.findDirectory("bin");
        
        ImageIOHelper.writeOutputImage(
            dirPath + "/transformed_edges_" + rotDegrees + ".png", img);
    }
    
    private void debugAddCurveToImage(PairIntArray edge, Image input, 
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        for (int i = 0; i < edge.getN(); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                float xx = x + dx;
                if ((xx > -1) && (xx < (input.getWidth() - 1))) {
                    for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); 
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
            
            CurvatureScaleSpaceInflectionMapperTest test = 
                new CurvatureScaleSpaceInflectionMapperTest();
            
            test.testMap();
            
            //test.testMap2();
            
            //test.testTransformerGreyscaleImage();
            
            //test.testTransformerImage();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
