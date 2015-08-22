package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import algorithms.util.PairIntArray;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class CurvatureScaleSpaceInflectionMapperTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());

    public CurvatureScaleSpaceInflectionMapperTest() {
    }
       
    public void testTransformerGreyscaleImage() throws Exception {
        
        String[] rotDegrees = new String[]{"225"};
        
        String fileName1 = "closed_curve.png";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        String fileName2 = "closed_curve_translate_scale_rotate225.png";

        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);;

        CurvatureScaleSpaceInflectionMapper mapper = new 
            CurvatureScaleSpaceInflectionMapper(img1, img2);

        mapper.useLineDrawingLineMode();
        
        mapper.useDebugMode();

        TransformationParameters transformationParams = 
            mapper.createEuclideanTransformation();

        assertNotNull(transformationParams);

        double rotDeg = transformationParams.getRotationInDegrees();
        
        assertTrue(Math.abs(rotDeg - 136.f) < 10.f);

        double scale = transformationParams.getScale();

        assertTrue(Math.abs(scale - 1.3) < 0.15);

        List<PairIntArray> edges1 = 
            mapper.getEdges1InOriginalReferenceFrame();
        int nPointsInEdges1 = 0;
        for (int i = 0; i < edges1.size(); i++) {
            nPointsInEdges1 += edges1.get(i).getN();
        }

        Transformer transformer = new Transformer();
        
        img2 = ImageIOHelper.readImageExt(filePath2);
        GreyscaleImage img2Masked = img2.copyToGreyscale();
        for (int i = 180; i < img2Masked.getWidth(); i++) {
            for (int j = 0; j < img2Masked.getHeight(); j++) {
                img2Masked.setValue(i, j, 0);
            }
        }
        
        img1 = ImageIOHelper.readImageExt(filePath1);
        
        GreyscaleImage img1Transformed = transformer.applyTransformation(
            img1.copyToGreyscale(), transformationParams, img2.getWidth(), img2.getHeight());
        
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
    
    public void testTransformerImage() throws Exception {
        
        String[] rotDegrees = new String[]{"225"};
        
        String fileName1 = "closed_curve.png";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        String fileName2 = "closed_curve_translate_scale_rotate225.png";

        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        CurvatureScaleSpaceInflectionMapper mapper = new 
            CurvatureScaleSpaceInflectionMapper(img1, img2);

        mapper.useLineDrawingLineMode();

        TransformationParameters transformationParams = 
            mapper.createEuclideanTransformation();

        assertNotNull(transformationParams);

        double rotDeg = transformationParams.getRotationInDegrees();
        
        assertTrue(Math.abs(rotDeg - 136.f) < 10.f);

        double scale = transformationParams.getScale();

        assertTrue(Math.abs(scale - 1.3) < 0.15);

        List<PairIntArray> edges1 = 
            mapper.getEdges1InOriginalReferenceFrame();
        int nPointsInEdges1 = 0;
        for (int i = 0; i < edges1.size(); i++) {
            nPointsInEdges1 += edges1.get(i).getN();
        }

        Transformer transformer = new Transformer();
        
        img2 = ImageIOHelper.readImageExt(filePath2);
            
        Image img2Masked = img2.copyImage();
        for (int i = 180; i < img2Masked.getWidth(); i++) {
            for (int j = 0; j < img2Masked.getHeight(); j++) {
                img2Masked.setRGB(i, j, 0, 0, 0);
            }
        }
        
        img1 = ImageIOHelper.readImageExt(filePath1);
        
        Image img1Transformed = transformer.applyTransformation(
            img1.copyImage(), transformationParams, 
            img2Masked.getWidth(), img2Masked.getHeight());
        
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
    
    public static void main(String[] args) {
        
        try {
            
            CurvatureScaleSpaceInflectionMapperTest test = 
                new CurvatureScaleSpaceInflectionMapperTest();
                                                
            test.testTransformerGreyscaleImage();
            
            test.testTransformerImage();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
