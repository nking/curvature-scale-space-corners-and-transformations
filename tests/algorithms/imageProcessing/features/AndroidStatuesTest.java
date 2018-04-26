package algorithms.imageProcessing.features;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.matching.ObjectMatcherWrapper;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class AndroidStatuesTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    private boolean debug = false;
    
    public AndroidStatuesTest() {
    }
    
    public void testObjectFinder() throws Exception {
        runMatcher_gingerbreadman();
        runMatcher_icecream();
        runMatcher_cupcake();
    }

    public void runMatcher_gingerbreadman() throws Exception {

        String[] fileNames0 = new String[]{
        //    "android_statues_03_sz1",
            "android_statues_03_sz3"
        };
        String[] fileNames1 = new String[]{
            "android_statues_01.jpg",
            "android_statues_02.jpg", 
            "android_statues_03.jpg",
            "android_statues_04.jpg", 
        };
        
        /* 
        ROUGH COORDS
        gbman template sz3    
            center
              60. 72  (sz1 has center 32,45)
        
        andr01 (scale 95/23 = 4.1   5:0   <=== 3rd
              78, 52
        
        andr02 (scale 95/70 = 1.35  1:0
            181, 57
            (seg: 175, 57)
        
        andr03 (scale 95/61 = 1.6   1:0
           38, 76
        
        andr04 (scale 95/46 = 2.1   2:0
           196, 87
        */
          
        int fn0 = 0;
        for (String fileNameRoot0 : fileNames0) {               
            fn0++;
           
            String filePathMask0 = ResourceFinder
                .findFileInTestResources(fileNameRoot0 + "_mask.png");
            String filePath0 = ResourceFinder
                .findFileInTestResources(fileNameRoot0 + ".jpg");
          
            for (int fIdx = 0; fIdx < fileNames1.length; ++fIdx) {               
        
                String fileName1 = fileNames1[fIdx];
                                
                String debugLabel = "gbman_" + fIdx + "_" + fn0;
                
                ObjectMatcherWrapper omw = new ObjectMatcherWrapper();

                if (debug) {
                    omw.setToDebug();
                }
                
                String filePath1
                    = ResourceFinder.findFileInTestResources(fileName1);

                List<CorrespondenceList> corresList = omw.find(filePath0, 
                    filePathMask0, filePath1, debugLabel);
                
                //omw.setToDebug();

                if (corresList == null || corresList.isEmpty()) {
                    return;
                }
                       
                plotCorrespondence(corresList, 
                    omw.getTemplateImage()[1].copyToImageExt(),
                    omw.getSearchImage().copyToImageExt(),
                    debugLabel);
            }
        }
    }

    public void runMatcher_cupcake() throws Exception {
        
        String[] fileNames0 = new String[]{
            "android_statues_04.jpg", 
            "android_statues_04_cupcake_mask.png"};

        String[] fileNames1 = new String[]{
            "android_statues_01.jpg",
            "android_statues_02.jpg",  
            "android_statues_04.jpg"
        };
        
        /*  ROUGH COORDS
        cupcake template
            top      bottom    approx center
          48, 32    44, 59         51, 54
        
        andr01  (scale 61/14 = 4.35  5:0    <==== not found
          15,48                                   needs to use template0 segmentation
        (seg: 14, 56)                             -> best cost for not found 0.44
        
        andr02  (scale 61/27 = 2.25  2:0
            52, 50                     
        (seg: 58, 60)
        
        andr04  (scale 61/58  = 1   0:0
         92, 74                   91, 103                
        */
        
        String filePath0 = ResourceFinder.findFileInTestResources(fileNames0[0]);
        String filePath0Mask = ResourceFinder.findFileInTestResources(fileNames0[1]);
        
        for (int fIdx = 0; fIdx < fileNames1.length; ++fIdx) {
            
            String fileName1 = fileNames1[fIdx];             

            String filePath1 = ResourceFinder.findFileInTestResources(
                fileName1);
            ImageExt img = ImageIOHelper.readImageExt(filePath1);
            
            if (fileName1.equals("android_statues_01.jpg")) {
                //A LOOK AT WHETHER zoom in helps find difficult shadowed patterns
                img = (ImageExt) img.copySubImage(0, img.getWidth()/5, 0, img.getHeight());
            } else if (fileName1.equals("android_statues_02.jpg")) {
                
            }

            String debugLabel = "cupcake_" + fIdx;
                
            ObjectMatcherWrapper omw = new ObjectMatcherWrapper();

            if (debug) {
                omw.setToDebug();
            }
            
            Set<PairInt> shape0 = new HashSet<PairInt>();

            ImageExt[] imgs0 = omw.maskAndBin2(filePath0, filePath0Mask,
                shape0);
        
            List<CorrespondenceList> corresList = omw.find(imgs0, shape0, img, 
                debugLabel);

            if (corresList == null || corresList.isEmpty()) {
                return;
            }

            plotCorrespondence(corresList, 
                omw.getTemplateImage()[1].copyToImageExt(),
                omw.getSearchImage().copyToImageExt(),
                debugLabel);
        }
    }
    
    public void runMatcher_icecream() throws Exception {

        String[] fileNames0 = new String[]{
            "android_statues_04.jpg",
            "android_statues_04_icecream_mask.png",
        };
        
        String[] fileNames1 = new String[]{
            "android_statues_01.jpg",  
            "android_statues_02.jpg",
            "android_statues_04.jpg", // descr are fine
        };
        
        /*  ROUGH COORDS
        icecream template
            top      bottom    approx center
           40, 40    38, 60       41, 51
           37,38
        
        andr01  (scale 53/17 = 3.1    3:0
                                  49, 53
        
        andr02  (scale 53/36 = 1.5    1:0
                                 108, 55,   98,63
        
        andr04  (scale 53/50 = 1      0:0
                                 153, 96
        */

        String filePath0 = ResourceFinder.findFileInTestResources(fileNames0[0]);
        String filePath0Mask = ResourceFinder.findFileInTestResources(fileNames0[1]);
        
        for (int fIdx = 0; fIdx < fileNames1.length; ++fIdx) {               

            String fileName1 = fileNames1[fIdx];
            
            String filePath1 = ResourceFinder.findFileInTestResources(
                fileName1);
            ImageExt img = ImageIOHelper.readImageExt(filePath1);
            
            if (fileName1.equals("android_statues_01.jpg")) {
                //A LOOK AT WHETHER zoom in helps find difficult shadowed patterns
                img = (ImageExt) img.copySubImage(0, img.getWidth()/5, 0, img.getHeight());
            } else if (fileName1.equals("android_statues_02.jpg")) {
                
            }

            String debugLabel = "icecream_" + fIdx;
                
            ObjectMatcherWrapper omw = new ObjectMatcherWrapper();
            if (debug) {
                omw.setToDebug();
            }

            List<CorrespondenceList> corresList = omw.find(filePath0, 
                filePath0Mask, filePath1, debugLabel);

            if (corresList == null || corresList.isEmpty()) {
                return;
            }

            plotCorrespondence(corresList, 
                omw.getTemplateImage()[1].copyToImageExt(),
                omw.getSearchImage().copyToImageExt(),
                debugLabel);
            
        }
    }
   
    public void estRot90() throws Exception {

        String fileName1, fileName2;

        FeatureMatcherSettings settings = new FeatureMatcherSettings();
        settings.setDebug(true);
        settings.setStartWithBinnedImages(true);

        for (int i = 0; i < 5; ++i) {
            fileName1 = null;
            fileName2 = null;
            switch(i) {
                case 0: {
                    fileName1 = "android_statues_02.jpg";
                    fileName2 = "android_statues_04.jpg";
                    //fileName1 = "campus_010.jpg";
                    //fileName2 = "campus_011.jpg";
                    settings.setUseNormalizedFeatures(true);
                    settings.setToUse2ndDerivCorners();
                    break;
                }
                default: {
                    fileName1 = "android_statues_01.jpg";
                    fileName2 = "android_statues_03.jpg";
                    settings.setUseNormalizedFeatures(true);
                    settings.setToUse2ndDerivCorners();
                    break;
                }
            }
            runCorrespondenceList(fileName1, fileName2, settings, true);
        }
    }

    private void runCorrespondenceList(String fileName1, String fileName2,
        FeatureMatcherSettings settings, boolean rotateBy90) throws Exception {

        if (fileName1 == null) {
            return;
        }

        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);
        idx = fileName2.lastIndexOf(".");
        String fileName2Root = fileName2.substring(0, idx);

        settings.setDebugTag(fileName1Root);

        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        ImageProcessor imageProcessor = new ImageProcessor();
        ImageSegmentation imageSegmentation = new ImageSegmentation();

        int w1 = img1.getWidth();
        int h1 = img1.getHeight();
        int w2 = img2.getWidth();
        int h2 = img2.getHeight();

        int maxDimension = 350;
        int binFactor1 = (int) Math.ceil(Math.max((float)w1/maxDimension,
            (float)h1/ maxDimension));
        int binFactor2 = (int) Math.ceil(Math.max((float)w2/maxDimension,
            (float)h2/ maxDimension));

        ImageExt img1Binned = imageProcessor.binImage(img1, binFactor1);
        ImageExt img2Binned = imageProcessor.binImage(img2, binFactor2);

        if (rotateBy90) {
            TransformationParameters params90 = new TransformationParameters();
            params90.setRotationInDegrees(90);
            params90.setOriginX(0);
            params90.setOriginY(0);
            params90.setTranslationX(0);
            params90.setTranslationY(img1.getWidth() - 1);
            Transformer transformer = new Transformer();
            img1 = (ImageExt) transformer.applyTransformation(img1,
                params90, img1.getHeight(), img1.getWidth());

            /*
            MatchedPointsTransformationCalculator tc =
                new MatchedPointsTransformationCalculator();

            TransformationParameters revParams = tc.swapReferenceFrames(params90);
            transformer.transformToOrigin(0, 0, revParams);
            revParams.setTranslationX(revParams.getTranslationX() + -74);
            revParams.setTranslationY(revParams.getTranslationY() + -0);
            revParams.setRotationInDegrees(revParams.getRotationInDegrees() - 0);
            log.info("revParams: " + revParams.toString());

            ImageExt img1RevTr = img1.copyToImageExt();
            img1RevTr = (ImageExt) transformer.applyTransformation(img1RevTr,
                revParams, img1RevTr.getHeight(), img1RevTr.getWidth());
            MiscDebug.writeImage(img1RevTr, "rot90_rev_trans");
            */
        }

        log.info("fileName1Root=" + fileName1Root);

        /*EpipolarColorSegmentedSolver solver = new EpipolarColorSegmentedSolver(img1, img2, settings);

        boolean solved = solver.solve();

        assertTrue(solved);

        //MiscDebug.writeImagesInAlternatingColor(img1, img2, stats,
        //    fileName1Root + "_matched_non_euclid", 2);
        */
    }

     public static void main(String[] args) {

        try {
            AndroidStatuesTest test = new AndroidStatuesTest();
            //test.test0();
            //test.testRot90();

        } catch(Exception e) {
            e.printStackTrace();
            System.err.println(e.getMessage());
            fail(e.getMessage());
        }
    }
 
    private void plotCorrespondence(List<CorrespondenceList> corresList, 
       ImageExt templateImage, ImageExt searchImage, 
       String debugLabel) throws IOException {
        
        //settings.setToExcludeColorFilter();

        CorrespondencePlotter plotter = new CorrespondencePlotter(templateImage,
            searchImage); 
        CorrespondenceList cor = corresList.get(0);
        for (int ii = 0; ii < cor.getPoints1().size(); ++ii) {
            PairInt p1 = cor.getPoints1().get(ii);
            PairInt p2 = cor.getPoints2().get(ii);

            //System.out.println("orb matched: " + p1 + " " + p2);
            //if (p2.getX() > 160)
            plotter.drawLineInAlternatingColors(p1.getX(), p1.getY(), 
                p2.getX(), p2.getY(), 1);
        }
        
        if (corresList.size() > 1) {
            cor = corresList.get(1);
            for (int ii = 0; ii < cor.getPoints1().size(); ++ii) {
                PairInt p1 = cor.getPoints1().get(ii);
                PairInt p2 = cor.getPoints2().get(ii);

                //System.out.println("orb matched: " + p1 + " " + p2);
                //if (p2.getX() > 160)
                plotter.drawDashedLine(p1.getX(), p1.getY(),
                    p2.getX(), p2.getY(), 255, 200, 200, 0, 7);
            }
        }

        plotter.writeImage("_orb_corres_final_" + 
            debugLabel);
        System.out.println(cor.getPoints1().size() + 
            " matches " + debugLabel);
        //MiscDebug.writeImage(img11, "_orb_matched_" + str
        //    + "_" + fileName1Root);
    }
  
}
