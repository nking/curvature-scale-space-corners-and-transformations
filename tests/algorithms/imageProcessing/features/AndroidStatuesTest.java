package algorithms.imageProcessing.features;

import algorithms.MultiArrayMergeSort;
import algorithms.imageProcessing.CannyEdgeFilterAdaptiveDeltaE2000;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.features.ObjectMatcher.Settings;
import algorithms.imageProcessing.matching.ObjectMatcherWrapper;
import algorithms.imageProcessing.segmentation.LabelToColorHelper;
import algorithms.imageProcessing.segmentation.NormalizedCuts;
import algorithms.imageProcessing.segmentation.SLICSuperPixels;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntPair;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
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

    public AndroidStatuesTest() {
    }
    
    public void testObjectFinder() throws Exception {
        boolean debug = true;
        runMatcher_gingerbreadman(debug);
        runMatcher_icecream(debug);
        runMatcher_cupcake(debug);
    }

    public void runMatcher_gingerbreadman(boolean debug) throws Exception {

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

                String filePath1
                    = ResourceFinder.findFileInTestResources(fileName1);

                List<CorrespondenceList> corresList = omw.find(filePath0, 
                    filePathMask0, filePath1, debugLabel);

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

    public void runMatcher_cupcake(boolean debug) throws Exception {
        
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
    
    public void runMatcher_icecream(boolean debug) throws Exception {

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
    
    public void estMkImgs() throws Exception {

        String fileName1 = "";

        for (int i = 0; i < 4; ++i) {

            switch(i) {
                case 0: {
                    fileName1
                        = "android_statues_01.jpg";
                    break;}
                case 1: {
                    fileName1 = "android_statues_02.jpg";
                    break;}
                case 2: {
                    fileName1 = "android_statues_03.jpg";
                    break;}
                case 3: {
                    fileName1 = "android_statues_04.jpg";
                    break;}
                default: {break;}
            }

            int idx = fileName1.lastIndexOf(".");
            String fileName1Root = fileName1.substring(0, idx);

            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            ImageExt img = ImageIOHelper.readImageExt(filePath1);

            int w1 = img.getWidth();
            int h1 = img.getHeight();

            int maxDimension = 512;
            int binFactor1 = (int) Math.ceil(Math.max((float) w1 / maxDimension,
                (float) h1 / maxDimension));

            ImageProcessor imageProcessor = new ImageProcessor();
            img = imageProcessor.binImage(img, binFactor1);

            {
                ImageExt imgCp = img.copyToImageExt();  
                CannyEdgeFilterAdaptiveDeltaE2000 canny = new 
                    CannyEdgeFilterAdaptiveDeltaE2000();
                //canny.setToDebug();
                canny.setToUseSingleThresholdIn2LayerFilter();
                canny.applyFilter(imgCp);
                MiscDebug.writeImage(canny.getFilterProducts().getGradientXY(),
                    "_GXY_" + fileName1Root);
            
            }
            
            ImageExt imgCp = img.copyToImageExt();

            int nClusters = 200;//100;
            //int clrNorm = 5;
            SLICSuperPixels slic = new SLICSuperPixels(img, nClusters);
            slic.calculate();
            int[] labels = slic.getLabels();
            ImageIOHelper.addAlternatingColorLabelsToRegion(
            //LabelToColorHelper.applyLabels(
                img, labels);
            MiscDebug.writeImage(img, "_slic_" + fileName1Root);

            NormalizedCuts normCuts = new NormalizedCuts();
            normCuts.setColorSpaceToHSV();
            int[] labels2 = normCuts.normalizedCut(imgCp, labels);
            labels = labels2;
            ImageIOHelper.addAlternatingColorLabelsToRegion(
                imgCp, labels);
            MiscDebug.writeImage(imgCp, "_norm_cuts_" + fileName1Root);


            //img = ImageIOHelper.readImageExt(filePath1);
            //img = imageProcessor.binImage(img, binFactor1);
            //MiscDebug.writeImage(img,  "_512_img_" + fileName1Root);
        }
    }

    public void estColorLayout() throws Exception {

        String[] fileNames = new String[]{
            "android_statues_01.jpg",
            "android_statues_02.jpg",
            "android_statues_03.jpg",
            "android_statues_04.jpg"
        };

        // paused here.  editing to use super pixels
        //   and a pattern of color aggregation
        //   w/ voronoi cells for comparison to model,
        //   then a partial shape matching algorithm.

        ImageProcessor imageProcessor = new ImageProcessor();

        ImageSegmentation imageSegmentation = new ImageSegmentation();

        List<ImageExt> images = new ArrayList<ImageExt>();

        for (int i = 0; i < fileNames.length; ++i) {

            String fileName = fileNames[i];
            String filePath = ResourceFinder.findFileInTestResources(fileName);
            String fileNameRoot = fileName.substring(0,
                fileName.lastIndexOf("."));
            ImageExt img = ImageIOHelper.readImageExt(filePath);
            images.add(img);

            /*
            wanting to look at color similarity of pixels
                holding known objects that have different
                orientation and lighting in different
                images.
                -- android
                -- ice cream
                -- eclair
                -- cupcake
            bigger goal is to use segmentation to find object
            outlines then identify the object in other images
            and follow that with detailed feature matching.

            the method has to work on objects that have changed
            position and possibly also lighting.

            the method may need some form of contour matching
            for pure sillouhette conditions
            but would need to allow for occlusion.

            deltaE for gingerbread man in the 4 images
            is at most 5, else 3.7 and 1.8.
            */
        }

    }

    private PairIntArray extractOrderedBoundary(ImageExt image) {
        return extractOrderedBoundary(image, SIGMA.TWO);
    }

    private PairIntArray extractOrderedBoundary(ImageExt image,
        SIGMA sigma) {

        GreyscaleImage img = image.copyToGreyscale();

        Set<PairInt> blob = new HashSet<PairInt>();
        for (int i = 0; i < img.getNPixels(); ++i) {
            if (img.getValue(i) > 0) {
                int x = img.getCol(i);
                int y = img.getRow(i);
                blob.add(new PairInt(x, y));
            }
        }

        ImageProcessor imageProcessor =
            new ImageProcessor();

        PairIntArray ordered =
            imageProcessor.extractSmoothedOrderedBoundary(
            blob, sigma, img.getWidth(), img.getHeight());

        return ordered;
    }

    private void sortByDecrSize(List<Set<PairInt>> clusterSets) {
        int n = clusterSets.size();
        int[] sizes = new int[n];
        int[] indexes = new int[n];
        for (int i = 0; i < n; ++i) {
            sizes[i] = clusterSets.get(i).size();
            indexes[i] = i;
        }
        MultiArrayMergeSort.sortByDecr(sizes, indexes);
        List<Set<PairInt>> out = new ArrayList<Set<PairInt>>();
        for (int i = 0; i < n; ++i) {
            int idx = indexes[i];
            out.add(clusterSets.get(idx));
        }
        clusterSets.clear();
        clusterSets.addAll(out);
    }

    private TIntList addIntersection(GreyscaleImage gradient,
        int[] labels) {

        assert(gradient.getNPixels() == labels.length);

        int maxLabel = MiscMath.findMax(labels) + 1;

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        int w = gradient.getWidth();
        int h = gradient.getHeight();

        TIntList change = new TIntArrayList();

        for (int i = 0; i < gradient.getNPixels(); ++i) {
            if (gradient.getValue(i) < 1) {
                continue;
            }
            int x = gradient.getCol(i);
            int y = gradient.getRow(i);
            int l0 = labels[i];
            int l1 = -1;
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if (x2 < 0 || y2 < 0 || (x2 > (w - 1))
                    || (y2 > (h - 1))) {
                    continue;
                }
                int j = gradient.getInternalIndex(x2, y2);
                int lt = labels[j];
                if (lt != l0) {
                    if (l1 == -1) {
                        l1 = lt;
                    }
                }
            }
            if (l1 != -1) {
                // gradient is on edge of superpixels
                change.add(i);
                System.out.println(
                    "x=" + x  + " y=" + y
                    + "  pixIdx=" + i);
            }
        }
        return change;
    }

    private int[] desegment(ImageExt img,
        TIntList gradSP, int[] labels, int[] labels2) {

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        int w = img.getWidth();
        int h = img.getHeight();

        TIntSet restore = new TIntHashSet();

        for (int i = 0; i < gradSP.size(); ++i) {
            int pixIdx = gradSP.get(i);

            int l0 = labels2[pixIdx];
            int l1 = -1;
            int x = img.getCol(pixIdx);
            int y = img.getRow(pixIdx);
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if (x2 < 0 || y2 < 0 || (x2 > (w - 1))
                    || (y2 > (h - 1))) {
                    continue;
                }
                int j = img.getInternalIndex(x2, y2);
                int lt = labels2[j];
                if (lt != l0) {
                    if (l1 == -1) {
                        l1 = lt;
                    }
                }
            }
            if (l1 == -1) {
                // these need boundaries restored
                ImageIOHelper.addPointToImage(x, y,
                    img, 1, 255, 0, 0);
                restore.add(labels[pixIdx]);
            }
        }
        MiscDebug.writeImage(img, "restore");

        if (restore.isEmpty()) {
            return labels2;
        }

        TIntObjectMap<Set<PairInt>> label1PointMap =
            LabelToColorHelper.extractLabelPoints(img,
                labels);

        int[] labels3 = Arrays.copyOf(labels2,
            labels2.length);

        int maxLabel = MiscMath.findMax(labels2);
        maxLabel++;
        TIntIterator iter = restore.iterator();
        while (iter.hasNext()) {
            int label = iter.next();
            Set<PairInt> pts = label1PointMap.get(label);
            for (PairInt pt : pts) {
                int x = pt.getX();
                int y = pt.getY();
                int pixIdx = img.getInternalIndex(x, y);
                labels3[pixIdx] = maxLabel;
            }
            maxLabel++;
        }
        return labels3;
    }

    public void estMatching() throws Exception {

        String fileName1, fileName2;

        FeatureMatcherSettings settings = new FeatureMatcherSettings();
        settings.setDebug(true);
        settings.setStartWithBinnedImages(true);
        settings.setToUse2ndDerivCorners();
        //for (int i = 0; i < 7; ++i) {
        for (int i = 0; i < 3; ++i) {
            switch(i) {
                case 0: {
                    fileName1 = "android_statues_02.jpg";
                    fileName2 = "android_statues_04.jpg";
                    //fileName1 = "android_statues_02_gingerbreadman.jpg";
                    //fileName2 = "android_statues_04_gingerbreadman.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 1: {
                    fileName1 = "android_statues_01.jpg";
                    fileName2 = "android_statues_03.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 2: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    fileName2 = "brown_lowe_2003_image2.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 3: {
                    fileName2 = "brown_lowe_2003_image1.jpg";
                    fileName1 = "brown_lowe_2003_image2.jpg";
                    break;
                }
                case 4: {
                    fileName1 = "venturi_mountain_j6_0001.png";
                    fileName2 = "venturi_mountain_j6_0010.png";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 5: {
                    fileName1 = "campus_010.jpg";
                    fileName2 = "campus_011.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 6: {
                    fileName1 = "merton_college_I_001.jpg";
                    fileName2 = "merton_college_I_002.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                default: {
                    fileName1 = "checkerboard_01.jpg";
                    fileName2 = "checkerboard_02.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
            }
            runCorrespondenceList(fileName1, fileName2, settings, false);
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

    private void populateClasses(Set<PairIntPair> similarClass,
        Set<PairIntPair> differentClass, String fileNameRoot) throws IOException {

        BufferedReader bReader = null;
        FileReader reader = null;

        String fileName = "label_" + fileNameRoot + "_coords.csv";

        String filePath = ResourceFinder.findFileInTestResources(fileName);

        try {
            reader = new FileReader(new File(filePath));

            bReader = new BufferedReader(reader);

            //read comment line and discard
            String line = bReader.readLine();
            line = bReader.readLine();

            while (line != null) {

                String[] items = line.split(",");
                if (items.length != 5) {
                    throw new IllegalStateException("Error while reading " +
                        fileName + " expecting 5 items in a line");
                }

                PairIntPair pp = new PairIntPair(
                    Integer.valueOf(items[0]).intValue(),
                    Integer.valueOf(items[1]).intValue(),
                    Integer.valueOf(items[2]).intValue(),
                    Integer.valueOf(items[3]).intValue());

                int classValue = Integer.valueOf(items[4]).intValue();

                if (classValue == 0) {
                    similarClass.add(pp);
                } else {
                    differentClass.add(pp);
                }

                line = bReader.readLine();
            }

        } catch (IOException e) {
            log.severe(e.getMessage());
        } finally {
            if (reader == null) {
                reader.close();
            }
            if (bReader == null) {
                bReader.close();
            }
        }
    }

    private void plot(PairIntArray p, int fn) throws Exception {

        float[] x = new float[p.getN()];
        float[] y = new float[p.getN()];

        for (int i = 0; i < x.length; ++i) {
            x[i] = p.getX(i);
            y[i] = p.getY(i);
        }

        float xMax = MiscMath.findMax(x) + 1;
        float yMax = MiscMath.findMax(y) + 1;

        PolygonAndPointPlotter plot = new PolygonAndPointPlotter();

        plot.addPlot(0, xMax, 0, yMax,
            x, y, x, y, "");

        plot.writeFile(fn);
    }
    
    private List<Set<PairInt>> extractNonZeros(ImageExt[] imgs) {
        
        List<Set<PairInt>> shape 
            = new ArrayList<Set<PairInt>>(imgs.length);
        
        for (int i = 0; i < imgs.length; ++i) {
            shape.add(extractNonZeros(imgs[i]));
        }
        
        return shape;
    }
    
    /**
     * extract all non-zero pixels
     * assuming all pixels with 0 value in r,g,b are mask pixels.
     * @param img
     * @return 
     */
    private Set<PairInt> extractNonZeros(ImageExt img) {
        
        Set<PairInt> set = new HashSet<PairInt>();
        for (int i = 0; i < img.getNPixels(); ++i) {
            if (img.getR(i) != 0 ||  img.getG(i) != 0 || img.getB(i) != 0) {
                set.add(new PairInt(img.getCol(i), img.getRow(i)));
            }
        }
        
        return set;
    }
        
    private ImageExt[] loadMaskedIcecream() throws IOException {
        
        // icecream images are x:0 to 255, y=0 to 63
        int d = 64;
        
        String filePath = ResourceFinder.findFileInTestResources(
            "android_statues_objects.png");
        ImageExt img = ImageIOHelper.readImageExt(filePath);

        ImageExt[] out = new ImageExt[3];
        
        out[0] = (ImageExt) img.copySubImage(0, d, 0, d);
        out[1] = (ImageExt) img.copySubImage(d, 2*d, 0, d);
        out[2] = (ImageExt) img.copySubImage(3*d, 4*d, 0, d);
        
        return out;
    }
    
    private ImageExt[] loadMaskedCupcake() throws IOException {
        
        // cupckes images are x:0 to 255, y=64 to 128
        int d = 64;
        
        String filePath = ResourceFinder.findFileInTestResources(
            "android_statues_objects.png");
        ImageExt img = ImageIOHelper.readImageExt(filePath);

        ImageExt[] out = new ImageExt[3];
        
        out[0] = (ImageExt) img.copySubImage(0, d, d, 2*d);
        out[1] = (ImageExt) img.copySubImage(d, 2*d, d, 2*d);
        out[2] = (ImageExt) img.copySubImage(3*d, 4*d, d, 2*d);
        
        return out;
    }
    
    private ImageExt[] loadMaskedGingerBreadMan() throws IOException {
        
        // cupckes images are x:0 to 255, y=128 to 192
        int d = 64;
        
        String filePath = ResourceFinder.findFileInTestResources(
            "android_statues_objects.png");
        ImageExt img = ImageIOHelper.readImageExt(filePath);

        ImageExt[] out = new ImageExt[4];
        
        out[0] = (ImageExt) img.copySubImage(0,     d, 2*d, 3*d);
        out[1] = (ImageExt) img.copySubImage(d,   2*d, 2*d, 3*d);
        out[2] = (ImageExt) img.copySubImage(2*d, 3*d, 2*d, 3*d);
        out[3] = (ImageExt) img.copySubImage(3*d, 4*d, 2*d, 3*d);
        
        return out;
    }
    
    private ImageExt[] loadMaskedEuclair() throws IOException {
        
        // cupckes images are x:0 to 255, y=192 to 256
        int d = 64;
        
        String filePath = ResourceFinder.findFileInTestResources(
            "android_statues_objects.png");
        ImageExt img = ImageIOHelper.readImageExt(filePath);

        ImageExt[] out = new ImageExt[4];
        
        out[0] = (ImageExt) img.copySubImage(0,     d, 3*d, 4*d);
        out[1] = (ImageExt) img.copySubImage(d,   2*d, 3*d, 4*d);
        out[2] = (ImageExt) img.copySubImage(2*d, 3*d, 3*d, 4*d);
        out[3] = (ImageExt) img.copySubImage(3*d, 4*d, 3*d, 4*d);
        
        return out;
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
                p2.getX(), p2.getY(), 0);
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
