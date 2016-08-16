package algorithms.imageProcessing.features;

import algorithms.compGeometry.PerimeterFinder2;
import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.DFSContiguousIntValueFinder;
import algorithms.imageProcessing.DFSContiguousValueFinder;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.ImageSegmentation.DecimatedData;
import algorithms.imageProcessing.matching.PartialShapeMatcher;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.SegmentationMergeThreshold;
import algorithms.imageProcessing.matching.Sequence;
import algorithms.imageProcessing.matching.Sequences;
import algorithms.imageProcessing.segmentation.LabelToColorHelper;
import algorithms.imageProcessing.segmentation.SLICSuperPixels;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.util.GroupAverageColors;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntPair;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.QuadInt;
import algorithms.util.ResourceFinder;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;
import org.ejml.data.Complex64F;
import org.ejml.simple.SimpleEVD;
import org.ejml.simple.SimpleMatrix;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class AndroidStatuesTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public AndroidStatuesTest() {
    }

     public void est0() throws Exception {

        String fileName1 = "";
        SegmentationMergeThreshold mt = SegmentationMergeThreshold.DEFAULT;

        //11, 12, 21, 30
        //for (int i = 11; i < 13; ++i) {
        for (int i = 0; i < 37; ++i) {
            
            /*
            if ( !( (i==8)||(i==12)||(i==21)||(i==30))) {
                continue;
            }
            */
            mt = SegmentationMergeThreshold.DEFAULT;
            switch(i) {
                case 0: {
                    fileName1 = "android_statues_01.jpg";
                    break;
                }
                case 1: {
                    fileName1 = "android_statues_02.jpg";
                    break;
                }
                case 2: {
                    fileName1 = "android_statues_03.jpg";
                    break;
                }
                case 3: {
                    fileName1 = "android_statues_04.jpg";
                    break;
                }
                case 4: {
                    fileName1 = "seattle.jpg";
                    break;
                }
                case 5: {
                    fileName1 = "stonehenge.jpg";
                    mt = SegmentationMergeThreshold.HIGH_CONTRAST_ONLY;
                    break;
                }
                case 6: {
                    fileName1 = "cloudy_san_jose.jpg";
                    mt = SegmentationMergeThreshold.HIGH_CONTRAST_ONLY;
                    break;
                }
                case 7: {
                    fileName1 = "patagonia_snowy_foreground.jpg";
                    break;
                }
                case 8: {
                    fileName1 = "mt_rainier_snowy_field.jpg";
                    mt = SegmentationMergeThreshold.EXTREMELY_LOW_CONTRAST;
                    break;
                }
                case 9: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    break;
                }
                case 10: {
                    fileName1 = "brown_lowe_2003_image2.jpg";
                    break;
                }
                case 11: {
                    fileName1 = "venturi_mountain_j6_0001.png";
                    //mt = SegmentationMergeThreshold.EXTREMELY_LOW_CONTRAST;
                    break;
                }
                case 12: {
                    fileName1 = "venturi_mountain_j6_0010.png";
                    //mt = SegmentationMergeThreshold.EXTREMELY_LOW_CONTRAST;
                    break;
                }
                case 13: {
                    fileName1 = "campus_010.jpg";
                    break;
                }
                case 14: {
                    fileName1 = "campus_011.jpg";
                    break;
                }
                case 15: {
                    fileName1 = "merton_college_I_001.jpg";
                    break;
                }
                case 16: {
                    fileName1 = "merton_college_I_002.jpg";
                    break;
                }
                case 17: {
                    fileName1 = "arches.jpg";
                    break;
                }
                case 18: {
                    fileName1 = "stinson_beach.jpg";
                    //mt = SegmentationMergeThreshold.HIGH_CONTRAST_ONLY;
                    break;
                }                
                case 19: {
                    fileName1 = "norwegian_mtn_range.jpg";
                    //mt = SegmentationMergeThreshold.HIGH_CONTRAST_ONLY;
                    break;
                }
                case 20: {
                    fileName1 = "halfdome.jpg";
                    break;
                }
                case 21: {
                    fileName1 = "halfdome2.jpg";
                    mt = SegmentationMergeThreshold.EXTREMELY_LOW_CONTRAST;
                    break;
                }
                case 22: {
                    fileName1 = "halfdome3.jpg";
                    mt = SegmentationMergeThreshold.HIGH_CONTRAST_ONLY;
                    break;
                }
                case 23: {
                    fileName1 = "costa_rica.jpg";
                    //mt = SegmentationMergeThreshold.HIGH_CONTRAST_ONLY;
                    break;
                }
                case 24: {
                    fileName1 = "new-mexico-sunrise_w725_h490.jpg";
                    mt = SegmentationMergeThreshold.HIGH_CONTRAST_ONLY;
                    break;
                }
                case 25: {
                    fileName1 = "arizona-sunrise-1342919937GHz.jpg";
                    mt = SegmentationMergeThreshold.HIGH_CONTRAST_ONLY;
                    break;
                }
                case 26: {
                    fileName1 = "sky_with_rainbow.jpg";
                    break;
                }
                case 27: {
                    fileName1 = "sky_with_rainbow2.jpg";
                    break;
                }
                case 28: {
                    fileName1 = "books_illum3_v0_695x555.png";
                    break;
                }
                case 29: {
                    fileName1 = "books_illum3_v6_695x555.png";
                    break;
                }
                case 30: {
                    fileName1 = "klein_matterhorn_snowy_foreground.jpg";
                    mt = SegmentationMergeThreshold.EXTREMELY_LOW_CONTRAST;
                    break;
                }
                case 31: {
                    fileName1 = "30.jpg";
                    break;
                }
                case 32: {
                    fileName1 = "arches_sun_01.jpg";
                    mt = SegmentationMergeThreshold.HIGH_CONTRAST_ONLY;
                    break;
                }
                case 33: {
                    fileName1 = "stlouis_arch.jpg";
                    mt = SegmentationMergeThreshold.HIGH_CONTRAST_ONLY;
                    break;
                }
                case 34: {
                    fileName1 = "contrail.jpg";
                    mt = SegmentationMergeThreshold.HIGH_CONTRAST_ONLY;
                    break;
                }
                case 35: {
                    fileName1 = "checkerboard_01.jpg";
                    //mt = SegmentationMergeThreshold.HIGH_CONTRAST_ONLY;
                    break;
                }
                default: {
                    fileName1 = "checkerboard_02.jpg";
                    //mt = SegmentationMergeThreshold.HIGH_CONTRAST_ONLY;
                    break;
                }
            }
       
            int idx = fileName1.lastIndexOf(".");
            String fileName1Root = fileName1.substring(0, idx);

            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            ImageExt img = ImageIOHelper.readImageExt(filePath1);

            ImageProcessor imageProcessor = new ImageProcessor();
            ImageSegmentation imageSegmentation = new ImageSegmentation();

            int w1 = img.getWidth();
            int h1 = img.getHeight();

            int maxDimension = 256;
            int binFactor1 = (int) Math.ceil(Math.max((float) w1 / maxDimension,
                (float) h1 / maxDimension));

            img = imageProcessor.binImage(img, binFactor1);

            List<Set<PairInt>> segmentedCellList
                = imageSegmentation.createColorEdgeSegmentation(
                    img, fileName1Root);
            
            MiscDebug.writeAlternatingColor(img, 
                segmentedCellList, "_final_" + fileName1Root);
                                   
            /*
            int nClusters = 200;
            //int clrNorm = 5;
            
            SLICSuperPixels slic 
                = new SLICSuperPixels(img, nClusters);

            slic.calculate();

            int[] labels = slic.getLabels();

            ImageIOHelper.addAlternatingColorLabelsToRegion(img, labels);
            MiscDebug.writeImage(img,  "_slic_" + fileName1Root);
            
            img = ImageIOHelper.readImageExt(filePath1);
            img = imageProcessor.binImage(img, binFactor1);
            LabelToColorHelper.applyLabels(img, labels);
            MiscDebug.writeImage(img,  "_slic_img_" + fileName1Root);
            */
        }
    }
     
    public void test1() throws Exception {

        String fileName0 
            = "android_statues_03_sz1_mask_small.png";
        int idx = fileName0.lastIndexOf(".");
        String fileName0Root = fileName0.substring(0, idx);
        String filePath0 = ResourceFinder
            .findFileInTestResources(fileName0);
        ImageExt img0 = ImageIOHelper.readImageExt(filePath0);

        PairIntArray p = extractOrderedBoundary(img0);
        plot(p, 100);
        
        String fileName1 = "";

        /*
        (required tol=0.25 or 0.35)
        0: offset approx 0  NOTE: p.n=187 q.n=231, needed topJ=5
        1: offset approx 0  NOTE: p.n=187 q.n=144, needs topJ=15
                         correct answer is trck 3,
                         sol= looking at tolerance and topK
                         to make it trck 0.
                         note that havne;t impl the append
                         be min offset near offse too/
                         nor did add back the higher error sequences/
        2: 
        3: offset approx 0  NOTE: p.n=187 q.n=181, needs topK=5
           
        -- the rigid model should solve this more easily
           as it will have an eval stage for transformation.
        -- meanwhile, could try two things with current model:
            -- lower threshold in search of diff matrix then
               filter to create a low and high threshold list
               and solve both, keeping best.
        -- when pattern is robust, witl consider a larger dp
           (which will make smaller number of points, hence
           speed up algorithm...larger dp's lose curvatur info though
        */
        for (int i = 0; i < 4; ++i) {
            
            switch(i) {
                case 0: {
                    fileName1 
                        = "android_statues_01_sz1_mask_small.png";
                    break;
                }
                case 1: {
                    fileName1 = "android_statues_02_sz1_mask_small.png";
                    break;
                }
                case 2: {
                    fileName1 = "android_statues_03_sz1_mask_small.png";
                    break;
                }
                case 3: {
                    fileName1 = "android_statues_04_sz1_mask_small.png";
                    break;
                }
                default: {
                    break;
                }
            }
       
            idx = fileName1.lastIndexOf(".");
            String fileName1Root = fileName1.substring(0, idx);

            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            ImageExt img = ImageIOHelper.readImageExt(filePath1);

            /*
            int w1 = img.getWidth();
            int h1 = img.getHeight();

            int maxDimension = 256;
            int binFactor1 = (int) Math.ceil(Math.max((float) w1 / maxDimension,
                (float) h1 / maxDimension));

            img = imageProcessor.binImage(img, binFactor1);
            */
            
            /*
            (1) make template objects to search for in other
                images. 
                - a color version for each image binned scale
                - a contour of the object at each image binned scale
            (2) use super pixels on segmented images.
            (3) for segment centroids, build k nearest neighbors.
            (4) search image for segments that match the search object
                in color.
                these form the core of what will be aggregated search.
            (5) create best pattern to search the candidate
                segmented cells.
                - determine different combinations of each with
                  its neighbors that might create the needed
                  area?  
                  - exclude adding neighbors that don't match any
                    color segment in the search object?
                    (might not be a good idea considering resolution
                     of super pixels)
                - smooth the bounds?
            ===> before that, extract the shapes 
                 aggregated (and smoothed if needed),
                 to make tests for the partial shape
                 matcher to adapt it for use.
                 (might need to add descriptors, the
                 color or gradient too, hopefully not)
            */
            
            PairIntArray q = extractOrderedBoundary(img);
            plot(q, (i+1)*100 + 1);

            System.out.println("matching " + fileName0Root
            + " to " + fileName1Root + " (" + p.getN() 
            + " points to " + q.getN() + " points");
            
            PartialShapeMatcher matcher = 
                new PartialShapeMatcher();
            matcher.overrideSamplingDistance(2);

            Sequences sequences = matcher.match(p, q);
            
            assertNotNull(sequences);
            
            System.out.println(
                "RESULTS=" + 
                fileName1Root + " : " +
                sequences.toString());
            
            CorrespondencePlotter plotter = new
                CorrespondencePlotter(p, q);
            
            for (Sequence s : sequences.getSequences()) {
                int len = s.getStopIdx2() - s.getStartIdx2()
                    + 1;
                for (int offset = 0; offset <= len; ++offset) {
                    int idx1 = s.getStartIdx1() + offset;
                    if (idx1 > (p.getN() - 1)) {
                        idx1 -= p.getN();
                    }
                    int idx2 = s.getStartIdx2() + offset;
                    if (idx2 > (q.getN() - 1)) {
                        idx2 -= q.getN();
                    }
                    int x1 = p.getX(idx1);
                    int y1 = p.getY(idx1);
                    int x2 = q.getX(idx2);
                    int y2 = q.getY(idx2);
                    //System.out.println(String.format(
                    //"(%d, %d) <=> (%d, %d)", x1, y1, x2, y2));
                
                    if ((offset % 5) == 0) {
                        plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 
                            0);
                    }
                }
                String filePath = plotter.writeImage("_" +
                    fileName1Root + "_corres");
            }
            
            /*
            MiscDebug.writeImage(img,  "_img_" + fileName1Root);
            int nClusters = 200;//100;
            //int clrNorm = 5;            
            SLICSuperPixels slic 
                = new SLICSuperPixels(img, nClusters);
            slic.calculate();
            int[] labels = slic.getLabels();
            ImageIOHelper.addAlternatingColorLabelsToRegion(img, labels);
            MiscDebug.writeImage(img,  "_slic_" + fileName1Root);
            img = ImageIOHelper.readImageExt(filePath1);
            //img = imageProcessor.binImage(img, binFactor1);
            LabelToColorHelper.applyLabels(img, labels);
            MiscDebug.writeImage(img,  "_slic_img_" + fileName1Root);
            */
            
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

        GreyscaleImage img = image.copyToGreyscale();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(img, 
            SIGMA.TWO);
        
        Set<PairInt> blob = new HashSet<PairInt>();
        for (int i = 0; i < img.getNPixels(); ++i) {
            if (img.getValue(i) > 0) {
                int x = img.getCol(i);
                int y = img.getRow(i);
                img.setValue(x, y, 255);
                blob.add(new PairInt(x, y));
            }
        }
        
        PerimeterFinder2 finder = new PerimeterFinder2();
        PairIntArray ordered = finder.extractOrderedBorder(blob);
    
        return ordered;
    }
    
    private class DeltaESim implements Comparable<DeltaESim> {

        private int x1;
        private int y1;
        private int x2;
        private int y2;
        private double deltaE;
        
        public DeltaESim(GroupAverageColors avg1,
            GroupAverageColors avg2) {
            this.x1 = avg1.getXCen();
            this.y1 = avg1.getYCen();
            this.x2 = avg2.getXCen();
            this.y2 = avg2.getYCen();
            CIEChromaticity cieC = new CIEChromaticity();
            this.deltaE = 
                Math.abs(cieC.calcDeltaECIE2000(
                avg1.getCIEL(), avg1.getCIEA(), avg1.getCIEB(), 
                avg2.getCIEL(), avg2.getCIEA(), avg2.getCIEB()));
        }
        
        @Override
        public int compareTo(DeltaESim other) {
            if (deltaE < other.deltaE) {
                return -1;
            } else if (deltaE > other.deltaE) {
                return 1;
            }
            return 0;
        }
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

        RotatedOffsets rotatedOffsets = RotatedOffsets.getInstance();

        log.info("fileName1Root=" + fileName1Root);

        EpipolarColorSegmentedSolver solver = new EpipolarColorSegmentedSolver(img1, img2, settings);

        boolean solved = solver.solve();

        assertTrue(solved);

        //MiscDebug.writeImagesInAlternatingColor(img1, img2, stats,
        //    fileName1Root + "_matched_non_euclid", 2);

    }

     public static void main(String[] args) {

        try {
            AndroidStatuesTest test = new AndroidStatuesTest();
            //test.test0();
            //test.testRot90();

        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            fail(e.getMessage());
        }
    }

    public void estColor() throws Exception {

        //String fileName1 = "android_statues_02.jpg";
        //String fileName2 = "android_statues_04.jpg";
        String fileName1 = "android_statues_02_gingerbreadman.jpg";
        String fileName2 = "android_statues_04_gingerbreadman.jpg";
        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);
        idx = fileName2.lastIndexOf(".");
        String fileName2Root = fileName2.substring(0, idx);

        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);

        CIEChromaticity cieC = new CIEChromaticity();
        int binnedImageMaxDimension = 512;
        int binFactor1 = (int) Math.ceil(
            Math.max((float) img1.getWidth() / (float)binnedImageMaxDimension,
            (float) img1.getHeight() / (float)binnedImageMaxDimension));
        int binFactor2 = (int) Math.ceil(
            Math.max((float) img2.getWidth() / (float)binnedImageMaxDimension,
            (float) img2.getHeight() / (float)binnedImageMaxDimension));
        ImageProcessor imageProcessor = new ImageProcessor();
        ImageExt imgBinned1 = imageProcessor.binImage(img1, binFactor1);
        ImageExt imgBinned2 = imageProcessor.binImage(img2, binFactor2);

        /*
        HistogramEqualizationForColor hEq = new HistogramEqualizationForColor(imgBinned1);
        hEq.applyFilter();
        hEq = new HistogramEqualizationForColor(imgBinned2);
        hEq.applyFilter();
        */
        RotatedOffsets rotatedOffsets = RotatedOffsets.getInstance();
        GreyscaleImage redBinnedImg1 = imgBinned1.copyRedToGreyscale();
        GreyscaleImage greenBinnedImg1 = imgBinned1.copyGreenToGreyscale();
        GreyscaleImage blueBinnedImg1 = imgBinned1.copyBlueToGreyscale();
        GreyscaleImage redBinnedImg2 = imgBinned2.copyRedToGreyscale();
        GreyscaleImage greenBinnedImg2 = imgBinned2.copyGreenToGreyscale();
        GreyscaleImage blueBinnedImg2 = imgBinned2.copyBlueToGreyscale();
        GreyscaleImage gsImg1 = imgBinned1.copyToGreyscale();
        GreyscaleImage gsImg2 = imgBinned2.copyToGreyscale();
        IntensityClrFeatures clrFeaturesBinned1 = new IntensityClrFeatures(gsImg1.copyImage(),
            5, rotatedOffsets);
        IntensityClrFeatures clrFeaturesBinned2 = new IntensityClrFeatures(gsImg2.copyImage(),
            5, rotatedOffsets);
        IntensityFeatures features1 = new IntensityFeatures(5, true, rotatedOffsets);
        IntensityFeatures features2 = new IntensityFeatures(5, true, rotatedOffsets);

        /*
        looking at trace/determinant of autocorrelation
        and eigenvalues of greyscale, autocorrelation, and lab colors for
        selected points in both images

        statues subsets:

        0 (64, 100) (96, 109)
        1 (67, 103) (103, 111)
        2 (68, 78)  (113, 86)
        3 (66, 49)  (106, 50)
        4 (92, 108) (157, 118)
        5 (92, 111) (160, 122)   delta e = 28.9
        6 (69, 129) (108, 142)   delta e = 26.4

        is the edelta for the gingerbread man's white stripes and shadows
        the same for shadow and higher illumination?
        */
        // single values of edelta
        List<PairInt> points1 = new ArrayList<PairInt>();
        List<PairInt> points2 = new ArrayList<PairInt>();
        points1.add(new PairInt(64, 100)); points2.add(new PairInt(96, 109));
        points1.add(new PairInt(67, 103)); points2.add(new PairInt(103, 111));
        points1.add(new PairInt(68, 78)); points2.add(new PairInt(113, 86));
        points1.add(new PairInt(66, 49)); points2.add(new PairInt(106, 50));
        points1.add(new PairInt(92, 108)); points2.add(new PairInt(157, 118));
        points1.add(new PairInt(92, 111)); points2.add(new PairInt(160, 122));
        points1.add(new PairInt(69, 129)); points2.add(new PairInt(108, 142));
        // this one is too look at localizability:
        points1.add(new PairInt(46, 65)); points2.add(new PairInt(6, 4));

        int n = points1.size();
        for (int i = 0; i < n; ++i) {

            StringBuilder sb = new StringBuilder();

            PairInt p1 = points1.get(i);
            PairInt p2 = points2.get(i);
            sb.append(p1.toString()).append(p2.toString());

            int r1 = redBinnedImg1.getValue(p1.getX(), p1.getY());
            int g1 = greenBinnedImg1.getValue(p1.getX(), p1.getY());
            int b1 = blueBinnedImg1.getValue(p1.getX(), p1.getY());
            int r2 = redBinnedImg2.getValue(p2.getX(), p2.getY());
            int g2 = greenBinnedImg2.getValue(p2.getX(), p2.getY());
            int b2 = blueBinnedImg2.getValue(p2.getX(), p2.getY());

            float[] lab1 = cieC.rgbToCIELAB(r1, g1, b1);
            float[] lab2 = cieC.rgbToCIELAB(r2, g2, b2);
            float[] cieXY1 = cieC.rgbToCIEXYZ(r1, g1, b1);
            float[] cieXY2 = cieC.rgbToCIEXYZ(r2, g2, b2);
            double deltaE = cieC.calcDeltaECIE94(lab1, lab2);

            sb.append(String.format("  dE=%.1f", (float)deltaE));

            int rot1 = clrFeaturesBinned1.calculateOrientation(p1.getX(), p1.getY());
            int rot2 = clrFeaturesBinned2.calculateOrientation(p2.getX(), p2.getY());

            IntensityDescriptor desc_l1 = clrFeaturesBinned1.extractIntensityLOfCIELAB(
                redBinnedImg1, greenBinnedImg1, blueBinnedImg1, p1.getX(), p1.getY(),
                rot1);

            IntensityDescriptor desc_a1 = clrFeaturesBinned1.extractIntensityAOfCIELAB(
                redBinnedImg1, greenBinnedImg1, blueBinnedImg1, p1.getX(), p1.getY(),
                rot1);

            IntensityDescriptor desc_b1 = clrFeaturesBinned1.extractIntensityBOfCIELAB(
                redBinnedImg1, greenBinnedImg1, blueBinnedImg1, p1.getX(), p1.getY(),
                rot1);

            IntensityDescriptor desc_l2 = clrFeaturesBinned2.extractIntensityLOfCIELAB(
                redBinnedImg2, greenBinnedImg2, blueBinnedImg2, p2.getX(), p2.getY(),
                rot2);

            IntensityDescriptor desc_a2 = clrFeaturesBinned2.extractIntensityAOfCIELAB(
                redBinnedImg2, greenBinnedImg2, blueBinnedImg2, p2.getX(), p2.getY(),
                rot2);

            IntensityDescriptor desc_b2 = clrFeaturesBinned2.extractIntensityBOfCIELAB(
                redBinnedImg2, greenBinnedImg2, blueBinnedImg2, p2.getX(), p2.getY(),
                rot2);

            IntensityDescriptor desc1 = features1.extractIntensity(gsImg1,
                p1.getX(), p1.getY(), rot1);

            IntensityDescriptor desc2 = features2.extractIntensity(gsImg2,
                p2.getX(), p2.getY(), rot2);

            double det, trace;
            SimpleMatrix a_l1 = clrFeaturesBinned1.createAutoCorrelationMatrix(desc_l1);
            det = a_l1.determinant();
            trace = a_l1.trace();
            sb.append(String.format("\n  L1_det(A)/trace=%.1f", (float)(det/trace)));
            SimpleMatrix a_l2 = clrFeaturesBinned2.createAutoCorrelationMatrix(desc_l2);
            det = a_l2.determinant();
            trace = a_l2.trace();
            sb.append(String.format("  L2_det(A)/trace=%.1f", (float)(det/trace)));

            try {
                sb.append("\n  eigen values:\n");
                SimpleEVD eigen1 = a_l1.eig();
                for (int j = 0; j < eigen1.getNumberOfEigenvalues(); ++j) {
                    Complex64F eigen = eigen1.getEigenvalue(j);
                    sb.append(String.format("    [1] %d %.1f %.1f\n", j,
                        (float)eigen.getReal(), (float)eigen.getMagnitude()));
                }
                sb.append("\n");
                SimpleEVD eigen2 = a_l2.eig();
                for (int j = 0; j < eigen2.getNumberOfEigenvalues(); ++j) {
                    Complex64F eigen = eigen2.getEigenvalue(j);
                    sb.append(String.format("    [1] %d %.1f %.1f\n", j,
                        (float)eigen.getReal(), (float)eigen.getMagnitude()));
                }
            } catch (Throwable t) {
            }

            if (desc1 != null && desc2 != null) {
                SimpleMatrix gs_l1 = features1.createAutoCorrelationMatrix(desc1);
                det = gs_l1.determinant();
                trace = gs_l1.trace();
                sb.append(String.format("\n  Grey_det(A)/trace=%.1f", (float)(det/trace)));
                SimpleMatrix gs_l2 = features2.createAutoCorrelationMatrix(desc2);
                det = gs_l2.determinant();
                trace = gs_l2.trace();
                sb.append(String.format("  Grey_det(A)/trace=%.1f", (float)(det/trace)));

                try {
                    sb.append("\n  eigen values:\n");
                    SimpleEVD eigen1 = gs_l1.eig();
                    for (int j = 0; j < eigen1.getNumberOfEigenvalues(); ++j) {
                        Complex64F eigen = eigen1.getEigenvalue(j);
                        sb.append(String.format("    [1] %d %.1f %.1f\n", j,
                            (float) eigen.getReal(), (float) eigen.getMagnitude()));
                    }
                    sb.append("\n");
                    SimpleEVD eigen2 = gs_l2.eig();
                    for (int j = 0; j < eigen2.getNumberOfEigenvalues(); ++j) {
                        Complex64F eigen = eigen2.getEigenvalue(j);
                        sb.append(String.format("    [2] %d %.1f %.1f\n", j,
                            (float) eigen.getReal(), (float) eigen.getMagnitude()));
                    }
                } catch (Throwable t) {
                }

            }

            log.info(sb.toString());
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

}
