package algorithms.imageProcessing.features;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.MedialAxis;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.compGeometry.RotatedOffsets;
import algorithms.compGeometry.clustering.KMeansHSV;
import algorithms.compGeometry.clustering.KMeansPlusPlus;
import algorithms.compGeometry.clustering.KMeansPlusPlusColor;
import algorithms.imageProcessing.AdaptiveThresholding;
import algorithms.imageProcessing.BetweenClassColorStats;
import algorithms.imageProcessing.BetweenClassColorStats.AllClassInterStats;
import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.CannyEdgeFilterAdaptive;
import algorithms.imageProcessing.CannyEdgeFilterAdaptiveDeltaE2000;
import algorithms.imageProcessing.ColorHistogram;
import algorithms.imageProcessing.DFSContiguousIntValueFinder;
import algorithms.imageProcessing.DFSContiguousValueFinder;
import algorithms.imageProcessing.EdgeFilterProducts;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.GroupPixelCIELAB;
import algorithms.imageProcessing.GroupPixelCIELAB1931;
import algorithms.imageProcessing.GroupPixelCIELCH;
import algorithms.imageProcessing.GroupPixelCIELUV;
import algorithms.imageProcessing.GroupPixelColors;
import algorithms.imageProcessing.GroupPixelHSV;
import algorithms.imageProcessing.GroupPixelRGB;
import algorithms.imageProcessing.GroupPixelRGB0;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageDisplayer;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageProcessor.Colors;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.ImageSegmentation.DecimatedData;
import algorithms.imageProcessing.IntraClassColorStats;
import algorithms.imageProcessing.PixelColors;
import algorithms.imageProcessing.matching.PartialShapeMatcher;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.SegmentationMergeThreshold;
import algorithms.imageProcessing.VanishingPoints;
import algorithms.imageProcessing.features.ORB.Descriptors;
import algorithms.imageProcessing.features.ObjectMatcher.Settings;
import algorithms.imageProcessing.matching.PartialShapeMatcher.Result;
import algorithms.imageProcessing.matching.SegmentedCellDescriptorMatcher;
import algorithms.imageProcessing.matching.ShapeFinder;
import algorithms.imageProcessing.segmentation.ColorSpace;
import algorithms.imageProcessing.segmentation.LabelToColorHelper;
import algorithms.imageProcessing.segmentation.NormalizedCuts;
import algorithms.imageProcessing.segmentation.SLICSuperPixels;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.GroupAverageColors;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.imageProcessing.util.RANSACAlgorithmIterations;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntPair;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.QuadInt;
import algorithms.util.ResourceFinder;
import algorithms.util.TwoDFloatArray;
import algorithms.util.TwoDIntArray;
import algorithms.util.VeryLongBitString;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.awt.image.ImageObserver;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
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

        int maxDimension = 256;//512;

        String fileName1 = "";

        //for (int i = 15; i < 16; ++i) {
        for (int i = 0; i < 37; ++i) {

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
                    break;
                }
                case 6: {
                    fileName1 = "cloudy_san_jose.jpg";
                    break;
                }
                case 7: {
                    fileName1 = "patagonia_snowy_foreground.jpg";
                    break;
                }
                case 8: {
                    fileName1 = "mt_rainier_snowy_field.jpg";
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
                    break;
                }
                case 12: {
                    fileName1 = "venturi_mountain_j6_0010.png";
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
                    break;
                }
                case 19: {
                    fileName1 = "norwegian_mtn_range.jpg";
                    break;
                }
                case 20: {
                    fileName1 = "halfdome.jpg";
                    break;
                }
                case 21: {
                    fileName1 = "halfdome2.jpg";
                    break;
                }
                case 22: {
                    fileName1 = "halfdome3.jpg";
                    break;
                }
                case 23: {
                    fileName1 = "costa_rica.jpg";
                    break;
                }
                case 24: {
                    fileName1 = "new-mexico-sunrise_w725_h490.jpg";
                    break;
                }
                case 25: {
                    fileName1 = "arizona-sunrise-1342919937GHz.jpg";
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
                    break;
                }
                case 31: {
                    fileName1 = "30.jpg";
                    break;
                }
                case 32: {
                    fileName1 = "arches_sun_01.jpg";
                    break;
                }
                case 33: {
                    fileName1 = "stlouis_arch.jpg";
                    break;
                }
                case 34: {
                    fileName1 = "contrail.jpg";
                    break;
                }
                case 35: {
                    fileName1 = "checkerboard_01.jpg";
                    break;
                }
                default: {
                    fileName1 = "checkerboard_02.jpg";
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

            int binFactor1 = (int) Math.ceil(Math.max(
                (float) w1 / maxDimension,
                (float) h1 / maxDimension));

            img = imageProcessor.binImage(img, binFactor1);

            CannyEdgeFilterAdaptiveDeltaE2000 canny =
                new CannyEdgeFilterAdaptiveDeltaE2000();
            //canny.overrideToUseAdaptiveThreshold();
            //canny.setToDebug();
            canny.applyFilter(img.copyToImageExt());
            EdgeFilterProducts products = canny.getFilterProducts();        
            
            int[] labels4 = imageSegmentation.objectSegmentation(img, products);

            ImageExt img11 = img.createWithDimensions();
            ImageIOHelper.addAlternatingColorLabelsToRegion(
                img11, labels4);
            MiscDebug.writeImage(img11, "_final_" + fileName1Root);
            //LabelToColorHelper.applyLabels(img, labels4);
            //MiscDebug.writeImage(img, "_final_" + fileName1Root);
             
            /*{// --- a look at the angles of phase and orientation plotted ----
                List<Set<PairInt>> contigSets = 
                    LabelToColorHelper.extractContiguousLabelPoints(
                    img, labels4);

                List<PairIntArray> orderedBoundaries = new ArrayList<PairIntArray>();

                int w = img.getWidth();
                int h = img.getHeight();
                SIGMA sigma = SIGMA.ZEROPOINTFIVE;//SIGMA.ONE;

                for (int ii = 0; ii < contigSets.size(); ++ii) {
                    Set<PairInt> set = contigSets.get(ii);
                    Set<PairInt> medialAxis = new HashSet<PairInt>();
                    PairIntArray p = imageProcessor.extractSmoothedOrderedBoundary(
                         set, sigma, w, h, medialAxis);
                    if (p.getN() > 24) {
                        orderedBoundaries.add(p);
                    }
                }

                EdgeFilterProducts products = imageSegmentation
                    .createPhaseCongruencyGradient(
                    img.copyToGreyscale());

                img11 = img.copyToImageExt();

                for (int ii = 0; ii < orderedBoundaries.size(); ++ii) {
                    PairIntArray a = orderedBoundaries.get(ii);
                    for (int j = 0; j < a.getN(); j += 10) {
                        int x = a.getX(j);
                        int y = a.getY(j);
                        double or = products.getTheta().getValue(x, y)
                            * Math.PI/180.;
                        double pa = products.getPhaseAngle().getValue(x, y)
                            * Math.PI/180.;
                        int dx0 = (int)Math.round(3. * Math.cos(or));
                        int dy0 = (int)Math.round(3. * Math.sin(or));
                        int dx1 = (int)Math.round(3. * Math.cos(pa));
                        int dy1 = (int)Math.round(3. * Math.sin(pa));

                        ImageIOHelper.addPointToImage(x, y, img11, 1, 255, 0, 0);
                        
                        int x2, y2;
                        
                        x2 = x + dx0;
                        y2 = y + dy0;
                        if (x2 >= 0 && x2 < w && y2 >= 0 && y2 < h) {
                            ImageIOHelper.drawLineInImage(x, y, 
                                x2, y2, img11, 0, 255, 255, 0);
                        }
                        x2 = x + dx1;
                        y2 = y + dy1;
                        if (x2 >= 0 && x2 < w && y2 >= 0 && y2 < h) {
                            ImageIOHelper.drawLineInImage(x, y, 
                                x2, y2, img11, 0, 0, 0, 255);
                        }
                    }
                }
                MiscDebug.writeImage(img11, "_aa_" + fileName1Root);
            }*/
        }
    }

    public void testORBMatcher_gingerbreadman() throws Exception {

        /*        
        this demonstrates ORB
            followed by filtering of search image keypoints by color.
            then matching by descriptors 
              and evaluation of pair combinations of best mathing keypoints
              from which euclidean transformaions are derived.
        */
    
        int maxDimension = 256;//512;
        SIGMA sigma = SIGMA.ZEROPOINTFIVE;//SIGMA.ONE;

        ImageProcessor imageProcessor = new ImageProcessor();
        ImageSegmentation imageSegmentation = new ImageSegmentation();

        String[] fileNames0 = new String[]{
            "android_statues_03_sz1",
        //    "android_statues_03_sz3"
        };

        /*
        NOTE: andr 03 and 04 have a very small number of keypoints, but
        just enough to match along w/ the very few number of
        sergmented cells remaining after filter which could possibly
        be matched.
        */
        String[] fileNames1 = new String[]{
        //    "android_statues_01.jpg", // no kp remaining
              "android_statues_02.jpg", // no kp remaining
        //    "android_statues_04.jpg", // OK: match0, matchSmall works
        //    "android_statues_03.jpg" // OK: match0, matchSmall works
        };
        
        int fn0 = 0;
        for (String fileNameRoot0 : fileNames0) {               
            fn0++;
            for (String fileName1 : fileNames1) {               
        
                long t0 = System.currentTimeMillis();
                
                Set<PairInt> shape0 = new HashSet<PairInt>();

                // to compare to "android_statues_01.jpg",
                //    set this to '2'
                int binFactor0 = 1;
               
                // 1st image is color image, 2nd is masked color image
                // 218 X 163... full size is 1280 X 960
                ImageExt[] imgs0 = maskAndBin(fileNameRoot0, 
                    binFactor0, shape0);
                int nShape0_0 = shape0.size();
               
                System.out.println("shape0 nPts=" + nShape0_0);
                
                
                String fileName1Root = fileName1.substring(0, 
                    fileName1.lastIndexOf("."));
                String filePath1 = ResourceFinder.findFileInTestResources(
                    fileName1);
                ImageExt img = ImageIOHelper.readImageExt(filePath1);
        
                // template img size is 218  163
                //img = (ImageExt) imageProcessor.bilinearDownSampling(
                //    img, 218, 163, 0, 255);
                
                long ts = MiscDebug.getCurrentTimeFormatted();
                
                int w1 = img.getWidth();
                int h1 = img.getHeight();
        
                int binFactor1 = (int) Math.ceil(Math.max(
                    (float) w1 / maxDimension,
                    (float) h1 / maxDimension));
        
                img = imageProcessor.binImage(img, binFactor1);
               
                int w = img.getWidth();
                int h = img.getHeight();

                /*RANSACAlgorithmIterations nsIter = new
                    RANSACAlgorithmIterations();
                long nnn = nsIter.estimateNIterFor99PercentConfidence(
                    300, 7, 20./300.);
                System.out.println("99 percent nIter for RANSAC=" 
                    + nnn);*/
                 
                Settings settings = new Settings();
                settings.setToFindVnishingPoints();
                           
                ObjectMatcher objMatcher = new ObjectMatcher();
                //if (fileName1Root.contains("_01")) {
                    settings.setToUseLargerPyramid0();
                    settings.setToUseLargerPyramid1();
                    settings.setToUseSmallObjectMethod();
                //}
                objMatcher.setToDebug();
                CorrespondenceList cor = objMatcher.findObject(imgs0[0], shape0, 
                    img, settings);
                
                long t1 = System.currentTimeMillis();
                System.out.println("matching took " + ((t1 - t0)/1000.) + " sec");
                         
                if (cor == null) {
                    continue;
                }
                
                CorrespondencePlotter plotter = new CorrespondencePlotter(
                    imgs0[1], img.copyImage());            
                for (int ii = 0; ii < cor.getPoints1().size(); ++ii) {
                    PairInt p1 = cor.getPoints1().get(ii);
                    PairInt p2 = cor.getPoints2().get(ii);

                    //System.out.println("orb matched: " + p1 + " " + p2);
                    //if (p2.getX() > 160)
                    plotter.drawLineInAlternatingColors(p1.getX(), p1.getY(), 
                        p2.getX(), p2.getY(), 0);
                }

                plotter.writeImage("_orb_corres_final_" + 
                    "_" + fileName1Root + "_" + fn0);
                System.out.println(cor.getPoints1().size() + 
                    " matches " + fileName1Root);
                //MiscDebug.writeImage(img11, "_orb_matched_" + str
                //    + "_" + fileName1Root);
            }
        }
    }

    public void estORBMatcher_cupcake() throws Exception {
        
        // TODO: this one either needs more keypoints across the cupcake in
        //       android statues 02 image,
        //    or it needs an algorithm similar to matchSmall which does not
        //    use descriptors, but instead uses color histograms and partial
        //    shape matching, but with the addition of aggregated shape
        //    comparisons (== the unfinished ShapeFinder)
        //
        //   and as always, improved segmentation would help, but the cupcake is
        //   in partly shaded locations.
        
        int maxDimension = 256;//512;
        SIGMA sigma = SIGMA.ZEROPOINTFIVE;//SIGMA.ONE;

        ImageProcessor imageProcessor = new ImageProcessor();
        ImageSegmentation imageSegmentation = new ImageSegmentation();

        String[] fileNames0 = new String[]{
            "android_statues_04.jpg",
            "android_statues_04_cupcake_mask.png",
        };

        String[] fileNames1 = new String[]{
             "android_statues_01.jpg",   // needs aggregated shape matching
             "android_statues_02.jpg", // needs aggregated shape matching 
             "android_statues_04.jpg", // descr are fine
        };

        for (String fileName1 : fileNames1) {               

            long t0 = System.currentTimeMillis();

            Set<PairInt> shape0 = new HashSet<PairInt>();

            ImageExt[] imgs0 = maskAndBin2(fileNames0, 
                maxDimension, shape0);
            int nShape0_0 = shape0.size();

            System.out.println("shape0 nPts=" + nShape0_0);


            String fileName1Root = fileName1.substring(0, 
                fileName1.lastIndexOf("."));
            String filePath1 = ResourceFinder.findFileInTestResources(
                fileName1);
            ImageExt img = ImageIOHelper.readImageExt(filePath1);

            long ts = MiscDebug.getCurrentTimeFormatted();

            int w1 = img.getWidth();
            int h1 = img.getHeight();

            int binFactor1 = (int) Math.ceil(Math.max(
                (float) w1 / maxDimension,
                (float) h1 / maxDimension));

            img = imageProcessor.binImage(img, binFactor1);

            int w = img.getWidth();
            int h = img.getHeight();

            /*RANSACAlgorithmIterations nsIter = new
                RANSACAlgorithmIterations();
            long nnn = nsIter.estimateNIterFor99PercentConfidence(
                300, 7, 20./300.);
            System.out.println("99 percent nIter for RANSAC=" 
                + nnn);*/

            //GreyscaleImage theta1 = imageProcessor.createCIELABTheta(imgs0[0], 255);
            //MiscDebug.writeImage(theta1, fileName1Root + "_theta_0");
            //theta1 = imageProcessor.createCIELABTheta(img, 255);
            //MiscDebug.writeImage(theta1, fileName1Root + "_theta_1");
        
            Settings settings = new Settings();
            settings.setToFindVnishingPoints();
            
            ObjectMatcher objMatcher = new ObjectMatcher();
          //  if (
          //      fileName1Root.contains("_01") ||
          //      fileName1Root.contains("_02")) {
          //      settings.setToUseShapeFinderMethod();
          //  }
            
            //settings.setToUseLargerPyramid0();
            objMatcher.setToDebug();
            CorrespondenceList cor = objMatcher.findObject(imgs0[0], shape0, 
                img, settings);

            long t1 = System.currentTimeMillis();
            System.out.println("matching took " + ((t1 - t0)/1000.) + " sec");

            if (cor == null) {
                continue;
            }
            
            CorrespondencePlotter plotter = new CorrespondencePlotter(
                imgs0[1], img.copyImage());            
            for (int ii = 0; ii < cor.getPoints1().size(); ++ii) {
                PairInt p1 = cor.getPoints1().get(ii);
                PairInt p2 = cor.getPoints2().get(ii);

                //System.out.println("orb matched: " + p1 + " " + p2);
                //if (p2.getX() > 160)
                plotter.drawLineInAlternatingColors(p1.getX(), p1.getY(), 
                    p2.getX(), p2.getY(), 0);
            }

            plotter.writeImage("_orb_corres_final_" + 
                "_" + fileName1Root);
            System.out.println(cor.getPoints1().size() + 
                " matches " + fileName1Root + " test2");
            //MiscDebug.writeImage(img11, "_orb_matched_" + str
            //    + "_" + fileName1Root);
        }
    }

    public void estORBMatcher_icecream() throws Exception {

        // TODO: needs better segmentation for the icecream in status 01 and 02
        //    AND/OR a different light source for polar theta CIE LAB
        
        int maxDimension = 256;//512;
        SIGMA sigma = SIGMA.ZEROPOINTFIVE;//SIGMA.ONE;

        ImageProcessor imageProcessor = new ImageProcessor();
        ImageSegmentation imageSegmentation = new ImageSegmentation();

        /*
        TODO: to match this:
           -- may need to use a later cielab for the polar theta image 
              for color space with corrections for neutral hues..
              colors close to white
           -- may need to define limits based upon segmentation,
              and/or return top results.
              the shapefinder requires good segmentation so need to look at
              the super pixels to see that the cupcake boundary is actually
              findable in the first place.
        */
        String[] fileNames0 = new String[]{
            "android_statues_04.jpg",
            "android_statues_04_icecream_mask.png",
        };

        String[] fileNames1 = new String[]{
            "android_statues_01.jpg",  
            "android_statues_02.jpg",
            "android_statues_04.jpg", // descr are fine
        };

        for (String fileName1 : fileNames1) {               

            long t0 = System.currentTimeMillis();

            Set<PairInt> shape0 = new HashSet<PairInt>();

            ImageExt[] imgs0 = maskAndBin2(fileNames0, 
                maxDimension, shape0);
            int nShape0_0 = shape0.size();

            System.out.println("shape0 nPts=" + nShape0_0);


            String fileName1Root = fileName1.substring(0, 
                fileName1.lastIndexOf("."));
            String filePath1 = ResourceFinder.findFileInTestResources(
                fileName1);
            ImageExt img = ImageIOHelper.readImageExt(filePath1);

            long ts = MiscDebug.getCurrentTimeFormatted();

            int w1 = img.getWidth();
            int h1 = img.getHeight();

            int binFactor1 = (int) Math.ceil(Math.max(
                (float) w1 / maxDimension,
                (float) h1 / maxDimension));

            img = imageProcessor.binImage(img, binFactor1);

            int w = img.getWidth();
            int h = img.getHeight();

            /*RANSACAlgorithmIterations nsIter = new
                RANSACAlgorithmIterations();
            long nnn = nsIter.estimateNIterFor99PercentConfidence(
                300, 7, 20./300.);
            System.out.println("99 percent nIter for RANSAC=" 
                + nnn);*/

            GreyscaleImage theta1 = imageProcessor.createCIELUVTheta(imgs0[0], 255);
            MiscDebug.writeImage(theta1, fileName1Root + "_theta_0");
            theta1 = imageProcessor.createCIELUVTheta(img, 255);
            MiscDebug.writeImage(theta1, fileName1Root + "_theta_1");
        
            Settings settings = new Settings();
            settings.setToFindVnishingPoints();
            
            ObjectMatcher objMatcher = new ObjectMatcher();
            //if (
            //    fileName1Root.contains("_01") ||
            //    fileName1Root.contains("_02")) {
            //    settings.setToUseShapeFinderMethod();
            //}
            //settings.setToUseLargerPyramid0();
            objMatcher.setToDebug();
            CorrespondenceList cor = objMatcher.findObject(imgs0[0], shape0, 
                img, settings);

            long t1 = System.currentTimeMillis();
            System.out.println("matching took " + ((t1 - t0)/1000.) + " sec");

            if (cor == null) {
                continue;
            }
            
            CorrespondencePlotter plotter = new CorrespondencePlotter(
                imgs0[1], img.copyImage());            
            for (int ii = 0; ii < cor.getPoints1().size(); ++ii) {
                PairInt p1 = cor.getPoints1().get(ii);
                PairInt p2 = cor.getPoints2().get(ii);

                //System.out.println("orb matched: " + p1 + " " + p2);
                //if (p2.getX() > 160)
                plotter.drawLineInAlternatingColors(p1.getX(), p1.getY(), 
                    p2.getX(), p2.getY(), 0);
            }

            plotter.writeImage("_orb_corres_final_" + 
                "_" + fileName1Root);
            System.out.println(cor.getPoints1().size() + 
                " matches " + fileName1Root + " test3");
            //MiscDebug.writeImage(img11, "_orb_matched_" + str
            //    + "_" + fileName1Root);
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

    private void printGradients(ImageExt img,
        String fileNameRoot) {

        GreyscaleImage gsImg = img.copyToGreyscale();
        GreyscaleImage gsImg1 = gsImg.copyImage();
        GreyscaleImage gsImg2 = gsImg.copyImage();

        CannyEdgeFilterAdaptive canny =
            new CannyEdgeFilterAdaptive();
        canny.overrideToNotUseLineThinner();
        canny.applyFilter(gsImg);
        for (int i = 0; i < gsImg.getNPixels(); ++i) {
            if (gsImg.getValue(i) > 0) {
                gsImg.setValue(i, 255);
            }
        }
        MiscDebug.writeImage(gsImg,
            "_canny_" + fileNameRoot);

        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.applyFirstDerivGaussian(gsImg1,
            SIGMA.ONE, 0, 255);
        for (int i = 0; i < gsImg1.getNPixels(); ++i) {
            if (gsImg1.getValue(i) > 0) {
                gsImg1.setValue(i, 255);
            }
        }
        MiscDebug.writeImage(gsImg1,
            "_firstderiv_" + fileNameRoot);

        PhaseCongruencyDetector pcd = new PhaseCongruencyDetector();
        PhaseCongruencyDetector.PhaseCongruencyProducts
            product = pcd.phaseCongMono(gsImg2);
        MiscDebug.writeImage(product.getThinned(),
            "_pcd_" + fileNameRoot);
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

    private ImageExt[] maskAndBin(String fileNamePrefix, int binFactor,
        Set<PairInt> outputShape) throws IOException, Exception {

        ImageProcessor imageProcessor = new ImageProcessor();

        String fileNameMask0 = fileNamePrefix + "_mask.png";
        String filePathMask0 = ResourceFinder
            .findFileInTestResources(fileNameMask0);
        ImageExt imgMask0 = ImageIOHelper.readImageExt(filePathMask0);

        String fileName0 = fileNamePrefix + ".jpg";
        String filePath0 = ResourceFinder
            .findFileInTestResources(fileName0);
        ImageExt img0 = ImageIOHelper.readImageExt(filePath0);
    
        if (binFactor != 1) {
            img0 = imageProcessor.binImage(img0, binFactor);
            imgMask0 = imageProcessor.binImage(imgMask0, binFactor);
        }
        
        ImageExt img0Masked = img0.copyToImageExt();
        
        assertEquals(imgMask0.getNPixels(), img0.getNPixels());

        for (int i = 0; i < imgMask0.getNPixels(); ++i) {
            if (imgMask0.getR(i) == 0) {
                img0Masked.setRGB(i, 0, 0, 0);
            } else {
                outputShape.add(new PairInt(imgMask0.getCol(i), imgMask0.getRow(i)));
            }
        }
        //MiscDebug.writeImage(img0Masked, "_MASKED");
   
        return new ImageExt[]{img0, img0Masked};
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
            System.err.println(e.getMessage());
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

    private void extractTemplateKeypoints(String fileNameRoot0,
        Set<PairInt> shape0, PairIntArray template,
        List<PairInt> templateKP, TDoubleList templateOrientations,
        Descriptors templateDescriptors) throws IOException, Exception {

        String fileName0 = fileNameRoot0 + ".jpg";
        String filePath0 = ResourceFinder
            .findFileInTestResources(fileName0);
        ImageExt img0 = ImageIOHelper.readImageExt(filePath0);

        int[] minMaxXY = MiscMath.findMinMaxXY(template);

        int w = img0.getWidth();
        int h = img0.getHeight();

        int xLL = minMaxXY[0] - 5;
        if (xLL < 0) {
            xLL = 0;
        }
        int yLL = minMaxXY[2] - 5;
        if (yLL < 0) {
            yLL = 0;
        }
        int xUR = minMaxXY[1] + 5;
        if (xUR > (w - 1)) {
            xUR = w - 1;
        }
        int yUR = minMaxXY[3] + 5;
        if (yUR > (h - 1)) {
            yUR = h - 1;
        }

        ORB.DescriptorDithers descrOffsets 
            = ORB.DescriptorDithers.NONE;
        //    = ORB.DescriptorDithers.FORTY_FIVE;
        //    = ORB.DescriptorDithers.FIFTEEN;
        
        ORBWrapper.extractKeypointsFromSubImage(
            img0, xLL, yLL, xUR, yUR,
            200, templateKP, templateOrientations, 
            templateDescriptors, 
            //0.01f,
            0.001f,
            true,
            descrOffsets);
        
        for (int i = 0; i < templateKP.size(); ++i) {
            PairInt p = templateKP.get(i);
            if (shape0.contains(p)) {
                ImageIOHelper.addPointToImage(p.getX(), p.getY(), img0, 1, 255, 0, 0);
            }
        }

        MiscDebug.writeImage(img0, "_template_orb");
    }

    private void extractTemplateORBKeypoints(ImageExt img,
        Set<PairInt> shape0, 
        List<PairInt> templateKP, TDoubleList templateOrientations,
        Descriptors templateDescriptorsH, Descriptors templateDescriptorsS,
        Descriptors templateDescriptorsV) throws IOException, Exception {

        int[] minMaxXY = MiscMath.findMinMaxXY(shape0);

        int w = img.getWidth();
        int h = img.getHeight();
        
        int buffer = 20;

        int xLL = minMaxXY[0] - buffer;
        if (xLL < 0) {
            xLL = 0;
        }
        int yLL = minMaxXY[2] - buffer;
        if (yLL < 0) {
            yLL = 0;
        }
        int xUR = minMaxXY[1] + buffer;
        if (xUR > (w - 1)) {
            xUR = w - 1;
        }
        int yUR = minMaxXY[3] + buffer;
        if (yUR > (h - 1)) {
            yUR = h - 1;
        }

        ORBWrapper.extractKeypointsFromSubImage(
            img, xLL, yLL, xUR, yUR,
            200,
            //100,
            templateKP, templateOrientations, 
            templateDescriptorsH, 
            templateDescriptorsS,
            templateDescriptorsV, 
            //0.01f, 
            0.001f,
            true);
        
        ImageExt imgCp = img.copyToImageExt();
        
        TIntList rm = new TIntArrayList();
        for (int i = 0; i < templateKP.size(); ++i) {
            PairInt p = templateKP.get(i);
            if (shape0.contains(p)) {
                ImageIOHelper.addPointToImage(p.getX(), p.getY(), 
                    imgCp, 1, 255, 0, 0);
            } else {
                rm.add(i);
                System.out.println("removing " + p);
            }
        }
        
        if (!rm.isEmpty()) {
            for (int i = (rm.size() - 1); i > -1; --i) {
                
                int rmIdx = rm.get(i);
                
                templateKP.remove(rmIdx);
                templateOrientations.removeAt(rmIdx);
                
                // move up operations.  everything with index > i moves up by 1
                for (int j = (i + 1); j < templateDescriptorsH.descriptors.length;
                    ++j) {
                    templateDescriptorsH.descriptors[j - 1] =
                        templateDescriptorsH.descriptors[j];
                    templateDescriptorsS.descriptors[j - 1] =
                        templateDescriptorsS.descriptors[j];
                    templateDescriptorsV.descriptors[j - 1] =
                        templateDescriptorsV.descriptors[j];
                }
            }
            
            int count = templateDescriptorsH.descriptors.length - rm.size();
            
            templateDescriptorsH.descriptors = 
                Arrays.copyOf(templateDescriptorsH.descriptors, count);
            
            templateDescriptorsS.descriptors = 
                Arrays.copyOf(templateDescriptorsS.descriptors, count);
            
            templateDescriptorsV.descriptors = 
                Arrays.copyOf(templateDescriptorsV.descriptors, count);
        }

        MiscDebug.writeImage(imgCp, "_template_orb");
    }

    private ORB extractTemplateORBKeypoints(ImageExt img,
        Set<PairInt> shape0, 
        int nKeypoints, float fastThresh,
        boolean useSmallPyramid,
        boolean createFirstDerivPts,
        boolean createCurvaturePts) throws IOException, Exception {

        int[] minMaxXY = MiscMath.findMinMaxXY(shape0);

        int w = img.getWidth();
        int h = img.getHeight();
        
        int buffer = 20;

        int xLL = minMaxXY[0] - buffer;
        if (xLL < 0) {
            xLL = 0;
        }
        int yLL = minMaxXY[2] - buffer;
        if (yLL < 0) {
            yLL = 0;
        }
        int xUR = minMaxXY[1] + buffer;
        if (xUR > (w - 1)) {
            xUR = w - 1;
        }
        int yUR = minMaxXY[3] + buffer;
        if (yUR > (h - 1)) {
            yUR = h - 1;
        }
        
        boolean overrideToCreateSmallestPyramid = useSmallPyramid;
        
        ORB orb = ORBWrapper.extractHSVKeypointsFromSubImage(
            img, xLL, yLL, xUR, yUR,
            nKeypoints,
            fastThresh, 
            createFirstDerivPts, createCurvaturePts,
            overrideToCreateSmallestPyramid);
                        
        // trim orb data that is outside of shape
        int ns = orb.getKeyPoint0List().size();
        
        for (int i = 0; i < ns; ++i) {
            TIntList kp0 = orb.getKeyPoint0List().get(i);
            TIntList kp1 = orb.getKeyPoint1List().get(i);
            TDoubleList or = orb.getOrientationsList().get(i);
            TFloatList s = orb.getScalesList().get(i);
            Descriptors dH = orb.getDescriptorsH().get(i);
            Descriptors dS = orb.getDescriptorsS().get(i);
            Descriptors dV = orb.getDescriptorsV().get(i);
            
            int n0 = kp0.size();
            
            TIntList rm = new TIntArrayList();
            for (int j = 0; j < n0; ++j) {
                PairInt p = new PairInt(kp1.get(j), kp0.get(j));
                if (!shape0.contains(p)) {
                    rm.add(j);
                }
            }
            if (!rm.isEmpty()) {
                int nb = n0 - rm.size();
                Descriptors dH2 = new Descriptors();
                dH2.descriptors = new VeryLongBitString[nb];
                Descriptors dS2 = new Descriptors();
                dS2.descriptors = new VeryLongBitString[nb];
                Descriptors dV2 = new Descriptors();
                dV2.descriptors = new VeryLongBitString[nb];
                
                TIntSet rmSet = new TIntHashSet(rm);
                for (int j = (rm.size() - 1); j > -1; --j) {
                    int idx = rm.get(j);
                    kp0.removeAt(idx);
                    kp1.removeAt(idx);
                    or.removeAt(idx);
                    s.removeAt(idx);
                }
                int count = 0;
                for (int j = 0; j < n0; ++j) {
                    if (rmSet.contains(j)) {
                        continue;
                    }
                    dH2.descriptors[count] = dH.descriptors[j];
                    dS2.descriptors[count] = dS.descriptors[j];
                    dV2.descriptors[count] = dV.descriptors[j];
                    count++;
                }
                assert(count == nb);
                dH.descriptors = dH2.descriptors;
                dS.descriptors = dS2.descriptors;
                dV.descriptors = dV2.descriptors;
            }
        }
        
        {// DEBUG print each pyramid to see if has matchable points
            // might need to change the ORb response filter to scale by scale level
            for (int i0 = 0; i0 < orb.getKeyPoint0List().size(); ++i0) {
                Image img0Cp = img.copyImage();
                float scale = orb.getScalesList().get(i0).get(0);
                for (int i = 0; i < orb.getKeyPoint0List().get(i0).size(); ++i) {
                    int y = orb.getKeyPoint0List().get(i0).get(i);
                    int x = orb.getKeyPoint1List().get(i0).get(i);
                    ImageIOHelper.addPointToImage(x, y, img0Cp, 
                        1, 255, 0, 0);
                }
                String str = Integer.toString(i0);
                if (str.length() < 2) {
                    str = "0" + str;
                }
                MiscDebug.writeImage(img0Cp, "_template_orb" + str);
            }
        }
        
        return orb;
    }
  
    private TDoubleList extractKeypoints(ImageExt img, 
        List<Set<PairInt>> listOfPointSets,
        List<PairInt> keypoints, Descriptors descriptors) throws IOException, Exception {

        // bins of size template size across image

        int w = img.getWidth();
        int h = img.getHeight();

        ORB orb = new ORB(1000);
        //orb.overrideFastThreshold(0.01f);
        orb.overrideFastThreshold(0.001f);
        orb.overrideToCreateHSVDescriptors();
        orb.overrideToAlsoCreate1stDerivKeypoints();
        orb.detectAndExtract(img);

        List<PairInt> kp = orb.getAllKeyPoints();
        
        Descriptors d = orb.getAllDescriptors();
        
        TDoubleList or = orb.getAllOrientations();

        ImageExt img0 = img.copyToImageExt();

        Set<PairInt> points = new HashSet<PairInt>();
        for (Set<PairInt> set : listOfPointSets) {
            points.addAll(set);
        }

        Set<PairInt> exists = new HashSet<PairInt>();
        TDoubleList orientations = new TDoubleArrayList();
        for (int i = 0; i < kp.size(); ++i) {
            PairInt p = kp.get(i);
            if (exists.contains(p) || !points.contains(p)) {
                continue;
            }
            exists.add(p);
            keypoints.add(p);
            orientations.add(or.get(i));
            ImageIOHelper.addPointToImage(p.getX(), p.getY(), img0, 1, 255, 0, 0);
        }
        
        exists.clear();
        
        VeryLongBitString[] outD = new VeryLongBitString[keypoints.size()];
        int count = 0;
        for (int i = 0; i < kp.size(); ++i) {
            PairInt p = kp.get(i);
            if (exists.contains(p) || !points.contains(p)) {
                continue;
            }
            exists.add(p);
            outD[count] = d.descriptors[i];
            count++;
        }
        descriptors.descriptors = outD;
        
        MiscDebug.writeImage(img0, "_srch_orb");

        return orientations;
    }

    private TDoubleList extractORBKeypoints(ImageExt img, 
        List<Set<PairInt>> listOfPointSets,
        List<PairInt> keypoints, Descriptors descriptorsH,
        Descriptors descriptorsS, Descriptors descriptorsV) 
        throws IOException, Exception {

        // bins of size template size across image

        int w = img.getWidth();
        int h = img.getHeight();

        ORB orb = new ORB(2000);//10000
        //orb.overrideFastThreshold(0.01f);
        orb.overrideFastThreshold(0.001f);
        orb.overrideToCreateHSVDescriptors();
        orb.overrideToAlsoCreate1stDerivKeypoints();
        orb.overrideToCreateCurvaturePoints();
        //orb.overrideToCreateOffsetsToDescriptors(ORB.DescriptorDithers.FIFTEEN);
        orb.detectAndExtract(img);

        List<PairInt> kp = orb.getAllKeyPoints();
        
        Descriptors[] dHSV = orb.getAllDescriptorsHSV();
        
        TDoubleList or = orb.getAllOrientations();

        ImageExt img0 = img.copyToImageExt();

        Set<PairInt> points = new HashSet<PairInt>();
        for (Set<PairInt> set : listOfPointSets) {
            points.addAll(set);
        }

        Set<PairInt> exists = new HashSet<PairInt>();
        TDoubleList orientations = new TDoubleArrayList();
        for (int i = 0; i < kp.size(); ++i) {
            PairInt p = kp.get(i);
            if (exists.contains(p) || !points.contains(p)) {
                continue;
            }
            exists.add(p);
            keypoints.add(p);
            orientations.add(or.get(i));
            ImageIOHelper.addPointToImage(p.getX(), p.getY(), img0, 1, 255, 0, 0);
        }
        
        exists.clear();
        
        VeryLongBitString[] outH = new VeryLongBitString[keypoints.size()];
        VeryLongBitString[] outS = new VeryLongBitString[keypoints.size()];
        VeryLongBitString[] outV = new VeryLongBitString[keypoints.size()];
        int count = 0;
        for (int i = 0; i < kp.size(); ++i) {
            PairInt p = kp.get(i);
            if (exists.contains(p) || !points.contains(p)) {
                continue;
            }
            exists.add(p);
            outH[count] = dHSV[0].descriptors[i];
            outS[count] = dHSV[1].descriptors[i];
            outV[count] = dHSV[2].descriptors[i];
            count++;
        }
        descriptorsH.descriptors = outH;
        descriptorsS.descriptors = outS;
        descriptorsV.descriptors = outV;
        
        MiscDebug.writeImage(img0, "_srch_orb");

        return orientations;
    }

    private Set<PairInt> createMedialAxis(Set<PairInt> points, 
        Set<PairInt> perimeter) {

        MedialAxis medAxis = new MedialAxis(points, perimeter);
        medAxis.fastFindMedialAxis();

        Set<PairInt> medAxisPts = medAxis.getMedialAxisPoints();

        return medAxisPts;
    }

    public void estCIETheta() throws Exception {
    
        int maxDimension = 256;//512;
        SIGMA sigma = SIGMA.ZEROPOINTFIVE;//SIGMA.ONE;

        ImageProcessor imageProcessor = new ImageProcessor();
        
        String[] fileNames1 = new String[]{
             "android_statues_01.jpg",
             "android_statues_02.jpg",
             "android_statues_04.jpg",
             "android_statues_03.jpg"
        };
        for (String fileName1 : fileNames1) {               
       
        String fileName1Root = fileName1.substring(0, fileName1.lastIndexOf("."));
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img = ImageIOHelper.readImageExt(filePath1);

        // template img size is 218  163
        //img = (ImageExt) imageProcessor.bilinearDownSampling(
        //    img, 218, 163, 0, 255);
        
        long ts = MiscDebug.getCurrentTimeFormatted();
        
        int w1 = img.getWidth();
        int h1 = img.getHeight();

        int binFactor1 = (int) Math.ceil(Math.max(
            (float) w1 / maxDimension,
            (float) h1 / maxDimension));

        img = imageProcessor.binImage(img, binFactor1);
       
        int w = img.getWidth();
        int h = img.getHeight();
        
        long t0 = System.currentTimeMillis();
    
        // descriptors w/ masks
        /*corList = ORB.match2(
            orb0, orb, tempListOfPointSets, listOfPointSets2,
            1.5f, 0.1f, false);
        */
       
        /*
        GreyscaleImage theta1 = imageProcessor.createCIELABTheta(img, 255);
        MiscDebug.writeImage(theta1, fileName1Root + "_theta_");
        GreyscaleImage theta_15 = 
            imageProcessor.createCIELABTheta(img, 255, 15);
        MiscDebug.writeImage(theta_15, fileName1Root + "_theta_15_");
        */
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        ImageExt imgCp = img.copyToImageExt();
        int[] labels = imageSegmentation.objectSegmentation3(imgCp);
        ImageIOHelper.addAlternatingColorLabelsToRegion(
            //LabelToColorHelper.applyLabels(
            imgCp, labels);
        MiscDebug.writeImage(imgCp, "_theta_segmentation_" + fileName1Root);

        imgCp = img.copyToImageExt();
        labels = imageSegmentation.objectSegmentation(imgCp);
        ImageIOHelper.addAlternatingColorLabelsToRegion(
            //LabelToColorHelper.applyLabels(
            imgCp, labels);
        MiscDebug.writeImage(imgCp, "_hsv_segmentation_" + fileName1Root);
        
        }
    }

    private ImageExt[] maskAndBin2(String[] fileNames, 
        int maxDimension, Set<PairInt> outputShape) throws IOException {
        
        ImageProcessor imageProcessor = new ImageProcessor();

        String fileNameMask0 = fileNames[1];
        String filePathMask0 = ResourceFinder
            .findFileInTestResources(fileNameMask0);
        ImageExt imgMask0 = ImageIOHelper.readImageExt(filePathMask0);

        String fileName0 = fileNames[0];
        String filePath0 = ResourceFinder
            .findFileInTestResources(fileName0);
        ImageExt img0 = ImageIOHelper.readImageExt(filePath0);
    
        int w0 = img0.getWidth();
        int h0 = img0.getHeight();

        int binFactor0 = (int) Math.ceil(Math.max(
             (float) w0 / maxDimension,
             (float) h0 / maxDimension));
        
        img0 = imageProcessor.binImage(img0, binFactor0);
        imgMask0 = imageProcessor.binImage(imgMask0, binFactor0);
        
        ImageExt img0Masked = img0.copyToImageExt();
        
        assertEquals(imgMask0.getNPixels(), img0.getNPixels());

        for (int i = 0; i < imgMask0.getNPixels(); ++i) {
            if (imgMask0.getR(i) == 0) {
                img0Masked.setRGB(i, 0, 0, 0);
            } else {
                outputShape.add(new PairInt(imgMask0.getCol(i), 
                    imgMask0.getRow(i)));
            }
        }
   
        return new ImageExt[]{img0, img0Masked};
    }
    
    public void estStatueColorMethods() throws Exception {
        
        ImageExt[] icecream = loadMaskedIcecream();
        ImageExt[] cupcake = loadMaskedCupcake();
        ImageExt[] gingerBreadman = loadMaskedGingerBreadMan();
        ImageExt[] euclair = loadMaskedEuclair();
        
        List<Set<PairInt>> iceCreamShape = extractNonZeros(icecream);
        List<Set<PairInt>> cupcakeShape = extractNonZeros(cupcake);
        List<Set<PairInt>> gingerBreadmanShape = extractNonZeros(gingerBreadman);
        List<Set<PairInt>> euclairShape = extractNonZeros(euclair);
        
        // -- compare average deltaE's within and between classes
        // -- compare average hsv's within and between classes
        // -- compare intersection of histograms of cieluv cyl angle and magnitude
        //      -- these are large so will look at avgs and stcev
        //         first to see if only 8 divisions in angle and 8 in
        //         magnitude (== hist size 64) would be fine enough
        //         to distinguish between objects.
        // -- compare intersection of hsv histograms
        
        // GroupPixelHSV, GroupPixelRGB, 
        // GroupPixelCIELAB, GroupPixelCIELUV, GroupPixelCIELAB1931
        // GroupPixelCIELCH 
        //    histogramCIECH64
       
        /*
        ----------
        icecream
                img1         img2          img3         (img4)
        rgb
        hsv  mean(+=srfv)
        lab
        lab31
        luv
        lch  mean(+=srdev)
        lch  inter w/ 2     inter w/3      (inter w/ 4)
             inter w/ 3     (inter w/ 4)
             (inter w/ 4)
        hsv  inter w/ 2     inter w/3      (inter w/ 4)
             inter w/ 3     (inter w/ 4)
             (inter w/ 4)
        
        averages of images:
        rgb       
        hsv
        lab
        lab31
        luv
        lch
        and average of all intersections
        -----------
        cupcake
        */
        
        IntraClassColorStats iccs = new IntraClassColorStats();
        
        float[][] deltaELAB = new float[4][2];
        float[][] deltaELAB31 = new float[4][2];
        float[][] deltaELUV = new float[4][2];
        
        float[][] diffHSV = new float[4][3];
        float[][] diffLAB = new float[4][3];
        float[][] diffLAB31 = new float[4][3];
        float[][] diffLUV = new float[4][3];
        float[][] diffLCH = new float[4][3];
        
        // storing average and st dev for each object class
        //    for each band in the color space
        float[][][] meanHSV = new float[4][3][3];
        float[][][] meanLAB = new float[4][3][3];
        float[][][] meanLAB31 = new float[4][3][3];
        float[][][] meanLUV = new float[4][3][3];
        float[][][] meanLCH = new float[4][3][3];
        
        for (int i = 0; i < 4; ++i) {
            ImageExt[] imgs = null;
            List<Set<PairInt>> shape = null;
            String lbl = "";
            switch(i) {
                case 0:
                    imgs = icecream;
                    shape = iceCreamShape;
                    lbl = "icecream";
                    break;
                case 1:
                    imgs = cupcake;
                    shape = cupcakeShape;
                    lbl = "cupcake";
                    break;
                case 2:
                    imgs = gingerBreadman;
                    shape = gingerBreadmanShape;
                    lbl = "gingerbreadman";
                    break;
                default:
                    imgs = euclair;
                    shape = euclairShape;
                    lbl = "euclair";
            }
            deltaELAB[i] = iccs.calculateWithinClassDeltaE(imgs, shape, 
                ColorSpace.CIELAB);
            deltaELAB31[i] = iccs.calculateWithinClassDeltaE(imgs, shape, 
                ColorSpace.CIELAB1931);
            deltaELUV[i] = iccs.calculateWithinClassDeltaE(imgs, shape, 
                ColorSpace.CIELUV);
            
            diffHSV[i] = iccs.calculateWithinClassNormalizedDiffererences(
                imgs, shape, ColorSpace.HSV);
            diffLAB[i] = iccs.calculateWithinClassNormalizedDiffererences(
                imgs, shape, ColorSpace.CIELAB);
            diffLAB31[i] = iccs.calculateWithinClassNormalizedDiffererences(
                imgs, shape, ColorSpace.CIELAB1931);
            diffLUV[i] = iccs.calculateWithinClassNormalizedDiffererences(
                imgs, shape, ColorSpace.CIELUV);
            diffLCH[i] = iccs.calculateWithinClassNormalizedDiffererences(
                imgs, shape, ColorSpace.CIELCH);
        
            meanHSV[i] = iccs.calculateWithinClassMeans(
                imgs, shape, ColorSpace.HSV);
            meanLAB[i] = iccs.calculateWithinClassMeans(
                imgs, shape, ColorSpace.CIELAB);
            meanLAB31[i] = iccs.calculateWithinClassMeans(
                imgs, shape, ColorSpace.CIELAB1931);
            meanLUV[i] = iccs.calculateWithinClassMeans(
                imgs, shape, ColorSpace.CIELUV);
            meanLCH[i] = iccs.calculateWithinClassMeans(
                imgs, shape, ColorSpace.CIELCH);   
        
            // -------- log intra-class means -------
            
            List<StringBuilder> lines = new ArrayList<StringBuilder>();
            lines.add(new StringBuilder("     "));
            lines.add(new StringBuilder(lbl)
                .append("   MEAN of class"));
            lines.add(new StringBuilder("  hsv"));//2
            lines.add(new StringBuilder("     "));//3
            lines.add(new StringBuilder("     "));//4
            lines.add(new StringBuilder("  lab"));//5
            lines.add(new StringBuilder("     "));//6
            lines.add(new StringBuilder("     "));//7
            lines.add(new StringBuilder("lab31"));//8
            lines.add(new StringBuilder("     "));//9
            lines.add(new StringBuilder("     "));//10
            lines.add(new StringBuilder(" luv:"));//11
            lines.add(new StringBuilder("     "));//12
            lines.add(new StringBuilder("     "));//13
            lines.add(new StringBuilder("  lch"));//14
            lines.add(new StringBuilder("     "));//15
            lines.add(new StringBuilder("     "));//16
            
            StringBuilder sb = lines.get(2);
            sb.append(String.format("  %.2f %.2f %.2f  ", 
                meanHSV[i][0][0],meanHSV[i][1][0], meanHSV[i][2][0]));
            sb = lines.get(3);
            sb.append(String.format(" (%.2f %.2f %.2f) ", 
                meanHSV[i][0][1],meanHSV[i][1][1], meanHSV[i][2][1]));
            sb = lines.get(4);
            sb.append(String.format(" (%.2f %.2f %.2f) ", 
                meanHSV[i][0][2],meanHSV[i][1][2], meanHSV[i][2][2]));

            sb = lines.get(5);
            sb.append(String.format("  %.2f %.2f %.2f ", 
                meanLAB[i][0][0],meanLAB[i][1][0], meanLAB[i][2][0]));
            sb = lines.get(6);
            sb.append(String.format(" (%.2f %.2f %.2f) ", 
                meanLAB[i][0][1],meanLAB[i][1][1], meanLAB[i][2][1]));
            sb = lines.get(7);
            sb.append(String.format(" (%.2f %.2f %.2f) ", 
                meanLAB[i][0][2],meanLAB[i][1][2], meanLAB[i][2][2]));

            sb = lines.get(8);
            sb.append(String.format("  %.2f %.2f %.2f ",
                meanLAB31[i][0][0],meanLAB31[i][1][0], meanLAB31[i][2][0]));
            sb = lines.get(9);
            sb.append(String.format(" (%.2f %.2f %.2f) ", 
                meanLAB31[i][0][1],meanLAB31[i][1][1], meanLAB31[i][2][1]));
            sb = lines.get(10);
            sb.append(String.format(" (%.2f %.2f %.2f) ", 
                meanLAB31[i][0][2],meanLAB31[i][1][2], meanLAB31[i][2][2]));

            sb = lines.get(11);
            sb.append(String.format("  %.2f %.2f %.2f ",
                meanLUV[i][0][0],meanLUV[i][1][0], meanLUV[i][2][0]));
            sb = lines.get(12);
            sb.append(String.format(" (%.2f %.2f %.2f) ", 
                meanLUV[i][0][1],meanLUV[i][1][1], meanLUV[i][2][1]));
            sb = lines.get(13);
            sb.append(String.format(" (%.2f %.2f %.2f) ", 
                meanLUV[i][0][2],meanLUV[i][1][2], meanLUV[i][2][2]));

            sb = lines.get(14);
            sb.append(String.format("  %.2f %.2f %.2f ",
                meanLCH[i][0][0],meanLCH[i][1][0], meanLCH[i][2][0]));
            sb = lines.get(15);
            sb.append(String.format(" (%.2f %.2f %.2f) ", 
                meanLCH[i][0][1],meanLCH[i][1][1], meanLCH[i][2][1]));
            sb = lines.get(16);
            sb.append(String.format(" (%.2f %.2f %.2f) ", 
                meanLCH[i][0][2],meanLCH[i][1][2], meanLCH[i][2][2]));
            
            for (int j = 0; j < lines.size(); ++j) {
                System.out.println(lines.get(j).toString());
            }
            
            // -------- log intra-class differences -------
            lines = new ArrayList<StringBuilder>();
            lines.add(new StringBuilder("     "));
            lines.add(new StringBuilder(lbl)
                .append("   DIFF within class"));
            lines.add(new StringBuilder("  hsv"));//2
            lines.add(new StringBuilder("     "));//3
            lines.add(new StringBuilder("  lab"));//4
            lines.add(new StringBuilder("     "));//5
            lines.add(new StringBuilder("lab31"));//6
            lines.add(new StringBuilder("     "));//7
            lines.add(new StringBuilder(" luv:"));//8
            lines.add(new StringBuilder("     "));//9
            lines.add(new StringBuilder("  lch"));//10
            lines.add(new StringBuilder("     "));//11
          
            sb = lines.get(2);
            sb.append(String.format("  %.2f  ", 
                diffHSV[i][0]));
            sb = lines.get(3);
            sb.append(String.format(" (%.2f, %.2f) ", 
                diffHSV[i][1], diffHSV[i][2]));

            sb = lines.get(4);
            sb.append(String.format("  %.2f  ", 
                diffLAB[i][0]));
            sb = lines.get(5);
            sb.append(String.format(" (%.2f, %.2f) ", 
                diffLAB[i][1], diffLAB[i][2]));;

            sb = lines.get(6);
            sb.append(String.format("  %.2f  ", 
                diffLAB31[i][0]));
            sb = lines.get(7);
            sb.append(String.format(" (%.2f, %.2f) ", 
                diffLAB31[i][1], diffLAB31[i][2]));;

            sb = lines.get(8);
            sb.append(String.format("  %.2f  ", 
                diffLUV[i][0]));
            sb = lines.get(9);
            sb.append(String.format(" (%.2f, %.2f) ", 
                diffLUV[i][1], diffLUV[i][2]));;

            sb = lines.get(10);
            sb.append(String.format("  %.2f  ", 
                diffLCH[i][0]));
            sb = lines.get(11);
            sb.append(String.format(" (%.2f, %.2f) ", 
                diffLCH[i][1], diffLCH[i][2]));;
            
            for (int j = 0; j < lines.size(); ++j) {
                System.out.println(lines.get(j).toString());
            }

            // ---- log the intra-class deltaEs
            
            lines = new ArrayList<StringBuilder>();
            lines.add(new StringBuilder("     "));
            lines.add(new StringBuilder(lbl)
               .append("   deltaE within class"));
            lines.add(new StringBuilder("  lab"));//2
            lines.add(new StringBuilder("     "));//3
            lines.add(new StringBuilder("lab31"));//4
            lines.add(new StringBuilder("     "));//5
            lines.add(new StringBuilder(" luv:"));//6
            lines.add(new StringBuilder("     "));//7
          
            sb = lines.get(2);
            sb.append(String.format("  %.2f  ", deltaELAB[i][0]));
            sb = lines.get(3);
            sb.append(String.format(" (%.2f) ",  deltaELAB[i][1]));

            sb = lines.get(4);
            sb.append(String.format("  %.2f  ", deltaELAB31[i][0]));
            sb = lines.get(5);
            sb.append(String.format(" (%.2f) ",  deltaELAB31[i][1]));

            sb = lines.get(6);
            sb.append(String.format("  %.2f  ", deltaELUV[i][0]));
            sb = lines.get(7);
            sb.append(String.format(" (%.2f) ",  deltaELUV[i][1]));
            
            for (int j = 0; j < lines.size(); ++j) {
                System.out.println(lines.get(j).toString());
            }            
            
            lines = null;
        }
            
        System.out.println("--- between class stats ----");
       
        BetweenClassColorStats bccs = new BetweenClassColorStats();
        AllClassInterStats hsvInterStats = 
            bccs.calculateDiffBetweenClasses(meanHSV);
        System.out.println("HSV min signif diff=" 
            + hsvInterStats.findMinSignficance()
            + "\n  min of avg sig="
            + hsvInterStats.findMinOfAverageSignficance()
        );
        
        AllClassInterStats labInterStats = 
            bccs.calculateDiffBetweenClasses(meanLAB);
        System.out.println("LAB min significant diff=" 
            + labInterStats.findMinSignficance()
            + "\n  min of avg sig="
            + labInterStats.findMinOfAverageSignficance()
        );
        
        AllClassInterStats lab31InterStats = 
            bccs.calculateDiffBetweenClasses(meanLAB31);
        System.out.println("LAB1931 min significant diff=" 
            + lab31InterStats.findMinSignficance()
            
            + "\n  min of avg sig="
            + lab31InterStats.findMinOfAverageSignficance()
        );
    
        AllClassInterStats luvInterStats = 
            bccs.calculateDiffBetweenClasses(meanLUV);
        System.out.println("LUV min significant diff=" 
            + luvInterStats.findMinSignficance()
            
            + "\n  min of avg sig="
            + luvInterStats.findMinOfAverageSignficance()
        );
        
        
        AllClassInterStats lchInterStats = 
            bccs.calculateDiffBetweenClasses(meanLCH);
        System.out.println("LCH min significant diff=" 
            + lchInterStats.findMinSignficance()
            
            + "\n  min of avg sig="
            + lchInterStats.findMinOfAverageSignficance()
        );
        
    
        AllClassInterStats labDeltaEInterStats = 
            bccs.calculateDeltaEBetweenClasses(meanLAB, deltaELAB);
        System.out.println("LAB DeltaE min significant diff=" 
            + labDeltaEInterStats.findMinSignficance()
            
            + "\n  min of avg sig="
            + labInterStats.findMinOfAverageSignficance()
        );
        
        
        AllClassInterStats lab31DeltaEInterStats = 
            bccs.calculateDeltaEBetweenClasses(meanLAB31, deltaELAB31);
        System.out.println("LAB31 DeltaE min significant diff=" 
            + lab31DeltaEInterStats.findMinSignficance()
            
            + "\n  min of avg sig="
            + lab31InterStats.findMinOfAverageSignficance()
        );
        
        
        AllClassInterStats luvDeltaEInterStats = 
            bccs.calculateDeltaEBetweenClasses(meanLUV, deltaELUV);
        System.out.println("LUV DeltaE min significant diff=" 
            + luvDeltaEInterStats.findMinSignficance()
            
            + "\n  min of avg sig="
            + luvInterStats.findMinOfAverageSignficance()
        );
        
    
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
}
