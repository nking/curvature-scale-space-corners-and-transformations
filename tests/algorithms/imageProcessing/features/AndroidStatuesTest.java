package algorithms.imageProcessing.features;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.compGeometry.RotatedOffsets;
import algorithms.compGeometry.clustering.KMeansPlusPlus;
import algorithms.compGeometry.clustering.KMeansPlusPlusColor;
import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.DFSContiguousIntValueFinder;
import algorithms.imageProcessing.DFSContiguousValueFinder;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.GroupPixelColors;
import algorithms.imageProcessing.GroupPixelRGB;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageProcessor.Colors;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.ImageSegmentation.DecimatedData;
import algorithms.imageProcessing.PixelColors;
import algorithms.imageProcessing.matching.PartialShapeMatcher;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.SegmentationMergeThreshold;
import algorithms.imageProcessing.segmentation.LabelToColorHelper;
import algorithms.imageProcessing.segmentation.NormalizedCuts;
import algorithms.imageProcessing.segmentation.SLICSuperPixels;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.util.AngleUtil;
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
import java.util.Arrays;
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

    public void testShapeMatcher() throws Exception {

        int maxDimension = 512;
        int nClusters = 200;
        //int clrNorm = 5;

        ImageProcessor imageProcessor = new ImageProcessor();
        ImageSegmentation imageSegmentation = new
            ImageSegmentation();
        CIEChromaticity cieC = new CIEChromaticity();

        String fileNameMask0
            = "android_statues_03_sz2_mask.png";
        String filePathMask0 = ResourceFinder
            .findFileInTestResources(fileNameMask0);
        ImageExt imgMask0 = ImageIOHelper.readImageExt(filePathMask0);

        String fileName0
            = "android_statues_03.jpg";
        int idx = fileName0.lastIndexOf(".");
        String fileName0Root = fileName0.substring(0, idx);
        String filePath0 = ResourceFinder
            .findFileInTestResources(fileName0);
        ImageExt img0 = ImageIOHelper.readImageExt(filePath0);
        int w0 = img0.getWidth();
        int h0 = img0.getHeight();
        int binFactor0 = (int) Math.ceil(Math.max(
            (float) w0 / maxDimension, (float) h0 / maxDimension));
        img0 = imageProcessor.binImage(img0, binFactor0);

        assert(img0.getNPixels() == imgMask0.getNPixels());

        Set<PairInt> shape0 = new HashSet<PairInt>();
        // multiply img0 by imgMask0 to leave only the pixels
        // of the model shape in the image.
        for (int i = 0; i < img0.getNPixels(); ++i) {
            if (imgMask0.getR(i) == 0) {
                img0.setRGB(i, 0, 0, 0);
            } else {
                shape0.add(new PairInt(img0.getCol(i),
                   img0.getRow(i)));
            }
        }

        //PairIntArray p = extractOrderedBoundary(imgMask0, SIGMA.ONE);
        //plot(p, 100);

        // extract color stats for img0
        GroupPixelRGB rgb0 = new GroupPixelRGB(shape0, img0, 0, 0);
        float[] lab0 = cieC.rgbToCIELAB(
            Math.round(rgb0.getAvgRed()),
            Math.round(rgb0.getAvgGreen()),
            Math.round(rgb0.getAvgBlue()));

        // -------

        String fileNameRoot1 = "";

        for (int i = 1; i < 2; ++i) {
            switch(i) {
                case 0: {
                    fileNameRoot1 = "android_statues_01";
                    break;}
                case 1: {
                    fileNameRoot1 = "android_statues_02";
                    break;}
                case 2: {
                    fileNameRoot1 = "android_statues_03";
                    break;}
                case 3: {
                    fileNameRoot1 = "android_statues_04";
                    break;}
                default: { break;}
            }

            String fileNameMask1 = fileNameRoot1 + "_sz1_mask.png";
            String fileName1 = fileNameRoot1 + ".jpg";

            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
            int w1 = img1.getWidth();
            int h1 = img1.getHeight();
            int binFactor1 = (int) Math.ceil(Math.max((float) w1 / maxDimension,
                (float) h1 / maxDimension));
            img1 = imageProcessor.binImage(img1, binFactor1);

            String filePathMask1 = ResourceFinder.findFileInTestResources(fileNameMask1);
            ImageExt imgMask1 = ImageIOHelper.readImageExt(filePathMask1);

            ImageExt img1Cp = img1.copyToImageExt();

            // -- use superpixels on img1
            // -- then cluster the super pixels by polar cie xy
            // -- filter superpixel labels for model colors

            SLICSuperPixels slic
                = new SLICSuperPixels(img1, nClusters);
            slic.calculate();
            int[] labels = slic.getLabels();

            ImageExt img1Labeled = img1Cp.copyToImageExt();
            ImageExt img1LabeledAlt = img1Cp.copyToImageExt();
            LabelToColorHelper.applyLabels(img1Labeled, labels);
            ImageIOHelper.addAlternatingColorLabelsToRegion(
                img1LabeledAlt, labels);
            MiscDebug.writeImage(img1Labeled,  "_slic_" + fileNameRoot1);
            MiscDebug.writeImage(img1LabeledAlt,  "_slic_alt_" + fileNameRoot1);

            ImageExt img1Cp2 = img1Cp.createWithDimensions();
            ImageExt img1Cp3 = img1Cp.copyToImageExt();
            ImageExt img1Cp4 = img1Cp.copyToImageExt();

            List<Set<PairInt>> filtered
                = new ArrayList<Set<PairInt>>();

            int nExtraForDot = 0;

            List<Set<PairInt>> contiguousSets = LabelToColorHelper
                .extractContiguousLabelPoints(img1Labeled, labels);
            List<TIntList> mergedContigIndexes =
                imageSegmentation.
                mergeUsingPolarCIEXYAndFrequency(img1Labeled,
                contiguousSets, 0.1f);
            List<Set<PairInt>> clusterSets1S = new
                ArrayList<Set<PairInt>>();
            for (TIntList list : mergedContigIndexes) {
                Set<PairInt> set = new HashSet<PairInt>();
                for (int ii = 0; ii < list.size(); ++ii) {
                    int cIdx = list.get(ii);
                    set.addAll(contiguousSets.get(cIdx));
                }
                clusterSets1S.add(set);
            }

            for (int j = 0; j < clusterSets1S.size(); ++j) {
                int[] rgb = ImageIOHelper.getNextRGB(j);
                Set<PairInt> set = clusterSets1S.get(j);

              //  imageSegmentation.erode(set, 2);

                ImageIOHelper.addToImage(set, 0, 0, img1Cp2,
                    nExtraForDot, rgb[0], rgb[1], rgb[2]);
            }
            MiscDebug.writeImage(img1Cp2,
                "_slic_pclstr_" + fileNameRoot1);

            // -- combine information in labels
            //    and clusterSets1S to make better segmentation
            //    (NOTE: the clustering algorithm could better
            //     incorporate this too...needs improvements one day)
            sortByDecrSize(clusterSets1S);

            TIntObjectMap<Set<PairInt>> labelMap =
                LabelToColorHelper.extractLabelPoints(
                img1Cp3,
                labels);

            List<Set<PairInt>> segmentedList = new ArrayList<Set<PairInt>>();
            TIntSet mergedIndexes = new TIntHashSet();
            for (Set<PairInt> set : clusterSets1S) {
                Set<PairInt> set2 = new HashSet<PairInt>();
                for (PairInt p : set) {
                    int pixIdx = img1.getInternalIndex(p.getX(), p.getY());
                    int lIdx = labels[pixIdx];
                    if (!mergedIndexes.contains(lIdx)) {
                        set2.addAll(labelMap.get(lIdx));
                        mergedIndexes.add(lIdx);
                    }
                }
                segmentedList.add(set2);
            }

            for (int j = 0; j < segmentedList.size(); ++j) {
                int[] rgb = ImageIOHelper.getNextRGB(j);
                Set<PairInt> set = segmentedList.get(j);
                ImageIOHelper.addToImage(set, 0, 0, img1Cp4,
                    nExtraForDot, rgb[0], rgb[1], rgb[2]);
            }
            MiscDebug.writeImage(img1Cp4,
                "_segmented_" + fileNameRoot1);

            for (int j = 0; j < segmentedList.size(); ++j) {

                Set<PairInt> points = segmentedList.get(j);

                GroupPixelRGB rgb1
                    = new GroupPixelRGB(points, img1, 0, 0);

                float[] lab1 = cieC.rgbToCIELAB(
                    Math.round(rgb1.getAvgRed()),
                    Math.round(rgb1.getAvgGreen()),
                    Math.round(rgb1.getAvgBlue()));

                double deltaE = cieC.calcDeltaECIE2000(
                    lab0[0], lab0[1], lab0[2],
                    lab1[0], lab1[1], lab1[2]);
                if (deltaE < 0) {
                    deltaE *= -1;
                }
                if (deltaE <= 8.0) {
                    filtered.add(points);
                }
            }
            for (Set<PairInt> set : filtered) {
                for (PairInt pt : set) {
                    img1.setRGB(pt.getX(), pt.getY(), 255, 255, 255);
                }
            }
            MiscDebug.writeImage(img1,  "_filtered_" + fileNameRoot1);

            // -- make one pass over filtered to fit against
            //    template.
            //    -- for the best fitting,
            //       try to aggregate adjacent sets to
            //       obtain a better fit to template
            //       (adjacent sets could be found in
            //        many different ways)

            /*
            PairIntArray q = extractOrderedBoundary(img);
            plot(q, fileNumber);

            int dp = 2;
            PartialShapeMatcher matcher =
                new PartialShapeMatcher();
            matcher.setToDebug();
            matcher.overrideSamplingDistance(dp);
            PartialShapeMatcher.Result result = matcher.match(p, q);
            */

            /*
            MiscDebug.writeImage(img,  "_img_" + fileName1Root);

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

            int nClusters = 200;//100;
            //int clrNorm = 5;
            SLICSuperPixels slic
                = new SLICSuperPixels(img, nClusters);
            slic.calculate();
            int[] labels = slic.getLabels();
            ImageIOHelper.addAlternatingColorLabelsToRegion(
            //LabelToColorHelper.applyLabels(
                img, labels);
            MiscDebug.writeImage(img, "_slic_" + fileName1Root);

            /*
            NormalizedCuts normCuts = new NormalizedCuts();
            normCuts.setColorSpaceToHSV();
            int[] labels2 = normCuts.normalizedCut(img, labels);
            labels = labels2;
            ImageIOHelper.addAlternatingColorLabelsToRegion(
                img, labels);
            MiscDebug.writeImage(img, "_norm_cuts_" + fileName1Root);
            */

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
            blob, sigma);

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

}
