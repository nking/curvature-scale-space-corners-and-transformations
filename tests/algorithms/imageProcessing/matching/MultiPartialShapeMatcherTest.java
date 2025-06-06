package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.*;
import algorithms.imageProcessing.features.PhaseCongruencyDetector2;
import algorithms.imageProcessing.features.mser.MSER;
import algorithms.imageProcessing.features.mser.MSEREdges;
import algorithms.imageProcessing.features.mser.MSEREdgesWrapper;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.imageProcessing.segmentation.MSEREdgesToLabelsHelper;
import algorithms.imageProcessing.segmentation.MergeLabels;
import algorithms.imageProcessing.segmentation.NormalizedCuts;
import algorithms.imageProcessing.segmentation.SLICSuperPixels;
import algorithms.misc.MiscDebug;
import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import gnu.trove.set.TIntSet;
import junit.framework.TestCase;

import java.io.IOException;
import java.util.*;

public class MultiPartialShapeMatcherTest extends TestCase {

    public void testAndroidStatues() throws IOException {

        String fileName0 = "android_statues_03_sz1_mask_small.png";
        int idx0 = fileName0.lastIndexOf(".");
        String fileName0Root = fileName0.substring(0, idx0);
        String filePath0 = ResourceFinder.findFileInTestResources(fileName0);
        ImageExt img0 = ImageIOHelper.readImageExt(filePath0);

        //PairFloatArray p = extractOrderedBoundary(img0);

        for (int i = 0; i < 4; ++i) {
            String fileName1;
            switch (i) {
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
                default: {
                    fileName1 = "android_statues_04.jpg";
                    break;
                }
            }

            //int minNumPts = p.getN();
            int minNumPts = 10;

            int idx = fileName1.lastIndexOf(".");
            String fileName1Root = fileName1.substring(0, idx);
            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            ImageExt img = ImageIOHelper.readImageExt(filePath1);

            int[][] labelsMserAndSlic = mserAndSLIC_2(img, fileName1Root);

            //TODO: get the centers of the labels.
            // will use them with yolov11 or yolov12


            //extractShapes(img, fileName1Root);

            /*
            List<Set<PairInt>> blobs = createBlobs(labels, img.getWidth(), minNumPts);
            List<PairFloatArray> cc = new ArrayList<>();
            for (Set<PairInt> item : blobs) {
                PairFloatArray q = extractOrderedBoundary(item, SIGMA.ONE, img.getWidth(), img.getHeight());
                if (q != null) {
                    cc.add(q);
                }
            }
            plot(img, cc, fileName1Root + "_closed_curves_" + nc + "_");
            */

            // use MultiPartialShapeMatcher to match these
        }

    }

    private int[][] mserAndSLIC(ImageExt img, String fileName1Root) throws IOException {

        ImageProcessor imageProcessor = new ImageProcessor();
        //imageProcessor.blur(img, SIGMA.getValue(SIGMA.TWO));

        MSEREdgesWrapper msew = new MSEREdgesWrapper();
        MSEREdges mserE = msew.extractAndMergeEdges(img, 512);
        //MSEREdges mserE = new MSEREdges(img);
        //mserE.extractAndMergeEdges();

        List<TIntSet> edgeList = mserE.getEdges();
        Image im = mserE.getGsImg().copyToColorGreyscale();
        int[] clr = new int[]{255, 0, 0};
        for (int ii = 0; ii < edgeList.size(); ++ii) {
            ImageIOHelper.addCurveToImage(edgeList.get(ii), im, 0, clr[0],clr[1], clr[2]);
        }
        MiscDebug.writeImage(im, "_" + fileName1Root + "_mser_edges_");

        //assign labels
        ImageExt im2 = msew.binImage(img);
        assertEquals(im.getNPixels(), im2.getNPixels());
        assertEquals(im.getWidth(), im2.getWidth());
        assertEquals(im.getHeight(), im2.getHeight());
        int[] labels = new int[im2.getNPixels()];
        int nComp = MSEREdgesToLabelsHelper.createLabels(im2, edgeList, labels);
        assertTrue(MSEREdgesToLabelsHelper.allAreNonNegative(labels));
        ImageIOHelper.addAlternatingColorLabelsToRegion(im2, labels);
        MiscDebug.writeImage(im2, "_" + fileName1Root + "_mser_segmentation_");


        int nc = 50;
        int[] labels2 = slic(msew.binImage(img), fileName1Root, nc);

        return new int[][]{labels, labels2};
    }

    private int[][] mserAndSLIC_2(ImageExt img, String fileName1Root) throws IOException {

        ImageProcessor imageProcessor = new ImageProcessor();
        //imageProcessor.blur(img, SIGMA.getValue(SIGMA.TWO));

        /*
        MSEREdges mserE = new MSEREdges(img);
        mserE.extractAndMergeEdges();

        List<TIntSet> edgeList = mserE.getEdges();
        Image im = mserE.getGsImg().copyToColorGreyscale();
        int[] clr = new int[]{255, 0, 0};
        for (int ii = 0; ii < edgeList.size(); ++ii) {
            ImageIOHelper.addCurveToImage(edgeList.get(ii), im, 0, clr[0],clr[1], clr[2]);
        }
        MiscDebug.writeImage(im, "_" + fileName1Root + "_mser_edges_");

        //assign labels
        ImageExt im2 = img.copyToImageExt();
        assertEquals(im.getNPixels(), im2.getNPixels());
        assertEquals(im.getWidth(), im2.getWidth());
        assertEquals(im.getHeight(), im2.getHeight());
        int[] labels = new int[im2.getNPixels()];
        int nComp = MSEREdgesToLabelsHelper.createLabels(im2, edgeList, labels);
        assertTrue(MSEREdgesToLabelsHelper.allAreNonNegative(labels));
        ImageIOHelper.addAlternatingColorLabelsToRegion(im2, labels);
        MiscDebug.writeImage(im2, "_" + fileName1Root + "_mser_segmentation_");
        */

        //int nc = 100; double thresh = 0.55;
        int nc = 70; double thresh = 0.20;

        int[] labels2 = slic(img.copyToImageExt(), fileName1Root, nc);

        int nLabels2 = MergeLabels.mergeUsingDeltaE2000(img, labels2, thresh);
        ImageExt im3 = img.copyToImageExt();
        ImageIOHelper.addAlternatingColorLabelsToRegion(im3, labels2);
        MiscDebug.writeImage(im3, "_" + fileName1Root + "_slic_merged_");

        //return new int[][]{labels, labels2};
        return new int[][]{labels2};
    }

    public void extractShapes(ImageExt img, String fileName1Root) throws IOException {

        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(img, SIGMA.getValue(SIGMA.ONE));

        MSEREdgesWrapper msew = new MSEREdgesWrapper();
        //msew.setToDebug();
        MSEREdges mserE = msew.extractAndMergeEdges(img);

        List<TIntSet> edgeList = mserE.getEdges();
        Image im = mserE.getGsImg().copyToColorGreyscale();
        int[] clr = new int[]{255, 0, 0};
        for (int ii = 0; ii < edgeList.size(); ++ii) {
            ImageIOHelper.addCurveToImage(edgeList.get(ii), im, 0, clr[0],
                    clr[1], clr[2]);
        }
        MiscDebug.writeImage(im, "_" + fileName1Root + "_mser_edges_");

        List<Region> regions = mserE.getRegions();
        List<List<Region>> regionsList = new ArrayList<>();
        regionsList.add(regions);
        regionsList.add(new ArrayList<Region>());
        Image img2 = msew.binImage(img).copyImage();
        Region.drawEllipses(img2, regionsList, 0);
        MiscDebug.writeImage(img2, "_" + fileName1Root + "_mser_regions_");
    }

    private void plot(ImageExt img, List<PairFloatArray> shapes,
                      String fileSuffix) throws IOException {

        Image im = img.copyToGreyscale().copyToColorGreyscale();
        int[] clr;
        int nDot = 0;
        for (int ii = 0; ii < shapes.size(); ++ii) {
            clr = ImageIOHelper.getNextRGB(ii);
            ImageIOHelper.addCurveToImage(shapes.get(ii), im, nDot, clr[0],
                    clr[1], clr[2]);
        }
        MiscDebug.writeImage(im, fileSuffix);

    }

    private List<Set<PairInt>> createBlobs(int[] labels, int width, int minNumPts) {
        Map<Integer, Set<Integer>> map = new HashMap<>();
        for (int i = 0; i < labels.length; ++i) {
            map.putIfAbsent(labels[i], new HashSet<>());
            map.get(labels[i]).add(i);
        }
        List<Set<PairInt>> out = new ArrayList<>();
        for (int label : map.keySet()) {
            if (map.get(label).size() < minNumPts) {
                continue;
            }
            Set<PairInt> set = new HashSet<>();
            for (int c : map.get(label)) {
                set.add(new PairInt(c % width, c / width));
            }
            out.add(set);
        }
        return out;
    }

    private int[] slic(ImageExt img, String fileName1Root, int nc) throws IOException {

        ImageSegmentation imS = new ImageSegmentation();

        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(img, SIGMA.getValue(SIGMA.TWO));

        int method = 0;

        EdgeFilterProducts edgeProducts = imS.createGradient(img, method, System.currentTimeMillis());

        SLICSuperPixels slic = new SLICSuperPixels(img, nc);
        slic.setGradient(edgeProducts.getGradientXY());
        slic.calculate();
        int[] labels = slic.getLabels();

        Image img3 = img.createWithDimensions();
        ImageIOHelper.addAlternatingColorLabelsToRegion(img3, labels);
        String str = Integer.toString(nc);
        while (str.length() < 3) {
            str = "0" + str;
        }
        String suffix = String.format("_slic_%s_%s", fileName1Root, str);
        MiscDebug.writeImage(img3, suffix);
        return labels;

    }

    private PairFloatArray extractOrderedBoundary(ImageExt image) {
        return extractOrderedBoundary(image, SIGMA.ONE);
    }

    private PairFloatArray extractOrderedBoundary(ImageExt image,
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

        return extractOrderedBoundary(blob, sigma, img.getWidth(), img.getHeight());
    }

    private PairFloatArray extractOrderedBoundary(Set<PairInt> blob,
                                                  SIGMA sigma, int width, int height) {

        try {
            ImageProcessor imageProcessor = new ImageProcessor();

            PairIntArray ordered =
                    imageProcessor.extractSmoothedOrderedBoundary(
                            blob, sigma, width, height);

            PairFloatArray out = new PairFloatArray();
            for (int i = 0; i < ordered.getN(); ++i) {
                out.add(ordered.getX(i), ordered.getY(i));
            }

            return out;
        } catch (Exception ex) {
            return null;
        }
    }

}
