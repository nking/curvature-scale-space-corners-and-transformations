package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.*;
import algorithms.imageProcessing.features.mser.MSEREdges;
import algorithms.imageProcessing.features.mser.MSEREdgesWrapper;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.imageProcessing.segmentation.MSEREdgesToLabelsHelper;
import algorithms.imageProcessing.segmentation.MergeLabels;
import algorithms.imageProcessing.segmentation.SLICSuperPixels;
import algorithms.imageProcessing.segmentation.MergeLabels.METHOD;
import algorithms.misc.MiscDebug;
import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import gnu.trove.set.TIntSet;
import junit.framework.TestCase;

import java.io.*;
import java.io.IOException;
import java.nio.file.FileSystems;
import java.util.*;
import java.util.logging.Logger;

public class MultiPartialShapeMatcherTest extends TestCase {

    static String sep = FileSystems.getDefault().getSeparator();
    static String eol = System.lineSeparator();
    Logger log = Logger.getLogger(MultiPartialShapeMatcherTest.class.getName());

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

            int[] labels = mergedSLIC(img, fileName1Root);



            //write_centroids(labels, img, fileName1Root);


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

    private int[] mergedMSER(ImageExt img, String fileName1Root) throws IOException {

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

        return labels2;
    }

    private int[] mergedSLIC(ImageExt img, String fileName1Root) throws IOException {

        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(img, SIGMA.getValue(SIGMA.TWO));

        METHOD method = METHOD.MIN_GRADIENT;
        //int nc = 100; double thresh = 2.5;
        //int nc = 70; double thresh = 2.5;
        //int nc = 50; double thresh = 3.0;//2nd best
        int nc = 40; double thresh = 1.0; // best results

        //METHOD method = METHOD.MODE;
        //int nc = 100; double thresh = 5.0;
        //int nc = 70; double thresh = 5.0;

        int[] labels2 = slic(img.copyToImageExt(), fileName1Root, nc);

        int nLabels2 = MergeLabels.mergeUsingDeltaE2000(img, labels2, thresh, method);
        System.out.printf("merged k=%d, nc=%d\n", nLabels2, nc);
        ImageExt im3 = img.copyToImageExt();
        ImageIOHelper.addAlternatingColorLabelsToRegion(im3, labels2);
        MiscDebug.writeImage(im3, "_" + fileName1Root + "_slic_merged_");

        //return new int[][]{labels, labels2};
        return labels2;
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

    private void write_centroids(int[] labels, ImageExt img, String fileName1Root) throws IOException {
        int[][] centroids = MergeLabels.calculateCentroids(labels, img.getWidth());

        // write to bin directory for use with python script to scale the
        // image and coordinates for input into a SAM model where SAM is
        // "Segment Anything Model" hosted on Google's Colab.
        String coordsFilePath = ResourceFinder.getBaseDir() + sep + "bin" + sep + fileName1Root + "_points.txt";
        BufferedWriter b = null;
        FileWriter w = null;
        try {
            w = new FileWriter(coordsFilePath);
            b = new BufferedWriter(w);
            for (int j = 0; j < centroids.length; ++j) {
                w.write(String.format("%d %d%s", centroids[j][0], centroids[j][1], eol));
            }
        } catch (IOException ex) {
            log.severe("Error: " + ex.getMessage());
        } finally {
            if (w != null) {
                w.close();
            }
            if (b != null) {
                b.close();
            }
        }
        // SAM model:
        //https://keras.io/keras_hub/guides/segment_anything_in_keras_hub/
        // B is the number of images in a batch
        //"points": A batch of point prompts. Each point is an (x, y) coordinate originating from the top-left
        // corner of the image. In other works, each point is of the form (r, c) where r and c are the row and
        // column of the pixel in the image. Must be of shape (B, N, 2).
        //    ==> presumably N can be different for each image


        //"images": A batch of images to segment. Must be of shape (B, 1024, 1024, 3).

        //Segment Anything allows prompting an image using points, boxes, and masks.
        //so will get points from the centroids of labeled regions

        //Point prompts are the most basic of all: the model tries to guess the object given a point on an image. The point can either be a foreground point (i.e. the desired segmentation mask contains the point in it) or a backround point (i.e. the point lies outside the desired mask).

        /*
        https://colab.research.google.com/drive/1ITKfey9d7MVQ4H32BbNl1jmCFGMyt--H
        inside the SAM code on colab:
        from PIL import Image

        !wget -q https://raw.githubusercontent.com/nking/curvature-scale-space-corners-and-transformations/main/testresources/android_statues_02.jpg
        !wget -q https://raw.githubusercontent.com/nking/curvature-scale-space-corners-and-transformations/main/testresources/android_statues_03.jpg
        !wget -q https://raw.githubusercontent.com/nking/curvature-scale-space-corners-and-transformations/main/testresources/android_statues_04.jpg
        !wget -q https://raw.githubusercontent.com/nking/curvature-scale-space-corners-and-transformations/main/testresources/android_statues_01.jpg

        TODO: same for test data

        for file_name in ["android_statues_02", "android_statues_03", "android_statues_04", "android_statues_01"]:
            image = np.array(keras.utils.load_img(f"{file_name}.jpg"))
            image = inference_resizing(image)

            #test point: col=533, row=193
            read the centroid array into format like this [x,y] w.r.t upper left corner of image

            input_point = np.array([[848, 270], [522, 118], [522,260], [209, 405], [278, 230]])
            input_label = np.array([1,1,1,1,1])

            #plt.figure(figsize=(10, 10))
            #plt.imshow(ops.convert_to_numpy(image) / 255.0)
            #show_points(input_point, input_label, plt.gca())
            #plt.axis("on")
            #plt.show()

            outputs = model.predict(
                {
                    "images": image[np.newaxis, ...],
                    "points": np.concatenate(
                        [input_point[np.newaxis, ...], np.zeros((1, 1, 2))], axis=1
                    ),
                    "labels": np.concatenate(
                        [input_label[np.newaxis, ...], np.full((1, 1), fill_value=-1)], axis=1
                    ),
                }
            )
            fig, ax = plt.subplots(1, 3, figsize=(20, 60))
            masks, scores = outputs["masks"][0][1:], outputs["iou_pred"][0][1:]
            for i, (mask, score) in enumerate(zip(masks, scores)):
                mask = inference_resizing(mask[..., None], pad=False)[..., 0]
                mask, score = map(ops.convert_to_numpy, (mask, score))
                mask = 1 * (mask > 0.0)
                #ax[i].imshow(ops.convert_to_numpy(image) / 255.0)
                #show_mask(mask, ax[i])
                #show_points(input_point, input_label, ax[i])
                #ax[i].set_title(f"Mask {i+1}, Score: {score:.3f}", fontsize=12)
                #ax[i].axis("off")
                #save masks as black and white, 1 bit images
                im_mask = Image.fromarray(mask.astype(np.uint8), mode='L').convert('1')
                im_mask.save(f'./{file_name}_mask_{i}.png')
                with open(f'./{file_name}_score_{i}.txt', "w") as file:
                    file.write(f"{score}\n")
            #plt.show()
         */
    }

}
