package algorithms.imageProcessing.matching;

import algorithms.compGeometry.PerimeterFinder3;
import algorithms.imageProcessing.*;
import algorithms.imageProcessing.features.mser.MSEREdges;
import algorithms.imageProcessing.features.mser.MSEREdgesWrapper;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.imageProcessing.segmentation.*;
import algorithms.imageProcessing.segmentation.MergeLabels.METHOD;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.*;
import com.spotify.voyager.jni.Index;
import gnu.trove.set.TIntSet;
import junit.framework.TestCase;

import java.io.*;
import java.io.IOException;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.logging.Logger;
import java.nio.file.DirectoryStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

public class MultiPartialShapeMatcherTest extends TestCase {

    static String sep = FileSystems.getDefault().getSeparator();
    static String eol = System.lineSeparator();
    Logger log = Logger.getLogger(MultiPartialShapeMatcherTest.class.getName());

    public void testDB() {
        // ns is length of an embedding
        int[] ns = new int[]{20, 15, 10};
        long seed = System.nanoTime();
        seed = 272458813734000L;
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);

        double max = 2*Math.PI;
        int nEmb = 10; // number of embeddings, that is, number of samples
        float tol = 0.001f;

        for (int i = 0; i < ns.length; ++i) {
            float[][] embeddings = new float[nEmb][ns[i]];
            float[] q = new float[ns[i]];
            float[][] qs = new float[nEmb][ns[i]];
            long[] ids = new long[nEmb];
            int topK = nEmb;
            for (int j = 0; j < nEmb; ++j) {
                embeddings[j] = new float[ns[i]];
                for (int k = 0; k < ns[i]; ++k) {
                    embeddings[j][k] = (float) rand.nextDouble(max);
                }
                ids[j] = j;
            }
            for (int k = 0; k < ns[i]; ++k) {
                q[k] = (float) rand.nextDouble(max);
            }

            Index index = new Index(Index.SpaceType.Euclidean, ns[i]);

            long[] outIds = index.addItems(embeddings, ids, -1);

            float[] expectedDist = new float[nEmb];
            for (int j = 0; j < nEmb; ++j) {
                expectedDist[j] = euclidSq(q, embeddings[j]);
            }
            Index.QueryResults res = index.query(q, topK);
            for (int j = 0; j < res.getLabels().length; ++j) {
                float d = res.getDistances()[j];
                int idx = (int)res.getLabels()[j];
                float expected = expectedDist[idx];
                assertTrue(Math.abs(expected - d) <= tol);
            }
        }
        //TODO: multiple queries
    }

    protected float euclidSq(float[] a, float[] b) {
        double sum = 0;
        double diff;
        for (int i = 0; i < a.length; ++i) {
            diff = a[i] - b[i];
            sum += (diff * diff);
        }
        return (float)sum;
    }

    public void testMPEG7_0() throws Exception {
        /** test using the MPEG-7 shape dataset

        https://dabi.temple.edu/external/shape/MPEG7/dataset.html

         each class (category) has several contours of the same object with small differences in articulation or occlusion,
         and more.

         testing whether one from each class can be randomly extracted and placed in a database, and then
         all other contours become queries to the database to see if the same class is found.

         the method tested returns the topk best hits even if it returns the same contour more than once,
         but different offsets or lengths.
        */

        long seed = System.nanoTime();
        seed = 211649731957458L;
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);

        // the randomly chosen single class member from each class:
        List<PairFloatArray> queries = new ArrayList<>();
        List<PairFloatArray> targets = new ArrayList<>();

        // key = target index, value = file_name
        Map<Integer, String> indexFileNameMap = new HashMap<>();
        Map<Integer, String> queryFileNameMap = new HashMap<>();

        // === load the contours ====
        String dataPath = ResourceFinder.findTestResourcesDirectory() + sep  +"mpeg7";
        File dataDir = new File(dataPath);

        String[] classes = Files.lines(Paths.get(dataPath + sep + "CLASSES.txt")).toArray(String[]::new);

        // for each class, load the contours into queries or targets, depending upon random choic of query
        for (int i = 0; i < classes.length; ++i) {
            final String category = classes[i];
            File[] contourFiles = dataDir.listFiles(new FilenameFilter() {
                @Override
                public boolean accept(File dir, String name) {
                    return name.startsWith(category + "-");
                }
            });
            if (contourFiles.length == 1) {
                continue;
            }
            int qIdx = rand.nextInt(contourFiles.length);

            for (int j = 0; j < contourFiles.length; ++j) {
                File file = contourFiles[j];
                String[] lines = Files.lines(file.toPath()).toArray(String[]::new);
                PairFloatArray p = new PairFloatArray(lines.length);
                for (String line : lines) {
                    String[] items = line.split(",");
                    p.add(Float.valueOf(items[0]), Float.valueOf(items[1]));
                }
                if (j == qIdx) {
                    queryFileNameMap.put(queries.size(), file.getName());
                    queries.add(p);
                } else {
                    indexFileNameMap.put(targets.size(), file.getName());
                    targets.add(p);
                }
            }
        }

        int topK = 10;
        int n = 100;
        int minBlockSize = (int)Math.round(0.25*n);
        System.out.printf("topK=%d, curveDim=%d, minBlockSize=%d\n", topK, n, minBlockSize);

        long start = System.nanoTime();
        MultiPartialShapeMatcher m = new MultiPartialShapeMatcher(n, minBlockSize, targets);
        long stop = System.nanoTime();
        double diffSec = (stop - start)/1.E9;
        System.out.printf("build index took %f sec (%d ns)\n", diffSec, stop-start);
        System.out.flush();

        float tol = 0.0000001f;

        for (int i = 1; i < classes.length; ++i) {
            String category = classes[i];
            System.out.printf("i=%d, category=%s\n", i, category);
            PairFloatArray q = queries.get(i);
            MultiPartialShapeMatcher.Results results = m.query(q, topK);
            assertEquals(topK, results.getDistances().size());
            for (int j = 0; j < results.getOriginalCurveIndexes().size(); ++j) {
                float dist = results.getDistances().get(j);
                int origIdx = results.getOriginalCurveIndexes().get(j);
                int tOffIdx = results.getOffsetsTargets().get(j);
                int qOffIdx = results.getOffsetsQuery().get(j);
                int len = results.getMatchingLengths().get(j);
                String targetFileName = indexFileNameMap.get(origIdx);
                //System.out.printf("  k=%d, res=%s (q=%s), len=%d\n",
                //        j, targetFileName, queryFileNameMap.get(i), len);
                PairFloatArray t, t2;
                t = results.getDBCurves().get(j);
                t2 = targets.get(origIdx);
                assertEquals(t, t2, tol);
                /*if (!(targetFileName.startsWith(category + "-"))) {
                    System.out.printf("ERROR: follow up: i=%d, targetOrigIdx=%d, qOffIdx=%d (qS=(%f,%f)), tOffIdx=%d (pS=(%f,%f)), len=%d\n",
                            i, origIdx, qOffIdx, q.getX(qOffIdx), q.getY(qOffIdx),
                            tOffIdx, t.getX(tOffIdx), t.getY(tOffIdx), len);
                    System.out.printf("plotting query %s\n", plot(q, i));
                    System.out.printf("plotting target %s\n", plot(t, j*targets.size()));
                    continue;
                    //System.out.printf("plotting target false positive %s\n", plot(tRes, j + queries.size()));
                }
                commenting out for now because class shape does match other class members sometimes.
                assertTrue(targetFileName.startsWith(category + "-"));

                 */
            }
        }
    }

    private void assertEquals(PairFloatArray t, PairFloatArray t2, float tol) {
        assertEquals(t.getN(), t2.getN());
        for (int i = 0; i < t.getN(); ++i) {
            assertTrue(Math.abs(t.getX(i) - t2.getX(i)) < tol);
            assertTrue(Math.abs(t.getY(i) - t2.getY(i)) < tol);
        }
    }

    public void testMPEG7_1() throws IOException {
        /** test using the MPEG-7 shape dataset

         https://dabi.temple.edu/external/shape/MPEG7/dataset.html

         each class (category) has several contours of the same object with small differences in articulation or occlusion,
         and more.

         randomly extract one contour from each class, then put all other contours in the database.
         the single contours are the queries.

         test the method for topK by unique contours.
         */

        //TODO: implement method for keeping topK matches, but unique contours among those

    }

    public void testSimple0() throws Exception {

        /* shape 1:

        7                    *
        6                 *     *
        5              *           *
        4           *                 *
        3        *                       *
        2     *  *  *  *  *  *  *  *  *  *  *
        1
        0
           0  1  2  3  4  5  6  7  8  9 10 11


        shape 2:
        7                    *
        6                 *  *
        5              *     *  *  *
        4           *              *
        3        *                 *  *  *
        2     *  *  *  *  *  *  *  *  *  *  *
        1
        0
           0  1  2  3  4  5  6  7  8  9 10 11

        shape 3:
        7
        6                 *  *  *
        5               *          *
        4               *          *
        3               *          *
        2                 *  *  *
        1
        0
           0  1  2  3  4  5  6  7  8  9 10 11
        */

        PairFloatArray shape1 = PartialShapeMatcherTest.getTriangle();
        PairFloatArray shape2 = PartialShapeMatcherTest.getShape2();
        PairFloatArray shape3 = PartialShapeMatcherTest.getShape3();

        List<PairFloatArray> curves = new ArrayList<>();
        curves.add(shape1);
        curves.add(shape2);
        curves.add(shape3);

        //System.out.printf("plotting %s\n", plot(shape1, 1));
        //System.out.printf("plotting %s\n", plot(shape2, 2));
        //System.out.printf("plotting %s\n", plot(shape3, 3));

        int n = shape1.getN();
        MultiPartialShapeMatcher m = new MultiPartialShapeMatcher(n, 3, curves);

        int topK = 10;
        MultiPartialShapeMatcher.Results results = m.query(shape1, topK);
        assertEquals(topK, results.distances.size());
        assertEquals(0, results.getOriginalCurveIndexes().get(0).intValue());
        assertEquals(n - 1, results.matchingLengths.get(0).intValue());
        assertEquals(0, results.offsetsTargets.get(0).intValue());
        assertTrue(Math.abs(results.getDistances().get(0).floatValue()) < 1E-6);
        assertEquals(0, results.getOffsetsQuery().get(0).intValue());

        /*
        results:
            curves = new ArrayList<>();
            offsetsTargets = new ArrayList<>();
            offsetsQuery = new ArrayList<>();
            matchingLengths = new ArrayList<>();
            distances = new ArrayList<>();
         */
    }

    public void testSimple1() throws Exception {

        /* shape 1:

        7                    *
        6                 *     *
        5              *           *
        4           *                 *
        3        *                       *
        2     *  *  *  *  *  *  *  *  *  *  *
        1
        0
           0  1  2  3  4  5  6  7  8  9 10 11


        shape 2:
        7                    *
        6                 *  *
        5              *     *  *  *
        4           *              *
        3        *                 *  *  *
        2     *  *  *  *  *  *  *  *  *  *  *
        1
        0
           0  1  2  3  4  5  6  7  8  9 10 11

        shape 3:
        7
        6                 *  *  *
        5               *          *
        4               *          *
        3               *          *
        2                 *  *  *
        1
        0
           0  1  2  3  4  5  6  7  8  9 10 11
        */

        PairFloatArray shape1 = PartialShapeMatcherTest.getTriangle();
        PairFloatArray shape2 = PartialShapeMatcherTest.getShape2();
        PairFloatArray shape3 = PartialShapeMatcherTest.getShape3();

        List<PairFloatArray> curves = new ArrayList<>();
        curves.add(shape3);
        curves.add(shape2);

        //System.out.printf("plotting %s\n", plot(shape1, 1));
        //System.out.printf("plotting %s\n", plot(shape2, 2));
        //System.out.printf("plotting %s\n", plot(shape3, 3));

        int n = shape1.getN();
        int minBlockSize = n/2;
        MultiPartialShapeMatcher m = new MultiPartialShapeMatcher(n, minBlockSize, curves);

        int topK = 10;
        MultiPartialShapeMatcher.Results results = m.query(shape1, topK);
        assertEquals(topK, results.distances.size());
        assertEquals(1, results.getOriginalCurveIndexes().get(0).intValue());
        assertTrue(results.offsetsTargets.get(0).intValue() >= 12 && results.offsetsTargets.get(0).intValue() <= 17);
        assertTrue( results.getOffsetsQuery().get(0).intValue() >= 10
                && results.getOffsetsQuery().get(0).intValue() <= 14);

    }

    public void testAndroidStatues_SAM2_contours() throws IOException {

        PairFloatArray queryCurve = getAndroidQueryCurve();

        // reducing curveDimension to 50, put android_statues_01 into top10, but not top1
        int curveDimension = queryCurve.getN();
        int minDim = (int)Math.round(0.2*curveDimension);

        String dirPath = ResourceFinder.findTestResourcesDirectory();
        String contourDirPath = dirPath + sep + "android_contours";
        for (int andI = 1; andI <= 4; ++andI) {
            String filenameRoot = "android_statues_0" + andI;

            String filePrefix = filenameRoot + "_contour";
            List<PairFloatArray> curves = readCurves(contourDirPath, filePrefix);

            ImageExt img = ImageIOHelper.readImageExt(dirPath + sep + filenameRoot + ".jpg");
            plot(img.copyToGreyscale().copyToColorGreyscaleExt(), curves,
                    filenameRoot + "_closed_curves_sam2_", 3);

            MultiPartialShapeMatcher matcher = new MultiPartialShapeMatcher(curveDimension, minDim, curves);

            int topK = 10;
            MultiPartialShapeMatcher.Results results = matcher.query(queryCurve, topK);

            List<PairFloatArray> top = new ArrayList<>();
            top.add(results.getDBCurves().getFirst());
            plot(img.copyToGreyscale().copyToColorGreyscaleExt(), top,
                    filenameRoot + "_found_top_sam2_", 1);
            plot(img.copyToGreyscale().copyToColorGreyscaleExt(), results.getDBCurves(),
                    filenameRoot + "_found_sam2_", 1);

        }

    }

    protected PairFloatArray getAndroidQueryCurve() throws IOException {
        String fileName0 = "android_statues_03_sz1_mask_small.png";
        int idx0 = fileName0.lastIndexOf(".");
        String fileName0Root = fileName0.substring(0, idx0);
        String filePath0 = ResourceFinder.findFileInTestResources(fileName0);
        GreyscaleImage img0 = ImageIOHelper.readImageAsBinary(filePath0);
        int[] labels0 = new int[img0.getNPixels()];
        for (int i = 0; i < img0.getNPixels(); ++i) {
            if (img0.getValue(i) > 0) {
                labels0[i] = 1;
            }
        }
        // template n = 174
        Map<Integer, PairIntArray> shapes0 = PerimeterFinder3.extractOrderedBorders(labels0, img0.getWidth(), img0.getHeight());
        PairFloatArray queryCurve = MultiPartialShapeMatcher.convert(shapes0.get(1));
        //plot(img0.copyToColorGreyscaleExt(), shapes0.get(1), 0,"_template_");
        return queryCurve;
    }

    public void testAndroidStatues() throws IOException {
        // this uses quick segmentation made from SLIC Super Pixels with polar angle theta of U and V in CIE LUV color space

        if (false) {
            calcAndWriteCurvesToFile();
        }

        PairFloatArray queryCurve = getAndroidQueryCurve();

        // reducing curveDimension to 50, put android_statues_01 into top10, but not top1
        int curveDimension = queryCurve.getN();
        int minDim = (int)Math.round(0.2*curveDimension);

        //PairFloatArray p = extractOrderedBoundary(img0);

        // create for each image, the db of closed curve shapes to search.

        // for this purpose, will use SLICSuperPixels to make the segmentation hence labeled pixels, though it produces
        // over segmented regions for the parameters I use.
        // then will use PerimeterFinder3 to extract the shapes.
        // then, we have curves to place in db (as outlined in MultiPartialShapeMatcher).

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

            int idx = fileName1.lastIndexOf(".");
            String fileName1Root = fileName1.substring(0, idx);
            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            ImageExt img = ImageIOHelper.readImageExt(filePath1);

            List<PairIntArray> curves = readInCurves(fileName1Root);
            System.out.printf("read %d curves\n", curves.size());

            List<PairFloatArray> dbCurves = MultiPartialShapeMatcher.convert(curves);

            /*
            int _n1 = curves.size()/2;
            int _n2 = curves.size();
            plot(img.copyToGreyscale().copyToColorGreyscaleExt(), dbCurves.subList(0, _n1),
                    fileName1Root + "_closed_curves_0");
            plot(img.copyToGreyscale().copyToColorGreyscaleExt(), dbCurves.subList(_n1, _n2),
                    fileName1Root + "_closed_curves_1");
            */
            plot(img.copyToGreyscale().copyToColorGreyscaleExt(), dbCurves,
                    fileName1Root + "_closed_curves_", 3);

            MultiPartialShapeMatcher matcher = new MultiPartialShapeMatcher(curveDimension, minDim, dbCurves);

            int topK = 10;

            MultiPartialShapeMatcher.Results results = matcher.query(queryCurve, topK);

            List<PairFloatArray> top = new ArrayList<>();
            top.add(results.getDBCurves().getFirst());
            plot(img.copyToGreyscale().copyToColorGreyscaleExt(), top,
                    fileName1Root + "_found_top_", 1);
            plot(img.copyToGreyscale().copyToColorGreyscaleExt(), results.getDBCurves(),
                    fileName1Root + "_found_", 1);

            //write_centroids(labels, img, fileName1Root);

            //extractShapes(img, fileName1Root);

            // use MultiPartialShapeMatcher to match these
        }

    }

    private List<PairIntArray> readInCurves(String fileNameRoot) throws IOException {
        List<PairIntArray> out = new ArrayList<>();

        String path = getCurveFilePath(fileNameRoot);
        FileReader rw = null;
        BufferedReader reader = null;
        try {
            rw = new FileReader(path);
            reader = new BufferedReader(rw);
            String line = reader.readLine();
            while (line != null) {
                String[] items = line.split(",");
                int n = items.length/2;
                PairIntArray p = new PairIntArray(n);
                for (int i = 0; i < n; i+=2) {
                    p.add(Integer.valueOf(items[2*i]), Integer.valueOf(items[2*i+1]));
                }
                out.add(p);
                line = reader.readLine();
            }
        } catch(IOException ex) {
            System.out.printf("ERROR: %s\n", ex.getMessage());
        } finally{
            if (reader != null) {
                reader.close();
            }
            if (rw != null) {
                rw.close();
            }
        }
        return out;
    }

    private void writeOutCurves(Map<Integer, PairIntArray> curves, String fileNameRoot, int minSz, int maxSz) throws IOException {
        String path = getCurveFilePath(fileNameRoot);
        FileWriter fw = null;
        BufferedWriter writer = null;
        try {
            fw = new FileWriter(path);
            writer = new BufferedWriter(fw);
            for (Map.Entry<Integer, PairIntArray> entry : curves.entrySet()) {
                PairIntArray s = entry.getValue();
                if (s.getN() < minSz || s.getN() > maxSz) {
                    continue;
                }
                StringBuilder sb = new StringBuilder();
                for (int k = 0; k < s.getN(); ++k) {
                    if (k > 0) {
                        sb.append(",");
                    }
                    sb.append(s.getX(k)).append(",").append(s.getY(k));
                }
                sb.append(eol);
                writer.write(sb.toString());
            }
            writer.flush();
        } catch(IOException ex) {
            System.out.printf("ERROR: %s\n", ex.getMessage());
        } finally{
            if (writer != null) {
                writer.close();
            }
            if (fw != null) {
                fw.close();
            }
        }
    }

    private String getCurveFilePath(String fileNameRoot) throws IOException {
        String dir = ResourceFinder.findTestResourcesDirectory();
        String path = dir + sep + "test_shapes" + sep + fileNameRoot + "_curves_.dat";
        return path;
    }

    private void calcAndWriteCurvesToFile() throws IOException {

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

            /*
            (1) closed curves from  over-segmented images
            (2) closed curves from under-segmented images that preserve a couple of larger compositions
             */

            //int minNumPts = p.getN();
            int minNumPts = 10;

            int idx = fileName1.lastIndexOf(".");
            String fileName1Root = fileName1.substring(0, idx);
            String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
            ImageExt img = ImageIOHelper.readImageExt(filePath1);

            ImageProcessor imgProc = new ImageProcessor();
            ImageSegmentation imgSeg = new ImageSegmentation();
            ImageExt img2 = img.copyToImageExt();
            for (int k = 0; k < 0; ++k) {
                imgProc.blur(img2, SIGMA.getValue(SIGMA.TWO));
            }
            GreyscaleImage ptImg = imgProc.createCIELUVTheta_WideRangeLightness(img2, 255);
            String fn = ResourceFinder.findOutputTestDirectory() + "/../" + fileName1Root + "_polar_theta_0_.png";
            ImageIOHelper.writeOutputImage(fn, ptImg);
            EdgeFilterProducts edgeProducts = imgSeg.createGradient(ptImg, 1, System.currentTimeMillis());
            fn = ResourceFinder.findOutputTestDirectory() + "/../" + fileName1Root + "_polar_theta_0_grad.png";
            ImageIOHelper.writeOutputImage(fn, edgeProducts.getGradientXY());

            System.out.printf("i=%d, min=%d, max=%d\n", i, ptImg.min(), ptImg.max());

            // over-segmented:
            int[] labels1 = mergedSLIC(img, fileName1Root, edgeProducts.getGradientXY(), ptImg);

            // under-segmented:
            //int[] labels2 = mergedSLIC(ptImg, fileName1Root, edgeProducts.getGradientXY());

            int[] labels3 = mergeUsingPT(ptImg, fileName1Root, labels1);

            int n = img.getNPixels();
            int minSz = 150;
            int maxSz = Math.round(0.7f * n);

            Map<Integer, PairIntArray> shapes1 = PerimeterFinder3.extractOrderedBorders(labels1,
                    img.getWidth(), img.getHeight(), minSz, maxSz);

            //Map<Integer, PairIntArray> shapes2 = PerimeterFinder3.extractOrderedBorders(labels2,
            //        img.getWidth(), img.getHeight(), minSz, maxSz);

            Map<Integer, PairIntArray> shapes3 = PerimeterFinder3.extractOrderedBorders(labels3,
                    img.getWidth(), img.getHeight(), minSz, maxSz);

            System.out.printf("i=%d, n_curves=%d, %d\n", i, shapes1.size(), shapes3.size());

            List<PairFloatArray> list1 = new ArrayList<>();
            for (Map.Entry<Integer, PairIntArray> entry : shapes1.entrySet()) {
                PairIntArray s = entry.getValue();
                if (s.getN() < minSz || s.getN() > maxSz) {
                    continue;
                }
                PairFloatArray f = new PairFloatArray(s.getN());
                for (int j = 0; j < s.getN(); ++j) {
                    f.add(s.getX(j), s.getY(j));
                }
                list1.add(f);
            }
            plot(img.copyToGreyscale().copyToColorGreyscaleExt(), list1, fileName1Root + "_closed_curves_under", 1);

            /*List<PairFloatArray> list2 = new ArrayList<>();
            for (Map.Entry<Integer, PairIntArray> entry : shapes2.entrySet()) {
                PairIntArray s = entry.getValue();
                PairFloatArray f = new PairFloatArray(s.getN());
                for (int j = 0; j < s.getN(); ++j) {
                    f.add(s.getX(j), s.getY(j));
                }
                list2.add(f);
            }
            plot(img.copyToGreyscale().copyToColorGreyscaleExt(), list2, fileName1Root + "_closed_curves_over");
            */

            List<PairFloatArray> list3 = new ArrayList<>();
            for (Map.Entry<Integer, PairIntArray> entry : shapes3.entrySet()) {
                PairIntArray s = entry.getValue();
                if (s.getN() < minSz || s.getN() > maxSz) {
                    continue;
                }
                PairFloatArray f = new PairFloatArray(s.getN());
                for (int j = 0; j < s.getN(); ++j) {
                    f.add(s.getX(j), s.getY(j));
                }
                list3.add(f);
            }
            plot(img.copyToGreyscale().copyToColorGreyscaleExt(), list3, fileName1Root + "_closed_curves_under_over", 1);

            writeOutCurves(shapes3, fileName1Root, minSz, maxSz);
        }
    }

    private int[] mergedSLIC(ImageExt img, String fileName1Root, GreyscaleImage gradImg, GreyscaleImage luvPolarThetaImg) throws IOException {

        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(img, SIGMA.getValue(SIGMA.TWO));

        METHOD method = METHOD.MIN_GRADIENT;
        //int nc = 100; double thresh = 2.5;
        //int nc = 70; double thresh = 2.5;
        //int nc = 50; double thresh = 3.0;//2nd best
        int nc = 40; double thresh = 1.0; // best results, but over-segmented

        //METHOD method = METHOD.MODE;
        //int nc = 100; double thresh = 5.0;
        //int nc = 70; double thresh = 5.0;

        int[] labels2 = slic(img.copyToImageExt(), fileName1Root, nc, gradImg);

        //int nLabels2 = MergeLabels.mergeWithMinGrad(img, labels2, thresh, gradImg);
        //System.out.printf("merged k=%d, nc=%d\n", nLabels2, nc);
        //int thresh2 = 2;
        //int nLabels3 = MergeLabels.mergeCIELUVPolarTheta(luvPolarThetaImg, labels2, thresh2);
        //System.out.printf("merged k=%d\n", nLabels3);
        //double thresh3 = 0.985;
        //int nLabels4 = MergeLabels.mergeCIELUVPolarThetaH(luvPolarThetaImg, labels2, thresh3);
        //System.out.printf("merged k=%d\n", nLabels4);

        ImageExt im3 = img.copyToImageExt();
        ImageIOHelper.addAlternatingColorLabelsToRegion(im3, labels2);
        MiscDebug.writeImage(im3, "_" + fileName1Root + "_slic_merged_");

        //return new int[][]{labels, labels2};
        return labels2;
    }

    private int[] mergeUsingPT(GreyscaleImage img, String fileName1Root, int[] labels) {
        ImageProcessor imageProcessor = new ImageProcessor();
        //imageProcessor.blur(img, SIGMA.getValue(SIGMA.TWO));

        labels = Arrays.copyOf(labels, labels.length);
        /*
        int thresh2 = 10;
        int nLabels3 = MergeLabels.mergeCIELUVPolarTheta(img, labels, thresh2);
        System.out.printf("merged k=%d\n", nLabels3);
        */
        ///*
        double thresh3 = 0.65;
        int nLabels4 = MergeLabels.mergeCIELUVPolarThetaH(img, labels, thresh3);
        System.out.printf("merged k=%d\n", nLabels4);
        //*/

        /*
        int thresh5 = 10;
        int nLabels5 = MergeLabels.mergeCIELUVPolarTheta(img, labels, thresh5);
        System.out.printf("merged k=%d\n", nLabels5);
        */

        ImageExt im3 = img.copyToColorGreyscaleExt();
        ImageIOHelper.addAlternatingColorLabelsToRegion(im3, labels);
        MiscDebug.writeImage(im3, "_" + fileName1Root + "_slic_merged_");

        //return new int[][]{labels, labels2};
        return labels;
    }

    private int[] mergedSLIC(GreyscaleImage img, String fileName1Root, GreyscaleImage gradImg) throws IOException {

        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(img, SIGMA.getValue(SIGMA.TWO));

        METHOD method = METHOD.MIN_GRADIENT;
        //int nc = 100; double thresh = 2.5;
        //int nc = 70; double thresh = 2.5;
        //int nc = 50; double thresh = 3.0;//2nd best
        int nc = 40;
        nc = 20;
        //nc = 30; // best
        //nc = 20;
        //nc = 10;
        //nc = 20; thresh = 0.;
        //METHOD method = METHOD.MODE;
        //int nc = 100; double thresh = 5.0;
        //int nc = 70; double thresh = 5.0;

        int[] labels2 = slic(img, fileName1Root, nc, gradImg);

        ///*
        int thresh2 = 10;
        int nLabels3 = MergeLabels.mergeCIELUVPolarTheta(img, labels2, thresh2);
        System.out.printf("merged k=%d\n", nLabels3);
        //*/
        ///*
        double thresh3 = 0.75;
        int nLabels4 = MergeLabels.mergeCIELUVPolarThetaH(img, labels2, thresh3);
        System.out.printf("merged k=%d\n", nLabels4);
        //*/

        /// *
        int thresh5 = 10;
        int nLabels5 = MergeLabels.mergeCIELUVPolarTheta(img, labels2, thresh5);
        System.out.printf("merged k=%d\n", nLabels5);
        //*/

        ImageExt im3 = img.copyToColorGreyscaleExt();
        ImageIOHelper.addAlternatingColorLabelsToRegion(im3, labels2);
        MiscDebug.writeImage(im3, "_" + fileName1Root + "_slic_merged_");

        //return new int[][]{labels, labels2};
        return labels2;
    }

    public void extractShapes(ImageExt img, String fileName1Root) throws IOException {

        ImageProcessor imageProcessor = new ImageProcessor();
        //imageProcessor.blur(img, SIGMA.getValue(SIGMA.ONE));

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
                      String fileSuffix, int nDot) throws IOException {

        Image im = img.copyToGreyscale().copyToColorGreyscale();
        int[] clr;
        for (int ii = 0; ii < shapes.size(); ++ii) {
            clr = ImageIOHelper.getNextRGB(ii);
            ImageIOHelper.addCurveToImage(shapes.get(ii), im, nDot, clr[0],
                    clr[1], clr[2]);
        }
        MiscDebug.writeImage(im, fileSuffix);

    }

    private void plot(ImageExt img, Set<Integer> points,
        String fileSuffix) throws IOException {
        Image im = img.copyToGreyscale().copyToColorGreyscale();
        int[] clr;
        int nDot = 0;
        clr = ImageIOHelper.getNextRGB(0);
        for (int idx : points) {
            ImageIOHelper.addPointToImage(idx % im.getWidth(),
               idx / im.getWidth(), im, nDot, clr[0], clr[1], clr[2]);
        }
        MiscDebug.writeImage(im, fileSuffix);
    }
    private void plot(ImageExt img, PairIntArray pts, int nDot,
                      String fileSuffix) throws IOException {
        Image im = img.copyToGreyscale().copyToColorGreyscale();
        int[] clr;
        clr = ImageIOHelper.getNextRGB(0);
        for (int i = 0; i < pts.getN(); ++i) {
            ImageIOHelper.addPointToImage(pts.getX(i),
                    pts.getY(i), im, nDot, clr[0], clr[1], clr[2]);
        }
        MiscDebug.writeImage(im, fileSuffix);
    }

    private int[] slic(ImageExt img, String fileName1Root, int nc, GreyscaleImage gradImg) throws IOException {

        ImageSegmentation imS = new ImageSegmentation();

        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(img, SIGMA.getValue(SIGMA.TWO));

        int method = 0;

        SLICSuperPixels slic = new SLICSuperPixels(img, nc);
        slic.setGradient(gradImg);
        slic.calculate();
        int[] labels = slic.getLabels();

        LabelHelper.resolveByConnectedness(labels, img.getWidth(), img.getHeight(), true);

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

    private int[] slic(GreyscaleImage img, String fileName1Root, int nc, GreyscaleImage gradImg) throws IOException {

        ImageSegmentation imS = new ImageSegmentation();

        ImageProcessor imageProcessor = new ImageProcessor();
        //imageProcessor.blur(img, SIGMA.getValue(SIGMA.TWO));

        int method = 0;

        SLICSuperPixelsGS slic = new SLICSuperPixelsGS(img, nc);
        slic.setGradient(gradImg);
        slic.calculate();
        int[] labels = slic.getLabels();

        LabelHelper.resolveByConnectedness(labels, img.getWidth(), img.getHeight(), true);

        Image img3 = img.copyToColorGreyscale();
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


    private String plot(PairFloatArray p, int fn) throws Exception {

        float[] x = Arrays.copyOf(p.getX(), p.getN());
        float[] y = Arrays.copyOf(p.getY(), p.getN());
        float xMax = MiscMath.findMax(x) + 1;
        float yMax = MiscMath.findMax(y) + 1;

        PolygonAndPointPlotter plot = new PolygonAndPointPlotter();

        plot.addPlot(0, xMax, 0, yMax,
                x, y, x, y, "");

        return plot.writeFile(fn);
    }

    private List<PairFloatArray> readCurves(String contourDirPath, String filePrefix) throws IOException {

        List<PairFloatArray> curves = new ArrayList<>();

        String globPattern = filePrefix + "_*.csv";
        try (DirectoryStream<Path> stream = Files.newDirectoryStream(Paths.get(contourDirPath), globPattern)) {
            for (Path entry : stream) {
                if (Files.isRegularFile(entry)) { // Ensure it's a file, not a directory
                    PairFloatArray curve = new PairFloatArray();
                    try (Stream<String> lines = Files.lines(entry)) {
                        lines.map(line -> line.split(","))
                        .forEach(parts -> {
                            curve.add(Float.valueOf(parts[0]), Float.valueOf(parts[1]));
                        });
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    curves.add(curve);
                }
            }
        }

        return curves;
    }

}
