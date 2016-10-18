package algorithms.imageProcessing.features;

import algorithms.compGeometry.clustering.KMeansHSV;
import algorithms.imageProcessing.ColorHistogram;
import algorithms.imageProcessing.FixedSizeSortedVector;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.GroupPixelHSV;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.MedianTransform;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.PeriodicFFT;
import algorithms.imageProcessing.SegmentationMergeThreshold;
import algorithms.imageProcessing.SummedAreaTable;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.Complex;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import com.climbwithyourfeet.clustering.DTClusterFinder;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PhaseCongruencyDetectorTest extends TestCase {

    public PhaseCongruencyDetectorTest() {
    }

    public void est0() throws Exception {

        String[] fileNames = new String[]{
           // "blox.gif", "lab.gif", "house.gif", "seattle.jpg", "merton_college_I_001.jpg",
           // "susan-in_plus.png", "lena.jpg",
           // "campus_010.jpg",
            "android_statues_01.jpg",
            "android_statues_02.jpg", "android_statues_03.jpg", "android_statues_04.jpg"
        };

        ImageProcessor imageProcessor = new ImageProcessor();

        for (String fileName : fileNames) {

            System.out.println("fileName=" + fileName);

            String filePath = ResourceFinder.findFileInTestResources(fileName);

            GreyscaleImage img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();

            PhaseCongruencyDetector phaseCDetector = new PhaseCongruencyDetector();
            phaseCDetector.setToCreateCorners();
            PhaseCongruencyDetector.PhaseCongruencyProducts products =
                phaseCDetector.phaseCongMono(img);

            assertNotNull(products);
            int[][] thinned = products.getThinned();

            GreyscaleImage pcImg = img.createWithDimensions();
            GreyscaleImage out2 = img.createWithDimensions();
            GreyscaleImage out = img.createWithDimensions();
            for (int i = 0; i < out.getWidth(); ++i) {
                for (int j = 0; j < out.getHeight(); ++j) {
                    int vPC = (int)Math.round(255. * products.getPhaseCongruency()[j][i]);
                    if (thinned[j][i] > 0) {
                        out.setValue(i, j, thinned[j][i]);
                        out2.setValue(i, j, vPC);
                    }
                    pcImg.setValue(i, j, vPC);
                }
            }
            MiscDebug.writeImage(out, "_thinned_" + fileName + "_");
            MiscDebug.writeImage(out2, "_pc_thinned_" + fileName + "_");
            MiscDebug.writeImage(pcImg, "_pc_" + fileName + "_");

            // ----- make O1 edges
            ImageExt imgClr = ImageIOHelper.readImageExt(filePath);

            GreyscaleImage o1 = imageProcessor.createO1(imgClr);

            products =
                phaseCDetector.phaseCongMono(o1);
            thinned = products.getThinned();
            out = img.createWithDimensions();
            pcImg = img.createWithDimensions();
            out2 = img.createWithDimensions();
            for (int i = 0; i < out.getWidth(); ++i) {
                for (int j = 0; j < out.getHeight(); ++j) {
                    int vPC = (int)Math.round(255. * products.getPhaseCongruency()[j][i]);
                    if (thinned[j][i] > 0) {
                        out.setValue(i, j, thinned[j][i]);
                        out2.setValue(i, j, vPC);
                    }
                    pcImg.setValue(i, j, vPC);
                }
            }

            MiscDebug.writeImage(out, "_thinned_o1_" + fileName);

            MiscDebug.writeImage(out2, "_pc_thinned_o1_" + fileName);

            MiscDebug.writeImage(pcImg, "_pc_o1_" + fileName);
        }
    }

    public void est1() throws Exception {

        String[] fileNames = new String[]{
            //"blox.gif",
            //"lab.gif",
            "house.gif",
            //"seattle.jpg",
            //"merton_college_I_001.jpg",
            // "susan-in_plus.png",
            //"lena.jpg",
            //"campus_010.jpg",
            //"android_statues_01.jpg",
            //"android_statues_02.jpg",
            //"android_statues_03.jpg",
            //"android_statues_04.jpg"
        };

        ImageSegmentation imageSegmentation = new ImageSegmentation();

        for (String fileName : fileNames) {

            System.out.println("fileName=" + fileName);

            String filePath = ResourceFinder.findFileInTestResources(fileName);

            ImageExt img = ImageIOHelper.readImageExt(filePath);

            GreyscaleImage edgeImage = imageSegmentation.createColorEdges(img);
        }
    }

    public void est2() throws Exception {

        String[] fileNames = new String[]{
            //"seattle.jpg",
            "merton_college_I_001.jpg",
            //"merton_college_I_002.jpg",
            //"lena.jpg",
            //"campus_010.jpg",
            //"android_statues_01.jpg",
            //"android_statues_02.jpg",
            //"android_statues_03.jpg",
            //"android_statues_04.jpg"
        };

        ImageSegmentation imageSegmentation = new ImageSegmentation();

        boolean doDecimate = true;
        int minDimension = 300;//512;//300;

        for (String fileName : fileNames) {
            try {
            System.out.println("fileName=" + fileName);

            String filePath = ResourceFinder.findFileInTestResources(fileName);

            ImageExt img = ImageIOHelper.readImageExt(filePath);
            List<ImageExt> transformed = null;

            int selectIdx = -1;

            if (doDecimate) {
                MedianTransform mt = new MedianTransform();
                transformed = new ArrayList<ImageExt>();
                //List<ImageExt> coeffs = new ArrayList<ImageExt>();
                //mt.multiscalePyramidalMedianTransform(img, transformed, coeffs);
                mt.<ImageExt>multiscalePyramidalMedianTransform2(img, transformed);
                //GreyscaleImage r = mt.reconstructPyramidalMultiscaleMedianTransform(
                //    transformed.get(transformed.size() - 1), coeffs);

                // choose the first image which is smaller than 300 x 300
                for (int j = 0; j < transformed.size(); ++j) {
                    ImageExt tr = transformed.get(j);
                    //MiscDebug.writeImage(tr, "_tr_" + j);
                    if (selectIdx == -1) {
                        if (tr.getWidth() <= minDimension && tr.getHeight() <= minDimension) {
                            selectIdx = j;
                            img = transformed.get(selectIdx);
                        }
                    }
                }
            }

            GreyscaleImage edgeImage = imageSegmentation.createColorEdges(img);

            edgeImage = imageSegmentation.fillInGapsOf1(edgeImage,
                new HashSet<PairInt>(), 255);

// TODO: edit performSegmentationWithColorEdges
            List<Set<PairInt>> segmentedPoints =
                imageSegmentation.performSegmentationWithColorEdges(img,
                edgeImage, SegmentationMergeThreshold.DEFAULT, fileName);

            List<PairIntArray> perimeters = BlobsAndPerimeters.extractBoundsOfBlobs(
                segmentedPoints, false, false, 1, img.getWidth(), img.getHeight());

            if (selectIdx > -1) {
                ImageProcessor imageProcessor = new ImageProcessor();
                for (int jj = (selectIdx - 1); jj > -1; --jj) {
                    perimeters = imageProcessor.unbinZeroPointLists(perimeters, 2);
                    //Image outImg = transformed.get(jj);
                    //ImageIOHelper.addAlternatingColorCurvesToImage(perimeters, outImg);
                    //MiscDebug.writeImage(outImg, "_final_edges_" + (jj + 1) + "_" + fileName);
                }
            }

            Image outImg = ImageIOHelper.readImage(filePath);
            ImageIOHelper.addAlternatingColorCurvesToImage(perimeters, outImg, 2);
            MiscDebug.writeImage(outImg, "_final_edges_" + fileName);

            } catch (Throwable t) {
                int z = 1;
            }
        }
    }

    // use of phase congruency on grey, r-g, g-b, and r-b then combining all results
    // the phase ongruency is performed on small overlapping regions.
    public void est3() throws Exception {

        String[] fileNames = new String[]{
            //"seattle.jpg",
            "merton_college_I_001.jpg",
            // "lena.jpg",
            // "campus_010.jpg",
            //"android_statues_01.jpg",
            //"android_statues_02.jpg",
            //"android_statues_03.jpg",
            //"android_statues_04.jpg"
        };

        ImageSegmentation imageSegmentation = new ImageSegmentation();

        ImageProcessor imageProcessor = new ImageProcessor();

        for (String fileName : fileNames) {

            System.out.println("fileName=" + fileName);

            String filePath = ResourceFinder.findFileInTestResources(fileName);

            ImageExt img = ImageIOHelper.readImageExt(filePath);

            GreyscaleImage combined = imageSegmentation.createColorEdges_1(img, 100);

            MiscDebug.writeImage(combined, "_combined_" + fileName);
        }
    }

    public void test4() throws Exception {

        String[] fileNames = new String[]{
            //"seattle.jpg",
            //"merton_college_I_001.jpg",
            //"house.gif",
            //"lena.jpg",
            //"campus_010.jpg",
            //"android_statues_01.jpg",
            "android_statues_02.jpg",
            //"android_statues_03.jpg",
            //"android_statues_04.jpg"
        };

        ImageSegmentation imageSegmentation = new ImageSegmentation();
        /*
        a look at O(N) patterns to make a single combined image for input to
        phase conguency that would result in closed curves for the main objects.
        */

        for (String fileName : fileNames) {

            System.out.println("fileName=" + fileName);

            String filePath = ResourceFinder.findFileInTestResources(fileName);

            ImageExt img = ImageIOHelper.readImageExt(filePath);

            //GreyscaleImage combined = imageSegmentation.createColorEdges_2(img);

            //MiscDebug.writeImage(combined, "_MAX_SOBEL_EDGES_");

            //GreyscaleImage combinedCopy = combined.copyImage();
            //ImageProcessor imageProcessor = new ImageProcessor();
            //imageProcessor.applyAdaptiveMeanThresholding(combinedCopy);
            //MiscDebug.writeImage(combinedCopy, "_MAX_SOBEL_EDGES__AMT_");

            GreyscaleImage img2 = img.copyToGreyscale2();

            PhaseCongruencyDetector phaseCDetector
                = new PhaseCongruencyDetector();
            phaseCDetector.setToExtractNoise();
            phaseCDetector.setToDebug();
            //phaseCDetector.setToCreateCorners();
            PhaseCongruencyDetector.PhaseCongruencyProducts products =
                phaseCDetector.phaseCongMono(img2);

            List<Set<PairInt>> subsetNoise = products.getSubsetNoise();
            ColorHistogram chhist = new ColorHistogram();
            int[][][] chs = new int[subsetNoise.size()][][];
            int[][] sum = chhist.createWithDefaultSize();
            int count = 0;
            for (int i = 0; i < subsetNoise.size(); ++i) {
                Set<PairInt> set = subsetNoise.get(i);
                chs[i] = chhist.histogramHSV(img, set);
                chhist.add2To1(sum, chs[i]);
                count += set.size();
            }
            final int nPoints = count;

            // multiply each histogram by count/set.size to normalize all to
            // same number of counts
            for (int i = 0; i < subsetNoise.size(); ++i) {
                int factor = Math.round((float)count/
                    (float)subsetNoise.get(i).size());
                int[][] ch = chs[i];
                for (int ii = 0; ii < sum.length; ++ii) {
                    for (int jj = 0; jj < sum[ii].length; ++jj) {
                        ch[ii][jj] *= factor;
                    }
                }
            }

            float[] valuesH = new float[subsetNoise.size()];
            float[] valuesS = new float[subsetNoise.size()];
            float[] valuesV = new float[subsetNoise.size()];
            float[] values = new float[subsetNoise.size()];
            float[][] x2 = new float[subsetNoise.size()][];
            // looking at the H, S, and V separately
            for (int i = 0; i < subsetNoise.size(); ++i) {
                x2[i] = new float[3];
                int[][] ch = chs[i];
                assert(ch.length == sum.length);
                chhist.chiSquaredSum(ch, sum, x2[i]);
                valuesH[i] = x2[i][0];
                valuesS[i] = x2[i][1];
                valuesV[i] = x2[i][2];
                values[i] = chhist.intersection(ch, sum);
                count++;
            }
            int nBins = 32;
            HistogramHolder histH = Histogram.createSimpleHistogram(
                nBins, valuesH, Errors.populateYErrorsBySqrt(valuesH));
            histH.plotHistogram("cch x2 H", "_color_hist_H");
            HistogramHolder histS = Histogram.createSimpleHistogram(
                nBins, valuesS, Errors.populateYErrorsBySqrt(valuesS));
            histS.plotHistogram("cch x2 S", "_color_hist_S");
            HistogramHolder histV = Histogram.createSimpleHistogram(
                nBins, valuesV, Errors.populateYErrorsBySqrt(valuesV));
            histV.plotHistogram("cch x2 V", "_color_hist_V");
            HistogramHolder hist = Histogram.createSimpleHistogram(
                nBins, values, Errors.populateYErrorsBySqrt(values));
            hist.plotHistogram("cch intersections", "_color_hist_");

            // peaks in h might be clusters when using nBins=16
            //   -- for each h peak,
            //      want to know whether their colors are similar, but
            //      avoid an n^2 comparison if possible,
            //      so make s and v histograms within that h range
            //      and hsv combined histograms
            //      -- look at whether
            List<Integer> indexes = MiscMath.findStrongPeakIndexes(histH,
                0.02f);
            int nGroups = indexes.size();
            TIntList[] groupIndexes = new TIntList[nGroups];
            for (int i = 0; i < nGroups; ++i) {
                groupIndexes[i] = new TIntArrayList();
            }
            Image dbg = img.createWithDimensions();
            for (int i = 0; i < subsetNoise.size(); ++i) {
                float x2H = x2[i][0];
                float minDiff = Float.MAX_VALUE;
                int minDiffIdx = -1;
                for (int j = 0; j < nGroups; ++j) {
                    float diff = Math.abs(histV.getXHist()[indexes.get(j).intValue()]
                        - x2H);
                    if (diff < minDiff) {
                        minDiff = diff;
                        minDiffIdx = j;
                    }
                }
                groupIndexes[minDiffIdx].add(i);
                int[] clr = ImageIOHelper.getNextRGB(minDiffIdx);
                ImageIOHelper.addCurveToImage(subsetNoise.get(i), dbg,
                    2, clr[0], clr[1], clr[2]);
                System.out.println("groupId=" + minDiffIdx + " "
                    + " r=" + clr[0] + " g=" + clr[1] + " b=" + clr[2]);
            }
            MiscDebug.writeImage(dbg, "_texture_clusters_2_" + fileName);

            List<GreyscaleImage> patterns = new ArrayList<GreyscaleImage>();
            List<GroupPixelHSV> patternColors = new ArrayList<GroupPixelHSV>();
            
            ImageProcessor imageProcessor = new ImageProcessor();
            
            MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

            // ---- a kmeans style iteration to re-org membership
            //      until no change
            int[][][] sumCHs = new int[nGroups][][];
            int nIter = 0;
            final int nIterMax = 10;
            double prevTotSim = Double.MAX_VALUE;
            while (nIter < nIterMax) {
                for (int i = 0; i < nGroups; ++i) {
                    sumCHs[i] = chhist.createWithDefaultSize();
                    TIntList snIndexes = groupIndexes[i];
                    for (int j = 0; j < snIndexes.size(); ++j) {
                        int idx = snIndexes.get(j);
                        chhist.add2To1(sumCHs[i], chs[idx]);
                    }
                    groupIndexes[i].clear();
                }
                double totSim = 0;
                for (int i = 0; i < subsetNoise.size(); ++i) {
                    int[][] ch = chs[i];
                    float maxSim = Float.MIN_VALUE;
                    int maxSimIdx = -1;
                    for (int ii = 0; ii < nGroups; ++ii) {
                        float intersection = chhist.intersection(
                            sumCHs[ii], ch);
                        if (intersection > maxSim) {
                            maxSim = intersection;
                            maxSimIdx = ii;
                        }
                    }
                    groupIndexes[maxSimIdx].add(i);
                    totSim += maxSim;
                }
                System.out.println("nIter=" + nIter + " totSim=" + totSim);
                // sim is 1.0 for perfect match
                if (nIter>0 &&
                    (Math.abs(totSim - prevTotSim) < (0.01 * nPoints))
                    || (totSim > (0.7 * nPoints))) {
                    break;
                }
                prevTotSim = totSim;
                nIter++;
            }

            dbg = img.createWithDimensions();
            for (int i = 0; i < nGroups; ++i) {
                int[] clr = ImageIOHelper.getNextRGB(i);
                TIntList snIndexes = groupIndexes[i];
                for (int j = 0; j < snIndexes.size(); ++j) {
                    int idx = snIndexes.get(j);
                    ImageIOHelper.addCurveToImage(
                        subsetNoise.get(idx), dbg, 2,
                        clr[0], clr[1], clr[2]);
                }
                System.out.println("groupId=" + i + " "
                    + " r=" + clr[0] + " g=" + clr[1] + " b=" + clr[2]);
            }
            MiscDebug.writeImage(dbg, "_texture_clusters_3_" + fileName);

            List<Set<PairInt>> rList = new ArrayList<Set<PairInt>>();
            //TFloatList critDists = new TFloatArrayList();
            TIntList groupIndexList = new TIntArrayList();
            List<GroupPixelHSV> colors = new ArrayList<GroupPixelHSV>();
            
            for (int i = 0; i < nGroups; ++i) {
                TIntList snIndexes = groupIndexes[i];
                if (snIndexes.size() == 0) {
                    continue;
                }
                Set<PairIntWithIndex> points2 = new HashSet<PairIntWithIndex>();
                int maxX = Integer.MIN_VALUE;
                int maxY = Integer.MIN_VALUE;
                for (int j = 0; j < snIndexes.size(); ++j) {
                    Set<PairInt> set = subsetNoise.get(snIndexes.get(j));
                    for (PairInt p : set) {
                        int x = p.getX();
                        int y = p.getY();
                        points2.add(new PairIntWithIndex(x, y, points2.size()));
                        if (x > maxX) {
                            maxX = x;
                        }
                        if (y > maxY) {
                            maxY = y;
                        }
                    }
                }

                DTClusterFinder<PairIntWithIndex> cFinder
                    = new DTClusterFinder<PairIntWithIndex>(points2,
                    maxX + 1, maxY + 1);

                cFinder.setMinimumNumberInCluster(1);
                cFinder.calculateCriticalDensity();
                cFinder.findClusters();
                final int n = cFinder.getNumberOfClusters();

                int maxN = Integer.MIN_VALUE;
                int maxNIdx = -1;
                for (int ii = 0; ii < n; ++ii) {
                    int sz = cFinder.getCluster(ii).size();
                    if (sz > maxN) {
                        maxN = sz;
                        maxNIdx = ii;
                    }
                }
                if (maxNIdx != -1) {

                    Set<PairInt> set = new HashSet<PairInt>();
                    for (PairIntWithIndex p : cFinder.getCluster(maxNIdx)) {
                        set.add(new PairInt(p.getX(), p.getY()));
                    }
                    
                    if (set.isEmpty()) {
                        continue;
                    }
                    
                    // 2 points within this distnce
                    float critDist = 1.f/cFinder.getCriticalDensity();
                    
                    //critDists.add(critDist);
  
                    int patchWidth = Math.round(critDist/2.f);
                    if ((patchWidth & 1) == 0) {
                        patchWidth++;
                    }
                    int patchHalf = patchWidth >> 1;
           
                    // --- make representative patches for this set, where
                    //     patch is the color image histograms centered at
                    //     point and window with width = patchWidth.
                    //     -- would like to pick the most different among them.
                    //        will approx that with sorted by intersection
                    //        with the sum and take the first, last and mid
                    //        set.
                    //        -- then the patches are the window of all pixels
                    //           around the centroid of pick
                    PairInt[] rPoints;
                    if (set.size() < 4) {
                        rPoints = new PairInt[set.size()];
                        int count2 = 0;
                        for (PairInt p : set) {
                            rPoints[count2] = p;
                            count2++;
                        }
                    } else {
                        int[][][] chs2 = new int[set.size()][][];
                        int[][] sum2 = chhist.createWithDefaultSize();
                        int count2 = 0;
                        for (PairInt p : set) {
                            int x = p.getX();
                            int y = p.getY();
                            Set<PairInt> rect = imageProcessor
                                .createWindowOfPoints(x, y, patchHalf,
                                img.getWidth(), img.getHeight());
                            chs2[count2] = chhist.histogramHSV(img, rect);
                            chhist.add2To1(sum2, chs2[count2]);
                            count2++;
                        }
                        FixedSizeSortedVector<IntegerPI> sortedInter = new
                            FixedSizeSortedVector<IntegerPI>(set.size(),
                            IntegerPI.class);
                        count2 = 0;
                        for (PairInt p : set) {
                            float intersection = chhist.intersection(
                                sum2, chs2[count2]);
                            IntegerPI ip = new IntegerPI(
                                Math.round(100 * intersection), p);
                            sortedInter.add(ip);
                            count2++;
                        }
                        rPoints = new PairInt[3];
                        rPoints[0] = sortedInter.getArray()[0].p;
                        rPoints[1] = sortedInter.getArray()[set.size()/2].p;
                        rPoints[2] = sortedInter.getArray()[set.size() - 1].p;
                    }
                   
                    for (PairInt r : rPoints) {
                        
                        int x = r.getX();
                        int y = r.getY();
                        
                        Set<PairInt> rect = imageProcessor.createWindowOfPoints(
                            x, y, patchHalf, img.getWidth(), img.getHeight());
                        
                        // determine avg and std of color to use as a filter for
                        // image. using HSV
                        GroupPixelHSV hsv = new GroupPixelHSV();
                        hsv.calculateColors(rect, img);

                        int[] minMaxXY = MiscMath.findMinMaxXY(rect);
                        int n0 = (minMaxXY[1] - minMaxXY[0]) + 1;
                        int n1 = (minMaxXY[3] - minMaxXY[2]) + 1;
                    
                        GreyscaleImage imagePattern = new GreyscaleImage(
                            patchWidth, patchWidth);
                        
                        for (PairInt p2 : rect) {
                            int x3 = p2.getX();
                            int y3 = p2.getY();
                            int v = img2.getValue(x3, y3);
                            x3 -= minMaxXY[0];
                            y3 -= minMaxXY[2];
                            imagePattern.setValue(x3, y3, v);
                        }
                        
                        patterns.add(imagePattern);
                        rList.add(rect);
                        colors.add(hsv);
                        groupIndexList.add(i);
                    }
                }
            }

            // plot rList
            for (int i = 0; i < rList.size(); ++i) {
                dbg = img.createWithDimensions();
                int[] clr = ImageIOHelper.getNextRGB(i);
                ImageIOHelper.addCurveToImage(rList.get(i), dbg, 2,
                    clr[0], clr[1], clr[2]);
                //System.out.println("**groupId=" + i + " "
                //    + " r=" + clr[0] + " g=" + clr[1] + " b=" + clr[2]);
                MiscDebug.writeImage(dbg, "_texture_clusters_4_" + i + "_" 
                    + fileName);
            }

            /*
            -- color histograms to look for color classes in the subset noise.
            -- look into the texture stats of Malik et al 2001.
               -- whether a sure edge lies along path in between two
                  texture xlusters.
                    along the line between the pixels:
                       W_IC_i_j = 1 - argmax(local maxima of pc perpendicular
                                             to curve)
                    to search for curve points, they use a radii of 30
                    around the  textons of interest.
                  -- this might be useful during some merging steps
                     after super pixels stage.
               -- different goal: adding other terms to the normalized cuts
                  weighting function is in this paper and the DNCuts paper.

                  The term they use for chi squared as the difference between
                  textons and intensities, could be used to make
                  a weighting function for the difference between color
                  histograms:
                   chi squared =
                      1/2 times sum over all bins of : (h1 - h2)^2/(h1 + h2)
                   W_TX_i_j = exp(- chi squared / sigma_TX)
                    (note, for the intensity weight component, sigma is 0.02,
                     and sigma_TX = 0.025,
                     so might expect similar for other weight sigmas.
                -- They also include suggestions for using position, that is
                   adjacency.

            -- explore making a frequency domain filter for spatial domain
               from a representative set of subset noise or centered on it.
               -- might be able to reduce computations due to sparse
                  data
               -- Malik et al. 2001 normalize their texton responses using:
                    F(x) = F(x) X log(1-(|F(x)|/0.03))/|F(x)|
            */

            // make hsv filtered copies of image, then extracy greuscale
            List<GreyscaleImage> filteredHSVImgs = new
                ArrayList<GreyscaleImage>();
            float eFactor = 1.f;
            for (int i = 0; i < colors.size(); ++i) {
                GroupPixelHSV hsv = colors.get(i);
                float h0 = hsv.getAvgH() - eFactor*hsv.getStdDevH();
                float h1 = hsv.getAvgH() + eFactor*hsv.getStdDevH();
                float s0 = hsv.getAvgS() - eFactor*hsv.getStdDevS();
                float s1 = hsv.getAvgS() + eFactor*hsv.getStdDevS();
                float v0 = hsv.getAvgV() - eFactor*hsv.getStdDevV();
                float v1 = hsv.getAvgV() + eFactor*hsv.getStdDevV();
                
                float l0 = hsv.getAvgL() - eFactor*hsv.getStdDevL();
                float l1 = hsv.getAvgL() + eFactor*hsv.getStdDevL();
                float a0 = hsv.getAvgA() - eFactor*hsv.getStdDevA();
                float a1 = hsv.getAvgA() + eFactor*hsv.getStdDevA();
                float b0 = hsv.getAvgB() - eFactor*hsv.getStdDevB();
                float b1 = hsv.getAvgB() + eFactor*hsv.getStdDevB();
                
                Image imgCp = img.copyImage();
                for (int j = 0; j < img.getNPixels(); ++j) {
                    float h = img.getHue(j);
                    float s = img.getSaturation(j);
                    float v = img.getBrightness(j);
                    if (h < h0 || h > h1 || s < s0 || s > s1 || v < v0 || v > v1) {
                        imgCp.setRGB(j, 0, 0, 0);
                        continue;
                    }
                    float[] lab = img.getCIELAB(j);
                    if (lab[0] < l0 || lab[0] > l1 || 
                        lab[1] < a0 || lab[1] > a1 || 
                        lab[2] < b0 || lab[2] > b1) {
                        imgCp.setRGB(j, 0, 0, 0);
                    }
                }
                
                filteredHSVImgs.add(imgCp.copyToGreyscale2());
                MiscDebug.writeImage(imgCp, "_hsv_filtered_" + i);
            }

            GreyscaleImage gsImg = img2;

            // use row major for FFTs
            int nCols = gsImg.getWidth();
            int nRows = gsImg.getHeight();
            PeriodicFFT perfft2 = new PeriodicFFT();
            Complex[][][] perfResults = perfft2.perfft2(gsImg, false);
            Complex[][] fftImage = perfResults[1];

            int buffer = 5;

            List<GreyscaleImage> spatialResponses
                = new ArrayList<GreyscaleImage>(rList.size());

            for (int i = 0; i < rList.size(); ++i) {
               
                Set<PairInt> set = rList.get(i);
                
                int[] minMaxXY = MiscMath.findMinMaxXY(set);
                                
                // imgPattern is in column major format
                GreyscaleImage imgPattern = patterns.get(i);

                Complex[][] fftPattern = PhaseCongruencyDetector
                    .createLowPassFreqDomainFilter(imgPattern);

                // --- image in the frequency domain convolved with texture patch ----
                Complex[][] freqDomainImageTimesPattern =
                    imageProcessor.convolveWithKernel(fftImage, fftPattern);

                // ----- transform that to spatial domain ----
                Complex[][] fComplex = imageProcessor.create2DFFT(
                    freqDomainImageTimesPattern, false, false);

                double[][] transformedReal = new double[nCols][];
                for (int i0 = 0; i0 < nCols; ++i0) {
                    transformedReal[i0] = new double[nRows];
                    for (int i1 = 0; i1 < nRows; ++i1) {
                        transformedReal[i0][i1] = fComplex[i1][i0].abs();
                    }
                }

                MiscMath.applyRescale(transformedReal, 0, 255);
                GreyscaleImage kpFreqR2Img = new GreyscaleImage(nCols, nRows);
                for (int ii = 0; ii < nCols; ++ii) {
                    for (int j = 0; j < nRows; ++j) {
                        kpFreqR2Img.setValue(ii, j,
                            (int)Math.round(transformedReal[ii][j]));
                    }
                }
                MiscDebug.writeImage(kpFreqR2Img, "_freq_spatial_" +
                    i + "_" + fileName);

                GreyscaleImage img3 = filteredHSVImgs.get(i);
                for (int j = 0; j < kpFreqR2Img.getNPixels(); ++j) {
                    if (img3.getValue(j) == 0) {
                        kpFreqR2Img.setValue(j, 0);
                    }
                }
                
                MiscDebug.writeImage(kpFreqR2Img, "_freq_spatial_filtered_" +
                    i + "_" + fileName);
                
                // -- threshold each freq response image
                //     -- the threshold might be dependent upon
                //        the intensity range in the template pattern.
                // -- combine thresholded responses that have the same
                //    groupIndexList value
                // -- store that as spatial responses with key = group index
                
            }

            /*
            assertNotNull(products);
            int[][] thinned = products.getThinned();

            GreyscaleImage pcImg = img2.createWithDimensions();
            GreyscaleImage out2 = img2.createWithDimensions();
            GreyscaleImage out = img2.createWithDimensions();
            for (int i = 0; i < out.getWidth(); ++i) {
                for (int j = 0; j < out.getHeight(); ++j) {
                    int vPC = (int)Math.round(255. * products.getPhaseCongruency()[j][i]);
                    if (thinned[j][i] > 0) {
                        out.setValue(i, j, thinned[j][i]);
                        out2.setValue(i, j, vPC);
                    }
                    pcImg.setValue(i, j, vPC);
                }
            }
            MiscDebug.writeImage(out, "_thinned_" + "_" + fileName + "_");
            MiscDebug.writeImage(out2, "_pc_thinned_"  + "_" + fileName + "_");
            MiscDebug.writeImage(pcImg, "_pc_" + "_" + fileName + "_");


            PhaseCongruencyDetectorPyramidal phaseCDetector0
                = new PhaseCongruencyDetectorPyramidal();
            phaseCDetector0.setToCreateCorners();
            PhaseCongruencyDetectorPyramidal.PhaseCongruencyProducts products0 =
                phaseCDetector0.phaseCongMono(img2);

            assertNotNull(products0);
            thinned = products0.getThinned();
            pcImg = img2.createWithDimensions();
            out2 = img2.createWithDimensions();
            out = img2.createWithDimensions();
            for (int i = 0; i < out.getWidth(); ++i) {
                for (int j = 0; j < out.getHeight(); ++j) {
                    int vPC = (int)Math.round(255.
                        * products0.getPhaseCongruency()[j][i]);
                    if (thinned[j][i] > 0) {
                        out.setValue(i, j, thinned[j][i]);
                        out2.setValue(i, j, vPC);
                    }
                    pcImg.setValue(i, j, vPC);
                }
            }
            MiscDebug.writeImage(out, "_thinned_" + "_" + fileName + "_0");
            MiscDebug.writeImage(out2, "_pc_thinned_"  + "_" + fileName + "_0");
            MiscDebug.writeImage(pcImg, "_pc_" + "_" + fileName + "_0");
            */
        }
    }

    public class IntegerPI implements Comparable<IntegerPI> {

        final Integer d;
        final PairInt p;

        public IntegerPI(Integer d, PairInt p) {
            this.d = d;
            this.p = p;
        }

        @Override
        public int compareTo(IntegerPI other) {
            return d.compareTo(other.d);
        }

    }
}
