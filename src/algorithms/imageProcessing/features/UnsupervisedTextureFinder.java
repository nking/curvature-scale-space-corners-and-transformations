package algorithms.imageProcessing.features;

import algorithms.imageProcessing.ColorHistogram;
import algorithms.imageProcessing.FFTUtil;
import algorithms.FixedSizeSortedVector;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.GroupPixelHSV;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.compGeometry.MiscellaneousCurveHelper;
import algorithms.imageProcessing.PeriodicFFT;
import algorithms.misc.Complex;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
import com.climbwithyourfeet.clustering.ClusterFinder;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.iterator.TLongIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TIntHashSet;
import gnu.trove.set.hash.TLongHashSet;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * class to search for textures within the noise
 * extracted from PhaseCongruencyDetector.java.
 * 
 * @author nichole
 */
public class UnsupervisedTextureFinder {
   
    public TexturePatchesAndResponse[] createTextureImages(ImageExt img,
        PhaseCongruencyDetector.PhaseCongruencyProducts products,
        String debugTag) {
        
        List<TexturePatchesAndResponse> output = new
            ArrayList<TexturePatchesAndResponse>();
        
        GreyscaleImage img2 = img.copyToGreyscale2();
        
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
            int factor = Math.round((float) count
                / (float) subsetNoise.get(i).size());
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
            assert (ch.length == sum.length);
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
        /*try {
            histH.plotHistogram("cch x2 H", "_color_hist_H");
        } catch (IOException ex) {
            Logger.getLogger(UnsupervisedTextureFinder.class.getName()).log(Level.SEVERE, null, ex);
        }*/
        HistogramHolder histS = Histogram.createSimpleHistogram(
            nBins, valuesS, Errors.populateYErrorsBySqrt(valuesS));
        /*try {
            histS.plotHistogram("cch x2 S", "_color_hist_S");
        } catch (IOException ex) {
            Logger.getLogger(UnsupervisedTextureFinder.class.getName()).log(Level.SEVERE, null, ex);
        }*/
        HistogramHolder histV = Histogram.createSimpleHistogram(
            nBins, valuesV, Errors.populateYErrorsBySqrt(valuesV));
        /*try {
            histV.plotHistogram("cch x2 V", "_color_hist_V");
        } catch (IOException ex) {
            Logger.getLogger(UnsupervisedTextureFinder.class.getName()).log(Level.SEVERE, null, ex);
        }*/
        HistogramHolder hist = Histogram.createSimpleHistogram(
            nBins, values, Errors.populateYErrorsBySqrt(values));
        /*try {
            hist.plotHistogram("cch intersections", "_color_hist_");
        } catch (IOException ex) {
            Logger.getLogger(UnsupervisedTextureFinder.class.getName()).log(Level.SEVERE, null, ex);
        }*/

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
        MiscDebug.writeImage(dbg, "_texture_clusters_2_" + debugTag);

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
                float maxSim = Float.NEGATIVE_INFINITY;
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
            if (nIter > 0
                && (Math.abs(totSim - prevTotSim) < (0.01 * nPoints))
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
        MiscDebug.writeImage(dbg, "_texture_clusters_3_" + debugTag);

        List<Set<PairInt>> rList = new ArrayList<Set<PairInt>>();
        //TFloatList critDists = new TFloatArrayList();
        TIntList groupIndexList = new TIntArrayList();
        List<GroupPixelHSV> colors = new ArrayList<GroupPixelHSV>();

        
        int w = img.getWidth();
        PixelHelper ph = new PixelHelper();
            
        for (int i = 0; i < nGroups; ++i) {
            TIntList snIndexes = groupIndexes[i];
            if (snIndexes.size() == 0) {
                continue;
            }
            TLongSet pixIdxs2 = new TLongHashSet();
            int maxX = Integer.MIN_VALUE;
            int maxY = Integer.MIN_VALUE;
            for (int j = 0; j < snIndexes.size(); ++j) {
                Set<PairInt> set = subsetNoise.get(snIndexes.get(j));
                for (PairInt p : set) {
                    int x = p.getX();
                    int y = p.getY();
                    long pixIdx = ph.toPixelIndex(p, w);
                    pixIdxs2.add(pixIdx);
                    if (x > maxX) {
                        maxX = x;
                    }
                    if (y > maxY) {
                        maxY = y;
                    }
                }
            }

            ClusterFinder cFinder = new ClusterFinder(pixIdxs2,
                img.getHeight(), w);
            cFinder.setThreshholdFactor(1.f);
            cFinder.setMinimumNumberInCluster(1);
            cFinder.calculateBackgroundSeparation();
            cFinder.findClusters();
            List<TLongSet> groupList = cFinder.getGroups();
            final int n = groupList.size();

            int maxN = Integer.MIN_VALUE;
            int maxNIdx = -1;
            for (int ii = 0; ii < n; ++ii) {
                int sz = groupList.get(ii).size();
                if (sz > maxN) {
                    maxN = sz;
                    maxNIdx = ii;
                }
            }
            if (maxNIdx != -1) {
                
                int[] xy = new int[2];

                Set<PairInt> set = new HashSet<PairInt>();
                TLongIterator iter3 = groupList.get(maxNIdx).iterator();
                while (iter3.hasNext()) {
                    long pixIdx = iter3.next();
                    ph.toPixelCoords(pixIdx, w, xy);
                    set.add(new PairInt(xy[0], xy[1]));
                }

                if (set.isEmpty()) {
                    continue;
                }

                int critSep = Math.round(cFinder.getBackgroundSeparationHolder()
                    .getXBackgroundSeparation());

                //critDists.add(critDist);
                int patchWidth = critSep;
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
                    FixedSizeSortedVector<IntegerPI> sortedInter = new FixedSizeSortedVector<IntegerPI>(set.size(),
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
                    rPoints[1] = sortedInter.getArray()[set.size() / 2].p;
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
                + debugTag);
        }

        // make hsv filtered copies of image, then extracy greuscale
        List<GreyscaleImage> filteredHSVImgs = new ArrayList<GreyscaleImage>();
        float eFactor = 1.f;//2.f;
        for (int i = 0; i < colors.size(); ++i) {
            GroupPixelHSV hsv = colors.get(i);
            float h0 = hsv.getAvgH() - eFactor * hsv.getStdDevH();
            float h1 = hsv.getAvgH() + eFactor * hsv.getStdDevH();
            float s0 = hsv.getAvgS() - eFactor * hsv.getStdDevS();
            float s1 = hsv.getAvgS() + eFactor * hsv.getStdDevS();
            float v0 = hsv.getAvgV() - eFactor * hsv.getStdDevV();
            float v1 = hsv.getAvgV() + eFactor * hsv.getStdDevV();

            float l0 = hsv.getAvgL() - eFactor * hsv.getStdDevL();
            float l1 = hsv.getAvgL() + eFactor * hsv.getStdDevL();
            float a0 = hsv.getAvgA() - eFactor * hsv.getStdDevA();
            float a1 = hsv.getAvgA() + eFactor * hsv.getStdDevA();
            float b0 = hsv.getAvgB() - eFactor * hsv.getStdDevB();
            float b1 = hsv.getAvgB() + eFactor * hsv.getStdDevB();

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
                if (lab[0] < l0 || lab[0] > l1
                    || lab[1] < a0 || lab[1] > a1
                    || lab[2] < b0 || lab[2] > b1) {
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

        TIntObjectMap<Set<GreyscaleImage>> groupResponseImages = new TIntObjectHashMap<Set<GreyscaleImage>>();
        TIntObjectMap<TexturePatchesAndResponse> groupPatches = 
            new TIntObjectHashMap<TexturePatchesAndResponse>();
        
        FFTUtil fftUtil = new FFTUtil();
        
        for (int i = 0; i < rList.size(); ++i) {

            Set<PairInt> set = rList.get(i);

            int[] minMaxXY = MiscMath.findMinMaxXY(set);

            GreyscaleImage imgPattern = patterns.get(i);

            Complex[][] fftPattern = PhaseCongruencyDetector
                .createLowPassFreqDomainFilter(imgPattern);

            // --- image in the frequency domain convolved with texture patch ----
            Complex[][] freqDomainImageTimesPattern
                = imageProcessor.convolveWithKernel(fftImage, fftPattern);

            // ----- transform that to spatial domain ----
            Complex[][] fComplex =  fftUtil.create2DFFT(
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
                        (int) Math.round(transformedReal[ii][j]));
                }
            }
            MiscDebug.writeImage(kpFreqR2Img, "_freq_spatial_"
                + i + "_" + debugTag);

            GreyscaleImage img3 = filteredHSVImgs.get(i);
            for (int j = 0; j < kpFreqR2Img.getNPixels(); ++j) {
                if (img3.getValue(j) == 0) {
                    kpFreqR2Img.setValue(j, 0);
                }
            }

            MiscDebug.writeImage(kpFreqR2Img, "_freq_spatial_filtered_"
                + i + "_" + debugTag);

            int groupIdx = groupIndexList.get(i);
            Set<GreyscaleImage> rImages = groupResponseImages.get(groupIdx);
            TexturePatchesAndResponse tpar = groupPatches.get(groupIdx);
            if (rImages == null) {
                rImages = new HashSet<GreyscaleImage>();
                groupResponseImages.put(groupIdx, rImages);
                
                tpar = new TexturePatchesAndResponse();
                groupPatches.put(groupIdx, tpar);
                tpar.patches = new ArrayList<GreyscaleImage>();
            }
            rImages.add(kpFreqR2Img);
            tpar.patches.add(imgPattern);
        }
     
        for (int i = 0; i < nGroups; ++i) {
            Set<GreyscaleImage> images = groupResponseImages.get(i);
            if (images == null) {
                continue;
            }
            GreyscaleImage gImage = new GreyscaleImage(img.getWidth(),
                img.getHeight());
            int thresh = 94;
            for (GreyscaleImage g : images) {
                for (int j = 0; j < g.getNPixels(); ++j) {
                    int v = g.getValue(j);
                    //choosing to keep strongest response if more than one
                    if ((v > gImage.getValue(j)) && v >= thresh) {
                        gImage.setValue(j, v);
                    }
                }
            }
            MiscDebug.writeImage(gImage, "_final_textures_" + i
                + "_" + debugTag);
            
            groupPatches.get(i).responseImage = gImage;
        }

        // --- merge any texture images/ with significnt overlap in set picels,
        //     and with similar colors.
        TIntObjectMap<TIntSet> pointsInImageMap = new
            TIntObjectHashMap<TIntSet>(groupPatches.size());
        
        TIntObjectIterator<TexturePatchesAndResponse> iter =
            groupPatches.iterator();
        
        for (int i = 0; i < groupPatches.size(); ++i) {
            
            iter.advance();
            
            int idx = iter.key();
            TexturePatchesAndResponse tPar = iter.value();
            
            GreyscaleImage gImage = tPar.responseImage;
            TIntSet ptSet = new TIntHashSet();
            for (int j = 0; j < gImage.getNPixels(); ++j) {
                if (gImage.getValue(j) > 0) {
                    ptSet.add(j);
                }
            }
            pointsInImageMap.put(idx, ptSet);
        }
        
        TIntSet merged = new TIntHashSet();
        
        int[] groupKeys = pointsInImageMap.keys();
        for (int j = 0; j < groupKeys.length; ++j) {
            int idx0 = groupKeys[j];
            if (merged.contains(idx0)) {
                continue;
            }
            for (int k = (j + 1); k < groupKeys.length; ++k) {
                int idx1 = groupKeys[k];
                TIntSet jIntersectionK = new TIntHashSet(pointsInImageMap.get(idx0));
                TIntSet setJ = new TIntHashSet(pointsInImageMap.get(idx0));
                // j - k subtracts the intersection, then J minu that leaves intersetion
                setJ.removeAll(pointsInImageMap.get(idx1));
                jIntersectionK.removeAll(setJ);
                int nI = jIntersectionK.size();
                float f0 = (float)nI/(float)pointsInImageMap.get(idx0).size();
                float f1 = (float)nI/(float)pointsInImageMap.get(idx1).size();
                System.out.println("j intersection=" + nI + 
                    " frac of J = " + f0 + " frac of k = " + f1);
                if (f0 > 0.05 || f1 > 0.05) {
                    int[][] ch0 = chhist.histogramCIELAB(img,
                        pointsInImageMap.get(idx0));
                    int[][] ch1 = chhist.histogramCIELAB(img,
                        pointsInImageMap.get(idx1));
                    float intersection = chhist.intersection(ch0, ch1);
                    
                    System.out.println("   CIELAB intersection=" + intersection);
                    
                    //TODO: revisit this
                    if (intersection < 0.45) {
                        continue;
                    }
                    
                    TexturePatchesAndResponse tPar0 = groupPatches.get(idx0);
                    TexturePatchesAndResponse tPar1 = groupPatches.remove(idx1);
                    
                    tPar0.patches.addAll(tPar1.patches);
                    for (int pixIdx = 0; pixIdx < tPar0.responseImage.getNPixels();
                        ++pixIdx) {
                        int v1 = tPar1.responseImage.getValue(pixIdx);
                        if (v1 > 0) {
                            int v0 = tPar0.responseImage.getValue(pixIdx);
                            if (v1 > v0) {
                                tPar0.responseImage.setValue(pixIdx, v1);
                            }
                        }
                    }
                    
                    merged.add(idx1);
                }
            }
        }
        
        return groupPatches.values(new TexturePatchesAndResponse[groupPatches.size()]);
    }
  
    public class TexturePatchesAndResponse {
        public List<GreyscaleImage> patches;
        public GreyscaleImage responseImage;
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
