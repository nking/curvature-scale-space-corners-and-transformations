package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.IntegralHistograms;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.MiscMath;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import algorithms.util.TwoDIntArray;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TIntHashSet;
import gnu.trove.set.hash.TLongHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class HOGRegionsManager {
    
    private static float eps = 0.000001f;
    
    // 9 is default
    private final int nAngleBins;
    
    private final int nHistBins;

    // 6 x 6 is recommended
    private final int N_PIX_PER_CELL_DIM;

    // 2x2 or 3x3 is recommended
    private final int N_CELLS_PER_BLOCK_DIM;

    /**
     * key = region index
     */
    private final TIntObjectMap<TwoDIntArray> regionIndexHistHOG = 
        new TIntObjectHashMap<TwoDIntArray>();
    private final TIntObjectMap<TwoDIntArray> regionIndexHistHCPT = 
        new TIntObjectHashMap<TwoDIntArray>();
    private final TIntObjectMap<TwoDIntArray> regionIndexHistHGS = 
        new TIntObjectHashMap<TwoDIntArray>();
    private final TIntObjectMap<OneDIntArray> regionIndexMinMaxXY =
        new TIntObjectHashMap<OneDIntArray>();
    private final TIntObjectMap<Canonicalizer.CRegion> regionIndexRegions;
    //In reference frame of subImage
    private final TIntObjectMap<TLongSet> regionCoords 
        = new TIntObjectHashMap<TLongSet>();
    
    private static int maskValue = 0;

    private final int w;
    private final int h;

    private boolean debug = false;
    
    public HOGRegionsManager(TIntObjectMap<Canonicalizer.CRegion> regionMap, 
        int imageWidth, int imageHeight, int nCellsPerDim, int nPixPerCellDim, 
        int nAngleBins, int nHCPTHGSBins) {

        regionIndexRegions = regionMap;
        this.nAngleBins = nAngleBins;
        nHistBins = nHCPTHGSBins;
        N_PIX_PER_CELL_DIM = nPixPerCellDim;
        N_CELLS_PER_BLOCK_DIM = nCellsPerDim;
        w = imageWidth;
        h = imageHeight;
    }

    public void setToDebug() {
        debug = true;
    }
    
    int[][] createHOGHistogram(GreyscaleImage gsImg, 
        Collection<PairInt> points,
        int[] outputMinMaxXY, TLongSet outputRefFramePixs) {

        if (w != gsImg.getWidth() || h != gsImg.getHeight()) {
            throw new IllegalArgumentException(
            "gsImg must have size same as constructor args");
        }

        GreyscaleImage gsImg2 = HOGUtil.createAndMaskSubImage(gsImg, 
            maskValue, points, outputMinMaxXY, outputRefFramePixs);
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        GreyscaleImage[] gXgY = 
            imageProcessor.createCentralDifferenceGradients(gsImg2);
        GreyscaleImage theta = imageProcessor.computeTheta180(gXgY[0], gXgY[1]);
        GreyscaleImage gXY = 
            imageProcessor.combineConvolvedImages(gXgY[0], gXgY[1]);
            
        return createHOGHistogram(gXY, theta);        
    }
    
    int[][] createHOGHistogram(GreyscaleImage gXY, GreyscaleImage theta) {
        
        GradientIntegralHistograms gh = new GradientIntegralHistograms();
        
        int[][] histograms = gh.createHistograms(gXY, theta, nAngleBins);

        //apply a windowed sum across the integral image.
        // result is that at each pixel is a histogram holding the sum of histograms
        //    from the surrounding N_PIX_PER_CELL_DIM window. 
        gh.applyWindowedSum(histograms, gXY.getWidth(), gXY.getHeight(), 
            N_PIX_PER_CELL_DIM);

        return histograms;
    }
    
    int[][] createHCPTHistogram(GreyscaleImage ptImg, TLongSet regionPixelCoords) {

        int w2 = ptImg.getWidth();
        int h2 = ptImg.getHeight();
        
        PolarThetaIntegralHistograms gh = new PolarThetaIntegralHistograms();
        
        int[][] histograms = gh.createHistograms(ptImg, 
            regionPixelCoords, nHistBins);

        //apply a windowed avg across the integral image
        gh.applyWindowedSum(histograms, w2, h2, N_PIX_PER_CELL_DIM);
        
        return histograms;
    }
    
    /**
     * add a region to the maps and construct a 2D histogram for it.
     * Note that the region index is taken from cRegion.dataIdx.
     * 
     * Also note that the method assumes that the sub-images have already 
     * been masked.
     * 
     * Also note that if the cRegion has fewer than 9 points, the region is not added.
     * 
     * @param subImageGSImg
     * @param subImagePTImg
     * @param subImageGradient
     * @param subImageTheta
     * @param cRegion
     * @param fullFrameMinMaxXY
     * @param subImagePixelCoords
     */
    public void addARegion(
        GreyscaleImage subImageGSImg, GreyscaleImage subImagePTImg,
        GreyscaleImage subImageGradient, GreyscaleImage subImageTheta,
        Canonicalizer.CRegion cRegion, 
        int[] fullFrameMinMaxXY, TLongSet subImagePixelCoords) {
        
        int rIndex = cRegion.dataIdx;
        if (rIndex == -1) {
            throw new IllegalArgumentException("cRegion.dataIdx must be >=");
        }
        
        if (cRegion.offsetsToOrigCoords.size() < 9) {
            return;
        }
        
        int w2 = subImageGSImg.getWidth();
        int h2 = subImageGSImg.getHeight();
        
        regionCoords.put(rIndex, subImagePixelCoords);
        regionIndexMinMaxXY.put(rIndex, new OneDIntArray(fullFrameMinMaxXY));        
        regionIndexRegions.putIfAbsent(rIndex, cRegion);
            
        // non-region pixels are excluded because magnitude is zero
        int[][] histogramsHOG = createHOGHistogram(subImageGradient, subImageTheta);        
        regionIndexHistHOG.put(rIndex, new TwoDIntArray(histogramsHOG));
        
        // exclude non-region pixels:
        int[][] histogramsHCPT = createHCPTHistogram(subImagePTImg, subImagePixelCoords);
        regionIndexHistHCPT.put(rIndex, new TwoDIntArray(histogramsHCPT));
        
        IntegralHistograms gh = new IntegralHistograms();
        int[][] histogramsHGS = gh.create(subImageGSImg, subImagePixelCoords, 
            0, 255, nHistBins);
        //apply a windowed avg across the integral image
        gh.applyWindowedSum(histogramsHGS, w2, h2, N_PIX_PER_CELL_DIM);
        regionIndexHistHGS.put(rIndex, new TwoDIntArray(histogramsHGS));        
    }
    
    private static enum TYPE {
        HOG, HCPT, HGS
    };
    
    /**
     * CAVEAT: small amount of testing done, not yet throughly tested.
     *
     * extract the block surrounding the coordinates.
     * the number of pixels in a cell and the number of cells in 
     * block were set during
     * construction.
     * 
     * The feature is nAngleBins in length for 180 degrees
     * and the bin with the largest value
     * is the bin holding the angle perpendicular to the windowed point.
     * (for example: a horizontal line, the feature of a point on the
     * line has largest bin being the 90 degrees bin).
     *
     * @param rIndex
     * @param x
     * @param y
     * @param outHist
     * @return true if method succeeded, else false.  can return false if the
     * addARegion failed due to having fewer than 9 pixels in the CRegion
     * for rIndex.
     */
    public boolean extractBlockHOG(int rIndex, int x, int y, int[] outHist) {
        return extractBlock(TYPE.HOG, rIndex, x, y, outHist);
    }
    
    /**
     * extract the block surrounding the feature.
     * the number of pixels in a cell and the number of cells in block were set during
     * construction.
     * @param rIndex
     * @param x
     * @param y
     * @param outHist
     */
    public boolean extractBlockHCPT(int rIndex, int x, int y, int[] outHist) {
        return extractBlock(TYPE.HCPT, rIndex, x, y, outHist);
    }
    
    /**
     * extract the block surrounding the feature.
     * the number of pixels in a cell and the number of cells in block were set during
     * construction.
     * @param rIndex
     * @param x
     * @param y
     * @param outHist
     */
    public boolean extractBlockHGS(int rIndex, int x, int y, int[] outHist) {
        return extractBlock(TYPE.HGS, rIndex, x, y, outHist);
    }
    
    /**
     * CAVEAT: small amount of testing done, not yet throughly tested.
     *
     * extract the block surrounding the feature.
     * the number of pixels in a cell and the number of cells in 
     * block were set during
     * construction.
     * 
     * The feature is nAngleBins in length for 180 degrees
     * and the bin with the largest value
     * is the bin holding the angle perpendicular to the windowed point.
     * (for example: a horizontal line, the feature of a point on the
     * line has largest bin being the 90 degrees bin).
     *
     * @param rIndex
     * @param x
     * @param y
     * @param outHist
     */
    private boolean extractBlock(TYPE type, int rIndex, int x, int y, int[] outHist) {

        if (!regionIndexRegions.containsKey(rIndex)) {
            return false;
        }
        
        if ((type.equals(TYPE.HOG) && outHist.length != nAngleBins) ||
            (!type.equals(TYPE.HOG) && outHist.length != nHistBins)) {
            throw new IllegalArgumentException("outHist.length != expected");
        }

        if (x < 0 || y < 0 || x >= w || y >= h) {
            throw new IllegalArgumentException("x or y is out of bounds of "
                + "original image");
        }
        
        if (!regionCoords.containsKey(rIndex) || 
            !regionIndexMinMaxXY.containsKey(rIndex)) {
            return false;
        }

        // uses the block normalization recomended by Dalal & Triggs,
        //   the summary of histogram counts over all cells
        //   is used to normalize each cell by that sum.

        int[][] hist = null;
        switch (type) {
            case HOG: hist = regionIndexHistHOG.get(rIndex).a; break;
            case HCPT: hist = regionIndexHistHCPT.get(rIndex).a; break;
            case HGS: hist = regionIndexHistHGS.get(rIndex).a; break;
            default: break;
        }
        
        int[] minMaxXY = regionIndexMinMaxXY.get(rIndex).a;
        TLongSet rCoords = regionCoords.get(rIndex);
        
        int width = minMaxXY[1] - minMaxXY[0] + 1;
        int height = minMaxXY[3] - minMaxXY[2] + 1;
        int xOffset = minMaxXY[0];
        int yOffset = minMaxXY[2];
        
        int nH = N_CELLS_PER_BLOCK_DIM * N_CELLS_PER_BLOCK_DIM;
        
        long[] tmp = new long[outHist.length];
        
        int count = 0;
        
        for (int cX = 0; cX < N_CELLS_PER_BLOCK_DIM; ++cX) {

            int cXOff = -(N_CELLS_PER_BLOCK_DIM/2) + cX;

            int x2 = x + (cXOff * N_PIX_PER_CELL_DIM);

            int xSub = x2 - xOffset;
            
            if ((xSub + N_PIX_PER_CELL_DIM - 1) < 0) {
                break;
            } else if (xSub < 0) {
                xSub = 0;
            } else if (xSub >= width) {
                break;
            }

            for (int cY = 0; cY < N_CELLS_PER_BLOCK_DIM; ++cY) {

                int cYOff = -(N_CELLS_PER_BLOCK_DIM/2) + cY;

                int y2 = y + (cYOff * N_PIX_PER_CELL_DIM);

                int ySub = y2 - yOffset;
                
                if ((ySub + N_PIX_PER_CELL_DIM - 1) < 0) {
                    break;
                } else if (ySub < 0) {
                    ySub = 0;
                } else if (ySub >= height) {
                    break;
                }
                
                int pixIdx = (ySub * width) + xSub;
                
                if (!rCoords.contains(pixIdx)) {
                    continue;
                }

                add(tmp, hist[pixIdx]);
            }
        }

        //System.out.println("  s=" + Arrays.toString(tmp));
        long blockTotal = 0;        
        for (int i = 0; i < tmp.length; ++i) {
            blockTotal += tmp[i];
        }
        double norm = 255./(blockTotal + eps);
        for (int i = 0; i < outHist.length; ++i) {
            outHist[i] = (int)(tmp[i]*norm);
        }
        //System.out.println("->out=" + Arrays.toString(outHist));

        return true;
    }
    
    public static void populateRegionsIfNeeded(
        TIntObjectMap<Canonicalizer.RegionPoints> regionPointsMap, float scale,
        HOGRegionsManager hogMgs,
        GreyscaleImage gsImg, GreyscaleImage ptImg) {

        TIntObjectMap<Canonicalizer.CRegion> cRegionMapReference 
            = hogMgs.getRegionIndexRegions();
        
        if (!cRegionMapReference.isEmpty()) {
            //System.out.println("*map is not empty, so not adding more.");
            return;
        }
        
        int w = gsImg.getWidth();
        int h = gsImg.getHeight();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Canonicalizer cn = new Canonicalizer();
        int[] outputMinMaxXY = new int[4];
        
        TIntObjectIterator<Canonicalizer.RegionPoints> iter 
            = regionPointsMap.iterator();
        
        int n = regionPointsMap.size();
        
        for (int i = 0; i < n; ++i) {
            iter.advance();
            
            Canonicalizer.RegionPoints regionPoints = iter.value();
            
            if (regionPoints.points.size() < 9) {
                continue;
            }
            
            List<Canonicalizer.CRegion> crs = cn.canonicalizeRegionsToFit(regionPoints, 
                outputMinMaxXY);
            
            //System.out.println("created " + crs.size() + " from RegionPoints " + i);
            
            for (Canonicalizer.CRegion cr : crs) {
                int rIdx = cRegionMapReference.size();
                cr.dataIdx = rIdx;
                
                if (cr.offsetsToOrigCoords.size() < 9) {
                    continue;
                }
                
                Canonicalizer.CRegion crScaled = cr.createNewDividedByScale(scale);
                
                if (crScaled.offsetsToOrigCoords.size() < 9) {
                    continue;
                }
                
                int[] minMaxXY = crScaled.getMinMaxXY();
                
                //cRegionMapReference.put(rIdx, crScaled);
                
                TLongSet outputRefFramePixs = new TLongHashSet();
                Set<PairInt> coords = crScaled.extractCoords();
                GreyscaleImage gsImg2 = HOGUtil.createAndMaskSubImage(gsImg, 
                    maskValue, coords, outputMinMaxXY, 
                    outputRefFramePixs);

                GreyscaleImage[] gXgY = 
                    imageProcessor.createCentralDifferenceGradients(gsImg2);
                GreyscaleImage theta2 = imageProcessor.computeTheta180(gXgY[0], gXgY[1]);
                GreyscaleImage gXY2 = 
                    imageProcessor.combineConvolvedImages(gXgY[0], gXgY[1]);
               
                GreyscaleImage ptImg2 = HOGUtil.createAndMaskSubImage2(ptImg, 
                    maskValue, coords, outputMinMaxXY, outputRefFramePixs);
                
                hogMgs.addARegion(gsImg2, ptImg2, gXY2, theta2, crScaled, 
                    outputMinMaxXY, outputRefFramePixs);
                assert(cRegionMapReference.get(rIdx).equals(crScaled));
            }
        }
    }
    
    /**
     * CAVEAT: small amount of testing done, not yet throughly tested.
     *
     * calculate the intersection of histA and histB which have already
     * been normalized to the same scale.
     * A result of 0 is maximally dissimilar and a result of 1 is maximally similar.
     *
     * The orientations are needed to compare the correct rotated bins to one another.
     * Internally, orientation of 90 leads to no shift for rotation,
     * and orientation near 0 results in rotation of nBins/2, etc...
     *
     * Note that an orientation of 90 is a unit vector from x,y=0,0 to
     * x,y=0,1.
     * 
     * Note also that you may want to try the rotation of oppossite direction.
     *
     * @param histA
     * @param orientationA
     * @param histB
     * @param orientationB
     * @return
     */
    public float intersection(int[] histA, int orientationA, int[] histB,
        int orientationB) {

        return HOGUtil.intersection(histA, orientationA, histB, orientationB);
    }
    
    /**
     * CAVEAT: small amount of testing done, not yet throughly tested.
     *
     * calculate the intersection of histA and histB which have already
     * been normalized to the same scale.
     * A result of 0 is maximally dissimilar and a result of 1 is maximally similar.
     *
     * @param histA
     * @param histB
     * @return
     */
    public float intersection(int[] histA, int[] histB) {

        return HOGUtil.intersection(histA, histB);
    }
    
    /**
     * CAVEAT: small amount of testing done, not yet throughly tested.
     * 
     * calculate the intersection of histA and histB which have already
     * been normalized to the same scale.
     * A result of 0 is maximally dissimilar and a result of 1 is maximally similar.
     * 
     The orientations are needed to compare the correct rotated bins to one another.
     * Internally, orientation of 90 leads to no shift for rotation,
     * and orientation near 0 results in rotation of nBins/2, etc...
     *
     * Note that an orientation of 90 is a unit vector from x,y=0,0 to
     * x,y=0,1.
     *
     * @param histA
     * @param orientationA
     * @param histB
     * @param orientationB
     * @return difference, error
     */
    public float[] diff(int[] histA, int orientationA, int[] histB,
        int orientationB) {

        return HOGUtil.diff(histA, orientationA, histB, orientationB);
    }
    
    private int sumCounts(int[] hist) {

        int sum = 0;
        for (int v : hist) {
            sum += v;
        }

        return sum;
    }

    public static void add(int[] addTo, int[] addFrom) {
        for (int i = 0; i < addTo.length; ++i) {
            addTo[i] += addFrom[i];
        }
    }

    public static void add(long[] addTo, int[] addFrom) {
        for (int i = 0; i < addTo.length; ++i) {
            addTo[i] += addFrom[i];
        }
    }
    
    public static void add(long[] addTo, long[] addFrom) {
        for (int i = 0; i < addTo.length; ++i) {
            addTo[i] += addFrom[i];
        }
    }

     
    /**
     * NOT READY FOR USE.
     *
     * @param rIndex
     * @return
     */
    public int calculateDominantOrientation(int rIndex) {

        int[][] hist = regionIndexHistHOG.get(rIndex).a;
        int[] minMaxXY = getRegionIndexMinMaxXY().get(rIndex).a;
        
        int width = minMaxXY[1] - minMaxXY[0] + 1;
        //int height = minMaxXY[3] - minMaxXY[2] + 1;
        int xOffset = minMaxXY[0];
        int yOffset = minMaxXY[2];
        
        Canonicalizer.CRegion r = regionIndexRegions.get(rIndex);
        
        long[] combined = new long[nAngleBins];

        for (Map.Entry<PairInt, PairInt> entry : r.offsetsToOrigCoords.entrySet()) {
            PairInt xy = entry.getValue();
            int pixIdx = ((xy.getY() - yOffset) * width) + (xy.getX() - xOffset);
            add(combined, hist[pixIdx]);
        }
        
        TIntList maxIdxs = new TIntArrayList();
        long maxValue = Long.MIN_VALUE;
        for (int i = 0; i < combined.length; ++i) {
            long v = combined[i];
            if (v > maxValue) {
                maxIdxs.clear();
                maxIdxs.add(i);
                maxValue = v;
            } else if (v == maxValue) {
                maxIdxs.add(i);
            }
        }
        
        int binWidth = 180 / nAngleBins;

        //TODO: revise to fit max peak using neighboring bins
        if (maxIdxs.size() == 1) {
            return Math.round((maxIdxs.get(0) + 0.5f) * binWidth);
        }

        //TODO: revise this
        if (maxIdxs.size() == 2) {

            // min of average with and without a wraparound phase
            float ang0 = (maxIdxs.get(0) + 0.5f) * binWidth;
            float ang1 = (maxIdxs.get(1) + 0.5f) * binWidth;

            if (ang0 < ang1) {
                float diff0 = 360 + ang0 - ang1;
                float diff1 = ang1 - ang0;
                if (diff0 < diff1) {
                    return Math.round(0.5f * (360 + ang0 + ang1));
                }

                return Math.round(0.5f * (ang0 + ang1));
            }

            float diff0 = 360 + ang1 - ang0;
            float diff1 = ang0 - ang1;
            if (diff0 < diff1) {
                return Math.round(0.5f * (360 + ang1 + ang0));
            }

            return Math.round(0.5f * (ang0 + ang1));

        } else {

            // average of all indexes

            double[] angles = new double[maxIdxs.size()];
            for (int i = 0; i < maxIdxs.size(); ++i) {
                angles[i] = ((maxIdxs.get(i) + 0.5f) * binWidth);
            }

            double angleAvg =
                AngleUtil.calculateAverageWithQuadrantCorrections(
                    angles, false);

            return (int)Math.round(angleAvg);
        }
    }

    public TIntSet calculateDominantOrientations(int rIndex) {

        int[][] hist = regionIndexHistHOG.get(rIndex).a;
        int[] minMaxXY = getRegionIndexMinMaxXY().get(rIndex).a;
        
        int width = minMaxXY[1] - minMaxXY[0] + 1;
        //int height = minMaxXY[3] - minMaxXY[2] + 1;
        int xOffset = minMaxXY[0];
        int yOffset = minMaxXY[2];
        
        Canonicalizer.CRegion r = regionIndexRegions.get(rIndex);
        
        Set<PairInt> xys = new HashSet<PairInt>();
        
        long[] combined = new long[nAngleBins];

        for (Map.Entry<PairInt, PairInt> entry : r.offsetsToOrigCoords.entrySet()) {
            PairInt xy = entry.getValue();
            int pixIdx = ((xy.getY() - yOffset) * width) + (xy.getX() - xOffset);
            add(combined, hist[pixIdx]);
            xys.add(xy);
        }

        int maxIdx = MiscMath.findYMaxIndex(combined);

        if (maxIdx == -1) {
            throw new IllegalArgumentException("histogram is full of "
                + " min value longs");
        }

        // if any bins have values within 80% of max, add to maxIdxs
        TIntList maxIdxs = new TIntArrayList();

        long max = combined[maxIdx];
        double limit = 0.8 * max;

        for (int i = 0; i < combined.length; ++i) {
            long v = combined[i];
            if (v >= limit) {
                maxIdxs.add(i);
            }
        }

        int binWidth = 180 / nAngleBins;

        TIntSet orientations = new TIntHashSet();

        for (int i = 0; i < maxIdxs.size(); ++i) {
            int idx = maxIdxs.get(i);
            int angle = Math.round((idx + 0.5f) * binWidth);
            orientations.add(angle);
        }
                        
        //NOTE: adding an orientation for the center of points
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        PairInt xyCen = ch.calculateXYCentroids2(xys);
        if (xys.contains(xyCen)) {
            int pixIdx = ((xyCen.getY() - yOffset) * width) + (xyCen.getX() - xOffset);
            maxIdx = MiscMath.findYMaxIndex(hist[pixIdx]);
            if (maxIdx > -1) {
                int angle = Math.round((maxIdx + 0.5f) * binWidth);
                orientations.add(angle);
            }
        }
        
        return orientations;
    }

    public int getImageWidth() {
        return w;
    }
    
    public int getImageHeight() {
        return h;
    }

    public int getNumberOfBins() {
        return nAngleBins;
    }

    public TIntObjectMap<Canonicalizer.CRegion> getRegionIndexRegions() {
        return regionIndexRegions;
    }

    /**
     * @return the regionIndexMinMaxXY
     */
    public TIntObjectMap<OneDIntArray> getRegionIndexMinMaxXY() {
        return regionIndexMinMaxXY;
    }
}
