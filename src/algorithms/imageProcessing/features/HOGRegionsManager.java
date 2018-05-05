package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.IntegralHistograms;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import algorithms.imageProcessing.features.mser.Canonicalizer.RegionPoints;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
import algorithms.util.TwoDIntArray;
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
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
    private final TIntObjectMap<OneDIntArray> regionIndexMinMaxXY =
        new TIntObjectHashMap<OneDIntArray>();
    private final TIntObjectMap<Canonicalizer.CRegion> regionIndexRegions;
    //In reference frame of subImage
    private final TIntObjectMap<TLongSet> regionCoords 
        = new TIntObjectHashMap<TLongSet>();
    
    /**
     * key = minmaxXY for region
     */
    private final Map<OneDIntArray, TwoDIntArray> histHOGMap = 
        new HashMap<OneDIntArray, TwoDIntArray>();
    private final Map<OneDIntArray,TwoDIntArray> histHCPTMap = 
        new HashMap<OneDIntArray, TwoDIntArray>();
    private final Map<OneDIntArray,TwoDIntArray> histHGSMap = 
        new HashMap<OneDIntArray, TwoDIntArray>();
    
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
        
        return HOGUtil.createHOGHistogram(gXY, theta, nAngleBins, N_PIX_PER_CELL_DIM);
    }
    
    int[][] createHCPTHistogram(GreyscaleImage ptImg, TLongSet regionPixelCoords) {

        return HOGUtil.createHCPTHistogram(ptImg, regionPixelCoords, nHistBins, 
            N_PIX_PER_CELL_DIM);
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
     * @param cRegion a region whose coordinates are in the reference frame of
     *    the uncropped image (that is, the sumImages).
     * @param minMaxXYUnCroppedRefFrame these are the bounds of the points in
     *    cRegion but the reference frame is the uncropped image (that is, the sumImages)
     *    to pass to this method the offsets from the uncropped image.
     * @param subImagePixelCoords
     */
    public void addARegion(
        GreyscaleImage subImageGSImg, GreyscaleImage subImagePTImg,
        GreyscaleImage subImageGradient, GreyscaleImage subImageTheta,
        Canonicalizer.CRegion cRegion, 
        int[] minMaxXYUnCroppedRefFrame, TLongSet subImagePixelCoords) {
        
        int rIndex = cRegion.dataIdx;
        if (rIndex == -1) {
            throw new IllegalArgumentException("cRegion.dataIdx must be >=");
        }
        
        if (cRegion.getOffsetsToOrigCoords().size() < 9) {
            return;
        }
        
        int w2 = subImageGSImg.getWidth();
        int h2 = subImageGSImg.getHeight();
        
        regionCoords.putIfAbsent(rIndex, subImagePixelCoords);
        regionIndexMinMaxXY.putIfAbsent(rIndex, new OneDIntArray(minMaxXYUnCroppedRefFrame));        
        regionIndexRegions.putIfAbsent(rIndex, cRegion);
        
        int[] minMaxXY = cRegion.getMinMaxXY();
        assert(minMaxXY[0] < Integer.MAX_VALUE);
        assert(minMaxXY[2] < Integer.MAX_VALUE);
        assert(minMaxXY[1] > Integer.MIN_VALUE);
        assert(minMaxXY[3] > Integer.MIN_VALUE);
        
        OneDIntArray key = new OneDIntArray(minMaxXY);
        
        if (!histHOGMap.containsKey(key)) {
            // non-region pixels are excluded because magnitude is zero
            int[][] histogramsHOG = createHOGHistogram(subImageGradient, subImageTheta);        
            histHOGMap.put(key, new TwoDIntArray(histogramsHOG));
        }     
        if (!histHCPTMap.containsKey(key)) {
            // exclude non-region pixels:
            int[][] histogramsHCPT = createHCPTHistogram(subImagePTImg, subImagePixelCoords);
            histHCPTMap.put(key, new TwoDIntArray(histogramsHCPT));
        }
        if (!histHGSMap.containsKey(key)) {
            int[][] histogramsHGS = HOGUtil.createHGSHistogram(
                subImageGSImg, subImagePixelCoords, nHistBins, N_PIX_PER_CELL_DIM);
            histHGSMap.put(key, new TwoDIntArray(histogramsHGS));
        }
        assert(histHOGMap.size() > 0);
        assert(histHCPTMap.size() > 0);
        assert(histHGSMap.size() > 0);
    }

    private TLongSet transformByOffsets(TLongSet pixs, int w, int w2, 
        int xOffset2, int yOffset2) {
        
        PixelHelper ph = new PixelHelper();
        int[] xy = new int[2];
        TLongSet pixs2 = new TLongHashSet();
        TLongIterator iter = pixs.iterator();
        while (iter.hasNext()) {
            long pix = iter.next();
            ph.toPixelCoords(pix, w, xy);
            long pix2 = ph.toPixelIndex(
                xy[0] - xOffset2, xy[1] - yOffset2, w2);
            pixs2.add(pix2);
        }
        return pixs2;
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
     * @return 
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
        
        OneDIntArray key = regionIndexMinMaxXY.get(rIndex);

        // uses the block normalization recomended by Dalal & Triggs,
        //   the summary of histogram counts over all cells
        //   is used to normalize each cell by that sum.

        int[][] hist = null;
        switch (type) {
            case HOG: hist = histHOGMap.get(key).a; break;
            case HCPT: hist = histHCPTMap.get(key).a; break;
            case HGS: hist = histHGSMap.get(key).a; break;
            default: break;
        }
        
        int[] minMaxXY = regionIndexMinMaxXY.get(rIndex).a;
        TLongSet rCoords = regionCoords.get(rIndex);
        
        return extractBlock(hist, rCoords, minMaxXY, x, y, outHist);
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
     * @param x
     * @param y
     * @param outHist
     */
    private boolean extractBlock(int[][] hist, TLongSet subImagePixelIndexes, 
        int[] minMaxXY,
        int x, int y, int[] outHist) {
        
        if (x < 0 || y < 0 || x >= w || y >= h) {
            throw new IllegalArgumentException("x or y is out of bounds of "
                + "original image");
        }
        
        // uses the block normalization recomended by Dalal & Triggs,
        //   the summary of histogram counts over all cells
        //   is used to normalize each cell by that sum.

        int width = minMaxXY[1] - minMaxXY[0] + 1;
        int height = minMaxXY[3] - minMaxXY[2] + 1;
        int xOffset = minMaxXY[0];
        int yOffset = minMaxXY[2];
        
        int nH = N_CELLS_PER_BLOCK_DIM * N_CELLS_PER_BLOCK_DIM;
        
        long[] tmp = new long[outHist.length];
                
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
                
                if (!subImagePixelIndexes.contains(pixIdx)) {
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
    
    /**
     * populate the instance variables if not already.
     * 
     * @param regionPointsMapNotScaled
     * @param scale
     * @param hogMgs
     * @param gsImgScaled
     * @param ptImgScaled 
     */
    public void populateRegionsIfNeeded(
        TIntObjectMap<Canonicalizer.RegionPoints> regionPointsMapNotScaled, 
        float scale, HOGRegionsManager hogMgs, GreyscaleImage gsImgScaled, 
        GreyscaleImage ptImgScaled) {

        TIntObjectMap<Canonicalizer.CRegion> cRegionMapReference 
            = hogMgs.getRegionIndexRegions();
        
        if (!cRegionMapReference.isEmpty()) {
            System.out.println("*map is not empty, so not adding more.");
            return;
        }
        
        int w = gsImgScaled.getWidth();
        int h = gsImgScaled.getHeight();
        
        PixelHelper ph = new PixelHelper();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        Canonicalizer cn = new Canonicalizer();
        
        TIntObjectIterator<Canonicalizer.RegionPoints> iter 
            = regionPointsMapNotScaled.iterator();
        
        int n = regionPointsMapNotScaled.size();
        
        int ts = MiscDebug.getCurrentTimeFormatted();
        
        int count = 0;
        
        for (int i = 0; i < n; ++i) {
            iter.advance();
            
            RegionPoints regionPointsNotScaled = iter.value();
            
            int nr = regionPointsNotScaled.accX.size();
            
            if (nr < 9) {
                continue;
            }
            
            assert(regionPointsNotScaled.accX.size() == 
                regionPointsNotScaled.accY.size());
            
            RegionPoints regionPointsScaled = 
                regionPointsNotScaled.createNewDividedByScale(
                scale, w - 1, h - 1);
            
            assert(!regionPointsScaled.accX.isEmpty());
            assert(!regionPointsScaled.accY.isEmpty());
            
            int[] minMaxXY = regionPointsScaled.getMinMaxXY();
            assert(minMaxXY[0] != Integer.MAX_VALUE);
            assert(minMaxXY[2] != Integer.MAX_VALUE);
            assert(minMaxXY[1] != Integer.MIN_VALUE);
            assert(minMaxXY[3] != Integer.MIN_VALUE);
            int w2 = minMaxXY[1] - minMaxXY[0] + 1;
            int h2 = minMaxXY[3] - minMaxXY[2] + 1;
            if (w2 < 3 || h2 < 3) {
                continue;
            }
            
            TLongSet refFramePixs = regionPointsScaled
                .createAccPixelCoords(w);
            assert(!refFramePixs.isEmpty());
            
            /*if (true) {
                Image tmp;
                int[] xy = new int[2];
                tmp = gsImgScaled.copyToColorGreyscale();
                for (int ii = 0; ii < regionPointsScaled.accX.size(); ++ii) {
                    ImageIOHelper.addPointToImage(
                        regionPointsScaled.accX.get(ii), 
                        regionPointsScaled.accY.get(ii), 
                        tmp, 0, 255, 0, 0);
                };
                MiscDebug.writeImage(tmp, "_" + ts + "_rpdebug0_");
                tmp = gsImgScaled.copyToColorGreyscale();
                TLongIterator iter2 = refFramePixs.iterator();
                while (iter2.hasNext()) {
                    long pixIdx = iter2.next();
                    ph.toPixelCoords(pixIdx, w, xy);
                    ImageIOHelper.addPointToImage(
                        xy[0], xy[1], 
                        tmp, 0, 255, 0, 0);
                };
                MiscDebug.writeImage(tmp, "_" +  "_rpdebug1_");
            }*/
            
            refFramePixs = transformByOffsets(refFramePixs, w, w2, minMaxXY[0], minMaxXY[2]);
                    
            GreyscaleImage gsImg2 = gsImgScaled.subImage2(
                minMaxXY[0], minMaxXY[1], minMaxXY[2], minMaxXY[3]);
            
            /*if (true) {
                Image tmp;
                int[] xy = new int[2];
                tmp = gsImg2.copyToColorGreyscale();
                TLongIterator iter2 = refFramePixs.iterator();
                while (iter2.hasNext()) {
                    long pixIdx = iter2.next();
                    ph.toPixelCoords(pixIdx, w2, xy);
                    ImageIOHelper.addPointToImage(
                        xy[0], xy[1], 
                        tmp, 0, 255, 0, 0);
                };
                MiscDebug.writeImage(tmp, "_" +  "_rpdebug2_");
            }*/
            
            GreyscaleImage[] gXgY = 
                imageProcessor.createCentralDifferenceGradients(gsImg2);
            GreyscaleImage theta2 = imageProcessor.computeTheta180(gXgY[0], gXgY[1]);
            GreyscaleImage gXY2 = 
                imageProcessor.combineConvolvedImages(gXgY[0], gXgY[1]);
            GreyscaleImage ptImg2 = HOGUtil.createAndMaskSubImage2(ptImgScaled, 
                maskValue, minMaxXY, refFramePixs);
            
            // the histograms are stored by bounds minMaxXY
            
            OneDIntArray key = new OneDIntArray(minMaxXY);
        
            // non-region pixels are excluded because magnitude is zero
            int[][] histogramsHOG = createHOGHistogram(gXY2, theta2);        
            histHOGMap.putIfAbsent(key, new TwoDIntArray(histogramsHOG));
        
            // exclude non-region pixels:
            int[][] histogramsHCPT = createHCPTHistogram(ptImg2, refFramePixs);
            histHCPTMap.putIfAbsent(key, new TwoDIntArray(histogramsHCPT));
        
            IntegralHistograms gh = new IntegralHistograms();
            int[][] histogramsHGS = gh.create(gsImg2, refFramePixs, 
                0, 255, nHistBins);
            //apply a windowed avg across the integral image
            gh.applyWindowedSum(histogramsHGS, w2, h2, N_PIX_PER_CELL_DIM);
            histHGSMap.putIfAbsent(key, new TwoDIntArray(histogramsHGS));
            
            int[] hogHistCenter = new int[nAngleBins];
            extractBlock(histogramsHOG, refFramePixs, minMaxXY, 
                regionPointsScaled.ellipseParams.xC, 
                regionPointsScaled.ellipseParams.yC, hogHistCenter);
            
            // need to create the HOG here to calculate the dominant orientations
            //   for multiple versions of CRegion if needed.
            TIntList dos = calculateDominantOrientations(hogHistCenter);
            if (false) {
                // use the MSER ellipse orientation and the top 80% HOG bins
                regionPointsScaled.hogOrientations.addAll(dos);
            } else {
                // only use hog orientation
                regionPointsScaled.hogOrientations.clear();
                regionPointsScaled.hogOrientations.add(dos.get(0));
            }
            
            List<CRegion> crsScaled = cn.canonicalizeRegions(
                w, h, regionPointsScaled);
            
            assert(!crsScaled.isEmpty());
          
            //DEBUG
            /*{
                Image tmp = gsImgScaled.copyToColorGreyscale();
                for (PairInt p : regionPointsScaled.points) {
                    ImageIOHelper.addPointToImage(
                        p.getX(), p.getY(), 
                        tmp,
                        0, 255, 0, 0);
                };
                MiscDebug.writeImage(tmp, "_DBG00_" + scale + "_" + ts); 
                int z = 0;
            }*/
            
            //System.out.println("created " + crs.size() + " from RegionPoints " + i);
            
            count++;
            
            for (Canonicalizer.CRegion crScaled : crsScaled) {

                int rIdx = cRegionMapReference.size();
                
                if (crScaled.getOffsetsToOrigCoords().size() < 9) {
                    continue;
                }
                
                crScaled.dataIdx = rIdx;
            
                hogMgs.addARegion(
                    gsImg2, ptImg2, gXY2, theta2, 
                    crScaled, minMaxXY, refFramePixs);
                
                assert(regionIndexRegions.containsKey(rIdx));
                assert(regionIndexMinMaxXY.containsKey(rIdx));
                assert(regionCoords.containsKey(rIdx));
                assert(cRegionMapReference.get(rIdx).equals(crScaled));
                assert(!histHOGMap.isEmpty());
                assert(!histHCPTMap.isEmpty());
                assert(!histHGSMap.isEmpty());
            }
        }
        
        assert(!histHOGMap.isEmpty());
        assert(!histHCPTMap.isEmpty());
        assert(!histHGSMap.isEmpty());
        assert(regionCoords.size() >= count);
        assert(regionIndexRegions.size() >= count);
        assert(regionIndexMinMaxXY.size() >= count);
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

    TIntList calculateDominantOrientations(int[] hogHist) {
        
        int maxIdx = MiscMath.findYMaxIndex(hogHist);

        if (maxIdx == -1) {
            throw new IllegalArgumentException("histogram is full of "
                + " min value longs");
        }        

        TIntSet set = new TIntHashSet();
        set.add(maxIdx);
        
        // if any bins have values within 80% of max, add to maxIdxs
        TIntList list = new TIntArrayList();
        list.add(maxIdx);

        long max = hogHist[maxIdx];
        double limit = 0.8 * max;

        for (int i = 0; i < hogHist.length; ++i) {
            long v = hogHist[i];
            if (v >= limit) {
                if (set.add(i)) {
                    list.add(i);
                }
            }
        }

        int binWidth = 180 / nAngleBins;
        
        for (int i = 0; i < list.size(); ++i) {
            int idx = list.get(i);
            int angle = Math.round((idx + 0.5f) * binWidth);
            list.set(i, angle);
        }
        
        return list;
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
