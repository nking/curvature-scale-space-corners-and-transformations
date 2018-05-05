package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TLongHashSet;
import java.util.Set;

/**
 * similar to HOGRegionsManager except that there is only one histogram of each
 * type for all regions.
 * 
 * @author nichole
 */
public class HOGsManager {
    
    private static float eps = 0.000001f;
    
    // 9 is default
    private final int nAngleBins;
    
    private final int nHistBins;

    // 6 x 6 is recommended
    private final int N_PIX_PER_CELL_DIM;

    // 2x2 or 3x3 is recommended
    private final int N_CELLS_PER_BLOCK_DIM;
    
    // histogram integral images with a windowed sum of N_PIX_PER_CELL_DIM
    private final int[][] histHOG2D;
    private final int[][] histHCPT2D;
    private final int[][] histHGS2D;

    private static int maskValue = 0;

    private final int w;
    private final int h;

    private boolean debug = false;
    
    public HOGsManager(GreyscaleImage gsImg, GreyscaleImage ptImg,
        TIntObjectMap<Canonicalizer.CRegion> regionMap, 
        int nCellsPerDim, int nPixPerCellDim, 
        int nAngleBins, int nHCPTHGSBins) {

        this.nAngleBins = nAngleBins;
        nHistBins = nHCPTHGSBins;
        N_PIX_PER_CELL_DIM = nPixPerCellDim;
        N_CELLS_PER_BLOCK_DIM = nCellsPerDim;
        w = gsImg.getWidth();
        h = gsImg.getHeight();

        ImageProcessor imageProcessor = new ImageProcessor();
        
        GreyscaleImage[] gXgY = 
            imageProcessor.createCentralDifferenceGradients(gsImg);
        GreyscaleImage theta = imageProcessor.computeTheta180(gXgY[0], gXgY[1]);
        GreyscaleImage gXY = 
            imageProcessor.combineConvolvedImages(gXgY[0], gXgY[1]);
        
        // non-region pixels are excluded because magnitude is zero
        histHOG2D = HOGUtil.createHOGHistogram(gXY, theta, 
            nAngleBins, N_PIX_PER_CELL_DIM);
        
        TLongSet regionPixelCoords = allPixelCoords(regionMap);
        
        histHCPT2D = HOGUtil.createHCPTHistogram(ptImg, regionPixelCoords, 
            nHistBins, N_PIX_PER_CELL_DIM);
        
        histHGS2D = HOGUtil.createHGSHistogram(gsImg, regionPixelCoords, 
            nHistBins, N_PIX_PER_CELL_DIM);
    }

    public HOGsManager(GreyscaleImage gsImg, GreyscaleImage ptImg,
        GreyscaleImage gradient, GreyscaleImage theta,
        TIntObjectMap<Canonicalizer.CRegion> regionMap, 
        int nCellsPerDim, int nPixPerCellDim, 
        int nAngleBins, int nHCPTHGSBins) {

        this.nAngleBins = nAngleBins;
        nHistBins = nHCPTHGSBins;
        N_PIX_PER_CELL_DIM = nPixPerCellDim;
        N_CELLS_PER_BLOCK_DIM = nCellsPerDim;
        w = gsImg.getWidth();
        h = gsImg.getHeight();

        ImageProcessor imageProcessor = new ImageProcessor();
        
        // non-region pixels are excluded because magnitude is zero
        histHOG2D = HOGUtil.createHOGHistogram(gradient, theta, 
            nAngleBins, N_PIX_PER_CELL_DIM);
        
        TLongSet regionPixelCoords = allPixelCoords(regionMap);
        
        histHCPT2D = HOGUtil.createHCPTHistogram(ptImg, regionPixelCoords, 
            nHistBins, N_PIX_PER_CELL_DIM);
        
        histHGS2D = HOGUtil.createHGSHistogram(gsImg, regionPixelCoords, 
            nHistBins, N_PIX_PER_CELL_DIM);
    }
    
    public void setToDebug() {
        debug = true;
    }
    
    int[][] createHCPTHistogram(GreyscaleImage ptImg, TLongSet regionPixelCoords) {
        return HOGUtil.createHCPTHistogram(ptImg, regionPixelCoords, nHistBins, 
            N_PIX_PER_CELL_DIM);
    }
    
    private TLongSet allPixelCoords(TIntObjectMap<CRegion> regionMap) {

        TLongSet pixs = new TLongHashSet();
        
        TIntObjectIterator<CRegion> iter = regionMap.iterator();
        for (int i = 0; i < regionMap.size(); ++i) {
            iter.advance();
            CRegion cRegion = iter.value();
            pixs.addAll(cRegion.getPixelCoords());
        }
        
        return pixs;
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
     * @param x
     * @param y
     * @param minMaxXY
     * @param rCoords
     * @param outHist
     * @return true if method succeeded, else false.  can return false if the
     * addARegion failed due to having fewer than 9 pixels in the CRegion
     * for rIndex.
     */
    public boolean extractBlockHOG(int x, int y, int[] minMaxXY, 
        TLongSet rCoords, int[] outHist) {
        return extractBlock(TYPE.HOG, x, y, minMaxXY, rCoords, outHist);
    }
    
    /**
     * extract the block surrounding the feature.
     * the number of pixels in a cell and the number of cells in block were set during
     * construction.
     * @param x
     * @param y
     * @param minMaxXY
     * @param rCoords
     * @param outHist
     * @return 
     */
    public boolean extractBlockHCPT(int x, int y, int[] minMaxXY, 
        TLongSet rCoords, int[] outHist) {
        return extractBlock(TYPE.HCPT, x, y, minMaxXY, rCoords, outHist);
    }
    
    /**
     * extract the block surrounding the feature.
     * the number of pixels in a cell and the number of cells in block were set during
     * construction.
     * 
     * @param x
     * @param y
     * @param minMaxXY
     * @param rCoords
     * @param outHist
     * @return 
     */
    public boolean extractBlockHGS(int x, int y, int[] minMaxXY, 
        TLongSet rCoords, int[] outHist) {
        return extractBlock(TYPE.HGS, x, y, minMaxXY, rCoords, outHist);
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
    private boolean extractBlock(TYPE type, int x, int y, 
        int[] minMaxXY, TLongSet rCoords, int[] outHist) {

        if ((type.equals(TYPE.HOG) && outHist.length != nAngleBins) ||
            (!type.equals(TYPE.HOG) && outHist.length != nHistBins)) {
            throw new IllegalArgumentException("outHist.length != expected");
        }

        if (x < 0 || y < 0 || x >= w || y >= h) {
            throw new IllegalArgumentException("x or y is out of bounds of "
                + "original image");
        }
        
        // uses the block normalization recomended by Dalal & Triggs,
        //   the summary of histogram counts over all cells
        //   is used to normalize each cell by that sum.

        int[][] hist = null;
        switch (type) {
            case HOG: hist = histHOG2D; break;
            case HCPT: hist = histHCPT2D; break;
            case HGS: hist = histHGS2D; break;
            default: break;
        }
        
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
    private boolean extractBlock(int[][] hist, TLongSet pixelIndexes, int[] minMaxXY,
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
                
                if (!pixelIndexes.contains(pixIdx)) {
                //    continue;
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

    public int getImageWidth() {
        return w;
    }
    
    public int getImageHeight() {
        return h;
    }

    public int getNumberOfBins() {
        return nAngleBins;
    }

}
