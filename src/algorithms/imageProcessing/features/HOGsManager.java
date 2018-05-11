package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TLongHashSet;
import java.util.Arrays;

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
     * @param outHist
     * @return true if method succeeded, else false.  can return false if the
     * addARegion failed due to having fewer than 9 pixels in the CRegion
     * for rIndex.
     */
    public boolean extractBlockHOG(int x, int y, int[] outHist) {
        return extractBlock(TYPE.HOG, x, y, outHist);
    }
    
    /**
     * extract the block surrounding the feature.
     * the number of pixels in a cell and the number of cells in block were set during
     * construction.
     * @param x
     * @param y
     * @param outHist
     * @return 
     */
    public boolean extractBlockHCPT(int x, int y, int[] outHist) {
        return extractBlock(TYPE.HCPT, x, y, outHist);
    }
    
    /**
     * extract the block surrounding the feature.
     * the number of pixels in a cell and the number of cells in block were set during
     * construction.
     * 
     * @param x
     * @param y
     * @param outHist
     * @return 
     */
    public boolean extractBlockHGS(int x, int y, int[] outHist) {
        return extractBlock(TYPE.HGS, x, y, outHist);
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
    private boolean extractBlock(TYPE type, int x, int y, int[] outHist) {

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
        
        return extractBlock(hist, x, y, outHist);
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
    private boolean extractBlock(int[][] hist, int x, int y, int[] outHist) {

        if (x < 0 || y < 0 || x >= w || y >= h) {
            throw new IllegalArgumentException("x or y is out of bounds of "
                + "original image");
        }
        
        // uses the block normalization recomended by Dalal & Triggs,
        //   the summary of histogram counts over all cells
        //   is used to normalize each cell by that sum.

        Arrays.fill(outHist, 0, outHist.length, 0);
                
        int r = N_CELLS_PER_BLOCK_DIM >> 1;
        int stopY = y + r;
        int stopX = x + r;
        int startX = x - r;
        int startY = y - r;
        if ((h & 1) == 0) {
            startX--;
            startY--;            
        }
        if (startX < 0) {
            startX = 0;
        }
        if (startY > 0) {
            startY = 0;
        }
        if (stopX >= w) {
            stopX = w - 1;
        }
        if (stopY >= h) {
            stopY = h - 1;
        }
        
        
        int[] outputN = new int[1];  
        
        HOGUtil.extractWindow(hist, startX, stopX, startY, stopY, w, h, 
            outHist, outputN);
        
        double blockTotal = HOGUtil.sumCounts(outHist);
        blockTotal *= blockTotal;

        double norm;
        if (blockTotal > 0) {
            blockTotal /= (double)outputN[0];
            blockTotal = Math.sqrt(blockTotal);
            norm = 255./blockTotal;
        } else {
            norm = 255.;
        }
          
        for (int i = 0; i < outHist.length; ++i) {
            outHist[i] = (int)Math.round(norm * outHist[i]);
        }
        
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
