package algorithms.imageProcessing;

import algorithms.QuickSort;
import algorithms.SubsetChooser;
import algorithms.compGeometry.ParabolaLeastSquares;
import algorithms.imageProcessing.Sky.SkyObject;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.RegionGeometry;
import algorithms.imageProcessing.features.mser.MSEREdges;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import algorithms.util.VeryLongBitString;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

/**
 * NOT READY FOR USE.
 * rewriting most of it.
 * 
 * class with methods to find a rainbow within an image and to create a hull
 * to encapsulate it for various methods.
 * 
 * @author nichole
 */
public class RainbowFinder {
    
    private Logger log = Logger.getLogger(this.getClass().getName());

    static class Hull {
        float[] xHull;
        float[] yHull;
    }
    
    private Hull rainbowHull = null;
    
    public RainbowFinder() {
    }
    
    /**
     * find rainbows in the blob detector results of MSER.
     * NOTE that the method needs more testing.
     * 
     * @param mserEdges
     * @return 
     */
    public List<SkyObject> findRainbows(MSEREdges mserEdges) {
        
        // there may be more than one rainbow in this set.
        // for example, a test image has 1 strong rainbow and 2 fainter ones.
            
        ImageExt img = mserEdges.getClrImg();
        
        List<Region> polarThetaPositive = mserEdges.getOrigGsPtRegions().get(2);
        
        /*
        segments of rainbows in two test images are found as elongated
        thin ellipses in the polar theta positive images,
        and those then have adjacent thin elongated regions in
        the polar theta negative images where the brighter arcs of the
        rainbow are.
        
        NOTE: once a rainbow is found, could look for others fainter in the
        image.
        */
        List<Set<PairInt>> listOfSets0 = new ArrayList<Set<PairInt>>();
        
        List<OneDIntArray> hists0 = new ArrayList<OneDIntArray>();
        
        List<RegionGeometry> rgs0 = new ArrayList<RegionGeometry>();

        //mserEdges._debugOrigRegions(2, "_PT_");
        findPositivePT(img, polarThetaPositive, listOfSets0, hists0, rgs0);
       
        List<Region> polarThetaNegative = mserEdges.getOrigGsPtRegions().get(3);
        
        List<Set<PairInt>> listOfSets1 = new ArrayList<Set<PairInt>>();
        
        List<OneDIntArray> hists1 = new ArrayList<OneDIntArray>();
        
        List<RegionGeometry> rgs1 = new ArrayList<RegionGeometry>();

        findNegativePT(img, polarThetaNegative, listOfSets1, hists1, rgs1);
       
        List<VeryLongBitString> listOfSetBits0 = makeBitStrings(listOfSets0, img);
        
        List<SkyObject> output = new ArrayList<SkyObject>();
        
        TIntList arcIdxs = findLargeArc(listOfSetBits0, listOfSets0, hists0, 
            rgs0, img);
        
        if (arcIdxs == null || arcIdxs.isEmpty()) {
            return null;
        }
            
        // if have found a large arc in the positive image, search for 
        //    adjacent arcs in the negative regions

        Set<PairInt> arcPoints = new HashSet<PairInt>();
        for (int j = 0; j < arcIdxs.size(); ++j) {
            arcPoints.addAll(listOfSets0.get(arcIdxs.get(j)));
        }

        NearestNeighbor2D nn = new NearestNeighbor2D(arcPoints, 
            img.getWidth(), img.getHeight());

        int[] negativeIdxs = findAdjacent(nn, listOfSets1, rgs1);

        if (negativeIdxs != null) {
            for (int idx1 : negativeIdxs) {
                arcPoints.addAll(listOfSets1.get(idx1));
            }

            Arrays.sort(negativeIdxs);
            for (int j = (negativeIdxs.length - 1); j > -1; --j) {
                int idx = negativeIdxs[j];
                listOfSets1.remove(idx);
                rgs1.remove(idx);
                hists1.remove(idx);
            }                
        }

        /*
        ptImg values for histogram bins:
         0:  red = 0 - 18
         1:  orange = 18 - 40
         2:  yellow = 41 - 60ish
         3:  green = 61 - 106
         4:  blue = 107 - 192
         5:  purple = 193 - 255
        */
        int[] ptCH = ColorHistogram.createPTHistogram(mserEdges.getPtImg(), 
            arcPoints);
        
        float[] normalizedHist = new float[ptCH.length];
        int tot = 0;
        for (int c : ptCH) {
            tot += c;
        }
        for (int j = 0; j < ptCH.length; ++j) {
            normalizedHist[j] = (float)ptCH[j]/(float)tot;
        }
        
        System.out.println("rainbow? " + Arrays.toString(ptCH) + 
            "\n   " + Arrays.toString(normalizedHist));
        
        boolean dark = normalizedHist[0] > 0.1 
            && normalizedHist[1] > 0.1;
        boolean bright = normalizedHist[1] > 0.01
            && normalizedHist[3] > 0.01 
            && normalizedHist[4] > 0.01;
        
        if (!dark && !bright) {
            return null;
        }
        
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        double[] xyCenter = ch.calculateXYCentroids(arcPoints);
        int x = (int)Math.round(xyCenter[0]);
        int y = (int)Math.round(xyCenter[1]);

        SkyObject obj = new SkyObject();
        obj.points = arcPoints;
        obj.xyCenter = new int[]{x, y};
        output.add(obj);

        arcIdxs.sort();
        for (int j = (arcIdxs.size() - 1); j > -1; --j) {
            int arcIdx = arcIdxs.get(j);
            //System.out.println("removing " + 
            //    rgs0.get(arcIdx).xC + ", " + rgs0.get(arcIdx).yC);
            listOfSets0.remove(arcIdx);
            listOfSetBits0.remove(arcIdx);
            hists0.remove(arcIdx);
            rgs0.remove(arcIdx);
        }

        return output;
    }
    
    /*
        rainbowCoeff = findRainbowPoints(skyPoints, 
            removedSets.getReflectedSunRemoved(), colorImg, 
            xOffset, yOffset, imageWidth, imageHeight,
            skyIsDarkGrey, allSkyColor, outputRainbowPoints);

        if (outputRainbowPoints.size() < 10) {
            outputRainbowPoints.clear();
            rainbowCoeff = null;
        }
        
        if (rainbowCoeff != null) {
            
            rainbowHull = createRainbowHull(rainbowCoeff, 
                outputRainbowPoints, colorImg, xOffset, yOffset);
            
            if (rainbowHull != null) {
                
                //TODO: may need adjustment for a boundary being an image boundary
                
                int minXHull = (int)MiscMath.findMin(rainbowHull.xHull);
                int maxXHull = (int)Math.ceil(MiscMath.findMax(rainbowHull.xHull));
                int minYHull = (int)MiscMath.findMin(rainbowHull.yHull);
                int maxYHull = (int)Math.ceil(MiscMath.findMax(rainbowHull.yHull));
                PointInPolygon p = new PointInPolygon();
                for (int col = minXHull; col <= maxXHull; col++) {
                    for (int row = minYHull; row <= maxYHull; row++) {
                        boolean in = p.isInSimpleCurve(col, row, 
                            rainbowHull.xHull, rainbowHull.yHull, 
                            rainbowHull.yHull.length);
                        if (in) {
                            excludePointsInRainbowHull.add(new PairInt(col, row));
                        }
                    }
                }
                
                if (!excludePointsInRainbowHull.isEmpty()) {
                    // addRainbow to Hull, but only if there are sky points adjacent to hull
                    addRainbowToPoints(skyPoints, imageWidth - 1, imageHeight - 1);
                }
            }
        }
    }
    
    public void addRainbowToSkyPoints(Set<PairInt> skyPoints, 
        int lastImgCol, int lastImgRow) {
        
        addRainbowToPoints(skyPoints, lastImgCol, lastImgRow);
    }
    
    private void addRainbowToPoints(Set<PairInt> skyPoints, 
        int lastImgCol, int lastImgRow) {
         
        
        int nHull = rainbowHull.xHull.length;
        int nHalf = nHull >> 1;
        for (int c = 0; c < nHalf; c++) {
            
            int count0 = c;
            int count1 = nHull - 2 - c;
            
            double dy0 = rainbowHull.yHull[count0 + 1] - rainbowHull.yHull[count0];
            double dx0 = rainbowHull.xHull[count0 + 1] - rainbowHull.xHull[count0];
            int dist0 = (int)Math.sqrt(dx0*dx0 + dy0*dy0);

            double dy1 = rainbowHull.yHull[count1 - 1] - rainbowHull.yHull[count1];
            double dx1 = rainbowHull.xHull[count1 - 1] - rainbowHull.xHull[count1];
            int dist1 = (int)Math.sqrt(dx1*dx1 + dy1*dy1);

            boolean removeSection = false;
            
            int dist = (dist0 < dist1) ? dist0 : dist1;
            for (int i = 0; i < dist; i++) {
                int x0 = (int)(rainbowHull.xHull[count0] + i*dx0);
                int y0 = (int)(rainbowHull.yHull[count0] + i*dy0);

                int x1 = (int)(rainbowHull.xHull[count1] + i*dx1);
                int y1 = (int)(rainbowHull.yHull[count1] + i*dy1);

                int n0Sky = 0;
                int n0SkyPossible = 0;
                int n1Sky = 0;
                int n1SkyPossible = 0;
                for (int type = 0; type < 2; type++) {
                    int x = (type == 0) ? x0 : x1;
                    int y = (type == 0) ? y0 : y1;
                    int n = 0;
                    int nPossible = 0;
                    for (int col = (x - 1); col <= (x + 1); col++) {
                        if ((col < 0) || (col > lastImgCol)) {
                            continue;
                        }
                        for (int row = (y - 1); row <= (y + 1); row++) {
                            if ((row < 0) || (row > lastImgRow)) {
                                continue;
                            }
                            PairInt p = new PairInt(col, row);
                            if (!excludePointsInRainbowHull.contains(p)) {
                                nPossible++;
                                if (skyPoints.contains(p)) {
                                    n++;
                                }
                            }
                        }
                    }
                    if (type == 0) {
                        n0Sky = n;
                        n0SkyPossible = nPossible;
                    } else {
                        n1Sky = n;
                        n1SkyPossible = nPossible;
                    }
                }
                // evaluate the n's and shorten rainbowHull values for count0 and count1 if needed
                float n0Div = (float)n0Sky/(float)n0SkyPossible;
                float n1Div = (float)n1Sky/(float)n1SkyPossible;
                if ((n0Div < 0.5) || (n1Div < 0.5)) {
                    removeSection = true;
                    break;
                }
            } // end for i
            if (removeSection) {
                
                // remove points from rainbowHull
                float[] xh = new float[]{
                    rainbowHull.xHull[count0], rainbowHull.xHull[count0 + 1],
                    rainbowHull.xHull[count1 - 1], rainbowHull.xHull[count1],
                    rainbowHull.xHull[count0]};
                float[] yh = new float[]{
                    rainbowHull.yHull[count0], rainbowHull.yHull[count0 + 1],
                    rainbowHull.yHull[count1 - 1], rainbowHull.yHull[count1],
                    rainbowHull.yHull[count0]};
                
                int minXHull = (int)MiscMath.findMin(xh);
                int maxXHull = (int)Math.ceil(MiscMath.findMax(xh));
                int minYHull = (int)MiscMath.findMin(yh);
                int maxYHull = (int)Math.ceil(MiscMath.findMax(yh));
                PointInPolygon p = new PointInPolygon();
                for (int col = minXHull; col <= maxXHull; col++) {
                    for (int row = minYHull; row <= maxYHull; row++) {
              
                        boolean in = p.isInSimpleCurve(col, row, xh, yh, 
                            xh.length);
                        if (in) {
                            excludePointsInRainbowHull.remove(new PairInt(col, row));
                        }
                    }
                }
            }
        }
        
        skyPoints.addAll(excludePointsInRainbowHull);
        
    }

     // @return polynomial fit coefficients to 
    // y[i] = c0*1 + c1*x[i] + c2*x[i]*x[i].  this may be null if a fit wasn't
     // possible.
     
    float[] findRainbowPoints(Set<PairInt> skyPoints, 
        Set<PairInt> reflectedSunRemoved,
        ImageExt colorImg, int xOffset, int yOffset, 
        int imageWidth, int imageHeight,
        boolean skyIsDarkGrey, GroupPixelColors allSkyColor,
        Set<PairInt> outputRainbowPoints) {

        Set<PairInt> rainbowPoints = findRainbowColoredPoints(colorImg, 
            reflectedSunRemoved, xOffset, yOffset, skyIsDarkGrey);

        if (rainbowPoints.isEmpty()) {
            return null;
        }
        
        if (rainbowPoints.size() < 12) {
            return null;
        }
        
        if (rainbowPoints.size() > 0.25*(colorImg.getWidth()*colorImg.getHeight())) {
            if ((allSkyColor.getAvgBlue()/allSkyColor.getAvgRed()) < 0.5) {
                return null;
            }
        }
        
        // fit a polynomial to rainbow points.  
        // would prefer a circle, but the optical depth of the dispersers and the
        // orientation of groups of them is not always a slab perpendicular to 
        // the camera

        int[] minMaxXY = MiscMath.findMinMaxXY(rainbowPoints);
        log.fine("rainbow range in x: " + minMaxXY[0] + " to " + minMaxXY[1]);
        
        //TODO: consider contiguous subsets at this point
        
        PolynomialFitter polyFitter = new PolynomialFitter();
        //y = c0*1 + c1*x[i] + c2*x[i]*x[i]
        float[] coef = polyFitter.solveAfterRandomSampling(rainbowPoints);
             
        if (coef == null) {
            return null;
        }
        
        log.info("rainbow polynomial coefficients = " + Arrays.toString(coef));
        log.info("image dimensions are " + colorImg.getWidth() + " X " + 
            colorImg.getHeight() + " pixels^2");
         
        polyFitter.plotFit(coef, rainbowPoints, colorImg.getWidth(),
            colorImg.getHeight(), 23, "rainbow points");
        
        double resid = polyFitter.calcResiduals(coef, rainbowPoints);

        //TODO: determine this more accurately:
        if (resid > 5) {
            
            Set<PairInt> bestFittingPoints = new HashSet<PairInt>();
            
            coef = polyFitter.solveForBestFittingContiguousSubSets(
                rainbowPoints, bestFittingPoints, colorImg.getWidth(), 
                colorImg.getHeight());
            
        }
    }
   
    private Hull createRainbowHull(float[] rainbowCoeff, 
        Set<PairInt> rainbowPoints, ImageExt originalColorImage, int xOffset, 
        int yOffset) {   
        
        int width = originalColorImage.getWidth() - xOffset;
        int height = originalColorImage.getHeight() - yOffset;
        
        rainbowSkeletonX = new float[10];
        rainbowSkeletonY = new float[10];
        generatePolynomialPoints(rainbowPoints, rainbowCoeff, rainbowSkeletonX, 
            rainbowSkeletonY);
       
        float maxOfPointMinDistances = maxOfPointMinDistances(rainbowPoints,
            rainbowSkeletonX, rainbowSkeletonY);
        
        float high = 2 * maxOfPointMinDistances;
        float low = maxOfPointMinDistances / 2;
        int nMatched = 0;
     
        float[] xPoly = new float[2 * rainbowSkeletonX.length + 1];
        float[] yPoly = new float[xPoly.length];
        int nMaxIter = 5;
        int nIter = 0;
        
        int eps = (int)(0.15f * rainbowPoints.size());
        
        while ((low < high) && (nIter < nMaxIter)) {
            
            float mid = (high + low)/2.f;
            
            hullHalfWidth = mid;
            
            populatePolygon(rainbowSkeletonX, rainbowSkeletonY, mid, 
                xPoly, yPoly, rainbowCoeff, width, height);
            
            nMatched = nPointsInPolygon(rainbowPoints, xPoly, yPoly);
            
            log.info("low=" + low + " high=" + high + " mid=" + mid 
                + " nMatched=" + nMatched + " out of " + rainbowPoints.size());

            if (Math.abs(nMatched - rainbowPoints.size()) < eps) {
                if (low < mid) {
                    // decrease high so next mid is lower
                    high = mid;
                } else {
                    break;
                }
            } else {
                // nMatched < rainbowPoints.size()
                // increase low so next mid is higher
                low = mid;
            }
            
            nIter++;
        }
        
        //TODO: once a reasonable hullHalfWidth has been determined,
        // may want to regenerate the hull with higher resolution,
        // that is 5 or 10 times the number of skeleton points on 
        // each side
         
        Hull hull = new Hull();
        hull.xHull = xPoly;
        hull.yHull = yPoly;

        return hull;
    }
    
    protected void generatePolynomialPoints(Set<PairInt> points, 
        float[] polyCoeff, float[] outputX, float[] outputY) {
       
        PolynomialFitter fitter = new PolynomialFitter();
        float[] minXYMaxXY = fitter.determineGoodEndPoints(polyCoeff, points);
        
        // find the furthest points that have the smallest residuals on each side
        int minX = (int)minXYMaxXY[0];
        int yForMinX = (int)minXYMaxXY[1];
        int maxX = (int)minXYMaxXY[2];
        int yForMaxX = (int)minXYMaxXY[3];
        
        log.info("polyCoeff=" + Arrays.toString(polyCoeff) 
            + " endpoints=(" + minX + "," + yForMinX + ") (" + maxX + "," 
            + yForMaxX + ")");
       
        int n = outputX.length;
        
        // max-min divided by 9 gives 8 more points
        float deltaX = (maxX - minX)/(float)(n - 1);
        
        outputX[0] = minX;
        outputY[0] = yForMinX;
        for (int i = 1; i < (n - 1); i++) {
            outputX[i] = outputX[i - 1] + deltaX;
            outputY[i] = polyCoeff[0] + polyCoeff[1] * outputX[i] 
                + polyCoeff[2] * outputX[i] * outputX[i];
            
        }
        
        outputX[n - 1] = maxX;
        outputY[n - 1] = yForMaxX;
        
    }
    
    protected float maxOfPointMinDistances(Set<PairInt> rainbowPoints, 
        float[] xc, float[] yc) {
        
        double maxDistSq = Double.MIN_VALUE;
        
        for (PairInt p : rainbowPoints) {
            int x = p.getX();
            int y = p.getY();
            double minDistSq = Double.MAX_VALUE;
            int minIdx = -1;
            for (int i = 0; i < xc.length; i++) {
                float diffX = xc[i] - x;
                float diffY = yc[i] - y;
                float dist = (diffX * diffX) + (diffY * diffY);
                if (dist < minDistSq) {
                    minDistSq = dist;
                    minIdx = i;
                }
            }
            if (minDistSq > maxDistSq) {
                maxDistSq = minDistSq;
            } 
        }
        
        return (float)Math.sqrt(maxDistSq);
    }

    protected void populatePolygon(float[] x, float[] y, float dist, 
        float[] outputXPoly, float[] outputYPoly, float[] polynomialCoeff,
        int imgWidth, int imgHeight) {
        
        if (x == null || y == null) {
            throw new IllegalArgumentException("neither x nor y can be null");
        }
        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y must be the same length");
        }
        if (outputXPoly == null || outputYPoly == null) {
            throw new IllegalArgumentException(
                "neither outputXPoly nor outputYPoly can be null");
        }
        if (outputXPoly.length != outputYPoly.length) {
            throw new IllegalArgumentException(
                "outputXPoly and outputYPoly must be the same length");
        }
        if (polynomialCoeff == null) {
            throw new IllegalArgumentException(
                "polynomialCoeff cannot be null");
        }
        if (polynomialCoeff.length != 3) {
            throw new IllegalArgumentException(
                "polynomialCoeff.length has to be 3");
        }
        if (outputXPoly.length != (2*x.length + 1)) {
            throw new IllegalArgumentException("outputXPoly.length must be " +
                " (2 * x.length) + 1");
        }
        
        //y = c0*1 + c1*x[i] + c2*x[i]*x[i]
        
        //dy/dx = c1 + 2*c2*x[i] = tan theta
        
        int n = outputXPoly.length;
        
       
        int count0 = 0;
        int count1 = n - 2;
        
        for (int i = 0; i < x.length; i++) {
            
            double dydx = polynomialCoeff[1] + (2. * polynomialCoeff[2] * x[i]);
            
            if (dydx == 0) {
                // same x, y's are +- dist
                outputXPoly[count0] = x[i];
                outputYPoly[count0] = y[i] + dist;
                count0++;
                outputXPoly[count1] = x[i];
                outputYPoly[count1] = y[i] - dist;
                count1--;
                continue;
            }
            
            double tangentSlope = -1./dydx;
            double theta = Math.atan(tangentSlope);
            double dy = dist * Math.sin(theta);    
            double dx = dist * Math.cos(theta);
            
            if ((count0 == 0) && (x[i] == 0)) {
                dy = dist;
                dx = 0;
            } else if (
                ((count0 == ((x.length/2) - 1) || (count0 == (x.length/2))))
                && (x[i] == (imgWidth - 1))) {
                dy = dist;
                dx = 0;
            }
            
            //System.out.println("i=" + i + " theta=" + theta
            //    + " x[i]=" + x[i] + " dx=" + dx + " dy=" + dy);
          
            float xHigh = (float)(x[i] + dx);
            float yHigh = (float)(y[i] + dy);
            
            float xLow = (float)(x[i] - dx);
            float yLow = (float)(y[i] - dy);
              
            if (theta < 0) {
                float tx = xHigh;
                float ty = yHigh;
                xHigh = xLow;
                yHigh = yLow;
                xLow = tx;
                yLow = ty;
            }
            
            if (xHigh < 0) {
                xHigh = 0;
            }
            if (yHigh < 0) {
                yHigh = 0;
            }
            if (xLow < 0) {
                xLow = 0;
            }
            if (yLow < 0) {
                yLow = 0;
            }
            if (xLow > (imgWidth - 1)) {
                xLow = (imgWidth - 1);
            }
            if (xHigh > (imgWidth - 1)) {
                xHigh = (imgWidth - 1);
            }
            if (yLow > (imgHeight - 1)) {
                yLow = (imgHeight - 1);
            }
            if (yHigh > (imgHeight - 1)) {
                yHigh = (imgHeight - 1);
            }
            
            outputXPoly[count0] = xHigh;
            outputYPoly[count0] = yHigh;
            count0++;
            outputXPoly[count1] = xLow;
            outputYPoly[count1] = yLow;
            count1--;
        }
       
        outputXPoly[outputXPoly.length - 1] = outputXPoly[0];
        outputYPoly[outputXPoly.length - 1] = outputYPoly[0];
    }
    */
    
    private void extractHues(ImageExt img, Set<PairInt> points, 
        int[] hues, float[] sbAvg) {
        
        float sum0 = 0;
        float sum1 = 0;
        
        float[] hsb = new float[3];
        for (PairInt p : points) {
            int r = img.getR(p);
            int g = img.getG(p);
            int b = img.getB(p);
            Color.RGBtoHSB(r, g, b, hsb);
            
            sum0 += hsb[1];
            sum1 += hsb[2];
            
            // TODO: consider satudation and brightness limits
            int bin = (int)(hsb[0]/0.1f);
            if (bin == 10) {
                bin--;
            }
            hues[bin]++;
        }
        
        sum0 /= (float)points.size();
        sum1 /= (float)points.size();
    
        sbAvg[0] = sum0;
        sbAvg[1] = sum1;
    }

    private boolean isWithinBounds(int x, int y, int width, int height) {
    
        if (x < 0 || y < 0 || (x >= width) || (y >= height)) {
            return false;
        }
        return true;
    }
    
    private boolean isEmpty(int[] hues) {
        for (int h : hues) {
            if (h != 0) {
                return false;
            }
        }
        return true;
    }

    private List<VeryLongBitString> clusterByIntersection(List<OneDIntArray> hists, 
        float minIntersection) {
        
        ColorHistogram ch = new ColorHistogram();
        
        List<VeryLongBitString> out = new ArrayList<VeryLongBitString>();
        
        for (int i = 0; i < hists.size(); ++i) {
            
            OneDIntArray h = hists.get(i);
            
            VeryLongBitString indexes = new VeryLongBitString(hists.size());
            boolean didSet = false;
            for (int j = (i + 1); j < hists.size(); ++j) {
                OneDIntArray h2 = hists.get(j);
                float intersection = ch.intersection(h.a, h2.a);
                if (intersection >= minIntersection) {
                    indexes.setBit(j);
                    didSet = true;
                }
            }
            if (didSet) {
                indexes.setBit(i);
                // if this bitstring is not a subset of an existing in out, add
                boolean doNotAdd = false;
                for (int j = 0; j < i; ++j) {
                    VeryLongBitString bs = out.get(j);
                    VeryLongBitString intersectionBS = bs.and(indexes);
                    if (intersectionBS.getNSetBits() == indexes.getNSetBits()) {
                        doNotAdd = true;
                        break;
                    }
                }
                if (!doNotAdd) {
                    out.add(indexes);
                }
            }
        }
        
        // add entries for the sets by themselves
        for (int i = 0; i < hists.size(); ++i) {            
            OneDIntArray h = hists.get(i);
            VeryLongBitString indexes = new VeryLongBitString(hists.size());
            indexes.setBit(i);
            out.add(indexes);
        }
       
        return out;
    }

    private void findPositivePT(ImageExt img, 
        List<Region> polarThetaPositive, List<Set<PairInt>> listOfSets, 
        List<OneDIntArray> hists, List<RegionGeometry> rgs) {
        
        for (int rIdx = 0; rIdx < polarThetaPositive.size(); ++rIdx) {
            
            Region r = polarThetaPositive.get(rIdx);
           
            Set<PairInt> set1 = null;
                        
            for (int i = 0; i < r.accX.size(); ++i) {
                int x = r.accX.get(i);
                int y = r.accY.get(i);
                PairInt p = new PairInt(x, y);
                int pixIdx = img.getInternalIndex(p);
               
                if (set1 == null) {
                    set1 = new HashSet<PairInt>();
                }
                set1.add(p);
            }
            
            if (set1 != null) {
                
                RegionGeometry rg = Canonicalizer.calculateEllipseParams(
                    r, img.getWidth(), img.getHeight());
                
                /*System.out.format(
                    "%d : xy=(%d,%d) angle=%.2f "
                        + " ecc=%.3f minor=%.3f major=%.3f n=%d\n",
                    rIdx, rg.xC, rg.yC,
                    (rg.orientation * 180. / Math.PI),
                    (float) rg.eccentricity, (float) rg.minor, 
                    (float) rg.major, set1.size());
                */
                
                if (rg.eccentricity >= 0.95) {
                    
                    int[] hues = new int[10];
                    float[] sbAvg = new float[2];
                    extractHues(img, set1, hues, sbAvg);
                    
                    //System.out.println("  " + rIdx + " hues=" + Arrays.toString(hues));
                    
                    if (isEmpty(hues) || Float.isNaN(sbAvg[0])) {
                        continue;
                    }
                                        
                    /*
                    bright sky:
                            sAvg=0.23  vAvg=0.75
                        rnbw hues peak is 2nd bin
                    dark sky:
                           sAvg=0.5   bck vAvg=0.43
                        rnbw hues peak is 1st peak
                    */
                    
                    int maxPeakIdx = MiscMath.findYMaxIndex(hues);
                    int hIdx = -1;
                    
                    //NOTE: this may need to be revised
                    if (sbAvg[1] < 0.55) {
                        if (maxPeakIdx == 0) {
                            hIdx = hists.size();
                            listOfSets.add(set1);
                            hists.add(new OneDIntArray(hues));
                            rgs.add(rg);
                        }
                    } else {
                        if (maxPeakIdx == 1) {
                            hIdx = hists.size();
                            listOfSets.add(set1);
                            hists.add(new OneDIntArray(hues));
                            rgs.add(rg);
                        }
                    }
                    
                    /*System.out.println(
                        "  " + rIdx + " xy=(" + rg.xC + "," + rg.yC + ") "
                        + " maxIdx=" + maxPeakIdx
                        + " hists.idx=" + hIdx
                        + " angle=" + (rg.orientation*180./Math.PI)
                        + " ecc=" + rg.eccentricity
                        + "\n  sAvg=" + sbAvg[0] 
                        + "    vAvg=" + sbAvg[1]
                        + " hues hist=" + Arrays.toString(hues)
                    );*/
                }
            }
        }        
    }
    
    private void findNegativePT(ImageExt img, 
        List<Region> polarThetaNegative, List<Set<PairInt>> listOfSets, 
        List<OneDIntArray> hists, List<RegionGeometry> rgs) {
        
        for (int rIdx = 0; rIdx < polarThetaNegative.size(); ++rIdx) {
            
            Region r = polarThetaNegative.get(rIdx);
           
            Set<PairInt> set1 = null;
                        
            for (int i = 0; i < r.accX.size(); ++i) {
                int x = r.accX.get(i);
                int y = r.accY.get(i);
                PairInt p = new PairInt(x, y);
                int pixIdx = img.getInternalIndex(p);
               
                if (set1 == null) {
                    set1 = new HashSet<PairInt>();
                }
                set1.add(p);
            }
            
            if (set1 != null) {
                
                RegionGeometry rg = Canonicalizer.calculateEllipseParams(
                    r, img.getWidth(), img.getHeight());

                //System.out.println("negative xy=" 
                //   + rg.xC + "," + rg.yC 
                //   + " ecc=" + rg.eccentricity);
                
                if (rg.eccentricity >= 0.95) {
                    
                    int[] hues = new int[10];
                    float[] sbAvg = new float[2];
                    extractHues(img, set1, hues, sbAvg);
                    
                    if (isEmpty(hues) || Float.isNaN(sbAvg[0])) {
                        continue;
                    }
                   
                    /*
                    bright sky:
                        sAvg=0.07   vAvg=0.86
                        hues hist=[0, 0, 0, 0, 3, 36, 0, 0, 0, 0]
                        rnbw hues peak is  bin
                    dark sky:
                        sAvg=  vAvg=
                        hues hist=
                        rnbw hues
                    */
                    
                    int maxPeakIdx = MiscMath.findYMaxIndex(hues);
                    
                    //NOTE: this may need to be revised
                    if (maxPeakIdx == 5) {
                        listOfSets.add(set1);
                        hists.add(new OneDIntArray(hues));
                        rgs.add(rg);
                    }
                    
                    /*System.out.println("negative xy=(" + rg.xC + "," + rg.yC + ") "
                        + " angle=" + (rg.orientation*180./Math.PI)
                        + " ecc=" + rg.eccentricity
                        + "\n   sAvg=" + sbAvg[0] 
                        + "\n   vAvg=" + sbAvg[1]
                        + "\n   hues hist=" + Arrays.toString(hues)
                    );*/
                   
                }
            }
        }        
    }
    
    private TIntList findLargeArc(List<VeryLongBitString> listOfSetBits, 
        List<Set<PairInt>> listOfSets, 
        List<OneDIntArray> hists0, List<RegionGeometry> rgs0, Image img) {
    
        int n = listOfSetBits.size();
        
        int[] sizes = new int[n];
        int[] indexes = new int[n];
        for (int i = 0; i < listOfSetBits.size(); ++i) {
            sizes[i] = (int)listOfSetBits.get(i).getNSetBits();
            indexes[i] = i;
        }
        QuickSort.sortBy1stArg(sizes, indexes);
        
        for (int i = n - 1; i > 0; --i) {
            int idx0 = indexes[i];
            VeryLongBitString bs = listOfSetBits.get(idx0);
            int nbs0 = (int)bs.getNSetBits();
            
            TIntList subsetIdxs = new TIntArrayList();
            int nInSubsets = 0;
            
            ParabolaLeastSquares polyFitter0 = null;
            float[] coeff0 = null;
            
            for (int j = i - 1; j > -1; --j) {
                int idx1 = indexes[j];
                VeryLongBitString bs2 = listOfSetBits.get(idx1);
                int nbs2 = (int)bs2.getNSetBits();
                
                VeryLongBitString inter = bs.and(bs2);
                int nbsi = (int)inter.getNSetBits();
                
                if ((nbs2 - nbsi) < 0.1*(nbs2)) {
                    subsetIdxs.add(idx1);
                    nInSubsets += nbs2;
                }
            }
            
            if (!subsetIdxs.isEmpty()) {
                
                // if not empty, check to see if set idx0 is an arc rather
                //    than a large level set which includes most of the image
                RegionGeometry rg0 = rgs0.get(idx0);
                double area = rg0.major * rg0.minor;
            
                double dens = ((double)listOfSets.get(idx0).size())/area;
                
                /*System.out.println("dens of x,y=" + rg0.xC + "," + rg0.yC
                    + " -> " + String.format("%.3f", dens)
                    + " nPts=" + nbs0 + " nInSubsets=" + nInSubsets    
                );*/
                
                if (dens > 0.1) {
                    
                    // fit a polynomial to rainbow points.  
                    // would prefer a circle, but the optical depth of the dispersers and the
                    // orientation of groups of them is not always a slab perpendicular to 
                    // the camera
                    
                    if (polyFitter0 == null) {
                        polyFitter0 = new ParabolaLeastSquares();
                        Set<PairInt> set = listOfSets.get(idx0);
                        
                        //y = c0*1 + c1*x[i] + c2*x[i]*x[i]
                        polyFitter0.accumulate(set);
                        coeff0 = polyFitter0.solve();
                        
                        polyFitter0.plotFit(coeff0, set, img.getWidth(),
                            img.getHeight(), i, 
                            "rainbow points");

                        double resid = polyFitter0.calcResiduals(coeff0, set);
                        
                        /*System.out.println("rainbow polynomial coefficients = " 
                            + Arrays.toString(coeff0));
                        System.out.println("image dimensions are " + img.getWidth() + " X "
                            + img.getHeight() + " pixels^2 " 
                            + " resid=" + resid
                            + " rg.coeff=" + Arrays.toString(rg0.m)
                        );*/

                        if (resid < 5) {
                            
                            TIntSet chk = new TIntHashSet(subsetIdxs);
                            
                            // since we have a polynomial now,
                            // look for other regions that may fit on the arc
                            for (int j = i - 1; j > -1; --j) {
                                int idx1 = indexes[j];
                                if (chk.contains(idx1)) {
                                    continue;
                                }
                                     
                                Set<PairInt> chkSet = listOfSets.get(idx1);
                                
                                ParabolaLeastSquares polyFitter2 = polyFitter0.copy();
                                polyFitter2.accumulate(chkSet);
                                
                                float[] coeff2 = polyFitter2.solve();
                                
                                Set<PairInt> tmp = new HashSet<PairInt>(set);
                                tmp.addAll(chkSet);
                                
                                double resid2 = polyFitter2.calcResiduals(coeff2, 
                                    tmp);
                                
                                //System.out.println("check " + " x=" + 
                                //    rgs0.get(idx1).xC + " y=" +
                                //    rgs0.get(idx1).yC + " resid=" + resid2);
                                
                                if (resid2 > 5) {
                                    continue;
                                }
                        
                                RegionGeometry chkRg = rgs0.get(idx1);
                                
                                //System.out.println("check " + " x=" + chkRg.xC + " y=" +
                                //    chkRg.yC + " resid=" + resid2);
                            
                                subsetIdxs.add(idx1);
                                //set.addAll(chkSet);
                            }
                            
                            subsetIdxs.add(idx0);
                            return subsetIdxs;
                        }
                    }
                }
            }
        }
        
        return null;
    }

    private List<VeryLongBitString> makeBitStrings(List<Set<PairInt>> 
        listOfSets, Image img) {
        
        List<VeryLongBitString> out = new ArrayList<VeryLongBitString>();
        for (Set<PairInt> set : listOfSets) {
            VeryLongBitString bs = new VeryLongBitString(img.getNPixels());
            for (PairInt p : set) {
                bs.setBit(img.getInternalIndex(p));
            }
            out.add(bs);
        }
        
        return out;
    }

    private int[] findAdjacent(NearestNeighbor2D nn, 
        List<Set<PairInt>> listOfSets, List<RegionGeometry> rgs) {

        TIntList idxs2 = new TIntArrayList();
        for (int i = 0; i < listOfSets.size(); ++i) {
            
            Set<PairInt> set = listOfSets.get(i);
            RegionGeometry rg = rgs.get(i);
            
            int n = set.size();
            int nAdj = 0;
            
            int d = 5;
            int m = 2 * (int)Math.round(rg.minor);
            
            for (PairInt p : set) {
                Set<PairInt> nearest = nn.findClosest(p.getX(), p.getY(), d);
                if (nearest != null && !nearest.isEmpty()) {
                    nAdj++;
                }
            }
            if (nAdj > 0) {
                float f = ((float)d)/(float)m;
                if (f > 1) {
                    f = 1;
                }
                //System.out.println("nAdj=" + nAdj 
                //    + " f*n=" + (f*n) + " setx=" + rg.xC + " y=" + rg.yC
                //    + " n=" + n);
                if (nAdj >= (f*n)) {
                    idxs2.add(i);
                }
            }
        }
        
        if (idxs2.isEmpty()) {
            return null;
        }
        
        return idxs2.toArray(new int[idxs2.size()]);
    }
    
    private TIntList findBestCombination(int[] idxs, 
        List<Set<PairInt>> listOfSets, List<RegionGeometry> rgs) {
        
        TIntList bestIdxs = new TIntArrayList();
        
        double minResid = Double.MAX_VALUE;
        
        for (int k = idxs.length; k > 0; k--) {
            
            int[] selectedIndexes = new int[k];
            
            SubsetChooser subsetChooser = new SubsetChooser(idxs.length, k);
            
            int nV = subsetChooser.getNextSubset(selectedIndexes);
            
            ParabolaLeastSquares polyFitter = new ParabolaLeastSquares();
        
            while (nV != -1) {
                
                Set<PairInt> subset = new HashSet<PairInt>();

                StringBuilder sb = new StringBuilder();
                
                for (int bitIndex : selectedIndexes) {

                    Set<PairInt> g = listOfSets.get(bitIndex);

                    subset.addAll(g);
                    
                    RegionGeometry rg = rgs.get(bitIndex);
                    
                    sb.append(String.format(" (%d,%d) ", rg.xC, rg.yC));
                }
                
                polyFitter.accumulate(subset);
                                
                float[] coeff = polyFitter.solve();
                
                if (coeff == null) {
                    continue;
                }
               
                double resid = polyFitter.calcResiduals(coeff, subset);
                
                if (resid > 5) {
                    nV = subsetChooser.getNextSubset(selectedIndexes);
                    continue;
                }
                
                sb.append(String.format(" resid=%.3f\n", resid));
                
                if (resid < minResid) {
                    minResid = resid;
                    bestIdxs.clear();
                    for (int bitIndex : selectedIndexes) {
                        bestIdxs.add(bitIndex);
                    }
                }
                
                System.out.println(sb.toString());
            
                nV = subsetChooser.getNextSubset(selectedIndexes);
            }
        }

        return bestIdxs;
    }

}
