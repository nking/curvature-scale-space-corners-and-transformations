package algorithms.imageProcessing;

import algorithms.imageProcessing.Sky.SkyObject;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.RegionGeometry;
import algorithms.imageProcessing.features.mser.EllipseHelper;
import algorithms.imageProcessing.features.mser.MSEREdges;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.VeryLongBitString;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
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
    public SkyObject[] findRainbows(MSEREdges mserEdges) {
        
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

        findPositivePT(img, polarThetaPositive, listOfSets0, hists0, rgs0);
        
        
        List<Region> polarThetaNegative = mserEdges.getOrigGsPtRegions().get(3);
        
        List<Set<PairInt>> listOfSets1 = new ArrayList<Set<PairInt>>();
        
        List<OneDIntArray> hists1 = new ArrayList<OneDIntArray>();
        
        List<RegionGeometry> rgs1 = new ArrayList<RegionGeometry>();

        findNegativePT(img, polarThetaNegative, listOfSets1, hists1, rgs1);
                
        // points from the same rainbow in the positive image
        // should have similar hue histograms.
        // gathered those into groups here.
        List<VeryLongBitString> lists = clusterByIntersection(hists0, 0.95f);
        for (int i = 0; i < lists.size(); ++i) {
            VeryLongBitString bs = lists.get(i);
            int[] histIdxs = bs.getSetBits();
            
            List<Set<PairInt>> setsI = new ArrayList<Set<PairInt>>(histIdxs.length);
            List<RegionGeometry> rgsI = new ArrayList<RegionGeometry>(histIdxs.length);
            for (int hIdx : histIdxs) {
                //int idx = indexes.get(hIdx);
                setsI.add(listOfSets0.get(hIdx));
                RegionGeometry rg = rgs0.get(hIdx);
                rgsI.add(rg);
            
                System.out.println(i + ") " 
                    + " xy=(" + rg.xC + "," + rg.yC + ") "
                    + " angle=" + (rg.orientation*180./Math.PI)
                    + " ecc=" + rg.eccentricity
                );
            }
            
            // there may be more than one rainbow in this set.
            // for example, a test image has 1 strong rainbow and 2 fainter ones.
            
            // combine these positive mser w/ the complementary MSER regions found in the
            //    negative polar theta image here. for the strong mser regions
            //    in the positive image that are rainbow arcs, there are
            //    complementary arcs which appear in the
            //    negative image near them.
            
            // TODO: paused here
        }
        
        /*can use a clustering based on intersection limit
        to get clusters of points and fit ellipses to them.
        because of the eccentricity limit, the individual point set
        mser regions already have elliptical fits.
        */  
        //SkyObject obj = new SkyObject();
        //obj.points = listOfSets2.get(0);
        //obj.xyCenter = ehs.get(0).getXYCenter();
        //return obj;
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    /*
    public void findRainbowInImage(Set<PairInt> skyPoints, 
        Set<PairInt> reflectedSunRemoved,
        ImageExt colorImg, int xOffset, int yOffset, 
        int imageWidth, int imageHeight,
        boolean skyIsDarkGrey, GroupPixelColors allSkyColor, 
        RemovedSets removedSets) {

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
    
    public Set<PairInt> getPointsToExcludeInHull() {
        return excludePointsInRainbowHull;
    }
    
    public Set<PairInt> getRainbowPoints() {
        return outputRainbowPoints;
    }
    
    public float[] getRainbowCoeff() {
        return rainbowCoeff;
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
            
            AbstractSkyRainbowColors colors = skyIsDarkGrey ?
                new DarkSkyRainbowColors() : new BrightSkyRainbowColors();
            
            Set<PairInt> bestFittingPoints = new HashSet<PairInt>();
            
            coef = polyFitter.solveForBestFittingContiguousSubSets(
                rainbowPoints, bestFittingPoints, colorImg.getWidth(), 
                colorImg.getHeight());
            
            if (coef == null) {
                return null;
            }
            
            int w = colorImg.getWidth() - xOffset;
            int h = colorImg.getHeight() - yOffset;
            if (bestFittingPoints.size() > (0.3f * (float)(w*h))) {
                return null;
            }
 
            // assert that colors are as expected, that is, that we don't
            // have only green and blue from rocks
            
            // filter to keep only those with a significant fraction with 
            //    cieX > 0.4 and cieY < 0.3
            // filter to keep only those with red
            
            int nGTX = 0;
            int nLTY = 0;
            int nPurpRed = 0;
            int nOranRed = 0;
            int nYellow = 0;
            int nGreen = 0;
            int nRed = 0;
            int nWhite = 0;
            int nBroadlyRed = 0;

            //MiscDebug.plotSkyColor(bestFittingPoints, colorImg, xOffset, yOffset);

            CIEChromaticity cieC = new CIEChromaticity();
            ArrayPair cPurpRed = cieC.getRedThroughPurplishRedPolynomial();
            ArrayPair cOranRed = cieC.getOrangePolynomial();
            ArrayPair cYellow = cieC.getGreenishYellowThroughOrangePolynomial();
            ArrayPair cGreen = cieC.getGreenPolynomial();
            ArrayPair cRed = cieC.getRedPolynomial();
            PointInPolygon pInPoly = new PointInPolygon();
            
            float minCIEX = Float.MAX_VALUE;
            float maxCIEX = Float.MIN_VALUE;
            float minCIEY = Float.MAX_VALUE;
            float maxCIEY = Float.MIN_VALUE;
            for (PairInt p : bestFittingPoints) {
                
                int x = p.getX();
                int y = p.getY();
                
                int idx = colorImg.getInternalIndex(x + xOffset, y + yOffset);
                
                int r = colorImg.getR(idx);
                int g = colorImg.getG(idx);
                int b = colorImg.getB(idx);
                float rDivB = (float)r/(float)b;
                float cieX = colorImg.getCIEX(idx);
                float cieY = colorImg.getCIEY(idx);

log.fine(String.format(
                    "(%d,%d) r=%d, g=%d, b=%d  rDivB=%f  cieX=%f  cieY=%f  hue=%f",
                    x, y, r, g, b, rDivB, cieX, cieY, colorImg.getHue(idx)));

                boolean isWhite = cieC.isWhite(cieX, cieY);
                
                if (cieX >= 0.35) {
                    nGTX++;
                }
                if ((cieY <= 0.3) || (rDivB < 0.89f)) {
                    nLTY++;
                }
                if (cieX < minCIEX) {
                    minCIEX = cieX;
                }
                if (cieX > maxCIEX) {
                    maxCIEX = cieX;
                }
                if (cieY < minCIEY) {
                    minCIEY = cieY;
                }
                if (cieY > maxCIEY) {
                    maxCIEY = cieY;
                }
                
                if (!isWhite) {
                    if (pInPoly.isInSimpleCurve(cieX,cieY, cPurpRed.getX(), 
                        cPurpRed.getY(), cPurpRed.getX().length)) {
                        nPurpRed++;
                    }
                    if (pInPoly.isInSimpleCurve(cieX, cieY, cOranRed.getX(), 
                        cOranRed.getY(), cOranRed.getX().length)) {
                        nOranRed++;
                    }
                    if (pInPoly.isInSimpleCurve(cieX, cieY, cYellow.getX(), 
                        cYellow.getY(), cYellow.getX().length)) {
                        nYellow++;
                    }
                    if (pInPoly.isInSimpleCurve(cieX, cieY, cGreen.getX(), 
                        cGreen.getY(), cGreen.getX().length)) {
                        nGreen++;
                    }
                    if (pInPoly.isInSimpleCurve(cieX, cieY, cRed.getX(), 
                        cRed.getY(), cRed.getX().length)) {
                        nRed++;
                    }
                }
                
                if (colors.isInBroadRedRainbowSpace(cieX, cieY)) {
                    nBroadlyRed++;
                }
            }
            
            float cieXRange = maxCIEX - minCIEX;
            float cieYRange = maxCIEY - minCIEY;
            
            float nTot = (float)bestFittingPoints.size();
            
            float nOrangeFrac = (float)nOranRed/nTot;
            
            float nPurpleRedFrac = (float)nPurpRed/nTot;
            
            float nGreenFrac = (float)nGreen/nTot;
            
            float nYellowFrac = (float)nYellow/nTot;
            
            float nBroadlyRedFrac = (float)nBroadlyRed/nTot;
    
            log.info("nGTX=" + nGTX + " nLTY=" + nLTY + " n=" 
                + bestFittingPoints.size() + " "
                + " CIE: minCIEX=" + minCIEX + " maxCIEX=" + maxCIEX
                + " minCIEY=" + minCIEY + " maxCIEY=" + maxCIEY
                + " range=(" + cieXRange + "," + cieYRange + ")"
                + "\n nPurpRed=" + nPurpRed + " nOranRed=" + nOranRed
                + " nYellow=" + nYellow + " nGreen=" + nGreen + " nRed=" + nRed
                + " nOrange/nPurpRed=" + ((float)nOranRed/(float)nPurpRed)
                + " nOrange/nTot=" + nOrangeFrac
                + " nPurpRed/nTot=" + nPurpleRedFrac
                + " nGreen/nTot=" + nGreenFrac
                + " nYellow/nTot=" + nYellowFrac
                + " nBroadlyRed/nTot=" + nBroadlyRedFrac
            );

            rainbowPoints.clear();

            //TODO: this entire section could use more testing
MiscDebug.plotPoints(bestFittingPoints, imageWidth, imageHeight, 
MiscDebug.getCurrentTimeFormatted());

            boolean colorsAreNotRainbow = false;
            
            if (skyIsDarkGrey) {
                //TODO: this needs to be tested on more night images.
                // rainbow: nOrange/nTot=0.36870417
                //          nPurpRed/nTot=0.2621027
                //          nBroadlyRed/nTot=0.65574574
                //          nGreen/nTot=0.19
                //          nYellow/nTot=0.44
                if (nPurpleRedFrac < 0.1 && nOrangeFrac < 0.1 && nGreenFrac < 0.1) {
                    colorsAreNotRainbow = true;
                }
            } else {
                if (
                    (nBroadlyRedFrac > 0.1 && nOrangeFrac > 0.05)
                    || (nPurpleRedFrac < 0.005)
                    || (nYellowFrac > 0.8)
                    ){
                    colorsAreNotRainbow = true;
                }
            }
            
            if (!colorsAreNotRainbow) {
                rainbowPoints.addAll(bestFittingPoints);
            }
        }
        
        outputRainbowPoints.addAll(rainbowPoints);
        
        if (!rainbowPoints.isEmpty()) {
            polyFitter.plotFit(coef, outputRainbowPoints, colorImg.getWidth(),
                colorImg.getHeight(), 234, "rainbow points");
        }
        
        log.info("rainbow points size=" + outputRainbowPoints.size());
 
        return coef;
    }
   
    Set<PairInt> findRainbowColoredPoints(ImageExt colorImg, 
        Set<PairInt> reflectedSunRemoved,
        int xOffset, int yOffset, boolean skyIsDarkGrey) {
        
        AbstractSkyRainbowColors colors = skyIsDarkGrey ?
            new DarkSkyRainbowColors() : new BrightSkyRainbowColors();
        
        Set<PairInt> set = new HashSet<PairInt>();
        
        for (int col = 0; col < colorImg.getWidth(); col++) {
            for (int row = 0; row < colorImg.getHeight(); row++) {
                
                PairInt p = new PairInt(col - xOffset, row - yOffset);
                
                if (reflectedSunRemoved.contains(p)) {
                    continue;
                }
                
                int idx = colorImg.getInternalIndex(col, row);
      
                float cieX = colorImg.getCIEX(idx);
                float cieY = colorImg.getCIEY(idx);

                int r = colorImg.getR(idx);
                int g = colorImg.getG(idx);
                int b = colorImg.getB(idx);

                float saturation = colorImg.getSaturation(idx);
                float brightness = colorImg.getBrightness(idx);
        
                if (colors.isInRedThroughPurplishRed(r, g, b, cieX, cieY, 
                    saturation, brightness)) {
                    
                    set.add(p);

                } else if (colors.isInOrangeRed(r, g, b, cieX, cieY, saturation, 
                    brightness)) {
                    
                    set.add(p);
                
                } else if (colors.isInGreenishYellowOrange(r, g, b, cieX, cieY,
                    saturation, brightness)) {
                 
                    set.add(p);
                    
                }
            }
        }
        
        //debugPlot(set, colorImg, xOffset, yOffset, "rainbow_colored_points");
        
        return set;
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
     
    protected int nPointsInPolygon(Set<PairInt> rainbowPoints, 
        float[] xPoly, float[] yPoly) {
        
        PointInPolygon pIn = new PointInPolygon();
        
        int nInside = 0;
        
        for (PairInt p : rainbowPoints) {
            float x = p.getX();
            float y = p.getY();
            
            boolean ans = pIn.isInSimpleCurve(x, y, xPoly, yPoly, xPoly.length);
            
            if (ans) {
                nInside++;
            }
        }
        
        return nInside;
    }
    
    void debugPlot(Set<PairInt> extSkyPoints, Image originalColorImage, 
        int xOffset, int yOffset, String outputPrefixForFileName) {
        
        //plot is made in aspects
        
    }
    
    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("Number of rainbowPoints=")
            .append(Integer.toString(outputRainbowPoints.size())).append("\n");
        
        if (rainbowCoeff != null) {
            sb.append("rainbow fit coeff=").append(Arrays.toString(rainbowCoeff))
                .append("\n");
        }
        
        return sb.toString();
    }
    */
    
    private EllipseHelper createRegion(Set<PairInt> points) {

        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        int[] xyCen = ch.calculateRoundedXYCentroids(points);
    
        PairIntArray xy = Misc.convertWithoutOrder(points);
        
        EllipseHelper eh = new EllipseHelper(xyCen[0], xyCen[1], xy);
        
        return eh;
    }
    
    private int[] extractHues(ImageExt img, Set<PairInt> set) {
        
        int[] hues = new int[10];
        
        float[] hsb = new float[3];
        for (PairInt p : set) {
            int r = img.getR(p);
            int g = img.getG(p);
            int b = img.getB(p);
            Color.RGBtoHSB(r, g, b, hsb);
            
            // TODO: consider satudation and brightness limits
            int bin = (int)(hsb[0]/0.1f);
            if (bin == 10) {
                bin--;
            }
            hues[bin]++;
        }
        
        return hues;
    }
    
    private Set<PairInt> extractOffRegion(ImageExt img, Region r, 
        RegionGeometry rg) {
        
        // region orientation is the direction that the minor axis points
        //  to, so is the direction needed for an offset patch
        // and so is 180 from that.
        
        int eAngle = (int)Math.round(rg.orientation * 180./Math.PI);
            // put into 0 to 180 ref frame
        if (eAngle > 179) {
            eAngle -= 180;
        }
        // put into ref frame of dominant orientations (major axis direction)
        eAngle -= 90;
        if (eAngle < 0) {
            eAngle += 180;
        }
        double d = eAngle * Math.PI/180.;
        
        int x = rg.xC;
        int y = rg.yC;
        
        int delta = (int)((rg.minor + rg.major)/2.);
        
        int x2 = (int)Math.round(x + 2.f * delta * Math.cos(d));
        int y2 = (int)Math.round(y + 2.f * delta * Math.sin(d));

        int w = img.getWidth();
        int h = img.getHeight();

        Set<PairInt> out = new HashSet<PairInt>();
        
        if (isWithinBounds(x2, y2, w, h)) {
        
            Set<PairInt> points = new HashSet<PairInt>();
            
            for (int i = (x2 - delta); i < (x2 + delta); ++i) {
               for (int j = (y2 - delta); j < (y2 + delta); ++j) {
                   if (isWithinBounds(i, j, w, h)) {
                       points.add(new PairInt(i, j));
                   }
               }
            }
            out.addAll(points);
        }

        int x3 = (int)Math.round(x - 2.f * delta * Math.cos(d));
        int y3 = (int)Math.round(y - 2.f * delta * Math.sin(d));
        if (isWithinBounds(x3, y3, w, h)) {
        
            Set<PairInt> points = new HashSet<PairInt>();
            
            for (int i = (x3 - delta); i < (x3 + delta); ++i) {
               for (int j = (y3 - delta); j < (y3 + delta); ++j) {
                   if (isWithinBounds(i, j, w, h)) {
                       points.add(new PairInt(i, j));
                   }
               }
            }
            out.addAll(points);
        }
        
        return out;
    }
    
    private void extractHSBProperties(ImageExt img, Set<PairInt> points, 
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
    
    private float[] normalize(int[] hues) {
        
        float[] norm0 = new float[hues.length];
        int n0 = 0;
        for (int h : hues) {
            n0 += h;
        }
        
        for (int i = 0; i < hues.length; ++i) {
             norm0[i] = (float)hues[i]/(float)n0;    
        }
        
        return norm0;
    }

    private int[] subtractNormalized(int[] hues, int[] offRegion) {
    
        float[] h = normalize(hues);
        float[] bck = normalize(offRegion);
        
        /*System.out.println(
            "hues=" + Arrays.toString(hues) + 
            " off=" + Arrays.toString(offRegion) + 
            " h=" + Arrays.toString(h) + 
            " bck=" + Arrays.toString(bck));
        */
        float tot = 0;
        for (int i = 0; i < h.length; ++i) {
            h[i] -= bck[i];
            if (h[i] < 0) {
                h[i] = 0;
            }
            tot += h[i];
        }
        
        float factor = 100.f/tot;
        int[] out = new int[hues.length];
        for (int i = 0; i < h.length; ++i) {
            out[i] = Math.round(factor * h[i]);
        }
        
        return out;
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
                
                if (rg.eccentricity >= 0.95) {
                    
                    int[] hues = extractHues(img, set1);
                    
                    if (isEmpty(hues)) {
                        continue;
                    }
                    
                    Set<PairInt> offRegionPoints = extractOffRegion(img, r, rg);
                    
                    int[] offRegion = new int[hues.length];
                    float[] sbAvg = new float[2];
                    extractHSBProperties(img, offRegionPoints, offRegion, sbAvg);
                    
                    if (Float.isNaN(sbAvg[0])) {
                        continue;
                    }
                    
                    int[] hues2 = subtractNormalized(hues, offRegion);
                    
                    /*
                    bright sky:
                            sAvg=0.23220387   bck vAvg=0.5838191
                        bck sAvg=0.16768077   bck vAvg=0.80450857
                        rnbw hues peak is 2nd bin
                    dark sky:
                           sAvg=0.26138106   bck vAvg=0.25395915
                        bck sAvg=0.3070194   bck vAvg=0.2282818
                        rnbw hues peak is 1st peak
                    */
                    
                    int maxPeakIdx = MiscMath.findYMaxIndex(hues2);
                    
                    //NOTE: this may need to be revised
                    if (sbAvg[1] < 0.4) {
                        if (maxPeakIdx == 0) {
                            listOfSets.add(set1);
                            hists.add(new OneDIntArray(hues2));
                            rgs.add(rg);
                        }
                    } else {
                        if (maxPeakIdx == 1) {
                            listOfSets.add(set1);
                            hists.add(new OneDIntArray(hues2));
                            rgs.add(rg);
                        }
                    }
                    
                    /*
                    System.out.println("xy=(" + rg.xC + "," + rg.yC + ") "
                        + " angle=" + (rg.orientation*180./Math.PI)
                        + " ecc=" + rg.eccentricity
                        + "\n   bck sAvg=" + sbAvg[0] 
                        + "   bck vAvg=" + sbAvg[1]
                        + " hues hist=" + Arrays.toString(hues2)
                        + "\n   bck=" + Arrays.toString(offRegion)
                    );
                    */
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
                    
                    int[] hues = extractHues(img, set1);
                    
                    if (isEmpty(hues)) {
                        continue;
                    }
                    
                    Set<PairInt> offRegionPoints = extractOffRegion(img, r, rg);
                    
                    int[] offRegion = new int[hues.length];
                    float[] sbAvg = new float[2];
                    extractHSBProperties(img, offRegionPoints, offRegion, sbAvg);
                    
                    if (Float.isNaN(sbAvg[0])) {
                        continue;
                    }
                    
                    int[] hues2 = subtractNormalized(hues, offRegion);
                    
                    /*
                    bright sky:
                        bck sAvg=0.07984637   bck vAvg=0.7707274
                        hues hist=[0, 0, 0, 0, 6, 94, 0, 0, 0, 0]
                        rnbw hues peak is  bin
                    dark sky:
                        bck sAvg=0.19074595   bck vAvg=0.44089583
                        hues hist=[0, 0, 0, 0, 0, 100, 0, 0, 0, 0]
                        rnbw hues peak is peak
                    */
                    
                    int maxPeakIdx = MiscMath.findYMaxIndex(hues2);
                    
                    //NOTE: this may need to be revised
                    if (maxPeakIdx == 5) {
                        listOfSets.add(set1);
                        hists.add(new OneDIntArray(hues2));
                        rgs.add(rg);
                    }
                    
                    System.out.println("negative xy=(" + rg.xC + "," + rg.yC + ") "
                        + " angle=" + (rg.orientation*180./Math.PI)
                        + " ecc=" + rg.eccentricity
                        + "\n   bck sAvg=" + sbAvg[0] 
                        + "   bck vAvg=" + sbAvg[1]
                        + "\n   hues hist=" + Arrays.toString(hues2)
                        + "\n   bck=" + Arrays.toString(offRegion)
                    );
                   
                }
            }
        }        
    }
}
