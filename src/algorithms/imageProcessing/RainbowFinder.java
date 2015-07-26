package algorithms.imageProcessing;

import algorithms.compGeometry.PointInPolygon;
import algorithms.imageProcessing.SkylineExtractor.RemovedSets;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.ArrayPair;
import algorithms.util.PairInt;
import algorithms.util.PolynomialFitter;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

/**
 * class with methods to find a rainbow within an image and to create a hull
 * to encapsulate it for various methods.
 * 
 *  TODO: Note, it may be necessary to build a hull from a spine of
    more than 10 points for complex images w/ rainbow intersecting
    multiple times with structured foreground and sky
    (see the end of method createRainbowHull).
        
 * @author nichole
 */
public class RainbowFinder {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    static class Hull {
        float[] xHull;
        float[] yHull;
    }
    
    private Set<PairInt> outputRainbowPoints = new HashSet<PairInt>();
    
    private Set<PairInt> excludePointsInRainbowHull = new HashSet<PairInt>();
    
    private Hull rainbowHull = null;
    
    private float[] rainbowCoeff = null;
    
    private float[] rainbowSkeletonX = null;
    
    private float[] rainbowSkeletonY = null;
    
    private float hullHalfWidth = 0;
    
    /**
     * search for a rainbow within the image, and if found, create a hull of
     * points around that encapsulates image points that are part of the
     * expanded rainbow.  The results are kept in member variables
     * rainbowCoeff, rainbowHull, outputRainbowPoints, and 
     * excludePointsInRainbowHull.
     * 
     * @param skyPoints
     * @param reflectedSunRemoved
     * @param colorImg
     * @param xOffset
     * @param yOffset
     * @param imageWidth
     * @param imageHeight
     * @param skyIsDarkGrey
     * @param removedSets 
     */
    public void findRainbowInImage(Set<PairInt> skyPoints, 
        Set<PairInt> reflectedSunRemoved,
        ImageExt colorImg, int xOffset, int yOffset, 
        int imageWidth, int imageHeight,
        boolean skyIsDarkGrey, RemovedSets removedSets) {
    
        rainbowCoeff = findRainbowPoints(skyPoints, 
            removedSets.getReflectedSunRemoved(), colorImg, 
            xOffset, yOffset, skyIsDarkGrey, outputRainbowPoints);

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
        
        /* n=21
         0,20    1    2    3    4    5    6    7    8    9
                                                  
         19   18   17   16   15   14   13   12   11   10
        
        Points 0 and 19 are opposite points on the hull, so are 1 and 18, etc.
        So can do quick check that for each segment such as 0->1->18->19->0
        that there are skyPoints surrounding them.
        If not, and if the coverage is partial, will reduce the segment polygon
        size and only add points to skypoints within the segment polygon.
        */ 
        
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

    /**
     * search for rainbow colored points over the entire image, fit an
     * ellipse to them, and assert that the points have certain colors in
     * them.  If the original fit to an ellipse is not good, the
     * method divides the rainbow points by contiguous subsets to find best
     * and similar fits.  The last step of color requirement helps to rule
     * out half ellipse patterns in rocks for instance that have only rock
     * colors in them. 
     * The return polynomial coefficients are float[]{c0, c1, c2}
     * where y[i] = c0 + c1 * x[i] + c2 * x[i] * x[i].
     * when c2 is negative, the parabola peak is upward (higher y).
       c2 also indicates the magnitude of the points in powers of 1/10.
     * @param skyPoints
     * @param reflectedSunRemoved
     * @param colorImg
     * @param xRelativeOffset
     * @param yRelativeOffset
     * @param outputRainbowPoints output variable to return found rainbow
     * points if any
     * @return polynomial fit coefficients to 
     * y[i] = c0*1 + c1*x[i] + c2*x[i]*x[i].  this may be null if a fit wasn't
     * possible.
     */
    float[] findRainbowPoints(Set<PairInt> skyPoints, 
        Set<PairInt> reflectedSunRemoved,
        ImageExt colorImg, int xOffset, int yOffset, boolean skyIsDarkGrey,
        Set<PairInt> outputRainbowPoints) {

        Set<PairInt> rainbowPoints = findRainbowColoredPoints(colorImg, 
            reflectedSunRemoved, xOffset, yOffset, skyIsDarkGrey);

        if (rainbowPoints.isEmpty()) {
            return null;
        }
        
        if (rainbowPoints.size() < 12) {
            return null;
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

                boolean isWhite = cieC.isCentralWhite(cieX, cieY);
                
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
            
            float nBroadlyRedFrac = 
                ((float)nBroadlyRed/(float)bestFittingPoints.size());
    
            log.info("nGTX=" + nGTX + " nLTY=" + nLTY + " n=" 
                + bestFittingPoints.size() + " "
                + " CIE: minCIEX=" + minCIEX + " maxCIEX=" + maxCIEX
                + " minCIEY=" + minCIEY + " maxCIEY=" + maxCIEY
                + " range=(" + cieXRange + "," + cieYRange + ")"
                + "\n nPurpRed=" + nPurpRed + " nOranRed=" + nOranRed
                + " nYellow=" + nYellow + " nGreen=" + nGreen + " nRed=" + nRed
                + " nOrange/nPurpRed=" + ((float)nOranRed/(float)nPurpRed)
                + " nOrange/nTot=" + ((float)nOranRed/(float)bestFittingPoints.size())
                + " nPurpRed/nTot=" + ((float)nPurpRed/(float)bestFittingPoints.size())
                + " nBroadlyRed/nTot=" + nBroadlyRedFrac
            );

            rainbowPoints.clear();

            //TODO: this entire section could use more testing
            
            /* 
            assert that orange and red are present
            */
            
            if (nRed == 0 && nOranRed == 0 && nPurpRed == 0) {
                return null;
            }

            /*
            would like to assert purple and green too.  difficult because
            green was not included in the rainbow colored point gathering 
            because so much of landscape is possibly green.
            so assert yellow for greeen.
            */
            float greenFraction = (float)nGreen/(float)bestFittingPoints.size();
            float yellowFraction = (float)nYellow/(float)bestFittingPoints.size();
            float purpleFraction = (float)nPurpRed/(float)bestFittingPoints.size();
            if ( 
                ((greenFraction < 0.01) && (yellowFraction < 0.05)) 
                //|| ((purpleFraction < 0.05) && (greenFraction < 0.03))
            ) {
                
                return null;
                
            }/* else if ((cieXRange < 0.08) && (cieYRange < 0.08)) {
                return null;
            }*/
 
            float rdr = ((float)(nOranRed + nYellow))/(float)nPurpRed;
            
            //TODO: this could use alot of tests and revision
            if ((nGTX > 10) && 
                ((nLTY > 10) || 
                  (!skyIsDarkGrey &&
                    (
                        (rdr > 1.1) && (nGreen > nPurpRed))
                  )
                  || (skyIsDarkGrey && (rdr > 1.1))
                )
                ) {
                
                float frac = (float)(nGTX + nLTY)/(float)bestFittingPoints.size();
                if (frac > 0.002) {
                    //NOTE: this doesn't include the one excluded thru random sampling
                    rainbowPoints.addAll(bestFittingPoints);
                }
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
   
     /**
     * search over the entire image for pixels that are rainbow colored.
     * 
     * @param colorImg
     * @param reflectedSunRemoved
     * @param xOffset
     * @param yOffset
     * @param skyIsDarkGrey
     * @return rainbowPoints pixels from the entire image containing rainbow 
     * colors.
     */
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
if (col==295 && row>= 173 && row <= 175) {
  int z = 1;
}      
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
                    
                } /*else if (colors.isInYellowishYellowGreen(r, g, b, cieX, cieY, 
                    saturation, brightness)) {
                    
                    // finds grass and trees
                
                    set.add(p);
                    
                }*/
            }
        }
        
        //debugPlot(set, colorImg, xOffset, yOffset, "rainbow_colored_points");
        
        return set;
    }
    
     /**
     * a method to roughly fit a hull around the rainbow polynomial described
     * by rainbowCoeff and populated by rainbowPoints.
     * 
     * @param skyPoints
     * @param rainbowCoeff coefficients for a 2nd order polynomial
     * @param rainbowPoints  rainbow colored points in the image that fit a polynomial.
     * there should be 10 or more points at least.
     * @param originalColorImage
     * @param xOffset
     * @param yOffset
     * @return 
     */
    private Hull createRainbowHull(float[] rainbowCoeff, 
        Set<PairInt> rainbowPoints, ImageExt originalColorImage, int xOffset, 
        int yOffset) {
        
        /*
        need to know the furthest closest distance to the polynomial, that is the
        furthest perpendicular point to the polynomial in order to expand the 
        region around the polynomial to become a hull that encloses all rainbow 
        points.
        
        (1) Could determine for every point, the distance to the polynomial:
        Robust and Efficient Computation of the Closest Point on a Spline Curve" 
        Hongling Wang, Joseph Kearney, and Kendall Atkinson
        http://homepage.cs.uiowa.edu/~kearney/pubs/CurvesAndSufacesClosestPoint.pdf
        implemented by:
        http://www.ogre3d.org/tikiwiki/tiki-index.php?page=Nearest+point+on+a+Spline
        
        -- The method for nearest point on a spline requires creating a
        spline out of the polynomial.  can assume that something like 1000
        points would be necessary, though maybe 100 would be fine.
        -- Then the method requires a good starting guess for 3 points on the
        polynomial that would be near the true perpendicular point.
        One could guess the first 3 splines by making the polynomial roughly
        10 splines separately and doing a distance test to each spline
        for each point (with shortcuts for not needing to check).
        -- Then about 10 iterations or less to find the best answer.
        
        The runtime might look roughly like
        N_poly +  N_points * 10 distance tests + N_points * 10 iterations of the 
        closestToSpline algorithm.
        
        (2) Or could instead, characterize the hull by 10 points and calculate 
        the coordinates of the points projected perpendicular to them above and 
        below at distances that are the hull of the rainbow.
        Then evaluation of the size would be "point in polygon" tests for
        each point.
        Improved estimates of the hull size can use a binomial search pattern
        to expand or shrink the distance estimate used for furthest extension
        of the hull.
        The first estimate of the hull size can use the rough distance point test
        just as above, by comparing each point to the polynomial as 10 segments.
        
        runtime is roughly:
            N_points * 10 distance tests +
            N_iter * N_points *
               "point in polygon" (which is roughly O(N_poly) where N_poly would be ~20)
            = N_points * 10 + N_iter * N_points * 20 for each iteration
        
        expect that N_iter is probably <= 5 because it's close already.
        
        So the 2nd method will be fine for the purposes here, but the first
        method is interesting for use cases which need precision.
        */
        
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
        
        /* 
        n=21
         0,20    1    2    3    4    5    6    7    8    9
                                                  
         19   18   17   16   15   14   13   12   11   10
        */
        
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
            
            /*log.info(String.format(
                "(%d,%d) is closest to poly point (%f,%f) dist=%f", x, y, 
                xc[minIdx], yc[minIdx], (float)Math.sqrt(minDistSq)));
            */
        }
        
        return (float)Math.sqrt(maxDistSq);
    }

    /**
    populate outputXPoly and outputYPoly with points perpendicular to x and y
    * at distances dist.  note that the lengths of outputXPoly and outputYPoly
    * should be 2*x.length+1.
    * Also note that one wants the separation between points in x and y
    * to be larger than dist (else, the concave portion of writing the hull
    * has a retrograde order and appearance).
    */
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
        
        /*
        y = c0*1 + c1*x[i] + c2*x[i]*x[i]
        
        dy/dx = c1 + 2*c2*x[i] = tan theta
        */
        
        int n = outputXPoly.length;
        
        /*
        want them in order so using count0 and count1
        n=21
        
        0,20    1    2    3    4    5    6    7    8    9
                                                  
         19    18   17   16   15   14   13   12   11   10
        */
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
    
}
