package algorithms.imageProcessing;

import algorithms.compGeometry.EllipseHelper;
import algorithms.misc.Histogram;
import algorithms.util.LinearRegression;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class SunFinder {
    
    private Set<PairInt> sunPoints = new HashSet<PairInt>();
    
    private double[] sunCoeff = null;
    
    private double pointDensity = 0;
    
    private Logger log = Logger.getLogger(this.getClass().getName());

    public SunFinder() {
    }
    
    /**
     * find sun colored points by searching entire image.
     * The sun's visible radius is the photosphere and that size is known so
     * if the points returned are well fit by a circle
     * or a fraction of a fittable circle, one gets a limit
     * on scale in the image (knowledge of the camera and lens are needed too
     * for projection).
     * 
     * @param colorImg
     * @param xOffset
     * @param yOffset
     * @param skyIsDarkGrey 
     */
    public void findSunPhotosphere(ImageExt colorImg, int xOffset, 
        int yOffset, boolean skyIsDarkGrey) {
        
        SunColors sunColors = new SunColors();
        
        if (skyIsDarkGrey) {
            sunColors.useDarkSkiesLogic();
        }
        
        for (int col = 0; col < colorImg.getWidth(); col++) {
            for (int row = 0; row < colorImg.getHeight(); row++) {
                
                int idx = colorImg.getInternalIndex(col, row);

                if (sunColors.isSunCenterColor(colorImg, idx)) {
                    sunPoints.add(new PairInt(col - xOffset, row - yOffset));
                }
            }
        }
        
        log.info("found " + sunPoints.size() 
            + " sun points from image of width " + colorImg.getWidth());

        if (sunPoints.size() < 6) {
            sunPoints.clear();
            return;
        }
        
        //fit ellipse to yellowPoints.  ellipse because of possible occlusion.
        EllipseHelper ellipseHelper = new EllipseHelper();
        double[] params = ellipseHelper.fitEllipseToPoints(sunPoints);
        
        if (params == null) {
            // not close to an ellipse
            sunPoints.clear();
            return;
        }
        
        float xc = (float)params[0];
        float yc = (float)params[1];
        float a = (float)params[2];
        float b = (float)params[3];
        float alpha = (float)params[4];
        
        log.info("elliptical fit to yellow points: " + Arrays.toString(params));
       
        //TODO:  this may need adjustments.  top of sun rising over mountains..
        /*if ((a/b) > 6) {
            return new HashSet<PairInt>();
        }*/
        
        sunCoeff = params;
        
        calculatePointDensity();
        
        evaluatePointDensity();
    }
    
    private void calculatePointDensity() {
        if (sunCoeff != null) {
            pointDensity = (sunPoints.size()/(2*Math.PI*sunCoeff[2]*sunCoeff[3]));
        }
    }
    
    private void evaluatePointDensity() {
        
        if (pointDensity < 0.5) {
            sunPoints.clear();
        } else if (sunCoeff == null) {
            sunPoints.clear();
        } else if ((sunCoeff[2]/sunCoeff[3])> 6) {
            sunPoints.clear();
        }
    }
    
    private void removePointsUnderSkyline(Set<PairInt> skyline, 
        Set<PairInt> sunPoints, Set<PairInt> skyPoints, 
        int width, int height) {
                
        //TODO:  needs alot of testing with images that are tilted and
        // have sun in them.
        
        int[] dSkylineXY = determineChangeTowardsSkyline(skyline, sunPoints,
            skyPoints);
        
        int dSkylineX = Math.round(dSkylineXY[0]);
        int dSkylineY = Math.round(dSkylineXY[1]);
        
        if ((dSkylineX == 0) && (dSkylineY == 0)) {
            // do nothing
            
        } else if ((dSkylineX >= 0) && (dSkylineY >= 0)) {
            
            for (PairInt s : skyline) {
                int vX = s.getX() + dSkylineX;
                int vY = s.getY() + dSkylineY;
                while ((vX < (width - 1)) && (vY < (height - 1))) {
                    PairInt p = new PairInt(vX, vY);
                    skyPoints.remove(p);
                    vX += dSkylineX;
                    vY += dSkylineY;
                }
            }
            
        } else if ((dSkylineX == -1) && (dSkylineY >= 0)) {
            
            for (PairInt s : skyline) {
                int vX = s.getX() + dSkylineX;
                int vY = s.getY() + dSkylineY;
                while ((vX > -1) && (vY < (height - 1))) {
                    PairInt p = new PairInt(vX, vY);
                    skyPoints.remove(p);
                    vX += dSkylineX;
                    vY += dSkylineY;
                }
            }
        
        } else if ((dSkylineX >= 0) && (dSkylineY == -1)) {
            
            for (PairInt s : skyline) {
                int vX = s.getX() + dSkylineX;
                int vY = s.getY() + dSkylineY;
                while ((vX < (width - 1)) && (vY > -1)) {
                    PairInt p = new PairInt(vX, vY);
                    skyPoints.remove(p);
                    vX += dSkylineX;
                    vY += dSkylineY;
                }
            }
            
        } else if ((dSkylineX == -1) && (dSkylineY == -1)) {
            
            for (PairInt s : skyline) {
                int vX = s.getX() + dSkylineX;
                int vY = s.getY() + dSkylineY;
                while ((vX > -1) && (vY > -1)) {
                    PairInt p = new PairInt(vX, vY);
                    skyPoints.remove(p);
                    vX += dSkylineX;
                    vY += dSkylineY;
                }
            }
        
        }
    }
    
    private int[] determineChangeTowardsSkyline(Set<PairInt> skyline, 
        Set<PairInt> sunPoints, Set<PairInt> skyPoints) {
        
        //TODO:  needs alot of testing with images that are tilted and
        // have sun in them.
        
        // need to fit a line to skyline to get the slope and hence
        // know what is perpendicular to the skyline near the sun.
        
        int[] xsl = new int[skyline.size()];
        int[] ysl = new int[xsl.length];
        int i = 0;
        for (PairInt p : skyline) {
            xsl[i] = p.getX();
            ysl[i] = p.getY();
            i++;
        }
        
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
            
        double[] skyXYCen = ch.calculateXYCentroids(skyPoints);
            
        double[] sunXYCen = ch.calculateXYCentroids(sunPoints);
        
        LinearRegression lr = new LinearRegression();
        float[] yInterceptAndSlope = lr.calculateTheilSenEstimatorParams(
            xsl, ysl);
        
        int dSkylineX = 0;
        int dSkylineY = 0;
        
        /*
        inf
           1.7
         |    1.0
         |   /    0.57
         | /     
        @ ----- 0        
        */
        
        if (Math.abs(yInterceptAndSlope[1]) < 0.57) {
            // ~ horizontal line, less than 30 degrees from horiz
            if (skyXYCen[1] > sunXYCen[1]) {
                dSkylineY = -1;
            } else {
                dSkylineY = 1;
            }
        } else if (Math.abs(yInterceptAndSlope[1]) < 1.73) {
            // ~ diagonal line, between 30 degrees and 60 degrees from horiz
            if (skyXYCen[0] > sunXYCen[0]) {
                dSkylineX = -1;
            } else {
                dSkylineX = 1;
            }
            if (skyXYCen[1] > sunXYCen[1]) {
                dSkylineY = -1;
            } else {
                dSkylineY = 1;
            }
        } else  {
            // approx w/ a vertical line
            if (skyXYCen[0] > sunXYCen[0]) {
                dSkylineX = -1;
            } else {
                dSkylineX = 1;
            }
        }
        
        return new int[]{dSkylineX, dSkylineY};
    }
    
    public void correctSkylineForSun( 
        Set<PairInt> skyPoints, Image colorImg, int xOffset, int yOffset, 
        GreyscaleImage gradientXY) {
        
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        double[] xySunCen = ch.calculateXYCentroids(sunPoints);
        
        int h = colorImg.getHeight();
        int w = colorImg.getWidth();
        
        int srchRadius = 125;
        
        int maxValue = Integer.MIN_VALUE;
        
        Set<PairInt> pointsNearSun = new HashSet<PairInt>();
                
        for (int col = ((int)xySunCen[0] - srchRadius); 
            col < ((int)xySunCen[0] + srchRadius); col++) {
        
            for (int row = ((int)xySunCen[1] - srchRadius); 
                row < ((int)xySunCen[1] + srchRadius); row++) {
                
                if ((col < 0) || (col > (w - 1)) || (row < 0) || (row > (h - 1))) {
                    continue;
                }
                
                int v = gradientXY.getValue(col, row);
                
                if (v > maxValue) {
                    maxValue = v;
                }
                
                if (v > 0) {
                    pointsNearSun.add(new PairInt(col, row));
                }
            }
        }
        
        PairIntArray valueCounts = 
            Histogram.createADescendingSortByKeyArray(pointsNearSun, gradientXY);
        int countIsDoubledValue = maxValue;
        int lastCount = Integer.MAX_VALUE;
        for (int i = 1; i < valueCounts.getN(); i++) {
            int count = valueCounts.getY(i);
            if ((count/lastCount) >= 2) {
                countIsDoubledValue = valueCounts.getX(i);
                break;
            }
            lastCount = count;
        }
        
        assert(maxValue == valueCounts.getX(0));
        
        int valueTolerance = 10;
        if ((maxValue - countIsDoubledValue) > valueTolerance) {
            valueTolerance = maxValue - countIsDoubledValue;
        }
        
        int minX = Integer.MAX_VALUE;
        int minY = Integer.MAX_VALUE;
        int maxX = Integer.MIN_VALUE;
        int maxY = Integer.MIN_VALUE;
        
        Set<PairInt> set = new HashSet<PairInt>();
        
        for (int col = ((int)xySunCen[0] - srchRadius); 
            col < ((int)xySunCen[0] + srchRadius); col++) {
        
            for (int row = ((int)xySunCen[1] - srchRadius); 
                row < ((int)xySunCen[1] + srchRadius); row++) {
                
                if ((col < 0) || (col > (w - 1)) || (row < 0) || (row > (h - 1))) {
                    continue;
                }
                
                int v = gradientXY.getValue(col, row);
                
                if (Math.abs(v - maxValue) < valueTolerance) {
                    
                    set.add(new PairInt(col, row));
                    
                    if (col < minX) {
                        minX = col;
                    }
                    if (col > maxX) {
                        maxX = col;
                    }
                    if (row < minY) {
                        minY = row;
                    }
                    if (row > maxY) {
                        maxY = row;
                    }
                }
            }
        }

        //debugPlot(set, colorImg, xOffset, yOffset, 
        //    "horizon_near_sun_before_thinning");
        
        // points in set now represent the skyline near the sun.
        // the points are a line widened by convolution so need to be thinned
        // to the line centroid.
        // ErosionFilter isn't currently able to provide a line that is centered,
        // it may be a line closer to one side than the other of the original
        // thick band of points.
        
        // so trying Zhang-Suen: 
        
        ZhangSuenLineThinner zsLT = new ZhangSuenLineThinner();
        zsLT.applyLineThinner(set, minX, maxX, minY, maxY);
        
        ch.populateGapsWithInterpolation(set);

        //debugPlot(set, colorImg, xOffset, yOffset,
        //    "horizon_near_sun");

        // if know skyline direction and any points are "below" set towards
        // the skyline, those should be removed from skyPoints
        
        removePointsUnderSkyline(set, sunPoints, skyPoints, w, h);

        skyPoints.addAll(set);
        skyPoints.addAll(sunPoints);
        
        // fill in the gaps
        //Set<PairInt> embeddedPoints = findEmbeddedNonPoints(skyPoints);
        //skyPoints.addAll(embeddedPoints);
        
        //debugPlot(skyPoints, colorImg, xOffset, yOffset,
        //    "horizon_near_sun_final");
    }
    
    
    public double getPointDensity() {
        return pointDensity;
    }
    
    public double[] getSunEllipseCoeff() {
        return sunCoeff;
    }
    
    public Set<PairInt> getSunPoints() {
        return sunPoints;
    }
}
