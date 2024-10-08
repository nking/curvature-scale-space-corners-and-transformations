package algorithms.imageProcessing.transform;

import algorithms.util.AngleUtil;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.list.TIntList;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import thirdparty.HungarianAlgorithm;

/**
 *
 * @author nichole
 */
public class MatchedPointsTransformationCalculator {

    protected transient Logger log = Logger.getLogger(
        MatchedPointsTransformationCalculator.class.getName());

    private boolean debug = false;

    public void useDebugMode() {
        debug = true;
    }

    /**
     * coordinate transformations from image 1 to image 2 are calculated from
     * matching lists of x, y coordinates, and given "scale" as a starting
     * parameter.  Scale is determined roughly from the contour matcher,
     * so can be used to get a rough first solution.
     * Note, the rotation, when applied, will result in a clockwise
      * direction (which is in the -z direction using right hand rule).
     *<pre>
     * positive Y is up
       positive X is right
       positive theta starts from Y=0, X>=0 and proceeds CW
             +Y 270
                 |
                 |
          180--------- 0   +X
                 |
                 |
                 90
                 -Y
     * </pre>
     * @param scale
     * @param matchedXY1
     * @param matchedXY2
     * @param centroidX1
     * @param centroidY1
     * @return 
     */
    public TransformationParameters calulateEuclideanGivenScale(
        double scale, PairIntArray matchedXY1, PairIntArray matchedXY2,
        double centroidX1, double centroidY1) {

        if (matchedXY1 == null) {
            throw new IllegalArgumentException("matchedXY1 cannot be null");
        }
        if (matchedXY2 == null) {
            throw new IllegalArgumentException("matchedXY2 cannot be null");
        }
        if (matchedXY1.getN() != matchedXY2.getN()) {
            throw new IllegalArgumentException(
                "matchedXY1 and matchedXY2 must have same number of points");
        }
        if (matchedXY1.getN() < 2) {
            return null;
        }

        log.info("start solution for " + matchedXY1.getN() + " points");

        /*
        solve for rotation.

        Take the same 2 pairs int both images and get the difference in their
        angles:
            tan(theta) = y / x

        For example:
            theta of pair in image1:
                theta = math.atan2( (y1-y0), (x1-x0) )
                      = 0.7853981633974483 radians
                      = 45 degrees

            theta of pair in image2:
                theta = math.atan2( (yt1-yt0), (xt1-xt0) )
                      = 0.3490522203358645
                      = 20.0

            rotation = theta_image1 - theta_image2 = 25 degrees
        */

        /*
        discard outside avg +- stdev
        */

        AngleUtil angleUtil = new AngleUtil();

        double[] thetas = new double[matchedXY1.getN()];
        double[] scales = new double[matchedXY1.getN()];
        double thetaSum = 0;
        double scaleSum = 0;
        for (int i = 0; i < matchedXY1.getN(); i++) {
            int x0im1 = matchedXY1.getX(i);
            int y0im1 = matchedXY1.getY(i);
            int x0im2 = matchedXY2.getX(i);
            int y0im2 = matchedXY2.getY(i);
            int x1im1, y1im1, x1im2, y1im2;
            if ((i + 1) == matchedXY1.getN()) {
                x1im1 = matchedXY1.getX(0);
                y1im1 = matchedXY1.getY(0);
                x1im2 = matchedXY2.getX(0);
                y1im2 = matchedXY2.getY(0);
            } else {
                x1im1 = matchedXY1.getX(i + 1);
                y1im1 = matchedXY1.getY(i + 1);
                x1im2 = matchedXY2.getX(i + 1);
                y1im2 = matchedXY2.getY(i + 1);
            }
            double diffX1 = (x1im1 - x0im1);
            double diffY1 = (y1im1 - y0im1);

            double diffX2 = (x1im2 - x0im2);
            double diffY2 = (y1im2 - y0im2);

            double t = angleUtil.subtract(diffX1, diffY1, diffX2, diffY2);

            t *= -1;
            
            if (t < 0) {
                t = 2 * Math.PI + t;
            }
            
            thetas[i] = t;

            thetaSum += thetas[i];

            double lenim1 = Math.sqrt(Math.pow(diffX1, 2)
                + Math.pow(diffY1, 2));
            double lenim2 = Math.sqrt(Math.pow(diffX2, 2)
                + Math.pow(diffY2, 2));
            scales[i] = lenim2/lenim1;
            scaleSum += scales[i];
        }

        double avgScale = scaleSum / (double)matchedXY1.getN();
        double avgTheta = thetaSum / (double)matchedXY1.getN();
        double stDevThetaSum = 0;
        double stDevScaleSum = 0;
        for (int i = 0; i < matchedXY1.getN(); i++) {
            stDevThetaSum += Math.pow(thetas[i] - avgTheta, 2);
            stDevScaleSum += Math.pow(scales[i] - avgScale, 2);
        }
        double stDevTheta= Math.sqrt(stDevThetaSum/(matchedXY1.getN() - 1));
        double stDevScale= Math.sqrt(stDevScaleSum/(matchedXY1.getN() - 1));

        double rotSum = 0;
        double rCount = 0;
        scaleSum = 0;
        double sCount = 0;

        List<Integer> rm = new ArrayList<Integer>();
        for (int i = 0; i < matchedXY1.getN(); i++) {
            double dss = Math.abs(scales[i] - avgScale);
            if (dss > 1.5*stDevScale) {
                rm.add(Integer.valueOf(i));
                continue;
            }
            double dtt = Math.abs(thetas[i] - avgTheta);
            if (dtt > 1.5*stDevTheta) {
                rm.add(Integer.valueOf(i));
                continue;
            }

log.info("scl=" + scales[i] + " stDevScale=" + stDevScale
+ " abs(scale-avg)=" + dss);

            scaleSum += scales[i];
            sCount++;

log.info("rot=" + thetas[i] + " stDevTheta=" + stDevTheta
+ " abs(theta-avg)=" + dtt);

            rotSum += thetas[i];
            rCount++;
        }

        if (!rm.isEmpty()) {

            PairIntArray xy1 = new PairIntArray();
            PairIntArray xy2 = new PairIntArray();

            for (int i = 0; i < matchedXY1.getN(); i++) {
                if (rm.contains(Integer.valueOf(i))) {
                    continue;
                }
                xy1.add(matchedXY1.getX(i), matchedXY1.getY(i));
                xy2.add(matchedXY2.getX(i), matchedXY2.getY(i));
            }

            return calulateEuclideanGivenScale(scale, xy1,
                xy2, centroidX1, centroidY1);
        }

        double theRotation = rotSum/rCount;
        double theScale = scaleSum/sCount;

        log.info("given scale=" + scale + " found scale=" + theScale);
        log.info("rotation = " + theRotation);

        if (Math.abs(theScale - scale) > scale*0.1) {
            log.warning("the differences in estimated scale and given scale are"
                + " large.  this can happen if the given scale value was "
                + " determined from contour matching."
                + " the estimate here uses pairs of points which may"
                + " be close to one another.  choosing the scale given to "
                + " the method and continuing.");
        }

        theScale = scale;

        /*
        estimate translation:

        transX = xt0 - 
            (xc*scale + (((x0-xc)*scale*math.cos(theta))
            + ((y0-yc)*scale*math.sin(theta)))

        transY = yt0 - 
            (yc*scale + ((-(x0-xc)*scale*math.sin(theta))
            + ((y0-yc)*scale*math.cos(theta)))
        */
        double mc = Math.cos(theRotation);
        double ms = Math.sin(theRotation);
        double transXSum = 0;
        double transYSum = 0;
        for (int i = 0; i < matchedXY1.getN(); i++) {
            int xim1 = matchedXY1.getX(i);
            int yim1 = matchedXY1.getY(i);
            int xim2 = matchedXY2.getX(i);
            int yim2 = matchedXY2.getY(i);

            double trX1 = centroidX1*scale + ((xim1 - centroidX1) * scale*mc)
                + ((yim1 - centroidY1) *scale*ms);

            double trY1 = centroidY1*scale + (-(xim1 - centroidX1) * scale*ms)
                + ((yim1 - centroidY1) * scale*mc);

            double transX = xim2 - trX1;

            double transY = yim2 - trY1;

            transXSum += transX;
            transYSum += transY;
        }
        double theTranslationX = transXSum/(double)matchedXY1.getN();
        double theTranslationY = transYSum/(double)matchedXY1.getN();

        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians((float)theRotation);
        params.setScale((float)theScale);
        params.setTranslationX((float)theTranslationX);
        params.setTranslationY((float)theTranslationY);
        params.setOriginX((float)centroidX1);
        params.setOriginY((float)centroidY1);
        params.setNumberOfPointsUsed(matchedXY1.getN());

        if (debug) {
            log.info("params: " + params.toString());
        }

        return params;
    }

     /**
     * coordinate transformations from pair 1 to pair 2 are calculated.
     *
     * positive Y is up
       positive X is right
       positive theta starts from Y=0, X>=0 and proceeds CW
                270
                 |
                 |
          180--------- 0   +X
                 |
                 |
                 90
                 -Y
     * </pre>
     * @param centroidX1
     * @param centroidY1
     * @return
     */
    public TransformationParameters calulateEuclidean(
        final int set1X1, final int set1Y1,
        final int set1X2, final int set1Y2,
        final int set2X1, final int set2Y1,
        final int set2X2, final int set2Y2,
        final double centroidX1, final double centroidY1) {

        /*
        solve for rotation.

        Take the same 2 pairs int both images and get the difference in their
        angles:
            tan(theta) = y / x

        For example:
            theta of pair in image1:
                theta = math.atan( (y1-y0)/(x1-x0) )
                      = 0.7853981633974483 radians
                      = 45 degrees

            theta of pair in image2:
                theta = math.atan( (yt1-yt0)/(xt1-xt0) )
                      = 0.3490522203358645
                      = 20.0

            rotation = theta_image1 - theta_image2 = 25 degrees
        */

        AngleUtil angleUtil = new AngleUtil();

        double dx1 = set1X1 - set1X2;
        double dy1 = set1Y1 - set1Y2;

        double dx2 = set2X1 - set2X2;
        double dy2 = set2Y1 - set2Y2;

        double theta = angleUtil.subtract(dx1, dy1, dx2, dy2);
        
        theta *= -1;
            
        if (theta < 0) {
            theta += (2 * Math.PI);
        }

        double sep1 = Math.sqrt((dx1*dx1) + (dy1*dy1));

        double sep2 = Math.sqrt((dx2*dx2) + (dy2*dy2));

        if (sep1 == 0 || sep2 == 0) {
            // the 2 points in set1 or set2 are the same point
            return null;
        }
        
        double scale = sep2/sep1;
        
        /*
        estimate translation:

        xr_0 = xc*scale + (((x0-xc)*scale*math.cos(theta)) - ((y0-yc)*scale*math.sin(theta)))

        xt_0 = xr_0 + transX = x1

        yr_0 = yc*scale + (((x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))

        yt_0 = yr_0 + transY = y1
        */
        double mc = Math.cos(theta);
        double ms = Math.sin(theta);

        double tr1X1 = centroidX1*scale + ((set1X1 - centroidX1) * scale * mc)
            + ((set1Y1 - centroidY1) *scale*ms);

        double tr1Y1 = centroidY1*scale + (-(set1X1 - centroidX1) * scale * ms)
            + ((set1Y1 - centroidY1) *scale*mc);

        double tr1X2 = centroidX1*scale + ((set1X2 - centroidX1) * scale*mc)
            + ((set1Y2 - centroidY1) *scale*ms);

        double tr1Y2 = centroidY1*scale + (-(set1X2 - centroidX1) *scale*ms)
            + ((set1Y2 - centroidY1) *scale*mc);

        double transX1 = (set2X1 - tr1X1);
        double transX2 = (set2X2 - tr1X2);
        double transY1 = (set2Y1 - tr1Y1);
        double transY2 = (set2Y2 - tr1Y2);

        double transX = 0.5 * (transX1 + transX2);

        double transY = 0.5 * (transY1 + transY2);

        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians((float)theta);
        params.setScale((float)scale);
        params.setTranslationX((float)transX);
        params.setTranslationY((float)transY);
        params.setOriginX((float)centroidX1);
        params.setOriginY((float)centroidY1);

        if (debug) {
            log.info("params: " + params.toString());
        }

        return params;
    }
    
    /**
     * coordinate transformations from pair 1 to pair 2 are calculated from
     * the widest pairings of the given matched points.
     *
     * positive Y is up
       positive X is right
       positive theta starts from Y=0, X>=0 and proceeds CW
                270
                 |
                 |
          180--------- 0   +X
                 |
                 |
                 90
                 -Y
     * </pre>
     * @param matchedXY1
     * @param matchedXY2
     * @param weights
     * @param centroidX1
     * @param centroidY1
     * @param outputScaleRotTransXYStDev output standard deviation of
     * the scale, rotation, and translations in X and Y. If null, calculations
     * are not performed for standard deviation from mean.
     * @return
     */
    public TransformationParameters calulateEuclidean(
        PairIntArray matchedXY1, PairIntArray matchedXY2, float[] weights,
        final double centroidX1, final double centroidY1,
        float[] outputScaleRotTransXYStDev) {
     
        TIntList inlierIndexes = null;
        
        return calulateEuclidean(matchedXY1, matchedXY2, weights,
            centroidX1, centroidY1, outputScaleRotTransXYStDev, inlierIndexes);
    }

    /**
     * coordinate transformations from pair 1 to pair 2 are calculated from
     * the widest pairings of the given matched points.
     *
     * positive Y is up
       positive X is right
       positive theta starts from Y=0, X>=0 and proceeds CW
                270
                 |
                 |
          180--------- 0   +X
                 |
                 |
                 90
                 -Y
     * </pre>
     * @param matchedXY1
     * @param matchedXY2
     * @param weights
     * @param centroidX1
     * @param centroidY1
     * @param outputScaleRotTransXYStDev output standard deviation of
     * the scale, rotation, and translations in X and Y. If null, calculations
     * are not performed for standard deviation from mean.
     * @return
     */
    public TransformationParameters calulateEuclidean(
        PairIntArray matchedXY1, PairIntArray matchedXY2, float[] weights,
        final double centroidX1, final double centroidY1,
        float[] outputScaleRotTransXYStDev, TIntList outputInlierIndexes) {
        
        if (matchedXY1 == null) {
            throw new IllegalArgumentException("matchedXY1 cannot be null");
        }
        if (matchedXY2 == null) {
            throw new IllegalArgumentException("matchedXY2 cannot be null");
        }
        if (weights == null) {
            throw new IllegalArgumentException("weights cannot be null");
        }
        if (matchedXY1.getN() != matchedXY2.getN()) {
            throw new IllegalArgumentException(
            "matchedXY1 must be the same length as matchedXY2");
        }
        if (matchedXY1.getN() != weights.length) {
            throw new IllegalArgumentException(
            "weights must be the same length as matchedXY1");
        }
        if (matchedXY1.getN() < 2) {
            throw new IllegalArgumentException(
            "matchedXY1 must have at least 2 points");
        }
        if ((outputScaleRotTransXYStDev != null) && (outputScaleRotTransXYStDev.length < 4)) {
            throw new IllegalArgumentException(
            "outputScaleRotTransXYStDev has to be at least 4 in length if not null");
        }
       
        /*
        solve for rotation.

        Take the same 2 pairs int both images and get the difference in their
        angles:
            tan(theta) = y / x

        For example:
            theta of pair in image1:
                theta = math.atan( (y1-y0)/(x1-x0) )
                      = 0.7853981633974483 radians
                      = 45 degrees

            theta of pair in image2:
                theta = math.atan( (yt1-yt0)/(xt1-xt0) )
                      = 0.3490522203358645
                      = 20.0

            rotation = theta_image1 - theta_image2 = 25 degrees
        */
            
        //choosing pairs of points by pairing for maximum distance.
        
        Set<PairInt> pairIndexesSet = new HashSet<PairInt>();
        
        if (matchedXY1.getN() < 100) {
            useBipartiteMatching(matchedXY1, pairIndexesSet);
        } else {
            useGreedyMatching(matchedXY1, pairIndexesSet);
        }
        
        List<PairInt> pairIndexes = new ArrayList<PairInt>(pairIndexesSet);
        
        AngleUtil angleUtil = new AngleUtil();

        float totW = 0;
        for (PairInt pairIndex : pairIndexes) {
            int idx1 = pairIndex.getX();
            int idx2 = pairIndex.getY();            
            float pairWeight = (weights[idx1] + weights[idx2]);
            totW += pairWeight;
        }
        Map<PairInt, Float> pairWeights = new HashMap<PairInt, Float>();
        for (PairInt pairIndex : pairIndexes) {
            int idx1 = pairIndex.getX();
            int idx2 = pairIndex.getY();            
            float pairWeight = (weights[idx1] + weights[idx2])/totW;
            pairWeights.put(pairIndex, Float.valueOf(pairWeight));
        }
        
        List<Double> thetas = new ArrayList<Double>();        
        List<Float> scales = new ArrayList<Float>();
        List<Float> transXs = new ArrayList<Float>();
        List<Float> transYs = new ArrayList<Float>();
        List<Float> pairWeights2 = new ArrayList<Float>();

        for (PairInt pairIndex : pairIndexes) {

            int idx1 = pairIndex.getX();
            int idx2 = pairIndex.getY();

            int set1X1 = matchedXY1.getX(idx1);
            int set1Y1 = matchedXY1.getY(idx1);
            int set1X2 = matchedXY1.getX(idx2);
            int set1Y2 = matchedXY1.getY(idx2);

            int set2X1 = matchedXY2.getX(idx1);
            int set2Y1 = matchedXY2.getY(idx1);
            int set2X2 = matchedXY2.getX(idx2);
            int set2Y2 = matchedXY2.getY(idx2);

            double dx1 = set1X1 - set1X2;
            double dy1 = set1Y1 - set1Y2;

            double dx2 = set2X1 - set2X2;
            double dy2 = set2Y1 - set2Y2;
            
            double theta = angleUtil.subtract(dx1, dy1, dx2, dy2);

            theta *= -1;
            
            if (theta < 0) {
                theta += 2 * Math.PI;
            }
            
            double sep1 = Math.sqrt((dx1*dx1) + (dy1*dy1));

            double sep2 = Math.sqrt((dx2*dx2) + (dy2*dy2));

            double scale = sep2/sep1;
            
            /*
            estimate translation:

            xr_0 = xc*scale + (((x0-xc)*scale*math.cos(theta)) - ((y0-yc)*scale*math.sin(theta)))

            xt_0 = xr_0 + transX = x1

            yr_0 = yc*scale + (((x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))

            yt_0 = yr_0 + transY = y1
            */
            double mc = Math.cos(theta);
            double ms = Math.sin(theta);

            double tr1X1 = centroidX1 * scale + (((set1X1 - centroidX1) * scale * mc)
                + ((set1Y1 - centroidY1) * scale * ms));

            double tr1Y1 = centroidY1 * scale + (-(set1X1 - centroidX1) * scale * ms)
                + ((set1Y1 - centroidY1) * scale * mc);

            double tr1X2 = centroidX1 * scale + (((set1X2 - centroidX1) * scale * mc)
                + ((set1Y2 - centroidY1) * scale * ms));

            double tr1Y2 = centroidY1 * scale + (-(set1X2 - centroidX1) * scale * ms)
                + ((set1Y2 - centroidY1) * scale * mc);

            double transX1 = (set2X1 - tr1X1);
            double transX2 = (set2X2 - tr1X2);
            double transY1 = (set2Y1 - tr1Y1);
            double transY2 = (set2Y2 - tr1Y2);

            double transX = 0.5 * (transX1 + transX2);

            double transY = 0.5 * (transY1 + transY2);
                        
            thetas.add(Double.valueOf(theta));
            pairWeights2.add(pairWeights.get(pairIndex));
            scales.add((float)scale);
            transXs.add((float)transX);
            transYs.add((float)transY);
        }
                
        boolean useRadians = true;
        List<Double> thetaCorr = new ArrayList<Double>(thetas); 
        
        AngleUtil.correctForQuadrants(thetaCorr, useRadians);
        
        double thetaAvg = 0;
        for (int i = 0; i < thetaCorr.size(); ++i) {
            thetaAvg += thetaCorr.get(i) * pairWeights2.get(i);
        }
        
        // remove discrepant scale ?  defeats use of weights?
        float unwScaleAvg = 0;
        for (Float s : scales) {
            unwScaleAvg += s.floatValue();
        }
        unwScaleAvg /= (float)scales.size();
        double unwScaleStDev = 0;
        for (Float s : scales) {
            float diff = s.floatValue() - unwScaleAvg;
            unwScaleStDev += (diff * diff);
        }
        unwScaleStDev = Math.sqrt(unwScaleStDev/((double)scales.size() - 1.));
        
        if (unwScaleStDev > (0.2*unwScaleAvg)) {
            // filter using the median of all possible differences of scale
            // this requires O(N^2/2) in space
            List<Integer> remove = findOutliersUsingMedianOfDifferences(scales);
          
            if (!remove.isEmpty()) {
                for (int i = (remove.size() - 1); i > -1; --i) {
                    int idx = remove.get(i).intValue();
                    scales.remove(idx);
                    thetaCorr.remove(idx);
                    transXs.remove(idx);
                    transYs.remove(idx);
                    pairWeights2.remove(idx);
                    pairIndexes.remove(idx);
                }
                // rebalance the remaining weights
                float ws = 0;
                for (Float pw : pairWeights2) {
                    ws += pw.floatValue();
                }
                for (int i = 0; i < pairWeights2.size(); ++i) {
                    float v = pairWeights2.get(i)/ws;
                    pairWeights2.set(i, v);
                }
                AngleUtil.correctForQuadrants(thetaCorr, useRadians);
            }
        }
        
        double tolTheta = 20. * Math.PI/180.;
        
        // remove those w/ discrepant thetas
        List<Integer> remove = new ArrayList<Integer>();
        for (int i = 0; i < thetaCorr.size(); ++i) {
            double diffT = Math.abs(thetaCorr.get(i) - thetaAvg);
            double diffS = Math.abs(scales.get(i) - unwScaleAvg);
            if (diffT > tolTheta) {
                remove.add(Integer.valueOf(i));
            } else if (diffS > (1.5*unwScaleStDev)) {
                remove.add(Integer.valueOf(i));
            }
        }
        
        if (remove.size() == scales.size()) {
            return null;
        }
        
        if (!remove.isEmpty()) {
            for (int i = (remove.size() - 1); i > -1; --i) {
                int idx = remove.get(i).intValue();
                scales.remove(idx);
                thetaCorr.remove(idx);
                transXs.remove(idx);
                transYs.remove(idx);
                pairWeights2.remove(idx);
                pairIndexes.remove(idx);
            }
            // rebalance the remaining weights
            float ws = 0;
            for (Float pw : pairWeights2) {
                ws += pw.floatValue();
            }
            for (int i = 0; i < pairWeights2.size(); ++i) {
                float v = pairWeights2.get(i)/ws;
                pairWeights2.set(i, v);
            }
            AngleUtil.correctForQuadrants(thetaCorr, useRadians);
        }
        
        float scaleAvg = 0;
        thetaAvg = 0;
        float transXAvg = 0;
        float transYAvg = 0;
        for (int i = 0; i < scales.size(); ++i) {
            float w = pairWeights2.get(i);
            scaleAvg += scales.get(i) * w;
            thetaAvg += thetaCorr.get(i) * w;
            transXAvg += transXs.get(i) * w;
            transYAvg += transYs.get(i) * w;
        }
            
        // ----- determine standard deviations ----
        if (outputScaleRotTransXYStDev != null) {
            Arrays.fill(outputScaleRotTransXYStDev, Float.MAX_VALUE);
            double scaleSum = 0;
            double thetaSum = 0;
            double transXSum = 0;
            double transYSum = 0;
            for (int i = 0; i < scales.size(); ++i) {
                double diff = scales.get(i) - scaleAvg;
                scaleSum += (diff * diff);
                
                diff = thetaCorr.get(i) - thetaAvg;
                thetaSum += (diff * diff);
                
                diff = transXs.get(i) - transXAvg;
                transXSum += (diff * diff);
                
                diff = transYs.get(i) - transYAvg;
                transYSum += (diff * diff);
            }
            double scaleStDv = Math.sqrt(scaleSum/((double)scales.size() - 1.));
            double thetaStDv = Math.sqrt(thetaSum/((double)scales.size() - 1.));
            double transXStDv = Math.sqrt(transXSum/((double)scales.size() - 1.));
            double transYStDv = Math.sqrt(transYSum/((double)scales.size() - 1.));
            
            outputScaleRotTransXYStDev[0] = (float)scaleStDv;
            outputScaleRotTransXYStDev[1] = (float)thetaStDv;
            outputScaleRotTransXYStDev[2] = (float)transXStDv;
            outputScaleRotTransXYStDev[3] = (float)transYStDv;
        }
        
        if (thetaAvg >= (2. * Math.PI)) {
            thetaAvg -= (2. * Math.PI);
        }          
                
        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians((float)thetaAvg);
        params.setScale(scaleAvg);
        params.setTranslationX(transXAvg);
        params.setTranslationY(transYAvg);
        params.setOriginX((float)centroidX1);
        params.setOriginY((float)centroidY1);
        if (outputScaleRotTransXYStDev != null) {
            params.setStandardDeviations(outputScaleRotTransXYStDev);
        }
        // this is the number of data points as the number of pairs
        params.setNumberOfPointsUsed(scales.size());
     
        if (outputInlierIndexes != null) {
            //remaining pairIndexes refer to the input matched correspondence.
            // x is with respect to matchedXY1
            // y is with respect to matchedXY2
            // and they're parallel lists, so they're the same.
            TIntSet idxs = new TIntHashSet();
            for (PairInt pIdxs : pairIndexes) {
                idxs.add(pIdxs.getX());
                idxs.add(pIdxs.getY());
            }
           
            outputInlierIndexes.clear();
            outputInlierIndexes.addAll(idxs);
            outputInlierIndexes.sort();
        }
        
        if (debug) {
            log.info("params: " + params.toString());
        }

        return params;
    }

    /**
     * coordinate transformations from pair 1 to pair 2 are calculated from
     * the widest pairings of the given matched points.  No standard 
     * deviations from mean values are calculated and no filtering is
     * performed.  The method is meant to be used with code which will use
     * other methods to filter the points.
     *
     * positive Y is up
       positive X is right
       positive theta starts from Y=0, X>=0 and proceeds CW
                270
                 |
                 |
          180--------- 0   +X
                 |
                 |
                 90
                 -Y
     * </pre>
     * @param matchedXY1
     * @param matchedXY2
     * @param centroidX1
     * @param centroidY1
     * @return
     */
    public TransformationParameters calulateEuclideanWithoutFilter(
        PairIntArray matchedXY1, PairIntArray matchedXY2,
        final double centroidX1, final double centroidY1) {
        
        if (matchedXY1 == null) {
            throw new IllegalArgumentException("matchedXY1 cannot be null");
        }
        if (matchedXY2 == null) {
            throw new IllegalArgumentException("matchedXY2 cannot be null");
        }
        if (matchedXY1.getN() != matchedXY2.getN()) {
            throw new IllegalArgumentException(
            "matchedXY1 must be the same length as matchedXY2");
        }
        if (matchedXY1.getN() < 2) {
            throw new IllegalArgumentException(
            "matchedXY1 must have at least 2 points");
        }
       
        /*
        solve for rotation.

        Take the same 2 pairs int both images and get the difference in their
        angles:
            tan(theta) = y / x

        For example:
            theta of pair in image1:
                theta = math.atan( (y1-y0)/(x1-x0) )
                      = 0.7853981633974483 radians
                      = 45 degrees

            theta of pair in image2:
                theta = math.atan( (yt1-yt0)/(xt1-xt0) )
                      = 0.3490522203358645
                      = 20.0

            rotation = theta_image1 - theta_image2 = 25 degrees
        */
            
        //choosing pairs of points by pairing for maximum distance.
        Set<PairInt> pairIndexes = new HashSet<PairInt>();     
        useGreedyMatching(matchedXY1, pairIndexes);
        
        AngleUtil angleUtil = new AngleUtil();
        
        double[] thetas = new double[pairIndexes.size()];

        float scaleSum = 0;
        float transXSum= 0;
        float transYSum = 0;
        int count = 0;
        for (PairInt pairIndex : pairIndexes) {

            int idx1 = pairIndex.getX();
            int idx2 = pairIndex.getY();

            int set1X1 = matchedXY1.getX(idx1);
            int set1Y1 = matchedXY1.getY(idx1);
            int set1X2 = matchedXY1.getX(idx2);
            int set1Y2 = matchedXY1.getY(idx2);

            int set2X1 = matchedXY2.getX(idx1);
            int set2Y1 = matchedXY2.getY(idx1);
            int set2X2 = matchedXY2.getX(idx2);
            int set2Y2 = matchedXY2.getY(idx2);

            double dx1 = set1X1 - set1X2;
            double dy1 = set1Y1 - set1Y2;

            double dx2 = set2X1 - set2X2;
            double dy2 = set2Y1 - set2Y2;
            
            double theta = angleUtil.subtract(dx1, dy1, dx2, dy2);

            theta *= -1;
            
            if (theta < 0) {
                theta += 2 * Math.PI;
            }
            
            double sep1 = Math.sqrt((dx1*dx1) + (dy1*dy1));

            double sep2 = Math.sqrt((dx2*dx2) + (dy2*dy2));

            double scale = sep2/sep1;
            
            /*
            estimate translation:

            xr_0 = xc*scale + (((x0-xc)*scale*math.cos(theta)) - ((y0-yc)*scale*math.sin(theta)))

            xt_0 = xr_0 + transX = x1

            yr_0 = yc*scale + (((x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))

            yt_0 = yr_0 + transY = y1
            */
            double mc = Math.cos(theta);
            double ms = Math.sin(theta);

            double tr1X1 = centroidX1 * scale + (((set1X1 - centroidX1) * scale * mc)
                + ((set1Y1 - centroidY1) * scale * ms));

            double tr1Y1 = centroidY1 * scale + (-(set1X1 - centroidX1) * scale * ms)
                + ((set1Y1 - centroidY1) * scale * mc);

            double tr1X2 = centroidX1 * scale + (((set1X2 - centroidX1) * scale * mc)
                + ((set1Y2 - centroidY1) * scale * ms));

            double tr1Y2 = centroidY1 * scale + (-(set1X2 - centroidX1) * scale * ms)
                + ((set1Y2 - centroidY1) * scale * mc);

            double transX1 = (set2X1 - tr1X1);
            double transX2 = (set2X2 - tr1X2);
            double transY1 = (set2Y1 - tr1Y1);
            double transY2 = (set2Y2 - tr1Y2);

            double transX = 0.5 * (transX1 + transX2);

            double transY = 0.5 * (transY1 + transY2);
                        
            thetas[count] = theta;
            
            scaleSum += scale;
            transXSum += transX;
            transYSum += transY;
            count++;
        }
        
        float scaleAvg = scaleSum/(float)pairIndexes.size();
        float transXAvg = transXSum/(float)pairIndexes.size();
        float transYAvg = transYSum/(float)pairIndexes.size();
                
        boolean useRadians = true;
        double thetaAvg = AngleUtil.calculateAverageWithQuadrantCorrections(
            thetas, useRadians);
        
        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians((float)thetaAvg);
        params.setScale(scaleAvg);
        params.setTranslationX(transXAvg);
        params.setTranslationY(transYAvg);
        params.setOriginX((float)centroidX1);
        params.setOriginY((float)centroidY1);
        // this is the number of data points as the number of pairs
        params.setNumberOfPointsUsed(matchedXY1.getN());
     
        if (debug) {
            log.info("params: " + params.toString());
        }

        return params;
    }

    public MatchedPointsTransformationCalculator() {
    }

    /**
     * from a set of transformation parameters params that transform
     * points in reference frame 1 into reference frame 2, create
     * a transformation that can transform points in reference frame 2
     * into reference frame 1.
     * @param params transformation parameters to apply to points in reference
     * frame 1 to put them in reference frame 2
     * @return transformation parameters that can transform points in reference
     * frame 2 into reference frame 1
     */
    public TransformationParameters swapReferenceFrames(TransformationParameters
        params) {

        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }

        /*
        xr_0 = xc*scale + (((x0-xc)*scale*math.cos(theta)) + ((y0-yc)*scale*math.sin(theta)))

        xt_0 = xr_0 + transX = x1

        yr_0 = yc*scale + (-((x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))

        yt_0 = yr_0 + transY = y1
        */

        double revRot = -1 * params.getRotationInRadians();
        if (revRot < 0) {
            revRot += 2. * Math.PI;
        }
        double revScale = 1. / params.getScale();

        double transformedXC =
            (params.getOriginX()*params.getScale()) + params.getTranslationX();

        double transformedYC =
            (params.getOriginY()*params.getScale()) + params.getTranslationY();

        TransformationParameters paramsRev = new TransformationParameters();
        paramsRev.setScale((float)revScale);
        paramsRev.setRotationInRadians((float)revRot);
        paramsRev.setTranslationX(0);
        paramsRev.setTranslationY(0);
        paramsRev.setOriginX((float)transformedXC);
        paramsRev.setOriginY((float)transformedYC);
        paramsRev.setNumberOfPointsUsed(params.getNumberOfPointsUsed());

        // transform the new origin, then the new translation is what is needed for it to equal the params origin
        Transformer transformer = new Transformer();

        double[] revTransformedXYOrigin = transformer.applyTransformation(
            paramsRev, transformedXC, transformedYC);

        double revTransX = params.getOriginX() - revTransformedXYOrigin[0];
        double revTransY = params.getOriginY() - revTransformedXYOrigin[1];

        paramsRev.setTranslationX((float)revTransX);
        paramsRev.setTranslationY((float)revTransY);

        return paramsRev;
    }
    
    public boolean areSimilarByScaleAndRotation(TransformationParameters
        params1, TransformationParameters params2) {
        
        float scale = params1.getScale();
        float rotation = params1.getRotationInDegrees();

        float scale2 = params2.getScale();
        float rotation2 = params2.getRotationInDegrees();

        // scale similar within 10%?
        float scaleAvg = (scale + scale2) / 2.f;
        float scaleTol = 0.1f * scaleAvg;

        if (Math.abs(scale - scaleAvg) > scaleTol) {
            return false;
        }
        if (Math.abs(scale2 - scaleAvg) > scaleTol) {
            return false;
        }

        float diffR = AngleUtil.getAngleDifference(rotation, rotation2);
        if (Math.abs(diffR) > 20) {
            return false;
        }
        
        return true;
    }
    
    public Map<Integer, Set<Integer>> findSimilarByScaleAndRotation(
        List<TransformationParameters> params) {
        
        Map<Integer, Set<Integer>> sMap = new HashMap<Integer, Set<Integer>>();
        
        if (params.isEmpty()) {
            return sMap;
        }
        
        for (int i = 0; i < params.size(); ++i) {
            
            TransformationParameters param = params.get(i);
                       
            Set<Integer> sim = null;
            
            for (int j = (i + 1); j < params.size(); ++j) {
                
                TransformationParameters param2 = params.get(j);
                
                if (!areSimilarByScaleAndRotation(param, param2)) {
                    continue;
                }
                
                if ((param.getOriginX() != param2.getOriginX() )
                    || (param.getOriginY() != param2.getOriginY())) {

                    throw new IllegalArgumentException(
                    "params must all be with respect to same origin (x,y)");
                }
        
                if (sim == null) {
                    sim = new HashSet<Integer>();
                    sim.add(Integer.valueOf(i));
                }
                sim.add(Integer.valueOf(j));
            }
            
            if (sim != null) {
                sMap.put(Integer.valueOf(i), sim);
            }
        }
        
        // remove similar solutions
        for (int i = 0; i < params.size(); ++i) {
            
            Integer key = Integer.valueOf(i);
            Set<Integer> sim = sMap.get(key);
            if (sim == null) {
                continue;
            }
            
            for (int j = (i + 1); j < params.size(); ++j) {
                Integer key2 = Integer.valueOf(j);
                Set<Integer> sim2 = sMap.get(key2);
                if (sim2 == null) {
                    continue;
                }
                if (sim.equals(sim2)) {
                    sMap.remove(key2);
                }
            }
        }
        
        return sMap;
    }
    
    /**
     * average the paramsList without using outlier removal.  Note, the
     * standard deviations are not recalculated (later, if that's useful, they
     * should be added in quadrature here and set in the result).
     * NOTE that if the origins are not zero, this method cannot be used
     * @param paramsList
     * @return 
     */
    public TransformationParameters averageWithoutRemoval(
        List<TransformationParameters> paramsList) {
        
        boolean useRadians = true;
        List<Double> thetaCorr = new ArrayList<Double>(paramsList.size());
        for (TransformationParameters params : paramsList) {
            thetaCorr.add(Double.valueOf(params.getRotationInRadians()));
            if (params.getOriginX() != 0) {
                throw new IllegalArgumentException("originX must be 0");
            }
            if (params.getOriginY() != 0) {
                throw new IllegalArgumentException("originY must be 0");
            }
        }

        AngleUtil.correctForQuadrants(thetaCorr, useRadians);

        double thetaAvg = 0;
        float scaleAvg = 0;
        float transXAvg = 0;
        float transYAvg = 0;
        int nStats = 0;
        for (int i = 0; i < thetaCorr.size(); ++i) {
            thetaAvg += thetaCorr.get(i);
            TransformationParameters params = paramsList.get(i);
            scaleAvg += params.getScale();
            transXAvg += params.getTranslationX();
            transYAvg += params.getTranslationY();
            if (params.getNumberOfPointsUsed() >0 -1) {
                nStats ++;
            }
        }
        float n = paramsList.size();
        thetaAvg /= n;
        scaleAvg /= n;
        transXAvg /= n;
        transYAvg /= n;
        
        double thetaStDev = 0;
        double scaleStDev = 0;
        double transXStDev = 0;
        double transYStDev = 0;
        for (int i = 0; i < thetaCorr.size(); ++i) {
            thetaStDev += Math.pow((thetaCorr.get(i) - thetaAvg), 2);
            TransformationParameters params = paramsList.get(i);
            scaleStDev += Math.pow((params.getScale() - scaleAvg), 2);
            transXStDev += Math.pow((params.getTranslationX() - transXAvg), 2);
            transYStDev += Math.pow((params.getTranslationY() - transYAvg), 2);
        }
        thetaStDev = Math.sqrt(thetaStDev/(n - 1.));
        scaleStDev = Math.sqrt(scaleStDev/(n - 1.));
        transXStDev = Math.sqrt(transXStDev/(n - 1.));
        transYStDev = Math.sqrt(transYStDev/(n - 1.));
        
        if (thetaAvg > (2. * Math.PI)) {
            thetaAvg -= 2. * Math.PI;
        }
        
        TransformationParameters params = new TransformationParameters();
        params.setOriginX(0);
        params.setOriginY(0);
        params.setRotationInRadians((float)thetaAvg);
        params.setScale(scaleAvg);
        params.setTranslationX(transXAvg);
        params.setTranslationY(transYAvg);
        params.setStandardDeviations(new float[]{(float)scaleStDev, 
            (float)thetaStDev, (float)transXStDev, (float)transYStDev});
        params.setNumberOfPointsUsed(nStats);
        
        return params;
    }

    private List<Integer> findOutliersUsingMedianOfDifferences(
        List<? extends Number> values) {
        
        if (values.size() < 3) {
            return new ArrayList<Integer>();
        }
        
        // number of permutations = n!/(2!*(n-2)! = n*(n-1)/2
        int n = (values.size() * (values.size() - 1))/2;
        
        double[] diffs = new double[n];
        
        int count = 0;
        for (int i = 0; i < values.size(); ++i) {
            double v1 = values.get(i).doubleValue();
            for (int j = (i + 1); j < values.size(); ++j) {
                double v2 = values.get(j).doubleValue();
                diffs[count] = v1 - v2;
                count++;
            }
        }
        
        Arrays.sort(diffs);
        
        int idx = diffs.length/2;
        double median;
        if ((idx & 1) == 0) {
            median = (diffs[idx] + diffs[idx - 1])/2.f;
        } else {
            median = diffs[idx];
        }
        median = Math.abs(median);
        if (median == 0) {
            // widen to central 2 or 3 values
            if ((idx & 1) == 0) {
                median = (Math.abs(diffs[idx]) + Math.abs(diffs[idx - 1]))/2.f;
            } else {
                median = (Math.abs(diffs[idx]) + Math.abs(diffs[idx - 1]) +
                    Math.abs(diffs[idx + 1]))/3.;
            }
        }
        
        Set<Integer> keepIndexes = new HashSet<Integer>();
        for (int i = 0; i < values.size(); ++i) {
            double v1 = values.get(i).doubleValue();
            for (int j = (i + 1); j < values.size(); ++j) {
                double v2 = values.get(j).doubleValue();
                double diff = Math.abs(v1 - v2);
                if (diff <= median) {
                    keepIndexes.add(Integer.valueOf(i));
                    keepIndexes.add(Integer.valueOf(j));
                }
            }
        }
       
        List<Integer> remove = new ArrayList<Integer>();
        for (int i = 0; i < values.size(); ++i) {
            Integer key = Integer.valueOf(i);
            if (!keepIndexes.contains(key)) {
                remove.add(key);
            }
        }
        
        return remove;
    }

    private void useBipartiteMatching(PairIntArray matchedXY1, 
        Set<PairInt> pairIndexes) {
        
        float[][] invDist = new float[matchedXY1.getN()][matchedXY1.getN()];

        for (int i1 = 0; i1 < matchedXY1.getN(); ++i1) {
            int x1 = matchedXY1.getX(i1);
            int y1 = matchedXY1.getY(i1);
            invDist[i1] = new float[matchedXY1.getN()];
            for (int i2 = 0; i2 < matchedXY1.getN(); ++i2) {
                if (i1 == i2) {
                    invDist[i1][i2] = Float.MAX_VALUE;
                    continue;
                }
                int x2 = matchedXY1.getX(i2);
                int y2 = matchedXY1.getY(i2);
                int diffX = x1 - x2;
                int diffY = y1 - y2;
                double dist = Math.sqrt(diffX*diffX + diffY*diffY);
                invDist[i1][i2] = (float)(1./dist);
            }
        }

        //TODO: need a version of this method that uses matching < O(N^3)
    
        HungarianAlgorithm b = new HungarianAlgorithm();
        int[][] match = b.computeAssignments(invDist);
        
        for (int i = 0; i < match.length; i++) {

            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            
            if (idx2 < idx1) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }
            
            pairIndexes.add(new PairInt(idx1, idx2));
        }        
    }
    
    private void useGreedyMatching(PairIntArray matchedXY1, 
        Set<PairInt> pairIndexes) {
        
        if (matchedXY1.getN() == 2) {
            pairIndexes.add(new PairInt(0, 1));
            return;
        }

        for (int i1 = 0; i1 < matchedXY1.getN(); ++i1) {
            
            int x1 = matchedXY1.getX(i1);
            int y1 = matchedXY1.getY(i1);
            
            double maxDist = Double.MIN_VALUE;
            int maxDistIdx = -1;
            
            for (int i2 = 0; i2 < matchedXY1.getN(); ++i2) {
                
                if (i1 == i2) {
                    continue;
                }
                int x2 = matchedXY1.getX(i2);
                int y2 = matchedXY1.getY(i2);
                int diffX = x1 - x2;
                int diffY = y1 - y2;
                double dist = Math.sqrt(diffX*diffX + diffY*diffY);
                if (dist > maxDist) {
                    maxDist = dist;
                    maxDistIdx = i2;
                }
            }
            
            if (maxDistIdx > -1) {
                pairIndexes.add(new PairInt(i1, maxDistIdx));
            }
        }
    }
}
