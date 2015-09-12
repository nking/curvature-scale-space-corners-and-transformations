package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class ScaleSpaceCurvature {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * compute scale space metrics of curve, given sigma
     * @param curve the x,y pairs of points for which to calculate the scale 
     * space metrics
     * @param kernelSigma the sigma of the Gaussian kernel to convolve the curve
     * with.
     * @param resultingSigma the equivalent sigma of the resulting curve.  This
     * is useful if feeding back a curve that is already convolved by a sigma=4
     * kernel.  The convolution of that with a sigma=2 kernel results in a
     * sigma=4*sqrt(2) kernel.
     * 
     * @return  
     */
    public ScaleSpaceCurve computeCurvature(PairIntArray curve, 
        SIGMA kernelSigma, float resultingSigma) {
        
        float[] gFirstDeriv = Gaussian1DFirstDeriv.getKernel(kernelSigma);

        float[] gSecondDeriv = Gaussian1DSecondDeriv.getKernel(kernelSigma);
 
        return computeCurvature(curve, resultingSigma, gFirstDeriv, 
            gSecondDeriv);
    }
    
    /**
     * compute scale space metrics of curve, given sigma
     * @param curve the x,y pairs of points for which to calculate the scale 
     * space metrics
     * @param kernelSigma the sigma of the Gaussian kernel to convolve the curve
     * with.
     * @param resultingSigma the equivalent sigma of the resulting curve.  This
     * is useful if feeding back a curve that is already convolved by a sigma=4
     * kernel.  The convolution of that with a sigma=2 kernel results in a
     * sigma=4*sqrt(2) kernel.
     * 
     * @return  
     */
    public ScaleSpaceCurve computeCurvature(PairIntArray curve, 
        float kernelSigma, float resultingSigma) {
        
        float[] gFirstDeriv = Gaussian1DFirstDeriv.getKernel(kernelSigma);

        float[] gSecondDeriv = Gaussian1DSecondDeriv.getKernel(kernelSigma);
 
        return computeCurvature(curve, resultingSigma, gFirstDeriv, 
            gSecondDeriv);
    }
    
    /**
     * compute scale space metrics of curve, given sigma
     * @param curve the x,y pairs of points for which to calculate the scale 
     * space metrics
     * @param resultingSigma the equivalent sigma of the resulting curve.  This
     * is useful if feeding back a curve that is already convolved by a sigma=4
     * kernel.  The convolution of that with a sigma=2 kernel results in a
     * sigma=4*sqrt(2) kernel.
     * 
     * @return  
     */
    public ScaleSpaceCurve computeCurvature(PairIntArray curve, 
        float resultingSigma, float[] gFirstDeriv, float[] gSecondDeriv) {

        int n = curve.getN();

        boolean isClosedCurved = (curve instanceof PairIntArrayWithColor)
            && (((PairIntArrayWithColor) curve).getColor() == 1);

        ScaleSpaceCurve scaleSpaceCurve = new ScaleSpaceCurve(resultingSigma,
            curve, isClosedCurved);

        /*
                  X_dot(t,o~) * Y_dot_dot(t,o~) - Y_dot(t,o~) * X_dot_dot(t,o~) 
        k(t,o~) = -------------------------------------------------------------
                                     (X^2(t,o~) + Y^2(t,o~))^1.5
        */
 
        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();

        for (int i = 0; i < n; i++) {

            double xFirstDerivInteg = kernel1DHelper.convolvePointWithKernel(
                curve, i, gFirstDeriv, true);
            
            double yFirstDerivInteg = kernel1DHelper.convolvePointWithKernel(
                curve, i, gFirstDeriv, false);
            
            double xSecondDerivInteg = kernel1DHelper.convolvePointWithKernel(
                curve, i, gSecondDeriv, true);
            
            double ySecondDerivInteg = kernel1DHelper.convolvePointWithKernel(
                curve, i, gSecondDeriv, false);
    
            double denominator = Math.pow(
                ((xFirstDerivInteg * xFirstDerivInteg) +
                (yFirstDerivInteg * yFirstDerivInteg)), 1.5);
            
            double numerator = ((xFirstDerivInteg * ySecondDerivInteg) -
                (yFirstDerivInteg * xSecondDerivInteg));
            
            double curvature = (denominator == 0)  ? 
                (numerator == 0) ? 0 : Double.POSITIVE_INFINITY
                : numerator / denominator;
           
            scaleSpaceCurve.setK(i, (float)curvature);
            
            /*
            if using ScaleSpaceCurve2 to capture the 2nd derivatives, set
            them here.
            */
        }
        
        calculateZeroCrossings(scaleSpaceCurve, curve, isClosedCurved);

        return scaleSpaceCurve;
    }
     
    private boolean isZeroCrossing(final float kPrev, final float k) {
        if (k <= 0) {
            if (kPrev >= 0) {
                return true;
            }
        }
        if (k >= 0) {
            if (kPrev <= 0) {
                return true;
            }
        }
        return false;
    }

    private void calculateZeroCrossings(ScaleSpaceCurve scaleSpaceCurve,
        PairIntArray curve, boolean isClosedCurved) {

        final int n = scaleSpaceCurve.getSize();
        
        if (n < 3) {
            scaleSpaceCurve.compressKIsZeroIdx();
            return;
        }
        
        /*
        // find the points where the curvature changes from + to -, but
        // discard such changes when they are noise.
        float[] x = new float[scaleSpaceCurve.getK().length];
        float[] y = new float[x.length];
        for (int i = 0; i < x.length; i++) {
            x[i] = i;
            y[i] = scaleSpaceCurve.getK(i);
        }
        float yMin = 1.2f*MiscMath.findMin(y);
        float yMax = 1.2f*MiscMath.findMax(y);
        try {
            String time = Long.toString((long)Math.floor(
                System.currentTimeMillis()));
            time = time.substring(4, time.length());
            int fileNumber = Long.valueOf(time).intValue();
            ScatterPointPlotterPNG plotter = new ScatterPointPlotterPNG();
            plotter.plot(0.0f, x[x.length - 1], yMin, yMax,
                x, y, "curvature", "t", "k");
            plotter.writeFile(fileNumber);

            StringBuilder sb = new StringBuilder(fileNumber + " coords:");
            for (int i = 0; i < curve.getN(); i++) {
                sb.append(String.format("t=%d (%d,%d)\n", i, curve.getX(i), curve.getY(i)));
            }
            log.info(sb.toString());
            
        } catch (IOException e) {
            log.severe(e.getMessage());
        }
        */
         
        // simple zero crossings works well on curves, but not on straight lines
        //log.info("new curve (" + scaleSpaceCurve.getSize() + ")");
        float maxInflectionK = Integer.MIN_VALUE;
        if (isClosedCurved) {
            if (isZeroCrossing(scaleSpaceCurve.getK(n - 1),
                scaleSpaceCurve.getK(0))) {
                float k = scaleSpaceCurve.getK(0);
                if (k < 0) {
                    k *= -1;
                }
                if (k > maxInflectionK) {
                    maxInflectionK = k;
                }
                //log.info("   i=" + 0 + " x=" + curve.getX(0) + " y=" + curve.getY(0) + " k=" + k);
                //scaleSpaceCurve.addKIsZeroIdx(0, curve.getX(0), curve.getY(0));
            }
        }

        for (int i = 1; i < scaleSpaceCurve.getK().length; ++i) {
            if (isZeroCrossing(scaleSpaceCurve.getK(i - 1),
                scaleSpaceCurve.getK(i))) {
                float k = scaleSpaceCurve.getK(i);
                if (k < 0) {
                    k *= -1;
                }
                if (k > maxInflectionK) {
                    maxInflectionK = k;
                }
                //log.info("   i=" + i + " x=" + curve.getX(i) + " y=" + curve.getY(i) + " k=" + k);
                //scaleSpaceCurve.addKIsZeroIdx(i, curve.getX(i), curve.getY(i));
            }
        }
        
        //adjusting for straight lines.  consecutive points ~0 should
        //be summarized to just the mid point as the inflection point
        
        int endingIdx = -1;
        int beginningIdx = -1;
        boolean addToEnd = false;
        int addIdx = -1;
        if (isClosedCurved) {
            // scan the wrap around of end to beginning
            if ((Math.abs(scaleSpaceCurve.getK()[n - 1]) <= maxInflectionK) && 
                (Math.abs(scaleSpaceCurve.getK()[0]) <= maxInflectionK)) {
                for (int i = (n - 1); i > -1; --i) {
                    if (Math.abs(scaleSpaceCurve.getK()[i]) <= maxInflectionK) {
                        endingIdx = i;
                    } else {
                        break;
                    }
                }
                for (int i = 0; i < n; ++i) {
                    if (Math.abs(scaleSpaceCurve.getK()[i]) <= maxInflectionK) {
                        beginningIdx = i;
                    } else {
                        break;
                    }
                }
                int midIdx = endingIdx + ((n - endingIdx + beginningIdx)/2);
                if (midIdx > (n - 1)) {
                    midIdx = midIdx - n;
                } else {
                    addToEnd = true;
                    addIdx = midIdx;
                }
                //log.info("   beginningIdx=" + beginningIdx + " endIdx=" + endingIdx + " midIdx=" + midIdx);
                //log.info("     " + beginningIdx + " (" + curve.getX(beginningIdx) + "," + curve.getY(beginningIdx) + ")"
                //    + "   " + endingIdx + " (" + curve.getX(endingIdx) + "," + curve.getY(endingIdx) + ")"
                //);
                if (!addToEnd) { 
                    scaleSpaceCurve.addKIsZeroIdx(midIdx, curve.getX(midIdx), curve.getY(midIdx));
                }
            }
        }
        
        int i0 = (beginningIdx > -1) ? (beginningIdx + 1) : 0;
        int i1 = (endingIdx > -1) ? endingIdx : n;
        int startIdx = -1;
        int endIdx = -1;
        for (int i = i0; i < i1; ++i) {
            float k = scaleSpaceCurve.getK()[i];
            if (Math.abs(k) <= maxInflectionK) {
                if (startIdx == -1) {
                    startIdx = i;
                    endIdx = i;
                } else {
                    endIdx = i;
                }
            } else {
                if (startIdx > -1) {
                    int midIdx = (startIdx + endIdx)/2;
                    //log.info("   startIdx=" + startIdx + " endIdx=" + (i-1) + " midIdx=" + midIdx);
                    //log.info("     " + startIdx + " (" + curve.getX(startIdx) + "," + curve.getY(startIdx) + ")"
                    //    + "   " + (i-1) + " (" + curve.getX(i-1) + "," + curve.getY(i-1) + ")"
                    //);
                    scaleSpaceCurve.addKIsZeroIdx(midIdx, curve.getX(midIdx), curve.getY(midIdx));
                }
                startIdx = -1;
                endIdx = -1;
            }
        }
        if (startIdx > -1) {
             int midIdx = (startIdx + endIdx)/2;
             //log.info("   startIdx=" + startIdx + " endIdx=" + endIdx + " midIdx=" + midIdx);
             //   log.info("     " + startIdx + " (" + curve.getX(startIdx) + "," + curve.getY(startIdx) + ")"
             //   + "   " + (endIdx) + " (" + curve.getX(endIdx) + "," + curve.getY(endIdx) + ")"
             //);
             scaleSpaceCurve.addKIsZeroIdx(midIdx, curve.getX(midIdx), curve.getY(midIdx));
        }
        if (addToEnd) { 
            scaleSpaceCurve.addKIsZeroIdx(addIdx, curve.getX(addIdx), curve.getY(addIdx));
        }
       
//log.info("k=" + java.util.Arrays.toString(scaleSpaceCurve.getK()));
//log.info("k=" + java.util.Arrays.toString(scaleSpaceCurve.getKIsZeroIdx()));

        scaleSpaceCurve.compressKIsZeroIdx();

    }

}
