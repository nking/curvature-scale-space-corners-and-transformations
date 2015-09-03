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
        
        calculateZeroCrossings(scaleSpaceCurve, curve);

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
        PairIntArray curve) {

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

        //TODO: adjust this to avoid noise
        for (int i = 1; i < scaleSpaceCurve.getK().length; ++i) {
            
            if (isZeroCrossing(scaleSpaceCurve.getK(i - 1),
                scaleSpaceCurve.getK(i))) {

                scaleSpaceCurve.addKIsZeroIdx(i, curve.getX(i), curve.getY(i));
            }
        }
        
//log.info("k=" + java.util.Arrays.toString(scaleSpaceCurve.getK()));

        scaleSpaceCurve.compressKIsZeroIdx();

    }

}
