package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;

/**
 *
 * @author nichole
 */
public class ScaleSpaceCurvature {
    
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
                     
            // set the zeros in curvature
            if (i > 0) {
                if (isZeroCrossing(scaleSpaceCurve.getK(i-1), 
                    scaleSpaceCurve.getK(i))) {
                    
                    float p = (float)i/(float)curve.getN();
                    
                    //System.out.println("sigma=" + resultingSigma 
                    //    + " t=" + p);
                    scaleSpaceCurve.addKIsZeroIdx(i, curve.getX(i), curve.getY(i));
                }
            }
        }

        scaleSpaceCurve.compressKIsZeroIdx();

        return scaleSpaceCurve;
    }
     
    private boolean isZeroCrossing(float kPrev, float k) {
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
    
}
