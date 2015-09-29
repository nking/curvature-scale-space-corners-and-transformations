package algorithms.curves;

import java.io.IOException;

import algorithms.curves.GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM;

/**
 *
 * @author nichole
 */
public interface ICurveFitter {

    /**
     *
     * @param weightMethod
     * @return
     * @throws FailedToConvergeException
     * @throws IOException
     */
    public GEVYFit fitCurveKGreaterThanZero(WEIGHTS_DURING_CHISQSUM weightMethod) throws FailedToConvergeException, IOException;

    /**
     *
     * @param doUseDebug
     */
    public void setDebug(boolean doUseDebug);
    
}
