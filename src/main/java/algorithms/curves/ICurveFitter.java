package algorithms.curves;

import java.io.IOException;

import algorithms.curves.GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM;

public interface ICurveFitter {

    public GEVYFit fitCurveKGreaterThanZero(WEIGHTS_DURING_CHISQSUM weightMethod) throws FailedToConvergeException, IOException;

    public void setDebug(boolean doUseDebug);
    
}
