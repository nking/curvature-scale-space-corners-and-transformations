package algorithms.curves;

import java.io.IOException;
import java.util.logging.Logger;

import algorithms.curves.GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM;
import algorithms.misc.MiscMath;
import algorithms.util.PolygonAndPointPlotter;

/**
 * a wrapper around NonQuadraticConjugateGradientSolver to use small ranges of parameters to
 * find best of local search results as the global result.
 * 
 * 
 * @author nichole
 *
 */
public class NQCGSWrapper implements ICurveFitter {
    
    protected NonQuadraticConjugateGradientSolver solver = null;
    
    public static final float kMinDefault = 0.00001f;
    public static final float kMaxDefault = 0.001f;
    public static final float kMinDefault2 = 0.001f;
    public static final float kMaxDefault2 = 0.1f;
    public static final float kMinDefault3 = 0.1f;
    public static final float kMaxDefault3 = 2.0f;
    public static final float sigmaMinDefault = 0.025f;
    public static final float sigmaMinDefault2 = 0.1f;
    public static final float sigmaMaxDefault = 0.5f;
    
    protected boolean debug = false;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public NQCGSWrapper(float[] xPoints, float[] yPoints,
        float[] xErrPoints, float[] yErrPoints) {

        solver = new NonQuadraticConjugateGradientSolver(xPoints, yPoints, xErrPoints, yErrPoints);
        
    } 
    
    public void setDebug(boolean dbg) {
        this.debug = dbg;
    }
    
    public GEVYFit fitCurveKGreaterThanZero(WEIGHTS_DURING_CHISQSUM weightMethod) throws FailedToConvergeException, IOException {

        GEVYFit bestFit = null;

        float yNorm = MiscMath.findMax(solver.y);

        int yMaxIndex = MiscMath.findYMaxIndex(solver.y);

        if (yMaxIndex == -1) {
            // all y's were zero in y
            return null;
        }
        float mu = solver.x[yMaxIndex];
        float yPeak = solver.y[yMaxIndex];

        int end = (yMaxIndex + 5);
        if (end > (solver.y.length - 1)) {
            end = solver.y.length - 1;
        }

        float nDmu = 2;
        for (int i = 0; i <= end; i++) {

            //float yNorm = y[i];
            float dmu;
            if (i < (solver.y.length - 1)) {
                dmu = (solver.x[i+1] - solver.x[i])/nDmu;
            } else {
                dmu = (solver.x[1] - solver.x[0])/nDmu;
            }
            
            for (int ii = 0; ii < nDmu; ii++) {
                
                mu = solver.x[i] + ii*dmu;

                //GEVYFit yfit = fitCurveKGreaterThanZero(weightMethod, mu, yNorm);
                GEVYFit yfit = fitCurveKGreaterThanZeroWithExtendedRanges(weightMethod, mu, dmu, yNorm);

                if (debug && (yfit != null)) {

                    String label = String.format("k=%.2e s=%.2e m=%.2e chisq=%.2f",
                        yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum());

                    plotFit(yfit, label);

                    log.fine(label);
                }

                if ( (bestFit == null) || (yfit.getChiSqSum() < bestFit.getChiSqSum()) ) {
                    bestFit = yfit;
                }
            }
        }

        return bestFit;
    }

    private GEVYFit fitCurveKGreaterThanZeroWithExtendedRanges(
        WEIGHTS_DURING_CHISQSUM weightMethod, float mu, float dmu, float yNorm) throws FailedToConvergeException {
        
        // start with ranges useful for EV Type I and II for the two-point correlation functions
        float kMin = kMinDefault;
        float kMax = kMaxDefault;
        float sigmaMin = sigmaMinDefault;
        float sigmaMax = sigmaMaxDefault;

        float chiSqStatisticLimit = 0f;//1000f;

        if (debug) {
            String str = String.format("kMin=%.7f kMax=%.7f sigmaMin=%.7f sigmaMax=%.7f", kMin, kMax, sigmaMin, sigmaMax);
            log.fine(str);
        }

        GEVYFit bestFit;

        bestFit = solver.fitCurveParametersSeparately(kMin, kMax, sigmaMin, sigmaMax, mu, mu + dmu);

        if (debug && (bestFit != null)) {

            String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f chistatistic=%.1f",
                bestFit.getK(), bestFit.getSigma(), bestFit.getMu(), bestFit.getChiSqSum(),
                bestFit.getChiSqStatistic());

            log.fine(str);
        }

        if ((bestFit == null) || (bestFit.getChiSqStatistic() > chiSqStatisticLimit)) {

            GEVYFit yfit;

            yfit = solver.fitCurveParametersSeparately(kMin, kMax, sigmaMin, sigmaMax, mu, mu + dmu);

            if (debug && (yfit != null)) {

                String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f chistatistic=%.1f",
                    yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum(),
                    yfit.getChiSqStatistic());

                log.fine(str);
            }

            if ((bestFit == null) || ((yfit != null) && (yfit.getChiSqStatistic() < bestFit.getChiSqStatistic()))) {
                bestFit = yfit;
            }

            if ((bestFit == null) || (bestFit.getChiSqStatistic() > chiSqStatisticLimit)) {

                kMin = kMinDefault2;
                kMax = kMaxDefault2;

                yfit = solver.fitCurveParametersSeparately(kMin, kMax, sigmaMin, sigmaMax, mu, mu + dmu);

                if (debug && (yfit != null)) {

                    String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f chistatistic=%.1f",
                    yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum(),
                    yfit.getChiSqStatistic());

                    log.fine(str);
                }

                if ((bestFit == null) || ((yfit != null) && (yfit.getChiSqStatistic() < bestFit.getChiSqStatistic()))) {
                    bestFit = yfit;
                }
            }

            if ((bestFit == null) || (bestFit.getChiSqStatistic() > chiSqStatisticLimit)) {

                kMin = kMinDefault3;
                kMax = kMaxDefault3;

                yfit = solver.fitCurveParametersSeparately(kMin, kMax, sigmaMin, sigmaMax, mu, mu + dmu);

                if (debug && (yfit != null)) {

                    String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f chistatistic=%.1f",
                    yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum(),
                    yfit.getChiSqStatistic());

                    log.fine(str);
                }

                if ((bestFit == null) || ((yfit != null) && (yfit.getChiSqStatistic() < bestFit.getChiSqStatistic()))) {
                    bestFit = yfit;
                }
            }

            if ((bestFit == null) || (bestFit.getChiSqStatistic() > chiSqStatisticLimit)) {

                kMin = kMinDefault3;
                kMax = kMaxDefault3;// might need kMinDefault3*10

                sigmaMin = sigmaMinDefault2;

                yfit = solver.fitCurveParametersSeparately(kMin, kMax, sigmaMin, sigmaMax, mu, mu + dmu);

                if (debug && (yfit != null)) {

                    String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f chistatistic=%.1f",
                    yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum(),
                    yfit.getChiSqStatistic());

                    log.fine(str);
                }

                if ((bestFit == null) || ((yfit != null) && (yfit.getChiSqStatistic() < bestFit.getChiSqStatistic()))) {
                    bestFit = yfit;
                }
            }

            if ((bestFit == null) || (bestFit.getChiSqStatistic() > chiSqStatisticLimit)) {

                kMin = kMinDefault3;
                kMax = 1.0f;

                sigmaMin = sigmaMinDefault;
                sigmaMax = 0.1f;

                yfit = solver.fitCurveParametersSeparately(kMin, kMax, sigmaMin, sigmaMax, mu, mu + dmu);

                if (debug && (yfit != null)) {

                    String str = String.format("fit: k=%.7f s=%.7f m=%.7f chisq=%.1f chistatistic=%.1f",
                    yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum(),
                    yfit.getChiSqStatistic());

                    log.fine(str);
                }

                if ((bestFit == null) || ((yfit != null) && (yfit.getChiSqStatistic() < bestFit.getChiSqStatistic()))) {
                    bestFit = yfit;
                }
            }
        }

        return bestFit;
    }

    protected void plotFit(GEVYFit yfit, String label) throws IOException {

        float xIntervalHalf = (yfit.getX()[1] - yfit.getX()[0]) / 2.0f;
        float xmin = yfit.getX()[0] - xIntervalHalf;
        float xmax = yfit.getX()[ yfit.getX().length - 1];
        float ymin = 0.0f;
        float ymax = MiscMath.findMax(solver.y);

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(xmin, xmax, ymin, ymax);

        try {
            plotter.addPlot(solver.x, solver.y, solver.xe, solver.ye, yfit.getX(), yfit.getYFit(), label);
            plotter.writeFile2();
        } catch (Exception e) {
            Logger.getLogger(this.getClass().getSimpleName()).severe(e.getMessage());
        }
    }

    public GEVYFit fitCurveKGreaterThanZero(WEIGHTS_DURING_CHISQSUM errors,
        float kMin, float kMax, float sigmaMin, float sigmaMax) throws FailedToConvergeException, IOException {
        
        GEVYFit bestFit = null;

        float yNorm = MiscMath.findMax(solver.y);

        int yMaxIndex = MiscMath.findYMaxIndex(solver.y);

        if (yMaxIndex == -1) {
            // all y's were zero in y
            return null;
        }
        float mu = solver.x[yMaxIndex];
        float yPeak = solver.y[yMaxIndex];

        int end = (yMaxIndex + 5);
        if (end > (solver.y.length - 1)) {
            end = solver.y.length - 1;
        }

        float nDmu = 2;
        for (int i = 0; i <= end; i++) {

            //float yNorm = y[i];
            float dmu;
            if (i < (solver.y.length - 1)) {
                dmu = (solver.x[i+1] - solver.x[i])/nDmu;
            } else {
                dmu = (solver.x[1] - solver.x[0])/nDmu;
            }
            
            for (int ii = 0; ii < nDmu; ii++) {
                
                mu = solver.x[i] + ii*dmu;

                //GEVYFit yfit = fitCurveKGreaterThanZero(weightMethod, mu, yNorm);
                GEVYFit yfit = solver.fitCurveParametersSeparately(kMin, kMax, sigmaMin, sigmaMax, mu, mu + dmu);
                
                if (debug && (yfit != null)) {

                    String label = String.format("k=%.2e s=%.2e m=%.2e chisq=%.2f",
                        yfit.getK(), yfit.getSigma(), yfit.getMu(), yfit.getChiSqSum());

                    plotFit(yfit, label);

                    log.fine(label);
                }

                if ( (bestFit == null) || (yfit.getChiSqSum() < bestFit.getChiSqSum()) ) {
                    bestFit = yfit;
                }
            }
        }

        return bestFit;
    }

}
