package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.PolygonPlotterPNG;
import java.io.IOException;

/**
 *
 * @author nichole
 */
public abstract class AbstractPointMatcher {

    protected boolean costIsNumAndDiff = false;
    
    public void setCostToNumMatchedAndDiffFromModel() {
        costIsNumAndDiff = true;
    }
    
    /**
     * sort the fits by descending number of matches. uses a slightly more 
     * iterative pattern adapted from "Programming Pearls" by Jon Bentley.
     * The average runtime is O(N*lgN) and the worse case runtime is O(N^2) 
     * but it is 3 times faster than sortByDescendingMatches0(...).
     * 
     * @param fits
     * @param idxLo
     * @param idxHi, upper index, inclusive
     */
    void sortByDescendingMatches(TransformationPointFit[] a, int idxLo,
        int idxHi) {

        if (idxLo < idxHi) {

            TransformationPointFit x = a[idxLo];
            int store = idxLo;
            int idxMid = idxHi + 1;
            
            //replace the recursive partition method frame load with iteration:
            while (true) {
                do {
                    store++;
                } while ((store <= idxHi) && (fitIsBetter(x, a[store])));
                
                do {
                    idxMid--;
                } while (fitIsBetter(a[idxMid], x));
                
                if (store > idxMid) {
                    break;
                }
                TransformationPointFit swap = a[store];
                a[store] = a[idxMid];
                a[idxMid] = swap;
            }
            TransformationPointFit swap = a[idxLo];
            a[idxLo] = a[idxMid];
            a[idxMid] = swap;

            sortByDescendingMatches(a, idxLo, idxMid - 1);

            sortByDescendingMatches(a, idxMid + 1, idxHi);
        }
    }


     /**
     * sort the fits by descending number of matches.
     * @param fits
     * @param idxLo
     * @param idxHi, upper index, inclusive
     */
    void sortByDescendingMatches0(TransformationPointFit[] fits, int idxLo,
        int idxHi) {

        if (idxLo < idxHi) {

            int idxMid = partition(fits, idxLo, idxHi);

            sortByDescendingMatches0(fits, idxLo, idxMid - 1);

            sortByDescendingMatches0(fits, idxMid + 1, idxHi);
        }
    }

    private int partition(TransformationPointFit[] fits, int idxLo, int idxHi) {

        TransformationPointFit x = fits[idxHi];

        int store = idxLo - 1;

        for (int i = idxLo; i < idxHi; ++i) {
            if (fitIsBetter(x, fits[i])) {
                store++;
                TransformationPointFit swap = fits[store];
                fits[store] = fits[i];
                fits[i] = swap;
            }
        }

        store++;
        TransformationPointFit swap = fits[store];
        fits[store] = fits[idxHi];
        fits[idxHi] = swap;

        return store;
    }

    protected boolean fitIsSimilar(TransformationPointFit fit1, 
        TransformationPointFit fit2, float scaleTolerance, 
        float rotInDegreesTolerance, float translationTolerance) {
        
        if (fit1 == null || fit2 == null) {
            return false;
        }
        TransformationParameters params1 = fit1.getParameters();
        TransformationParameters params2 = fit2.getParameters();
        float rotA = params1.getRotationInDegrees();
        float rotB = params2.getRotationInDegrees();
        if (rotA > rotB) {
            float swap = rotA;
            rotA = rotB;
            rotB = swap;
        }
        if (((rotB - rotA) > rotInDegreesTolerance) && (((rotA + 360) - rotB) > rotInDegreesTolerance)) {
            return false;
        }
        if ((Math.abs(params1.getScale() - params2.getScale()) <= scaleTolerance) && (Math.abs(params1.getTranslationX() - params2.getTranslationX()) <= translationTolerance) && (Math.abs(params1.getTranslationY() - params2.getTranslationY()) <= translationTolerance)) {
            return true;
        }
        return false;
    }

    public boolean fitIsBetter(TransformationPointFit bestFit, TransformationPointFit compareFit) {
        if (costIsNumAndDiff) {
            return fitIsBetterUseNumAndDiff(bestFit, compareFit);
        }
        if (compareFit == null) {
            return false;
        }
        if (bestFit == null) {
            return true;
        }
        int comp = bestFit.compareTo(compareFit);
        if (comp == 1) {
            return true;
        }
        return false;
    }
    
    public boolean fitIsBetterNormalized(TransformationPointFit bestFit, 
        TransformationPointFit compareFit) {
        
        if (compareFit == null) {
            return false;
        }
        if (bestFit == null) {
            return true;
        }
        
        boolean ans = false;
        int comp = bestFit.compareToUsingNormalizedMatches(compareFit);
        if (comp == 1) {
            ans = true;
        }
        
        return ans;
    }

    /**
     * compare bestFit to compareFit and return
     * -1 if bestFit is better
     * 0 if they are equal
     * 1 if compareFit is better
     * @param bestFit
     * @param compareFit
     * @return
     */
    public int compare(TransformationPointFit bestFit, TransformationPointFit compareFit) {
        /*
        if (costIsNumAndDiff) {
        return fitIsBetterUseNumAndDiff(bestFit, compareFit);
        }
         */
        if (compareFit == null && bestFit == null) {
            return 0;
        } else if (compareFit == null) {
            return -1;
        } else if (bestFit == null) {
            return 1;
        }
        int comp = bestFit.compareTo(compareFit);
        return comp;
    }

    public boolean fitIsBetterUseNumAndDiff(TransformationPointFit bestFit, TransformationPointFit compareFit) {
        if (compareFit == null) {
            return false;
        }
        if (bestFit == null) {
            if (compareFit.getNumberOfMatchedPoints() > 0) {
                return true;
            } else {
                return false;
            }
        }
        float compN = compareFit.getNumberOfMatchedPoints();
        float bestN = bestFit.getNumberOfMatchedPoints();
        double compAvg = compareFit.getMeanDistFromModel();
        double compS = compareFit.getStDevFromMean();
        double compAvgS = compAvg + compS;
        double bestAvg = bestFit.getMeanDistFromModel();
        double bestS = bestFit.getStDevFromMean();
        double bestAvgS = bestAvg + bestS;
        float f = 1.5f;
        if (!Double.isNaN(compAvg)) {
            if ((compN / bestN) >= f) {
                return true;
            } else if ((compN >= bestN) && (compAvg < bestAvg) && (compS < bestS)) {
                return true;
            }
        }
        return false;
    }

    public boolean fitIsBetter(ProjectiveFit bestFit, ProjectiveFit compareFit) {
        if (compareFit == null) {
            return false;
        }
        if (bestFit == null) {
            return true;
        }
        int nMatches = compareFit.getNumberOfPoints();
        if (nMatches > bestFit.getNumberOfPoints()) {
            return true;
        } else if (nMatches == bestFit.getNumberOfPoints()) {
            if (!Double.isNaN(compareFit.getMeanDistFromModel()) && (compareFit.getMeanDistFromModel() < bestFit.getMeanDistFromModel())) {
                return true;
            } else if (compareFit.getMeanDistFromModel() == bestFit.getMeanDistFromModel()) {
                if (compareFit.getStdDevOfMean() < bestFit.getStdDevOfMean()) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * a fitness function that tries to allow a smaller number of points roughly
     * fit to be seen in contrast to a larger number of points that
     * are not a better fit, but have better stats due to matching alot of
     * scattered points.
     * if numberOfMatched/maxMatchable is not infinite:
     * if the mean/10 of the comparison fit is better than the best mean/10,
     * a true is returned, else
     * compares the standard
     * deviation from the mean difference to the model fit and returns true
     * if compareFit has a smaller value.
     * @param bestFit
     * @param compareFit
     * @return
     */
    protected boolean fitIsBetter2(TransformationPointFit bestFit, TransformationPointFit compareFit) {
        if (compareFit == null) {
            return false;
        }
        if (bestFit == null) {
            return true;
        }
        float compareNStat = (float) compareFit.getNumberOfMatchedPoints() / (float) compareFit.getNMaxMatchable();
        if (Float.isInfinite(compareNStat)) {
            return false;
        }
        float bestNStat = (float) bestFit.getNumberOfMatchedPoints() / (float) bestFit.getNMaxMatchable();
        if (Float.isInfinite(bestNStat)) {
            return true;
        }
        if ((bestFit.getNumberOfMatchedPoints() == 0) && (compareFit.getNumberOfMatchedPoints() > 0)) {
            return true;
        } else if (compareFit.getNumberOfMatchedPoints() == 0) {
            return false;
        }
        double bestMean = bestFit.getMeanDistFromModel();
        double compareMean = compareFit.getMeanDistFromModel();
        int comp = Double.compare(compareMean, bestMean);
        if (comp < 0) {
            return true;
        } else if (comp > 0) {
            return false;
        }
        double bestStdDevMean = bestFit.getStDevFromMean();
        double compareStdDevMean = compareFit.getStDevFromMean();
        // a smaller std dev from mean is a better fit
        if (Double.compare(compareStdDevMean, bestStdDevMean) < 0) {
            return true;
        } else if (compareStdDevMean == bestStdDevMean) {
            return Double.compare(compareMean, bestMean) < 0;
        }
        return false;
    }

    protected boolean areEqual(TransformationParameters[] lastParams, TransformationParameters[] currentParams) {
        if (lastParams.length != currentParams.length) {
            throw new IllegalArgumentException("lastParams.length must be equal to currentParams.length");
        }
        for (int i = 0; i < lastParams.length; ++i) {
            TransformationParameters p0 = lastParams[i];
            TransformationParameters p1 = currentParams[i];
            if ((p0 == null) && (p1 != null)) {
                return false;
            } else if ((p0 != null) && (p1 == null)) {
                return false;
            } else if (p0 == null && p1 == null) {
                continue;
            } else if (!p0.equals(p1)) {
                return false;
            }
        }
        return true;
    }

    /**
     * compare the fields mean distance from model and standard deviation
     * from mean to find if they are similar within a tolerance, then
     * find if the parameters are the same and return
     * <pre>
     * 0==same fit;  1==similar fits;  -1==different fits
     * </pre>
     * @param bestFit
     * @param fit
     * @return comparisonResult 0==same fit;  1==similar fits;  -1==different fits
     */
    protected int fitsAreSimilarWithDiffParameters(TransformationPointFit bestFit, TransformationPointFit fit) {
        if (bestFit == null || fit == null) {
            return -1;
        }
        int comp = bestFit.isSimilarWithDiffParameters(fit);
        return comp;
    }

    /** compare the bestFit and fit tolerances and return
    <pre>
    -1 : both are not null and bestFit tolerances are smaller
    0 : both are not null and tolerances are same.
    1 : both are not null and fit tolerances are smaller
    2 : both are not null and the x and y fits and smaller and larger in a mix
    3 : either bestFit or fit is null
    </pre>
     */
    protected int compareTolerance(TransformationPointFit bestFit, TransformationPointFit fit) {
        if (bestFit == null || fit == null) {
            return 3;
        }
        float diffTolX = bestFit.getTranslationXTolerance() - fit.getTranslationXTolerance();
        float diffTolY = bestFit.getTranslationYTolerance() - fit.getTranslationYTolerance();
        if ((Math.abs(diffTolX) < 1) && (Math.abs(diffTolY) < 1)) {
            return 0;
        } else if ((diffTolX > 0) && (diffTolY > 0)) {
            return 1;
        } else if ((diffTolX < 0) && (diffTolY < 0)) {
            return -1;
        } else {
            return 2;
        }
    }
    
    /**
     * estimate whether fit has converged by the mean distance from model and 
     * the standard deviation of a point from that mean and by the
     * number of points matched.  note that a false return may still be a 
     * converged solution for some cases.
     * @param fit
     * @param nMaxMatchable
     * @return 
     */
    public boolean hasConverged(TransformationPointFit fit, int nMaxMatchable) {
    
        if (fit == null) {
            return false;
        }
        
        int bestNMatches = fit.getNumberOfMatchedPoints();

        double bestAvg = fit.getMeanDistFromModel();

        double bestS = fit.getStDevFromMean();

        float fracMatched = (float)bestNMatches/(float)nMaxMatchable;

        if ((bestAvg < 1) && (bestS < 1)) {
            if (fracMatched > 0.9) {
                return true;
            }
        } else if ((bestAvg < 0.5) && (bestS < 0.5)) {
            if (nMaxMatchable > 10 && bestNMatches > 10) {
                return true;
            }
        }

        return false;
    }
    
    private PolygonPlotterPNG plotter = null;
    //TODO: put this in an aspect
    private void plotTranslationSimplex(TransformationPointFit[] fits,
        float minX, float maxX, float minY, float maxY, java.awt.Color clr) {

        try {
            if (plotter == null) {
                plotter = new PolygonPlotterPNG(minX, maxX, minY, maxY,
                    "translation simplex", "transX", "transY");
            }
            int count = 0;
            for (TransformationPointFit fit : fits) {
                if (fit == null) {
                    continue;
                }
                count++;
            }
            double[] x = new double[count];
            double[] y = new double[x.length];
            count = 0;
            for (TransformationPointFit fit : fits) {
                if (fit == null) {
                    continue;
                }
                x[count] = fit.getTranslationX();
                y[count] = fit.getTranslationY();
                count++;
            }

            if (clr == null) {
                plotter.addPolygon(x, y);
            } else {
                plotter.addPolygon(x, y, clr);
            }

        } catch (IOException e) {
            System.err.println(e.getMessage());
        }
    }

    private void writeTranslationSimplexPlot() {

        if (plotter == null) {
            return;
        }

        try {
            plotter.writeFile(MiscDebug.getCurrentTimeFormatted());
        } catch (IOException e) {
            System.err.println(e.getMessage());
        }

    }
}
