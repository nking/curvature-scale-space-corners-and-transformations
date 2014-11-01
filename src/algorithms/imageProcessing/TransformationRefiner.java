package algorithms.imageProcessing;

import java.util.logging.Logger;

/**
 A minimization routine to refine the given transformation parameters.
 
 Internally, it's using the Nelder-Mead Downhill Simplex method,
but that could be improved with a non-linear conjugate gradient solver
(unless skew is included... then the first derivative used in the 
conjugate gradient solver would not proceed in one direction locally).
 
 @author nichole
 */
public class TransformationRefiner {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    public TransformationRefiner() {
        
    }
    
    private SearchableCurve[] createSearchableCurves(PairIntArray[] edges) {
        
        if (edges == null) {
            return null;
        }
        if (edges.length == 0) {
            return new SearchableCurve[0];
        }
        
        SearchableCurve[] curves = new SearchableCurve[edges.length];
        for (int i = 0; i < edges.length; i++) {
            curves[i] = new SearchableCurve(edges[i]);
        }
        
        return curves;
    }
    
    public TransformationParameters refineTransformation(PairIntArray[] edges1, 
        PairIntArray[] edges2, TransformationParameters params) {
        
        SearchableCurve[] searchableEdges2 = createSearchableCurves(edges2);
            
        //TODO: set this empirically from tests
        double convergence = 0;
                
        double r = params.getRotationInRadians();
        double s = params.getScale();
        
        double dR = (1.0 * Math.PI/180.);
        double dS = 0.1;
        
        double rMin = r - (10 * Math.PI/180);
        double rMax = r + (10 * Math.PI/180);
        double sMin = s - 0.5;
        double sMax = s + 0.5;
        
        //TODO: these starting points could be improved.
        
        // start with simplex for at least 3 points (fitting 2 parameters)
        TransformationFit[] fits = new TransformationFit[9];
        fits[0] = transformEdges(r, s, edges1, searchableEdges2);
        fits[1] = transformEdges(r + dR, s, edges1, searchableEdges2);
        fits[2] = transformEdges(r - dR, s, edges1, searchableEdges2);
        // adding additional points.  need small changes more than changes
        // near the size of the value
        fits[3] = transformEdges(r, s + dS, edges1, searchableEdges2);
        fits[4] = transformEdges(r, s - dS, edges1, searchableEdges2);
        fits[5] = transformEdges(r + dR, s + dS, edges1, searchableEdges2);
        fits[6] = transformEdges(r + dR, s - dS, edges1, searchableEdges2);
        fits[7] = transformEdges(r - dR, s + dS, edges1, searchableEdges2);
        fits[8] = transformEdges(r - dR, s - dS, edges1, searchableEdges2);

        float alpha = 1;   // > 0
        float gamma = 2;   // > 1
        float beta = -0.5f; 
        float tau = 0.5f;
            
        boolean go = true;

        int nMaxIter = 100;
        int nIter = 0;
        
        int bestFitIdx = 0;
        int secondWorstFitIdx = fits.length - 2;
        int worstFitIdx = fits.length - 1;

        while (go && (nIter < nMaxIter)) {

            sortFromMinToMax(fits, 0, (fits.length - 1));

            // determine center for all points excepting the worse fit
            double rSum = 0.0;
            double sSum = 0.0;
            for (int i = 0; i < (fits.length - 1); i++) {
                rSum += fits[i].getRotationInRadians();
                sSum += fits[i].getScale();
            }
            r = rSum / (double)(fits.length - 1);
            s = sSum / (double)(fits.length - 1);

            // "Reflection"
            double rReflect = r + (alpha * 
                (r - fits[worstFitIdx].getRotationInRadians()));
            double sReflect = s +
                (s - alpha * (fits[worstFitIdx].getScale()));
            TransformationFit fitReflected = transformEdges(rReflect, sReflect, 
                edges1, searchableEdges2);
            
            boolean relectIsWithinBounds = 
                (rReflect >= rMin) && (rReflect <= rMax) 
                && (sReflect >= sMin) && (sReflect <= sMax);

            if ((fitReflected.getChiSqSum() < fits[secondWorstFitIdx].getChiSqSum())
                && relectIsWithinBounds &&
                (fitReflected.getChiSqSum() > fits[bestFitIdx].getChiSqSum())
                ) {
                
                fits[worstFitIdx] = fitReflected;
                
            } else {
                
                if ((fitReflected.getChiSqSum() < fits[bestFitIdx].getChiSqSum())
                && relectIsWithinBounds) {
                    
                    // "Expansion"
                    double rExpansion = r + (gamma * 
                        (r - fits[worstFitIdx].getRotationInRadians()));
                    double sExpansion = s + (gamma * 
                        (s - fits[worstFitIdx].getScale()));
                    TransformationFit fitExpansion = transformEdges(rExpansion, 
                        sExpansion, edges1, searchableEdges2);
                    
                    if ((fitExpansion.getChiSqSum() < fitReflected.getChiSqSum())
                    && ((rExpansion >= rMin) && (rExpansion <= rMax) 
                    && (sExpansion >= sMin) && (sExpansion <= sMax))
                    ) {

                        fits[worstFitIdx] = fitExpansion;
                        
                    } else {
                        
                        fits[worstFitIdx] = fitReflected;
                    }
                    
                } else {
                
                    // we know that the reflection fit is worse than the 2nd worse

                    // "Contraction"
                    double rContraction = r + (beta * 
                        (r - fits[worstFitIdx].getRotationInRadians()));
                    double sContraction = s + (beta * 
                        (s - fits[worstFitIdx].getScale()));
                    TransformationFit fitContraction = transformEdges(
                        rContraction, sContraction, edges1, searchableEdges2);
                
                    if (fitContraction.getChiSqSum() > fits[worstFitIdx].getChiSqSum()
                        && (
                        (rContraction >= rMin) && (rContraction <= rMax) 
                        && (sContraction >= sMin) && (sContraction <= sMax)
                        )
                    ) {

                        fits[worstFitIdx] = fitContraction;
                        
                    } else {
                        
                        // "Reduction"
                        for (int i = 1; i < fits.length; i++) {
                            
                            double rReduction = 
                                fits[bestFitIdx].getRotationInRadians()
                                + (tau * 
                                (fits[i].getRotationInRadians()
                                - fits[bestFitIdx].getRotationInRadians()));
                            double sReduction = 
                                fits[bestFitIdx].getScale()
                                + (tau * 
                                (fits[i].getScale()
                                - fits[bestFitIdx].getScale()));
                            fits[i] = transformEdges(rReduction, sReduction, 
                                edges1, searchableEdges2);
                        }                        
                    }
                }
            }

            log.finest("best fit so far: " + fits[bestFitIdx].getChiSqSum());
            
            nIter++;

            if (fits[bestFitIdx].getChiSqSum() < convergence) {
                go = false;
            } else if ((r > rMax) || (r < rMin)) {
                go = false;
            } else if ((s > sMax) || (s < sMin)) {
                go = false;
            }
        }
        
        return fits[0].getParameters();
    }
    
    private TransformationFit transformEdges(double rotInRad, double scl, 
        PairIntArray[] edges1, SearchableCurve[] edges2) {
        
        Transformer transformer = new Transformer();
        
        edges1 = transformer.applyTransformation(rotInRad, scl, 0, 0, edges1);
        
        TransformationParameters params = new TransformationParameters();
        params.setScale((float)scl);
        params.setRotationInRadians((float)rotInRad);
                
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        double offsetX = 0;
        double offsetY = 0;
        for (int i = 0; i < edges1.length; i++) {
            
            double[] centroidXY1 = curveHelper.calculateXYCentroids(edges1[i]);
            double centroidX1 = centroidXY1[0];
            double centroidY1 = centroidXY1[1];
            
            double[] centroidXY2 = curveHelper.calculateXYCentroids(edges2[i]);
            double centroidX2 = centroidXY2[0];
            double centroidY2 = centroidXY2[1];
            
            offsetX += (centroidX2 - centroidX1);
            offsetY += (centroidY2 - centroidY1);
        }
        offsetX /= (double)edges1.length;
        offsetY /= (double)edges1.length;
        
        // ==== apply translation to points in edges1 =====
        for (int i = 0; i < edges1.length; i++) {            
            PairIntArray edge1 = edges1[i];
            for (int j = 0; j < edge1.getN(); j++) {
                int x = edge1.getX(j) + (int)offsetX;
                int y = edge1.getY(j) + (int)offsetY;
                edge1.set(j, x, y);
            }
        }
        //Image img1 = new Image(300, 300);
        //debugDisplay(edges1, edges2, img1, "after translation");
        
        //==== write the curve and compare it to the model, edge2 ====
        
        TransformationFit fit = new TransformationFit(params, edges1);
        
        double diffSum = differenceOfCurves(edges1, edges2);
        
        fit.setChiSqSum(diffSum);
        
        return fit;
    }
    
    /**
     * calculate the differences of the points in edges1 to their closest 
     * points in edges2.
     * 
     * @param edges1
     * @param edges2
     * @return 
     */
    protected double differenceOfCurves(PairIntArray[] edges1, 
        SearchableCurve[] edges2) {
        
        double sumDiff = 0;
        
        for (int i = 0; i < edges1.length; i++) {
            
            double diff = differenceOfCurves(edges1[i], edges2[i]);
                
            sumDiff += diff;
        }
        
        return sumDiff;
    }
    
    /**
     * calculate the differences of the points in edge1 to their closest 
     * points in edge2.
     * 
     * @param edge1
     * @param edge2
     * @return 
     */
    protected double differenceOfCurves(PairIntArray edge1, SearchableCurve 
        edge2) {
        
        double sumDiff = 0;
        for (int i = 0; i < edge1.getN(); i++) {
            
            int x1 = edge1.getX(i);
            int y1 = edge1.getY(i);
            
            PairInt xy2 = edge2.findClosestMatch(x1, y1);
            
            double diff;
            
            if (xy2 != null) {
                int x2 = xy2.getX();
                int y2 = xy2.getY();
                diff = ((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2));
            } else {
                // penalty is full cost of point set
                diff = (x1 * x1) + (y1 * y1);
            }
            
            sumDiff += diff;
        }
        
        return Math.sqrt(sumDiff);
    }
    
    /**
     * sort the array fits by ascending chisquare sum using
     * the quick sort algorithm.
     *
     * @param yfits
     * @param p first index of the array fits to be sorted
     * @param r the last index of the array fits to be sorted, inclusive
     */
    void sortFromMinToMax(TransformationFit[] fits, int p, int r) {

        if (p < r) {

            int q = partition(fits, p, r);

            sortFromMinToMax(fits, p, q - 1);

            sortFromMinToMax(fits, q + 1, r);
        }
    }
    
    /**
     * the partition function of the TransformationFit quick sort method.
     *
     * @param yfits
     * @param p
     * @param r
     * @return
     */
    int partition(TransformationFit[] fits, int p, int r) {

        double xxp = fits[r].getChiSqSum();

        int i = p - 1;

        for (int j = p; j < r ; j++ ) {
            if (fits[j].getChiSqSum() <= xxp) {

                i++;

                TransformationFit swap = fits[i];
                fits[i] = fits[j];
                fits[j] = swap;
            }
        }

        i++;
        TransformationFit swap = fits[i];
        fits[i] = fits[r];
        fits[r] = swap;

        return i;
    }
}
