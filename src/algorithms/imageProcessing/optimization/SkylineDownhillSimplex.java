package algorithms.imageProcessing.optimization;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.SkylineExtractor;
import algorithms.util.PairInt;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class SkylineDownhillSimplex {
    
    private final List<ImageExt> images;
    private final List<GreyscaleImage> thetaImages;
    private final List<Set<PairInt>> seedPoints;
    private final List<Set<PairInt>> excludePoints;
    private final List<Set<PairInt>> expectedPoints;
    private final ANDedClauses[] clauses;
    private final float[][] coeffLowerLimits;
    private final float[][] coeffUpperLimits;
    private final float[][] coeffStepSizes;
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public SkylineDownhillSimplex(List<ImageExt> images, 
        List<GreyscaleImage> thetaImages, 
        List<Set<PairInt>> seedPoints, List<Set<PairInt>> excludePoints, 
        List<Set<PairInt>> expectedPoints, ANDedClauses[] clauses,
        float[][] coeffLowerLimits, float[][] coeffUpperLimits,
        float[][] coeffStepSizes) {
        
        if (images == null) {
            throw new IllegalArgumentException("images cannot be null");
        }
        if (thetaImages == null) {
            throw new IllegalArgumentException("thetaImages cannot be null");
        }
        if (seedPoints == null) {
            throw new IllegalArgumentException("seedPoints cannot be null");
        }
        if (excludePoints == null) {
            throw new IllegalArgumentException("excludePoints cannot be null");
        }
        if (expectedPoints == null) {
            throw new IllegalArgumentException("expectedPoints cannot be null");
        }
        if (clauses == null) {
            throw new IllegalArgumentException("clauses cannot be null");
        }
        
        this.images = images; 
        this.thetaImages = thetaImages;
        this.seedPoints = seedPoints;
        this.excludePoints = excludePoints;
        this.expectedPoints = expectedPoints;
        this.clauses = clauses;
        this.coeffLowerLimits = coeffLowerLimits;
        this.coeffUpperLimits = coeffUpperLimits;
        this.coeffStepSizes = coeffStepSizes;
    }
    
    protected SkylineFits process(SkylineExtractor skylineExtractor,
        ANDedClauses[] c) {

        SkylineFits fit = new SkylineFits();
        List<SetComparisonResults> results = new ArrayList<SetComparisonResults>();
        for (int i = 0; i < images.size(); i++) {
            SkylineFits tFit = process(skylineExtractor, c, i);
            SetComparisonResults result = tFit.results;
            results.add(result);
        }
        SetComparisonResults result = new SetComparisonResults(results);
        fit.clauses = c;
        fit.results = result;

        return fit;
    }
    
    protected SkylineFits process(SkylineExtractor skylineExtractor,
        ANDedClauses[] clauses, int imageIndex) {

        ImageExt img = images.get(imageIndex);
        GreyscaleImage thetaImg = thetaImages.get(imageIndex);
        Set<PairInt> skyPoints = new HashSet<PairInt>(seedPoints.get(imageIndex));
        Set<PairInt> exclude = excludePoints.get(imageIndex);
        Set<PairInt> expected = expectedPoints.get(imageIndex);

        skylineExtractor.findClouds(skyPoints, exclude, img, thetaImg,
            clauses);

        SetCompare setCompare = new SetCompare();
        
        SetComparisonResults results = setCompare.compare(expected, skyPoints);

        SkylineFits fit = new SkylineFits();
        fit.results = results;
        fit.clauses = clauses;
        
        return fit;
    }
    
    public SkylineFits fit() throws NoSuchAlgorithmException {
                
        SkylineExtractor skylineExtractor = new SkylineExtractor();
                
        //TODO: edit convergence.  it's the fraction of matched to expected matches.
        float convergence = 0.98f;
      
        int nStarterPoints = 10;
        
        SkylineFits[] fits = createStarterPoints(skylineExtractor,
            nStarterPoints);
        
        int nMaxIter = 100;
        int nIter = 0;
        
        int bestFitIdx = 0;
        int secondWorstFitIdx = fits.length - 2;
        int worstFitIdx = fits.length - 1;

        SetComparisonResults lastBest = null;
        int nIterSameMin = 0;
        
        float alpha = 1;   // > 0
        float gamma = 2;   // > 1
        float beta = -0.5f; 
        float tau = 0.5f;
        
        boolean go = true;
        
        while (go && (nIter < nMaxIter)) {

            // best matches should be at smaller indexes
            Arrays.sort(fits, 0, (fits.length - 1));
            
            if ((lastBest != null) && 
                (lastBest.compareTo(fits[bestFitIdx].results) == 0)) {
                
                nIterSameMin++;
                
                if (nIterSameMin >= 5) {
                    break;
                }
            } else {
                nIterSameMin = 0;
            }
            lastBest = fits[bestFitIdx].results;
            
            float[][] summedCoeff = sumAllButLastCoefficients(fits);
            
            //"Reflection"
            ANDedClauses[] tClauses = reflect(fits[worstFitIdx], summedCoeff, 
                alpha);
            
            SkylineFits fitReflected = process(skylineExtractor, tClauses);
                        
            boolean relectIsWithinBounds = isWithinBounds(fitReflected);
            
            if (fitIsBetter(fits[secondWorstFitIdx], fitReflected) 
                && relectIsWithinBounds
                && !fitIsBetter(fits[bestFitIdx], fitReflected)) {
                
                fits[worstFitIdx] = fitReflected;
                
            } else {
                
                if (fitIsBetter(fits[bestFitIdx], fitReflected)
                    && relectIsWithinBounds) {
                    
                    // "Expansion"
                    tClauses = expand(fits[worstFitIdx], summedCoeff, gamma);
                    
                    SkylineFits fitExpansion = process(skylineExtractor, tClauses);
                    
                    boolean expandIsWithinBounds = isWithinBounds(fitExpansion);
                    
                    if (fitIsBetter(fitReflected, fitExpansion)
                        && expandIsWithinBounds) {

                        fits[worstFitIdx] = fitExpansion;
                        
                    } else {
                        
                        fits[worstFitIdx] = fitReflected;
                    }
                    
                } else {
                
                    // we know that the reflection fit is worse than the 2nd worse

                    // "Contraction"
                    tClauses = contract(fits[worstFitIdx], summedCoeff, beta);
                    
                    SkylineFits fitContraction = process(skylineExtractor, tClauses);
                    
                    boolean IsWithinBounds = isWithinBounds(fitContraction);
                    
                    if (fitIsBetter(fits[worstFitIdx], fitContraction)
                        && IsWithinBounds) {

                        fits[worstFitIdx] = fitContraction;
                        
                    } else {
                                                
                        // "Reduction"
                        for (int i = 1; i < fits.length; i++) {
                            
                            tClauses = reduce(fits[bestFitIdx], fits[i], tau);
                            
                            SkylineFits fitReduction = process(skylineExtractor, 
                                tClauses);
                            
                            //TODO: consider what to do when solution is out of
                            // bounds
                            fits[i] = fitReduction;
                        }
                    }
                }
            }

            log.finest("best fit so far: " + 
                fits[bestFitIdx].results.numberMatchedDivExpected);
            
            nIter++;
            
            //TODO: test the bounds of the best answer

            // convergence is when overrun is zero and number matched is within
            // some limit
            if ((fits[bestFitIdx].results.numberOverrunDivExpected == 0) &&
                (fits[bestFitIdx].results.numberMatchedDivExpected == convergence)) {
                go = false;
            }
        }
        
        return fits[bestFitIdx];
    }

    private SkylineFits[] createStarterPoints(SkylineExtractor skylineExtractor,
        int nStarterPoints) throws NoSuchAlgorithmException {
                          
        /*
        what pattern to use to alter ANDedClauses[]?
        -- unaltered ANDedClauses[]
        -- ANDedClauses[] adjusted to all lower limits.
        -- ANDedClauses[] adjusted all all higher limits.
        -- 4 divisions between the higher and lower limits?
        
        that's 7 "starter points", but the values are not representing the
        possible mix of ranges very well, so might need to sample that randomly
        by mixing the combinations of the modified coefficients and increase
        the number of starter points to 20 or so.  
        for example, to avoid all of the lowest bounds being present in one
        set of clauses, a random mix would put some of them in each possibly.
        
        the reason for mixing this small set of coefficients is that the 
        simplex, when it alters the coefficients, performs the same changes
        on all of the coefficients in that starter point's list of coefficients, 
        so an "all low bounds" when altered
        might not find the local best for all 40 something coefficients.
        
        With so many parameters, the Nelder-Maede seems like one of the only
        feasible choices for optimization, but it requires more starter points 
        and good bounds to succeed for this use.
        */
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(System.currentTimeMillis());
        
        ANDedClauses[][] starterPointClauses = new ANDedClauses[nStarterPoints][];
        starterPointClauses[0] = Arrays.copyOf(clauses, clauses.length);
        
        for (int i = 1; i < nStarterPoints; i++) {
            starterPointClauses[i] = Arrays.copyOf(clauses, clauses.length);
            for (int ii = 0; ii < starterPointClauses.length; ii++) {
                ANDedClauses clause = starterPointClauses[i][ii];
                for (int jj = 0; jj < clause.coefficients.length; jj++) {
                    float coeff = clause.coefficients[jj];
                    // pick randomly between low bounds, high bounds, center,
                    // or random between bounds
                    int type = sr.nextInt(4);
                    switch(type) {
                        case 0:
                            clause.coefficients[jj] = coeffLowerLimits[ii][jj];
                            break;
                        case 1:
                            clause.coefficients[jj] = coeffUpperLimits[ii][jj];
                            break;
                        case 2:
                            //remain the same
                            break;
                        default: {
                            float range = coeffUpperLimits[ii][jj] -
                                coeffLowerLimits[ii][jj];
                            // dividing the range by 3 and choosing randomly
                            // between those 2 middle marks
                            int d = sr.nextInt(3);
                            clause.coefficients[jj] = (float)(coeffLowerLimits[ii][jj] +
                                ((float)(d + 1.))*(range/3.));
                            break;
                        }
                    }
                }
            }
        }        
        
        SkylineFits[] fits = new SkylineFits[nStarterPoints];
        for (int idx = 0; idx < fits.length; idx++) { 
            SkylineFits fit = process(skylineExtractor, starterPointClauses[idx]);
            fits[idx] = fit;
        }
        
        return fits;
    }

    private float[][] sumAllButLastCoefficients(SkylineFits[] fits) {
        
        int n = 0;
        SkylineFits fit = fits[0];
        int nClauses = fit.clauses.length;
        float[][] avgCoeffs = new float[nClauses][];
        for (int clauseIdx = 0; clauseIdx < nClauses; clauseIdx++) {
            
            int nCoefficients = clauses[clauseIdx].coefficients.length;
            
            avgCoeffs[clauseIdx] = new float[nCoefficients];
            
            for (int clauseCoeffIdx = 0; clauseCoeffIdx < nCoefficients; 
                clauseCoeffIdx++) {
                
                // iterate over each clauseIdx and each clauseCoeffIdx
                // for each fit in fits.
                double sumCoeff = 0;
                for (int fitIdx = 0; fitIdx < fits.length; fitIdx++) {
                    fit = fits[fitIdx];
                    float coeff = fit.clauses[clauseIdx].coefficients[clauseCoeffIdx];
                    sumCoeff += coeff;
                }
                avgCoeffs[clauseIdx][clauseCoeffIdx] = (float)(sumCoeff/(float)fits.length);
            }
        }
        
        return avgCoeffs;
    }

    /**
     * perform "reflection" to change each coefficient using pattern
     *     averagedCoeff[i] + (alpha * (averagedCoeff[i] - fit.coefficient)).
     *
     * @param skylineFits
     * @param summedCoeff
     * @param alpha
     * @return 
     */
    private ANDedClauses[] reflect(SkylineFits fit, float[][] averagedCoeff, 
        float alpha) {
        
        ANDedClauses[] tClauses = performAction(fit, averagedCoeff, alpha);
        
        return tClauses;
    }

    /**
     * perform "expansion" to change each coefficient using pattern
     *     averagedCoeff[i] + (gamma * (averagedCoeff[i] - fit.coefficient)).
     *
     * @param skylineFits
     * @param averagedCoeff
     * @param alpha
     * @return 
     */
    private ANDedClauses[] expand(SkylineFits fit, 
        float[][] averagedCoeff, float gamma) {
        
        ANDedClauses[] tClauses = performAction(fit, averagedCoeff, gamma);
        
        return tClauses;
    }
    
    /**
     * perform "contraction" to change each coefficient using pattern
     *     averagedCoeff[i] + (beta * (averagedCoeff[i] - fit.coefficient)).
     *
     * @param skylineFits
     * @param averagedCoeff
     * @param alpha
     * @return 
     */
    private ANDedClauses[] contract(SkylineFits fit, 
        float[][] averagedCoeff, float beta) {
        
        ANDedClauses[] tClauses = performAction(fit, averagedCoeff, beta);
        
        return tClauses;
    }
    
    /**
     * perform "contraction" to change each coefficient using pattern
     *    bestFit + (tau * (indivFit - bestFit)).
     * @param skylineFits
     * @param averagedCoeff
     * @param alpha
     * @return 
     */
    private ANDedClauses[] reduce(SkylineFits bestFit, SkylineFits fitI, 
        float tau) {
        
        ANDedClauses[] bestClauses = Arrays.copyOf(bestFit.clauses, 
            bestFit.clauses.length);
         
        for (int clauseIdx = 0; clauseIdx < bestClauses.length; clauseIdx++) {
            
            ANDedClauses clauseBest = bestClauses[clauseIdx];
            
            ANDedClauses clauseI = fitI.clauses[clauseIdx];
            
            for (int clauseCoeffIdx = 0; clauseCoeffIdx < 
                clauseBest.coefficients.length; clauseCoeffIdx++) {
                
                float coeffBest = clauseBest.coefficients[clauseCoeffIdx];
                
                float coeffI = clauseI.coefficients[clauseCoeffIdx];
                
                float coeffM = coeffBest + (tau * (coeffI - coeffBest));
                
                bestClauses[clauseIdx].coefficients[clauseCoeffIdx] = coeffM;
            }
        }
        
        return bestClauses;
    }
    
    private ANDedClauses[] performAction(SkylineFits fit, float[][] averagedCoeff, 
        float factor) {
        
        ANDedClauses[] tClauses = Arrays.copyOf(fit.clauses, fit.clauses.length);
         
        for (int clauseIdx = 0; clauseIdx < tClauses.length; clauseIdx++) {
            
            ANDedClauses clause = tClauses[clauseIdx];
            
            for (int clauseCoeffIdx = 0; clauseCoeffIdx < 
                clause.coefficients.length; clauseCoeffIdx++) {
                
                float coeff = averagedCoeff[clauseIdx][clauseCoeffIdx];
                
                float coeffM = coeff + (factor * 
                    (coeff - fit.clauses[clauseIdx].coefficients[clauseCoeffIdx]));
                
                tClauses[clauseIdx].coefficients[clauseCoeffIdx] = coeffM;
            }
        }
        
        return tClauses;
    }
    
    private boolean isWithinBounds(SkylineFits fit) {
        
        ANDedClauses[] tClauses = fit.clauses;
         
        for (int clauseIdx = 0; clauseIdx < tClauses.length; clauseIdx++) {
            
            ANDedClauses clause = tClauses[clauseIdx];
            
            for (int clauseCoeffIdx = 0; clauseCoeffIdx < 
                clause.coefficients.length; clauseCoeffIdx++) {
                
                float coeff = clause.coefficients[clauseCoeffIdx];
                
                float upper = this.coeffUpperLimits[clauseIdx][clauseCoeffIdx];
                
                float lower = this.coeffLowerLimits[clauseIdx][clauseCoeffIdx];
                
                if ((coeff < lower) || (coeff > upper)) {
                    return false;
                }
            }
        }
        
       return true;
    }

    private boolean fitIsBetter(SkylineFits fit, SkylineFits compareToFit) {
        
        // -1 is fit is bettern than compareToFit
        int comp = fit.compareTo(compareToFit);
        
        return (comp > 0);
    }

}
