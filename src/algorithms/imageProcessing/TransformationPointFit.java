package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class TransformationPointFit implements Comparable<TransformationPointFit> {
    
    private final TransformationParameters parameters;
    
    private final int nMatchedPoints;
    
    private int nMaxMatchable = 0;
    
    private final double meanDistFromModel;
    
    private final double stDevFromMean;
    
    private float transTolX = Float.MAX_VALUE;
    
    private float transTolY = Float.MAX_VALUE;
    
    public TransformationPointFit(
        TransformationParameters theParameters, 
        int numberOfMatchedPoints, double theMeanDistFromModel, 
        double theStDevFromMean, float theTranslationXTolerance,
        float theTranslationYTolerance) {
        
        parameters = theParameters;
        
        nMatchedPoints = numberOfMatchedPoints;
        
        meanDistFromModel = theMeanDistFromModel;
        
        stDevFromMean = theStDevFromMean;
        
        transTolX = theTranslationXTolerance;
        
        transTolY = theTranslationYTolerance;
    }

    /**
     * @return the parameters
     */
    public TransformationParameters getParameters() {
        return parameters;
    }
    
    public int getNumberOfMatchedPoints() {
        return nMatchedPoints;
    }
    
    public float getRotationInRadians() {
        if (parameters == null) {
            return 0;
        }
        
        return parameters.getRotationInRadians();
    }
    
    public float getScale() {
        if (parameters == null) {
            return 0;
        }
        
        return parameters.getScale();
    }
    
    public float getTranslationX() {
        if (parameters == null) {
            return 0;
        }
        
        return parameters.getTranslationX();
    }
    
    public float getTranslationY() {
        if (parameters == null) {
            return 0;
        }
        
        return parameters.getTranslationY();
    }
    
    /**
     * tolerance used when including only residuals whose absolute value
     * is less than tolerance.
     * @return 
     */
    public float getTranslationXTolerance() {
        return transTolX;
    }
    
    /**
     * tolerance used when including only residuals whose absolute value
     * is less than tolerance.
     * @return 
     */
    public float getTranslationYTolerance() {
        return transTolY;
    }

    /**
     * @return the meanDistFromModel
     */
    public double getMeanDistFromModel() {
        return meanDistFromModel;
    }

    /**
     * @return the stDevFromMean
     */
    public double getStDevFromMean() {
        return stDevFromMean;
    }
    
    /**
     * @return the nMaxMatchable
     */
    public int getNMaxMatchable() {
        return nMaxMatchable;
    }

    public void setTranslationXTolerance(float theTolerance) {
        transTolX = theTolerance;
    }
    
    public void setTranslationYTolerance(float theTolerance) {
        transTolY = theTolerance;
    }
    
    /**
     * @param maximumNumberMatchable the nMaxMatchable to set
     */
    public void setMaximumNumberMatchable(int maximumNumberMatchable) {
        this.nMaxMatchable = maximumNumberMatchable;
    }
    
    /**
     * compare object other to this instance and return -1 if this instance
     * is "better", 0 if they are equal, else 1 if other is "better".
     * @param other
     * @return 
     */
    @Override
    public int compareTo(TransformationPointFit other) {
    
        if (other == null) {
            return -1;
        }
        if ((other.getNumberOfMatchedPoints() == 0) && (nMatchedPoints > 0)) {
            return -1;
        }
        if (Double.isNaN(other.getStDevFromMean()) && !Double.isNaN(stDevFromMean)) {
            return -1;
        }
        if ((other.getStDevFromMean() == Double.MAX_VALUE) && 
            (stDevFromMean != Double.MAX_VALUE)) {
            return -1;
        }
        if ((other.getNumberOfMatchedPoints() > 0) && (nMatchedPoints == 0)) {
            return 1;
        }
        if (!Double.isNaN(other.getStDevFromMean()) && Double.isNaN(stDevFromMean)) {
            return 1;
        }
        if ((other.getStDevFromMean() != Double.MAX_VALUE) && 
            (stDevFromMean == Double.MAX_VALUE)) {
            return 1;
        }
        
        TransformationPointFit compareFit = other;
        
        int compNMatches = compareFit.getNumberOfMatchedPoints();
        int bestNMatches = nMatchedPoints;

        double compAvg = compareFit.getMeanDistFromModel();
        double bestAvg = meanDistFromModel;
        double diffAvg = Math.abs(compAvg - bestAvg);

        double compS = compareFit.getStDevFromMean();
        double bestS = stDevFromMean;
        double diffS = Math.abs(compS - bestS);

        double avgDiv = bestAvg/compAvg;

        int diffEps = (int)Math.ceil(Math.max(bestNMatches, compNMatches)/10.);
        if (diffEps == 0) {
            diffEps = 1;
        }

        if (Math.abs(bestNMatches - compNMatches) <= diffEps) {

            if (Double.isNaN(compareFit.getMeanDistFromModel())) {
                return -1;
            }
            
            if (compAvg < bestAvg) {
                return 1;
            } else if (compAvg > bestAvg) {
                return -1;
            }

            if (compS < bestS) {
                return 1;
            } else if (compS > bestS) {
                return -1;
            }
            
        } else if (compNMatches > bestNMatches) {

            return 1;
        }
        
        return -1;
    }

    /**
     * compare object other to this instance and return -1 if this instance
     * is "better", 0 if they are equal, else 1 if other is "better".
     * Uses the same logic as compareTo, except replaces nMatched with
     * nMatched/nMaxMatchable.
     * @param other
     * @return 
     */
    public int compareToUsingNormalizedMatches(TransformationPointFit other) {
    
        if (other == null) {
            return -1;
        }
        if ((other.getNumberOfMatchedPoints() == 0) && (nMatchedPoints > 0)) {
            return -1;
        }
        if (Double.isNaN(other.getStDevFromMean()) && !Double.isNaN(stDevFromMean)) {
            return -1;
        }
        if ((other.getStDevFromMean() == Double.MAX_VALUE) && 
            (stDevFromMean != Double.MAX_VALUE)) {
            return -1;
        }
        if ((other.getNumberOfMatchedPoints() > 0) && (nMatchedPoints == 0)) {
            return 1;
        }
        if (!Double.isNaN(other.getStDevFromMean()) && Double.isNaN(stDevFromMean)) {
            return 1;
        }
        if ((other.getStDevFromMean() != Double.MAX_VALUE) && 
            (stDevFromMean == Double.MAX_VALUE)) {
            return 1;
        }
        
        TransformationPointFit compareFit = other;
        
        float compNMatches = (float)compareFit.getNumberOfMatchedPoints()/
            (float)compareFit.getNMaxMatchable();
        float bestNMatches = (float)nMatchedPoints/(float)nMaxMatchable;

        double compAvg = compareFit.getMeanDistFromModel();
        double bestAvg = meanDistFromModel;

        double compS = compareFit.getStDevFromMean();
        double bestS = stDevFromMean;

        double diffEps = Math.sqrt(Math.min(compareFit.getNMaxMatchable(),
            nMaxMatchable));
        diffEps = 1./(10.*diffEps);

        if (Math.abs(bestNMatches - compNMatches) <= diffEps) {

            if (Double.isNaN(compareFit.getMeanDistFromModel())) {
                return -1;
            }
            
            if (compAvg < bestAvg) {
                return 1;
            } else if (compAvg > bestAvg) {
                return -1;
            }

            if (compS < bestS) {
                return 1;
            } else if (compS > bestS) {
                return -1;
            }
            
        } else if (compNMatches > bestNMatches) {

            return 1;
        }
        
        return -1;
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("\nnMatchedPoints=").append(Integer.toString(nMatchedPoints))
            .append(" nMaxMatchable=")
            .append(Double.toString(nMaxMatchable))
            .append("\nmeanDistFromModel=")
            .append(Double.toString(meanDistFromModel))
            .append(" stDevFromMean=")
            .append(Double.toString(stDevFromMean))
            .append("\ntranslationXTolerance=")
            .append(Float.toString(transTolX))
            .append(" translationYTolerance=")
            .append(Float.toString(transTolY))
            .append("\n")
            .append(parameters.toString());
        
        return sb.toString();
    }

    public int isSimilarWithDiffParameters(TransformationPointFit compareFit) {
        
        /*if the fit and bestFit is very close,
        need an infrastructure to return more than one
        (and to reset it when another fit not similar to bestFit is found)

        bestFit:
            fit=nMatchedPoints=63 nMaxMatchable=63.0
            meanDistFromModel=37.57523582095192
            stDevFromMean=22.616214121394517
            tolerance=144.2497833620557
            rotationInRadians=0.34906584
            rotationInDegrees=19.99999941818584
            scale=1.0
            translationX=86.61143 translationY=-18.330833

        fit:
            fit=nMatchedPoints=63 nMaxMatchable=63.0
            meanDistFromModel=36.79744166041177            <-- mean dist and stdev are similar
            stDevFromMean=23.167906020816925
            tolerance=144.2497833620557
            rotationInRadians=0.5235988
            rotationInDegrees=30.000000834826057            <--- rot is very different
            scale=1.0
            translationX=91.6094 translationY=-72.9244      <--- transY is very different

        the correct answer is closer to "bestFit" but is currently
        then replaced with fit, so need to preserve both.

        correct answer:
            rotationInDegrees=18.000000500895634
            scale=1.0
            translationX=7.0 translationY=33.0
            number of vertical partitions=3

        ----
        in contrast, these 2 sets of meanDistFromModel and stDevFromMean
        are significantly different:
            bestFit: meanDistFromModel=0.3658992320831333 stDevFromMean=0.15709856816653883
            fit:     meanDistFromModel=2.267763360270432 stDevFromMean=1.0703222291650063
        -----
        so ratios are needed rather than differences
        similar fits:    divMean = 1.02
                         divStDev = 0.976
        different fits:  divMean = 0.161
                         divStDev = 0.147
        */

        //TODO: this may need to be adjusted and in the 2 fitness
        // functions which use it... should be it's own method...
        int nEps = (int)(1.5*Math.ceil(nMatchedPoints/10.));
        if (nEps == 0) {
            nEps = 1;
        }

        int diffNMatched = Math.abs(nMatchedPoints -
            compareFit.getNumberOfMatchedPoints());

        if (diffNMatched > nEps) {
            return -1;
        }

        double divMean = Math.abs(meanDistFromModel/
            compareFit.getMeanDistFromModel());

        double divStDev = Math.abs(stDevFromMean/compareFit.getStDevFromMean());

        if ((Math.abs(1 - divMean) < 0.05) && (Math.abs(1 - divStDev) < 0.3)) {
            if (parameters.equals(compareFit.getParameters())) {
                return 0;
            } else {
                return 1;
            }
        }

        return -1;
    }

    public int isSimilarNormalizedWithDiffParameters(TransformationPointFit compareFit) {
        
        float compNMatches = (float)compareFit.getNumberOfMatchedPoints()/
            (float)compareFit.getNMaxMatchable();
        float bestNMatches = (float)nMatchedPoints/(float)nMaxMatchable;
        
        double diffEps = Math.sqrt(Math.min(compareFit.getNMaxMatchable(),
            nMaxMatchable));
        diffEps = 1./(10.*diffEps);

        if (Math.abs(bestNMatches - compNMatches) > diffEps) {
            return -1;
        }

        double divMean = Math.abs(meanDistFromModel/
            compareFit.getMeanDistFromModel());

        double divStDev = Math.abs(stDevFromMean/compareFit.getStDevFromMean());

        if ((Math.abs(1 - divMean) < 0.05) && (Math.abs(1 - divStDev) < 0.3)) {
            if (parameters.equals(compareFit.getParameters())) {
                return 0;
            } else {
                return 1;
            }
        }

        return -1;
    }
}
