package algorithms.imageProcessing.optimization;

import java.util.List;

/**
 *
 * @author nichole
 */
public class SetComparisonResults implements Comparable<SetComparisonResults> {

    protected final int nExpectedPoints;
    protected final float numberOverrunDivExpectedMatchedPoints;
    protected final float numberMatchedDivExpected;    
    protected final int nExpectedBorderPoints;
    protected final float numberMatchedBorderDivExpected;
    public final double eps0;
    
    private final float eps0Factor = 1000.f;
    //private final float eps0Factor = 10000.f;
    
    public SetComparisonResults(int nExpected, int nOverrun, int nMatched,
        int nExpectedBorder, int nMatchedBorder) {
        
        this.nExpectedPoints = nExpected;
        this.numberOverrunDivExpectedMatchedPoints = (nExpected > 0) ?
            (float) nOverrun / (float) nExpected : Float.POSITIVE_INFINITY;
        this.numberMatchedDivExpected = (nExpected > 0) ? 
            (float) nMatched / (float) nExpected : Float.POSITIVE_INFINITY;
        this.nExpectedBorderPoints = nExpectedBorder;
        this.numberMatchedBorderDivExpected = (nExpected > 0) ? 
            (float) nMatchedBorder / (float) nExpectedBorder 
            : Float.POSITIVE_INFINITY;
        
        System.out.println("eps0 is " 
            + (eps0Factor/(float)nExpectedPoints) + " fraction of total");
        
        this.eps0 = eps0Factor * (1./(double)nExpectedPoints);
    }
    
    private double numberOverrunDivExpectedMatchedPointsStDev = 0;
    private double numberMatchedDivExpectedStDev = 0;
    private double numberMatchedBorderDivExpectedStDev = 0;
    public SetComparisonResults(List<SetComparisonResults> results) {
        
        int nTotal = 0;
        for (SetComparisonResults result : results) {
            nTotal += result.nExpectedPoints;
        }
        this.nExpectedPoints = nTotal;
        
        float sumOverrunFraction = 0;
        for (SetComparisonResults result : results) {
            double fraction = result.numberOverrunDivExpectedMatchedPoints;
            sumOverrunFraction += fraction;
        }
        this.numberOverrunDivExpectedMatchedPoints = 
            sumOverrunFraction/(float)results.size();
        
        double sumstdev = 0;
        for (SetComparisonResults result : results) {
            double fraction = result.numberOverrunDivExpectedMatchedPoints;
            double diff = fraction - numberOverrunDivExpectedMatchedPoints;
            sumstdev += (diff * diff);
        }
        numberOverrunDivExpectedMatchedPointsStDev = (Math.sqrt(sumstdev/
            (float)(results.size() - 1.0)));
        
        
        float sumnMatchedFraction = 0;
        for (SetComparisonResults result : results) {
            double fraction = result.numberMatchedDivExpected;
            sumnMatchedFraction += fraction;
        }
        this.numberMatchedDivExpected = sumnMatchedFraction/(float)results.size();
        
        sumstdev = 0;
        for (SetComparisonResults result : results) {
            double fraction = result.numberMatchedDivExpected;
            double diff = fraction - numberMatchedDivExpected;
            sumstdev += (diff * diff);
        }
        this.numberMatchedDivExpectedStDev = (Math.sqrt(sumstdev/
            (float)(results.size() - 1.0)));
        
        
        nTotal = 0;
        for (SetComparisonResults result : results) {
            nTotal += result.nExpectedBorderPoints;
        }
        this.nExpectedBorderPoints = nTotal;
        sumnMatchedFraction = 0;
        for (SetComparisonResults result : results) {
            double fraction = result.numberMatchedBorderDivExpected;
            sumnMatchedFraction += fraction;
        }
        this.numberMatchedBorderDivExpected = sumnMatchedFraction/(float)results.size();
        
        sumstdev = 0;
        for (SetComparisonResults result : results) {
            double fraction = result.numberMatchedBorderDivExpected;
            double diff = fraction - numberMatchedDivExpected;
            sumstdev += (diff * diff);
        }
        this.numberMatchedBorderDivExpectedStDev = (Math.sqrt(sumstdev/
            (float)(results.size() - 1.0)));
        
        this.eps0 = eps0Factor * (1./(double)nExpectedPoints);
    }
    
    public SetComparisonResults(List<Integer> nExpected, List<Integer> nOverrun, 
        List<Integer> nMatched, List<Integer> nExpectedBorder, 
        List<Integer> nMatchedBorder) {
        
        if (nExpected.size() != nOverrun.size() 
            || nOverrun.size() != nMatched.size() 
            || nMatched.size() != nExpectedBorder.size()
            || nExpectedBorder.size() != nMatchedBorder.size()) {
            throw new IllegalArgumentException("lists must be same size");
        }
        
        int nTotal = 0;
        for (Integer n : nExpected) {
            nTotal += n.intValue();
        }
        this.nExpectedPoints = nTotal;
        
        float sumOverrunFraction = 0;
        for (int i = 0; i < nOverrun.size(); i++) {
            double fraction = (float)nOverrun.get(i)/(float)nExpected.get(i);
            sumOverrunFraction += fraction;
        }
        this.numberOverrunDivExpectedMatchedPoints = sumOverrunFraction/(float)nOverrun.size();
        
        float sumnMatchedFraction = 0;
        for (int i = 0; i < nMatched.size(); i++) {
            double fraction = (float)nMatched.get(i)/(float)nExpected.get(i);
            sumnMatchedFraction += fraction;
        }
        this.numberMatchedDivExpected = sumnMatchedFraction/(float)nMatched.size();
        
        nTotal = 0;
        for (Integer n : nExpectedBorder) {
            nTotal += n.intValue();
        }
        this.nExpectedBorderPoints = nTotal;
        
        sumnMatchedFraction = 0;
        for (int i = 0; i < nMatchedBorder.size(); i++) {
            double fraction = (float)nMatchedBorder.get(i)/(float)nExpectedBorder.get(i);
            sumnMatchedFraction += fraction;
        }
        this.numberMatchedBorderDivExpected = 
            sumnMatchedFraction/(float)nMatchedBorder.size();
        
        
        double sumstdev = 0;
        for (int i = 0; i < nOverrun.size(); i++) {
            double fraction = (float)nOverrun.get(i)/(float)nExpected.get(i);
            double diff = fraction - numberMatchedBorderDivExpected;
            sumstdev += (diff * diff);
        }
        this.numberOverrunDivExpectedMatchedPointsStDev = (Math.sqrt(sumstdev/
            (float)(nOverrun.size() - 1.0)));
        
        sumstdev = 0;
        for (int i = 0; i < nMatched.size(); i++) {
            double fraction = (float)nMatched.get(i)/(float)nExpected.get(i);
            double diff = fraction - numberMatchedDivExpected;
            sumstdev += (diff * diff);
        }
        this.numberMatchedDivExpectedStDev = (Math.sqrt(sumstdev/
            (float)(nMatched.size() - 1.0)));
        
        sumstdev = 0;
        for (int i = 0; i < nMatchedBorder.size(); i++) {
            double fraction = (float)nMatchedBorder.get(i)/(float)nExpectedBorder.get(i);
            double diff = fraction - numberMatchedBorderDivExpected;
            sumstdev += (diff * diff);
        }
        this.numberMatchedBorderDivExpectedStDev = (Math.sqrt(sumstdev/
            (float)(nMatchedBorder.size() - 1.0)));
        
        
        this.eps0 = eps0Factor * (1./(double)nExpectedPoints);
    }

    @Override
    public boolean equals(Object obj) {

        if (!(obj instanceof SetComparisonResults)) {
            return false;
        }
        SetComparisonResults other = (SetComparisonResults) obj;
        
        double eps1 = 5.0E-4;//0.536e-7; // 1 pixel in 1024x104:
        
        if (Math.abs(numberOverrunDivExpectedMatchedPoints - other.numberOverrunDivExpectedMatchedPoints)
            > eps0) {
            return false;
        }
        if (Math.abs(numberMatchedBorderDivExpected - other.numberMatchedBorderDivExpected)
            > eps1) {
            return false;
        }
        if (Math.abs(numberMatchedDivExpected - other.numberMatchedDivExpected)
            > eps1) {
            return false;
        }
        return true;
    }

    /**
     * Compare this instance's numberOverrunDivExpected, 
     * numberMatchedBorderDivExpected, numberMatchedDivExpected to other and 
     * return -1 if this instance has better fields as a result, else 0 if they 
     * are equal, else +1 if other has better results. The sign convention is 
     * meant to put the best answers at the top of a list sorted using this 
     * comparison function, where the default behavior by java framework 
     * algorithms is to make an ascending sort. 
     * Best is defined as having the smallest numberOverrunDivExpected and
     * ties are broken by having the largest numberMatchedBorderDivExpected,
     * else if there are no overruns and no border points matched, the
     * instance with the largest numberMatchedDivExpected is the best.
     *
     * @param other
     * @return
     */
    @Override
    public int compareTo(SetComparisonResults other) {
        
        //return compareToMaximizeBorderMatches(other);
        //return compareToMaximizeMatches(other);
        return compareToMinimizeOverruns(other);
        //return compareToMinimizeOverruns2(other);
    }
    
    protected int compareToMinimizeOverruns(SetComparisonResults other) {
        
        // this one is close, but needs to finish the last approach towards skyline
        
        double eps1 = 5.0E-4;//0.536e-7; // 1 pixel in 1024x104:
        
        float diff0 = Math.abs(numberOverrunDivExpectedMatchedPoints
            - other.numberOverrunDivExpectedMatchedPoints);

        if (diff0 <= eps0) {

            float diff1 = Math.abs(numberMatchedBorderDivExpected
                - other.numberMatchedBorderDivExpected);
            
            if (diff1 <= eps1) {
                
                float diff2 = Math.abs(numberMatchedDivExpected
                    - other.numberMatchedDivExpected);

                if (diff2 <= eps1) {
                    return 0;
                } else {
                    if (numberMatchedDivExpected < other.numberMatchedDivExpected) {
                        return 1;
                    } else {
                        return -1;
                    }
                }
            } else {
                if (numberMatchedBorderDivExpected < other.numberMatchedBorderDivExpected) {
                    return 1;
                } else {
                    return -1;
                }
            }
            
        } else {
            if (numberOverrunDivExpectedMatchedPoints < other.numberOverrunDivExpectedMatchedPoints) {
                return -1;
            } else {
                return 1;
            }
        }
    }
    
    protected int compareToMinimizeOverruns2(SetComparisonResults other) {
        
        // this one is close, but needs to finish the last approach towards skyline
        
        double eps1 = 5.0E-4;//0.536e-7; // 1 pixel in 1024x104:
        
        float diff0 = Math.abs(numberOverrunDivExpectedMatchedPoints
            - other.numberOverrunDivExpectedMatchedPoints);

        if (diff0 <= eps0) {
            
            // closest to fraction 1 for matching expected
            
            if (Math.abs(numberMatchedDivExpected - 1.0) < 
                Math.abs(other.numberMatchedDivExpected - 1.0)) {
                return -1;
            } else {
                return 1;
            }
            
        } else {
            if (numberOverrunDivExpectedMatchedPoints < other.numberOverrunDivExpectedMatchedPoints) {
                return -1;
            } else {
                return 1;
            }
        }
    }

    protected int compareToMaximizeBorderMatches(SetComparisonResults other) {
        
        double eps1 = 5.0E-4;//0.536e-7; // 1 pixel in 1024x104:
        
        float diff1 = Math.abs(numberMatchedBorderDivExpected
                - other.numberMatchedBorderDivExpected);
            
        if ((diff1 <= 0 /*eps1*/) || 
            ((numberMatchedBorderDivExpected < eps1) && 
            (other.numberMatchedBorderDivExpected < eps1))) {

            // either the number of border matches is equal or both are zero
            float diff2 = Math.abs(numberMatchedDivExpected
                - other.numberMatchedDivExpected);

            if (diff2 <= eps1) {
                
                float diff0 = Math.abs(numberOverrunDivExpectedMatchedPoints
                    - other.numberOverrunDivExpectedMatchedPoints);
                
                if (diff0 <= eps0) {
                    
                    return 0;
                    
                } else {
                    
                    // minimize the number of overruns
                    if (numberOverrunDivExpectedMatchedPoints < other.numberOverrunDivExpectedMatchedPoints) {
                        return -1;
                    } else {
                        return 1;
                    }
                }
                                
            } else {
                
                if (numberMatchedDivExpected < other.numberMatchedDivExpected) {
                    return 1;
                } else {
                    return -1;
                }
            }
            
        } else {
            
            if (numberMatchedBorderDivExpected < other.numberMatchedBorderDivExpected) {
                return 1;
            } else {
                return -1;
            }
        } 
    }

    protected int compareToMaximizeMatches(SetComparisonResults other) {
        
        double eps1 = 5.0E-4;//0.536e-7; // 1 pixel in 1024x104:
        
        float diff0 = Math.abs(numberOverrunDivExpectedMatchedPoints
            - other.numberOverrunDivExpectedMatchedPoints);

        float diff1 = Math.abs(numberMatchedDivExpected
            - other.numberMatchedDivExpected);
        
        if (diff1 <= eps1) {
            
            // the number matched are roughly equal

            if (diff0 <= eps0) {
                
                // maximize the number of border points matched
                if (numberMatchedBorderDivExpected < other.numberMatchedBorderDivExpected) {
                    return 1;
                } else {
                    return -1;
                }
                
            } else {
                
                // number of matched were equal, return the one with fewest
                // overruns
                
                if (numberOverrunDivExpectedMatchedPoints < other.numberOverrunDivExpectedMatchedPoints) {
                    return -1;
                } else {
                    return 1;
                }
            }
            
        } else {
            
            // want the largest number of matches closest to the fraction of 1.0
            
            if ((numberMatchedDivExpected < 1.0) && (other.numberMatchedDivExpected < 1.0)) {
                if (numberMatchedDivExpected > other.numberMatchedDivExpected) {
                    return -1;
                } 
                return 1;
            } else if ((numberMatchedDivExpected < 1.0) && (other.numberMatchedDivExpected > 1.0)) {
                return -1;
            }
            
            return 1;
            
            /*
            if (Math.abs(numberMatchedDivExpected - 1.0) < 
                Math.abs(other.numberMatchedDivExpected - 1.0)) {
                return -1;
            } else {
                return 1;
            }*/
        } 
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("nExpectedPoints=").append(Integer.toString(nExpectedPoints));
        
        sb.append("\n  numberOverrunDivExpectedMatchedPoints=")
            .append(Float.toString(numberOverrunDivExpectedMatchedPoints))
            .append(" (+-")
            .append(Double.toString(numberOverrunDivExpectedMatchedPointsStDev))
            .append(")");
        
        sb.append("\n  numberMatchedDivExpected=")
            .append(Float.toString(numberMatchedDivExpected))
            .append(" (+-")
            .append(Double.toString(numberMatchedDivExpectedStDev))
            .append(")");
        
        sb.append("\n  nExpectedBorderPoints=")
            .append(Integer.toString(nExpectedBorderPoints));
        
        sb.append("\n  numberMatchedBorderDivExpected=")
            .append(Float.toString(numberMatchedBorderDivExpected))
            .append(" (+-")
            .append(Double.toString(numberMatchedBorderDivExpectedStDev))
            .append(")");
        
        return sb.toString();
    }
    
}
