package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.List;

/**
 * a class to hold methods for calculating the differences
 * between color properties of labeled classes.
 * 
 * note that there is a method regarding LDA in MatrixUtil that, but it was
 * written quickly and not used or tested much since.
 * 
 * @author nichole
 */
public class BetweenClassColorStats {
    
    /**
     * 
     * @param meanLAB [class idx][color idx][avg, std dev]
     * @param deltaELAB [class idx][deltaE avg, std dev]
     * @return 
     */
    public AllClassInterStats calculateDeltaEBetweenClasses(
        float[][][] meanLAB, float[][] deltaELAB) {
        
        AllClassInterStats allStats = new AllClassInterStats();
        
        int nClasses = meanLAB.length;
        int nBands = meanLAB[0].length;
        
        if (nClasses != deltaELAB.length) {
            throw new IllegalArgumentException("meanLAB and deltaELAB"
                + " must be same length");
        }

        CIEChromaticity cieC = new CIEChromaticity();

        int nComb = nClasses - 1;
        
        for (int i = 0; i < nClasses; ++i) {
            
            OneClassInterStats stats = new OneClassInterStats();
            stats.totDiffs = new float[nComb];
            stats.signifOfDiffs = new float[nComb];
            allStats.stats.add(stats);
            
            float stdvDeltaE = deltaELAB[i][1];
            float avgDeltaE = deltaELAB[i][0];
            
            int count = 0;
            for (int j = 0; j < nClasses; ++j) {
                if (i == j) { continue;}
                float totDiff = (float)cieC.calcDeltaECIE2000(
                        meanLAB[i][0][0], meanLAB[i][1][0], meanLAB[i][2][0], 
                        meanLAB[j][0][0], meanLAB[j][1][0], meanLAB[j][2][0]);
                stats.totDiffs[count] = totDiff;
                stats.signifOfDiffs[count] = totDiff / stdvDeltaE;
                count++;
            }
            assert(count == nComb);
            stats.calcSignficanceStats();
        }
       
        return allStats;
    }
    
    /**
     * 
     * @param meanClassColors [class idx][color idx][avg, std dev]
     * @return 
     */
    public AllClassInterStats calculateDiffBetweenClasses(
        float[][][] meanClassColors) {
        
        AllClassInterStats allStats = new AllClassInterStats();
        
        int nClasses = meanClassColors.length;
        int nBands = meanClassColors[0].length;
        
        int nComb = nClasses - 1;
        
        for (int i = 0; i < nClasses; ++i) {
            
            OneClassInterStats stats = new OneClassInterStats();
            stats.totDiffs = new float[nComb];
            stats.signifOfDiffs = new float[nComb];
            allStats.stats.add(stats);
            
            int count = 0;
            for (int j = 0; j < nClasses; ++j) {
                if (i == j) { continue;}
                float totDiff = 0;
                float totStdv = 0;
                for (int k = 0; k < nBands; ++k) {
                    totDiff += Math.abs(meanClassColors[i][k][0] -
                        meanClassColors[j][k][0]);
                    totStdv += (meanClassColors[i][k][1] * 
                        meanClassColors[i][k][1]);
                }
                stats.totDiffs[count] = totDiff;
                stats.signifOfDiffs[count] = totDiff /
                    (float)Math.sqrt(totStdv);
                count++;
            }
            assert(count == nComb);
            stats.calcSignficanceStats();
        }
       
        return allStats;
    }
    
    public static class AllClassInterStats {
        List<OneClassInterStats> stats = new ArrayList<OneClassInterStats>();
    
        public float findMinSignficance() {
            float min = Float.MAX_VALUE;
            for (OneClassInterStats s : stats) {
                if (s.minSignificance < min) {
                    min = s.minSignificance;
                }
            }
            return min;
        }
        
        public float findMinOfAverageSignficance() {
            float min = Float.MAX_VALUE;
            for (OneClassInterStats s : stats) {
                if (s.avgSignificance < min) {
                    min = s.avgSignificance;
                }
            }
            return min;
        }
    }
    
    public static class OneClassInterStats {
        float[] totDiffs;
        float[] signifOfDiffs;
        float minSignificance;
        float avgSignificance;

        void calcSignficanceStats() {
            
            float min = Float.MAX_VALUE;
            double sum = 0;
            for (float s : signifOfDiffs) {
                if (s < min) {
                    min = s;
                }
                sum += s;
            }
            sum /= (double)signifOfDiffs.length;
            
            this.avgSignificance = (float)sum;
            this.minSignificance = min;
        }
    }
    
}
