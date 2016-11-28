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
    
    //TODO: need a method for deltaE
    
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
