package algorithms.imageProcessing.features;

/**
 *
 * @author nichole
 */
public class HOGUtil {
    
    private static float eps = 0.000001f;
    
    /**
     * CAVEAT: small amount of testing done, not yet throughly tested.
     *
     * calculate the intersection of histA and histB which have already
     * been normalized to the same scale.
     * A result of 0 is maximally dissimilar and a result of 1 is maximally similar.
     *
     * The orientations are needed to compare the correct rotated bins to one another.
     * Internally, orientation of 90 leads to no shift for rotation,
     * and orientation near 0 results in rotation of nBins/2, etc...
     *
     * Note that an orientation of 90 is a unit vector from x,y=0,0 to
     * x,y=0,1.
     * 
     * Note also that you may want to try the rotation of oppossite direction.
     *
     * @param histA
     * @param orientationA
     * @param histB
     * @param orientationB
     * @return
     */
    public static float intersection(int[] histA, int orientationA, int[] histB,
        int orientationB) {

        if ((histA.length != histB.length)) {
            throw new IllegalArgumentException(
                "histA and histB must be same dimensions");
        }

        if (orientationA < 0 || orientationA > 180 || orientationB < 0 ||
            orientationB > 180) {
            throw new IllegalArgumentException("orientations must be in range 0 to 180,"
                + "  inclusive,  or!=" + orientationA + " orB=" + orientationB);
        }
        if (orientationA == 180) {
            orientationA = 0;
        }
        if (orientationB == 180) {
            orientationB = 0;
        }

        int nBins = histA.length;

        int binWidth = 180/nBins;

        int shiftA = (orientationA - 90)/binWidth;
        int shiftB = (orientationB - 90)/binWidth;

        /*
        histograms are already normalized

        K(a,b) =
            (summation_over_i_from_1_to_n( min(a_i, b_i))
             /
            (min(summation_over_i(a_i), summation_over_i(b_i))
        */

        float sum = 0;
        float sumA = 0;
        float sumB = 0;
        for (int j = 0; j < nBins; ++j) {

            int idxA = j + shiftA;
            if (idxA < 0) {
                idxA += nBins;
            } else if (idxA > (nBins - 1 )) {
                idxA -= nBins;
            }

            int idxB = j + shiftB;
            if (idxB < 0) {
                idxB += nBins;
            } else if (idxB > (nBins - 1 )) {
                idxB -= nBins;
            }
            
            float yA = histA[idxA];
            float yB = histB[idxB];

            sum += Math.min(yA, yB);
            sumA += yA;
            sumB += yB;

            //System.out.println(" " + yA + " -- " + yB + " sum="+sum + ", " + sumA + "," + sumB);
        }

        float d = eps + Math.min(sumA, sumB);
        float sim = sum/d;

        return sim;
    }
    
    /**
     * CAVEAT: small amount of testing done, not yet throughly tested.
     * 
     * calculate the intersection of histA and histB which have already
     * been normalized to the same scale.
     * A result of 0 is maximally dissimilar and a result of 1 is maximally similar.
     * 
     The orientations are needed to compare the correct rotated bins to one another.
     * Internally, orientation of 90 leads to no shift for rotation,
     * and orientation near 0 results in rotation of nBins/2, etc...
     *
     * Note that an orientation of 90 is a unit vector from x,y=0,0 to
     * x,y=0,1.
     *
     * @param histA
     * @param orientationA
     * @param histB
     * @param orientationB
     * @return difference, error
     */
    public static float[] diff(int[] histA, int orientationA, int[] histB,
        int orientationB) {

        if ((histA.length != histB.length)) {
            throw new IllegalArgumentException(
                "histA and histB must be same dimensions");
        }

        if (orientationA < 0 || orientationA > 180 || orientationB < 0 ||
            orientationB > 180) {
            throw new IllegalArgumentException("orientations must be in range 0 to 180,"
                + "  inclusive,  or!=" + orientationA + " orB=" + orientationB);
        }
        if (orientationA == 180) {
            orientationA = 0;
        }
        if (orientationB == 180) {
            orientationB = 0;
        }

        int nBins = histA.length;

        int binWidth = 180/nBins;

        int shiftA = (orientationA - 90)/binWidth;
        int shiftB = (orientationB - 90)/binWidth;

        /*
        histograms are already normalized

        K(a,b) =
            (summation_over_i_from_1_to_n( min(a_i, b_i))
             /
            (min(summation_over_i(a_i), summation_over_i(b_i))
        */

        double sumDiff = 0;
        
        double err = 0;
        for (int j = 0; j < nBins; ++j) {

            int idxA = j + shiftA;
            if (idxA < 0) {
                idxA += nBins;
            } else if (idxA > (nBins - 1 )) {
                idxA -= nBins;
            }

            int idxB = j + shiftB;
            if (idxB < 0) {
                idxB += nBins;
            } else if (idxB > (nBins - 1 )) {
                idxB -= nBins;
            }

            float yA = histA[idxA];
            float yB = histB[idxB];

            float maxValue = Math.max(yA, yB) + eps;

            float diff = Math.abs((yA - yB)/maxValue);
            
            //sumDiff += (diff * diff);
            sumDiff += diff;

            //      already squared
            err += (diff/maxValue);
        }
        
        sumDiff /= (double)nBins;
        
        //sumDiff = Math.sqrt(sumDiff);
        
        err /= (double)nBins;
        err = Math.sqrt(err);
        
        return new float[]{(float)sumDiff, (float)err};
    }
}
