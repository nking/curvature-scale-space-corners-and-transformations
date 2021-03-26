package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import no.uib.cipr.matrix.NotConvergedException;

/**
 * various 3D alignment methods
 * 
 * @author nichole
 */
public class Alignment {
    
    /**
     * determine the rotation between measurements x1 and x2 when both datasets
     * have the same center, that is, there is no translation between them,
     * only rotation.
     * (see Golub & van Loan "Matrix Computations" 11.12.4,
     * Szeliski 2010, Sect 6.1.5).
     * @param x1 a set of measurements having same center as x2, that is, 
     * there is no translation between them, only rotation.
     * the expected format is nData X nDimensions.
     * @param x2 another set of measurements having same center as x1, that is, 
     * there is no translation between them, only rotation.
     * the expected format is nData X nDimensions.
     * @return 
     */
    public static double[][] procrustesAlgorithmForRotation(double[][] x1,
        double[][] x2) throws NotConvergedException {
        
        int m = x1.length;
        int p = x1[0].length;
        
        if (x2.length != m || x2[0].length != p) {
            throw new IllegalArgumentException("x1 and x2 must have same sizes");
        }
        
        // minimize || x1 - x2*Q ||_F
        //    subject to Q^T * Q = I_P
        double[][] c = MatrixUtil.multiply(MatrixUtil.transpose(x2), x1);
        
        MatrixUtil.SVDProducts svdC = MatrixUtil.performSVD(c);

        double[][] q = MatrixUtil.multiply(svdC.u, svdC.vT);
        
        return q;
    }
}
