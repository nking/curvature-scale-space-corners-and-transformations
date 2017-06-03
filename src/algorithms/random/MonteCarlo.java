package algorithms.random;

import algorithms.misc.Misc;
import java.util.Random;

/**
 * random sampling of a system
 * 
 * @author nichole
 */
public class MonteCarlo {
    
    public double f(AFunction function, int nPoints) {
        
        double sum = 0;
        int nParams = function.getNumberOfCoeffs1();
        double[] c = new double[nParams];
        
        Random rng = Misc.getSecureRandom();
        
        for (int i = 0; i < nPoints; ++i) {
            for (int j = 0; j < nParams; ++j) {
                c[j] = rng.nextDouble();
            }
            sum += function.f(c);
        }
        
        sum /= (double)nPoints;
        
        return sum;
    }
    
}
