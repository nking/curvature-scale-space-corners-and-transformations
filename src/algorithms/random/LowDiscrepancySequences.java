package algorithms.random;

/**
 * LDS are quasi-random or sub-random sequences, 
 * and are commonly used as a replacement of uniformly 
 * distributed random numbers.
 * 
 * https://en.wikipedia.org/wiki/Low-discrepancy_sequence
 * 
 * D_N(x_1, ..., x_N) .lte. C*(((ln N)^s)/N)
 *    where C is a const dep upon seq
 * 
 * Common LDS algorithms are
 *  van der Corput, Halton, Sobol, and Tezuka
 * 
 * @author nichole
 */
public class LowDiscrepancySequences {
    
    /**
     * caveat is correlated points for some specific
     * primes.  TODO: add known corrections for bases
     * which need them.
     * @param index
     * @param outputPair 
     */
    public void haltonPoint(int index, double[] outputPair) {
        
        // base 2 for x, base 3 for y
        outputPair[0] = halton(index, 2);
        outputPair[1] = halton(index, 3);
    }
    
    /**
     * caveat is correlated points for some specific
     * primes.  TODO: add known corrections for bases
     * which need them.
     * @param index
     * @param base
     * @return 
     */
    public double halton(int index, int base) {
        
        double f = 1;
        double r = 0;
        
        while (index > 0) {
            f /= (double)base;
            r += (f * (index % base));
            index = (index/base);
        }
        
        return r;
    }
    
}
