package algorithms.imageProcessing;

import algorithms.misc.Complex;
import algorithms.misc.MiscMath;

public class FFT {
    
    private boolean performNormalization = true;
    
    public void setToNotNormalize() {
        performNormalization = false;
    }
    
    /**
     * compute the FFT of x[], assuming its length is a power of 2.
     * runtime complexity is O(N).
     * (adapted from  Cormen et al. pseudocode)
     * ability added
     */ 
    public Complex[] fft(Complex[] x) {
        
        return fft(x, true);
    }

    /**
     * FFT from  Cormen et al. pseudocode for iterative FFT w/ an inverse
     * ability added.  Note that the length of x must be a power of 2.
     * runtime complexity is O(N).
     * 
     * @param x
     * @param forward run the transform in forward if true, else perform inverse 
     * transform.
     * @return 
     */
    protected Complex[] fft(Complex[] x, boolean forward) {
       
        if (x == null || x.length == 0) {
            throw new IllegalArgumentException("xReal cannot be null or empty");
        }
     
        int n = x.length;

        if (n == 1) {
            return x;
        }
        
        if (!MiscMath.isAPowerOf2(n)) {
            throw new IllegalArgumentException("x's length has to be a power of 2");
        }
        
        Complex[] a = bitReverseCopy(x);
        
        //TODO: could improve the speed at expense of space by caching
        // norm, end, m, eCoeff, and wn
        
        double norm = 1./Math.sqrt(n);
        
        int end = (int)(Math.log(n)/Math.log(2));
        
        for (int s = 1; s <= end; s++) {
            
            int m = 1 << s;
            
            double eCoeff = 2. * Math.PI/(double)m;
            
            Complex wn = forward ? 
                new Complex(Math.cos(eCoeff), -Math.sin(eCoeff)) :
                new Complex(Math.cos(eCoeff), Math.sin(eCoeff));
            
            for (int k = 0; k < n; k+=m) {
                
                Complex w = new Complex(1, 0);
                
                for (int j = 0; j < (m/2); j++) {
                    
                    Complex t = w.times(a[k + j + (m/2)]);
                    
                    Complex u = a[k + j];
                    a[k + j + (m/2)] = u.minus(t);
                    a[k + j] = u.plus(t);
                    
                    w = w.times(wn);
                }
            }
        }
        
        if (performNormalization) {
            for (int i = 0; i < a.length; i++) {
                a[i] = a[i].times(norm);
            }
        }
        
        /*
        bit-reverse-copy(a,A)
        for s=1 to lg n {
            m = 2^s
            wm = exp^(i*2*PI/m)
            for k=0 to n-1 by m {
                w=1
                for j=0 to ((m/2)-1) {
                    do t=w*A[k + j + (m/2)]
                    u = A[k + j]
                    A[k + j] = u + t
                    A[k + j + (m/2)] = u - t
                    w = w * wm
                }
            }
        }
        */
        
        return a;
    }
    
    /**
     * compute the inverse FFT of x[], assuming its length is a power of 2
     * runtime complexity is O(N).
     * (adapted from  Cormen et al. pseudocode for forward transform)
     */
    public Complex[] ifft(Complex[] x) {
        
       return fft(x, false);
        
    }

    protected double[] bitReverseCopy(double[] x) {
        
        int n = x.length;
        
        int nBits = MiscMath.numberOfBits(n - 1);
                        
        double[] r = new double[n];
        
        for (int k = 0; k < n; k++) {
            
            int idx = MiscMath.bitReverse(k, nBits);
            
            r[idx] = x[k];
        }
        
        return r;
    }
    
    protected double[] bitReverseCopy(int[] x) {
        
        int n = x.length;
        
        int nBits = MiscMath.numberOfBits(n - 1);
                        
        double[] r = new double[n];
        
        for (int k = 0; k < n; k++) {
            
            int idx = MiscMath.bitReverse(k, nBits);
            
            r[idx] = x[k];
        }
        
        return r;
    }
    
    protected Complex[] bitReverseCopy(Complex[] x) {
        
        int n = x.length;
        
        int nBits = MiscMath.numberOfBits(n - 1);
                        
        Complex[] r = new Complex[n];
        
        for (int k = 0; k < n; k++) {
            
            int idx = MiscMath.bitReverse(k, nBits);
            
            r[idx] = x[k].copy();
        }
        
        return r;
    }
    
    /**
     * perform FFT on x
     * runtime complexity is O(N).
     * (adapted from  Cormen et al. pseudocode for forward transform).
     * Note that the result cannot be inverted because only the real portion is 
     * returned and inverse FFT needs the real and complex components.
     */
    public double[] fft(int[] x) {

        if (x == null || x.length == 0) {
            throw new IllegalArgumentException("xReal cannot be null or empty");
        }
     
        int n = x.length;

        if (n == 1) {
            return new double[]{x[0]};
        }
        
        if (!MiscMath.isAPowerOf2(n)) {
            throw new IllegalArgumentException("xReal's length has to be a power of 2");
        }
        
        double[] a = bitReverseCopy(x);
        
        double norm = 1./Math.sqrt(n);
        
        int end = (int)(Math.log(n)/Math.log(2));
        
        for (int s = 1; s <= end; s++) {
            
            int m = 1 << s;
            
            double eCoeff = 2. * Math.PI/(double)m;
            double wnReal = Math.cos(eCoeff);
            double wnImag = Math.sin(eCoeff);
        
            for (int k = 0; k < n; k+=m) {
                
                double wReal = 1;
                double wImag = 0;
                
                for (int j = 0; j < (m/2); j++) {
                    
                    //complex multiplication:
                    double tReal = wReal * a[k + j + (m/2)];
                    double tImag = wImag * a[k + j + (m/2)];
                    double tAbs = Math.hypot(tReal, tImag);;
                    
                    double u = a[k + j];
                    a[k + j] = (u + tAbs);
                    a[k + j + (m/2)] = (u - tAbs);
                    
                    //complex multiplication:
                    wReal = wReal * wnReal - (wImag * wnImag);
                    wImag = wReal * wnImag + (wImag * wnReal);
                }
            }
        }
        
        if (performNormalization) {
            for (int i = 0; i < a.length; i++) {
                a[i] *= norm;
            }
        }
        
        /*
        bit-reverse-copy(a,A)
        for s=1 to lg n {
            m = 2^s
            wm = exp^(i*2*PI/m)
            for k=0 to n-1 by m {
                w=1
                for j=0 to ((m/2)-1) {
                    do t=w*A[k + j + (m/2)]
                    u = A[k + j]
                    A[k + j] = u + t
                    A[k + j + (m/2)] = u - t
                    w = w * wm
                }
            }
        }
        */
        
        return a;
    }
}
