package algorithms.misc;

/**
 * adapted from:
 
***********************************************************************
 *  Compilation:  javac Complex.java
 *  Execution:    java Complex
 *
 * Code is from http://introcs.cs.princeton.edu/java/stdlib/
 * "An Introduction to Programming in Java".
 * 
 * The code is released under the GNU General Public License, version 3 (GPLv3). 
 * If you wish to license the code under different terms, please contact our 
 * publisher to discuss.
 * 
 *  Data type for complex numbers.
 *
 *  The data type is "immutable" so once you create and initialize
 *  a Complex object, you cannot change it. The "final" keyword
 *  when declaring re and im enforces this rule, making it a
 *  compile-time error to change the .re or .im fields after
 *  they've been initialized.
 *
 *  % java Complex
 *  a            = 5.0 + 6.0i
 *  b            = -3.0 + 4.0i
 *  Re(a)        = 5.0
 *  Im(a)        = 6.0
 *  b + a        = 2.0 + 10.0i
 *  a - b        = 8.0 + 2.0i
 *  a * b        = -39.0 + 2.0i
 *  b * a        = -39.0 + 2.0i
 *  a / b        = 0.36 - 1.52i
 *  (a / b) * b  = 5.0 + 6.0i
 *  conj(a)      = 5.0 - 6.0i
 *  |a|          = 7.810249675906654
 *  tan(a)       = -6.685231390246571E-6 + 1.0000103108981198i
 *
 *************************************************************************/

public class ComplexModifiable {
    private double re;   // the real part
    private double im;   // the imaginary part

    // create a new object with the given real and imaginary parts
    public ComplexModifiable(double real, double imag) {
        re = real;
        im = imag;
    }
    
    public void setReal(double real) {
        re = real;
    }
    public void setImag(double imag) {
        im = imag;
    }
    
    public void resetTo(ComplexModifiable m) {
        this.re = m.re;
        this.im = m.im;
    }

    // return a string representation of the invoking Complex object
    public String toString() {
        if (im == 0) return re + "";
        if (re == 0) return im + "i";
        if (im <  0) return re + " - " + (-im) + "i";
        return re + " + " + im + "i";
    }

    // return abs/modulus/magnitude and angle/phase/argument
    public double abs()   { return Math.hypot(re, im); }  // Math.sqrt(re*re + im*im)
    public double phase() { return Math.atan2(im, re); }  // between -pi and pi

    public double squareSum() {
        return (re * re) + (im * im);
    }
    
    /**
     * add the value of b to this instance
     * @param b
     */
    public void plus(ComplexModifiable b) {
        re += b.re;
        im += b.im;
    }
    
    /**
     * add the value of b to this instance
     * @param b
     */
    public void plus(Complex b) {
        re += b.re();
        im += b.im();
    }

    /**
     * subtract the value of b from this instance
     * @param b
     */
    public void minus(ComplexModifiable b) {
        re -= b.re;
        im -= b.im;
    }
    
    /**
     * subtract the value of b from this instance
     * @param b
     */
    public void minus(Complex b) {
        re -= b.re();
        im -= b.im();
    }
    
    public ComplexModifiable copy() {
        return new ComplexModifiable(this.re, this.im);
    }

    /**
     * multiply the value of b with this instance
     * @param b
     */
    public void times(ComplexModifiable b) {
        double r = this.re;
        double i = this.im;
        re = r * b.re - i * b.im;
        im = r * b.im + i * b.re;
    }
    
    /**
     * multiply the value of b with this instance
     * @param b
     */
    public void times(Complex b) {
        double r = this.re;
        double i = this.im;
        re = r * b.re() - i * b.im();
        im = r * b.im() + i * b.re();
    }

    /**
     * multiply the value of b with this instance
     * @param b
     */
    public void times(double b) {
        re *= b;
        im *= b;
    }

    /**
     * turn this instance into the conjugate of itself
     */
    public void conjugate() {  
        im *= -1; 
    }
    
    /**
     * turn this instance into the reciprocal of itself
     */
    public void reciprocal() {  
        double scale = re*re + im*im;
        re /= scale;
        im /= (-scale);
    }

    // return the real or imaginary part
    public double re() { return re; }
    public double im() { return im; }

    /**
     * divide this instance by the value of b 
     * @param b
     */
    public void divided(ComplexModifiable b) {
        // reciprocal of b:
        double bR = b.re();
        double bI = b.im();
        double scaleB = bR*bR + bI*bI;
        bR /= scaleB;
        bI /= (-scaleB);
        
        re *= bR;
        im *= bI;
    }
    
    /**
     * divide this instance by the value of b 
     * @param b
     */
    public void divided(Complex b) {
        // reciprocal of b:
        double bR = b.re();
        double bI = b.im();
        double scaleB = bR*bR + bI*bI;
        bR /= scaleB;
        bI /= (-scaleB);
        
        re *= bR;
        im *= bI;
    }

    // return a new Complex object whose value is the complex exponential of this
    public ComplexModifiable exp() {
        return new ComplexModifiable(
            Math.exp(re) * Math.cos(im), 
            Math.exp(re) * Math.sin(im));
    }

    // return a new Complex object whose value is the complex sine of this
    public ComplexModifiable sin() {
        return new ComplexModifiable(Math.sin(re) * Math.cosh(im), Math.cos(re) * Math.sinh(im));
    }

    // return a new Complex object whose value is the complex cosine of this
    public ComplexModifiable cos() {
        return new ComplexModifiable(Math.cos(re) * Math.cosh(im), 
            -Math.sin(re) * Math.sinh(im));
    }

    // return a new Complex object whose value is the complex tangent of this
    public ComplexModifiable tan() {
        ComplexModifiable c = cos();
        ComplexModifiable s = sin();
        s.divided(c);
        return s;
    }
    
    // a static version of plus
    public static ComplexModifiable plus(ComplexModifiable a, ComplexModifiable b) {
        double real = a.re + b.re;
        double imag = a.im + b.im;
        ComplexModifiable sum = new ComplexModifiable(real, imag);
        return sum;
    }
}
