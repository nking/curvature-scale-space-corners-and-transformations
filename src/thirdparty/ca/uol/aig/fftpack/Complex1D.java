package thirdparty.ca.uol.aig.fftpack;
/**
  * Construct a 1-D complex data sequence.

    jfftpack is a Java version of fftpack. jfftpack is based
on Paul N. Swarztraubre's Fortran code and Pekka Janhuen's
C code. It is developed as part of my official duties as
lead software engineer for SCUBA-2 FTS projects
(www.roe.ac.uk/ukatc/projects/scubatwo/)

The original fftpack was public domain, so jfftpack is
public domain too. As such, this package is also released
under GPL. This software is in no way certified or guaranteed.

Notes:
Please read the following documents for FFT formula if necessary
http://www.netlib.org/fftpack/doc

*/
public class Complex1D
{
/**
  * <em>x</em>[<em>i</em>] is the real part of <em>i</em>-th complex data.
*/
    public double x[];
/**
  * <em>y</em>[<em>i</em>] is the imaginary part of <em>i</em>-th complex data.
*/
    public double y[];
}
