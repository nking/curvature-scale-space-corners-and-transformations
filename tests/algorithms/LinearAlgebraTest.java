package algorithms;

import java.util.Map;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseVectorSub;
import no.uib.cipr.matrix.sparse.ArpackSym;
import no.uib.cipr.matrix.sparse.FlexCompColMatrix;
import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.simple.SimpleEVD;
import org.ejml.simple.SimpleMatrix;

/**
 * a general location for the matrix linear algebra tests.
 * 
 * for the eigen solver tests,
    adapting from:
       https://github.com/hughperkins/jeigen/blob/master/src/java/jeigen/TestEigenvalues.java
    which uses Mozilla public license
    http://mozilla.org/MPL/2.0/
     
    one test is from NIST's JAMA and is referenced below too.
 */
public class LinearAlgebraTest extends TestCase {
    
    /*
    for the eigen solver tests,
    adapting from:
       https://github.com/hughperkins/jeigen/blob/master/src/java/jeigen/TestEigenvalues.java
    which uses Mozilla public license
    http://mozilla.org/MPL/2.0/
    */    
    public void testEigenvaluesEigenExample() {
        // based on example at top-ish of http://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver.html#a4140972e2b45343d1ef1793c2824159c
		
        double[][] a = getData0();
        
        double[][] expectedEigValues = getExpectedValues0();
        
        String expectedEigVec = getExpectedVextors0();
        
        SimpleMatrix m = new SimpleMatrix(a);
        SimpleEVD evd = m.eig();
        
        EigenDecomposition evdD = evd.getEVD();
        SimpleMatrix d = new SimpleMatrix(org.ejml.ops.EigenOps.createMatrixD(evdD));
        SimpleMatrix v = new SimpleMatrix(org.ejml.ops.EigenOps.createMatrixV(evdD));
        
        System.out.println("expected values=");
        (new SimpleMatrix(expectedEigValues)).print();;
        System.out.println("result values=");
        d.print();
        
        System.out.println("\nexpected vector=");
        System.out.println(expectedEigVec);
        System.out.println("result vector=");
        v.print();
        
        // A * V
        SimpleMatrix check0 = m.mult(v);
        // V * D
        SimpleMatrix check1 = v.mult(d);
        
        try {
            check(check0, check1);
        } catch (java.lang.RuntimeException e) {
            fail(e.getMessage());
        }
    }

    /*
    adapting from:
       https://github.com/hughperkins/jeigen/blob/master/src/java/jeigen/TestEigenvalues.java
    which uses Mozilla public license
    http://mozilla.org/MPL/2.0/
    */
    // from http://lpsa.swarthmore.edu/MtrxVibe/EigMat/MatrixEigen.html
    public void testEigenvaluesLpsa() {
        
        double[][] a = getData1();
        
        double[][] expectedValues = getExpectedValues1();
        
        double[][] expectedVectors = getExpectedVectors1();
        
        SimpleMatrix m = new SimpleMatrix(a);
        SimpleEVD evd = m.eig();
        
        EigenDecomposition evdD = evd.getEVD();
        SimpleMatrix d = new SimpleMatrix(org.ejml.ops.EigenOps.createMatrixD(evdD));
        SimpleMatrix v = new SimpleMatrix(org.ejml.ops.EigenOps.createMatrixV(evdD));
        
        System.out.println("expected values=");
        (new SimpleMatrix(expectedValues)).print();
        System.out.println("result values=");
        d.print();
        
        System.out.println("\nexpected vector=");
        (new SimpleMatrix(expectedVectors)).print();;
        System.out.println("result vector=");
        v.print();
        
        // A * V
        SimpleMatrix check0 = m.mult(v);
        // V * D
        SimpleMatrix check1 = v.mult(d);
        
        try {
            check(check0, check1);
        } catch (java.lang.RuntimeException e) {
            fail(e.getMessage());
        }
    }
    
    /*
    adapting from:
       https://github.com/hughperkins/jeigen/blob/master/src/java/jeigen/TestEigenvalues.java
    which uses Mozilla public license
    http://mozilla.org/MPL/2.0/
    */
    public void testEigenvaluesCircul() { // from http://www.mathworks.com/help/matlab/ref/eig.html
		
        double[][] a = getData2();
        
        // results differ from matlabs results.  Just another equivalent reuslt? or a bug?
    }
   
    /*
    adapting from:
       https://github.com/hughperkins/jeigen/blob/master/src/java/jeigen/TestEigenvalues.java
    which uses Mozilla public license
    http://mozilla.org/MPL/2.0/
    */
    public void testPseudoEigenvaluesComplexCase() {
        // based on example at bottom of http://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver.html#a4140972e2b45343d1ef1793c2824159c
		double[][] a = getData3();

        double[][] expectedValues = getExpectedValues3();
        
        double[][] expectedVectors = getExpectedVectors3();
    }
    
    // this is from JAMA matrix test http://math.nist.gov/javanumerics/jama/
    // Jama-1.0.3.zip
    public void testEig2() {
        
        double[][] a = {{4.,1.,1.},{1.,2.,3.},{1.,3.,6.}};
        
        SimpleMatrix m = new SimpleMatrix(a);
        SimpleEVD evd = m.eig();
        
        EigenDecomposition evdD = evd.getEVD();
        SimpleMatrix d = new SimpleMatrix(org.ejml.ops.EigenOps.createMatrixD(evdD));
        SimpleMatrix v = new SimpleMatrix(org.ejml.ops.EigenOps.createMatrixV(evdD));
        
        System.out.println("result values=");
        d.print();
        
        System.out.println("\nresult vectors=");
        v.print();
        
        // A * V
        SimpleMatrix check0 = m.mult(v);
        // V * D
        SimpleMatrix check1 = v.mult(d);
        
        try {
            check(check0, check1);
        } catch (java.lang.RuntimeException e) {
            fail(e.getMessage());
        }
    }
    
    public void testEig1_MTJ() {
        
        // javadocs are at http://www.javadoc.io/doc/com.googlecode.matrix-toolkits-java/mtj/1.0.4
        
        /*
        ArpackSym 
           ARPACK is designed to compute a subset of eigenvalues/eigenvectors
               ArpackSym(Matrix matrix) 
               Map<Double,DenseVectorSub>	solve(int eigenvalues, ArpackSym.Ritz ritz)
                   Solve the eigensystem for the number of eigenvalues requested.
            The two relevant Ritz enums for normalized cuts are
                SA - compute the NEV smallest (algebraic) eigenvalues.
                SM - compute the NEV smallest (in magnitude) eigenvalues.
                     *this is what scipy is using 'SM'
        
        CompColMatrix uses compressed column storage format
        
        CompDiagMatrix uses Compressed diagonal storage (CDS) matrix
        
        FlexCompRowMatrix is a Matrix stored row-wise into sparse vectors
        
        */
                
        // results in faster solution, but are lookup times increased?
        //LinkedSparseMatrix a = new LinkedSparseMatrix(data.length, data[0].length);
        
        FlexCompColMatrix colM = new FlexCompColMatrix(2, 2);
        colM.add(0, 0, 2);
        colM.add(1, 1, 2);
        
        ArpackSym arpackSym = new ArpackSym(colM);
        Map<Double, DenseVectorSub> rMap = arpackSym.solve(1, ArpackSym.Ritz.SM);
        
        for (Map.Entry<Double, DenseVectorSub> result : rMap.entrySet()) {           
            System.out.println("resulting eigenvalue=" + result.getKey().toString());
            System.out.println("resulting eigenvector=" + result.getValue().toString());
        }        
    }
    
    // this is from JAMA matrix test http://math.nist.gov/javanumerics/jama/
    // Jama-1.0.3.zip
    private static void check(SimpleMatrix X, SimpleMatrix Y) {
        double eps = Math.pow(2.0, -52.0);
        if (X.normF() == 0. & Y.normF() < 10 * eps) {
            return;
        }
        if (Y.normF() == 0. & X.normF() < 10 * eps) {
            return;
        }
        if (X.minus(Y).normF() > 1000 * eps * Math.max(X.normF(), Y.normF())) {
            throw new RuntimeException("The norm of (X-Y) is too large: " + 
                Double.toString(X.minus(Y).normF()));
        }
    }
    
    private double[][] getData0() {
        double[][] a = new double[6][];
        a[0] = new double[]{0.68, -0.33, -0.27, -0.717, -0.687, 0.0259};
        a[1] = new double[]{-0.211, 0.536, 0.0268, 0.214, -0.198, 0.678};
        a[2] = new double[]{0.566, -0.444, 0.904, -0.967, -0.74, 0.225};
        a[3] = new double[]{0.597, 0.108, 0.832, -0.514, -0.782, -0.408};
        a[4] = new double[]{0.823, -0.0452, 0.271, -0.726, 0.998, 0.275};
        a[5] = new double[]{-0.605, 0.258, 0.435, 0.608, -0.563, 0.0486};
        return a;
    }
    
    private double[][] getExpectedValues0() {
        double[][] expectedEigValues = new double[6][];
        expectedEigValues[0] = new double[]{0.049, 1.06};
        expectedEigValues[1] = new double[]{0.049, -1.06};
        expectedEigValues[2] = new double[]{0.967, 0};
        expectedEigValues[3] = new double[]{0.353, 0};
        expectedEigValues[4] = new double[]{0.618, 0.129};
        expectedEigValues[5] = new double[]{0.618, -0.129};
        return expectedEigValues;
    }
    
    private String getExpectedVextors0() {
        String str = " (-0.292,-0.454)   (-0.292,0.454)      (-0.0607,0)       (-0.733,0)    (0.59,-0.122)     (0.59,0.122)\n"
            + "  (0.134,-0.104)    (0.134,0.104)       (-0.799,0)        (0.136,0)    (0.335,0.368)   (0.335,-0.368)\n"
            + "  (-0.422,-0.18)    (-0.422,0.18)        (0.192,0)       (0.0563,0)  (-0.335,-0.143)   (-0.335,0.143)\n"
            + " (-0.589,0.0274) (-0.589,-0.0274)      (-0.0788,0)       (-0.627,0)   (0.322,-0.156)    (0.322,0.156)\n"
            + "  (-0.248,0.132)  (-0.248,-0.132)        (0.401,0)        (0.218,0)  (-0.335,-0.076)   (-0.335,0.076)\n"
            + "    (0.105,0.18)    (0.105,-0.18)       (-0.392,0)     (-0.00564,0)  (-0.0324,0.103) (-0.0324,-0.103\n";
        
        return str;
    }

    private double[][] getData1() {
        double[][] a = new double[2][];
        a[0] = new double[]{0, 1};
        a[1] = new double[]{-2, -3};
        return a;
    }
    
    private double[][] getExpectedValues1() {
        double[][] expectedValues = new double[2][];
        expectedValues[0] = new double[]{-1, 0};
        expectedValues[1] = new double[]{0, -2};
        return expectedValues;
    }
    
    private double[][] getExpectedVectors1() {
        double[][] expectedVectors = new double[2][];
        expectedVectors[0] = new double[]{0.7071, -0.4472};
        expectedVectors[1] = new double[]{-0.7071, 0.8944};
        return expectedVectors;
    }
    
    public double[][] getData2() {
        double[][] a = new double[3][];
        a[0] = new double[]{1, 2, 3};
        a[1] = new double[]{1, 1, 2};
        a[2] = new double[]{2, 3, 1};
        return a;
    }
    
    public double[][] getData3() {
        double[][] a = new double[][]{
            {0.68, -0.33, -0.27, -0.717, -0.687, 0.0259},
            {-0.211, 0.536, 0.0268, 0.214, -0.198, 0.678},
            {0.566, -0.444, 0.904, -0.967, -0.74, 0.225},
            {0.597, 0.108, 0.832, -0.514, -0.782, -0.408},
            {0.823, -0.0452, 0.271, -0.726, 0.998, 0.275},
            {-0.605, 0.258, 0.435, 0.608, -0.563, 0.0486}};
        return a;
    }
    
    public double[][] getExpectedValues3() {
        double[][] expectedValues = new double[][]{
            {0.049, 1.06, 0, 0, 0, 0},
            {-1.06, 0.049, 0, 0, 0, 0},
            {0, 0, 0.967, 0, 0, 0},
            {0, 0, 0, 0.353, 0, 0},
            {0, 0, 0, 0, 0.618, 0.129},
            {0, 0, 0, 0, -0.129, 0.618}};
        return expectedValues;
    }
    
    public double[][] getExpectedVectors3() {
        double[][] expectedVectors = new double[][]{
            {-0.571, -0.888, -0.066, -1.13, 17.2, -3.54},
            {0.263, -0.204, -0.869, 0.21, 9.73, 10.7},
            {-0.827, -0.352, 0.209, 0.0871, -9.75, -4.17},
            {-1.15, 0.0535, -0.0857, -0.971, 9.36, -4.53},
            {-0.485, 0.258, 0.436, 0.337, -9.74, -2.21},
            {0.206, 0.353, -0.426, -0.00873, -0.942, 2.98}};
        return expectedVectors;
    }
    
    
}
