package algorithms;

import algorithms.matrix.MatrixUtil;
import gnu.trove.list.TDoubleList;
import java.util.Arrays;
import java.util.Map;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVectorSub;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.sparse.ArpackSym;
import no.uib.cipr.matrix.sparse.FlexCompRowMatrix;

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
    public void testEigenvaluesEigenExample() throws NotConvergedException {
        // based on example at top-ish of http://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver.html#a4140972e2b45343d1ef1793c2824159c
		
        double[][] a = getData0();
        
        double[][] expectedEigValues = getExpectedValues0();
        
        String expectedEigVec = getExpectedVextors0();
        
        /*EVD is LAPACK dgeev
        
        DGEEV computes for an N-by-N real nonsymmetric matrix A, the
         eigenvalues and, optionally, the left and/or right eigenvectors.

         The right eigenvector v(j) of A satisfies
                          A * v(j) = lambda(j) * v(j)
         where lambda(j) is its eigenvalue.
         
        The left eigenvector u(j) of A satisfies
                       u(j)**H * A = lambda(j) * u(j)**H
         where u(j)**H denotes the conjugate-transpose of u(j).

         The computed eigenvectors are normalized to have Euclidean norm
         equal to 1 and largest component real.

        */
        
        DenseMatrix m = new DenseMatrix(a);
        EVD evd = no.uib.cipr.matrix.EVD.factorize(m);
        // the diagonal of D:
        double[] eigenValues = evd.getRealEigenvalues();
        double[] eigenValuesI = evd.getImaginaryEigenvalues();
        
        DenseMatrix leftEigenVectors = evd.getLeftEigenvectors();
        
        // is column 1 incorrect??
        // python bindings to LAPACK -evd and even -evr produces
        // column 1 same as column 0, which is what the expected test
        // values are.
        // Because the test for right eigen vector below fails for the
        //    first 2 columns of the resulting A * right eigen vector,
        //    that would suggest an error in the calculation of the 
        //    2nd eigen vector at least.
        // may need to add to MtJ, access to error bounds for the eigen values
        // (e.g. http://www.netlib.org/lapack/lug/node89.html)
        /*
          http://www.netlib.org/lapack/lug/node90.html
        ..."if [an eigenvalue] is close to other eigenvalues, 
            its corresponding eigenvector zi may be inaccurate."

        */
        
        DenseMatrix rightEigenVectors = evd.getRightEigenvectors();
        
        System.out.println("a=\n" + m.toString());
        System.out.println("left e=\n" + leftEigenVectors.toString());
        System.out.println("right e=\n" + rightEigenVectors.toString());
        
        System.out.println("leftT * M * right=\n" +
            MatrixUtil.multiply(leftEigenVectors.transpose(),
                MatrixUtil.multiply(m, rightEigenVectors)));
        
        System.out.println("expected values=\n" +
            new DenseMatrix(expectedEigValues).toString());
        System.out.println("result values=\n" + Arrays.toString(eigenValues));
        
        System.out.println("\nexpected vector=");
        System.out.println(expectedEigVec);
        //System.out.println("result vector=" + v.toString();
        
        DenseMatrix d = new DenseMatrix(eigenValues.length, eigenValues.length);
        for (int i = 0; i < eigenValues.length; ++i) {
            d.set(i, i, eigenValues[i]);
        }
        
        /*
        The right eigenvector v(j) of A satisfies
                          A * v(j) = lambda(j) * v(j)
         where lambda(j) is its eigenvalue.
        */
        System.out.println("a.norm=" + m.norm(Matrix.Norm.One));
        System.out.println("a.norm=" + m.norm(Matrix.Norm.Frobenius));
        // A * V
        DenseMatrix check0_right = MatrixUtil.multiply(m, rightEigenVectors);
        // D * V
        DenseMatrix check1_right = MatrixUtil.multiply(rightEigenVectors, d);
        
        // NOTE!!  The first 2 columns of check0_right and check1_right
        //   do not match...see close eigenvalue errors comments above...
                
        System.out.println("check0_right=\n" + check0_right);
        System.out.println("check1_right=\n" + check1_right);

        Object[] vecAndValues = algorithms.imageProcessing.util.MatrixUtil.eigenWithErrorFilter(m);
        assertNotNull(vecAndValues);
        DenseMatrix filteredVecs = (DenseMatrix) vecAndValues[0];
        TDoubleList filteredValues = (TDoubleList) vecAndValues[1];
        assertEquals(2, filteredValues.size());
    }

    /*
    adapting from:
       https://github.com/hughperkins/jeigen/blob/master/src/java/jeigen/TestEigenvalues.java
    which uses Mozilla public license
    http://mozilla.org/MPL/2.0/
    */
    // from http://lpsa.swarthmore.edu/MtrxVibe/EigMat/MatrixEigen.html
    public void testEigenvaluesLpsa() throws NotConvergedException {
        
        double[][] a = getData1();
        
        double[][] expectedValues = getExpectedValues1();
        
        double[][] expectedVectors = getExpectedVectors1();
        
        DenseMatrix m2 = new DenseMatrix(a);
        EVD evd2 = no.uib.cipr.matrix.EVD.factorize(m2);
        // the diagonal of D:
        double[] eigenValues2 = evd2.getRealEigenvalues();
        double[] eigenValuesI2 = evd2.getImaginaryEigenvalues();
        
        DenseMatrix leftEigenVectors2 = evd2.getLeftEigenvectors();
        DenseMatrix rightEigenVectors2 = evd2.getRightEigenvectors();

        System.out.println("\nexpected vector=" +
            (new DenseMatrix(expectedVectors)).toString());
        
        System.out.println("left e=\n" + leftEigenVectors2.toString());
        System.out.println("right e=\n" + rightEigenVectors2.toString());
        
        System.out.println("leftT * M * right=\n" +
            MatrixUtil.multiply(leftEigenVectors2.transpose(),
                MatrixUtil.multiply(m2, rightEigenVectors2)));
        
        System.out.println("expected values=" +
        (new DenseMatrix(expectedValues)).toString());
        
        System.out.println("LAPACK result values=\n" + 
            Arrays.toString(eigenValues2));
        
        DenseMatrix d2 = new DenseMatrix(eigenValues2.length, eigenValues2.length);
        for (int i = 0; i < eigenValues2.length; ++i) {
            d2.set(i, i, eigenValues2[i]);
        }
        
        /*
        The right eigenvector v(j) of A satisfies
                          A * v(j) = lambda(j) * v(j)
         where lambda(j) is its eigenvalue.
        */
        // A * V
        DenseMatrix check0_right = MatrixUtil.multiply(m2, rightEigenVectors2);
        // D * V
        DenseMatrix check1_right = MatrixUtil.multiply(rightEigenVectors2, d2);
        
        // NOTE:  The first 2 columns of check0_right and check1_right
        //   do not match.  possibly an error in the eigenvector
        //   calculations.
                
        System.out.println("check0_right=\n" + check0_right);
        System.out.println("check1_right=\n" + check1_right);

        try {
            check(check0_right, check1_right);
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
    public void testEig2() throws NotConvergedException {
        
        double[][] a = {{4.,1.,1.},{1.,2.,3.},{1.,3.,6.}};
        
        DenseMatrix m = new DenseMatrix(a);
        EVD evd = no.uib.cipr.matrix.EVD.factorize(m);
        // the diagonal of D:
        double[] eigenValues = evd.getRealEigenvalues();
        double[] eigenValuesI = evd.getImaginaryEigenvalues();
        
        DenseMatrix leftEigenVectors = evd.getLeftEigenvectors();
        
        DenseMatrix rightEigenVectors = evd.getRightEigenvectors();
        
        System.out.println("a=\n" + m.toString());
        System.out.println("left e=\n" + leftEigenVectors.toString());
        System.out.println("right e=\n" + rightEigenVectors.toString());
         
        System.out.println("leftT * M * right=\n" +
            MatrixUtil.multiply(leftEigenVectors.transpose(),
                MatrixUtil.multiply(m, rightEigenVectors)));
       
        DenseMatrix d = new DenseMatrix(eigenValues.length, eigenValues.length);
        for (int i = 0; i < eigenValues.length; ++i) {
            d.set(i, i, eigenValues[i]);
        }
        
        
        /*
        The right eigenvector v(j) of A satisfies
                          A * v(j) = lambda(j) * v(j)
         where lambda(j) is its eigenvalue.
        */
        // A * V
        DenseMatrix check0_right = MatrixUtil.multiply(m, rightEigenVectors);
        // D * V
        DenseMatrix check1_right = MatrixUtil.multiply(rightEigenVectors, d);
        
        try {
            check(check0_right, check1_right);
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
        
        FlexCompRowMatrix colM = new FlexCompRowMatrix(2, 2);
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
    private static void check(DenseMatrix X, DenseMatrix Y) {
    
        double eps = Math.pow(2.0, -52.0);
        
        if (X.norm(Matrix.Norm.Frobenius) == 0. && 
            Y.norm(Matrix.Norm.Frobenius) < 10 * eps) {
            return;
        }
        if (Y.norm(Matrix.Norm.Frobenius) == 0. && 
            X.norm(Matrix.Norm.Frobenius) < 10 * eps) {
            return;
        }
        DenseMatrix diff = MatrixUtil.subtract(X, Y);
        System.out.println("diff=\n" + diff.toString());
        if (diff.norm(Matrix.Norm.One) >
            1000 * eps * Math.max(X.norm(Matrix.Norm.One), 
                Y.norm(Matrix.Norm.Frobenius))) {
            throw new RuntimeException("The norm of (X-Y) is too large: " + 
                Double.toString(
                    MatrixUtil.subtract(X, Y).norm(Matrix.Norm.Frobenius)));
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
        String str 
            = " (-0.292,-0.454)   (-0.292,0.454)      (-0.0607,0)       (-0.733,0)    (0.59,-0.122)     (0.59,0.122)\n"
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
    
    private double[][] getExpectedRealVectors0() {
        double[][] a = new double[6][];
        a[0] = new double[]{-0.292, -0.292, -0.0607, -0.733, 0.59, 0.59};
        a[1] = new double[]{0.134, 0.134, -0.799, 0.136, 0.335, 0.335};
        a[2] = new double[]{-0.422, -0.422, 0.192, 0.0563, -0.335, -0.335};
        a[3] = new double[]{-0.589, -0.589, -0.0788, -0.627, 0.322, 0.322};
        a[4] = new double[]{-0.248, -0.248, 0.401, 0.218, -0.335, -0.335};
        a[5] = new double[]{0.105, 0.105, -0.392, -0.00564, -0.0324, -0.0324};
       
        return a;
    }

    private double[][] getExpectedImagVectors0() {
        double[][] a = new double[6][];
        a[0] = new double[]{-0.454,0.454,0,0,-0.122,0.122};
        a[1] = new double[]{-0.104,0.104,0,0,0.368,-0.368};
        a[2] = new double[]{-0.18,0.18,0,0,-0.143,0.143};
        a[3] = new double[]{0.0274,-0.0274,0,0,-0.156,0.156};
        a[4] = new double[]{ 0.132,-0.132,0,0,-0.076,0.076};
        a[5] = new double[]{0.18,-0.18,0,0,0.103,-0.103};
        
        return a;
    }

}
