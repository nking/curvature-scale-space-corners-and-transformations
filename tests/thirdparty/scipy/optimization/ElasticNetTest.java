package thirdparty.scipy.optimization;

import algorithms.misc.Misc;
import gnu.trove.list.array.TIntArrayList;
import java.util.Arrays;
import java.util.Random;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

/**
 *
 * @author nichole
 */
public class ElasticNetTest extends TestCase {
    
    public ElasticNetTest() {
    }
    
    boolean debug = false;
    
    /**
     * 
     * @throws NotConvergedException 
     */
    public void testPoly() throws NotConvergedException {
        
        //test data from: http://rosettacode.org/wiki/Polynomial_Fitting
        double[] x = new double[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        double[] y = new double[]{1, 6, 17, 34, 57, 86, 121, 162, 209, 262, 321};
       
        double[][] x2 = new double[x.length][];
        for (int i = 0; i < x.length; i++) {
            x2[i] = new double[]{x[i], x[i] * x[i]};
        }
        
        double alpha = 0.01;
        double l1_ratio = 0.7;
        ElasticNet en = new ElasticNet(alpha, l1_ratio);
        en.fit(x2, y);
        
        double[] coef = en.getCoef();
        
        if (debug) {
        System.out.println("coef=" + Arrays.toString(coef));
        }
        
        double[] expected = new double[]{2., 3.};
        assertEquals(expected.length, coef.length);
        for (int i = 0; i < coef.length; ++i) {
            assertTrue(Math.abs(coef[i] - expected[i]) < 0.05);
        }
    }
    
    /**
     test from scipy
     https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/linear_model/tests/test_coordinate_descent.py
       
     scikit-learn/scikit-learn is licensed under the
        BSD 3-clause "New" or "Revised" License

     A permissive license similar to the BSD 2-Clause License, but with a 
       3rd clause that prohibits others from using the name of the project 
       or its contributors to promote derived products without written consent.
    */
    public void testENetToy() {
        
        /*
        Test ElasticNet for various parameters of alpha and l1_ratio.
        Actually, the parameters alpha = 0 should not be allowed. However,
        we test it as a border case.
        ElasticNet is tested with and without precomputed Gram matrix
        */
        double[][] x = new double[3][1];
        x[0] = new double[]{-1};
        x[1] = new double[]{0};
        x[2] = new double[]{1};
        
        double[] y = new double[]{-1, 0, 1};
        
        // test sample
        double[][] t = new double[3][1];
        t[0] = new double[]{2};
        t[1] = new double[]{3};
        t[2] = new double[]{4};
        
        // this should be the same as lasso
        double alpha = 1e-8;
        double l1_ratio = 1.0;
        ElasticNet clf = new ElasticNet(alpha, l1_ratio);
        clf.fit(x, y);
        double[] pred = clf.predict(t);
        double[] coef = clf.getCoef();
        assertEquals(1, coef.length);
        assertEquals(3, pred.length);
        double[] expected = new double[]{1};
        
        if (debug) {
        System.out.println("0 coef=" + Arrays.toString(coef));
        System.out.println("0 pred=" + Arrays.toString(pred));
        System.out.println("0 expected pred=" + Arrays.toString(expected));
        System.out.println("0 dualGaps=" + clf.getDualGap());
        }
        
        for (int i = 0; i < coef.length; ++i) {
            assertTrue(Math.abs(coef[i] - 1.) < 0.1);
        }
        
        expected = new double[]{2, 3, 4};
        for (int i = 0; i < pred.length; ++i) {
            assertTrue(Math.abs(pred[i] - expected[i]) < 0.1);
        }
        //assert_almost_equal(clf.dual_gap_, 0)

        alpha = 0.5;
        l1_ratio = 0.3;
        int max_iter = 100;
        boolean precompute = false;
        clf = new ElasticNet(alpha, l1_ratio, max_iter);
        clf.fit(x, y);
        pred = clf.predict(t);
        expected = new double[]{0.50819};
        coef = clf.getCoef();
        
        if (debug) {
        System.out.println("1 coef=" + Arrays.toString(coef));
        System.out.println("1 pred=" + Arrays.toString(pred));
        System.out.println("1 expected pred=" + Arrays.toString(expected));
        System.out.println("1 dualGaps=" + clf.getDualGap());
        }
        
        assertEquals(expected.length, coef.length);
        for (int i = 0; i < coef.length; ++i) {
            assertTrue(Math.abs(coef[i] - expected[i]) < 0.001);
        }
        expected = new double[]{1.0163, 1.5245, 2.0327};
        for (int i = 0; i < pred.length; ++i) {
            assertTrue(Math.abs(pred[i] - expected[i]) < 0.001);
        }
        //assert_almost_equal(clf.dual_gap_, 0)

        alpha = 0.5;
        l1_ratio = 0.5;
        clf = new ElasticNet(alpha, l1_ratio);
        clf.fit(x, y);
        pred = clf.predict(t);
        expected = new double[]{0.45454};
        coef = clf.getCoef();
        
        if (debug) {
        System.out.println("2 coef=" + Arrays.toString(coef));
        System.out.println("2 pred=" + Arrays.toString(pred));
        System.out.println("2 expected pred=" + Arrays.toString(expected));
        System.out.println("2 dualGaps=" + clf.getDualGap());
        }
        
        assertEquals(expected.length, coef.length);
        for (int i = 0; i < coef.length; ++i) {
            assertTrue(Math.abs(coef[i] - expected[i]) < 0.001);
        }
        expected = new double[]{0.9090, 1.3636, 1.8181};
        for (int i = 0; i < pred.length; ++i) {
            assertTrue(Math.abs(pred[i] - expected[i]) < 0.001);
        }        
        //assert_almost_equal(clf.dual_gap_, 0)        
    }
    
    public void test_warm_start_convergence_overconstrined() {
        
        Data data = build_dataset();
        data = build_dataset(200, 180, 5);
        
        double[][] X = data.X;
        double[] y = data.y;
        
        ElasticNet model = new ElasticNet(1e-3, 1e-3);
        model.fit(X, y);
        
        if (debug) {
        System.out.println("2 coef=" + Arrays.toString(model.getCoef()));
        System.out.println("2 dualGaps=" + model.getDualGap());
        }
        
        TIntArrayList nIters = model._nIters();
        if (debug) {
        for (int i = 0; i < nIters.size(); ++i) {
            System.out.println("i=" + i + " nIter=" + nIters.get(i));
        }
        }
        int n0 = nIters.get(0);

        //# This dataset is not trivial enough for the model to 
        // converge in one pass.
        assertTrue(n0 > 2);

        //# Check that n_iter_ is invariant to multiple calls to fit
        //# when warm_start=False, all else being equal.
        model.fit(X, y);
        nIters = model._nIters();
        if (debug) {
        for (int i = 0; i < nIters.size(); ++i) {
            System.out.println("i=" + i + " nIter=" + nIters.get(i));
        }
        }
        int n1 = nIters.get(0);
        assertEquals(n0, n1);

        //# Fit the same model again, using a warm start: the optimizer 
        //just performs
        //# a single pass before checking that it has already converged
        model.setToUseWarmStart();
        model.fit(X, y);
        nIters = model._nIters();
        if (debug) {
        for (int i = 0; i < nIters.size(); ++i) {
            System.out.println("i=" + i + " nIter=" + nIters.get(i));
        }
        }
        int n2 = nIters.get(0);
        assertEquals(1, n2);
    }
    
    public void test_warm_start_convergence_underconstrined() {
        
        Data data = build_dataset();
        data = build_dataset(100, 200, 10);
        
        double[][] X = data.X;
        double[] y = data.y;
        
        ElasticNet model = new ElasticNet(1e-3, 1e-3);
        model.fit(X, y);
        
        if (debug) {
        System.out.println("2 coef=" + Arrays.toString(model.getCoef()));
        System.out.println("2 dualGaps=" + model.getDualGap());
        }
        
        TIntArrayList nIters = model._nIters();
        if (debug) {
        for (int i = 0; i < nIters.size(); ++i) {
            System.out.println("i=" + i + " nIter=" + nIters.get(i));
        }
        }
        int n0 = nIters.get(0);

        //# This dataset is not trivial enough for the model to 
        // converge in one pass.
        assertTrue(n0 > 2);

        //# Check that n_iter_ is invariant to multiple calls to fit
        //# when warm_start=False, all else being equal.
        model.fit(X, y);
        nIters = model._nIters();
        if (debug) {
        for (int i = 0; i < nIters.size(); ++i) {
            System.out.println("i=" + i + " nIter=" + nIters.get(i));
        }
        }
        int n1 = nIters.get(0);
        assertEquals(n0, n1);

        //# Fit the same model again, using a warm start: the optimizer 
        //just performs
        //# a single pass before checking that it has already converged
        model.setToUseWarmStart();
        model.fit(X, y);
        nIters = model._nIters();
        if (debug) {
        for (int i = 0; i < nIters.size(); ++i) {
            System.out.println("i=" + i + " nIter=" + nIters.get(i));
        }
        }
        int n2 = nIters.get(0);
        assertEquals(1, n2);
    }
    
    public void test_random_descent() {
        
        //# Test that both random and cyclic selection give the same results.
        //# Ensure that the test models fully converge and check a wide
        //# range of conditions.

        //# This uses the coordinate descent algo using the gram trick.
        Data data = build_dataset();
        data = build_dataset(50, 20, 20);
        
        double[][] X = data.X;
        double[] y = data.y;
        
        ElasticNet model = new ElasticNet(1e-3, 1e-8);
        model.setSeletion(ElasticNet.Selection.CYCLIC);
        model.fit(X, y);
        
        ElasticNet rand = new ElasticNet(1e-3, 1e-8);
        rand.setSeletion(ElasticNet.Selection.RANDOM);
        rand.fit(X, y);
        
        TIntArrayList nIters = model._nIters();
        if (debug) {
        for (int i = 0; i < nIters.size(); ++i) {
            System.out.println("i=" + i + " nIter=" + nIters.get(i));
        }
        }
        int n0 = nIters.get(0);
        
        double[] c0 = model.getCoef();
        double[] c1 = rand.getCoef();
        
        if (debug) {
        System.out.println("sequential coef=" + Arrays.toString(c0));
        System.out.println("random coef=" + Arrays.toString(c1));
        }
        
        assertTrue(ElasticNet.allClose(c0, c1));

        double inter0 = model.getIntercept();
        double inter1 = rand.getIntercept();
        assertTrue(Math.abs(inter0 - inter1) < 0.01);
    }
    
    private Data build_dataset() {
        int nSamples = 50;
        int nFeatures = 200;
        int nInformative = 10;
        
        return build_dataset(nSamples, nFeatures, nInformative);
    }
    
    private Data build_dataset(int nSamples, int nFeatures, 
        int nInformative) {
       
        //build an ill-posed linear regression problem with many noisy 
        //features and comparatively few samples
    
        int nTargets = 1;
        
        Random rng = Misc.getSecureRandom();
        long seed = System.currentTimeMillis();
        if (debug) {
        System.out.println("SEED=" + seed);
        }
        rng.setSeed(seed);
        
        double[] w = new double[nFeatures];
        for (int i = 0; i < nFeatures; ++i) {
            w[i] = rng.nextDouble();
        }
        for (int i = nInformative; i < nFeatures; ++i) {
            w[i] = 0.;
        }
        
        //X = random_state.randn(n_samples, n_features)
        double[][] X = new double[nSamples][nFeatures];
        for (int i = 0; i < nSamples; ++i) {
            X[i] = new double[nFeatures];
            for (int j = 0; j < nFeatures; ++j) {
                X[i][j] = rng.nextDouble();
            }
        }
        
        double[] y = new double[nSamples];
        for (int row = 0; row < nSamples; ++row) {
            double sum = 0;
            for (int col = 0; col < nFeatures; ++col) {
                sum += X[row][col] * w[col];
            }
            y[row] = sum;
        }
        
        //X_test = random_state.randn(n_samples, n_features)
        double[][] Xtest = new double[nSamples][nFeatures];
        for (int i = 0; i < nSamples; ++i) {
            Xtest[i] = new double[nFeatures];
            for (int j = 0; j < nFeatures; ++j) {
                Xtest[i][j] = rng.nextDouble();
            }
        }
        
        //y_test = np.dot(X_test, w)
        double[] ytest = new double[nSamples];
        for (int row = 0; row < nSamples; ++row) {
            double sum = 0;
            for (int col = 0; col < nFeatures; ++col) {
                sum += X[row][col] * w[col];
            }
            ytest[row] = sum;
        }
        
        Data data = new Data();
        data.X = X;
        data.Xtest = Xtest;
        data.y = y;
        data.ytest = ytest;
        data.w = w;
                
        return data;
    }
    
    private static class Data {
        double[][] X;
        double[] y;
        double[][] Xtest;
        double[] ytest;
        double[] w;
    }
   
}
