package thirdparty.scipy.optimization;

import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TIntArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.logging.Logger;

/**
 * A port and adaptation of the scipy scikit version of the elastic-net algorithm
  which has copyright:

   scikit-learn/scikit-learn is licensed under the
   BSD 3-clause "New" or "Revised" License

   A permissive license similar to the BSD 2-Clause License, but with a
   3rd clause that prohibits others from using the name of the project
   or its contributors to promote derived products without written consent.

   https://github.com/scikit-learn/scikit-learn/blob/master/COPYING
 <pre>
   Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
          Fabian Pedregosa <fabian.pedregosa@inria.fr>
          Olivier Grisel <olivier.grisel@ensta.org>
          Gael Varoquaux <gael.varoquaux@inria.fr>

  License: BSD 3 clause
 </pre>

 NOT READY FOR USE
 */
public class ElasticNet {

    //have removed precompute logic to simplify the code for a "toy" version.
    // NOTE that a post suggets that the resulting gram version of
    //  corrdinate descent was slower to update in any case.

    private final double eps;
    private int nMaxIter;
    private boolean debug = true;

    /*
    alpha = a + b and l1_ratio = a / (a + b).
    when l1_ratio==1, the penalty is the lasso penalty.
    note, should not use l1_ratio <= 0.01
    */
    private final double alpha;
    private final double l1Ratio;

    private double tol = 1.e-4;

    private boolean doFitIntercept = true;
    private boolean doNormalize = false;
    
    // if doFitIntercept is true, this is populated
    private double intercept = 0.;

    private double[] coef;

    private double[][] X = null;
    private double[] y = null;

    // X.T * y, so is one dimension also
    private double[] Xy = null;

    // coef selection = 0 for 'cyclic' else '1' for random.
    public static enum Selection {
        CYCLIC, RANDOM;
    }
    private Selection selection = Selection.CYCLIC;

    /**
     When set to ``True``, reuse the solution of the previous call to fit as
        initialization, otherwise, just erase the previous solution.
    */
    private boolean warmStart = false;

    // dualGaps is a scalar because y is one dimensional in this edted class
    private double dualGap;
    private TIntArrayList nIter = new TIntArrayList();
    private boolean positive = false;

    /*
    random_state : int, RandomState instance, or None (default)
                The seed of the pseudo random number generator that selects
                a random feature to update. Useful only when selection is set to
                'random'.
    */
    private Random rng = null;
    private long ranSeed = 0;

    private int verbose = 0;

    private Logger log = Logger.getLogger(this.getClass().getName());

    /**
     alpha = a + b and l1_ratio = a / (a + b).
      when l1_ratio .eq. 1, the penalty is the lasso penalty.
      note, should not use l1_ratio .lte. 0.01.
     NOTE also, that for alpha .eq. 0, the algorithm does not converge well.
     * @param alpha
     * @param l1Ratio
     */
    public ElasticNet(double alpha, double l1Ratio) {

        /*
        default scipy values:
         alpha=1.0, l1_ratio=0.5, fit_intercept=True,
         normalize=False, max_iter=1000,
         copy_X=True, tol=1e-4, warm_start=False, positive=False,
         random_state=None, selection='cyclic'):
        */

        this.alpha = alpha;
        this.l1Ratio = l1Ratio;

        this.eps = 1.e-3;

        this.nMaxIter = 1000;
    }

    /**
      alpha = a + b and l1_ratio = a / (a + b).
      when l1_ratio .eq. 1, the penalty is the lasso penalty.
      note, should not use l1_ratio .lte. 0.01.
      NOTE also, that for alpha .eq. 0, the algorithm does not converge well.
     * @param alpha
     * @param l1Ratio
     * @param maxIter
     */
    public ElasticNet(double alpha, double l1Ratio, int maxIter) {

        /*
        default scipy values:
         alpha=1.0, l1_ratio=0.5, fit_intercept=True,
         normalize=False, max_iter=1000,
         copy_X=True, tol=1e-4, warm_start=False, positive=False,
         random_state=None, selection='cyclic'):
        */

        this.alpha = alpha;
        this.l1Ratio = l1Ratio;

        this.eps = 1.e-3;

        this.nMaxIter = maxIter;
    }

    /**
     * default is 0 produces no extra logging, v=1 results in more logging,
     * while v >= 2 produces most logging.
     * @param v
     */
    public void setVerbosity(int v) {
        this.verbose = v;
    }

    /**
     *
     * @param x each column of x is the data and any operation performed
     * upon it.  for example, to fit a 2nd order polynomial,
     * x[*][0] is xdata[*] and x[*][1] is xdata[*]^2.
     * @param y
     */
    public void fit(double[][] x, double[] y) {

        //NOTE: should consider changing the format of X,
        //   would be more convenient to access all data
        //   for a single coefficient as x[coefIdx],
        //   that is, the data for a single coefficient
        //   would be present in a row

        this.X = Arrays.copyOf(x, x.length);
        this.y = Arrays.copyOf(y, y.length);

        /*
        Coordinate descent:
        each column of data is solved one at a time (in orther words,
        each coefficient is estimated separately in a random or cyclical
        pattern).
        */

        int nXData = this.X.length;
        int nCoef = this.X[0].length;

        PreFitResults preFitResults = preFit(this.X, this.y, null,
            doNormalize, doFitIntercept, true);

        this.X = preFitResults.X;
        this.y = preFitResults.y;
        double[] XMean = preFitResults.XMean;
        double[] XStdv = preFitResults.XStdv;
        this.Xy = preFitResults.Xy;
        double yMean = preFitResults.yMean;

        /*
        leaving this here to edit in future if want to change to allow y
           to be 2 dimensional rather than a vector.
        if y.ndim == 1:
            y = y[:, np.newaxis]
        if Xy is not None and Xy.ndim == 1:
            Xy = Xy[:, np.newaxis]
        */

        int nSamples = nXData;
        int nFeatures = nCoef;

        // length of 2nd dimension of y, which is restricted to 1 dimension
        //    for now
        final int nTargets = 1;

        if (!warmStart || coef == null) {
            coef = new double[nCoef];
        } else {
            // copy to detach from instance present elsewhere?
            coef = Arrays.copyOf(coef, coef.length);
        }

        this.nIter = new TIntArrayList(nTargets);

        double[] thisXy = null;
        double[][]thisCoef = null;
        double[] thisDualGaps = null;
        TIntArrayList thisNIter = null;

        int nAlphas = 0;
        double[] alphas = new double[]{alpha};

        // edited for y being one dimension
        if (Xy != null) {
            thisXy = Arrays.copyOf(Xy, Xy.length);
        } else {
            // call to enet_path on line 272 of coordinate_descent.py
            PathResults pathResults =
                enetPath(
                    X, y,
                    l1Ratio,
                    0, //eps=0,
                    nAlphas,
                    alphas,
                    thisXy,
                    false, //fit_intercept=False,
                    false, //normalize=False,
                    true, // copy_X=True,
                    false,//verbose=False,
                    this.tol,
                    this.positive,
                    XMean, XStdv,
                    true, //return_n_iter=True,
                    this.coef, //coef_init=coef_[k],
                    this.nMaxIter, //max_iter=self.max_iter,
                    //random_state=self.random_state,
                    this.selection);

            //_, this_coef, this_dual_gap, this_iter = \

            thisCoef = pathResults.coefs;
            thisDualGaps = pathResults.dualGaps;
            thisNIter = pathResults.nIters;
        }

        //coef_[k] = this_coef[:, 0]
        for (int ii = 0; ii < coef.length; ++ii) {
            this.coef[ii] = thisCoef[ii][0];
        }
        this.dualGap = thisDualGaps[0];
        this.nIter.add(thisNIter.get(0));

        if (nTargets == 1) {
            //self.n_iter_ = self.n_iter_[0]
            this.nIter.remove(1, this.nIter.size() - 1);
        }

        setIntercept(XMean, yMean, XStdv);
    }

    public void setSeletion(Selection sel) {
        if (sel == null) {
            throw new IllegalArgumentException("sel cannot be null");
        }
        this.selection = sel;
    }

    /**
     * Predict using the linear model
     * 
     * @param testData 2 dimensional data of format [nSamples[nFeatures]
     * where features is the coefficients for polynomial regression, for
     * example.
     * @return 
     */
    public double[] predict(double[][] testData) {
        return decision_function(testData);            
    }
    
    /**
     * Predict confidence scores for samples.

        The confidence score for a sample is the signed distance of that
        sample to the hyperplane.
      
     * @param testData
     * @return 
     */
    private double[] decision_function(double[][] XD) {
    
        int nFeatures = coef.length;
        if (XD[0].length != nFeatures) {
            throw new IllegalArgumentException(
                "X has " + XD[0].length + 
                " features per sample; expecting " + nFeatures);
        }
        
        double[] scores = new double[XD.length];
        for (int row = 0; row < XD.length; ++row) {
            double sum = 0;
            for (int col = 0; col < nFeatures; ++col) {
                sum += XD[row][col] * coef[col];
            }
            scores[row] = sum + intercept;
        }

        return scores;
    }

    public double[] getCoef() {
        return coef;
    }
    
    public double getDualGap() {
        return dualGap;
    }

    public void setNMaxIter(int n) {
        this.nMaxIter = n;
    }

    public void setTol(double tolerance) {
        this.tol = tolerance;
    }

    public void setPositiveParam(boolean pos) {
        this.positive = pos;
    }
    
    public void setToUseWarmStart() {
        this.warmStart = true;
    }

    private CenResults _centerData(double[][] X2, double[] y2,
        boolean fitIntercept2, boolean normalize2,
        boolean doCopy) {

        boolean doSampleWeight = false;

        /*
        Centers data to have mean zero along axis 0. This is here because
        nearly all linear models will want their data to be centered.

        If sample_weight is not None, then the weighted mean of X and y
        is zero, and not the mean itself
        */

        int nXData = X2.length;
        int nCoef = X2[0].length;

        CenResults results = new CenResults();
        results.XMean = new double[nCoef];
        results.XStdv = new double[nCoef];

        if (doCopy) {
            results.X = MatrixUtil.copy(X2);
            results.y = Arrays.copyOf(y2, y2.length);
        } else {
           results.X = X2;
           results.y = y2;
        }

        if (fitIntercept2) {
            for (int i = 0; i < nCoef; ++i) {
                double sum = 0;
                for (int i0 = 0; i0 < nXData; ++i0) {
                    sum += results.X[i0][i];
                }
                double mean = sum / (double)nXData;
                results.XMean[i] = mean;
                for (int i0 = 0; i0 < nXData; ++i0) {
                    results.X[i0][i] -= mean;
                }
                if (normalize2) {
                    sum = 0;
                    for (int i0 = 0; i0 < nXData; ++i0) {
                        sum += results.X[i0][i] * results.X[i0][i];
                    }
                    sum = Math.sqrt(sum);
                    results.XStdv[i] = sum;
                    if (sum != 0) {
                        for (int i0 = 0; i0 < nXData; ++i0) {
                            results.X[i0][i] /= sum;
                        }
                    }
                } else {
                    results.XStdv[i] = 1;
                }
            }
            double sum = 0;
            for (int i0 = 0; i0 < nXData; ++i0) {
                sum += results.y[i0];
            }
            sum /= (double)nXData;
            results.yMean = sum;
            for (int i0 = 0; i0 < nXData; ++i0) {
                results.y[i0] -= sum;
            }
        } else {
            Arrays.fill(results.XMean, 0);
            Arrays.fill(results.XStdv, 1);
            results.yMean = 0;
        }

        return results;
    }

    private boolean allClose(double[] a, double value) {
        double rtol=1.e-5;
        double atol=1.e-8;
        for (double a0 : a) {
           if (!(Math.abs(a0 - value) <= (atol + rtol * Math.abs(value)))) {
               return false;
           }
        }
        return true;
    }
    private boolean allClose(double[] a, double[] b) {
        double rtol=1.e-5;
        double atol=1.e-8;
        for (int i = 0; i < a.length; ++i) {
            double a0 = a[i];
            double b0 = b[i];
           if (!(Math.abs(a0 - b0) <= (atol + rtol * Math.abs(b0)))) {
               return false;
           }
        }
        return true;
    }

    /**
     * NOTE that for now, y2 and this.y are restricted to being a single
     * dimension, that is, a vector.
     *
     * @param X2
     * @param y2
     * @param l1Ratio2
     * @param eps2
     * @param nAlphas2
     * @param alphas2
     * @param thisXy
     * @param fitIntercept2
     * @param normalize2
     * @param copyX2
     * @param verbose2
     * @param tol2
     * @param positive2
     * @param XMean2
     * @param XStdv2
     * @param returnNIter
     * @param coef2
     * @param nMaxIter2
     * @param selection2
     * @return
     */
    private PathResults enetPath(
        double[][] X2, double[] y2,
        double l1Ratio2, double eps2,
        int nAlphas2, double[] alphas2,
        double[] Xy2,
        boolean fitIntercept2, boolean normalize2, boolean copyX2,
        boolean verbose2, double tol2,
        boolean positive2,
        double[] XMean2, double[] XStdv2,
        boolean returnNIter, double[] coefInit,
        int nMaxIter2, Selection selection2) {

        /*
        Compute elastic net path with coordinate descent

        The elastic net optimization function varies for mono and
            multi-outputs.

        For mono-output tasks it is::
            1 / (2 * n_samples) * ||y - Xw||^2_2 +
            + alpha * l1_ratio * ||w||_1
            + 0.5 * alpha * (1 - l1_ratio) * ||w||^2_2

        For multi-output tasks it is::
            (1 / (2 * n_samples)) * ||Y - XW||^Fro_2
            + alpha * l1_ratio * ||W||_21
            + 0.5 * alpha * (1 - l1_ratio) * ||W||_Fro^2

        Where::

            ||W||_21 = \sum_i \sqrt{\sum_j w_{ij}^2}

            i.e. the sum of norm of each row.
        */

        //NOTE: have removed sparse matrix support

        if (fitIntercept2) {
            log.warning("scipy deprecates this feature, so will consider removing it"
                + "after a look at associcated logic");
        }

        if (copyX2) {
            X2 = Arrays.copyOf(X2, X2.length);
        }

        int nSamples = X2.length;
        int nFeatures = X2[0].length;

        //NOTE: in future, if edit for y dimension > 1, will need to make
        // changes here

        PreFitResults preFitResults =
            preFit(X2, y2, Xy2, normalize2, fitIntercept2, false);

        X2 = preFitResults.X;
        y2 = preFitResults.y;
        XMean2 = preFitResults.XMean;
        XStdv2 = preFitResults.XStdv;
        Xy2 = preFitResults.Xy;
        double yMean2 = preFitResults.yMean;

        if (alphas2 == null) {
            alphas2 = _alphaGrid(X2, y2, Xy2, l1Ratio2, false, eps2, nAlphas2,
                false, false);
        } else {
            Arrays.sort(alphas2);
        }

        nAlphas2 = alphas2.length;

        //tol = params.get('tol', 1e-4)
        //positive = params.get('positive', False)
        //max_iter = params.get('max_iter', 1000)
        TIntArrayList nIters2 = new TIntArrayList();
        double[] dualGaps2 = new double[nAlphas2];

        boolean useRandom = false;
        if (selection2.equals(Selection.RANDOM)) {
            if (rng == null) {
                this.rng = Misc.getSecureRandom();
                this.ranSeed = System.nanoTime();
                this.rng.setSeed(ranSeed);
            }
            useRandom = true;
        }

        //TODO: consider rewriting arrays so that a row
        //   holds the coefficients for a given alpha
        double[][] coefs = new double[nFeatures][nAlphas2];
        for (int i = 0; i < nFeatures; ++i) {
            coefs[i] = new double[nAlphas2];
        }

        if (coefInit == null) {
            //coef_ = np.asfortranarray(np.zeros(coefs.shape[:-1]))
            if (this.coef == null) {
                this.coef = new double[nFeatures];
            }
            for (int i = 0; i < nFeatures; ++i) {
                this.coef[i] = coefs[i][coefs.length - 1];
            }
        } else {
            //coef_ = np.asfortranarray(coef_init)
            if (this.coef == null) {
                this.coef = Arrays.copyOf(coefInit, coefInit.length);
            } else {
                System.arraycopy(coefInit, 0, this.coef, 0, coefInit.length);
            }
        }

        CDescResults model;

        //NOTE: have removed logic for multiple output, y dimension > 1
        //NOTE: have removed logic for return_models
        for (int i = 0; i < alphas2.length; ++i) {
            double a = alphas2[i];
            double l1Reg = a * l1Ratio2 * nSamples;
            double l2Reg = a * (1.0 - l1Ratio2) * nSamples;

            model = enet_coordinate_descent(
                this.coef, l1Reg, l2Reg, X2, y2, nMaxIter2, tol2,
                useRandom, positive2);

            this.coef = model.w;
            double dualGap = model.gap;

            //NOTE: assign tol to this.eps??
            tol2 = model.tol;
            int nIter2 = model.nIter;

            //n_iters.append(n_iter_)
            nIters2.add(nIter2);

            //coefs[..., i] = coef_
            for (int ii = 0; ii < coef.length; ++ii) {
                coefs[ii][i] = coef[ii];
            }

            //dual_gaps[i] = dual_gap_
            dualGaps2[i] = dualGap;

            if (dualGap > tol2) {
                log.warning("Objective did not converge." +
                    " You might want" +
                    " to increase the number of iterations");
            }

            if (verbose2 && verbose > 0) {
                if (verbose > 2) {
                    log.info("i=" + i + " nAlphas=" + nAlphas2);
                    log.info(model.toString());
                } else if (verbose > 1) {
                    log.info("i=" + i + " nAlphas=" + nAlphas2);
                } else {
                    log.info(".");
                }
            }
        }

        /*
        if return_n_iter:
            return alphas, coefs, dual_gaps, n_iters
        else:
            return alphas, coefs, dual_gaps
        */

        PathResults results = new PathResults();
        results.alphas = alphas2;
        results.coefs = coefs;
        results.dualGaps = dualGaps2;
        results.nIters = nIters2;

        return results;
    }

    /**
     *
     * @param X2
     * @param y2
     * @param Xy2  can be null
     * @param normalize2
     * @param fitIntercept2
     * @param doCopy
     * @return
     */
    private PreFitResults preFit(double[][] X2, double[] y2, double[] Xy2,
        boolean normalize2,
        boolean fitIntercept2, boolean doCopy) {

        //NOTE: should change local variable names to ones
        //   different from instance var names

        CenResults cr = _centerData(X2, y2, fitIntercept2, normalize2, doCopy);

        int nXData = X2.length;
        int nCoef = X2[0].length;

        X2 = cr.X;
        y2 = cr.y;
        double[] XMean2 = cr.XMean;
        double[] XStdv2 = cr.XStdv;
        double yMean2 = cr.yMean;

        //cannot use Xy if precompute is not Gram
        Xy2 = null;

        PreFitResults results = new PreFitResults();
        results.X = X2;
        results.XMean = XMean2;
        results.XStdv = XStdv2;
        results.Xy = Xy2;
        results.y = y2;
        results.yMean = yMean2;

        return results;
    }

    /**
     * Compute the grid of alpha values for elastic net parameter search
     * @param X2
     * @param Y2
     * @param Xy2 can be null
     * @param l1Ratio2
     * @param fitIntercept2
     * @param eps2
     * @param nAlphas2
     * @param normalize2
     * @param copyX2
     * @return
     */
    private double[] _alphaGrid(double[][] X2, double[] y2,
        double[] Xy2, double l1Ratio2, boolean fitIntercept2,
        double eps2, int nAlphas2, boolean normalize2,
        boolean copyX2) {

        /*
        def _alpha_grid(X, y, Xy=None, l1_ratio=1.0, fit_intercept=True,
            eps=1e-3, n_alphas=100, normalize=False, copy_X=True):
        */

        int nXData = X.length;
        int nSamples = y2.length;

        //NOTE: have removed sparse matrix possibility from class
        if (Xy2 == null) {
            if (copyX2 && fitIntercept2) {
                X2 = Arrays.copyOf(X2, X2.length);
            }
            CenResults cr = _centerData(X2, y2, fitIntercept2, normalize2, false);
            X2 = cr.X;
            y2 = cr.y;

            Xy2 = MatrixUtil.multiply(
                MatrixUtil.transpose(X2), y2);
        } else {
            Xy2 = Arrays.copyOf(Xy2, Xy2.length);
        }

        double alphaMax = Math.sqrt(MatrixUtil.dot(Xy2, Xy2)) /
            (nSamples * l1Ratio2);

        double start = Math.log(alphaMax * eps)/Math.log(10);

        double stop = Math.log(alphaMax)/Math.log(10);

        TDoubleList a2 =
            MiscMath.logspace((float)start, (float)stop, nAlphas2, true);

        return a2.toArray(new double[a2.size()]);
    }

    private static class CenResults {
        protected double[][] X;
        protected double[] y;
        protected double[] XMean;
        protected double yMean;
        protected double[] XStdv;
    }

    private static class PreFitResults extends CenResults {
        protected double[] Xy;
    }

    private static class PathResults {
        double[] alphas;
        double[][] coefs;
        double[] dualGaps;
        TIntArrayList nIters;
    }

    /**
    from scipy/scikit, method for searching solutions over coordinate
    (a.k.a. coefficient, a.k.a. feature) space.

    <pre>
    coordinate descent algorithm for Elastic-Net regression

    We minimize
        (1/2) * norm(y - X w, 2)^2 + alpha norm(w, 1) + (beta/2) norm(w, 2)^2
    </pre>
    https://github.com/scikit-learn/scikit-learn/blob/4d9a12d175a38f2bcb720389ad2213f71a3d7697/sklearn/linear_model/cd_fast.pyx

    Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
              Fabian Pedregosa <fabian.pedregosa@inria.fr>
              Olivier Grisel <olivier.grisel@ensta.org>
              Alexis Mignon <alexis.mignon@gmail.com>
              Manoj Kumar <manojkumarsivaraj334@gmail.com>

     License: BSD 3 clause

    */
    private CDescResults enet_coordinate_descent(
        double[] w,
        double alpha2, double beta2,
        double[][] X2, double[] y2,
        int nMaxIter2, double tol2,
        boolean useRandom, boolean positive2
        ) {

        int nSamples = X2.length;
        int nFeatures = X2[0].length; // nCoef
        int nTasks = y2.length;

        double d_w_tol = tol2;

        //cdef np.ndarray[floating, ndim=1] norm_cols_X = (X**2).sum(axis=0)
        double[] normColsX = new double[nFeatures];
        for (int i = 0; i < nFeatures; ++i) {
            double sum = 0;
            for (int i0 = 0; i0 < nSamples; ++i0) {
                sum += (X2[i0][i] * X2[i0][i]);
            }
            normColsX[i] = sum;
        }

        // init values of the residuals
        double[] r = new double[nSamples];
        double[] XtA = new double[nFeatures];

        // R = y - np.dot(X, w)
        for (int i = 0; i < nSamples; ++i) {

            //double v = dot(nFeatures,  & X_data[i], nSamples, w_data, 1);
            // X2[i] is column of length nFeatures, that is coefficients
            // w[i] is the coefficient
            double dot = 0;
            for (int j = 0; j < nFeatures; ++j) {
                dot += (X2[i][j] * w[j]);
            }

            r[i] = y2[i] - dot;
        }

        // tol *= np.dot(y, y)
        tol2 *= (MatrixUtil.dot(y2, y2));

        int ii, nIter2;
        double w_ii;
        double t;
        double gap = Double.NEGATIVE_INFINITY;

        for (nIter2 = 0; nIter2 < nMaxIter2; ++nIter2) {
            double w_max = 0.0;
            double d_w_max = 0.0;

            for (int fIter = 0; fIter < nFeatures; ++fIter) {

                if (useRandom) {
                    //ii = rand_int(n_features, rand_r_state)
                    ii = rng.nextInt(nFeatures);
                } else {
                    ii = fIter;
                }

                if (normColsX[ii] == 0.0) {
                    continue;
                }

                // Store previous value
                w_ii = w[ii];

                if (w_ii != 0.0) {
                    // R += w_ii * X[:,ii]
                    //axpy(nSamples, w_ii, X2[ii * nSamples], 1, r, 1)
                    for (int j = 0; j < nSamples; ++j) {
                        r[j] += (X2[j][ii] * w_ii);
                    }
                }

                //# tmp = (X[:,ii]*R).sum()
                //tmp = dot(nSamples,  & X_data[ii * nSamples], 1, R_data, 1)
                double tmp = 0;
                for (int j = 0; j < nSamples; ++j) {
                    tmp += (X2[j][ii] * r[j]);
                }

                if (positive2 && tmp < 0) {
                    w[ii] = 0.0;
                } else {
                    w[ii] =
                        (fsign(tmp) * fmax(Math.abs(tmp) - alpha2, 0)
                        / (normColsX[ii] + beta2));
                }

                if (w[ii] != 0.0) {
                    //# R -= w[ii] * X[:,ii] # Update residual
                    //axpy(nSamples, -w[ii],  & X_data[ii * nSamples], 1,
                    //    R_data, 1)
                    for (int j = 0; j < nSamples; ++j) {
                        r[j] -= (X2[j][ii] * w[ii]);
                    }
                }

                //# update the maximum absolute coefficient update
                double d_w_ii = Math.abs(w[ii] - w_ii);
                if (d_w_ii > d_w_max) {
                    d_w_max = d_w_ii;
                }

                if (Math.abs(w[ii]) > w_max) {
                    w_max = Math.abs(w[ii]);
                }
            }

            if (w_max == 0.0 ||
                d_w_max / w_max < d_w_tol ||
                nIter2 == nMaxIter2 - 1) {
                //# the biggest coordinate update of this iteration was smaller
                //# than the tolerance:
                //check the duality gap as ultimate
                //# stopping criterion

                //# XtA = np.dot(X.T, R) - beta * w
                for (int i = 0; i < nFeatures; ++i) {
                    //XtA[i] = dot(nSamples,  & X_data[i * nSamples],
                    //    1, R_data, 1)
                    // - beta * w[i];
                    double dot = 0;
                    for (int j = 0; j < nSamples; ++j) {
                        dot += (X2[j][i] * r[j]);
                    }
                    XtA[i] = dot - beta2 * w[i];
                }

                double dual_norm_XtA;
                if (positive2) {
                    dual_norm_XtA = max(XtA);
                } else {
                   dual_norm_XtA = absMax(XtA);
                }

                //# R_norm2 = np.dot(R, R)
                double rNorm2 = MatrixUtil.dot(r, r);

                //# w_norm2 = np.dot(w, w)
                double wNorm2 = MatrixUtil.dot(w, w);


                if (dual_norm_XtA > alpha2) {
                    t = alpha2 / dual_norm_XtA;
                    double aNorm2 = rNorm2 * (t * t);
                    gap = 0.5 * (rNorm2 + aNorm2);
                } else {
                    t = 1.0;
                    gap = rNorm2;
                }

                double l1Norm = asum(w);

                double rDotY = MatrixUtil.dot(r, y2);
                //# np.dot(R.T, y)
                gap += (alpha2 * l1Norm
                    - t * rDotY + 0.5 * beta2 * (1 + t*t) * (wNorm2));

                if (gap < tol2) {
                    //# return if we reached desired tolerance
                    break;
                }
            }
        }

        CDescResults results = new CDescResults();
        results.w = w;
        results.gap = gap;
        results.tol = tol2;
        results.nIter = nIter2 + 1;

        return results;
    }

    private static class CDescResults {
        double[] w = null;
        double gap = Double.NEGATIVE_INFINITY;
        double tol = Double.POSITIVE_INFINITY;
        int nIter = Integer.MIN_VALUE;

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("w=");
            if (w != null) {
                for (int i = 0; i < w.length; ++i) {
                    sb.append(String.format(" %.4f", (float)w[i]));
                }
            }
            sb.append(" gap=").append(String.format(" %.4f", (float)gap));
            sb.append(" tol=").append(String.format(" %.4f", (float)tol));
            sb.append(" nIter=").append(Integer.toString(nIter));

            return sb.toString();
        }
    }

    private double fsign(double f) {
        if (f == 0) {
            return 0;
        } else if (f > 0) {
            return 1.0;
        } else {
            return -1.0;
        }
    }

    private double fmax(double x, double y) {
        if (x > y) {
            return x;
        } else {
            return y;
        }
    }

    private double max(double[] a) {
        int i;
        double m = a[0];
        double d;
        for (i = 1; i < a.length; ++i) {
            d = a[i];
            if (d > m) {
                m = d;
            }
        }
        return m;
    }

    private double absMax(double[] a) {
        int i;
        double m = Math.abs(a[0]);
        double d;
        for (i = 1; i < a.length; ++i) {
            d = Math.abs(a[i]);
            if (d > m) {
                m = d;
            }
        }
        return m;
    }

    private double asum(double[] a) {
        double sum = 0;
        for (int i = 0; i < a.length; ++i) {
            sum += Math.abs(a[i]);
        }
        return sum;
    }

    private void setIntercept(double[] XMean, double yMean, double[] XStdv) {
        if (doFitIntercept) {
            for (int i = 0; i < coef.length; ++i) {
                coef[i] /= XStdv[i];
            }
            this.intercept = yMean - MatrixUtil.dot(XMean, coef);
        } else {
            this.intercept = 0.;
        }
    }
    
    public TIntArrayList _nIters() {
        return nIter;
    }
}
