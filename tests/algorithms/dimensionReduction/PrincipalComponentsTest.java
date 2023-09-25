package algorithms.dimensionReduction;

import algorithms.correlation.BruteForce;
import algorithms.correlation.MultivariateDistance;
import algorithms.correlation.UnivariateDistance;
import algorithms.dimensionReduction.PrincipalComponents.PCAStats;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath0;
import algorithms.statistics.Standardization;
import algorithms.util.FormatArray;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

/**
 *
 * @author nichole
 */
public class PrincipalComponentsTest extends TestCase {
    
    public PrincipalComponentsTest(String testName) {
        super(testName);
    }
    
    public void testPCA() throws Exception {
        
        // from:
        // https://online.stat.psu.edu/stat505/book/export/html/670
        
        double[][] x = readPlaces();
        
        int i, j;
        
        //NOTE: to match the results of the psu tutorial, follow section
        //    Example 11-3: Place Rated (after Standardization)

        /*
        double[] mean = new double[x[0].length];
        double[] stDev = new double[x[0].length];
        x = Standardization.standardUnitNormalization(x, mean, stDev);
        System.out.printf("mean x=\n%s\n", FormatArray.toString(mean, "%.5e"));
        System.out.flush();
         */
        double[] mean = MatrixUtil.columnMeans(x);

        x = Standardization.zeroCenterMean(x);

        /*
        Step 1: Examine the eigenvalues to determine how many principal 
        components should be considered:
        
        Table 1. Eigenvalues and the proportion of variation explained by the 
        principal components.

            Component	Eigenvalue	Proportion	Cumulative
                1	3.2978	0.3664	0.3664
                2	1.2136	0.1348	0.5013
                3	1.1055	0.1228	0.6241
                4	0.9073	0.1008	0.7249
                5	0.8606	0.0956	0.8205
                6	0.5622	0.0625	0.8830
                7	0.4838	0.0538	0.9368
                8	0.3181	0.0353	0.9721
                9	0.2511	0.0279	1.0000	
            The first principal component explains about 37% of the variation. 
            Furthermore, the first four principal components explain 72%, while 
            the first five principal components explain 82% of the variation.         
        */

        System.out.printf("x is [%d X %x]\n", x.length, x[0].length);
        
        PCAStats stats = PrincipalComponents.calcPrincipalComponents(x, 3);

        System.out.printf("mean of A =\n%s\n", FormatArray.toString(mean, "%.4e"));
        System.out.printf("pA=\n%s\n", FormatArray.toString(stats.principalAxes, "%.4e"));
        System.out.printf("pC=\n%s\n", FormatArray.toString(stats.principalComponents, "%.4e"));

        double[] expectedEig = new double[]{0.3775, 0.0511, 0.0279, 0.0230, 0.0168, 0.0120, 0.0085, 0.0039, 0.0018};
        double[] expectedFracEig = new double[]{0.7227, 0.0977, 0.0535, 0.0440, 0.0321, 0.0229, 0.0162, 0.0075, 0.0034};
        double[] expectedCumulativeFracEig = new double[]{0.7227, 0.8204, 0.8739, 0.9178, 0.9500, 0.9728,
            0.9890, 0.9966, 1.00};
        double[] expectedPA0 = new double[]{
                0.03551, 0.0933, 0.4078, 0.1004, 0.1501, 0.0321, 0.8743, 0.1590, 0.0195
        };

        double tol = 1E-3;
        assertEquals(expectedEig.length, stats.eigenValues.length);
        for (i = 0; i < expectedEig.length; ++i) {
            assertTrue(Math.abs(expectedEig[i] - stats.eigenValues[i]) < tol);
        }


        /*
        correlation between stats.pC and the original variables:

        	            Principal Component
            Variable	1	2	3
            --------   --  --  --
            Climate	0.190	0.017	0.207
            Housing	0.544	0.020	0.204
            Health	0.782	-0.605	0.144
            Crime	0.365	0.294	0.585
            Transportation	0.585	0.085	0.234
            Education	0.394	-0.273	0.027
            Arts	0.985	0.126	-0.111
            Recreation	0.520	0.402	0.519
            Economy	0.142	0.150	0.239
         */
        // projected [297 X 9]
        // X [297 X 9]
        // correlation between projected and X
        int k;
        double[][] dCor = MultivariateDistance.fastDCor(x, stats.principalComponents);

        // correlation and fast distance correlation are similar, except dCor ~ abs(cor)
        System.out.printf("dCor=\n%s\n", FormatArray.toString(dCor, "%.4e"));

        double[][] expectedCorr = new double[9][];
        expectedCorr[0] = new double[]{0.190,	0.017,	0.207};
        expectedCorr[1] = new double[]{0.544,	0.020,	0.204};
        expectedCorr[2] = new double[]{0.782,	-0.605,	0.144};
        expectedCorr[3] = new double[]{0.365,	0.294,	0.585};
        expectedCorr[4] = new double[]{0.585,	0.085,	0.234};
        expectedCorr[5] = new double[]{0.394,	-0.273,	0.027};
        expectedCorr[6] = new double[]{0.985,	0.126,	-0.111};
        expectedCorr[7] = new double[]{0.520,	0.402,	0.519};
        expectedCorr[8] = new double[]{0.142,	0.150,	0.239};

        double diff, e;
        assertEquals(expectedCorr.length, dCor.length);
        int n10 = 0;
        int n20 = 0;
        int n30 = 0;
        int n40 = 0;
        for (i = 0; i < expectedCorr.length; ++i) {
            assertEquals(expectedCorr[i].length, dCor[i].length);
            for (j = 0; j < expectedCorr[i].length; ++j) {
                // these are roughly equivalent... diff is 10-20% of expected
                e = Math.abs(expectedCorr[i][j]);
                diff = Math.abs(e - dCor[i][j]);
                if (diff <= 0.1*e) {
                    ++n10;
                } else
                if (diff <= 0.2*e) {
                    ++n20;
                } else
                if (diff <= 0.3*e) {
                    ++n30;
                } else
                if (diff <= 0.4*e) {
                    ++n40;
                }
                //System.out.printf("%.3e %.3e => %.3e => %.3e\n", expectedCorr[i][j], dCor[i][j], diff, 0.1*e);
            }
        }
        System.out.printf("out of %d: %d within 10%% of expected, then %d within 20%% of expected, " +
                        "then %d within 30%% of expected, then %d within 40%% of expected\n",
                expectedCorr.length*expectedCorr[0].length, n10, n20, n30, n40);

        // correlation values furthest from 0 (positive or negative)
        // are the most strongly correlated.
        // for this table just printed, |correlation| >= 0.5 are significant.

        // principal component scores = principal components
        // X * pa^T
        //  [nx9] * [9x3] *= [nX3]
        //double[][] y = MatrixUtil.multiply(x, MatrixUtil.transpose(stats.principalAxes));
        //System.out.printf("principal component scores=\n%s\n", FormatArray.toString(y, "%.4e"));

        // variance-covariance(standardized data) == correlation(unstandardized data).
        //                                        == correlation(zero mean centered data).
        // therefore, pca using the standardized data == pca using the correlation matrix.
        // eigen of cov(standardized) == eigen of cor(unstandardized)

        double[][] x2 = readPlaces();
        double[][] cor = BruteForce.correlation(x);

        double[] outMean = new double[x2[0].length];
        double[] outS = new double[x2[0].length];
        double[][] z = Standardization.standardUnitNormalization(x2, outMean, outS);
        double[][] cov = BruteForce.covariance(z);

        EVD evdCov = EVD.factorize(new DenseMatrix(cov));
        EVD evdCor = EVD.factorize(new DenseMatrix(cor));
        System.out.printf("EVD(cov(z))=\n%s\n", FormatArray.toString(evdCov.getRealEigenvalues(), "%.3e"));
        System.out.printf("EVD(cor(x2s))=\n%s\n", FormatArray.toString(evdCor.getRealEigenvalues(), "%.3e"));

        /*
        cov = (1/(n-1)) * X^T*X where X has been zero-centered
        cor_i_j = cov(X_i, X_j) / (sigma_i * sigma_j)
        */

        System.out.println("unit standardized Z:");
        PCAStats stats2 = PrincipalComponents.calcPrincipalComponents(z, 2);
        System.out.printf("Z pA=\n%s\n", FormatArray.toString(stats2.principalAxes, "%.4e"));
        System.out.printf("Z pC=\n%s\n", FormatArray.toString(stats2.principalComponents, "%.4e"));

        /*
        following Jepson lecture to find coefficients a_j to minimize the objective:
        solve for a_j, knowing x_j and assuming B is stats.U_P
        // x_j is [n X 1]
        // B is [n X p]
        // a_j is [p X 1]
        objective is (x_j - B*a_j)^2
                  = (x_j)^2 - 2*(x_j)*(B*a_j) + (B*a_j)^2

        He proves that B is the first p columns of SVD(X).U, but doesn't solve for a_j, so here is my attempt
        to solve for a_j, then a comparison of the residuals using B*a_j to those using B*ones(p,1).

        looking for a closed form solution:
          set deriv of objective to 0:
            d/da_j(objective) = - 2*(x_j).*(B*ones(p, 1)) + 2*(B*a_j).*(B*ones(p, 1))

            in the factored, if the polynomial order of a_j was larger than 1 and linear in the orders,
            could group by polynomial powers of a_j
            and use Vandermonde matrix (V = matrix with columns of x^0  x^1  x^2).

            since the factored only has polynomial order of 1 for a_j and all other variables are known,
            can use DLT to form the solution in terms of M * a = 0

            deriv =  - 2*(x_j).*(B*ones(p, 1)) + 2*(B*a_j)*(B*ones(p, 1))

            2*B*a_j = |(2*b00*a0 + 2*b01*a1 + 2*b02*a2)|
                      |(2*b10*a0 + 2*b11*a1 + 2*b12*a2)|
                                ...
            B*ones(p, 1) = |(b00 + b01 + b02 +...)|
                           |(b10 + b11 + b12 +...)|
                               ...
         =>   -2*(x_j).*(B*ones(p, 1)) = |-2*x0*(b00 + b01 + b02 +...)|
                                       |-2*x1*(b10 + b11 + b12 +...)|
                                              ...
         =>   2*(B*a_j).*(B*ones(p, 1)) = |(2*b00*a0 + 2*b01*a1 + 2*b02*a2).*(b00 + b01 + b02 +...)|
                                        |(2*b10*a0 + 2*b11*a1 + 2*b12*a2).*(b10 + b11 + b12 +...)|
                                          ...
                                      = |( a0*(2*b00*(b00 + b01 + b02 +...)) + a1*(2*b01*(b00 + b01 + b02 +...))
                                           + a2*(2*b02*(b00 + b01 + b02 +...)) + ...) |
                                           ...
          grouping by terms of a:
          M rows are | (-2*x0*(b00 + b01 + b02 +...))    ( a0*(2*b00*(b00 + b01 + b02 +...))    a1*(2*b01*(b00 + b01 + b02 +...)) ...
             note that M will be [n x (p + 1)]

          then can use SVD(M).VT[n-1]
          or for an exact solution, can place the 0 order terms on the left as y in for M * a = y and solve a = pInv(M)*y

          a_j is one column and is [pX1]
          A is j columns of a_j and is [p X k] where k is x[0].length.
          then B*A = [nXp][pXk]=[nXk] which is same size as x

          evaluate the objective for each column, summing the square of the differences
         */
        double[][] b = stats.uP; // [n X p]
        int n = x.length;
        int p = stats.nComponents;;

        //double[][] a = MatrixUtil.zeros(p, x[0].length);
        double[][] aT = MatrixUtil.zeros(x[0].length, p);
        int ii;
        double sumB;
        double[] xj, baj;
        double[] ones = new double[p];
        Arrays.fill(ones, 1);
        double[] ssd = new double[n];
        double[] ssd1 = new double[n];
        double ssdTotal = 0;
        double[] b1 = MatrixUtil.multiplyMatrixByColumnVector(b, ones);
        for (j = 0; j < x[0].length; ++j) {
            // form the design matrix for DLT:
            double[][] mj = MatrixUtil.zeros(n, p+1);
            for (i = 0; i < x.length; ++i) {
                sumB = 0;
                for (ii = 0; ii < p; ++ii) {
                    sumB += b[i][ii];
                }
                sumB *= 2;
                //mj = | (-2*x0*(b00 + b01 + b02 +...))    ( a0*(2*b00*(b00 + b01 + b02 +...))    a1*(2*b01*(b00 + b01 + b02 +...))   ...
                mj[i][0] = -sumB * x[i][j];
                for (ii = 0; ii < p; ++ii) {
                    mj[i][ii+1] = sumB * b[i][ii];
                }
            }
            // mj is populated now, so solve for a_j
            SVD svd = SVD.factorize(new DenseMatrix(mj));
            //vT is p X p
            // set aj into column j of matrix a, which is row j of aT
            aT[j] = Arrays.copyOf(MatrixUtil.convertToRowMajor(svd.getVt())[p - 1], p);
            System.out.printf("a_%d = %s\n", j, FormatArray.toString(aT[j], "%.4e"));

            // evaluate the objective for this column j
            // x_j - B*a_j
            xj = MatrixUtil.extractColumn(x, j);
            baj = MatrixUtil.multiplyMatrixByColumnVector(b, aT[j]);
            for (i = 0; i < x.length; ++i) {
                diff = xj[i] - baj[i];
                ssd[j] += (diff * diff);
                diff = xj[i] - b1[i];
                ssd1[j] += (diff * diff);
            }
            // ssd[j] should be similar to ssd_p
            ssdTotal += ssd[j];
            // ssd1[j] does not use a_j.  can see that ssd1 is larger than ssd
            System.out.printf("ssd[%d]=%.4e  ssd1=%.4e total=%.4e\n", j, ssd[j], ssd1[j], ssdTotal);
        }
    }
    
    public void estPCA2() throws Exception {
        
        double[][] x = new double[4][2];
        x[0] = new double[]{1, 2};
        x[1] = new double[]{2, 1};
        x[2] = new double[]{3, 4};
        x[3] = new double[]{4, 3};

        double[][] xT = MatrixUtil.transpose(x);

        PolygonAndPointPlotter plot = new PolygonAndPointPlotter();
        plot.addPlot(xT[0], xT[1], null, null, "X");
        
        int i, j;
        
        double[][] xZC = Standardization.zeroCenterMean(x);

        double[][] xZCT = MatrixUtil.transpose(xZC);
        plot.addPlot(xZCT[0], xZCT[1], null, null, "X-m");

        double n = x.length;
        int p = 1;
        // x is 4 X 2.
        // p = 1
        // U is 4x4,  U_p is 4x1
        // V is 2x2, V_P is 2x1, V^T_p is 1x2
        // U_p * diag(s) = [4x1][1x1]=[4x1]
        // X * V_p = [4x2]*[2x1]=[4x1]
        PCAStats stats = PrincipalComponents.calcPrincipalComponents(xZC, p);

        plot.writeFile("_pca_");
    }

    public void testEigenFaces() throws NotConvergedException {

        // following Jepson for eigen eyes

        // images 25 x 20 vectorized to [1 X 500]

        //Jepson references:
        // Matthew Turk and Alex Pentland, Eigenfaces for Recognition, Journal
        // of Cognitive Neuroscience, 3(1), 1991, pp.71-86.
        // and
        // Brett Allen, Brian Curless, and Zoran Popovic,
        // The space of human body shapes: reconstruction and parameterization from range scans,
        // in ACM SIGGRAPH 2003, pp.27-31.

        /*
        a couple of datasets from kaggle with open or no licenses and somewhat small
        in size for use here in testing use of projection components and reconstruction:

        test dataset w/ public license
            https://www.kaggle.com/datasets/vikrantrajput/lungs-disease-data


        test dataset license Data files © Original Authors
        reference:
        @misc{pavel biza_2021,
                title={Female and male eyes},
                url={https://www.kaggle.com/ds/1438879},
            DOI={10.34740/KAGGLE/DS/1438879},
                    publisher={Kaggle},
                    author={Pavel Biza},
                    year={2021}
        */

        // read in images as greyscale, then downscale, then vectorize
        // store each vectorized image in a row
        // imgs is [nImages X imgArea]
        double[][] imgs = null;

        // mean of each column  np.mean(imgs, axis=0).  length = imgArea
        double[] meanCols = null;

        // reshape meanCols into an image and plot it

        // x = imgs - mean
        double[][] x = null;

        // use principal component analysis
        //reshape and plot the pA

        // reconstruction
        //x[0] is the first face.  its number of columns is [imgArea]
        // w = x[0] * U_p = [1 X imgArea] * [imgArea X p] = [1 X p]
        // projected face = w * U_p^T = [1 X p] * [p X imgArea] = [1 X imgArea]
        // reshape to image size and plot it
    }
    
    public void estReconstruction() throws Exception {
        
        double d2r = Math.PI/180.;
        double angle, dx, dy;
        double[][] x = new double[4][2];
       
        angle = -30;//-45;
        double x1=3./Math.sqrt(2); double y1=1./Math.sqrt(2);
        dx= x1*Math.cos(angle*d2r) + y1*Math.sin(angle*d2r);
        dy= -x1*Math.sin(angle*d2r) + y1*Math.cos(angle*d2r);
        double xt1 = dx;
        double yt1 = dy;
        
        y1*=-1.;
        dx= x1*Math.cos(angle*d2r) + y1*Math.sin(angle*d2r);
        dy= -x1*Math.sin(angle*d2r) + y1*Math.cos(angle*d2r);
        double xt2 = dx;
        double yt2 = dy;
        
        x1=7./Math.sqrt(2); y1*=-1.;
        dx= x1*Math.cos(angle*d2r) + y1*Math.sin(angle*d2r);
        dy= -x1*Math.sin(angle*d2r) + y1*Math.cos(angle*d2r);
        double xt3 = dx;
        double yt3 = dy;
        
        y1*=-1.;
        dx= x1*Math.cos(angle*d2r) + y1*Math.sin(angle*d2r);
        dy= -x1*Math.sin(angle*d2r) + y1*Math.cos(angle*d2r);
        double xt4 = dx;
        double yt4 = dy;
       
        x[0] = new double[]{xt1, yt1};
        x[1] = new double[]{xt2, yt2};
        x[2] = new double[]{xt3, yt3};
        x[3] = new double[]{xt4, yt4};
        
        int i, j;
        
        System.out.println("x0:");
        for (i = 0; i < x.length; ++i) {
            for (j = 0; j < x[i].length; ++j) {
                System.out.printf("%11.3e  ", x[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        
        double[] mean = new double[x[0].length];
        double[] stDev = new double[x[0].length];
        x = Standardization.standardUnitNormalization(x, mean, stDev);
        
        System.out.println("mean x=");
        for (i = 0; i < mean.length; ++i) {
            System.out.printf("%11.3e  ", mean[i]);
        }
        System.out.println();
        System.out.println("stdev x=");
        for (i = 0; i < stDev.length; ++i) {
            System.out.printf("%11.3e  ", stDev[i]);
        }
        System.out.println();
        System.out.flush();
        
        double n = x.length;
        
        // U is 2x2,  U_p is 2x1
        // V is 2x2, V_P is 2x1, V^T_p is 1x2
        // x is nx2
        PCAStats stats = PrincipalComponents.calcPrincipalComponents(x, 1);
                
        double[][] b = stats.principalComponents;
        System.out.printf("b = x * v^T_p=\n%s\n", FormatArray.toString(b, "%.5e"));
        System.out.flush();

        double[][] xMinusB = MatrixUtil.copy(x);
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                xMinusB[i][j] -= b[i][0];
            }
        }
        System.out.printf("x - Ba=\n%s\n", FormatArray.toString(xMinusB, "%11.3e"));
        System.out.flush();
        double xMBSum = 0;
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                xMBSum += xMinusB[i][j]*xMinusB[i][j];
            }
        }
        System.out.printf("sum (x - (m + B*a))^2 = %.5e\n", xMBSum);

        // =====================
        System.out.println("=== new points in same reference frame ===");
        
        angle = -30;//-45;
        x1=10./Math.sqrt(2); y1=1./Math.sqrt(2);
        dx= x1*Math.cos(angle*d2r) + y1*Math.sin(angle*d2r);
        dy= -x1*Math.sin(angle*d2r) + y1*Math.cos(angle*d2r);
        xt1 = dx;
        yt1 = dy;
        
        x1=10./Math.sqrt(2); y1=-1./Math.sqrt(2);
        dx= x1*Math.cos(angle*d2r) + y1*Math.sin(angle*d2r);
        dy= -x1*Math.sin(angle*d2r) + y1*Math.cos(angle*d2r);
        xt2 = dx;
        yt2 = dy;
        
        x1=15./Math.sqrt(2); y1=1./Math.sqrt(2);
        dx= x1*Math.cos(angle*d2r) + y1*Math.sin(angle*d2r);
        dy= -x1*Math.sin(angle*d2r) + y1*Math.cos(angle*d2r);
        xt3 = dx;
        yt3 = dy;
        
        x1=15./Math.sqrt(2); y1=-1./Math.sqrt(2);
        dx= x1*Math.cos(angle*d2r) + y1*Math.sin(angle*d2r);
        dy= -x1*Math.sin(angle*d2r) + y1*Math.cos(angle*d2r);
        xt4 = dx;
        yt4 = dy;
        
        x[0] = new double[]{xt1, yt1};
        x[1] = new double[]{xt2, yt2};
        x[2] = new double[]{xt3, yt3};
        x[3] = new double[]{xt4, yt4};
        
        System.out.println("*x0:");
        for (i = 0; i < x.length; ++i) {
            for (j = 0; j < x[i].length; ++j) {
                System.out.printf("%11.3e  ", x[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        
        mean = new double[x[0].length];
        stDev = new double[x[0].length];
        x = Standardization.standardUnitNormalization(x, mean, stDev);
        
        b = MatrixUtil.multiply(x, MatrixUtil.transpose(stats.principalAxes));

        double[][] bReconstruction = PrincipalComponents.reconstruct(mean, stats);
        for (i = 0; i < bReconstruction.length; ++i) {
            for (j = 0; j < bReconstruction[i].length; ++j) {
                bReconstruction[i][j] *= stDev[j];
                bReconstruction[i][j] += mean[j];
            }
        }
        
        System.out.println("*b = x * pr.dir * v^T_p=");
        for (i = 0; i < b.length; ++i) {
            for (j = 0; j < b[i].length; ++j) {
                System.out.printf("%11.3e  ", b[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        
        System.out.println("*B =");
        for (i = 0; i < bReconstruction.length; ++i) {
            for (j = 0; j < bReconstruction[i].length; ++j) {
                System.out.printf("%11.3e  ", bReconstruction[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        
        xMinusB = MatrixUtil.copy(x);
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                //xMinusB[i][j] -= mean[j];
                xMinusB[i][j] -= b[i][0];
            }
        }
        System.out.println("*x:");
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                System.out.printf("%11.3e  ", x[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        System.out.println("*(x - (m + B*a))^2:");
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                System.out.printf("%11.3e  ", xMinusB[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        
    }
    
    /**
     * get places data from a PSU statistics tutorial.  the
     * array has n=329 and dimensions=9.
     * the dimensions are 
     *     climate housing health crime trans educate arts recreate econ
     * @return
     * @throws IOException 
     */
    private double[][] readPlaces() throws IOException {
        // from:
        // https://online.stat.psu.edu/stat505/book/export/html/670
        
        double[][] x = new double[329][9];
        for (int i = 0; i < 329; ++i) {
            x[i] = new double[9];
        }
        
        String path = ResourceFinder.findFileInTestResources("places.txt");

        FileReader reader = null;
        BufferedReader in = null;
        int count = 0;
        try {
            in = new BufferedReader(new FileReader(new File(path)));
            
            String line = null;

            String pattern = "^(\\d+)";
            for (int i = 0; i < 9; ++i) {
                pattern = pattern + "\\s+(\\d+)";
            }
            pattern = pattern + "$";
            Pattern p = Pattern.compile(pattern);
            Matcher m = null;
            line = in.readLine();
            do {
                line = line.trim();
                //521  6200  237  923 4031 2757   996 1405 7633   1
                //521  6200  237  923 4031 2757   996 1405 7633   1
                //575  8138 1656  886 4883 2438  5564 2632 4350   2
                //468  7339  618  970 2531 2560   237  859 5250   3
                m = p.matcher(line);
                if (m.matches()) {
                    for (int c = 1; c <= 9; ++c) {
                        String s = m.group(c);
                        x[count][c-1] = Integer.valueOf(s);
                        x[count][c-1] = Math.log10(x[count][c-1]);
                    }
                    count++;
                }
                line = in.readLine();
            } while (line != null);
        } finally {
            if (in != null) {
                in.close();
            }
            if (reader != null) {
                reader.close();
            }
        }

        return x;
    }
}
