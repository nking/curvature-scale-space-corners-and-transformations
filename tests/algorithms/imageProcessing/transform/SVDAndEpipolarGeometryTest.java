package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.transform.EpipolarTransformer.FundamentalMatrixStructures;
import algorithms.matrix.MatrixUtil;
import algorithms.statistics.Standardization;
import algorithms.util.FormatArray;
import algorithms.util.LinearRegression;
import algorithms.util.ResourceFinder;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import java.awt.Color;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class SVDAndEpipolarGeometryTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public SVDAndEpipolarGeometryTest() {
    }
    
    public void testHartley_house() throws Exception {
        
        int m = 3;
        int n = 16;//8;
        
        double[][] x1 = new double[m][n];
        double[][] x2 = new double[m][n];
        
        x1[0] = new double[]{282, 178, 265, 177, 161, 282, 49
                , 440, 249, 109, 196
                , 244, 171, 215, 144
                , 345
        };
        x1[1] = new double[]{40, 208, 228, 236, 248, 260, 281
                , 384, 458, 380, 72
                , 88, 122, 145, 155
                , 194
        };
        x1[2] = new double[n];  Arrays.fill(x1[2], 1.0);
        
        x2[0] = new double[]{281, 172, 265, 169, 150, 281, 37
                , 433, 270, 121, 183
                , 243, 157, 215, 129
                , 339
        };
        x2[1] = new double[]{40, 210, 226, 238, 250, 255, 288
                , 372, 455, 385, 72
                , 88, 122, 145, 156
                , 189
        };
        x2[2] = new double[n];  Arrays.fill(x2[2], 1.0);
        
        
        String fileName1 = "Hartley1997_house_01.png";
        String fileName2 = "Hartley1997_house_02.png";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);

        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        int image1Width = img1.getWidth();
        int image1Height = img1.getHeight();
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        int image2Width = img2.getWidth();
        int image2Height = img2.getHeight();
        
        //overplotSVD(x1, img1, "tmp_nc_01_svd");
        
        //overplotSVD(x2, img2, "tmp_nc_02_svd");
        
        DenseMatrix expectedFM  = new DenseMatrix(3, 3);
        expectedFM.set(0, 0, -9.796e-8); expectedFM.set(0, 1, 1.473e-6); expectedFM.set(0, 2, -6.660e-4);
        expectedFM.set(1, 0, -6.346e-7); expectedFM.set(1, 1, 1.049e-8); expectedFM.set(1, 2, 7.536e-3);
        expectedFM.set(2, 0,  9.107e-4); expectedFM.set(2, 1, -7.739e-3); expectedFM.set(2, 2, -2.364e-2);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                expectedFM.set(i, j, expectedFM.get(i, j)/expectedFM.get(2, 2));
            }
        }
        System.out.printf("expected FM=\n%s\n", FormatArray.toString(expectedFM, "%.3e"));
            
        EpipolarTransformer tr = new EpipolarTransformer();
        
        for (int tst = 0; tst < 2; tst++) {
            
            if (tst == 1) {
                // redo for only the first 8 points which are not as "generally distributed":
                x1 = MatrixUtil.copySubMatrix(x1, 0, 2, 0, 7);
                x2 = MatrixUtil.copySubMatrix(x2, 0, 2, 0, 7);
            }
            
            System.out.printf("%d points in Hartley 1997 house\n============\n", x1[0].length);
        
            DenseMatrix x1M = new DenseMatrix(x1);
            DenseMatrix x2M = new DenseMatrix(x2);
            DenseMatrix fm = null;

            //fm = tr.calculateEpipolarProjection(x1M, x2M);
            
            //List<DenseMatrix> fms = tr.calculateEpipolarProjectionFor7Points(x1M, x2M);
            //fm = fms.get(0);

            System.out.printf("de-normalized FM=\n%s\n", 
                FormatArray.toString(fm, "%.3e"));

            overplotEpipolarLines(fm, x1M, x2M, img1, img2, 
                image1Width, image1Height, image2Width, image2Height, 
                "hartley97_house__" + tst);

            Distances distances = new Distances();

            double tolerance = 2;

            EpipolarTransformationFit fit = distances.calculateError(fm, x1M, 
                x2M, ErrorType.DIST_TO_EPIPOLAR_LINE, tolerance);
            System.out.println("epipolar fit=" + fit.toString());

            fit = distances.calculateError(fm, x1M, 
                x2M, ErrorType.SAMPSONS, tolerance);
            System.out.println("samspon errors=" + fit.toString());
        }
    }
    
    public void _testNC() throws Exception {
        
        int m = 3;
        int n = 8;
        
        double[][] x1 = new double[m][n];
        double[][] x2 = new double[m][n];
        
        /*
        SVD(A).U == SVD(AA^T).U == SVD(AA^T).V
        SVD(A).V == SVD(A^TA).V == SVD(A^TA).U  
        */
        
        x1[0] = new double[]{262, 316, 260, 284, 234, 177, 216
           , 220
        };
        x1[1] = new double[]{356, 342, 305, 279, 217, 76, 63
            , 158
        };
        x1[2] = new double[n];  Arrays.fill(x1[2], 1.0);
        
        x2[0] = new double[]{156, 191, 153, 167, 135, 97, 119
                , 125
        };
        x2[1] = new double[]{308, 301, 270, 249, 202, 97, 83
                , 156
        };
        x2[2] = new double[n];  Arrays.fill(x2[2], 1.0);
        
        
        String fileName1 = "nc_book_01.png";
        String fileName2 = "nc_book_02.png";
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);

        Image img1 = ImageIOHelper.readImageAsGrayScale(filePath1);
        int image1Width = img1.getWidth();
        int image1Height = img1.getHeight();
        Image img2 = ImageIOHelper.readImageAsGrayScale(filePath2);
        int image2Width = img2.getWidth();
        int image2Height = img2.getHeight();
        
        //overplotSVD(x1, img1, "tmp_nc_01_svd");
        
        //overplotSVD(x2, img2, "tmp_nc_02_svd");
        
        DenseMatrix x1M = new DenseMatrix(x1);
        DenseMatrix x2M = new DenseMatrix(x2);
        
        DenseMatrix fm = null;
        EpipolarTransformer tr = new EpipolarTransformer();
        //fm = tr.calculateEpipolarProjection(x1M, x2M);
        List<DenseMatrix> fms = tr.calculateEpipolarProjectionFor7Points(
           x1M, x2M);
        fm = fms.get(0);
        
        System.out.printf("de-normalized FM=\n%s\n", 
            FormatArray.toString(fm, "%.3e"));
        
        overplotEpipolarLines(fm, x1M, x2M, img1, img2, 
            image1Width, image1Height, image2Width, image2Height, 
            "nc_book"); 
        
        //Distances distances = new Distances();
        
        //DenseMatrix matchedLeftXY = Util.rewriteInto3ColumnMatrix(leftTrueMatches);
        
        //EpipolarTransformationFit fit = distances.calculateError(fm, matchedLeftXY, 
        //    matchedRightXY, ErrorType.DIST_TO_EPIPOLAR_LINE, tolerance);
        
        /*PolygonAndPointPlotter plotter = 
            new PolygonAndPointPlotter(0, 512, 0, 512);
        plotter.addPlotWithLines(float[] xPoints, float[] yPoints, 
            float[] xLinePairs, float[] yLinePairs, String plotLabel);
        */
        
        /*PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(0, 512, 0, image1Height,
            xPoints, yPoints, null, null, xPoints, yPoints,
            "X fft_" + fileNameRoot);
        */
        
    }
    
    public void _testMisc() throws Exception {
        
        //http://www.cs.cmu.edu/~16385/s17/Slides/12.4_8Point_Algorithm.pdf
        double[][] F = new double[3][3];
        F[0] = new double[]{-0.00310695, 0.0025646, 2.96584};
        F[1] = new double[]{-0.028094, 0.00771621, 56.3813};
        F[2] = new double[]{13.1905, -29.2007, -9999.79};
        
        double[][] conjugateTransposeF = MatrixUtil.transpose(F);
        
        double[] x = new double[]{343.53, 221.7, 0};
        
        // epipolar line l' = F * x
        double[] ellPrime = new double[]{0.0295,0.9996, -265.1531};
        
        double[][] fTf_U = new double[3][3];
        fTf_U[0] = new double[]{-0.0013, 0.2586, -0.9660};
        fTf_U[1] = new double[]{0.0029, -0.9660, -0.2586};
        fTf_U[2] = new double[]{1.0000, 0.0032,  -0.0005};
        
        double[] rightEpipole = new double[]{1861.02, 498.21, 1.0};
        //EpipolarTransformer sTransformer = new EpipolarTransformer();
        
        //SVD(A).U is the same as SVD(AA^T).U
        //SVD(A).V is the same as SVD(A^TA).V
        
        double[][] ffT = MatrixUtil.multiply(F, MatrixUtil.transpose(F));
        double[][] fTf = MatrixUtil.multiply(MatrixUtil.transpose(F), F);
        
        SVD svd_f = SVD.factorize(new DenseMatrix(F));
        SVD svd_fT = SVD.factorize(new DenseMatrix(conjugateTransposeF));
        SVD svd_ffT = SVD.factorize(new DenseMatrix(ffT));
        SVD svd_fTf = SVD.factorize(new DenseMatrix(fTf));
        
        System.out.printf("SVD(F):\nU=%s\nV^T=%s\n\n", svd_f.getU().toString(), svd_f.getVt().toString());
        System.out.printf("SVD(F^T):\nU=%s\nV^T=%s\n\n", svd_fT.getU().toString(), svd_fT.getVt().toString());
        System.out.printf("SVD(FF^T):\nU=%s\nV^T=%s\n\n", svd_ffT.getU().toString(), svd_ffT.getVt().toString());
        System.out.printf("SVD(F^TF):\nU=%s\nV^T=%s\n\n", svd_fTf.getU().toString(), svd_fTf.getVt().toString());
        System.out.flush();
    }
   
    /*
    for more datasets:
    http://www.robots.ox.ac.uk/~vgg/data/
    */
    
    private void overplotSVD(double[][] input, Image img, String outfileName) 
        throws IOException, NotConvergedException {

        int m = input.length;
        int n = input[0].length;
        
        // input1 reformatted to hold all dimensions of point_0, all dimensions of point 1,
        //     through all dimensions of pointn-1)
        int nDimensions = m;
        double[] data = new double[nDimensions * n];
        int c = 0;
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < m; ++i) {
                data[c] = input[i][j];
                c++;
            }
        }
        double[] outputMean = new double[nDimensions];
        double[] outputStandardDeviation = new double[nDimensions];
        double[] normalizedData = Standardization.standardUnitNormalization(data, nDimensions, outputMean, outputStandardDeviation);
        
        double[][] normalizedInput = new double[m][n];
        c = 0;
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < m; ++i) {
                normalizedInput[i][j] = normalizedData[c];
                c++;
            }
        }
        System.out.printf("normalized mean=%s  stDev=%s\n", FormatArray.toString(outputMean, "%.4f"),
            FormatArray.toString(outputStandardDeviation, "%.4f"));
        
        {
            // expand scale just for plotting a point per pixel
            double factor = 100;
            double[][] qn1 = MatrixUtil.copy(normalizedInput);
            //double[][] qn1 = MatrixUtil.copy(input);
            //factor = 1.0;
            float[] qn1X = new float[n];
            float[] qn1Y = new float[n];
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    qn1[i][j] *= factor;
                }
            }
            for (int j = 0; j < n; ++j) {
               qn1X[j] = (float)qn1[0][j];
               qn1Y[j] = (float)qn1[1][j];
            }
            
            int xMin = -(int)Math.floor(factor * 2);
            int xMax = (int)Math.ceil(factor * 2);
            int yMin = xMin;
            int yMax = xMax;
            
            LinearRegression lr = new LinearRegression();
            String str = lr.plotTheLinearRegression(qn1X, qn1Y, xMin, xMax, yMin, yMax);
            float[] yInterceptAndSlope = lr.calculateTheilSenEstimatorParams(qn1X, qn1Y);
            System.out.printf("normalized yIntercept, slope=%s\n", Arrays.toString(yInterceptAndSlope));
            System.out.println("wrote to file: " + str);
            System.out.flush();
        }
        
        //input1 = MatrixUtil.multiply(input1, MatrixUtil.transpose(input1));
        
        for (int ii = 0; ii < n; ii++) {
            double x = input[0][ii];
            double y = input[1][ii];
            ImageIOHelper.addPointToImage((float) x, (float) y, img, 3,
                255, 0, 0);
        }

        /*
        SVD(A).U == SVD(AA^T).U == SVD(AA^T).V
        SVD(A).V == SVD(A^TA).V == SVD(A^TA).U
        */
        
        //double[][] input1Sq = MatrixUtil.multiply(input1, MatrixUtil.transpose(input1));
        //System.out.println("A is input^2");
        double[][] a = input;
        SVD svd = SVD.factorize(new DenseMatrix(a));
        System.out.printf("SVD(X):\nU=\n%s\nS=%s\n\nV^T=\n%s\n\n", 
            FormatArray.toString(svd.getU(), "%.4f"), 
            Arrays.toString(svd.getS()), 
            FormatArray.toString(svd.getVt(), "%.4f")
        );
        // principal directions = U
        double[][] u = MatrixUtil.convertToRowMajor(svd.getU());
        double[][] vT = MatrixUtil.convertToRowMajor(svd.getVt());
        double s0 = svd.getS()[0];
        double s1 = svd.getS()[1];
        double s2 = svd.getS()[2];
        double[] u0s0 = new double[]{u[0][0]*s0/u[2][0], u[1][0]*s0/u[2][0], u[2][0]*s0/u[2][0]};
        double[] u1s1 = new double[]{u[0][1]*s1/u[2][1], u[1][1]*s1/u[2][1], u[2][1]*s1/u[2][1]};
        double[] u2s2 = new double[]{u[0][2]*s2/u[2][2], u[1][2]*s2/u[2][2], u[2][2]*s2/u[2][2]};        
        double[][] av = MatrixUtil.multiply(a, vT);
                
        System.out.printf("x =\n%s\n", FormatArray.toString(a, "%.3f  "));
        System.out.printf("u0 * s0 =\n%s\n", FormatArray.toString(u0s0, "%.3f  "));
        System.out.printf("u1 * s1 =\n%s\n", FormatArray.toString(u1s1, "%.3f  "));
        System.out.printf("u2 * s2 =\n%s\n", FormatArray.toString(u2s2, "%.3f  "));
        System.out.printf("av =\n%s\n", FormatArray.toString(av, "%.3f  "));
        System.out.flush();  
        
        System.out.println("normalized: \n=============\n");
        {
            a = normalizedInput;
            
            //reduce dimensions to 2
            
            a = MatrixUtil.copySubMatrix(a, 0, a.length-2, 0, a[0].length-1);
            svd = SVD.factorize(new DenseMatrix(a));
            System.out.printf("SVD(X):\nU=\n%s\nS=%s\n\nV^T=\n%s\n\n",
                FormatArray.toString(svd.getU(), "%.4f"), 
                Arrays.toString(svd.getS()), 
                FormatArray.toString(svd.getVt(), "%.4f")
            );
            // principal directions = U
            u = MatrixUtil.convertToRowMajor(svd.getU());
            vT = MatrixUtil.convertToRowMajor(svd.getVt());
            
            s0 = svd.getS()[0];
            s1 = svd.getS()[1];
            u0s0 = new double[]{u[0][0]*s0, u[1][0]*s0};
            u1s1 = new double[]{u[0][1]*s1, u[1][1]*s1};
            av = MatrixUtil.multiply(a, MatrixUtil.transpose(vT));
                
            //MatrixUtil.transpose(vT[0]);
            double[][] v0 = new double[vT[0].length][1];
            for (int i = 0; i < vT[0].length; ++i) {
                v0[i][0] = vT[0][i];
            }
            double[][] av0 = MatrixUtil.multiply(a, v0);
            
            System.out.printf("a = x =\n%s\n", FormatArray.toString(a, "%.3f  "));
            System.out.printf("u0 * s0 =\n%s\n", FormatArray.toString(u0s0, "%.3f  "));
            System.out.printf("u1 * s1 =\n%s\n", FormatArray.toString(u1s1, "%.3f  "));
            System.out.printf("a * v =\n%s\n", FormatArray.toString(av, "%.3f  "));
            System.out.printf("a * v0 =\n%s\n", FormatArray.toString(av0, "%.3f  "));
            //System.out.printf("u * x2 =\n%s\n", toString(ux2, "%.3f  "));
            System.out.flush();  
            
            double[][] aT = MatrixUtil.transpose(a);
            double[][] u0 = new double[2][1];
            u0[0] = new double[]{u[0][0]}; u0[1] = new double[]{u[1][0]};
            double[][] u1 = new double[2][1];
            u1[0] = new double[]{u[0][1]}; u1[1] = new double[]{u[1][1]};
            
            // a^T * u0 then u1
            double[][] aTu0 = MatrixUtil.multiply(aT, u0);
            double[][] aTu1 = MatrixUtil.multiply(aT, u1);
            System.out.printf("a^T * u0 =\n%s\n", FormatArray.toString(aTu0, "%.3f  "));
            System.out.printf("a^T * u1 =\n%s\n", FormatArray.toString(aTu1, "%.3f  "));
            
        }
        // plot the columns of U
        /*The representation of lines in homogeneous projective coordinates
        is:   line a*x + b*y + c = 0
            | a |
            | b |
            | c |
        The line can be rewritten in slope, intercept form:
            y = intercept + slope * x
              = -(c/b) - slope*(a/b)*x

        written as homogenization form of lines:
            | -a/b |
        */
        
        /*
        Color clr = null;
        for (int ii = 0; ii < input2.numColumns(); ii++) {
            clr = getColor(clr);
            DenseMatrix epipolarLinesInLeft = 
                MatrixUtil.multiply(
                algorithms.matrix.MatrixUtil.transpose(fm), input2);
            PairIntArray leftLine = spTransformer.getEpipolarLine(
                epipolarLinesInLeft, image1Width, image1Height, ii);
            ImageIOHelper.addCurveToImage(leftLine, img1, 0,
                clr.getRed(), clr.getGreen(), clr.getBlue());
        }

        clr = null;
        for (int ii = 0; ii < input1.numColumns(); ii++) {
            clr = getColor(clr);
            DenseMatrix epipolarLinesInRight = MatrixUtil.multiply(fm, input1);
            PairIntArray rightLine = spTransformer.getEpipolarLine(
                epipolarLinesInRight, img2.getWidth(), img2.getHeight(), ii);
            ImageIOHelper.addCurveToImage(rightLine, img2, 0,
                clr.getRed(), clr.getGreen(), clr.getBlue());
        }
        */
        String dirPath = ResourceFinder.findDirectory("bin");
        String path = ImageIOHelper.writeOutputImage(dirPath + "/" + outfileName + ".png", img);
        
        System.out.println("writing to: " + path);
    }
    
    private void overplotEpipolarLines(DenseMatrix fm, PairFloatArray set1,
        PairFloatArray set2, Image img1, Image img2, int image1Width,
        int image1Height, int image2Width, int image2Height, String outfileNumber) 
        throws IOException {
        
        DenseMatrix input1 =
            Util.rewriteInto3ColumnMatrix(set1);

        DenseMatrix input2 =
            Util.rewriteInto3ColumnMatrix(set2);

        overplotEpipolarLines(fm, input1, input2, img1, img2, 
            image1Width, image1Height, image2Width, image2Height, outfileNumber);
    }
    
    private void overplotEpipolarLines(DenseMatrix fm, DenseMatrix input1,
        DenseMatrix input2, Image img1, Image img2, int image1Width,
        int image1Height, int image2Width, int image2Height, String outfileNumber) 
        throws IOException {
        
        for (int ii = 0; ii < input1.numColumns(); ii++) {
            double x = input1.get(0, ii);
            double y = input1.get(1, ii);
            ImageIOHelper.addPointToImage((float) x, (float) y, img1, 3,
                255, 0, 0);
        }
        for (int ii = 0; ii < input2.numColumns(); ii++) {
            double x2 = input2.get(0, ii);
            double y2 = input2.get(1, ii);
            ImageIOHelper.addPointToImage((float) x2, (float) y2, img2, 3,
                255, 0, 0);
        }

        EpipolarTransformer spTransformer = new EpipolarTransformer();

        Color clr = null;
        for (int ii = 0; ii < input2.numColumns(); ii++) {
            clr = getColor(clr);
            
            // epipolar lines in left image are the projections of the
            //    right image points.
            //    = F^T * x2
            DenseMatrix epipolarLinesInLeft = 
                MatrixUtil.multiply(
                algorithms.matrix.MatrixUtil.transpose(fm), input2);
            
            PairIntArray leftLine = spTransformer.getEpipolarLine(
                epipolarLinesInLeft, image1Width, image1Height, ii);
            
            ImageIOHelper.addCurveToImage(leftLine, img1, 0,
                clr.getRed(), clr.getGreen(), clr.getBlue());
        }

        clr = null;
        for (int ii = 0; ii < input1.numColumns(); ii++) {
            clr = getColor(clr);
            
            // epipolar lines in right image are the projections of the
            //    left image points.
            //    = F * x1
            
            DenseMatrix epipolarLinesInRight = MatrixUtil.multiply(fm, input1);
            
            PairIntArray rightLine = spTransformer.getEpipolarLine(
                epipolarLinesInRight, img2.getWidth(), img2.getHeight(), ii);
            
            ImageIOHelper.addCurveToImage(rightLine, img2, 0,
                clr.getRed(), clr.getGreen(), clr.getBlue());
        }

        String dirPath = ResourceFinder.findDirectory("bin");
        String fl1 =  dirPath + "/tmp_m_1_" + outfileNumber + ".png";
        String fl2 =  dirPath + "/tmp_m_2_" + outfileNumber + ".png";
        ImageIOHelper.writeOutputImage(fl1, img1);
        ImageIOHelper.writeOutputImage(fl2, img2);
        System.out.printf("writing files:\n  %s\n  %s\n", fl1, fl2);
    }
    
    private Color getColor(Color clr) {
        if ((clr == null) || clr.equals(Color.MAGENTA)) {
            return Color.BLUE;
        }
        if (clr.equals(Color.BLUE)) {
            return Color.PINK;
        } else if (clr.equals(Color.PINK)) {
            return Color.GREEN;
        } else if (clr.equals(Color.GREEN)) {
            return Color.RED;
        } else if (clr.equals(Color.RED)) {
            return Color.CYAN;
        } else if (clr.equals(Color.CYAN)) {
            return Color.MAGENTA;
        } else if (clr.equals(Color.MAGENTA)) {
            return Color.LIGHT_GRAY;
        } else {
            return Color.ORANGE;
        }
    }
    
    private double[][] dot(double[] v, double[][] m) {

        if (v == null || v.length == 0) {
            throw new IllegalArgumentException("v cannot be null or empty");
        }
        if (m == null || m.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        
        int mRows = m.length;

        int mCols = m[0].length;
        
        if (v.length != mRows) {
            throw new IllegalArgumentException(
                "the lenght of v must equal the number of rows of m");
        }
        
        double[][] c = new double[mRows][mCols];
        for (int row = 0; row < mRows; row++) {
            c[row] = new double[mCols];
        }

        for (int row = 0; row < mRows; row++) {
            for (int col = 0; col < mCols; col++) {
                c[row][col] = (m[row][col] * v[row]);
            }
        }

        return c;
    }
    
    public static void main(String[] args) {
        
        try {
            EpipolarTransformer test = new EpipolarTransformer();
            
            //test.testRANSAC();
                        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }
    }

}
