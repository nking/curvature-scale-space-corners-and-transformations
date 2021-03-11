package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.features.RANSACSolver;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.scaleSpace.CurvatureScaleSpaceCornerDetector;
import algorithms.matrix.MatrixUtil;
import algorithms.matrix.MatrixUtil.QAndR;
import algorithms.matrix.MatrixUtil.SVDProducts;
import algorithms.statistics.Standardization;
import algorithms.util.FormatArray;
import algorithms.util.LinearRegression;
import algorithms.util.ResourceFinder;
import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import gnu.trove.list.TIntList;
import java.awt.Color;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.Matrices;
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
    
    public void estHartley_house() throws Exception {
        
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
        Image _img1 = img1.copyImage();
        Image _img2 = img2.copyImage();
        
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
        
        
        // ----- ransac
        DenseMatrix x1M = new DenseMatrix(x1);
        DenseMatrix x2M = new DenseMatrix(x2);
        
        EpipolarTransformer.NormalizedXY normXY1 = EpipolarTransformer.normalize(x1M);
        EpipolarTransformer.NormalizedXY normXY2 = EpipolarTransformer.normalize(x2M);
        DenseMatrix leftM = normXY1.getXy();
        DenseMatrix rightM = normXY2.getXy();
        
        double tolerance = 3.84;// 5.99 7.82        
        boolean useToleranceAsStatFactor = true;
        ErrorType errorType = ErrorType.DIST_TO_EPIPOLAR_LINE;
        boolean reCalcIterations = false;
        
        RANSACSolver solver = new RANSACSolver();
        EpipolarTransformationFit fitR = solver.calculateEpipolarProjection(
            leftM, rightM, errorType, useToleranceAsStatFactor, tolerance,
                reCalcIterations);
        
        DenseMatrix fm = EpipolarTransformer.denormalizeTheFundamentalMatrix(
            fitR.getFundamentalMatrix(), normXY1.getNormalizationMatrix(),
            normXY2.getNormalizationMatrix());
                
        x1M = extractIndices(x1M, fitR.inlierIndexes);
        x2M = extractIndices(x2M, fitR.inlierIndexes);
        overplotEpipolarLines(fm, x1M, x2M, img1, img2, 
            image1Width, image1Height, image2Width, image2Height, 
            "hartley97_house__RANSAC");
        System.out.println("RANSAC fit=" + fitR.toString());
        System.out.println("RANSAC inlier indexes=" + 
            Arrays.toString(fitR.inlierIndexes.toArray()));
        System.out.println("RANSAC errors=" + 
            Arrays.toString(fitR.errors.toArray()));
        // -- end ransac
        
        Distances distances = new Distances();
        EpipolarTransformer tr = new EpipolarTransformer();
        EpipolarTransformationFit fit = null;
        int[] indexesToUse = null;
        
        for (int tst = 0; tst < 2; tst++) {
            img1 = _img1.copyImage();
            img2 = _img2.copyImage();
            
            //in un-normalized units:
            x1M = new DenseMatrix(x1);
            x2M = new DenseMatrix(x2);
            
            if (tst == 0) {
                indexesToUse = new int[leftM.numColumns()];
                for (int i = 0; i < indexesToUse.length; ++i) {
                    indexesToUse[i] = i;
                }
            } else if (tst == 1) {
                // redo for only the first 8 points which are not as "generally distributed":
                indexesToUse = new int[8];
                for (int i = 0; i < indexesToUse.length; ++i) {
                    indexesToUse[i] = i;
                }
            }
            
            System.out.printf("%d points in Hartley 1997 house\n============\n", x1[0].length);
        
            //un-normalized coordinates
            x1M = new DenseMatrix(Matrices.getSubMatrix(x1M, new int[]{0,1,2}, indexesToUse));
            x2M = new DenseMatrix(Matrices.getSubMatrix(x2M, new int[]{0,1,2}, indexesToUse));

            //noremalized coordinates
            DenseMatrix normX1M = new DenseMatrix(
                Matrices.getSubMatrix(normXY1.getXy(), new int[]{0,1,2}, indexesToUse));
            DenseMatrix normX2M = new DenseMatrix(
                Matrices.getSubMatrix(normXY2.getXy(), new int[]{0,1,2}, indexesToUse));
            
            useToleranceAsStatFactor = false;
            tolerance = 2;
            
            DenseMatrix normalizedFM = tr.calculateEpipolarProjection(normX1M, normX2M);

            assertNotNull(normalizedFM);
                
            fm = EpipolarTransformer.denormalizeTheFundamentalMatrix(
                normalizedFM, normXY1.getNormalizationMatrix(),
                normXY2.getNormalizationMatrix());
            
            //List<DenseMatrix> fms = tr.calculateEpipolarProjectionFor7Points(x1M, x2M);
            //fm = fms.get(0);

            System.out.printf("de-normalized FM=\n%s\n", 
                FormatArray.toString(fm, "%.3e"));

            overplotEpipolarLines(fm, x1M, x2M, img1, img2, 
                image1Width, image1Height, image2Width, image2Height, 
                "hartley97_house__" + tst);

            if (useToleranceAsStatFactor) {
                fit = distances.calculateError2(normalizedFM,
                    normXY1.getXy(), normXY2.getXy(),
                    errorType, tolerance);
            } else {
                fit = distances.calculateError(normalizedFM,
                    normXY1.getXy(), normXY2.getXy(),
                    errorType, tolerance);
            }
        
            System.out.println("errors=" + fit.toString());
        }
    }
    
    
    public void testNC() throws Exception {
        
        int m = 3;
        int n = 24;
        
        double[][] x1 = new double[m][n];
        double[][] x2 = new double[m][n];
        
        /*
        NOTE: because we know ahead of time that the objects to match in the
        images have corners (unlike the android statues), one could use a 
        corner detector and then HOGs to make a correspondence list.
        
        If the objects do not necessarily have corners, this project has
        code to use MSER to find candidate objects, then HOGs to match
        those.   One could make a correspondence list after that if still
        useful.
        */
        
        /*
        SVD(A).U == SVD(A^T).V == SVD(AA^T).U == SVD(AA^T).V
        SVD(A).V == SVD(A^T).U == SVD(A^TA).V == SVD(A^TA).U  
        */
        
        x1[0] = new double[]{262, 316, 260, 284, 234, 177, 216
           , 220, 248, 248, 319, 159, 176, 407, 393, 119, 117, 428, 427, 112, 109, 425, 256, 411
        };
        x1[1] = new double[]{356, 342, 305, 279, 217, 76, 63
            , 158, 27, 46, 36, 54, 76, 64, 85,       115, 141, 120, 147, 320, 375, 333, 192, 213
        };
        x1[2] = new double[n];  Arrays.fill(x1[2], 1.0);
        
        x2[0] = new double[]{156, 191, 153, 167, 135, 97, 119
            , 125, 137, 137, 183, 87, 97, 248, 238, 71, 70, 267, 267, 73, 74, 272,    147, 256
        };
        x2[1] = new double[]{308, 301, 270, 249, 202, 97, 83
            , 156, 51,  65,   49,  83, 97, 61,  81, 130, 149, 110, 134, 276, 315, 299, 180, 192
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
        
        EpipolarTransformer.NormalizedXY normXY1 = EpipolarTransformer.normalize(x1M);
        EpipolarTransformer.NormalizedXY normXY2 = EpipolarTransformer.normalize(x2M);
        DenseMatrix leftM = normXY1.getXy();
        DenseMatrix rightM = normXY2.getXy();
        
        double tolerance = 3.84; //3.84 5.99 7.82        
        boolean useToleranceAsStatFactor = true;
        ErrorType errorType = ErrorType.SAMPSONS;
        EpipolarTransformationFit fitR = null;
        boolean reCalcIterations = false;
        EpipolarTransformer tr = new EpipolarTransformer();
        
        /*DenseMatrix normalizedFM = tr.calculateEpipolarProjection(leftM, rightM);
        DenseMatrix vNFM = tr.validateSolution(normalizedFM, leftM, rightM);
        
        Distances distances = new Distances();
        if (useToleranceAsStatFactor) {
            fitR = distances.calculateError2(normalizedFM, leftM, rightM,
                    errorType, tolerance);
        } else {
            fitR = distances.calculateError(normalizedFM, leftM, rightM,
                    errorType, tolerance);
        }*/
        
        RANSACSolver solver = new RANSACSolver();
        fitR = solver.calculateEpipolarProjection(
            leftM, rightM, errorType, useToleranceAsStatFactor, tolerance,
                reCalcIterations);
        
        DenseMatrix fm = EpipolarTransformer.denormalizeTheFundamentalMatrix(
            fitR.getFundamentalMatrix(), 
            normXY1.getNormalizationMatrix(),
            normXY2.getNormalizationMatrix());
                
        x1M = extractIndices(x1M, fitR.inlierIndexes);
        x2M = extractIndices(x2M, fitR.inlierIndexes);
        
        System.out.println("RANSAC fit=" + fitR.toString());
        
        //System.out.println("RANSAC inlier indexes=" + 
        //    Arrays.toString(fitR.inlierIndexes.toArray()));
        //System.out.println("RANSAC errors=" + 
        //    Arrays.toString(fitR.errors.toArray()));
        System.out.println("inliers=");
        for (int i = 0; i < x1M.numColumns();++i) {
            if (i > 0) {
                System.out.printf(", ");
            }
            System.out.printf(" (%.0f, %.0f; %.0f, %.0f)", 
                x1M.get(0, i), x1M.get(1, i), x2M.get(0, i), x2M.get(1, i)); 
        }
        System.out.println();
        
        double[][] _fm = MatrixUtil.convertToRowMajor(fm);
        
        if (true) {
            // fix the solution to examine K
            _fm[0] = new double[]{5.8661e-07, -3.9359e-05, 9.2096e-03};
            _fm[1] = new double[]{3.6772e-05, 1.6925e-06, -1.8747e-02};
            _fm[2] = new double[]{-8.1094e-03, 1.2621e-02, 1.0000e+00};
            fm = new DenseMatrix(_fm);
        }
        
        overplotEpipolarLines(fm, x1M, x2M, img1, img2, 
            image1Width, image1Height, image2Width, image2Height, 
            "nc__RANSAC");
        
        System.out.printf("de-normalized FM=\n%s\n", 
            FormatArray.toString(fm, "%.4e"));
        double[][] fEpipoles = tr.calculateEpipoles(fm);
        
        //512x384
        System.out.printf("FM e1 = %s\n", 
                FormatArray.toString(fEpipoles[0], "%.4e"));
        System.out.printf("FM e2 = %s\n", 
                FormatArray.toString(fEpipoles[1], "%.4e"));
        
        //TODO: implement algorithm in section 3.3 of Hartley 1992 
        // "Estimation of Relative Camera Positions for Uncalibrated Cameras"
        
        /*
        pixel width = 1.4e-3mm
        image dimensions originally 4032x3024, 
            then divided by 7.875 to 512x384, 300 dpi
        FOV = 77 degrees = 1.344 radians
        focal length = (4032/2.) / tan(1.344/2.)
                     ~ 1604 pixels = 2.245 mm
           //NOTE: a webstie mentioned a focal length of 27 mm whichis 19286 pixels
        intrinsic camera matrix 
        K = | 1604     0  2016 |
            |    0  1604  1512 |
            |    0     0     1 |
        */
        // 
        double focalLength = 1604;//1604;  19286
        System.out.printf("focal Length=%.1f pixels\n", focalLength);
        double xC = 4032./2.;
        double yC = 3024./2.;
        double scale1 = 512./4032.;// 1./7.875;
        double[][] k = new double[3][3];
        k[0] = new double[]{-focalLength, 0, xC};
        k[1] = new double[]{0, -focalLength, yC};
        k[2]= new double[]{0, 0, 1};
       
        double[][] kScaled = MatrixUtil.copy(k);
        MatrixUtil.multiply(kScaled, scale1);
        kScaled[2][2] = 1.;
        System.out.printf("K/7.875=\n%s\n", FormatArray.toString(kScaled, "%.3e"));
        
        printMallonWhelanRectification(fm, fEpipoles, x1M, x2M, focalLength, xC, yC, scale1);
        
        printMonasseRectification(fm, fEpipoles, x1M, x2M, focalLength, xC, yC, scale1);
        
        // a quick look at
        // http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
        
        double[][] kInv = MatrixUtil.pseudoinverseFullRank(kScaled);
        //double[][] kInv = MatrixUtil.pseudoinverseRankDeficient(kScaled);
        
        double[][] k2TF = MatrixUtil.multiply(kInv, _fm);
        double[][] k2T_F_K1 = MatrixUtil.multiply(k2TF, kScaled);
        double[][] _essentialM = k2T_F_K1;
        System.out.printf("E=\n%s\n", FormatArray.toString(_essentialM, "%.3e"));
        
        /*
        3 by 4 matrix P = K * [R | T]
        
        M_intr is K with negative focal lengths
                     | -f_x      0  o_x |
            M_intr = |    0   -f_y  o_y |
                     |    0      0    1 |
        
        M_extr = [R | T ].  can see DOF=11
                  | r00  r01  r02  t_x |     | R_0^T  t_x |
               =  | r10  r11  r12  t_y |  =  | R_1^T  t_y |
                  | r20  r21  r22  t_z |     | R_2^T  t_z |
        
        projective matrix is P, called M here
            M = M_intr * M_extr
                
        using projective space:
           | u |                      | X_w |
           | v | =  M_intr * M_extr * | Y_w |
           | w |                      | Z_w |
                                      |   1 |
        where x_img = u/w
              y_img = v/w
        */

        /*
        from Szeliski Sect 7.2 Structure From Motion:
        Note that the absolute distance between the two cameras can never be recovered
        from pure image measurements alone, regardless of how many cameras or points
        are used. Knowledge about absolute camera and point positions or distances,
        often called ground control points in photogrammetry, is always required to
        establish the final scale, position, and orientation.
                */
        
        // roughly, would expect for the projection, 
        //    a translation in x of about 100-150 pixels at image plane
        //   and a rotation around y between 30 and 45 degrees.
        // After have applied a scaling of the image by factor 1./7.875
        
        
        SVDProducts svdE = MatrixUtil.performSVD(_essentialM);
        
        assert(svdE.u[0].length == 3 && svdE.u.length == 3);

        System.out.printf("SVD.u=\n%s", FormatArray.toString(svdE.u, "%.3e"));
        System.out.printf("det(SVD.u)=%.2f\n", MatrixUtil.determinant(svdE.u));
        System.out.printf("SVD.vT=\n%s", FormatArray.toString(svdE.vT, "%.3e"));
        System.out.printf("det(SVD.vT)=%.2f\n", MatrixUtil.determinant(svdE.vT));
        
        // same as Hartley1992's E^T
        double[][] w = new double[3][3];
        w[0] = new double[]{0, -1, 0};
        w[1] = new double[]{1, 0, 0};
        w[2] = new double[]{0, 0, 1};
        double[][] eH92 = MatrixUtil.transpose(w);
        double[][] zH92 = new double[3][3];
        zH92[0] = new double[]{0, -1, 0};
        zH92[1] = new double[]{1, 0, 0};
        zH92[2] = new double[]{0, 0, 0};

        double[][] R1 = MatrixUtil.multiply(svdE.u, w);
        R1 = MatrixUtil.multiply(R1, svdE.vT);

        double[][] R2 = MatrixUtil.multiply(svdE.u, MatrixUtil.transpose(w));
        R2 = MatrixUtil.multiply(R2, svdE.vT);

        DenseMatrix uM = new DenseMatrix(svdE.u);
        double[] t1 = Matrices.getColumn(uM, 2).getData();

        double[] t2 = Matrices.getColumn(uM, 2).getData();
        MatrixUtil.multiply(t2, -1);
        
        // R1Hartley92 is the same as R2 of Kitani's lecture notes
        // R2Hartley92 is the same as R1 of Kitani's lecture notes
        double[][] R1Hartley92 = MatrixUtil.multiply(svdE.u, eH92);
        R1Hartley92 = MatrixUtil.multiply(R1Hartley92, svdE.vT);
        double[][] R2Hartley92 = MatrixUtil.multiply(svdE.u, w);
        R2Hartley92 = MatrixUtil.multiply(R2Hartley92, svdE.vT);
        double[][] SHartley92 = MatrixUtil.multiply(MatrixUtil.transpose(svdE.vT), zH92);
        SHartley92 = MatrixUtil.multiply(SHartley92, svdE.vT);
        double[][] Q1Hartley92 = MatrixUtil.multiply(R1Hartley92, SHartley92);
        double[][] Q2Hartley92 = MatrixUtil.multiply(R2Hartley92, SHartley92);
        
        // solution 1:  R1 and T1
        // solution 2:  R1 and T2
        // solution 3:  R2 and T2
        // solution 4:  R2 and T1

        // valid equation has det(R) = 1 (rotation and reflection)
        double detR1 = MatrixUtil.determinant(R1);
        double detR2 = MatrixUtil.determinant(R2);

        //Compute 3D point using triangulation, valid solution has positive Z value
        // (Note: negative Z means point is behind the camera )

        System.out.printf("R1=\n%s\n", FormatArray.toString(R1, "%.3e"));
        System.out.printf("R2=\n%s\n", FormatArray.toString(R2, "%.3e"));
        System.out.printf("t1=\n%s\n", FormatArray.toString(t1, "%.3e"));
        System.out.printf("t2=\n%s\n", FormatArray.toString(t2, "%.3e"));
        System.out.printf("det(R1)=%.3e\n", detR1);
        System.out.printf("det(R2)=%.3e\n\n", detR2);
        System.out.printf("R1_Hartly92=\n%s\n", FormatArray.toString(R1Hartley92, "%.3e"));
        System.out.printf("R2_Hartely92=\n%s\n", FormatArray.toString(R2Hartley92, "%.3e"));
        System.out.printf("S_Hartly92=\n%s\n", FormatArray.toString(SHartley92, "%.3e"));
        System.out.printf("Q1_Hartly92=\n%s\n", FormatArray.toString(Q1Hartley92, "%.3e"));
        System.out.printf("Q2_Hartly92=\n%s\n", FormatArray.toString(Q2Hartley92, "%.3e"));
        System.out.flush();
        
        
        // editing
        //(256, 192) -> (147, 180)
        //(411, 213) -> (256, 192)
        double[][] a = new double[3][2];      
        a[0] = new double[]{256, 411};    
        a[1] = new double[]{192, 213};  
        a[2] = new double[]{1, 1};
        
        PairIntArray xy1 = new PairIntArray();
        xy1.add(256, 192);
        xy1.add(411, 213);
        PairIntArray xy2 = new PairIntArray();
        xy2.add(147, 180);
        xy2.add(256, 192);
        printTransformation(xy1, xy2, 1.2); // -35 degrees is -0.61 radians
                
        
        /* euler transformations
        
        about z-axis (yaw):           about x-axis (roll):       about the y-axis (pitch):
            | cos φ   -sin φ    0 |    |    1       0       0 |  |  cos ψ    0  sin ψ |
            | sin φ    cos φ    0 |    |    0   cos θ   sin θ |  |      0    1      0 |
            |     0       0    1 |    |     0  -sin θ   cos θ |  | -sin ψ    0  cos ψ |        
        */
        
        
        //DenseMatrix p2 = tr.pFromF(fm);
        DenseMatrix p2 = tr.pFromF(new DenseMatrix(_essentialM));
        System.out.printf("P2=\n%s\n",  FormatArray.toString(p2, "%.3e"));
        
        {
            SVDProducts svdP = MatrixUtil.performSVD(p2);
            
            System.out.printf("SVD(p2).U==\n%s\n", FormatArray.toString(
                svdP.u, " %.3e"));
            System.out.printf("SVD(p2).V^T==\n%s\n", FormatArray.toString(
                svdP.vT, " %.3e"));
            
            // smallest eigenvalue's eigenvector in svdP is the
            //   camera center c (or rather, c*R I think...)
            
        }
                
        // p1 = [I | 0]; p2 = [ [e2]_x * F | e2 ]
        // x1 = P1*X
        // x2 = P2*X
        //  ====> X = pseudoInv(p2) * x2
        //            4X3 3XN --> 4XN
        
        //pseudoinverse(A) = inverse(A^T*A) * A^T
        double[][] p2Inv = MatrixUtil.pseudoinverseFullRank(MatrixUtil.convertToRowMajor(p2));
        
        double[][] XW = MatrixUtil.multiply(p2Inv, 
            MatrixUtil.convertToRowMajor(x2M));
        System.out.printf("X from inv(P2)*x2=\n%s\n", FormatArray.toString(XW, " %.0f"));
        
        // save the first that pass the tests for Z>=0.
        double[][] R = null;
        double[] T = null;
        String goodSolnLabel = null;
        
        double eps = 1.e-5;
        boolean goodSoln, allZsPos;
        XW = transformx(R1, t1, x2M);
        allZsPos = allZsArePositive(XW);
        goodSoln = ((Math.abs(detR1) - 1.) < 1e-5) && allZsPos;
        System.out.printf("Good Solution=%b\nX from inv(|[R1|t1])*x2=\n%s\n", goodSoln, FormatArray.toString(XW, " %.0f"));
        if (goodSoln && R == null) {
            goodSolnLabel = "R1, T1";
            R = R1;
            T = t1;
        } 
        
        XW = transformx(R1, t2, x2M);
        allZsPos = allZsArePositive(XW);
        goodSoln = ((Math.abs(detR1) - 1.) < 1e-5) && allZsPos;
        System.out.printf("Good Solution=%b\nX from inv(|[R1|t2])*x2=\n%s\n", goodSoln, FormatArray.toString(XW, " %.0f"));
        if (goodSoln && R == null) {
            goodSolnLabel = "R1, T2";
            R = R1;
            T = t2;
        } 
        
        XW = transformx(R2, t1, x2M);
        allZsPos = allZsArePositive(XW);
        goodSoln = ((Math.abs(detR1) - 1.) < 1e-5) && allZsPos;
        System.out.printf("Good Solution=%b\nX from inv(|[R2|t1])*x2=\n%s\n", goodSoln, FormatArray.toString(XW, " %.0f"));
        if (goodSoln && R == null) {
            goodSolnLabel = "R2, T1";
            R = R2;
            T = t1;
        } 
        
        XW = transformx(R2, t2, x2M);
        allZsPos = allZsArePositive(XW);
        goodSoln = ((Math.abs(detR1) - 1.) < 1e-5) && allZsPos;
        System.out.printf("Good Solution=%b\nX from inv(|[R2|t2])*x2=\n%s\n", goodSoln, FormatArray.toString(XW, " %.0f"));
        if (goodSoln && R == null) {
            goodSolnLabel = "R2, T2";
            R = R2;
            T = t2;
        } 
        
        // check the components of this one
        double estimatedRotP2 = Math.atan(p2.get(0, 2)/p2.get(0, 0)) * (180./Math.PI);
        System.out.printf("estimated rotation about y axis from P2=%.2f\n", estimatedRotP2);
        
        if (R == null) {
            //TODO: wrap the solution finiding in a retry.  may need a retry loop for adjusting f
            // once a good solution is found too.
            return;
        }
        System.out.println("choosing solution: " + goodSolnLabel);
        //double estimatedRotY = Math.atan(R[0][2]/R[0][0]) * (180./Math.PI);
        double estimatedRotY = Math.atan(-R[2][0]/R[2][2]) * (180./Math.PI);
        System.out.printf("estimated rotation about y axis from R=%.2f\n", estimatedRotY);
        System.out.flush();
        
                
        /* http://www.cs.cmu.edu/~16385/s17/Slides/13.1_Stereo_Rectification.pdf
        1. Compute E to get R
              Estimate E using the 8 point algorithm (SVD)
              Estimate the epipole e (SVD of E)
              Build Rrect from e
              Decompose E into R and T
              Set R1 = R_rect and R2 = R * R_rect
        2. Rotate right image by R
        3. Rotate both images by Rrect
        4. Scale both images by H
              ?? Rectified points as p = f/z’[x’ y’ z’]
        */
        
        double[][] eEpipoles = tr.calculateEpipoles(new DenseMatrix(_essentialM));
        
        double normForT = Math.sqrt(T[0]*T[0] + T[1]*T[1] + T[2]*T[2]);
        double rRect2Denom = Math.sqrt(T[0]*T[0] + T[1]*T[1]);
        double[] rRect1 = new double[]{T[0]/normForT, T[1]/normForT, T[2]/normForT};
        double[] rRect2 = new double[]{-T[1]/rRect2Denom, T[0]/rRect2Denom, 0};
        
        //double[] rRect1 = eEpipoles[0];
        //double rRect2Denom = Math.sqrt(eEpipoles[0][0]*eEpipoles[0][0] + 
        //    eEpipoles[0][1]*eEpipoles[0][1]);
        //double[] rRect2 = new double[]{-eEpipoles[0][1]/rRect2Denom, eEpipoles[0][0]/rRect2Denom, 0};
        double[] rRect3 = MatrixUtil.crossProduct(rRect1, rRect2);
        double[][] rRect = new double[3][3];
        rRect[0] = rRect1;
        rRect[1] = rRect2;
        rRect[2] = rRect3;
        
        System.out.printf("e1=\n%s\n", FormatArray.toString(eEpipoles[0], "%.3f"));
        System.out.printf("e2=\n%s\n", FormatArray.toString(eEpipoles[1], "%.3f"));
                
        System.out.printf("rRect=\n%s\n", FormatArray.toString(rRect, "%.3f"));
        
        // this should be [1, 0, 0]
        double[] rRectE1 = MatrixUtil.multiplyMatrixByColumnVector(rRect, 
            rRect1);
        System.out.printf("Expecting [1, 0, 0]:\n");
        System.out.printf("rRect*rRect[0]=\n%s\n", FormatArray.toString(rRectE1, "%.3f"));
        rRectE1 = MatrixUtil.multiplyMatrixByColumnVector(rRect, 
            eEpipoles[0]);
        System.out.printf("rRect*e1=\n%s\n", FormatArray.toString(rRectE1, "%.3f"));
        
        
        // [x’ y’ z’] = rRect * [x y z]
        //   (  f/z’ * ’[x’ y’ z’]  )
        // *may need to alter the focal length (inside K) to keep points within the original image size
        Image img2Rect = img2.createWithDimensions();
        double[] coords = new double[]{0, 0, 1};
        double[] coordsRect;
        int i2, j2;
        
        ///bounds (0,0)  (0, imageWidth-1, imageHeight-1)
        coords[0] = image2Width-1; coords[1] = image2Height-1;
        coordsRect = MatrixUtil.multiplyMatrixByColumnVector(rRect, coords);
        System.out.printf("\ncoordsR=%s\n",
            FormatArray.toString(coordsRect, "%.3e"));
        coords[0] = 1; coords[1] = 1;
        coordsRect = MatrixUtil.multiplyMatrixByColumnVector(rRect, coords);
        
        System.out.printf("\ncoordsR=%s\n",
            FormatArray.toString(coordsRect, "%.3e"));
        /* paused here
        for (int i = 1; i < image2Width; ++i) {
            coords[0] = i;
            for (int j = 1; j < image2Height; ++j) {
                coords[1] = j;
                coordsR = MatrixUtil.multiplyMatrixByColumnVector(rRect, coords);
                i2 = (int)Math.round(coordsR[0]);
                j2 = (int)Math.round(coordsR[1]);
                if (i2>0 && j2>0 && i2 < image2Width && j2 < image2Height) {
                    img2Rect.setRGB(i2, j2, img2.getRGB(i, j));
                }
            }
        }
        */
        // [x  y  z] = R * rRect * [x'  y'  z'] where R is the correct solution 
        //     of the 4 previous R1 and R2 combinations with T1 and T2 above
        
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
        
        System.out.printf("SVD(F):\nU=\n%s\nV^T=\n%s\n\n", 
            FormatArray.toString(svd_f.getU(), "%.4f"), 
            FormatArray.toString(svd_f.getVt(), "%.4f"));
        System.out.printf("SVD(F^T):\nU=\n%s\nV^T=\n%s\n\n", 
            FormatArray.toString(svd_fT.getU(), "%.4f"), 
            FormatArray.toString(svd_fT.getVt(), "%.4f"));
        System.out.printf("SVD(FF^T):\nU=\n%s\nV^T=\n%s\n\n", 
            FormatArray.toString(svd_ffT.getU(), "%.4f"), 
            FormatArray.toString(svd_ffT.getVt(), "%.4f"));
        System.out.printf("SVD(F^TF):\nU=\n%s\nV^T=\n%s\nD=\n%s\n\n", 
            FormatArray.toString(svd_fTf.getU(), "%.4f"), 
            FormatArray.toString(svd_fTf.getVt(), "%.4f"),
            FormatArray.toString(svd_fTf.getS(), "%.4f"));
        EVD evd_fTF = EVD.factorize(new DenseMatrix(fTf));
        System.out.printf("EVD(F^TF):\nU=\n%s\nV^T=\n%s\nD=\n%s\n\n", 
            FormatArray.toString(evd_fTF.getLeftEigenvectors(), "%.4f"), 
            FormatArray.toString(evd_fTF.getRightEigenvectors(), "%.4f"),
            FormatArray.toString(evd_fTF.getRealEigenvalues(), "%.4f"));
        
        EpipolarTransformer tr = new EpipolarTransformer();
        double[][] leftRightE = tr.calculateEpipoles(new DenseMatrix(F));
        System.out.printf("e1 from F = \n%s\n", FormatArray.toString(leftRightE[0], "%.4f"));
        System.out.printf("e2 from F = \n%s\n", FormatArray.toString(leftRightE[1], "%.4f"));
        System.out.flush();
        
        /*
        NOTE: when images are finally rectified, one should find
            FM_rect = 0   0   0
                      0   0  -1
                      0   1   0
        
            using H = 
        */
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
        SVD(A).U == SVD(A^T).V == SVD(AA^T).U == SVD(AA^T).V
        SVD(A).V == SVD(A^T).U == SVD(A^TA).V == SVD(A^TA).U
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

    private DenseMatrix extractIndices(DenseMatrix m, List<Integer> inlierIndexes) {
        DenseMatrix out = new DenseMatrix(m.numRows(), inlierIndexes.size());
        int r = 0;
        for (int i = 0; i < inlierIndexes.size(); ++i) {
            int idx = inlierIndexes.get(i);
            for (int j = 0; j < m.numRows(); ++j) {
                out.add(j, r, m.get(j, idx));
            }
            r++;
        }
        return out;
    }

    private void printTransformation(PairIntArray xy1, PairIntArray xy2, double scale) {
        
        /* euler transformations
        
        about z-axis (yaw):           about x-axis (roll):       about the y-axis (pitch):
            | cos φ   -sin φ    0 |    |    1       0       0 |  |  cos ψ    0  sin ψ |
            | sin φ    cos φ    0 |    |    0   cos θ   sin θ |  |      0    1      0 |
            |     0       0    1 |    |     0  -sin θ   cos θ |  | -sin ψ    0  cos ψ |
        
        */
        MatchedPointsTransformationCalculator c = new MatchedPointsTransformationCalculator();
        TransformationParameters params =  c.calulateEuclideanGivenScale(
            scale, xy1, xy2, xy1.getX(0), xy1.getY(0));
        
        // roation about z-axis:
        System.out.printf("transformation params (rot is for z-axis): \n%s\n", params.toString());
        
    }

    private double[][] transformx(double[][] R, double[] t, DenseMatrix x) throws NotConvergedException {
        double[][] P = concatenateColumns(R, t);
        double[][] pInv = MatrixUtil.pseudoinverseFullRank(P);
        
        double[][] XW = MatrixUtil.multiply(pInv, 
            MatrixUtil.convertToRowMajor(x));
        return XW;
    }
    
    private double[][] concatenateColumns(double[][] a, double[] b) {
        assertTrue(a.length == b.length);
        double[][] out = new double[a.length][a[0].length + 1];
        for (int i = 0; i < a.length; ++i) {
            out[i] = new double[a[0].length + 1];
            for (int j = 0; j < a[i].length; ++j) {
                out[i][j] = a[i][j];
            }
            out[i][a[0].length] = b[i];
        }
        return out;
    }

    private boolean allZsArePositive(double[][] XW) {
        return allOfDimensionArePositive(XW, 2);
    }
    private boolean allOfDimensionArePositive(double[][] XW, final int dimension) {
        for (int i = 0; i < XW[0].length; ++i) {
            if (XW[dimension][i] < 0) {
                return false;
            }
        }
        return true;
    }

    private void printMonasseRectification(DenseMatrix fm, double[][] fEpipoles, 
        DenseMatrix x1M, DenseMatrix x2M, 
        double focalLength, double xC, double yC, double scale1) throws NotConvergedException {
        
        double[][] k = new double[3][3];
        k[0] = new double[]{focalLength, 0, xC};
        k[1] = new double[]{0, focalLength, yC};
        k[2]= new double[]{0, 0, 1};
               
        double[][] kScaled = MatrixUtil.copy(k);
        MatrixUtil.multiply(kScaled, scale1);
        kScaled[2][2] = 1.;
                
        double[][] kInv = MatrixUtil.pseudoinverseFullRank(kScaled);
        //double[][] kInv = MatrixUtil.pseudoinverseRankDeficient(kScaled);
        
        double[][] fEpipolesGoal1 = MatrixUtil.copy(fEpipoles);
        fEpipolesGoal1[0][2] = 0;
        fEpipolesGoal1[1][2] = 0;
        
        double[][] fEpipolesGoal2 = MatrixUtil.copy(fEpipolesGoal1);
        for (int i = 0; i < 2; ++i) {
           fEpipolesGoal2[i][0] = 1;
           fEpipolesGoal2[i][1] = 0;
        }
        
        double[] a = MatrixUtil.multiplyMatrixByColumnVector(kInv, fEpipoles[0]);
        double[] b = MatrixUtil.multiplyMatrixByColumnVector(kInv, fEpipolesGoal1[0]);
        
        double aLength = MatrixUtil.lPSum(a, 2);
        double bLength = MatrixUtil.lPSum(b, 2);
        double ab = aLength*bLength;
        double aDotB = MatrixUtil.innerProduct(a, b);
        double[] aCrossB = MatrixUtil.crossProduct(a, b);
        double theta = Math.acos(aDotB/ab);
        double sinTheta = Math.sin(theta);
        double oneMinusCosTheta = 1. - Math.cos(theta);
        double[] rotAxis = new double[3];
        for (int i = 0; i < rotAxis.length; ++i) {
            rotAxis[i] = Math.acos(aCrossB[i]/ab);
        }
        System.out.printf("t rotation axis=\n%s\n", 
            FormatArray.toString(rotAxis, "%.3e"));
        
        // Rodrigues formula for small rotations:
        //   R(θ,t) = I + sinθ * [t]_× + (1−cosθ)*([t]_x)^2
        //      is [t]_x)^2 defined by power series P^(2*k+ 1) = ((-1)^k) *P?
        //      still reading...
        double[][] skewSymT = MatrixUtil.skewSymmetric(rotAxis);
        double[][] skewSymTT;// = skewSymT;//MatrixUtil.multiply(skewSymT, skewSymT);
        // Dai 2015, "Euler–Rodrigues formula variations, quaternion conjugation and intrinsic connections"
        // https://www.sciencedirect.com/science/article/pii/S0094114X15000415
        // the squared symmetric matrix given by eqn 14 is
        // s*s^T - I where s is the original vector
        
        // elsewhere, a definition with vector length squared times identity:
        //    [u×]2 = u ⊗ u − ||u||2 * I
        //       where ⊗ is the outer product
        double rotAxisLength = MatrixUtil.lPSum(rotAxis, 2);
        skewSymTT = MatrixUtil.outerProduct(rotAxis, rotAxis);
        for (int i = 0; i < skewSymTT.length; ++i) {
            //skewSymTT[i][i] -= 1;
            skewSymTT[i][i] -= (rotAxisLength*rotAxisLength);
        }
                
        System.out.printf("acos(aXb/(|a||b|))=\n%s\n", FormatArray.toString(rotAxis, "%.4e"));
        System.out.printf("skewSymT=\n%s\n", FormatArray.toString(skewSymT, "%.4e"));
        System.out.printf("skewSymTT=\n%s\n", FormatArray.toString(skewSymTT, "%.4e"));
        
        double[][] R_theta_t = MatrixUtil.zeros(3, 3);
        for (int i = 0; i < R_theta_t.length; ++i) {
            R_theta_t[i][i] = 1;
        }
        for (int i = 0; i < R_theta_t.length; ++i) {
            for (int j = 0; j < R_theta_t[i].length; ++j) {
                R_theta_t[i][j] += (sinTheta*skewSymT[i][j] + oneMinusCosTheta*skewSymTT[i][j]);
            }
        }
        System.out.printf("R_theta_t=\n%s\n", FormatArray.toString(R_theta_t, "%.4e"));
        
 
        double[] aRight = MatrixUtil.multiplyMatrixByColumnVector(kInv, fEpipoles[1]);
        double[] bRight = MatrixUtil.multiplyMatrixByColumnVector(kInv, fEpipolesGoal1[1]);
        
        double aLengthRight = MatrixUtil.lPSum(aRight, 2);
        double bLengthRight = MatrixUtil.lPSum(bRight, 2);
        double abRight = aLengthRight * bLengthRight;
        double aDotBRight = MatrixUtil.innerProduct(aRight, bRight);
        double[] aCrossBRight = MatrixUtil.crossProduct(aRight, bRight);
        double thetaRight = Math.acos(aDotBRight/abRight);
        double sinThetaRight = Math.sin(thetaRight);
        double oneMinusCosThetaRight = 1. - Math.cos(thetaRight);
        double[] rotAxisRight = new double[3];
        for (int i = 0; i < rotAxisRight.length; ++i) {
            rotAxisRight[i] = Math.acos(aCrossBRight[i]/abRight);
        }

        // Rodrigues formula for small rotations, but this angle is  large
        //   R(θ,t) = I + sinθ * [t]_× + (1−cosθ)*([t]_x)^2
        //      is [t]_x)^2 defined by power series P^(2*k+ 1) = ((-1)^k) *P?
        //      still reading...
        double rotAxisLengthRight = MatrixUtil.lPSum(rotAxisRight, 2);
        double[][] skewSymTRight = MatrixUtil.skewSymmetric(rotAxisRight);
        double[][] skewSymTTRight;// = skewSymTRight;//MatrixUtil.multiply(skewSymTRight, skewSymTRight);
        skewSymTTRight = MatrixUtil.outerProduct(rotAxisRight, rotAxisRight);
        for (int i = 0; i < skewSymTT.length; ++i) {
            //skewSymTTRight[i][i] -= 1;
            skewSymTTRight[i][i] -= (rotAxisLengthRight*rotAxisLengthRight);
        }
        
        System.out.printf("acos(aXb/(|a||b|)) Right =\n%s\n", FormatArray.toString(rotAxisRight, "%.4e"));
        System.out.printf("skewSymTRight=\n%s\n", FormatArray.toString(skewSymTRight, "%.4e"));
        System.out.printf("skewSymTTRight=\n%s\n", FormatArray.toString(skewSymTTRight, "%.4e"));
        
        double[][] R_theta_tRight = MatrixUtil.zeros(3, 3);
        for (int i = 0; i < R_theta_tRight.length; ++i) {
            R_theta_tRight[i][i] = 1;
        }
        for (int i = 0; i < R_theta_tRight.length; ++i) {
            for (int j = 0; j < R_theta_tRight[i].length; ++j) {
                R_theta_tRight[i][j] += (sinThetaRight * skewSymTRight[i][j] 
                    + oneMinusCosThetaRight * skewSymTTRight[i][j]);
            }
        }
        System.out.printf("R_theta_t Right=\n%s\n", FormatArray.toString(R_theta_tRight, "%.4e"));

        System.out.printf("K*K^T=\n%s\n", 
            FormatArray.toString(MatrixUtil.multiply(kScaled, kInv), "%.4e"));

        double[][] H1Left = MatrixUtil.multiply(kScaled, R_theta_t);
        H1Left = MatrixUtil.multiply(H1Left, kInv);
        
        double[][] H1Right = MatrixUtil.multiply(kScaled, R_theta_tRight);
        H1Right = MatrixUtil.multiply(H1Right, kInv);
        
        System.out.printf("H1_left=\n%s\n", FormatArray.toString(H1Left, "%.4e"));
        
        System.out.printf("H1_right=\n%s\n", FormatArray.toString(H1Right, "%.4e"));
     
        System.out.printf("H1_left * eLeft =\n%s\n", 
            FormatArray.toString(
                MatrixUtil.multiplyMatrixByColumnVector(H1Left, fEpipoles[0]),
            "%.4e"));
        
        System.out.println("-----------");
    }

    private void printMallonWhelanRectification(
        DenseMatrix fm, double[][] fEpipoles,
        DenseMatrix x1M, DenseMatrix x2M, 
        double focalLength, double xC, double yC, double scale1) {
        
        /*
        Projective Rectification from the
Fundamental Matrix
John Mallon ∗ Paul F. Whelan
        */
    }
}
