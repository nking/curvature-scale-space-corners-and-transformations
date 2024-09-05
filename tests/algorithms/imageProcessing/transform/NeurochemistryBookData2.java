package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.*;
import algorithms.imageProcessing.features.orb.ORB;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.matching.ORBMatcher;
import algorithms.matrix.MatrixUtil;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.FormatArray;
import algorithms.util.ResourceFinder;
import no.uib.cipr.matrix.NotConvergedException;

import java.io.File;
import java.io.IOException;

/**
<pre>
 4 images in this set.

 * images are 3024 X 4032
    the x-y-z directions are considered to be relative to the book cover bottom center,
        and from our perspective.
        x-axis is toward the right, y-axis is up and z-axis is out of the picture, towards
        the centered camera.

    all pictures are from the same camera, but will use the convention left camera for the
        image taken to the left of the book center, right camera for the image taken
        from the right (and centered in the image).
        the rightUp picture has more rotation in it for the unit tests.

    the centered camera is 41.25 cm from the book center.
    when camera is to the left of the book, it is 45 cm from the book center
        and 18 cm from the centered camera location.
    same for the camera when to the left of the book.
    3 images have 0 rotation about x-axis and the rightUp rotation should be about 6 or 7 degrees around
    the x-axis.
    the right triangle formed by the book's face cover bottom center as 1 point, the centered camera
        as another point and the left or right cameras as the 3rd point, has angles
        0.41146 radians (= 23.575 degrees), 1.1593400833811942 radians (= 66.425 degrees), and
        1.571 radians (= 90 degrees).
        1 degree is 0.0175 radians.
        errors in those angle estimates is possibly a few degrees which is 0.05 radians, and the
        distances might have measurement errors of 0.2 cm each.

    image exif info:
        aperture |  shutter |          |  focal  |  subject
        value    |  speed    | f-number  | length  |  distance
        1.58      1/120       1.73        4.38 mm    (0.371, 0.456, 0.391, 0.547)

        3024 X 4032

        expecting
             focalLength ~ 1604 pixels = 2.245 mm
             no skew
             xc=1512
             yc=2016
             little to no radial distortion (if was present, it is already removed)
             rotation between images = 23.575 degrees
             translation between images = 18 cm

          other information:
            pixel width = 1.4e-3mm
            FOV = 77 degrees = 1.344 radians
 </pre>
 * @author nichole
 */
public class NeurochemistryBookData2 {
 
    public static final int nFeatures = 0;
    public static final int mImages = 6;

    public static String DIR = "nc_book";
    protected final static String sep = System.getProperty("file.separator");


    /**
     * get the features in world scene coordinates.
     * @return double array of
     */
    public static double[][] getFeatureWCS()  {

        // book width is 22.2 cm.
        // height is 28.3 cm.
        // origin is the centered camera center, so z coordinate is 41.25 cm.

        double[][] wcs = new double[nFeatures][];
        /*wcs[0] = new double[]{, 41.25};
        wcs[1] = new double[]{, 41.25};
        wcs[2] = new double[]{, 41.25};
        wcs[3] = new double[]{, 41.25};
        wcs[4] = new double[]{, 41.25};
        wcs[5] = new double[]{, 41.25};
        wcs[6] = new double[]{, 41.25};
        wcs[7] = new double[]{, 41.25};*/
        // mean=
        // stdev=

        //double[] mn = new double[3];
        //double[] sd = new double[3];
        //double[][] wcs2 = Standardization.standardUnitNormalization(
        //    wcs, mn, sd);
        // mean=-1.4625, 0.4625, 41.5
        // stdev=8.82, 11.303, 0.

        wcs = MatrixUtil.transpose(wcs);        
        return wcs;
    }

    // size is 3 X (nFeatures*mImages)
    public static double[][] getObservedFeaturesInImage(int idx) {
        if (idx < 0 || idx > mImages-1) {
            throw new IllegalArgumentException("idx must be 0 through 3, inclusive");
        }

        double[][] uv = new double[nFeatures][];
        switch(idx) {
            case 0:
                //uv[0] = new double[]{608, 530, 1};
                break;
            case 1:
                //uv[0] = new double[]{678, 718, 1};
                break;
            case 2:
                //uv[0] = new double[]{678, 718, 1};
                break;
            case 3:
                //uv[0] = new double[]{678, 718, 1};
                break;
            case 4:
                //uv[0] = new double[]{678, 718, 1};
                break;
            case 5:
                //uv[0] = new double[]{678, 718, 1};
                break;
            default:
                //uv[0] = new double[]{744, 806, 1};
                break;
        }        
        
        return MatrixUtil.transpose(uv);
    }

    public static double[][] getRotation(int idx) {
        if (idx < 0 || idx > mImages-1) {
            throw new IllegalArgumentException("idx must be 0 through 3, inclusive");
        }
        double[] thetas;
        switch (idx) {
            // cc along y is + for right-handed system
            case 0: // centered
                thetas = new double[]{0, 0, 0};
                break;
            case 1: // to the right
                thetas = new double[]{0, 23.575*Math.PI/180., 0};
                break;
            case 2: // to the left
                thetas = new double[]{0, -23.575*Math.PI/180., 0};
                break;
            case 3: // to the right and rotated up too
                thetas = new double[]{6.5*Math.PI/180., 23.575*Math.PI/180., 0};
                break;
            case 4: // centered, facing the book. the z-axis tilt was about 17 degrees.
                thetas = new double[]{6.5*Math.PI/180., 23.575*Math.PI/180., 17*Math.PI/180.};
                break;
            case 5: // to the right of book, rotated up along x-axis too, and also about the z-axis
                // did not have a way to measure the x-axis tilt.  the z-axis tilt was about 17 degrees.
                thetas = new double[]{6.5*Math.PI/180., 23.575*Math.PI/180., 17*Math.PI/180.};
                break;
            default:
                throw new IllegalArgumentException("idx out of range");
        }
        return Rotation.createRotationZYX(thetas);
    }
    
    /**
     * translation is in WCS units of cm.
     * @param idx
     * @return 
     */
    public static double[] getTranslation(int idx) {
        if (idx < 0 || idx > mImages-1) {
            throw new IllegalArgumentException("idx must be 0 through 3, inclusive");
        }
        // offset of image center of book and center of image
        //    gives the offset in Y of the WCS reference frame origin to the image center (==camera origin).
        // that's (2016-1678)/81=4.2
        switch (idx) {
            case 0:
                return new double[]{0, 0, 0};
            case 1:
                return new double[]{+18, 0, 0};
            case 2:
                return new double[]{-18, 0, 0};
            case 3:
                return new double[]{+18, 0, 0};
            case 4:
                return new double[]{0, 0, 0};
            case 5:
                return new double[]{+18, 0, 0};
            default:
                //-7.376e+00, -4.029e+00, 4.009e+01 
                return new double[]{18, -4.2, 41.5};
        }        
    }
    
    public static double[][] getIntrinsicCameraMatrix() {
        double[][] intr = MatrixUtil.zeros(3, 3);
        intr[0][0] = 1600; //2.189e+03
        intr[1][1] = 1600; //2.886e+03
        intr[0][2] = 1512; //1500
        intr[1][2] = 2016; //2000
        intr[2][2] = 1;
        return intr;
    }
    
    /**
     * the radial distortion coefficients are for the distortion polynomial k1*r^2 + k2*r^4.
     * @return 
     */
    public static double[] getRadialDistortionR2R4() {
        return new double[]{0, 0};
    }

    // 6 images
    public static ImageExt getImage(int idx) throws IOException {
        if (idx < 0 || idx > 5) {
            throw new IllegalArgumentException("idx must be 1 through 5, inclusive");
        }
        ++idx;

        String path = ResourceFinder.findTestResourcesDirectory() + sep
                + DIR
                + sep + "nc_book_0" + Integer.toString(idx) + ".png";

        File f = new File(path);
        if (!f.exists()) {
            throw new IOException("could not find file at " + path);
        }

        ImageExt image = ImageIOHelper.readImageExt(path);

        return image;
    }

    public static void runCorresMaker() throws IOException {

        boolean useToleranceAsStatFactor = true;
        boolean recalcIterations = false;// possibly faster if set to true
        double tolerance = 2;
        ErrorType errorType = ErrorType.SAMPSONS;

        int np = 500;

        int f = 4;
        int w2 = 3024/f;  // 378
        int h2 = 4032/f; // 540

        ImageProcessor ip = new ImageProcessor();
        GreyscaleImage img1, img2;
        ORB orb1, orb2;
        ORB.Descriptors d1, d2;
        double[][] xKP1, xKP2;
        int nKP1, nKP2;
        double[][] xKP1n, xKP2n;
        double[][] t1, t2;

        int i, j, ii;

        for (i = 0; i < mImages; ++i) {
            img1 = getImage(i).copyToGreyscale2();
            img1 = ip.downSample(img1, w2, h2, 0, 255);

            //TODO: scale the resulting coordinate lists

            orb1 = new ORB(img1, np);
            orb1.detectAndExtract();
            d1 = orb1.getAllDescriptors();
            xKP1 = orb1.getAllKeyPointsHomogenous();
            nKP1 = xKP1[0].length;

            int x1, x2, y1, y2;
            Image tmp1 = img1.copyToColorGreyscale();
            for (ii = 0; ii < nKP1; ++ii) {
                y1 = (int)Math.round(xKP1[1][ii]);
                x1 = (int)Math.round(xKP1[0][ii]);
                ImageIOHelper.addPointToImage(x1, y1, tmp1, 2, 255, 0, 0);
            }
            xKP1n = MatrixUtil.copy(xKP1);
            t1 = EpipolarNormalizationHelper.unitStandardNormalize(xKP1n);

            for (j = 1000/*i+1*/; j < mImages; ++j) {

                img2 = getImage(j).copyToGreyscale2();
                orb2 = new ORB(img2, np);
                orb2.detectAndExtract();
                d2 = orb2.getAllDescriptors();
                xKP2 = orb2.getAllKeyPointsHomogenous();
                nKP2 = xKP2[0].length;
                xKP2n = MatrixUtil.copy(xKP2);
                t2 = EpipolarNormalizationHelper.unitStandardNormalize(xKP2n);

                ORBMatcher.FitAndCorres fitAndCorres = ORBMatcher.matchDescriptors(
                        d1, d2, xKP1n, xKP2n, useToleranceAsStatFactor,
                            tolerance, errorType, recalcIterations, false);

                int idx1, idx2;
                Image tmp2 = img2.copyToColorGreyscale();
                CorrespondencePlotter plotter = new CorrespondencePlotter(tmp1, tmp2);
                if (fitAndCorres.mIF != null) {
                    System.out.printf("%d, %d) #keypoints=(%d,%d) #matched=%d, #after filter=%d\n", i, j,
                            xKP1n[0].length, xKP2n[0].length, fitAndCorres.mI.length,
                            fitAndCorres.mIF.length);
                    for (ii = 0; ii < fitAndCorres.mIF.length; ++ii) {
                        idx1 = fitAndCorres.mIF[ii][0];
                        idx2 = fitAndCorres.mIF[ii][1];
                        x1 = (int)xKP1[0][idx1];
                        y1 = (int)xKP1[1][idx1];
                        x2 = (int)xKP2[0][idx2];
                        y2 = (int)xKP2[1][idx2];
                        plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 1);
                    }
                } else {
                    System.out.printf("%d, %d)  #keypoints=(%d,%d)  #matched=%d, but RANSAC fits failed\n", i, j,
                            xKP1n[0].length, xKP2n[0].length, fitAndCorres.mI.length);
                    for (ii = 0; ii < fitAndCorres.mI.length; ++ii) {
                        idx1 = fitAndCorres.mI[ii][0];
                        idx2 = fitAndCorres.mI[ii][1];
                        x1 = (int)xKP1[0][idx1];
                        y1 = (int)xKP1[1][idx1];
                        x2 = (int)xKP2[0][idx2];
                        y2 = (int)xKP2[1][idx2];
                        plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 1);
                    }
                }
                plotter.writeImage("_corres_orb_nc_" + Integer.toString(i) + "_" + Integer.toString(j) + "_");
            }
        }
    }
}
