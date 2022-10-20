package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.transform.Rectification.RectifiedPoints;
import algorithms.matrix.MatrixUtil;

import algorithms.misc.MiscDebug;
import algorithms.util.FormatArray;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertNotNull;

import no.uib.cipr.matrix.DenseMatrix;

/**
 *
 * @author nichole
 */
public class RectificationTest extends TestCase {
    
    public RectificationTest() {
    }

    public void testRectification0() throws Exception {
        System.out.println("testRectification oldhouse");

        ImageExt im1 = MASKSData.getOldhouse2A2000000();
        ImageExt im2 = MASKSData.getOldhouse2A2000084();

        double[][] x1 = MASKSData.getXOldhouse2A2000000();
        double[][] x2 = MASKSData.getXOldhouse2A2000084();

        double[][] fm = MASKSData.getFMOldhouse2A2_0_84();

        // these epipoles are both outside the image boundaries so the
        // MASKS rectify algorithm can be used on these data.
        double[][] e1e2 = EpipolarTransformer.calculateEpipoles(new DenseMatrix(fm));
        System.out.printf("oldhouse e1e2=\n%s\n", FormatArray.toString(e1e2, "%.3e"));

        EpipolarTransformer tr = new EpipolarTransformer();
        double[][] fm2 = tr.calculateEpipolarProjection2(x1, x2, false);

        RectifiedPoints rP = Rectification.rectify(fm, x1, x2, im1.getWidth()/2, im1.getHeight()/2);

        int i;
        double[][] x1R = rP.getX1();
        double[][] x2R = rP.getX2();
        int n = x1R[0].length;
        System.out.println("oldhouse keypoints rectified, y's should be similar");
        for (i = 0; i < n; ++i) {
            System.out.printf("%d) (%.1f, %.1f, %.1f)  (%.1f, %.1f, %.1f)\n",
                i, x1R[0][i], x1R[1][i], x1R[2][i],
                x2R[0][i], x2R[1][i], x2R[2][i]);
        }

        Rectification.RectifiedImage imgW1 = Rectification.hWarp(im1.copyToGreyscale2(), rP.getH1());

        Rectification.RectifiedImage imgW2 = Rectification.hWarp(im2.copyToGreyscale2(), rP.getH2());

        MiscDebug.writeImage(imgW1.copyToGreyscale2(), "_rect_oldhouse_00_");
        MiscDebug.writeImage(imgW2.copyToGreyscale2(), "_rect_oldhouse_84_");

    }

    public void testRectification1() throws Exception {
        System.out.println("testRectification Zhang");

        int idx1 = 1;
        int idx2 = 2;//5;

        double[][] k1 = Zhang98Data.getIntrinsicCameraMatrix();
        //x1, x2 size is 3 X 256
        double[][] x1 = Zhang98Data.getObservedFeaturesInImage(idx1);
        double[][] x2 = Zhang98Data.getObservedFeaturesInImage(idx2);

        boolean useR2R4 = true;
        double[] rd = Zhang98Data.getRadialDistortionR2R4();
        double[][] x1c = Camera.pixelToCameraCoordinates(x1, k1, null /*rd*/, useR2R4);
        double[][] x2c = Camera.pixelToCameraCoordinates(x2, k1, null /*rd*/, useR2R4);

        ImageExt im1 = Zhang98Data.getImage(idx1);
        ImageExt im2 = Zhang98Data.getImage(idx2);

        double[][] x1n = MatrixUtil.copy(x1);
        double[][] t1 = EpipolarNormalizationHelper.unitStandardNormalize(x1n);
        double[][] x2n = MatrixUtil.copy(x2);
        double[][] t2 = EpipolarNormalizationHelper.unitStandardNormalize(x2n);

        EpipolarTransformer tr = new EpipolarTransformer();
        double[][] fmN = tr.calculateEpipolarProjection2(x1n, x2n, false);

        //de-normalize to work with image coordinates
        double[][] fm = MatrixUtil.copy(fmN);
        EpipolarNormalizationHelper.denormalizeFM(fm, t1, t2);

        // when using idx1 = 1 and idx2 = 5, BOTH epipoles are within image bounds so we cannot use the
        //  MASKS uncalibrated method on them.
        // TODO: implement the non-linear polar rectification algorithm of Pollefys

        //NOTE: if idx1=1 and idx2=2, can see from the epipoles that the rectifications will
        // rotate the image around z-axis by a large amount

        // for idx1 = 1 and idx2 = 2, there must be another problem

        double[][] e1e2 = EpipolarTransformer.calculateEpipoles(new DenseMatrix(fm));

        System.out.printf("e1e2=\n%s\n", FormatArray.toString(e1e2, "%.3e"));

        RectifiedPoints rP = Rectification.rectify(fm, x1, x2, im1.getWidth()/2,
            im1.getHeight()/2);

        int i;
        double[][] x1R = rP.getX1();
        double[][] x2R = rP.getX2();
        int n = x1R[0].length;
        System.out.println("zhang checkerboard keypoints rectified, y's should be similar");
        for (i = 0; i < n; ++i) {
            System.out.printf("%d) (%.1f, %.1f, %.1f)  (%.1f, %.1f, %.1f)\n",
                    i, x1R[0][i], x1R[1][i], x1R[2][i],
                    x2R[0][i], x2R[1][i], x2R[2][i]);
        }

        Rectification.RectifiedImage imgW1 = Rectification.hWarp(im1.copyToGreyscale2(), rP.getH1());

        Rectification.RectifiedImage imgW2 = Rectification.hWarp(im2.copyToGreyscale2(), rP.getH2());

        MiscDebug.writeImage(imgW1.copyToGreyscale2(), "_rect_zhang98_01_");
        MiscDebug.writeImage(imgW2.copyToGreyscale2(), "_rect_zhang98_02_");

        Rectification.RectifiedImage[] imgR2 = Rectification.hWarpUnionSize(
                im1.copyToGreyscale2(), im2.copyToGreyscale2(), rP.getH1(), rP.getH2());
        MiscDebug.writeImage(imgR2[0].copyToGreyscale2(), "_union_rect_zhang98_01_");
        MiscDebug.writeImage(imgR2[1].copyToGreyscale2(), "_union_rect_zhang98_02_");

    }
    
}
