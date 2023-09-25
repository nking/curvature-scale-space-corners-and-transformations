package algorithms.compGeometry;

import algorithms.matrix.MatrixUtil;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

public class Misc {

    /**
     * given a set of 3D points (format [4Xn] where n is the number of points),
     * calculate the plane normal to the points and the distance from the
     * origin to the plane.
     * <pre>
     * extracted from  homography2Motion.m from the example code from the book
     * Chapter 5, "An introduction to 3-D Vision"
     * by Y. Ma, S. Soatto, J. Kosecka, S. Sastry (MASKS)
     * Code distributed free for non-commercial use
     * Copyright (c) MASKS, 2003
     * </pre>
     * @param xWCS
     * @return an array holding the plane normal to the points xWCS and the distance of the
     * plane from the origin.
     * the first 3 elements are the plane perpendicular to the points.
     * the last element is the distance of the place from the origin.
     */
    public static double[] orthogonalPlaneAndDistance(double[][] xWCS) throws NotConvergedException {
        if (xWCS.length != 4) {
            throw new IllegalArgumentException("xWCS length must be 4");
        }
        double[][] p = MatrixUtil.copy(xWCS);
        int i, j;
        for (i = 0; i < p[0].length; ++i) {
            for (j = 0; j < 4; ++j) {
                p[j][i] /= p[3][i];
            }
        }
        // [nX4]
        p = MatrixUtil.transpose(p);

        SVD svd = SVD.factorize(new DenseMatrix(p));

        // V is [4X4]
        // the first 3 elements are the plane perpendicular to the points.
        // the last element is the distance of the place from the origin.
        double[] orth = new double[4];
        for (i = 0; i < 4; ++i) {
            orth[i] = svd.getVt().get(3, i);
        }
        return orth;
    }
}
