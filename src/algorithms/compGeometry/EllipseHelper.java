package algorithms.compGeometry;

import algorithms.misc.MiscMath;
import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import java.io.IOException;
import java.util.Set;
import java.util.logging.Logger;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.EigenOps;
import org.ejml.simple.*;

/**
 * 
 * @author nichole
 */
public class EllipseHelper {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     <pre>
     adapted from:
        "Fitting of Circles and Ellipses Least Squares Solution" by Gander, 
        Golub, &amp; Strebel, 1994 is available from anonymous ftp.inf.ethz.ch 
        in directory doc/tech-reports/2xx
        The files begin with 217.* and the list includes matlab source code.
        files algellipse.m and ellipse_params.m are adapted here.
     
       currently uses the algebraic, least square ellipse fit
       
     The parameters returned can be used as:
     x(t) = xCenter + aParam*cos(alpha)*cos(t) − bParam*sin(alpha)*sin(t)
     y(t) = yCenter + aParam*sin(alpha)*cos(t) + bParam*cos(alpha)*sin(t)
     where t range is 0 &lt; == t &lt; == 2π
      </pre>
     * @param xyPoints
     * @return new double[]{xCenter, yCenter, aParam, bParam, alpha};
     */
    public double[] fitEllipseToPoints(PairFloatArray xyPoints) {
        
        return fitEllipseToPointsWithAlgLSQ(xyPoints);
    }
    
    /**
     <pre>
     adapted from:
        "Fitting of Circles and Ellipses Least Squares Solution" by Gander, 
        Golub, &amp; Strebel, 1994 is available from anonymous ftp.inf.ethz.ch 
        in directory doc/tech-reports/2xx
        The files begin with 217.* and the list includes matlab source code.
        files algellipse.m and ellipse_params.m are adapted here.
     
       currently uses the algebraic, least square ellipse fit
       
     The parameters returned can be used as:
     x(t) = xCenter + aParam*cos(alpha)*cos(t) − bParam*sin(alpha)*sin(t)
     y(t) = yCenter + aParam*sin(alpha)*cos(t) + bParam*cos(alpha)*sin(t)
     where t range is 0 &lt; == t &lt; == 2π
      </pre>
     * @param xyPoints
     * @return new double[]{xCenter, yCenter, aParam, bParam, alpha};
     */
    public double[] fitEllipseToPoints(PairIntArray xyPoints) {
        
        int nPoints = xyPoints.getN();
        
        //x*x, x*y, y*y, x, y, 1
        SimpleMatrix a = new SimpleMatrix(nPoints, 6);
        for (int row = 0; row < nPoints; row++) {
            float x = xyPoints.getX(row);
            float y = xyPoints.getY(row);
            a.setRow(row, 0, x*x, x*y, y*y, x, y, 1);
        }
       
        try {
            
            return fitEllipseToPointsWithAlgLSQ(a);
                        
        } catch(RuntimeException t) {
            log.warning(t.getMessage());
        }
        
        return null;
    }
    
    /**
     <pre>
     adapted from:
        "Fitting of Circles and Ellipses Least Squares Solution" by Gander, 
        Golub, &amp; Strebel, 1994 is available from anonymous ftp.inf.ethz.ch 
        in directory doc/tech-reports/2xx
        The files begin with 217.* and the list includes matlab source code.
        files algellipse.m and ellipse_params.m are adapted here.
     
       currently uses the algebraic, least square ellipse fit
       
     The parameters returned can be used as:
     x(t) = xCenter + aParam*cos(alpha)*cos(t) − bParam*sin(alpha)*sin(t)
     y(t) = yCenter + aParam*sin(alpha)*cos(t) + bParam*cos(alpha)*sin(t)
     where t range is 0 &lt; == t &lt; == 2π
      </pre>
     * @param xyPoints
     * @return new double[]{xCenter, yCenter, aParam, bParam, alpha};
     */
    public double[] fitEllipseToPoints(Set<PairInt> xyPoints) {
        
        int nPoints = xyPoints.size();
        
        //x*x, x*y, y*y, x, y, 1
        SimpleMatrix a = new SimpleMatrix(nPoints, 6);
        int row = 0;
        for (PairInt p : xyPoints) {
            float x = p.getX();
            float y = p.getY();
            a.setRow(row, 0, x*x, x*y, y*y, x, y, 1);
            row++;
        }
       
        try {
            
            return fitEllipseToPointsWithAlgLSQ(a);
                        
        } catch(RuntimeException t) {
            log.warning(t.getMessage());
        }
        
        return null;
    }
    
    /**
     <pre>
     adapted from:
        "Fitting of Circles and Ellipses Least Squares Solution" by Gander, 
        Golub, &amp; Strebel, 1994 is available from anonymous ftp.inf.ethz.ch 
        in directory doc/tech-reports/2xx
        The files begin with 217.* and the list includes matlab source code.
        files algellipse.m and ellipse_params.m are adapted here.
     
       algebraic, least square ellipse fit
      </pre>
     * @param xyPoints
     * @return new double[]{xCenter, yCenter, aParam, bParam, alpha};
     */
    protected double[] fitEllipseToPointsWithAlgLSQ(PairFloatArray xyPoints) {
        
        int nPoints = xyPoints.getN();
        
        //x*x, x*y, y*y, x, y, 1
        SimpleMatrix a = new SimpleMatrix(nPoints, 6);
        for (int row = 0; row < nPoints; row++) {
            float x = xyPoints.getX(row);
            float y = xyPoints.getY(row);
            a.setRow(row, 0, x*x, x*y, y*y, x, y, 1);
        }
       
        try {
            
            return fitEllipseToPointsWithAlgLSQ(a);
                        
        } catch(RuntimeException t) {
            log.warning(t.getMessage());
        }
        
        return null;
    }
    
    /**
     <pre>
     adapted from:
        "Fitting of Circles and Ellipses Least Squares Solution" by Gander, 
        Golub, &amp; Strebel, 1994 is available from anonymous ftp.inf.ethz.ch 
        in directory doc/tech-reports/2xx
        The files begin with 217.* and the list includes matlab source code.
        files algellipse.m and ellipse_params.m are adapted here.
     
       algebraic, least square ellipse fit
      
     NOTE: SimpleMatrix.solve(...) method may throw:
          org.ejml.factory.SingularMatrixException: Solution contains 
          uncountable numbers
     </pre> 
     * @param a
     * @return new double[]{xCenter, yCenter, aParam, bParam, alpha};
     */
    @SuppressWarnings({"unchecked"})
    protected double[] fitEllipseToPointsWithAlgLSQ(SimpleMatrix a) {
        
        int nPoints = a.numRows();
        
        /*
        adapted from:
        "Fitting of Circles and Ellipses Least Squares Solution" by Gander, 
        Golub, & Strebel, 1994 is available from anonymous ftp.inf.ethz.ch 
        in directory doc/tech-reports/2xx
        The files begin with 217.* and the list includes matlab source code.

        algebraic, least square ellipse fit
        
        this is an adaptation to the code they provide on the ftp site
        referenced in their paper, algellipse.m
        
        input:
            X: given points Pi = [X(i,1), X(i,2)]
            W: weight W(i) for the i-th equation
        
        fits an ellipse by minimizing the "algebraic distance"
        in the least squares sense x^T A x + b^T x + c = 0
        weighting the i-th data by W(i)
     
        solves for z, a, b, alpha parameters of the ellipse
        
        x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0
        
        subject to constraints:  
            B^2 - 4*C < 0
        
            D^2 + E^2
            ---   --- - F &gt; 0
             4    4*C
        */
   
        SimpleMatrix s = a.transpose().mult(a);
       
        //T = S(1:3,1:3) - S(1:3,4:6)*(S(4:6,4:6)'\S(4:6,1:3));
        //      0:2, 0:2     0:2,3:5     3:5,3:5     3:5,0:2
        
        // S(4:6,4:6) is a reduction of the matrix to S(rows 4 to 6, columns 4 to 6)
        // ' is transpose
        // \ is Matrix left division where x = A\B is the solution to the 
        //          equation xA = B. 
        //          Matrices A and B must have the same number of rows.
        
        SimpleMatrix matrixA = new SimpleMatrix(3, 3);
        for (int row = 3; row <= 5; row++) {
            for (int col = 3; col <= 5; col++) {
                matrixA.set(row-3, col-3, s.get(row, col));
            }
        }
        matrixA = matrixA.transpose();
        
        SimpleMatrix matrixB = new SimpleMatrix(3, 3);
        for (int row = 3; row <= 5; row++) {
            for (int col = 0; col <= 2; col++) {
                matrixB.set(row-3, col, s.get(row, col));
            }
        }
        
        //xA = B.   x is [3x3]
        SimpleMatrix x = matrixA.solve(matrixB);
        
        //T = S(1:3,1:3) - S(1:3,4:6)*(x);
        //      0:2, 0:2     0:2,3:5   
        
        SimpleMatrix matrixC = new SimpleMatrix(3, 3);
        for (int row = 0; row <= 2; row++) {
            for (int col = 3; col <= 5; col++) {
                matrixC.set(row, col-3, s.get(row, col));
            }
        }
        
        matrixC = matrixC.mult(x);
        
        SimpleMatrix matrixD = new SimpleMatrix(3, 3);
        for (int row = 0; row <= 2; row++) {
            for (int col = 0; col <= 2; col++) {
                matrixD.set(row, col, s.get(row, col));
            }
        }
        
        //T = matrixD - matrixC
        SimpleMatrix matrixT = matrixD.minus(matrixC);
        
        //T = diag([1,2,1])*T;
        SimpleMatrix diag = new SimpleMatrix(3, 3);
        diag.set(0, 0, 1);
        diag.set(1, 1, 2);
        diag.set(2, 2, 1);
        
        matrixT = diag.mult(matrixT);
        
        /*
        [V, D] = eig(T);
        returns two optional outputs for any of the previous input syntaxes. 
        D is a diagonal matrix containing the eigenvalues. 
        V is a matrix whose columns are the corresponding right eigenvectors.
        */
        SimpleEVD<SimpleMatrix> evd = matrixT.eig();
        
        //3X3
        DenseMatrix64F v = EigenOps.createMatrixV(evd.getEVD());
        //v.print();
        SimpleMatrix vMatrix = new SimpleMatrix(v);
      
        double emin = 0;
        int kmin = -1;
        for (int k = 0; k < 3; k++) {
            double aa = v.get(0, k);
            double bb = v.get(1, k);
            double cc = v.get(2, k);
            double i0 = aa + cc;
            double i1 = (aa * cc) - ((bb*bb)/4.);
            if (i1 <= 0) {
                // not an ellipse
            } else {
                double val = (i0*i0 - 4*i1)/(i0*i0 - 2*i1);
                if ((emin == 0) || (val < emin)) {
                    emin = val;
                    kmin = k;
                }
            }
        }
        if (kmin == -1) {
            // not an ellipse, or my port to java is wrong
            log.severe("not an ellipse, or error in my port to java?  need more tests for this");
            return null;
        }
        
        SimpleMatrix y1 = vMatrix.extractVector(false, kmin);
        
        //-(S(4:6,4:6)')
        SimpleMatrix numer = new SimpleMatrix(3, 3);
        for (int row = 3; row <= 5; row++) {
            for (int col = 3; col <= 5; col++) {
                numer.set(row - 3, col - 3, s.get(row, col));
            }
        }
        numer = numer.transpose();
        
        //(S(1:3,4:6)'*y1);
        SimpleMatrix denom = new SimpleMatrix(3, 3);
        for (int row = 0; row <= 2; row++) {
            for (int col = 3; col <= 5; col++) {
                denom.set(row, col - 3, s.get(row, col));
            }
        }
        denom = denom.transpose();
        denom = denom.mult(y1);
        
        // y2 = -(S(4:6,4:6)')\(S(1:3,4:6)'*y1);
        SimpleMatrix y2 = numer.solve(denom);
        y2 = y2.scale(-1);
        
        // u  = [y1; y2]; where ';' is matlab notation for row separator
        // 
        int nRows = y1.numRows() + y2.numRows();
        SimpleMatrix uMatrix = new SimpleMatrix(nRows, y1.numCols());
        
        int row = 0;
        for (int j = 0; j < y1.numRows(); j++) {
            for (int col = 0; col < y1.numCols(); col++) {
                uMatrix.set(row, col, y1.get(j, col));
            }
            row++;
        }
        for (int j = 0; j < y2.numRows(); j++) {
            for (int col = 0; col < y2.numCols(); col++) {
                uMatrix.set(row, col, y2.get(j, col));
            }
            row++;
        }
        
        //get the ellipse parameters from algebraic equation 
        //    u(0)x^2 + u(1)xy + u(2)y^2 + u(3)x + u(4)y + u(5) = 0.
        
        return extractEllipseParams(uMatrix);
    }
    
    /**
    <pre>
    adapted from:
        "Fitting of Circles and Ellipses Least Squares Solution" by Gander, 
        Golub, &amp; Strebel, 1994 is available from anonymous ftp.inf.ethz.ch 
        in directory doc/tech-reports/2xx
        file: ellipse_params.m in compressed archive for 217.*
       
     The parameters returned can be used as:
     x(t) = xCenter + aParam*cos(alpha)*cos(t) − bParam*sin(alpha)*sin(t)
     y(t) = yCenter + aParam*sin(alpha)*cos(t) + bParam*cos(alpha)*sin(t)
     where t range is 0 &lt; == t &lt; == 2π
     </pre>  
     * @param u
     * @return new double[]{xCenter, yCenter, aParam, bParam, alpha};
     */
    @SuppressWarnings({"unchecked"})
    private double[] extractEllipseParams(SimpleMatrix u) {
        
        if (u == null) {
            throw new IllegalArgumentException("u cannot be null");
        }
        if (u.numCols() != 1) {
            throw new IllegalArgumentException("u must have 1 column only");
        }
        if (u.numRows() != 6) {
            throw new IllegalArgumentException("u must have 6 rows");
        }
        
        SimpleMatrix a = new SimpleMatrix(2, 2);
        a.setRow(0, 0, u.get(0, 0), u.get(1, 0)/2.);
        a.setRow(1, 0, u.get(1, 0)/2., u.get(2, 0));
        
        SimpleMatrix bb = new SimpleMatrix(2, 1);
        bb.set(0, 0, u.get(3, 0));
        bb.set(1, 0, u.get(4, 0));
        
        SimpleMatrix c = new SimpleMatrix(1, 1);
        c.set(0, 0, u.get(5, 0));
        
        SimpleEVD<SimpleMatrix> qd = a.eig();
        
        DenseMatrix64F q = EigenOps.createMatrixV(qd.getEVD());
        //q.print();
        
        DenseMatrix64F d = EigenOps.createMatrixD(qd.getEVD());
        //d.print();
        
        double det = d.get(0, 0) * d.get(1, 1);
        
        if (det <= 0) {
            log.severe("error while extracting ellipse parameters");
            return null;
        }
        
        SimpleMatrix qMatrix = new SimpleMatrix(q);
        
        SimpleMatrix bs = (qMatrix.transpose()).mult(bb);
        
        double alpha = Math.atan2(q.get(1, 0), q.get(0, 0));
        
        SimpleMatrix numer = new SimpleMatrix(d);
        numer = numer.scale(-2.);
        SimpleMatrix zs = numer.solve(bs);
        
        SimpleMatrix z = qMatrix.mult(zs);
        
        SimpleMatrix h = bs.transpose().mult(zs);
        h = h.divide(-2.);
        h = h.minus(c);

        double aParam = Math.sqrt(h.get(0, 0)/d.get(0, 0));
        double bParam = Math.sqrt(h.get(0, 0)/d.get(1, 1));
        
        double z0 = z.get(0, 0);
        double z1 = z.get(1, 0);
        
        /*
        TODO: check this
        making a correction for aParam < bParam
        */
        if (aParam < bParam) {
            double swap = aParam;
            aParam = bParam;
            bParam = swap;
            alpha += (Math.PI/2.);
        }
        
        return new double[]{z0, z1, aParam, bParam, alpha};
    }
    
    /**
     * compute the statistics of the residuals of (xP,yP) from an ellipse 
     * described by the given parameters.
    
     * @param xP
     * @param yP
     * @param xCenterParam
     * @param yCenterParam
     * @param aParam
     * @param bParam
     * @param alphaParam
     * @return new double[avgResid, stDevResid]
     */
    @SuppressWarnings({"unchecked"})
    public double[] calculateEllipseResidualStats(float[] xP, float[] yP, 
        float xCenterParam, float yCenterParam, float aParam, float bParam, 
        float alphaParam) {
        
        if (bParam >= aParam) {
            throw new IllegalArgumentException("for an ellipse, a < b");
        }
        
        double[][] translated = new double[xP.length][];
        for (int row = 0; row < xP.length; row++) {
            translated[row] = new double[]{xP[row] - xCenterParam, 
                yP[row] - yCenterParam};
        }
  
        double sine = Math.sin(alphaParam);
        double cosine = Math.cos(alphaParam);
        
        SimpleMatrix q = new SimpleMatrix(2, 2);
        q.setRow(0, 0, cosine, -1*sine);
        q.setRow(1, 0, sine, cosine);
  
        SimpleMatrix transformed = (new SimpleMatrix(translated)).mult(q);
                
        // c^2 = a^2 - b^2 and focii are at (+-c, 0)
        double c = Math.sqrt(aParam*aParam - bParam*bParam);
        
        double[] resid = new double[xP.length];
                
        for (int row = 0; row < xP.length; row++) {
            
            double xPoint = transformed.get(row, 0);
            double yPoint = transformed.get(row, 1);
            
            /*
            compute the residual as the difference between 
            the point to the 2 foci (+-c, 0) and the value 2*a.
            2*a = dist to F1 + dist to F2 for any point on the ellipse
            */
            
            double diffX1 = xPoint - c;
            double diffX0 = xPoint - -c;
                        
            double distF1 = Math.sqrt(diffX1*diffX1 + yPoint*yPoint);
            
            double distF0 = Math.sqrt(diffX0*diffX0 + yPoint*yPoint);
            
            double residSum = 2.*aParam - (distF1 + distF0);
                        
            resid[row] = residSum;
        }
        
        double[] avgStdDev = MiscMath.getAvgAndStDev(resid);
                
        return new double[]{avgStdDev[0], avgStdDev[1]};
    }
    
    public void plotEllipseAndPoints(Set<PairInt> points, 
        float xCenterParam, float yCenterParam, float aParam, float bParam, 
        float alphaParam, int xMin, int xMax, int yMin, int yMax,
        int plotNumber, String plotLabel) throws IOException {
        
        float[] ellipseX = new float[360];
        float[] ellipseY = new float[360];
        double ca = Math.cos(alphaParam);
        double sa = Math.sin(alphaParam);
            
        for (int angle = 0; angle < 360; angle++) {
            double theta = ((double)angle) * Math.PI/180.;
            double ct = Math.cos(theta);
            double st = Math.sin(theta);
            double g = xCenterParam + (aParam * ca * ct) - (bParam * sa * st);
            ellipseX[angle] = (float)g;
            double h = yCenterParam + (aParam * sa * ct) + (bParam * ca * st);
            ellipseY[angle] = (float)h;
        }
        
        float[] xp = new float[points.size()];
        float[] yp = new float[xp.length];
        int i = 0;
        for (PairInt p : points) {
            xp[i] = p.getX();
            yp[i] = p.getY();
            i++;
        }
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(xMin, xMax, 
            yMin, yMax);
        
        plotter.addPlot(xp, yp, ellipseX, ellipseY, plotLabel);
        
        String fileName = plotter.writeFile(Integer.valueOf(plotNumber));
        
        log.info("wrote to file: " + fileName);
    }
        
}
