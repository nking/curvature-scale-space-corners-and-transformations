package algorithms.util;

import algorithms.misc.MiscMath;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.QRDecomposition;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author nichole
 */
public class PolynomialFitter {
    
    /**
     * solve for 2nd order curve for a random sample of 1000 points from
     * (dataX, dataY).
     * 
     * @param points
     * @return 
     */
    public float[] solveAfterRandomSampling(Set<PairInt> points) {
        
        int n = (points.size() > 1000) ? 1000 : points.size();
        
        List<PairInt> tmp = new ArrayList<PairInt>(points);
        
        SecureRandom sr = new SecureRandom();
        sr.setSeed(System.currentTimeMillis());
        
        int[] indexes = new int[n];
        MiscMath.chooseRandomly(sr, indexes, points.size() - 1);
        
        float[] xP = new float[n];
        float[] yP = new float[xP.length];
        int i = 0;
        for (int idx : indexes) {
            PairInt p = tmp.get(idx);
            xP[i] = p.getX();
            yP[i] = p.getY();
            i++;
        }
        tmp = null;
        
        return solve(xP, yP);
    }
    
    /**
     * solve for 2nd order curves.
     * 
     * @param dataX
     * @param dataY
     * @return 
     */
    public float[] solve(float[] dataX, float[] dataY) {
        
        // from http://rosettacode.org/wiki/Polynomial_Fitting
        // adapted from the Go solution
        
        int m = dataX.length;
        
        int polyDegree = 2;
        
        int n = polyDegree + 1;
        
        SimpleMatrix y = new SimpleMatrix(m, 1);
        
        SimpleMatrix x = new SimpleMatrix(m, n);
        
        for (int i = 0; i < m; i++) {
            y.set(i, 0, dataY[i]);
            float ip = 1.f;
            for (int j = 0; j < n; j++) {
                x.set(i, j, ip);
                ip *= dataX[i];
            }
        }
        QRDecomposition<DenseMatrix64F> qrDecomp = DecompositionFactory.qr(m, n);
        if (!qrDecomp.decompose(x.getMatrix())) {
            return null; 
        }
       
        DenseMatrix64F q = qrDecomp.getQ(null, false);
        DenseMatrix64F r = qrDecomp.getR(null, false);
        
        SimpleMatrix qq = new SimpleMatrix(q);
        
        SimpleMatrix qty = qq.transpose().mult(y);
        
        float[] c = new float[n];

        for (int i = (n - 1); i >= 0; i--) {
            c[i] = (float)qty.get(i, 0);
            for (int j = (i + 1); j < n; j++) {
                c[i] -= c[j] * r.get(i, j);
            }
            c[i] /= r.get(i, i);
        }
        
        return c;
    }

    /**
     * plot the points and 2nd order curve given the coefficients.
     * note that if there are more than 1000 points, only 1000 of
     * them are 
     * @param coefficients
     * @param points
     * @param plotXMax
     * @param plotYMax
     * @param plotNumber
     * @param plotLabel
     * @return 
     */
    public String plotFit(float[] coefficients, Set<PairInt> points, int plotXMax,
    int plotYMax, int plotNumber, String plotLabel) {
        
        // shape the rainbow points into a more even ribbon
        float xMin = Float.MAX_VALUE;
        float xMax = Float.MIN_VALUE;
        float[] xP = new float[points.size()];
        float[] yP = new float[xP.length];
        int i = 0;
        for (PairInt p : points) {
            float x = p.getX();
            if (x < xMin) {
                xMin = x;
            }
            if (x > xMax) {
                xMax = x;
            }
            xP[i] = x;
            yP[i] = p.getY();
            i++;
        }
        
        int nCurve = 100;
        float dx = (xMax - xMin)/(float)nCurve;
        float[] xPoly = new float[nCurve];
        float[] yPoly = new float[xPoly.length];
        for (i = 0; i < nCurve; i++) {
            xPoly[i] = xMin + i*dx;
            yPoly[i] = coefficients[0] + coefficients[1]*xPoly[i] + coefficients[2]*xPoly[i]*xPoly[i];
        }
        
        try {
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0, 
                plotXMax, 0, plotYMax);
            plotter.addPlot(xP, yP, xPoly, yPoly, plotLabel);
            
            String fileName = plotter.writeFile(Integer.valueOf(plotNumber));
            
            return fileName;
            
        } catch (IOException e) {
            Logger.getLogger(this.getClass().getName()).severe(e.getMessage());
        }
        
        return null;
    }
}
