package algorithms.util;

import algorithms.MultiArrayMergeSort;
import algorithms.imageProcessing.DFSConnectedGroupsFinder;
import algorithms.misc.MiscMath;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.QRDecomposition;
import org.ejml.simple.SimpleMatrix;

/**
 * methods associated with fitting a 2nd order polynomial curve.
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
        
        SecureRandom sr = new SecureRandom();
        
        sr.setSeed(System.currentTimeMillis());
        
        return solveAfterRandomSampling(points, sr);
    }
    
    /**
     * solve for 2nd order curve for a random sample of 1000 points from
     * (dataX, dataY).
     * 
     * @param points
     * @return 
     */
    protected float[] solveAfterRandomSampling(Set<PairInt> points,
        SecureRandom sr) {
        
        int n = (points.size() > 1000) ? 1000 : points.size();
        
        List<PairInt> tmp = new ArrayList<PairInt>(points);
        
        int[] indexes = new int[n];
        MiscMath.chooseRandomly(sr, indexes, points.size());
        
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
    
    /**
     * calculate the square root of the sum of the squared differences between 
     * a 2nd order polygon defined by the given coefficients and the given 
     * points.
     * Note that if coefficients or points are null or empty, it returns
     * a result of infinity.
     * 
     * @param coefficients
     * @param points
     * @return 
     */
    public double calcResiduals(float[] coefficients, Set<PairInt> points) {
        
        if (points == null || points.isEmpty()) {
            return Double.POSITIVE_INFINITY;
        }
        
        if (coefficients == null || (coefficients.length != 3)) {
            return Double.POSITIVE_INFINITY;
        }
                
        int[] minMaxXY = MiscMath.findMinMaxXY(points);
        
        double xMin = minMaxXY[0];
        
        // these can be large, so use abs value instead of sum of squares
        
        double sum = 0;
        
        for (PairInt p : points) {
            
            double x = p.getX();
            
            double y = p.getY();            
            
            double yPoly = xMin + coefficients[0] + (coefficients[1]*x) 
                + (coefficients[2]*x*x);
            
            double diff = y - yPoly;
            
            sum += Math.abs(diff);
        }
                
        return sum;
    }

    public float[] solveForBestFittingContiguousSubSamples(Set<PairInt> points, 
        Set<PairInt> outputPoints, int imageWidth, int imageHeight) {
        
        DFSConnectedGroupsFinder groupsFinder = new DFSConnectedGroupsFinder();
        groupsFinder.setMinimumNumberInCluster(100);
        groupsFinder.findConnectedPointGroups(points, imageWidth, imageHeight);
        
        int n = groupsFinder.getNumberOfGroups();
        
        List<Set<PairInt>> contigList = new ArrayList<Set<PairInt>>();
        
        if (n == 1) {
            
            outputPoints.addAll(groupsFinder.getXY(0));
            
            return solveAfterRandomSampling(outputPoints);
            
        } else if (n > 7) {
            
            int[] groupIndexes = new int[n];
            int[] groupN = new int[n];
            for (int gId = 0; gId < n; gId++) {
                int n2 = groupsFinder.getNumberofGroupMembers(gId);
                groupIndexes[gId] = gId;
                groupN[gId] = n2;
            }
            MultiArrayMergeSort.sortByDecr(groupN, groupIndexes);
            
            n = 5;
            
            for (int i = 0; i < n; i++) {
                int gId = groupIndexes[i];
                Set<PairInt> s = groupsFinder.getXY(gId);
                contigList.add(s);
            }
            
        } else {
            
            for (int gId = 0; gId < n; gId++) {
                Set<PairInt> s = groupsFinder.getXY(gId);
                contigList.add(s);
            }
        }
        
        groupsFinder = null;
        
        SecureRandom sr = new SecureRandom();
        
        sr.setSeed(System.currentTimeMillis());
        
        
        double bestSubsetResiduals = Double.MAX_VALUE;
        
        float[] bestSubsetCoeff = null;
                
        Set<PairInt> bestSubsetPoints = null;
        
        // TODO: consider storing all coeff and bitstrings, so that best match
        // can aggregate points from similar fits
        
        List<Integer> setBits = new ArrayList<Integer>();
        
        for (int k = 1; k <= n; k++) {
            
            Long bitstring = MiscMath.getNextSubsetBitstring(n, k, null);
            
            long nExpected = MiscMath.computeNDivNMinusK(n, k)
                /MiscMath.factorial(k);
                        
            int count = 0;
            
            while ((bitstring != null) && (count < nExpected)) {

                setBits.clear();
                
                System.out.println("n=" + n + " k=" + k + " " +
                    Long.toBinaryString(bitstring.longValue()));
                 
                MiscMath.readSetBits(bitstring, setBits);
                
                Set<PairInt> subset = new HashSet<PairInt>();

                for (Integer bitIndex : setBits) {

                    int groupId = bitIndex.intValue();

                    Set<PairInt> g = contigList.get(groupId);

                    subset.addAll(g);
                }
                
                float[] coeff = solveAfterRandomSampling(subset, sr);
                
                double resid = calcResiduals(coeff, subset);
                
                if (resid < bestSubsetResiduals) {
                    
                    if (resid < (subset.size() * 5)) {
                        bestSubsetResiduals = resid;
                        bestSubsetCoeff = coeff;
                        bestSubsetPoints = subset;
                    }
                }
                
                bitstring = MiscMath.getNextSubsetBitstring(n, k, bitstring);
                
                count++;
            }
        }
        
        if (bestSubsetPoints != null) {
            outputPoints.clear();
            outputPoints.addAll(bestSubsetPoints);
        }
        
        return bestSubsetCoeff;
    }   
    
}
