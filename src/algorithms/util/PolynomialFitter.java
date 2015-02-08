package algorithms.util;

import algorithms.MultiArrayMergeSort;
import algorithms.imageProcessing.DFSConnectedGroupsFinder;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageDisplayer;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.misc.MiscMath;
import java.io.IOException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
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
     * @return 2nd order polynomial coefficients if solved, else null
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
     * @param sr instance of secure random to use for generating random numbers
     * @return 2nd order polynomial coefficients if solved, else null
     */
    protected float[] solveAfterRandomSampling(Set<PairInt> points,
        SecureRandom sr) {
        
        int n = (points.size() > 2500) ? 2500 : points.size();
        
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
     * @return 2nd order polynomial coefficients if solved, else null
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
    public double[] calcResidualsForAvg(float[] coefficients, Set<PairInt> points) {
        
        if (points == null || points.isEmpty()) {
            return new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        }
        
        if (coefficients == null || (coefficients.length != 3)) {
            return new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
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
            
            sum += diff;
        }
        
        double avg = sum/(double)points.size();
        
        sum = 0;
        
        for (PairInt p : points) {
            
            double x = p.getX();
            
            double y = p.getY();
            
            double yPoly = xMin + coefficients[0] + (coefficients[1]*x) 
                + (coefficients[2]*x*x);
            
            double diff = (y - yPoly) - avg;
            
            sum += (diff * diff);
        }
        
        double stDev = Math.sqrt(sum/((double)points.size() - 1.0));
                
        return new double[]{avg, stDev};
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
                
        // key = coefficients for a subset's fit, value = bitstring representing
        //   indexes to extract from contigList
        Map<CoefficientWrapper, Long> subsetCoeffMap = new 
            HashMap<CoefficientWrapper, Long>();
        
        List<Integer> setBits = new ArrayList<Integer>();
        
        //TODO: simplify the iteration thru subsets
        
        int nIter = 0;
        
        for (int k = n; k > 0; k--) {
            
            Long bitstring = MiscMath.getNextSubsetBitstring(n, k, null);
            
            long nExpected = MiscMath.computeNDivNMinusK(n, k)
                /MiscMath.factorial(k);
                        
            int count = 0;
            
            while ((bitstring != null) && ((nIter == 0) || (count < nExpected))) {

                setBits.clear();
                
                //System.out.println("n=" + n + " k=" + k + " " +
                //    Long.toBinaryString(bitstring.longValue())
                //    + " nIter=" + nIter);
                
                MiscMath.readSetBits(bitstring, setBits);
                
                Set<PairInt> subset = new HashSet<PairInt>();

                for (Integer bitIndex : setBits) {

                    int groupId = bitIndex.intValue();

                    Set<PairInt> g = contigList.get(groupId);

                    subset.addAll(g);
                }
                
                float[] coeff = solveAfterRandomSampling(subset, sr);
                
                if (coeff != null) {
                    subsetCoeffMap.put(new CoefficientWrapper(coeff), bitstring);
                }
                
                double resid = calcResiduals(coeff, subset);
/*                
String label = Long.toBinaryString(bitstring.longValue()) + " " 
+ Double.toString(resid) + " "
+ Arrays.toString(coeff);                
Image img = new Image(imageWidth, imageHeight);
try {
ImageIOHelper.addToImage(subset, 0, 0, img);
ImageDisplayer.displayImage(label, img);
} catch(IOException e) {
}
*/
                if (resid < bestSubsetResiduals) {
                    
                    //if (resid < (subset.size() * 5)) {
                        bestSubsetResiduals = resid;
                        bestSubsetCoeff = coeff;
                    //}
                }
                
                bitstring = MiscMath.getNextSubsetBitstring(n, k, bitstring);
                
                count++;
                
                nIter++;
            }
        }
        
        if (bestSubsetCoeff != null) {
            
            // find all fits similar to best fit and add subsets to outputPoints
            
            Iterator<Entry<CoefficientWrapper, Long>> iter = 
                subsetCoeffMap.entrySet().iterator();
            
            while (iter.hasNext()) {
                
                Entry<CoefficientWrapper, Long> entry = iter.next();
                
                float[] c = entry.getKey().getCoefficients();
                
                float diffC0 = Math.abs(bestSubsetCoeff[0] - c[0]);
                float divC1 = bestSubsetCoeff[1]/c[1];
                float divC2 = bestSubsetCoeff[2]/c[2];
                
                //TODO: this probably needs adjustment. 
                // diffC0 needs real world scale or relative size knowledge
                // divC1 comparison might be reduced to 0.1
                if ((Math.abs(divC1 - 1) < 0.2) && (Math.abs(divC2 - 1) < 0.6)
                    && (diffC0 < 20)) {
                    
                    setBits.clear();
                    
                    Long bitstring = entry.getValue();
                    
                    System.out.println("adding subsets: "
                        + Long.toBinaryString(bitstring.longValue()));
                    
                    MiscMath.readSetBits(bitstring, setBits);
                    for (Integer bitIndex : setBits) {
                        int groupId = bitIndex.intValue();
                        Set<PairInt> g = contigList.get(groupId);
                        outputPoints.addAll(g);
                    }
                }
            }
            
            // redo fit for outputPoints
            float[] coeff = solveAfterRandomSampling(outputPoints, sr);
            
            if (coeff == null) {
                return null;
            }
                        
            /*
            double[] avgAndStDevDiffs = calcResidualsForAvg(coeff, outputPoints);
            
            System.out.println("avgAndStDevDiffs[0]=" + avgAndStDevDiffs[0] +
                " avgAndStDevDiffs[1]=" + avgAndStDevDiffs[1]);
            */
            
            return coeff;
        }
        
        return bestSubsetCoeff;
    }   
    
}
