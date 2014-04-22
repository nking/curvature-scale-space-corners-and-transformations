package algorithms.curves;

import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.GEVSimilarityParametersPlotter;
import java.util.Arrays;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 * a tool to find common curves among the large range of possible curves
 * in the GEV function which depends upon k, sigma, and (x-mu).
 * 
 * @author nichole
 */
public class GEVSimilarityToolTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    protected boolean debug = true;

    protected boolean enable = true;

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void test0() throws Exception {

        log.info("test0()");

        if (!enable) {
            return;
        }

        int nx = 20;
        float[] x = new float[nx];
        for (int i = 0; i < nx; i++) {
            x[i] = 1.f/(float)i;
        }
        
        PolygonAndPointPlotter plotter = null;

        // for nWithinTen=4 and starts at 1e-3 with powersOfTen=5, test0    took 232.141 sec 
        /* int powersOfTenK = 8; int powersOfTenS = 5; int powersOfTenM = 4;
        float k0 = 1e-5f; float sigma0 = 1e-2f; float xMinusMu0 = 1e-3f;
        float nWithinTen = 4.f;
            ==> took 408 sec
        */
        
        /*
        int powersOfTenK = 6;
        int powersOfTenS = 5;
        int powersOfTenM = 4;
        float k0 = 1e-3f;
        float sigma0 = 1e-2f;
        float xMinusMu0 = 1e-3f;        
        float nWithinTen = 10.f;
            ===> 519 min = 8.65 h
        */
        
        int powersOfTenK = 6;
        int powersOfTenS = 5;
        int powersOfTenM = 4;
        
        float k0 = 1e-3f;
        float sigma0 = 1e-2f;
        float xMinusMu0 = 1e-3f;
        
        float nWithinTen = 1.f;
        
        float nKPermutations = nWithinTen * powersOfTenK;
        float nSigmaPermutations = nWithinTen * powersOfTenS;
        float nXMinusMuPermutations = nWithinTen * powersOfTenM;
        float fctr = 10.f/nWithinTen;
        
        int nCurves = (int)(nKPermutations * nSigmaPermutations * 
            nXMinusMuPermutations * nx);
        
        float[][] curves = new float[nCurves][];
        for (int i = 0; i < nCurves; i++) {
            curves[i] = new float[nx];
        }
        float[] ks = new float[nCurves];
        float[] sigmas = new float[nCurves];
        float[] xminusmus = new float[nCurves];
        float[] mus = new float[nCurves];
        
        int nc = 0;
        
        float k00 = k0;
        for (int i = 0; i < nKPermutations; i++) {
            float nk = (i % nWithinTen);
            if (nk == 0.f) {
                k00 *= 10;
            }
            float k = (nk == 0) ? k00 : nk*fctr*k00;
            
            float s00 = sigma0;
            for (int ii = 0; ii < nSigmaPermutations; ii++) { 
                float ns = (ii % nWithinTen);
                if (ns == 0.f) {
                    s00 *= 10;
                }
                float sigma = (ns == 0) ? s00 : ns*fctr*s00;
                
                float m00 = xMinusMu0;
                for (int iii = 0; iii < nXMinusMuPermutations; iii++) {
                    float nxm = (iii % nWithinTen);
                    if (nxm == 0.f) {
                        m00 *= 10;
                    }
                    float xMinusMu = (nxm == 0) ? m00 : nxm*fctr*m00;
                    
                    for (int iiii = 0; iiii < nx; iiii++) {
                        float x0 = x[iiii];
                        float mu = x0 - xMinusMu;
                        
                        float[] y = GeneralizedExtremeValue.generateNormalizedCurve(x, k, sigma, mu);
                        
                        if (!isNearlyAllZeros(y) && !isNearlyAllOnes(y)) {
                            curves[nc] = y;
                            ks[nc] = k;
                            sigmas[nc] = sigma;
                            xminusmus[nc] = xMinusMu;
                            mus[nc] = mu;
                            nc++;
                        }
                    }
                }
            }
        }
                
        // capture similarity of a spectra compared to all others as a set of
        //    similar.
        double epsDiff = 1e-2;
        boolean[] inASetAlready = new boolean[nCurves];
        
        // store diff and each member of similarity set
        // float[] diff;  char[][]indexes
        log.info("nCurves=" + nCurves);
        float[] diff = new float[nCurves];
        Arrays.fill(diff, Float.MAX_VALUE);
        char[][] diffIndexes = new char[nCurves][];
        
        StringBuffer sb = new StringBuffer();
        
        nc = 0;
        for (int i = 0; i < nCurves; i++) {
            if (inASetAlready[i]) {
                continue;
            }
            
            float[] y0 = curves[i];
            
            if (sb.length() > 0) {
                sb = sb.delete(0, sb.length());
            }
            
            for (int ii = i; ii < nCurves; ii++) {
                if (inASetAlready[ii]) {
                    continue;
                }
                
                float[] y1 = curves[ii];
                
                double diffSumSSq = 0.f;
                for (int iii = 0; iii < y0.length; iii++) {
                    float d = y0[iii] - y1[iii];
                    diffSumSSq += (d * d);
                }
                
                if (diffSumSSq > 0) {
                                        
                    float d = (float) Math.sqrt(diffSumSSq);
                    
                    if (d < epsDiff) {
                        if (sb.length() == 0) {
                            diff[nc] = d;
                            sb.append(Integer.toString(i));
                            inASetAlready[i] = true;
                        }
                        sb.append(",").append(Integer.toString(ii));
                        inASetAlready[ii] = true;
                    }
                }
            }
            if (sb.length() > 0) {
                diffIndexes[nc] = sb.toString().toCharArray();
                nc++;
            }
        }
        
        // sort by similarity
        //sort(diff, diffIndexes, 0 , nc - 1);
        
        log.info("nc=" + nc);
        
        int end = (nc > 4000) ? 4000 : nc;
        //end = nc;
        int nci = -1;
        int nci2 = -1;
        
        
        float k1 = (float) (k0 * Math.pow(10, powersOfTenK));
        float sigma1 = (float) (sigma0 * Math.pow(10, powersOfTenS));
        float m0 = -10;
        float m1 = 10;
        GEVSimilarityParametersPlotter paramsPlotter = null;
                
        for (int i = 0; i < end; i++) {
            String[] indexesStr = new String(diffIndexes[i]).split(",");
            int[] indexes = new int[indexesStr.length];
            for (int j = 0; j < indexesStr.length; j++) {
                indexes[j] = Integer.valueOf(indexesStr[j]).intValue();
            }
            
            float[] y = curves[indexes[0]];
            
            String lbl = String.format("%d)", i);
            
            if ((i == 0) || (i == 1000) || (i == 2000)){
                nci++;
                plotter = new PolygonAndPointPlotter();
            }
            plotter.addPlot(x, y, x, y, lbl);
            switch (nci) {
                case 0:
                    plotter.writeFile();
                    break;
                case 1:
                    plotter.writeFile2();
                    break;
                default:
                    plotter.writeFile3();
                    break;
            } 
            
            StringBuffer sb2 = new StringBuffer(100);
            sb2.append(Integer.toString(i)).append(") [").append(new String(diffIndexes[i])).append("]")
                .append("  diff=").append(Float.toString(diff[i]));
            sb2.append("\n   k=");
            for (int j = 0; j < indexes.length; j++) {
                float k = ks[indexes[j]];
                sb2.append(k);
                if (j < (indexes.length - 1)) {
                    sb2.append(",");
                }
            }
            sb2.append("  sigma=");
            for (int j = 0; j < indexes.length; j++) {
                float s = sigmas[indexes[j]];
                sb2.append(s);
                if (j < (indexes.length - 1)) {
                    sb2.append(",");
                }
            }
            sb2.append("  mu=");
            for (int j = 0; j < indexes.length; j++) {
                float m = mus[indexes[j]];
                sb2.append(m);
                if (j < (indexes.length - 1)) {
                    sb2.append(",");
                }
            }
            sb2.append("  x-mu=");
            for (int j = 0; j < indexes.length; j++) {
                float m = xminusmus[indexes[j]];
                sb2.append(m);
                if (j < (indexes.length - 1)) {
                    sb2.append(",");
                }
            }

            log.info(sb2.toString());            

            // plot the parameters
            if ((i == 0) || ((i % 100) == 0) ){
                nci2++;
                int ymin = i;
                paramsPlotter = new GEVSimilarityParametersPlotter(
                    ymin, ymin + 100, k0, k1, sigma0, sigma1, m0, m1);
            }
            
            float[] yy = new float[indexes.length];
            Arrays.fill(yy, i);
            float[] kPoints = new float[indexes.length];
            float[] sPoints = new float[indexes.length];
            float[] mPoints = new float[indexes.length];
            for (int j = 0; j < indexes.length; j++) {
                kPoints[j] = ks[ indexes[j] ];
                sPoints[j] = sigmas[ indexes[j] ];
                mPoints[j] = mus[ indexes[j] ];
            }
                        
            paramsPlotter.addPlot(yy, kPoints, sPoints, mPoints);
                        
            paramsPlotter.writeFile(nci2);
        }        
    }
  
    protected boolean isNearlyAllZeros(float[] y) {
        int n = 0;
        for (int i = 0; i < y.length; i++) {
            if (y[i] < 1E-2) {
                n++;
            }
        }
        return (n > (y.length - 2));
    }
    
    protected boolean isNearlyAllOnes(float[] y) {
        int n = 0;
        for (int i = 0; i < y.length; i++) {
            if (y[i] > 0.98) {
                n++;
            }
        }
        return (n > (y.length - 2));
    }

    private void sort(float[] a, char[][] b) {
        if (a == null || b == null) {
            return;
        }
        log.info("sort " + a.length + " items");
        sort(a, b, 0, a.length - 1);
    }

    /** 
     * use quick sort to sort a by increasing value and perform same swaps on b
     * @param a
     * @param b
     * @param idxLo
     * @param idxHi last index to be sorted in a, inclusive
     */
    private void sort(float[] a, char[][] b, int idxLo, int idxHi) {
        if (idxLo < idxHi) {
            int idxMid = partition(a, b, idxLo, idxHi);
            sort(a, b, idxLo, idxMid - 1);
            sort(a, b, idxMid + 1, idxHi);
        }
    }

    private int partition(float[] a, char[][] b, int idxLo, int idxHi) {
        float x = a[idxHi];
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            if (a[i] <= x) {
                store++;
                float swap = a[store];
                a[store] = a[i];
                a[i] = swap;
                char[] swap2 = b[store];
                b[store] = b[i];
                b[i] = swap2;
            }
        }
        
        float swap = a[store + 1];
        a[store + 1] = a[idxHi];
        a[idxHi] = swap;
        
        char[] swap2 = b[store + 1];
        b[store + 1] = b[idxHi];
        b[idxHi] = swap2;
        
        return store + 1;
    }
}
