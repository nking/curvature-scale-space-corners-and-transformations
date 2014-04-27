package algorithms.curves;

import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.GEVSimilarityParametersPlotter;
import algorithms.util.ResourceFinder;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.logging.Logger;

/**
 * tool to create a text file and images of the calculated similarities between
 * GEV curves with the given set parameter ranges and steps.
 *
 * @author nichole
 */
public class GEVSimilarityTool {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    protected int powersOfTenK;
    protected int powersOfTenS;
    protected int powersOfTenM;
    protected int nx;
    protected float k0;
    protected float sigma0;
    protected float xMinusMu0;
    protected float nWithinTen; // when = 10.f, runtime = 8 hrs
    protected float nKPermutations;
    protected float nSigmaPermutations;
    protected float nXMinusMuPermutations;
    protected float fctr;
    protected int nCurves;
    
    protected float[][] curves;
    
    protected float[] x;
    protected float[] ks;
    protected float[] sigmas;
    protected float[] xminusmus;
    protected float[] mus;
    protected int nc;
        
    protected float[] diff;
    protected char[][] diffIndexes;
        
    protected boolean calculatedCurveDiffs;
    
    public GEVSimilarityTool() {
        init();
    }
    
    private void init() {        
        powersOfTenK = 6;
        powersOfTenS = 5;
        powersOfTenM = 4;
        nx = 20;
        k0 = 1e-3f;
        sigma0 = 1e-2f;
        xMinusMu0 = 1e-3f;
        
        nc = 0;
        
        x = new float[nx];
        x[0] = 0.0f;
        for (int i = 1; i < nx; i++) {
            x[i] = 1.f/(float)i;
        }

        calculatedCurveDiffs = false;
        
        initNWithinTenVars(10.f);
    }
    
    private void initNWithinTenVars(float nwt) {
        
        nWithinTen = nwt;
        nKPermutations = nWithinTen * powersOfTenK;
        nSigmaPermutations = nWithinTen * powersOfTenS;
        nXMinusMuPermutations = nWithinTen * powersOfTenM;
        fctr = 10.f/nWithinTen;
    
        nCurves = (int)(nKPermutations * nSigmaPermutations * 
            nXMinusMuPermutations * nx);
            
        curves = new float[nCurves][];
        for (int i = 0; i < nCurves; i++) {
            curves[i] = new float[nx];
        }        
    
        ks = new float[nCurves];
        sigmas = new float[nCurves];
        xminusmus = new float[nCurves];
        mus = new float[nCurves];
        
        diff = new float[nCurves];
        diffIndexes = new char[nCurves][];
        
        fctr = 10.f/nWithinTen;
    }
    
    protected void resetForNWithinTen(float nwt) {
        
        log.info("resetForNWithinTen");
        
        if (calculatedCurveDiffs) {
            log.severe("resetForNWithinTen can only be used before calculating cureves or readPersisted");
        }
        
        initNWithinTenVars(nwt);
    }
    
    public void calculateCurveDiffs() {
        
        log.info("calculateCurveDiffs nCurves=" + nCurves);

        long t0 = System.nanoTime();
        
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

        // expecting runtime complexity to be O(nCurves) for creating curves
        long t1 = System.nanoTime();
        long time = (t1 - t0)*1000000000;
        log.info("runtime for creating curves = " + time + " seconds");
                
        // capture similarity of a spectra compared to all others as a set of
        //    similar.
        double epsDiff = 1e-2;
        boolean[] inASetAlready = new boolean[nCurves];
        
        // store diff and each member of similarity set
        log.info("nCurves=" + nCurves);
        Arrays.fill(diff, Float.MAX_VALUE);
        
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
            
            for (int ii = i; ii < nCurves; ii++) {  // cost  is less than  nCurves * nCurves * nx
                if (inASetAlready[ii]) {
                    continue;
                }
                
                float[] y1 = curves[ii];
                
                double diffSumSSq = 0.f;
                for (int iii = 0; iii < y0.length; iii++) {
                    float d = y0[iii] - y1[iii];
                    diffSumSSq += (d * d);
                }
                
                //if (diffSumSSq > 0) {
                                        
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
                //}
            }
            if (sb.length() > 0) {
                diffIndexes[nc] = sb.toString().toCharArray();
                nc++;
            }
            if (!inASetAlready[i]) {
                // this is unique
                log.info(String.format("UNIQUE:  k=%e  sigma=%e  x-mu=%e  mu=%e", ks[i], sigmas[i], xminusmus[i], mus[i]));
            }
        }
     
        calculatedCurveDiffs = true;

        // expected runtime complexity for diff finding is less than nCurves * nCurves * nx
        // for default params, nCurves=(10*6 + 10*5 + 10*4)*20=3000 so complexity is less than O(180_000_000)
        long t2 = System.nanoTime();
        time = (t2 - t1)*1000000000;
        log.info("runtime for finding similar among curves = " + time + " seconds");
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

    public void plotResults() throws IOException {
        
        log.info("plotResults");
        
        // sort by similarity
        //sort(diff, diffIndexes, 0 , nc - 1);
        
        if (!calculatedCurveDiffs) {
            log.severe("need to calculate curve diffs before can plot results");
        }
        
        log.info("nc=" + nc);
        
        if (nc == 0) {
            return;
        }
        
        int end = (nc > 4000) ? 4000 : nc;
        //end = nc;
        int nci = -1;
        int nci2 = -1;
        
        PolygonAndPointPlotter plotter = null;
        
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
            if ((i == 0) || ((i % 100) == 0) ) {
                                
                if (paramsPlotter != null) {
System.out.println("i=" + i + " n=" + end + " nci2=" + nci2);
                    paramsPlotter.writeFile(nci2);
                }
                
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
            
            paramsPlotter.addToPlot(yy, kPoints, sPoints, mPoints);
            
        }
        
        paramsPlotter.writeFile(nci2);
    }

    public void readPersisted() throws IOException {
        readPersisted(0);
    }
    
    public void readPersisted(int fileNumber) throws IOException {
        
        log.info("readPersisted");
        
        String fileName = "similarity_" + Integer.toString(fileNumber) + ".dat";
        String filePath = ResourceFinder.getAFilePathInTmpData(fileName);
        
        FileInputStream fis = null;
        ObjectInputStream ois = null;
        
        try {
            fis = new FileInputStream(filePath);
            ois = new ObjectInputStream(fis);
            
            powersOfTenK = ois.readInt();
            powersOfTenS = ois.readInt();
            powersOfTenM = ois.readInt();
            nx = ois.readInt();
            k0 = ois.readFloat();
            sigma0 = ois.readFloat();
            xMinusMu0 = ois.readFloat();
            nWithinTen = ois.readFloat();
            nc = ois.readInt();
            
            nKPermutations = nWithinTen * powersOfTenK;
            nSigmaPermutations = nWithinTen * powersOfTenS;
            nXMinusMuPermutations = nWithinTen * powersOfTenM;
            fctr = 10.f/nWithinTen;
            nCurves = (int)(nKPermutations * nSigmaPermutations * 
                nXMinusMuPermutations * nx);
            
            curves = new float[nCurves][];
            for (int i = 0; i < nCurves; i++) {
                curves[i] = new float[nx];
                for (int j = 0; j < nx; j++) {
                    curves[i][j] = ois.readFloat();
                }
            }
            
            x = new float[nx];
            for (int j = 0; j < nx; j++) {
                x[j] = ois.readFloat();
            }
            
            ks = new float[nCurves];
            for (int i = 0; i < nCurves; i++) {
                ks[i] = ois.readFloat();
            }
            
            sigmas = new float[nCurves];
            for (int i = 0; i < nCurves; i++) {
                sigmas[i] = ois.readFloat();
            }
            
            xminusmus = new float[nCurves];
            for (int i = 0; i < nCurves; i++) {
                xminusmus[i] = ois.readFloat();
            }
            
            mus = new float[nCurves];
            for (int i = 0; i < nCurves; i++) {
                mus[i] = ois.readFloat();
            }
            
            diff = new float[nc];
            for (int i = 0; i < nc; i++) {
                diff[i] = ois.readFloat();
            }
                        
            for (int i = 0; i < nc; i++) {
                int len = ois.readInt();
                diffIndexes[i] = new char[len];
                for (int j = 0; j < len; j++) {
                    diffIndexes[i][j] = ois.readChar();
                }
            }
                        
            calculatedCurveDiffs = true;
            
        } catch (FileNotFoundException e) {
            log.severe(e.getMessage());
        } catch (IOException e) {
            log.severe(e.getMessage());
        } finally {
            if (fis != null) {
                try { fis.close();} catch (IOException e2){ }
            }
            if (ois != null) {
                try { ois.close();} catch (IOException e2){ }
            }
        }
    }
    
    public void persist() throws IOException {
        persist(0);
    }
    
    public void persist(int fileNumber) throws IOException {
        
        log.info("persist");
        
        String fileName = "similarity_" + Integer.toString(fileNumber) + ".dat";
        String filePath = ResourceFinder.getAFilePathInTmpData(fileName);
        
        FileOutputStream fos = null;
        ObjectOutputStream oos = null;
        
        try {
            fos = new FileOutputStream(filePath);
            oos = new ObjectOutputStream(fos);
            
            oos.writeInt(powersOfTenK);
            oos.writeInt(powersOfTenS);
            oos.writeInt(powersOfTenM);
            oos.writeInt(nx);
            oos.writeFloat(k0);
            oos.writeFloat(sigma0);
            oos.writeFloat(xMinusMu0);
            oos.writeFloat(nWithinTen);
            oos.writeInt(nc);
            
            for (int i = 0; i < nCurves; i++) {
                for (int j = 0; j < nx; j++) {
                    oos.writeFloat(curves[i][j]);
                }
            }
            
            for (int j = 0; j < nx; j++) {
                oos.writeFloat(x[j]);
            }
            
            for (int i = 0; i < nCurves; i++) {
                oos.writeFloat(ks[i]);
            }
            
            for (int i = 0; i < nCurves; i++) {
                oos.writeFloat(sigmas[i]);
            }
            
            for (int i = 0; i < nCurves; i++) {
                oos.writeFloat(xminusmus[i]);
            }
            
            for (int i = 0; i < nCurves; i++) {
                oos.writeFloat(mus[i]);
            }
        
            for (int i = 0; i < nc; i++) {
                oos.writeFloat(diff[i]);
            }
            
            for (int i = 0; i < nc; i++) {
                int len = (diffIndexes[i] == null) ? 0 : diffIndexes[i].length;
                oos.writeInt(len);
                for (int j = 0; j < len; j++) {
                    oos.writeChar(diffIndexes[i][j]);
                }
                oos.flush();
            }
            
            oos.flush();
            
        } catch (FileNotFoundException e) {
            log.severe(e.getMessage());
        } catch (IOException e) {
            log.severe(e.getMessage());
        } finally {
            if (oos != null) {
                oos.close();
            }
            if (fos != null) {
                fos.close();
            }
        }
        
    }
}
