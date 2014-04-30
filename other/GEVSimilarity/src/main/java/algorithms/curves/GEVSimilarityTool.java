package algorithms.curves;

import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.GEVSimilarityParametersPlotter;
import algorithms.util.ResourceFinder;
import java.io.File;
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
    protected float mu0;
    protected float nWithinTen; 
    protected float nKPermutations;
    protected float nSigmaPermutations;
    protected float nMuPermutations;
    protected float fctr;
    protected int nCurves;
    
    protected float[][] curves;
    
    protected float[] x;
    protected float[] ks;
    protected float[] sigmas;
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
        mu0 = 1e-3f;
        
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
        nMuPermutations = nWithinTen * powersOfTenM;
        fctr = 10.f/nWithinTen;
    
        nCurves = (int)(nKPermutations * nSigmaPermutations * 
            nMuPermutations);
            
        curves = new float[nCurves][];
        for (int i = 0; i < nCurves; i++) {
            curves[i] = new float[nx];
        }        
    
        ks = new float[nCurves];
        sigmas = new float[nCurves];
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
    
    protected void generateCurves() {
        
        log.info("generateCurves nCurves=" + nCurves + " nWithinTen=" + nWithinTen + " fctr=" + fctr);

        long t0 = System.nanoTime();
        
        float k00 = k0;
        for (int i = 0; i < nKPermutations; i++) {
            float nk = (i % nWithinTen);
            if ((i > 0) && (nk == 0.f)) {
                k00 *= 10;
            }
            float k = (nk == 0) ? k00 : nk*fctr*k00;
            float s00 = sigma0;
            for (int ii = 0; ii < nSigmaPermutations; ii++) { 
                float ns = (ii % nWithinTen);
                if ((ii > 0) && (ns == 0.f)) {
                    s00 *= 10;
                }
                float sigma = (ns == 0) ? s00 : ns*fctr*s00;
                
                float m00 = mu0;
                for (int iii = 0; iii < nMuPermutations; iii++) {
                    float nxm = (iii % nWithinTen);
                    if ((iii > 0) && (nxm == 0.f)) {
                        m00 *= 10;
                    }
                    float mu = (nxm == 0) ? m00 : nxm*fctr*m00;
                         
                    float[] y = GeneralizedExtremeValue.generateNormalizedCurve(x, k, sigma, mu);
                        
                    if (!isNearlyAllZeros(y) && !isNearlyAllOnes(y)) {
                        if (nc == 76) {
                            System.out.println("*****k=" + k + " sigma=" + sigma + " mu" + mu);
                        }
                        curves[nc] = y;
                        ks[nc] = k;
                        sigmas[nc] = sigma;
                        mus[nc] = mu;
                        nc++;
                    }
                }
            }
        }
        
        log.info("generated nc=" + nc + " curves" );
        
        // expecting runtime complexity to be O(nCurves) for creating curves
        long t1 = System.nanoTime();
        long time = (t1 - t0)/1000000000;
        log.info("runtime for creating curves = " + time + " seconds");
    }

    public void calculateCurveDiffs() {
        
        log.info("calculateCurveDiffs");

        generateCurves();
        

        long t1 = System.nanoTime();
                
        // capture similarity of a spectra compared to all others as a set of
        //    similar.
        double epsDiff = 0.1;
        boolean[] inASetAlready = new boolean[nCurves];
        
        // store diff and each member of similarity set
        log.info("looking for similarities between " + nc + " curves");
        Arrays.fill(diff, Float.MAX_VALUE);
        
        StringBuffer sb = new StringBuffer();
        
        int end = nc;
        
        nc = 0;
        for (int i = 0; i < end; i++) {
            if (inASetAlready[i]) {
                continue;
            }
            
            float[] y0 = curves[i];
            
            if (sb.length() > 0) {
                sb = sb.delete(0, sb.length());
            }
            
            for (int ii = (i + 1); ii < end; ii++) {  // cost  is less than  nCurves * nCurves
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
                log.info(String.format("UNIQUE:  k=%e  sigma=%e  mu=%e", ks[i], sigmas[i], mus[i]));
            }
        }
     
        calculatedCurveDiffs = true;

        // nc is < nCurves
        // expected runtime complexity for diff finding is less than nc * nc *
        long t2 = System.nanoTime();
        long time = (t2 - t1)/1000000000;
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
            
            for (int nv = 0; nv < 3; nv++) {
                float[] a;
                switch(nv) {
                    case 0:
                        sb2.append("\n   k=");
                        a = ks;
                        break;
                    case 1:
                        sb2.append("  sigma=");
                        a = sigmas;
                        break;
                    default:
                        sb2.append("  mu=");
                        a = mus;
                        break;
                }
                for (int j = 0; j < indexes.length; j++) {
                    float v = a[indexes[j]];
                    sb2.append(v);
                    if (j < (indexes.length - 1)) {
                        sb2.append(",");
                    }
                }
            }

            log.info(sb2.toString());            

            // plot the parameters
            if ((i == 0) || ((i % 100) == 0) ) {
                                
                if (paramsPlotter != null) {
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
    
    private String getFilePath(int fileNumber) throws IOException {
        String fileName = "similarity_" + Integer.toString(fileNumber) + ".dat";
        String filePath = ResourceFinder.findDirectory("tmpdata2");
        filePath = filePath + System.getProperty("file.separator") + fileName;
        return filePath;
    }
    
    public void readPersisted(int fileNumber) throws IOException {
        
        log.info("readPersisted");
        
        String filePath = getFilePath(fileNumber);   
        
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
            mu0 = ois.readFloat();
            nWithinTen = ois.readFloat();
            nc = ois.readInt();
            nKPermutations = nWithinTen * powersOfTenK;
            nSigmaPermutations = nWithinTen * powersOfTenS;
            nMuPermutations = nWithinTen * powersOfTenM;
            fctr = 10.f/nWithinTen;
            nCurves = (int)(nKPermutations * nSigmaPermutations * 
                nMuPermutations);
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
        
        String filePath = getFilePath(fileNumber);
        
        FileOutputStream fos = null;
        ObjectOutputStream oos = null;
        
        try {
            File fl = new File(filePath);
            fl.delete();
            fl.createNewFile();
            fos = new FileOutputStream(filePath);
            oos = new ObjectOutputStream(fos);
            
            oos.writeInt(powersOfTenK);
            oos.writeInt(powersOfTenS);
            oos.writeInt(powersOfTenM);
            oos.writeInt(nx);
            oos.writeFloat(k0);
            oos.writeFloat(sigma0);
            oos.writeFloat(mu0);
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
