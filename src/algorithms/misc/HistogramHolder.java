package algorithms.misc;

import algorithms.util.PolygonAndPointPlotter;
import java.io.IOException;

/**
 * @author nichole
 */
public class HistogramHolder {

    protected float[] xHist = null;
    protected int[] yHist = null;
    protected float[] yHistFloat = null;
    protected float[] yErrors = null;
    protected float[] xErrors = null;
    
    public float[] getHistArea(float maxXToUse, int nPartitions) {
        
        if (yHistFloat == null) {
            return null;
        }
        
        double[] area = new double[nPartitions];
        
        float binSize = maxXToUse/(float)nPartitions;
        
        // trapezoidal rule for area under the curve
        
        for (int i = 0; i < (xHist.length - 1); ++i) {
            
            float yTerm = yHistFloat[i + 1] + yHistFloat[i];
            float xLen = xHist[i + 1] - xHist[i];
            if (xLen < 0) {
                xLen *= -1;
            }
            
            float x = xHist[i];
            
            int partition = (int)(x/binSize);
            
            if (partition > (nPartitions - 1)) {
                partition = nPartitions - 1;
            }
            
            area[partition] += (yTerm * xLen);
        }
                
        double sum = 0;
        for (int i = 0; i < nPartitions; ++i) {
            sum += area[i]; // area should be multiplied by 0.5, but that's not needed for normalization
        }
        
        float[] frac = new float[nPartitions];
        for (int i = 0; i < nPartitions; ++i) {
            frac[i] = (float)(area[i]/sum);
        }
        
        return frac;
    }

    public int calculateHalfYMaxIndexPastYMax() {

        if (yHistFloat == null) {
            return -1;
        }

        int yMaxIndex = MiscMath.findYMaxIndex(yHistFloat);

        int halfMaxIndex = -1;
        float halfMax = yHistFloat[yMaxIndex]/2.0f;

        for (int i = yMaxIndex; i < yHistFloat.length; i++) {
            if (halfMax <= yHistFloat[i]) {
                halfMaxIndex = i;
            }
        }
        return halfMaxIndex;
    }
    
    public String plotHistogram(String label, 
        long outputFileNumber) throws IOException {
                
        float[] xh = xHist;
        float[] yh = yHistFloat;
        
        float yMin = MiscMath.findMin(yh);
        int yMaxIdx = MiscMath.findYMaxIndex(yh);
        float yMax = yh[yMaxIdx];
        
        float xMin = MiscMath.findMin(xh);
        float xMax = MiscMath.findMax(xh);        
                
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        plotter.addPlot(
            xMin, xMax, yMin, yMax,
            xh, yh, xh, yh, label);

        return plotter.writeFile(outputFileNumber);
    }
    
    public String plotHistogram(String label, 
        String outputFileSuffix) throws IOException {
                
        float[] xh = xHist;
        float[] yh = yHistFloat;
        
        float yMin = MiscMath.findMin(yh);
        int yMaxIdx = MiscMath.findYMaxIndex(yh);
        if (yMaxIdx == -1) {
            return null;
        }
        float yMax = yh[yMaxIdx];
        
        float xMin = MiscMath.findMin(xh);
        float xMax = MiscMath.findMax(xh);        
                
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        plotter.addPlot(
            xMin, xMax, yMin, yMax,
            xh, yh, xh, yh, label);

        return plotter.writeFile(outputFileSuffix);
    }
    
    public String plotHistogram(float xMin, float xMax, String label, 
        String outputFileSuffix) throws IOException {
                
        float[] xh = xHist;
        float[] yh = yHistFloat;
        
        float yMin = MiscMath.findMin(yh);
        int yMaxIdx = MiscMath.findYMaxIndex(yh);
        if (yMaxIdx == -1) {
            return null;
        }
        float yMax = yh[yMaxIdx];
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        plotter.addPlot(
            xMin, xMax, yMin, yMax,
            xh, yh, xh, yh, label);

        return plotter.writeFile(outputFileSuffix);
    }
    
    public String plotLogHistogram(String label, 
        String outputFileSuffix) throws IOException {
                
        float[] xh = xHist;
        float[] yh = yHistFloat;
        
        float[] yLogH = new float[yh.length];
        for (int i = 0; i < yh.length; ++i) {
            yLogH[i] = (float)Math.log(yh[i]/Math.log(10));
        }
        
        float yMin = MiscMath.findMin(yLogH);
        int yMaxIdx = MiscMath.findYMaxIndex(yLogH);
        float yMax = yLogH[yMaxIdx];
        
        float xMin = MiscMath.findMin(xh);
        float xMax = MiscMath.findMax(xh);
                        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        plotter.addPlot(
            xMin, xMax, yMin, yMax,
            xh, yLogH, xh, yLogH, label);

        return plotter.writeFile(outputFileSuffix);
    }

    /**
     * @return the xHist
     */
    public float[] getXHist() {
        return xHist;
    }

    /**
     * @return the yHist
     */
    public int[] getYHist() {
        return yHist;
    }

    /**
     * @return the yHistFloat
     */
    public float[] getYHistFloat() {
        return yHistFloat;
    }

    /**
     * @return the yErrors
     */
    public float[] getYErrors() {
        return yErrors;
    }

    /**
     * @return the xErrors
     */
    public float[] getXErrors() {
        return xErrors;
    }

    /**
     * @param xHist the xHist to set
     */
    public void setXHist(float[] xHist) {
        this.xHist = xHist;
    }

    /**
     * @param yHist the yHist to set
     */
    public void setYHist(int[] yHist) {
        this.yHist = yHist;
    }

    /**
     * @param yHistFloat the yHistFloat to set
     */
    public void setYHistFloat(float[] yHistFloat) {
        this.yHistFloat = yHistFloat;
    }

    /**
     * @param yErrors the yErrors to set
     */
    public void setYErrors(float[] yErrors) {
        this.yErrors = yErrors;
    }

    /**
     * @param xErrors the xErrors to set
     */
    public void setXErrors(float[] xErrors) {
        this.xErrors = xErrors;
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("histogram=[");
        
        for (int i = 0; i < xHist.length; i++) {
            
            sb.append("(").append(xHist[i]).append(", ");
            
            if (yHist != null) {
                sb.append(yHist[i]);
            } else {
                sb.append(yHistFloat[i]);
            }
            sb.append(") ");
            
        }
        
        sb.append("]\n");
        
        if (xErrors != null) {
            
            sb.append("histogram errors=[");
            
            for (int i = 0; i < xErrors.length; i++) {

                sb.append("(").append(xErrors[i]).append(", ");

                sb.append(yErrors[i]).append(") ");
            }
            
            sb.append("]\n");
        }      
        
        return sb.toString();
    }
}
