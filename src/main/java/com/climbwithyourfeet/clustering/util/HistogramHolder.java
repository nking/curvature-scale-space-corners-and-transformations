package com.climbwithyourfeet.clustering.util;

/**
 *
 * @author nichole
 */
public class HistogramHolder {
    
    private static final long serialVersionUID = -7105371701539621066L;

    protected float[] xHist = null;
    protected int[] yHist = null;
    protected float[] yHistFloat = null;
    protected float[] yErrors = null;
    protected float[] xErrors = null;

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
