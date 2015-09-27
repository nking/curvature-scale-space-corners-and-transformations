package com.climbwithyourfeet.clustering.util;

import com.climbwithyourfeet.clustering.CriticalDensitySolver;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.util.logging.Level;
import java.util.logging.Logger;

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

    /**
     * if a specific class exists in the classpath, this will plot the
     * histogram to an html file.
     * @param label
     * @param outputFileSuffix
     * @return
     */
    public String plotHistogram(String label, String outputFileSuffix) {
                
        try {
            // only using a class if it exists already in classpath. 
            // class not imported to avoid a dependency in the packaged jar
            ClassLoader cls = this.getClass().getClassLoader();
            Class<?> plotClass = cls.loadClass("algorithms.util.PolygonAndPointPlotter");
         
            Class<?>[] argTypes0 = new Class<?>[]{float.class, float.class, 
                float.class, float.class, float[].class, float[].class, 
                float[].class, float[].class, String.class};
            
            Method method0 = plotClass.getMethod("addPlot", argTypes0);
            
            Constructor constructor = null;
            for (Constructor c : plotClass.getConstructors()) {
                if (c.getParameterCount() == 0) {
                    constructor = c;
                    break;
                }
            }
            if (constructor == null) {
                return null;
            }
            
            Object plotterObj = constructor.newInstance();

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
            
            Object[] args0 = new Object[]{xMin, xMax, yMin, yMax, xh, yh, xh, yh, label};

            method0.invoke(plotterObj, args0);
            
            Class<?>[] argTypes1 = new Class<?>[]{String.class};
            
            Method method1 = plotClass.getMethod("writeFile", argTypes1);
            
            Object[] args1 = new Object[]{outputFileSuffix};
            
            String filePath = (String) method1.invoke(plotterObj, args1);
        
            return filePath;
            
        } catch (Exception ex) {
            Logger.getLogger(CriticalDensitySolver.class.getName()).
                log(Level.SEVERE, null, ex);
        }
        
        return null;
    }
    
}
