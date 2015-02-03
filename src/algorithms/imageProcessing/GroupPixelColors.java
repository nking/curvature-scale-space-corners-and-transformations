package algorithms.imageProcessing;

import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class GroupPixelColors {
    
    private float avgContrast;
    private float avgColorDiff;
    private double stdDevContrast;
    private double stdDevColorDiff;
    private float avgRed;
    private float avgGreen;
    private float avgBlue;
    private double stdDevRed;
    private double stdDevGreen;
    private double stdDevBlue;
    
    private float[] avgYUV;
    
    public GroupPixelColors(final Set<PixelColors> colors) {
        
        if (colors == null) {
            throw new IllegalArgumentException("colors cannot be null");
        }
        
        initialize(colors);
    }
    
    public GroupPixelColors(final Set<PairInt> points, Image colorImage,
        int xOffset, int yOffset) {
        
        if (points == null) {
            throw new IllegalArgumentException("points cannot be null");
        }
        if (colorImage == null) {
            throw new IllegalArgumentException("colorImage cannot be null");
        }
        
        Set<PixelColors> colors = new HashSet<PixelColors>();
        
        for (PairInt p : points) {
            
            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            
            int r = colorImage.getR(x, y);
            int g = colorImage.getG(x, y);
            int b = colorImage.getB(x, y);
            
            PixelColors pixClr = new PixelColors(r, g, b);
            
            colors.add(pixClr);
        }
        
        initialize(colors);
    }
    
    private void initialize(final Set<PixelColors> colors) {
        
        if (colors == null) {
            throw new IllegalArgumentException("colors cannot be null");
        }
        
        float n = colors.size();
        
        float sumRed = 0;
        float sumGreen = 0;
        float sumBlue = 0;
        
        int i = 0;
        for (PixelColors pixColor : colors) {
            int r = pixColor.getRed();
            int g = pixColor.getGreen();
            int b = pixColor.getBlue();
            sumRed += r;
            sumGreen += g;
            sumBlue += b;
            i++;
        }
        
        avgRed = sumRed/n;
        avgGreen = sumGreen/n;
        avgBlue = sumBlue/n;
         
        float[][] m = new float[3][];
        m[0] = new float[]{0.256f, 0.504f, 0.098f};
        m[1] = new float[]{-0.148f, -0.291f, 0.439f};
        m[2] = new float[]{0.439f, -0.368f, -0.072f};
        
        avgYUV = MatrixUtil.multiply(m, new float[]{avgRed, avgGreen, avgBlue});
        
        float sumContrast = 0;
        float sumColorDiff = 0;
        double sumStdDevRed = 0;
        double sumStdDevGreen = 0;
        double sumStdDevBlue = 0;
        float[] contrast = new float[colors.size()];
        float[] colorDiff = new float[colors.size()];
       
        i = 0;
        for (PixelColors pixColor : colors) {
            int r = pixColor.getRed();
            int g = pixColor.getGreen();
            int b = pixColor.getBlue();
            
            float[] yuv = MatrixUtil.multiply(m, new float[]{r, g, b});
            
            float diffR = r - avgRed;
            float diffG = g - avgGreen;
            float diffB = b - avgBlue;
            
            sumStdDevRed += (diffR * diffR);
            sumStdDevGreen += (diffG * diffG);
            sumStdDevBlue += (diffB * diffB);
            
            contrast[i] = (avgYUV[0] - yuv[0])/yuv[0];
            sumContrast += contrast[i];
            
            colorDiff[i] = (float)Math.sqrt(diffR*diffR + diffG*diffG + 
                diffB*diffB);
            sumColorDiff += colorDiff[i];
            i++;
        }
        
        avgContrast = sumContrast/n;
        avgColorDiff = sumColorDiff/n;
        
        stdDevRed = (n > 1) ? 
            Math.sqrt(sumStdDevRed/(n - 1)) : Double.POSITIVE_INFINITY;
        
        stdDevGreen = (n > 1) ? 
            Math.sqrt(sumStdDevGreen/(n - 1)) : Double.POSITIVE_INFINITY;
        stdDevBlue = (n > 1) ? 
            Math.sqrt(sumStdDevBlue/(n - 1)) : Double.POSITIVE_INFINITY;
        
        double sumStdDevContrast = 0;
        double sumStdDevColorDiff = 0;
        
        for (i = 0; i < colors.size(); i++) {
            float diffContrast = contrast[i] - avgContrast;
            sumStdDevContrast += (diffContrast * diffContrast);
            
            float diffColorDiff = colorDiff[i] - avgColorDiff;
            sumStdDevColorDiff += (diffColorDiff * diffColorDiff);            
        }
        
        stdDevContrast = (n > 1) ? 
            Math.sqrt(sumStdDevContrast/(n - 1)) : Double.POSITIVE_INFINITY;
        stdDevColorDiff = (n > 1) ? 
            Math.sqrt(sumStdDevColorDiff/(n - 1)) : Double.POSITIVE_INFINITY;
    }
    
    public float calcContrastToOther(int otherRed, int otherGreen, 
        int otherBlue) {
        
        float[][] m = new float[3][];
        m[0] = new float[]{0.256f, 0.504f, 0.098f};
        m[1] = new float[]{-0.148f, -0.291f, 0.439f};
        m[2] = new float[]{0.439f, -0.368f, -0.072f};
        
        float[] yuvOther = MatrixUtil.multiply(m, new float[]{otherRed, 
            otherGreen, otherBlue});
        
        float contrast = (avgYUV[0] - yuvOther[0])/yuvOther[0];
        
        return contrast;
    }
    
    public float calcColorDiffToOther(int otherRed, int otherGreen, 
        int otherBlue) {
        
        float rDiff = otherRed - avgRed;
        float gDiff = otherGreen - avgGreen;
        float bDiff = otherBlue - avgBlue;
        
        float dist = (float)Math.sqrt(rDiff*rDiff + gDiff*gDiff + bDiff*bDiff);
        
        return dist;
    }    

    /**
     * @return the averageContrast
     */
    public float getAvgContrast() {
        return avgContrast;
    }

    /**
     * @return the averageColorDifference
     */
    public float getAvgColorDiff() {
        return avgColorDiff;
    }

    /**
     * @return the standardDeviationContrast
     */
    public double getStdDevContrast() {
        return stdDevContrast;
    }

    /**
     * @return the standardDeviationColorDifference
     */
    public double getStdDevColorDiff() {
        return stdDevColorDiff;
    }

    /**
     * @return the averageRed
     */
    public float getAvgRed() {
        return avgRed;
    }

    /**
     * @return the averageGreen
     */
    public float getAvgGreen() {
        return avgGreen;
    }

    /**
     * @return the averageBlue
     */
    public float getAvgBlue() {
        return avgBlue;
    }

    /**
     * @return the standardDeviationRed
     */
    public double getStdDevRed() {
        return stdDevRed;
    }

    /**
     * @return the standardDeviationGreen
     */
    public double getStdDevGreen() {
        return stdDevGreen;
    }

    /**
     * @return the standardDeviationBlue
     */
    public double getStdDevBlue() {
        return stdDevBlue;
    }

    /**
     * @return the averageYUV
     */
    public float[] getAverageYUV() {
        return avgYUV;
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        sb.append(" avgRed=").append(avgRed).append("\n")
        .append(" avgGreen=").append(avgGreen).append("\n")
        .append(" avgBlue=").append(avgBlue).append("\n")
        .append("avgContrast=").append(avgContrast).append("\n")
        .append(" avgColorDiff=").append(avgColorDiff).append("\n")
        .append(" stdDevRed=").append(stdDevRed).append("\n")
        .append(" stdDevGreen=").append(stdDevGreen).append("\n")
        .append(" stdDevBlue=").append(stdDevBlue).append("\n")
        .append(" stdDevContrast=").append(stdDevContrast).append("\n")
        .append(" stdDevColorDiff=").append(stdDevColorDiff).append("\n");
       
        return sb.toString();
    }
    
}
