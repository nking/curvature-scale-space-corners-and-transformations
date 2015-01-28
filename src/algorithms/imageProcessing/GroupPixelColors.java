package algorithms.imageProcessing;

import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class GroupPixelColors {
    
    private float averageContrast;
    private float averageColorDifference;
    private double standardDeviationContrast;
    private double standardDeviationColorDifference;
    private float averageRed;
    private float averageGreen;
    private float averageBlue;
    private double standardDeviationRed;
    private double standardDeviationGreen;
    private double standardDeviationBlue;
    
    private float[] averageYUV;
    
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
        
        averageRed = sumRed/n;
        averageGreen = sumGreen/n;
        averageBlue = sumBlue/n;
         
        float[][] m = new float[3][];
        m[0] = new float[]{0.256f, 0.504f, 0.098f};
        m[1] = new float[]{-0.148f, -0.291f, 0.439f};
        m[2] = new float[]{0.439f, -0.368f, -0.072f};
        
        averageYUV = MatrixUtil.multiply(m, new float[]{averageRed, 
            averageGreen, averageBlue});
        
        float sumContrast = 0;
        float sumColorDifference = 0;
        double sumStandardDeviationRed = 0;
        double sumStandardDeviationGreen = 0;
        double sumStandardDeviationBlue = 0;
        float[] contrast = new float[colors.size()];
        float[] colorDiff = new float[colors.size()];
       
        i = 0;
        for (PixelColors pixColor : colors) {
            int r = pixColor.getRed();
            int g = pixColor.getGreen();
            int b = pixColor.getBlue();
            
            float[] yuv = MatrixUtil.multiply(m, new float[]{r, g, b});
            
            float diffR = r - averageRed;
            float diffG = g - averageGreen;
            float diffB = b - averageBlue;
            
            sumStandardDeviationRed += (diffR * diffR);
            sumStandardDeviationGreen += (diffG * diffG);
            sumStandardDeviationBlue += (diffB * diffB);
            
            contrast[i] = (averageYUV[0] - yuv[0])/yuv[0];
            sumContrast += contrast[i];
            
            colorDiff[i] = (float)Math.sqrt(diffR*diffR + diffG*diffG + diffB*diffB);
            sumColorDifference += colorDiff[i];
            i++;
        }
        
        averageContrast = sumContrast/n;
        averageColorDifference = sumColorDifference/n;
        
        standardDeviationRed = (n > 1) ? 
            Math.sqrt(sumStandardDeviationRed/(n - 1)) : Double.POSITIVE_INFINITY;
        
        standardDeviationGreen = (n > 1) ? 
            Math.sqrt(sumStandardDeviationGreen/(n - 1)) : Double.POSITIVE_INFINITY;
        standardDeviationBlue = (n > 1) ? 
            Math.sqrt(sumStandardDeviationBlue/(n - 1)) : Double.POSITIVE_INFINITY;
        
        double sumStandardDeviationContrast = 0;
        double sumStandardDeviationColorDifference = 0;
        
        for (i = 0; i < colors.size(); i++) {
            float diffContrast = contrast[i] - averageContrast;
            sumStandardDeviationContrast += (diffContrast * diffContrast);
            
            float diffColorDiff = colorDiff[i] - averageColorDifference;
            sumStandardDeviationColorDifference += (diffColorDiff * diffColorDiff);            
        }
        
        standardDeviationContrast = (n > 1) ? 
            Math.sqrt(sumStandardDeviationContrast/(n - 1)) : Double.POSITIVE_INFINITY;
        standardDeviationColorDifference = (n > 1) ? 
            Math.sqrt(sumStandardDeviationColorDifference/(n - 1)) : Double.POSITIVE_INFINITY;
    }
    
    public float calculateContrastToOther(int otherRed, int otherGreen, int otherBlue) {
        
        float[][] m = new float[3][];
        m[0] = new float[]{0.256f, 0.504f, 0.098f};
        m[1] = new float[]{-0.148f, -0.291f, 0.439f};
        m[2] = new float[]{0.439f, -0.368f, -0.072f};
        
        float[] yuvOther = MatrixUtil.multiply(m, new float[]{otherRed, otherGreen, otherBlue});
        
        float contrast = (averageYUV[0] - yuvOther[0])/yuvOther[0];
        
        return contrast;
    }
    
    public float calculateColorDifferenceToOther(int otherRed, int otherGreen, int otherBlue) {
        
        float rDiff = otherRed - averageRed;
        float gDiff = otherGreen - averageGreen;
        float bDiff = otherBlue - averageBlue;
        
        float dist = (float)Math.sqrt(rDiff*rDiff + gDiff*gDiff + bDiff*bDiff);
        
        return dist;
    }    

    /**
     * @return the averageContrast
     */
    public float getAverageContrast() {
        return averageContrast;
    }

    /**
     * @return the averageColorDifference
     */
    public float getAverageColorDifference() {
        return averageColorDifference;
    }

    /**
     * @return the standardDeviationContrast
     */
    public double getStandardDeviationContrast() {
        return standardDeviationContrast;
    }

    /**
     * @return the standardDeviationColorDifference
     */
    public double getStandardDeviationColorDifference() {
        return standardDeviationColorDifference;
    }

    /**
     * @return the averageRed
     */
    public float getAverageRed() {
        return averageRed;
    }

    /**
     * @return the averageGreen
     */
    public float getAverageGreen() {
        return averageGreen;
    }

    /**
     * @return the averageBlue
     */
    public float getAverageBlue() {
        return averageBlue;
    }

    /**
     * @return the standardDeviationRed
     */
    public double getStandardDeviationRed() {
        return standardDeviationRed;
    }

    /**
     * @return the standardDeviationGreen
     */
    public double getStandardDeviationGreen() {
        return standardDeviationGreen;
    }

    /**
     * @return the standardDeviationBlue
     */
    public double getStandardDeviationBlue() {
        return standardDeviationBlue;
    }

    /**
     * @return the averageYUV
     */
    public float[] getAverageYUV() {
        return averageYUV;
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        sb.append(" averageRed=").append(averageRed).append("\n")
        .append(" averageGreen=").append(averageGreen).append("\n")
        .append(" averageBlue=").append(averageBlue).append("\n")
        .append("averageContrast=").append(averageContrast).append("\n")
        .append(" averageColorDifference=").append(averageColorDifference).append("\n")
        .append(" standardDeviationRed=").append(standardDeviationRed).append("\n")
        .append(" standardDeviationGreen=").append(standardDeviationGreen).append("\n")
        .append(" standardDeviationBlue=").append(standardDeviationBlue).append("\n")
        .append(" standardDeviationContrast=").append(standardDeviationContrast).append("\n")
        .append(" standardDeviationColorDifference=").append(standardDeviationColorDifference).append("\n");
       
        return sb.toString();
    }
    
}
