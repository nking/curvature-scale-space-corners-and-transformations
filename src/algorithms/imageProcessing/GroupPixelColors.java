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
    
    private double averageContrast;
    private double averageColorDifference;
    private double standardDeviationContrast;
    private double standardDeviationColorDifference;
    private double averageRed;
    private double averageGreen;
    private double averageBlue;
    private double standardDeviationRed;
    private double standardDeviationGreen;
    private double standardDeviationBlue;
    
    private double[] averageYUV;
    
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
        
        double n = colors.size();
        
        double sumRed = 0;
        double sumGreen = 0;
        double sumBlue = 0;
        
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
         
        double[][] m = new double[3][];
        m[0] = new double[]{0.256, 0.504, 0.098};
        m[1] = new double[]{-0.148, -0.291, 0.439};
        m[2] = new double[]{0.439, -0.368, -0.072};
        
        averageYUV = MatrixUtil.multiply(m, new double[]{averageRed, 
            averageGreen, averageBlue});
        
        double sumContrast = 0;
        double sumColorDifference = 0;
        double sumStandardDeviationRed = 0;
        double sumStandardDeviationGreen = 0;
        double sumStandardDeviationBlue = 0;
        double[] contrast = new double[colors.size()];
        double[] colorDiff = new double[colors.size()];
       
        i = 0;
        for (PixelColors pixColor : colors) {
            int r = pixColor.getRed();
            int g = pixColor.getGreen();
            int b = pixColor.getBlue();
            
            double[] yuv = MatrixUtil.multiply(m, new double[]{r, g, b});
            
            double diffR = r - averageRed;
            double diffG = g - averageGreen;
            double diffB = b - averageBlue;
            
            sumStandardDeviationRed += (diffR * diffR);
            sumStandardDeviationGreen += (diffG * diffG);
            sumStandardDeviationBlue += (diffB * diffB);
            
            contrast[i] = (averageYUV[0] - yuv[0])/yuv[0];
            sumContrast += contrast[i];
            
            colorDiff[i] = Math.sqrt(diffR*diffR + diffG*diffG + diffB*diffB);
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
            double diffContrast = contrast[i] - averageContrast;
            sumStandardDeviationContrast += (diffContrast * diffContrast);
            
            double diffColorDiff = colorDiff[i] - averageColorDifference;
            sumStandardDeviationColorDifference += (diffColorDiff * diffColorDiff);            
        }
        
        standardDeviationContrast = (n > 1) ? 
            Math.sqrt(sumStandardDeviationContrast/(n - 1)) : Double.POSITIVE_INFINITY;
        standardDeviationColorDifference = (n > 1) ? 
            Math.sqrt(sumStandardDeviationColorDifference/(n - 1)) : Double.POSITIVE_INFINITY;
    }
    
    public double calculateContrastToOther(int otherRed, int otherGreen, int otherBlue) {
        
        double[][] m = new double[3][];
        m[0] = new double[]{0.256, 0.504, 0.098};
        m[1] = new double[]{-0.148, -0.291, 0.439};
        m[2] = new double[]{0.439, -0.368, -0.072};
        
        double[] yuvOther = MatrixUtil.multiply(m, new double[]{otherRed, otherGreen, otherBlue});
        
        double contrast = (averageYUV[0] - yuvOther[0])/yuvOther[0];
        
        return contrast;
    }
    
    public double calculateColorDifferenceToOther(int otherRed, int otherGreen, int otherBlue) {
        
        double rDiff = otherRed - averageRed;
        double gDiff = otherGreen - averageGreen;
        double bDiff = otherBlue - averageBlue;
        
        double dist = Math.sqrt(rDiff*rDiff + gDiff*gDiff + bDiff*bDiff);
        
        return dist;
    }    

    /**
     * @return the averageContrast
     */
    public double getAverageContrast() {
        return averageContrast;
    }

    /**
     * @return the averageColorDifference
     */
    public double getAverageColorDifference() {
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
    public double getAverageRed() {
        return averageRed;
    }

    /**
     * @return the averageGreen
     */
    public double getAverageGreen() {
        return averageGreen;
    }

    /**
     * @return the averageBlue
     */
    public double getAverageBlue() {
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
    public double[] getAverageYUV() {
        return averageYUV;
    }
}
