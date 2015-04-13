package algorithms.imageProcessing;

import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.util.PairInt;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class GroupPixelColors extends GroupPixelRGB {
    
    protected float avgContrast;
    protected float stdDevContrast;
    
    protected float avgLuma;
    protected float stdDevLuma;
    
    protected float avgCIEX;
    protected float stdDevCIEX;
    protected float avgCIEY;
    protected float stdDevCIEY;
    
    /**
     * calculate average colors and standard deviations of the average for
     * the points in the colorImage.  The xOffset and yOffset are the offsets
     * of the points from the colorImage reference frame.
     * @param points
     * @param colorImage
     * @param xOffset
     * @param yOffset 
     */
    public GroupPixelColors(final Set<PairInt> points, ImageExt colorImage,
        int xOffset, int yOffset) {
        
        super(points, colorImage, xOffset, yOffset);
    }
    
    @Override
    protected void calculateColors(final Set<PairInt> points, ImageExt 
        colorImage, int xOffset, int yOffset) {
        
        super.calculateColors(points, colorImage, xOffset, yOffset);
        
        float n = points.size();
        
        float sumCIEX = 0;
        float sumCIEY = 0;
        
        float sumLuma = 0;
        
        int i = 0;
        
        for (PairInt p : points) {
            
            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            int idx = colorImage.getInternalIndex(x, y);
                        
            sumCIEX += colorImage.getCIEX(idx);
            sumCIEY += colorImage.getCIEY(idx);
            
            sumLuma += colorImage.getLuma(idx);
            
            i++;
        }
    
        avgCIEX = sumCIEX/n;
        avgCIEY = sumCIEY/n;
        avgLuma = sumLuma/n;
        
        float sumContrast = 0;
        float[] contrast = new float[points.size()];
       
        double sumStdDevCIEX = 0;
        double sumStdDevCIEY = 0;
        double sumStdDevLuma = 0;
        
        i = 0;
        for (PairInt p : points) {
            
            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            int idx = colorImage.getInternalIndex(x, y);
            
            float cieX = colorImage.getCIEX(idx);
            float cieY = colorImage.getCIEY(idx);
            float luma = colorImage.getLuma(idx);
                        
            float diffCIEX = cieX - avgCIEX;
            float diffCIEY = cieY - avgCIEY;
            float diffLuma = luma - avgLuma;
            sumStdDevCIEX += (diffCIEX * diffCIEX);
            sumStdDevCIEY += (diffCIEY * diffCIEY);
            sumStdDevLuma += (diffLuma * diffLuma);
            
            contrast[i] = (avgLuma - luma)/luma;
            sumContrast += contrast[i];
            
            i++;
        }
        
        avgContrast = sumContrast/n;
       
        stdDevCIEX = (n > 1) ? 
            (float)Math.sqrt(sumStdDevCIEX/(n - 1)) : Float.POSITIVE_INFINITY;
        stdDevCIEY = (n > 1) ? 
            (float)Math.sqrt(sumStdDevCIEY/(n - 1)) : Float.POSITIVE_INFINITY;
        
        stdDevLuma = (n > 1) ? 
            (float)Math.sqrt(sumStdDevLuma/(n - 1)) : Float.POSITIVE_INFINITY;
        
        double sumStdDevContrast = 0;
        
        for (i = 0; i < points.size(); i++) {
            float diffContrast = contrast[i] - avgContrast;
            sumStdDevContrast += (diffContrast * diffContrast);            
        }
        
        stdDevContrast = (n > 1) ? 
            (float)Math.sqrt(sumStdDevContrast/(n - 1)) : Float.POSITIVE_INFINITY;
    }
    
    public float calcContrastToOther(int otherRed, int otherGreen, 
        int otherBlue) {
        
        float[][] m = new float[3][];
        m[0] = new float[]{0.256f, 0.504f, 0.098f};
        m[1] = new float[]{-0.148f, -0.291f, 0.439f};
        m[2] = new float[]{0.439f, -0.368f, -0.072f};
        
        float[] yuvOther = MatrixUtil.multiply(m, new float[]{otherRed, 
            otherGreen, otherBlue});
        
        float contrast = (avgLuma - yuvOther[0])/yuvOther[0];
        
        return contrast;
    }
    
    public float calcContrastToOther(float otherLuma) {
        
        float contrast = (avgLuma - otherLuma)/otherLuma;
        
        return contrast;
    }
    
    /**
     * @return the averageContrast
     */
    public float getAvgContrast() {
        return avgContrast;
    }

    /**
     * @return the standardDeviationContrast
     */
    public double getStdDevContrast() {
        return stdDevContrast;
    }

    /**
     * @return the averageLuma
     */
    public float getAverageLuma() {
        return avgLuma;
    }
    
    /**
     * @return the avgCIEX
     */
    public double getAverageCIEX() {
        return avgCIEX;
    }
    /**
     * @return the avgCIEY
     */
    public double getAverageCIEY() {
        return avgCIEY;
    }
    /**
     * @return the stdDevCIEX
     */
    public double getStdDevCIEX() {
        return stdDevCIEX;
    }
    /**
     * @return the stdDevCIEY
     */
    public double getStdDevCIEY() {
        return stdDevCIEY;
    }
    /**
     * @return the stdDevLuma
     */
    public double getStdDevLuma() {
        return stdDevLuma;
    }

    @Override
    public String toString() {
        
        String str = super.toString();
        
        StringBuilder sb = new StringBuilder(str);
        sb.append("avgContrast=").append(avgContrast).append("\n")
        .append(" avgCIEX=").append(avgCIEX).append("\n")
        .append(" avgCIEY=").append(avgCIEY).append("\n")
        .append(" avgLuma=").append(avgCIEY).append("\n")
        .append(" stdDevContrast=").append(stdDevContrast).append("\n")
        .append(" stdDevCIEX=").append(stdDevCIEX).append("\n")
        .append(" stdDevCIEY=").append(stdDevCIEY).append("\n")
        .append(" stdDevLuma=").append(stdDevLuma).append("\n")
        ;
       
        return sb.toString();
    }
    
}
