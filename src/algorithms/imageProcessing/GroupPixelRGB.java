package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class GroupPixelRGB extends GroupPixelRGB0 {
    
    protected float stdDevRed;
    protected float stdDevGreen;
    protected float stdDevBlue;
    protected float avgColorDiff;
    protected float stdDevColorDiff;
    
    public GroupPixelRGB(final Set<PairInt> points, ImageExt colorImage,
        int xOffset, int yOffset) {
        
        if (points == null) {
            throw new IllegalArgumentException("points cannot be null");
        }
        if (colorImage == null) {
            throw new IllegalArgumentException("colorImage cannot be null");
        }
        
        calculateColors(points, colorImage, xOffset, yOffset);
    }
    
    protected void calculateColors(final Set<PairInt> points, 
        ImageExt colorImage, int xOffset, int yOffset) {
        
        super.calculateColors(points, colorImage, xOffset, yOffset);
        
        float n = points.size();
        int i = 0;
        
        double sumStdDevRed = 0;
        double sumStdDevGreen = 0;
        double sumStdDevBlue = 0;
        double sumColorDiff = 0;
        float[] colorDiff = new float[points.size()];
        for (PairInt p : points) {
            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            int idx = colorImage.getInternalIndex(x, y);
            int r = colorImage.getR(idx);
            int g = colorImage.getG(idx);
            int b = colorImage.getB(idx);
            float diffR = r - avgRed;
            float diffG = g - avgGreen;
            float diffB = b - avgBlue;
            sumStdDevRed += (diffR * diffR);
            sumStdDevGreen += (diffG * diffG);
            sumStdDevBlue += (diffB * diffB);
            colorDiff[i] = (float) Math.sqrt(diffR * diffR + diffG * diffG + diffB * diffB);
            sumColorDiff += colorDiff[i];
            i++;
        }
        avgColorDiff = (float) sumColorDiff / n;
        stdDevRed = (n > 1) ? (float) Math.sqrt(sumStdDevRed / (n - 1)) : Float.POSITIVE_INFINITY;
        stdDevGreen = (n > 1) ? (float) Math.sqrt(sumStdDevGreen / (n - 1)) : Float.POSITIVE_INFINITY;
        stdDevBlue = (n > 1) ? (float) Math.sqrt(sumStdDevBlue / (n - 1)) : Float.POSITIVE_INFINITY;
        double sumStdDevColorDiff = 0;
        for (i = 0; i < points.size(); i++) {
            float diffColorDiff = colorDiff[i] - avgColorDiff;
            sumStdDevColorDiff += (diffColorDiff * diffColorDiff);
        }
        stdDevColorDiff = (n > 1) ? (float) Math.sqrt(sumStdDevColorDiff / (n - 1)) : Float.POSITIVE_INFINITY;
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
     * @return the standardDeviationColorDifference
     */
    public double getStdDevColorDiff() {
        return stdDevColorDiff;
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
     * @return the averageColorDifference
     */
    public float getAvgColorDiff() {
        return avgColorDiff;
    }
    
    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        sb.append(" avgRed=").append(avgRed).append("\n")
        .append(" avgGreen=").append(avgGreen).append("\n")
        .append(" avgBlue=").append(avgBlue).append("\n")
        .append(" avgColorDiff=").append(avgColorDiff).append("\n")
        .append(" stdDevRed=").append(stdDevRed).append("\n")
        .append(" stdDevGreen=").append(stdDevGreen).append("\n")
        .append(" stdDevBlue=").append(stdDevBlue).append("\n")
        .append(" stdDevColorDiff=").append(stdDevColorDiff).append("\n")
        ;
       
        return sb.toString();
    }
    
}
