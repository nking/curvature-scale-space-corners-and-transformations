package algorithms.imageProcessing;

import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class GroupPixelRGB0 {
    
    protected float avgRed;
    protected float avgGreen;
    protected float avgBlue;
    

    public GroupPixelRGB0() {
    }

    protected void calculateColors(final Set<PairInt> points, ImageExt colorImage, int xOffset, int yOffset) {
        float n = points.size();
        float sumRed = 0;
        float sumGreen = 0;
        float sumBlue = 0;
        int i = 0;
        for (PairInt p : points) {
            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            int idx = colorImage.getInternalIndex(x, y);
            sumRed += colorImage.getR(idx);
            sumGreen += colorImage.getG(idx);
            sumBlue += colorImage.getB(idx);
            i++;
        }
        avgRed = sumRed / n;
        avgGreen = sumGreen / n;
        avgBlue = sumBlue / n;
    }
    
    protected void calculateColors(final TIntSet pointIndexes, 
        ImageExt colorImage, int xOffset, int yOffset) {
        
        float n = pointIndexes.size();
        float sumRed = 0;
        float sumGreen = 0;
        float sumBlue = 0;
        
        TIntIterator iter = pointIndexes.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int x = colorImage.getCol(pixIdx) + xOffset;
            int y = colorImage.getRow(pixIdx) + yOffset;
            sumRed += colorImage.getR(pixIdx);
            sumGreen += colorImage.getG(pixIdx);
            sumBlue += colorImage.getB(pixIdx);
        }
        avgRed = sumRed / n;
        avgGreen = sumGreen / n;
        avgBlue = sumBlue / n;
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

    public int getAvgRGB() {
        float avg = (avgRed + avgGreen + avgBlue) / 3.f;
        return Math.round(avg);
    }
   
    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        sb.append(" avgRed=").append(avgRed).append("\n")
        .append(" avgGreen=").append(avgGreen).append("\n")
        .append(" avgBlue=").append(avgBlue).append("\n")
        ;
       
        return sb.toString();
    }
}
