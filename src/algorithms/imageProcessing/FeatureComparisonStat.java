package algorithms.imageProcessing;

import algorithms.util.PairInt;

/**
 *
 * @author nichole
 */
public class FeatureComparisonStat implements Comparable<FeatureComparisonStat> {

    private PairInt img1Point = null;
    private PairInt img2Point = null;
    private float sumIntensitySqDiff = Float.POSITIVE_INFINITY;
    private float sumGradientSqDiff = Float.POSITIVE_INFINITY;
    private float img2PointIntensityErr = Float.POSITIVE_INFINITY;
    private float img2PointGradientErr = Float.POSITIVE_INFINITY;
    private float img1RotInDegrees = Float.POSITIVE_INFINITY;
    private float img2RotInDegrees = Float.POSITIVE_INFINITY;

    /**
     * @return the img1Point
     */
    public PairInt getImg1Point() {
        return img1Point;
    }

    /**
     * @param img1Point the img1Point to set
     */
    public void setImg1Point(PairInt img1Point) {        
        this.img1Point = img1Point;
    }

    /**
     * @return the img2Point
     */
    public PairInt getImg2Point() {
        return img2Point;
    }

    /**
     * @param img2Point the img2Point to set
     */
    public void setImg2Point(PairInt img2Point) {        
        this.img2Point = img2Point;
    }

    /**
     * @return the sum of the square of differences of intensity
     */
    public float getSumIntensitySqDiff() {
        return sumIntensitySqDiff;
    }

    /**
     * @param sumSqDiff the sum of the square of differences of intensity
     */
    public void setSumIntensitySqDiff(float sumSqDiff) {
        this.sumIntensitySqDiff = sumSqDiff;
    }

    /**
     * @return the sum of the square of differences of intensity
     */
    public float getSumGradientSqDiff() {
        return sumGradientSqDiff;
    }

    /**
     * @param sumSqDiff the sum of the square of differences of intensity
     */
    public void setSumGradientSqDiff(float sumSqDiff) {
        this.sumGradientSqDiff = sumSqDiff;
    }

    /**
     * @return the img2PointIntensityErr
     */
    public float getImg2PointIntensityErr() {
        return img2PointIntensityErr;
    }

    /**
     * @param img2PointIntensityErr the img2PointIntensityErr to set
     */
    public void setImg2PointIntensityErr(float img2PointIntensityErr) {
        this.img2PointIntensityErr = img2PointIntensityErr;
    }
    
    /**
     * @return the img2PointGradientErr
     */
    public float getImg2PointGradientErr() {
        return img2PointGradientErr;
    }

    /**
     * @param img2PointGradientErr the img2PointIntensityErr to set
     */
    public void setImg2PointGradientErr(float img2PointGradientErr) {
        this.img2PointGradientErr = img2PointGradientErr;
    }
    
    /**
     * compare object other to this instance and return -1 if this instance
     * is "better", 0 if they are equal, else 1 if other is "better".
     * @param other
     * @return 
     */
    @Override
    public int compareTo(FeatureComparisonStat other) {
    
        if (other == null) {
            return -1;
        }
        if (Float.isInfinite(other.getSumIntensitySqDiff())) {
            return -1;
        }
        
        throw new UnsupportedOperationException("not yet implemented");
        
    }

    /**
     * @return the img1RotInDegrees
     */
    public float getImg1RotInDegrees() {
        return img1RotInDegrees;
    }

    /**
     * @param img1RotInDegrees the img1RotInDegrees to set
     */
    public void setImg1RotInDegrees(float img1RotInDegrees) {
        this.img1RotInDegrees = img1RotInDegrees;
    }
    
    /**
     * @return the img2RotInDegrees
     */
    public float getImg2RotInDegrees() {
        return img2RotInDegrees;
    }

    /**
     * @param img2RotInDegrees the img2RotInDegrees to set
     */
    public void setImg2RotInDegrees(float img2RotInDegrees) {
        this.img2RotInDegrees = img2RotInDegrees;
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        sb.append("p1=").append(img1Point.toString()).append(" ")
            .append("p2=").append(img2Point.toString())
            .append(String.format(" ssdIntensity=%.4f", sumIntensitySqDiff))
            .append(String.format(" err2Intensity=%.4f", img2PointIntensityErr))
            .append(String.format(" ssdGradient=%.4f", sumGradientSqDiff))
            .append(String.format(" err2Gradient=%.4f", img2PointGradientErr))
            .append(String.format(" rot1D=%.4f", img1RotInDegrees))
            .append(String.format(" rot2D=%.4f", img2RotInDegrees))
        ;
        
        return sb.toString();
    }

}
