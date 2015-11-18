package algorithms.imageProcessing;

import algorithms.util.PairInt;

/**
 *
 * @author nichole
 */
public class FeatureComparisonStat implements Comparable<FeatureComparisonStat> {

    private int index1 = -1;
    private int index2 = -1;
    
    private PairInt img1Point = null;
    private PairInt img2Point = null;
    private float sumIntensitySqDiff = Float.POSITIVE_INFINITY;
    private float sumGradientSqDiff = Float.POSITIVE_INFINITY;
    private float sumThetaSqDiff = Float.POSITIVE_INFINITY;
    private float img2PointIntensityErr = Float.POSITIVE_INFINITY;
    private float img2PointGradientErr = Float.POSITIVE_INFINITY;
    private float img2PointThetaErr = Float.POSITIVE_INFINITY;
    private float img1PointRotInDegrees = Float.POSITIVE_INFINITY;
    private float img2PointRotInDegrees = Float.POSITIVE_INFINITY;
    
    /**
     * bin factors for images which these properties were measured from.
     * Note that these fields are recent and not always set by using
     * class so check before using.
     */
    private int binFactor1 = 1;
    private int binFactor2 = 1;

    /**
     * @return the img1Point
     */
    public PairInt getImg1Point() {
        return img1Point;
    }
    
    public void setBinFactor1(int factor) {
        binFactor1 = factor;
    }
    public void setBinFactor2(int factor) {
        binFactor2 = factor;
    }
    public int getBinFactor1() {
        return binFactor1;
    }
    public int getBinFactor2() {
        return binFactor2;
    }

    public void setIndex1(int theIndex) {
        index1 = theIndex;
    }
    
    public void setIndex2(int theIndex) {
        index2 = theIndex;
    }
    
    public int getIndex1() {
        return index1;
    }
    
    public int getIndex2() {
        return index2;
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
        
        boolean hasIntensity = Float.isFinite(sumIntensitySqDiff);
        boolean hasTheta = Float.isFinite(sumThetaSqDiff);
        boolean hasGradient = Float.isFinite(sumGradientSqDiff);
        
        if (hasIntensity && !hasTheta && !hasGradient) {
            return compareToUsingIntensity(other);
        }
        
        throw new UnsupportedOperationException("not yet implemented");
        
    }
    
    public int compareToUsingIntensity(FeatureComparisonStat other) {
        
        // prefer smaller intensities
        
        if (sumIntensitySqDiff < other.sumIntensitySqDiff) {
            return -1;
        } else if (sumIntensitySqDiff > other.sumIntensitySqDiff) {
            return 1;
        }
        
        return 0;
    }

    /**
     * @return the img1PointRotInDegrees
     */
    public float getImg1PointRotInDegrees() {
        return img1PointRotInDegrees;
    }

    /**
     * @param rotationInDegrees the orientation rotation for point 1 region
     */
    public void setImg1PointRotInDegrees(float rotationInDegrees) {
        this.img1PointRotInDegrees = rotationInDegrees;
    }
    
    /**
     * @return the img2PointRotInDegrees
     */
    public float getImg2PointRotInDegrees() {
        return img2PointRotInDegrees;
    }

    /**
     * @param rotationInDegrees the orientation rotation for point 2 region
     */
    public void setImg2PointRotInDegrees(float rotationInDegrees) {
        this.img2PointRotInDegrees = rotationInDegrees;
    }

    /**
     * @return the sumThetaSqDiff
     */
    public float getSumThetaSqDiff() {
        return sumThetaSqDiff;
    }

    /**
     * @param sumSquareDiffs the sum of the square of the differences
     */
    public void setSumThetaSqDiff(float sumSquareDiffs) {
        this.sumThetaSqDiff = sumSquareDiffs;
    }

    /**
     * @return the img2PointThetaErr
     */
    public float getImg2PointThetaErr() {
        return img2PointThetaErr;
    }

    /**
     * @param img2PointThetaErr the img2PointThetaErr to set
     */
    public void setImg2PointThetaErr(float img2PointThetaErr) {
        this.img2PointThetaErr = img2PointThetaErr;
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        sb.append("p1=");
        if (img1Point != null) {
            sb.append(img1Point.toString()).append(" ");
        }
        sb.append("p2=");
        if (img2Point != null) {
            sb.append(img2Point.toString()).append(" ");
        }
        if (Float.isFinite(sumIntensitySqDiff)) {
            sb.append(String.format(" ssdInt=%.4f", sumIntensitySqDiff));
        }
        if (Float.isFinite(img2PointIntensityErr)) {
            sb.append(String.format(" err2Int=%.4f", img2PointIntensityErr));
        }
        if (Float.isFinite(sumGradientSqDiff)) {
            sb.append(String.format(" ssdGrd=%.4f", sumGradientSqDiff));
        }
        if (Float.isFinite(img2PointGradientErr)) {
            sb.append(String.format(" err2Grd=%.4f", img2PointGradientErr));
        }
        if (Float.isFinite(sumThetaSqDiff)) {
            sb.append(String.format(" thtDf=%.4f", sumThetaSqDiff));
        }
        if (Float.isFinite(img2PointThetaErr)) {
            sb.append(String.format(" thtErr=%.4f", img2PointThetaErr));
        }
        if (Float.isFinite(img1PointRotInDegrees)) {
            sb.append(String.format(" rot1D=%.4f", img1PointRotInDegrees));
        } 
        if (Float.isFinite(img2PointRotInDegrees)) {
            sb.append(String.format(" rot2D=%.4f", img2PointRotInDegrees));
        }
        
        return sb.toString();
    }

}