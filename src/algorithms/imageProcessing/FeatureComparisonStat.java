package algorithms.imageProcessing;

import algorithms.util.PairInt;

/**
 *
 * @author nichole
 */
public class FeatureComparisonStat implements Comparable<FeatureComparisonStat> {

    private PairInt img1Point;
    private PairInt img2Point;
    private float sumSqDiff;
    private float img2PointErr;
    private float img1RotInDegrees = Float.POSITIVE_INFINITY;

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
     * @return the sumSqDiff
     */
    public float getSumSqDiff() {
        return sumSqDiff;
    }

    /**
     * @param sumSqDiff the sumSqDiff to set
     */
    public void setSumSqDiff(float sumSqDiff) {
        this.sumSqDiff = sumSqDiff;
    }

    /**
     * @return the img2PointErr
     */
    public float getImg2PointErr() {
        return img2PointErr;
    }

    /**
     * @param img2PointErr the img2PointErr to set
     */
    public void setImg2PointErr(float img2PointErr) {
        this.img2PointErr = img2PointErr;
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
        if (Float.isInfinite(other.getSumSqDiff())) {
            return -1;
        }
        
        //TODO: consider an eps within which to use sumSqDiff/err
        
        float div = sumSqDiff/img2PointErr;
        float divOther = other.getSumSqDiff()/other.getImg2PointErr();
        
        if (other.getSumSqDiff() > other.getImg2PointErr()) {
            if (sumSqDiff > img2PointErr) {
                if (div < divOther) {
                    return -1;
                } else if (div > divOther) {
                    return 1;
                }
                return 0;
            }
            return -1;
        } else if (sumSqDiff > img2PointErr) {
            return 1;
        }
        
        if (sumSqDiff < other.getSumSqDiff()) {
            return -1;
        } else if (sumSqDiff == other.getSumSqDiff()) {
            return 0;
        } else {
            return 1;
        }
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

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        sb.append("p1=").append(img1Point.toString()).append(" ")
            .append("p2=").append(img2Point.toString()).append("\n")
            .append("sumSqDiff=").append(sumSqDiff)
            .append(" err2=").append(img2PointErr)
            .append(" img1RotInDegrees=").append(img1RotInDegrees);
        
        return sb.toString();
    }

}
