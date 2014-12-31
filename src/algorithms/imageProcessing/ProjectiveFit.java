package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class ProjectiveFit {
   
    private int numberOfPoints = 0;
    private double meanDistFromModel = Double.MAX_VALUE;
    private double stdDevOfMean = Double.MAX_VALUE;
    private double[][] projection = null;
    private double tolerance = Double.MAX_VALUE;
    private Object data = null;

    /**
     * @return the numberOfPoints
     */
    public int getNumberOfPoints() {
        return numberOfPoints;
    }

    /**
     * @param numberOfPoints the numberOfPoints to set
     */
    public void setNumberOfPoints(int numberOfPoints) {
        this.numberOfPoints = numberOfPoints;
    }

    /**
     * @return the meanDistFromModel
     */
    public double getMeanDistFromModel() {
        return meanDistFromModel;
    }

    /**
     * @param meanDistFromModel the meanDistFromModel to set
     */
    public void setMeanDistFromModel(double meanDistFromModel) {
        this.meanDistFromModel = meanDistFromModel;
    }

    /**
     * @return the stdDevOfMean
     */
    public double getStdDevOfMean() {
        return stdDevOfMean;
    }

    /**
     * @param stdDevOfMean the stdDevOfMean to set
     */
    public void setStdDevOfMean(double stdDevOfMean) {
        this.stdDevOfMean = stdDevOfMean;
    }

    /**
     * @return the projection
     */
    public double[][] getProjection() {
        return projection;
    }

    /**
     * @param projection the projection to set
     */
    public void setProjection(double[][] projection) {
        this.projection = projection;
    }

    /**
     * @return the tolerance
     */
    public double getTolerance() {
        return tolerance;
    }

    /**
     * @param tolerance the tolerance to set
     */
    public void setTolerance(double tolerance) {
        this.tolerance = tolerance;
    }

    /**
     * @return the data
     */
    public Object getData() {
        return data;
    }

    /**
     * @param data the data to set
     */
    public void setData(Object data) {
        this.data = data;
    }
}
