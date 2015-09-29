package algorithms.curves;

import algorithms.misc.MiscMath;

/**
 *
 * @author nichole
 */
public class GEVYFit implements IYFit {

    /**
     *
     */
    protected float[] yfit;

    /**
     *
     */
    protected float[] x;

    /**
     *
     */
    protected float xScale = 1;

    /**
     *
     */
    protected float yScale = 1;

    /**
     *
     */
    protected int xPeakIndex = -1;

    /**
     *
     */
    protected float k;

    /**
     *
     */
    protected float sigma;

    /**
     *
     */
    protected float mu;

    /**
     *
     */
    protected float chiSqSum = Float.MAX_VALUE;

    /**
     *
     */
    protected float chiSqStatistic = Float.MAX_VALUE;
    float kSolutionResolution;
    float sigmaSolutionResolution;
    float muSolutionResolution;

    /**
     *
     */
    protected float yDataErrSq;

    /**
     *
     */
    protected String[] parameterNames = new String[]{
        "k", "sigma", "mu"
    };

    /**
     *
     * @return
     */
    public long approximateMemoryUsed() {

        String arch = System.getProperty("sun.arch.data.model");

        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;

        int nbits = (is32Bit) ? 32 : 64;

        int overheadBytes = 16;

        int intBytes = (is32Bit) ? 4 : 8;
        int arrayBytes = 32/8;

        long sumBytes = 0;

        if (yfit != null) {
            sumBytes += (arrayBytes + (yfit.length*arrayBytes));
        }

        if (x != null) {
            sumBytes += (arrayBytes + (x.length*arrayBytes));
        }

        // 17 variables on the stack, each of size stack word size
        sumBytes += (17 * intBytes);

        // String size on the heap = reference size + content size?
        // parameterNames
        sumBytes += (arrayBytes + (3*(nbits/8) + (intBytes*1 + intBytes*5 + intBytes*2)));

        sumBytes += overheadBytes;

        // amount of padding needed to make it a round 8 bytes
        long padding = (sumBytes % 8);

        sumBytes += padding;

        return sumBytes;
    }

    /**
     *
     * @return
     */
    public String[] getParameterNames() {
        return parameterNames;
    }

    /**
     *
     * @return
     */
    public float[] getParameters() {
        return new float[]{k, sigma, mu};
    }

    public String toString() {

        StringBuffer sb = new StringBuffer();
        sb.append(yfit.length).append(" points, k=").append(k).append(" sigma=").append(sigma)
            .append(" mu=").append(mu).append(" chiSqSum=").append(chiSqSum)
            .append(" chiSqStatistic=").append(chiSqStatistic);

        return sb.toString();
    }

    /**
     *
     * @return
     */
    public float getXPeak() {
        if (xPeakIndex == -1) {
            xPeakIndex = MiscMath.findYMaxIndex(yfit);
        }
        return xScale*x[xPeakIndex];
    }
    
    /**
     *
     * @return
     */
    public int getXPeakIndex() {
        if (xPeakIndex == -1) {
            xPeakIndex = MiscMath.findYMaxIndex(yfit);
        }
        return xPeakIndex;
    }
    
    /**
     *
     * @return
     */
    public float[] getYFit() {
        return yfit;
    }

    /**
     *
     * @return
     */
    public float getK() {
        return k;
    }

    /**
     *
     * @return
     */
    public float getSigma() {
        return sigma;
    }

    /**
     *
     * @return
     */
    public float getMu() {
        return mu;
    }

    /**
     *
     * @return
     */
    public float getKResolution() {
        return kSolutionResolution;
    }

    /**
     *
     * @return
     */
    public float getSigmaResolution() {
        return sigmaSolutionResolution;
    }
    
    /**
     *
     * @return
     */
    public float getMuSolutionResolution() {
        return muSolutionResolution;
    }

    /**
     *
     * @return
     */
    public float getChiSqSum() {
        return chiSqSum;
    }

    /**
     *
     * @return
     */
    public float getYDataErrSq() {
        return yDataErrSq;
    }

    /**
     *
     * @return
     */
    public float getChiSqStatistic() {
        // minus one if mean was computed from the data
        if (yfit == null) {
            return chiSqStatistic;
        }
        return chiSqSum / (yfit.length - 3 - 1);
    }

    /**
     *
     * @param index
     * @param isStepFunction
     * @return
     */
    protected float calculateArea(int index, boolean isStepFunction) {

        return CurveMisc.calculateArea(x, yfit, index, isStepFunction, xScale, yScale);
    }

    /**
     *
     * @param index
     * @return
     */
    public float getX(int index) {
        return xScale*x[index];
    }

    /**
     * @param yfit the yfit to set
     */
    public void setYFit(float[] yfit) {
        this.yfit = yfit;
    }

    /**
     *
     * @param scale
     */
    public void setYScale(float scale) {
        this.yScale = scale;
    }

    /**
     * @return the x array of the fit
     */
    public float[] getX() {
        return x;
    }

    /**
     * @param x array of the fit
     */
    public void setX(float[] x) {
        this.x = x;
    }

    /**
     *
     * @param scale
     */
    public void setXScale(float scale) {
        this.xScale = scale;
    }

    /**
     * @param k the k to set
     */
    public void setK(float k) {
        this.k = k;
    }

    /**
     * @param sigma the sigma to set
     */
    public void setSigma(float sigma) {
        this.sigma = sigma;
    }

    /**
     * @param mu the mu to set
     */
    public void setMu(float mu) {
        this.mu = mu;
    }

    /**
     * @param chiSq the chiSqSum to set
     */
    public void setChiSqSum(float chiSq) {
        this.chiSqSum = chiSq;
    }

    /**
     * @param chiSqStat the chiSqStatistic to set
     */
    public void setChiSqStatistic(float chiSqStat) {
        this.chiSqStatistic = chiSqStat;
    }

    /**
     * @param yErrSq the yDataErrSq to set
     */
    public void setYDataErrSq(float yErrSq) {
        this.yDataErrSq = yErrSq;
    }

    /**
     *
     * @return
     */
    public float[] getOriginalScaleX() {
        if (x == null) {
            return null;
        }
        float[] xsc = new float[x.length];
        for (int i = 0; i < xsc.length; i++) {
            xsc[i] = x[i] * xScale;
        }
        return xsc;
    }

    /**
     *
     * @return
     */
    public float[] getOriginalScaleYFit() {
        if (yfit == null) {
            return null;
        }
        float[] ysc = new float[yfit.length];
        for (int i = 0; i < ysc.length; i++) {
            ysc[i] = yfit[i] * yScale;
        }
        return ysc;
    }

    /**
     *
     * @return
     */
    public float getXScale() {
        return xScale;
    }

    /**
     *
     * @return
     */
    public float getYScale() {
        return yScale;
    }

}
