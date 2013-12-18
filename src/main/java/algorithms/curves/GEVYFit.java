package algorithms.curves;

public class GEVYFit implements IYFit {

    protected float[] yfit;
    protected float[] x;
    protected float xScale = 1;
    protected float yScale = 1;
    protected float k;
    protected float sigma;
    protected float mu;
    protected float chiSqSum = Float.MAX_VALUE;
    protected float chiSqStatistic = Float.MAX_VALUE;
    float kSolutionResolution;
    float sigmaSolutionResolution;
    protected float yDataErrSq;

    // these are only set internally.  values of -1 mean the stats haven't been calculated
    protected int xMeanIndex = -1;
    protected int xPeakIndex = -1;
    protected int xMedianIndex = -1;
    protected int x05PercentIndex = -1;
    protected int x10PercentIndex = -1;
    protected int x80PercentIndex = -1;
    protected int x95PercentIndex = -1;

    protected String[] parameterNames = new String[]{
        "k", "sigma", "mu"
    };

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

    public String[] getParameterNames() {
        return parameterNames;
    }

    public float[] getParameters() {
        return new float[]{k, sigma, mu};
    }

    public String toString() {

        StringBuffer sb = new StringBuffer();
        sb.append(yfit.length).append(" points, k=").append(k).append(" sigma=").append(sigma)
            .append(" mu=").append(mu).append(" chiSqSum=").append(chiSqSum)
            .append(" chiSqStatistic=").append(chiSqStatistic);

        if (statsHaveBeenCalculated()) {
            sb.append(" xMean=").append(getXMean())
                .append(" xPeak=").append(getXPeak())
                .append(" xMedian=").append(getXMedian())
                .append(" x05Percent=").append(getX05Percent())
                .append(" x10Percent=").append(getX10Percent())
                .append(" x80Percent=").append(getX80Percent())
                .append(" x95Percent=").append(getX95Percent());
        }

        return sb.toString();
    }

    public float[] getYFit() {
        return yfit;
    }

    public float getK() {
        return k;
    }

    public float getSigma() {
        return sigma;
    }

    public float getMu() {
        return mu;
    }

    public float getKResolution() {
        return kSolutionResolution;
    }

    public float getSigmaResolution() {
        return sigmaSolutionResolution;
    }

    public float getChiSqSum() {
        return chiSqSum;
    }

    public float getYDataErrSq() {
        return yDataErrSq;
    }

    public float getChiSqStatistic() {
        // minus one if mean was computed from the data
        if (yfit == null) {
            return chiSqStatistic;
        }
        return chiSqSum / (yfit.length - 3 - 1);
    }

    protected float calculateArea(int index, boolean isStepFunction) {

        return CurveMisc.calculateArea(x, yfit, index, isStepFunction, xScale, yScale);
    }

    /**
     * compute the stats using area under the curve and assumption that the
     * curve is a step function (histogram)
     * @return
     */
    public boolean calculateStatsAsStepFunction() {
        return calcStats(true);
    }

    public boolean calculateStatsAsCurve() {
        return calcStats(false);
    }

    public int findIndexForPeakFraction(float peakFraction, boolean isStepFunction) {
        if (yfit == null) {
            return -1;
        }
        if (!statsHaveBeenCalculated()) {
            calcStats(isStepFunction);;
        }

        float findY = peakFraction * yScale * yfit[xPeakIndex];
        int index = -1;
        for (int i = (xPeakIndex + 1); i < yfit.length; i++) {
            float y = yScale * yfit[i];
            if ( y > findY) {
                index = i;
            }
        }
        return index;
    }

    protected boolean calcStats(boolean isStepFunction) {

        float ymax = Float.MIN_VALUE;

        float total = 0;

        float[] binAreas = new float[yfit.length];

        for (int i = 0; i < yfit.length; i++) {

            float binArea = calculateArea(i, isStepFunction);

            binAreas[i] = binArea;

            total += binArea;

            if (binArea > ymax) {
                ymax = binArea;
                xPeakIndex = i;
            }
        }

        if (Float.isNaN(total) || total == 0) {
            return false;
        }


        float limit05 = 0.05f * total;
        float limit10 = 0.10f * total;
        float limit80 = 0.80f * total;
        float limit95 = 0.95f * total;

        float mean = total/yfit.length;

        float total2 = 0;
        for (int i = 0; i < yfit.length; i++) {
            //String str = String.format("[%d] %.7e, %.1f => %.1f %.3f", i, histogram.xHist[i], ygev[i], total2, (total2/total));
            //log.info(str);
            if (total2 < limit05) {
                x05PercentIndex = i;
            }
            if (total2 < limit10) {
                x10PercentIndex = i;
            }
            if (total2 < limit80) {
                x80PercentIndex = i;
            }
            if (total2 < limit95) {
                x95PercentIndex = i;
            }
            if (total2 < mean) {
                xMeanIndex = i;
            }

            total2 += binAreas[i];
        }

        xMedianIndex = (yfit.length/2);

        return true;
    }

    protected boolean statsHaveBeenCalculated() {
        if (yfit == null) {
            return false;
        }
        if (this.xMeanIndex == -1) {
            return false;
        }
        return true;
    }

    protected void checkStateOfCalculationsWithException() {
        if (!statsHaveBeenCalculated()) {
            throw new IllegalStateException("cannot return statistics until calc stats has been run");
        }
    }

    /**
     * @return the xMean
     */
    public float getXMean() {
        checkStateOfCalculationsWithException();
        return xScale*x[xMeanIndex];
    }

    public float getX(int index) {
        checkStateOfCalculationsWithException();
        return xScale*x[index];
    }

    /**
     * @return the xPeak
     */
    public float getXPeak() {
        checkStateOfCalculationsWithException();
        return xScale*x[xPeakIndex];
    }

    /**
     * @return the xMedian
     */
    public float getXMedian() {
        checkStateOfCalculationsWithException();
        return xScale*x[xMedianIndex];
    }

    /**
     * @return the x10Percent
     */
    public float getX05Percent() {
        checkStateOfCalculationsWithException();
        return xScale*x[x05PercentIndex];
    }

    /**
     * @return the x10Percent
     */
    public float getX10Percent() {
        checkStateOfCalculationsWithException();
        return xScale*x[x10PercentIndex];
    }

    /**
     * @return the x80Percent
     */
    public float getX80Percent() {
        checkStateOfCalculationsWithException();
        return xScale*x[x80PercentIndex];
    }

    /**
     * @return the x95Percent
     */
    public float getX95Percent() {
        checkStateOfCalculationsWithException();
        return xScale*x[x95PercentIndex];
    }

    /**
     * @return the xMeanIndex
     */
    public int getxMeanIndex() {
        if (!statsHaveBeenCalculated()) {
            calcStats(true);;
        }
        return xMeanIndex;
    }

    /**
     * @return the xPeakIndex
     */
    public int getXPeakIndex() {
        if (!statsHaveBeenCalculated()) {
            calcStats(true);;
        }
        return xPeakIndex;
    }

    /**
     * @return the xMedianIndex
     */
    public int getXMedianIndex() {
        return xMedianIndex;
    }

    /**
     * @return the x10PercentIndex
     */
    public int getX10PercentIndex() {
        return x10PercentIndex;
    }

    /**
     * @return the x05PercentIndex
     */
    public int getX05PercentIndex() {
        return x05PercentIndex;
    }

    /**
     * @return the x80PercentIndex
     */
    public int getX80PercentIndex() {
        return x80PercentIndex;
    }

    /**
     * @return the x95PercentIndex
     */
    public int getX95PercentIndex() {
        return x95PercentIndex;
    }

    /**
     * @param yfit the yfit to set
     */
    public void setYFit(float[] yfit) {
        this.yfit = yfit;
    }

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

    public float getXScale() {
        return xScale;
    }

    public float getYScale() {
        return yScale;
    }

}
