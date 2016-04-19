package algorithms.imageProcessing.features;

/**
 *
 * @author nichole
 */
 public class PhaseCongruencyParameters {

     private int nScale = 5;
     private int minWavelength = 3;
     private float mult = 2.1f;
     private float sigmaOnf = 0.55f;
     private int k = 5;//2;
     private float cutOff = 0.5f;
     private float g = 10;
     private float deviationGain = 1.5f;
     private int noiseMethod = -1;
     private double tLow = 0.0001;//0.1;
     private double tHigh = 0.1;//0.3;
     private boolean increaseKIfNeeded = true;

     public void setParameters(int nScale, int minWavelength, float mult,
         float sigmaOnf, int k, float cutOff, float g, float deviationGain,
         int noiseMethod, double tLow, double tHigh, boolean increaseKIfNeeded) {
         
         this.nScale = nScale;
         this.minWavelength = minWavelength;
         this.mult = mult;
         this.sigmaOnf = sigmaOnf;
         this.k = k;
         this.cutOff = cutOff;
         this.g = g;
         this.deviationGain = deviationGain;
         this.noiseMethod = noiseMethod;
         this.tLow = tLow;
         this.tHigh = tHigh;
         this.increaseKIfNeeded = increaseKIfNeeded;
     }

    /**
     * @return the nScale
     */
    public int getNScale() {
        return nScale;
    }

    /**
     * @return the minWavelength
     */
    public int getMinWavelength() {
        return minWavelength;
    }

    /**
     * @return the mult
     */
    public float getMult() {
        return mult;
    }

    /**
     * @return the sigmaOnf
     */
    public float getSigmaOnf() {
        return sigmaOnf;
    }

    /**
     * @return the k
     */
    public int getK() {
        return k;
    }

    /**
     * @return the cutOff
     */
    public float getCutOff() {
        return cutOff;
    }

    /**
     * @return the g
     */
    public float getG() {
        return g;
    }

    /**
     * @return the deviationGain
     */
    public float getDeviationGain() {
        return deviationGain;
    }

    /**
     * @return the noiseMethod
     */
    public int getNoiseMethod() {
        return noiseMethod;
    }

    /**
     * @return the tLow
     */
    public double gettLow() {
        return tLow;
    }

    /**
     * @return the tHigh
     */
    public double gettHigh() {
        return tHigh;
    }

    /**
     * @return the increaseKIfNeeded
     */
    public boolean doIncreaseKIfNeeded() {
        return increaseKIfNeeded;
    }
    
    @Override
    public String toString() {
         
         StringBuilder sb = new StringBuilder();
         
         sb.append("nScale=").append(nScale)
             .append(" minWavelength=").append(minWavelength)
             .append(" mult=").append(mult)
             .append(" sigmaOnf=").append(sigmaOnf)
             .append(" k=").append(k)
             .append(" cutOff=").append(cutOff)
             .append(" g=").append(g)
             .append(" deviationGain=").append(deviationGain)
             .append(" noiseMethod=").append(noiseMethod)
             .append(" tLow=").append(tLow)
             .append(" tHigh=").append(tHigh)
             .append(" increaseKIfNeeded=").append(increaseKIfNeeded);
         
         return sb.toString();
     }
}