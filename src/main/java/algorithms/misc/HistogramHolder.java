package algorithms.misc;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * @author nichole
 */
public class HistogramHolder implements Externalizable {

    private static final long serialVersionUID = -7105371701539621066L;

    protected float[] xHist = null;
    protected int[] yHist = null;
    protected float[] yHistFloat = null;
    protected float[] yErrors = null;
    protected float[] xErrors = null;

    public long approximateMemoryUsed() {

        String arch = System.getProperty("sun.arch.data.model");

        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;

        int nbits = (is32Bit) ? 32 : 64;

        int overheadBytes = 16;

        int intBytes = (is32Bit) ? 4 : 8;
        int arrayBytes = 32/8;

        long sumBytes = 0;

        if (xHist != null) {

            sumBytes += (arrayBytes + (xHist.length*arrayBytes));

            sumBytes += (arrayBytes + (yHist.length*arrayBytes));
        }

        if (yHistFloat != null) {
            sumBytes += (arrayBytes + (yHistFloat.length*arrayBytes));
        }

        if (yErrors != null) {

            sumBytes += (arrayBytes + (xErrors.length*arrayBytes));

            sumBytes += (arrayBytes + (yErrors.length*arrayBytes));
        }

        sumBytes += overheadBytes;

        // amount of padding needed to make it a round 8 bytes
        long padding = (sumBytes % 8);

        sumBytes += padding;

        return sumBytes;
    }

    /**
     * @return the xHist
     */
    public float[] getXHist() {
        return xHist;
    }

    /**
     * @return the yHist
     */
    public int[] getYHist() {
        return yHist;
    }

    /**
     * @return the yHistFloat
     */
    public float[] getYHistFloat() {
        return yHistFloat;
    }

    /**
     * @return the yErrors
     */
    public float[] getYErrors() {
        return yErrors;
    }

    /**
     * @return the xErrors
     */
    public float[] getXErrors() {
        return xErrors;
    }

    /**
     * @param xHist the xHist to set
     */
    public void setXHist(float[] xHist) {
        this.xHist = xHist;
    }

    /**
     * @param yHist the yHist to set
     */
    public void setYHist(int[] yHist) {
        this.yHist = yHist;
    }

    /**
     * @param yHistFloat the yHistFloat to set
     */
    public void setYHistFloat(float[] yHistFloat) {
        this.yHistFloat = yHistFloat;
    }

    /**
     * @param yErrors the yErrors to set
     */
    public void setYErrors(float[] yErrors) {
        this.yErrors = yErrors;
    }

    /**
     * @param xErrors the xErrors to set
     */
    public void setXErrors(float[] xErrors) {
        this.xErrors = xErrors;
    }

    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("histogram=[");
        
        for (int i = 0; i < xHist.length; i++) {
            
            sb.append("(").append(xHist[i]).append(", ");
            
            if (yHist != null) {
                sb.append(yHist[i]);
            } else {
                sb.append(yHistFloat[i]);
            }
            sb.append(") ");
            
        }
        
        sb.append("]\n");
        
        if (xErrors != null) {
            
            sb.append("histogram errors=[");
            
            for (int i = 0; i < xErrors.length; i++) {

                sb.append("(").append(xErrors[i]).append(", ");

                sb.append(yErrors[i]).append(") ");
            }
            
            sb.append("]\n");
        }      
        
        return sb.toString();
    }

    @Override
    public void writeExternal(ObjectOutput out) throws IOException {

        if (xHist == null) {
            return;
        }

        out.writeInt(xHist.length);

        for (int i = 0; i < xHist.length; i++) {

            out.writeFloat(xHist[i]);
            out.writeInt(yHist[i]);
            out.writeFloat(yHistFloat[i]);

            out.flush();
        }

        if (xErrors == null) {
            out.writeInt(0);
        } else {
            out.writeInt(xErrors.length);
        
            for (int i = 0; i < xErrors.length; i++) {

                out.writeFloat(xErrors[i]);
                out.writeFloat(yErrors[i]);

                out.flush();
            }
        }
    }

    @Override
    public void readExternal(ObjectInput in) throws IOException {

        int n = in.readInt();

        xHist = new float[n];
        yHist = new int[n];
        yHistFloat = new float[n];

        for (int i = 0; i < n; i++) {
            xHist[i] = in.readFloat();
            yHist[i] = in.readInt();
            yHistFloat[i] = in.readFloat();
        }

        int ne = in.readInt();

        if (ne == 0) {
            return;
        }

        xErrors = new float[ne];
        yErrors = new float[ne];

        for (int i = 0; i < ne; i++) {
            xErrors[i] = in.readFloat();
            yErrors[i] = in.readFloat();
        }
    }
}
