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

    public int calculateHalfYMaxIndexPastYMax() {

        if (yHistFloat == null) {
            return -1;
        }

        int yMaxIndex = MiscMath.findYMaxIndex(yHistFloat);

        int halfMaxIndex = -1;
        float halfMax = yHistFloat[yMaxIndex]/2.0f;

        for (int i = yMaxIndex; i < yHistFloat.length; i++) {
            if (halfMax <= yHistFloat[i]) {
                halfMaxIndex = i;
            }
        }
        return halfMaxIndex;
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
        }
        for (int i = 0; i < xErrors.length; i++) {

            out.writeFloat(xErrors[i]);
            out.writeFloat(yErrors[i]);

            out.flush();
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
