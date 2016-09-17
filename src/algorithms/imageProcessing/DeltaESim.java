package algorithms.imageProcessing;

import algorithms.imageProcessing.util.GroupAverageColors;

/**
 *
 * @author nichole
 */
public class DeltaESim implements Comparable<DeltaESim> {
    
    private int x1;
    private int y1;
    private int x2;
    private int y2;
    private double deltaE;

    public DeltaESim(GroupAverageColors avg1, GroupAverageColors avg2) {
        
        this.x1 = avg1.getXCen();
        this.y1 = avg1.getYCen();
        this.x2 = avg2.getXCen();
        this.y2 = avg2.getYCen();
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        this.deltaE = Math.abs(cieC.calcDeltaECIE2000(
            avg1.getCIEL(), avg1.getCIEA(), avg1.getCIEB(),
            avg2.getCIEL(), avg2.getCIEA(), avg2.getCIEB()));
    }

    @Override
    public int compareTo(DeltaESim other) {
        if (deltaE < other.deltaE) {
            return -1;
        } else if (deltaE > other.deltaE) {
            return 1;
        }
        return 0;
    }

    /**
     * @return the x1
     */
    public int getX1() {
        return x1;
    }

    /**
     * @return the y1
     */
    public int getY1() {
        return y1;
    }

    /**
     * @return the x2
     */
    public int getX2() {
        return x2;
    }

    /**
     * @return the y2
     */
    public int getY2() {
        return y2;
    }

    /**
     * @return the deltaE
     */
    public double getDeltaE() {
        return deltaE;
    }
}
