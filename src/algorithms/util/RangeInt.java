package algorithms.util;

/**
 *
 * @author nichole
 */
public class RangeInt {
    
    protected int start;
    protected int stop;
    
    public RangeInt(final int theStart, final int theStop) {
        start = theStart;
        stop = theStop;
    }
    
    public RangeInt(final RangeInt r) {
        this.start = r.getStart();
        this.stop = r.getStop();
    }
    
    public void setStart(final int theStart) {
        start = theStart;
    }
    
    public void setStop(final int theStop) {
        stop = theStop;
    }
    
    public int getStart() {
        return start;
    }
    
    public int getStop() {
        return stop;
    }
    /**
     * check that start >= minValue and if not, reset it to minValue, and check
     * if stop is less than or equal to maxValue and if not, reset it.
     * @param minValue
     * @param maxValue 
     */
    public void resetBoundsIfNeeded(final int minValue, final int maxValue) {
        if (start < minValue) {
            start = minValue;
        }
        if (start > maxValue) {
            start = maxValue;
        }
        if (stop < minValue) {
            stop = minValue;
        }
        if (stop > maxValue) {
            stop = maxValue;
        }
    }
}
