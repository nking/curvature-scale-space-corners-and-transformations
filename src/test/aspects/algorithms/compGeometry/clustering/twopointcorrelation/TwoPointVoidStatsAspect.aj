package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

public aspect TwoPointVoidStatsAspect {

    protected int x0, x1, y0, y1;
    
    protected boolean testWithinBounds = false;

	before(int xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi, int nDiv, float bFactor) : 
	call(void *.findVoidsRoughRangeSearch( int, int, int, int, int, float ) ) && args(xIndexLo, xIndexHi, yIndexLo, yIndexHi, nDiv, bFactor) 
	&& target(algorithms.compGeometry.clustering.twopointcorrelation.TwoPointVoidStats) { 
						
	    TwoPointVoidStats cls = (TwoPointVoidStats)thisJoinPoint.getThis();
	    
		Logger log = Logger.getLogger(cls.getClass().getName());
		
		this.x0 = xIndexLo;
		
		this.x1 = xIndexHi;
		
		this.y0 = yIndexLo;
		
		this.y1 = yIndexHi;
		
		this.testWithinBounds = true;
		
		log.info("***ASPECT storing index ranges in aspect to assert in methods it invokes");
	}

    after(int xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi, int nDiv, float bFactor) : 
    call(void *.findVoidsRoughRangeSearch( int, int, int, int, int, float ) ) && args(xIndexLo, xIndexHi, yIndexLo, yIndexHi, nDiv, bFactor) 
    && target(algorithms.compGeometry.clustering.twopointcorrelation.TwoPointVoidStats) { 
                
        TwoPointVoidStats cls = (TwoPointVoidStats)thisJoinPoint.getThis();
        
        Logger log = Logger.getLogger(cls.getClass().getName());
      
        this.testWithinBounds = false;
        
        log.info("***ASPECT finished asserts for findVoidsRoughRangeSearch");
    }
    
	before(int xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi) : 
	call(void *.processIndexedRegion(int, int, int, int) ) && args(xIndexLo, xIndexHi, yIndexLo, yIndexHi) 
	&& target(algorithms.compGeometry.clustering.twopointcorrelation.TwoPointVoidStats) { 
		 
	    // cannot use logger here unless as an instance variable
				
	    if (!testWithinBounds) {
	        return;
	    }
		                
        if (!isWithinRange(xIndexLo, x0, x1) ) {
            String err = "ERROR: " + xIndexLo + " is not within bounds of [" + x0 + ":" + x1 + "]";
            System.err.println(err);
            throw new IllegalStateException(err);
        }
        
        if (!isWithinRange(xIndexHi, x0, x1) ) {
            String err = "ERROR: " + xIndexHi + " is not within bounds of [" + x0 + ":" + x1 + "]";
            System.err.println(err);
            throw new IllegalStateException(err);
        }
        
        if (!isWithinRange(yIndexLo, y0, y1) ) {
            String err = "ERROR: " + yIndexLo + " is not within bounds of [" + y0 + ":" + y1 + "]";
            System.err.println(err);
            throw new IllegalStateException(err);
        }
        
        if (!isWithinRange(yIndexHi, y0, y1) ) {
            String err = "ERROR: " + yIndexHi + " is not within bounds of [" + y0 + ":" + y1 + "]";
            System.err.println(err);
            throw new IllegalStateException(err);
        }
	}
	
	public boolean isWithinRange(int testVar, int var0, int var1) {
	    return (  (testVar >= var0) && (testVar <= var1) );
	}
}
