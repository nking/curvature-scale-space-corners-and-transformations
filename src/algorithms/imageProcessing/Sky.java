package algorithms.imageProcessing;

import algorithms.imageProcessing.features.mser.MSEREdges;
import algorithms.util.PairInt;
import java.util.Set;

/**
 * a class with methods related to sky specific tasks such as finding the sun,
 * rainbow, sky, and skyline.
 * 
 * NOT READY FOR USE.
 * still in design phase
 * 
 * @author nichole
 */
public class Sky {
    
    private final MSEREdges mserEdges;
    
    private final ImageExt img;
    
    private boolean debug = false;
    
    public Sky(ImageExt img) {
        
        this.img = img.copyToImageExt();
        
        mserEdges = new MSEREdges(this.img);
        mserEdges.setToDebug();
        mserEdges.mergeAndExtractEdges();
    }
    
    public SkyObject findSun() {
        SunFinder finder = new SunFinder();
        return finder.findSun(mserEdges);
    }
    
    public SkyObject findRainbows() {
        throw new UnsupportedOperationException("not ready for use");
    }
    
    public SkyObject findMoonDogs() {
        throw new UnsupportedOperationException("not ready for use");
    }
    
    public void setToDebug() {
        debug = true;
    }
    
    public void _printGsRegions0() {
        mserEdges._debugOrigRegions(0, "gs");
    }
    public void _printGsRegions1() {
        mserEdges._debugOrigRegions(1, "gs");
    }
    public void _printPtRegions0() {
        mserEdges._debugOrigRegions(2, "pt");
    }
    public void _printPtRegions1() {
        mserEdges._debugOrigRegions(3, "pt");
    }
    
    public static class SkyObject {
        Set<PairInt> points;
        int[] xyCenter;
    }
}
