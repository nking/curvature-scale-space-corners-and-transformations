package algorithms.imageProcessing;

import algorithms.imageProcessing.features.mser.MSEREdges;
import algorithms.util.PairInt;
import java.util.List;
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
    
    private RainbowFinder rFinder = null;
    
    private SunFinder sFinder = null;
    
    public Sky(ImageExt img) {
        
        this.img = img.copyToImageExt();
        
        mserEdges = new MSEREdges(this.img);
        mserEdges.setToDebug();
        mserEdges.setToLowerContrast();
        mserEdges.mergeAndExtractEdges();
    }
    
    public SkyObject findSun() {
        if (sFinder == null) {
            sFinder = new SunFinder();
        }
        return sFinder.findSun(mserEdges);
    }
    
    public List<SkyObject> findRainbows() {
        if (rFinder == null) {
            rFinder = new RainbowFinder();
        }
        return rFinder.findRainbows(mserEdges);
    }
    
    public SkyObject findMoonDogs() {
        throw new UnsupportedOperationException("not ready for use");
    }
    
    public GreyscaleImage extractSkyMask() {
        
        // (1) choosing labeled regions as candidates that
        //  have constancy in color but may have a gradient in illumination.
        //  blue - can use the polar theta images or mser edges (ptRegions[1]) 
        //  then the lch, c image, which is the magnitude of the
        //   LUV u and v, can find the gradual change if any
        //  dark red - same except ptRegions[0]
        
       
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
