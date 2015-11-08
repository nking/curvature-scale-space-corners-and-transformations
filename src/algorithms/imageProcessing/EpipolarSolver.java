package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * class encapsulating the steps from scale calculation to matching corners
 * to making correspondence lists to solving for epipolar projection.
 * 
 * @author nichole
 */
public class EpipolarSolver {
    
    private final ImageExt img1;
    private final ImageExt img2;
    
    private GreyscaleImage gsImg1 = null;
    private GreyscaleImage gsImg2 = null;
        
    private Set<CornerRegion> cornerRegions1 = null;
    private Set<CornerRegion> cornerRegions2 = null;
    
    private enum State {
        DID_APPLY_HIST_EQ, COULD_NOT_DETERMINE_SCALE
    }
    private Set<State> stateSet = new HashSet<State>();
    
    private final boolean doDetermineScale;
    
    private final boolean debug;
    
    private final String debugTagPrefix;
    
    private TransformationParameters params = null;
            
    private float scaleTol = 0.2f;
    
    private float rotationInRadiansTol = (float)(20. * Math.PI/180.);
    
    //TODO: revise this...
    private int transXYTol = 20;
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public EpipolarSolver(ImageExt image1, ImageExt image2) {
        img1 = image1;
        img2 = image2;
        doDetermineScale = true;
        debug = false;
        debugTagPrefix = "";
    }
    
    public EpipolarSolver(ImageExt image1, ImageExt image2, 
        String debugTagPrefix) {
        img1 = image1;
        img2 = image2;
        doDetermineScale = true;
        debug = true;
        this.debugTagPrefix = debugTagPrefix;
    }
    
    /**
     * constructor accepting transformation parameters.  Note, for best results,
     * the standard deviations within parameters should be populated because they
     * are used as tolerances in matching.
     * @param image1
     * @param image2
     * @param parameters 
     */
    public EpipolarSolver(ImageExt image1, ImageExt image2, 
        TransformationParameters parameters) {
        img1 = image1;
        img2 = image2;
        doDetermineScale = false;
        params = parameters;
        debug = false;
        debugTagPrefix = "";
    }
    
    /**
     * constructor accepting transformation parameters and a debugging tag for
     * image names.  Note, for best results, the standard deviations within 
     * parameters should be populated because they are used as tolerances in 
     * matching.
     * @param image1
     * @param image2
     * @param parameters
     * @param debugTagPrefix 
     */
    public EpipolarSolver(ImageExt image1, ImageExt image2, 
        TransformationParameters parameters, String debugTagPrefix) {
        img1 = image1;
        img2 = image2;
        doDetermineScale = false;
        params = parameters;
        debug = true;
        this.debugTagPrefix = debugTagPrefix;
    }
    
}
