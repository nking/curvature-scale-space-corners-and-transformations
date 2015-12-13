package algorithms.imageProcessing;

import algorithms.compGeometry.NearestPoints;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import thirdparty.HungarianAlgorithm;

/**
 *  NOT READY FOR USE
 * @author nichole
 */
public class BlobContoursScaleFinder extends AbstractBlobScaleFinder {
    
    /*
    this needs to be revised to use corner contours only.
    need a matcher for sillhouetted objects, for example.
    
    possibly based upon CurvatureScaleSpaceInflectionMapper
    
    */

    public MatchingSolution solveForScale(
        BlobContourHelper img1Helper, IntensityFeatures features1,
        SegmentationType type1, boolean useBinned1,
        BlobContourHelper img2Helper, IntensityFeatures features2,
        SegmentationType type2, boolean useBinned2) {

        throw new UnsupportedOperationException("not implemented");
    }
}
