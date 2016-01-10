package algorithms.imageProcessing;

import algorithms.compGeometry.NearestPoints;
import algorithms.compGeometry.clustering.FixedDistanceGroupFinder;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * given a map of lists of blob matched points, calculate euclidean
 * transformations and evaluate them against all points.
 *
 * @author nichole
 */
public class BlobCornersEuclideanCalculator2 {

    public MatchingSolution solveTransformation(
        GreyscaleImage image1, GreyscaleImage image2,
        IntensityFeatures features1, IntensityFeatures features2, int dither,
        Map<PairInt, List<FeatureComparisonStat>> matchedLists, 
        List<CornerRegion> allCorners1, List<CornerRegion> allCorners2) {
        
        // calculate all transformations
        //List<TransformationParameters> params = calculateTransformations(
        //    matchedLists);
        
        // evaluate against all points or against matched points only
        
        throw new UnsupportedOperationException("not yet implemented");
        
    }
}
