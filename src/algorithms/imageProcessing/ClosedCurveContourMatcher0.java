package algorithms.imageProcessing;

import java.util.List;
import java.util.ArrayList;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class ClosedCurveContourMatcher0 {
    
    /*
    O(N) pattern which assumes that points in contour1[i] match points
    in contour2[i + deltaIdx] in an ordered manner.

    matching points in the contour lists:
       contour1=contourList1[0]        curve2=contourList2[0]
          tries deltaIdx=0 with contour1[i] : contour2[i] for all i in contour1
          tries deltaIdx=1 with contour1[i] : contour2[i + 1] for all i in contour1
             ...
          tries deltaIdx=n-1 with contour1[i] : contour2[i + (n-1)] for all i in contour1
          and the best becomes the cost for contour1=contourList1[0]    contour2=contourList2[0]
       Note that the matching is not tried if contour1 and contour2 sizes are different
       proceeds with attempt for next in contour list matches,
       contour1=contourList1[0]        curve2=contourList2[1]
       to find the best solution if any for contourList1[0] and stores it.
       Then does the same for contourList1[1]...
       And uses various statistics to find the best solution and
       combine similar with it to make a large point set to solve
       the euclidean transformation with.
    */
    
    protected final List<CurvatureScaleSpaceContour> c1;

    protected final List<CurvatureScaleSpaceContour> c2;
    
    protected final List<BlobPerimeterRegion> cr1;
    
    protected final List<BlobPerimeterRegion> cr2;

    protected final IntensityFeatures features1;

    protected final IntensityFeatures features2;

    protected final int dither = 3;

    protected final int rotationTolerance = 20;

    //TODO: tune this
    private int tolerance = 2;

    private TransformationParameters solutionParameters = null;

    private double solutionCost = Double.MAX_VALUE;

    private List<BlobPerimeterRegion> solutionMatchedContourRegions1 = null;

    private List<BlobPerimeterRegion> solutionMatchedContourRegions2 = null;

    private List<FeatureComparisonStat> solutionMatchedCompStats = null;

    private final Logger log = Logger.getLogger(this.getClass().getName());

    private boolean solverHasFinished = false;

    private boolean hasBeenInitialized = false;

    private boolean solutionHasSomeScalesSmallerThanOne = false;

    public ClosedCurveContourMatcher0(final IntensityFeatures features1,
        final IntensityFeatures features2, 
        final List<CurvatureScaleSpaceContour> contours1,
        final List<CurvatureScaleSpaceContour> contours2,
        final List<BlobPerimeterRegion> regions1,
        final List<BlobPerimeterRegion> regions2) {
        
        if (contours1.size() != contours2.size()) {
            throw new IllegalArgumentException("contours must be same size");
        }

        c1 = new ArrayList<CurvatureScaleSpaceContour>(contours1.size());

        c2 = new ArrayList<CurvatureScaleSpaceContour>(contours2.size());
        
        cr1 = new ArrayList<BlobPerimeterRegion>(regions1.size());
        
        cr2 = new ArrayList<BlobPerimeterRegion>(regions2.size());

        solutionMatchedCompStats = new ArrayList<FeatureComparisonStat>();

        this.features1 = features1;

        this.features2 = features2;

        c1.addAll(contours1);
        c2.addAll(contours2);
        cr1.addAll(regions1);
        cr2.addAll(regions2);

        hasBeenInitialized = true;
    }

    public boolean matchCorners() {

        if (solverHasFinished) {
            throw new IllegalStateException(
            "matchContours cannot be invoked more than once");
        }

        solutionParameters = null;
        solutionCost = Double.MAX_VALUE;
        
        int deltaIdx = 0;
        
        FeatureMatcher featureMatcher = new FeatureMatcher();
        
        List<TransformationPair3> trList = new ArrayList<TransformationPair3>();
        
        while (deltaIdx < c1.size()) {
            
            TransformationPair3 transformationPair 
                = new TransformationPair3(0, deltaIdx, c1.size());
            
            // store matched and cost
            for (int i = 0; i < c1.size(); ++i) {
                
                int j = i + deltaIdx;
                if (j > (c2.size() - 1)) {
                    j = j - c2.size();
                }
                
                BlobPerimeterRegion region1 = cr1.get(i);
                BlobPerimeterRegion region2 = cr2.get(j);
               
                FeatureComparisonStat compStat = 
                    featureMatcher.ditherAndRotateForBestLocation(
                    features1, features2, region1, region2, dither);
                
                if (compStat != null) {
                    transformationPair.addMatched(region1, region2, compStat);
                }                
            }
            
            // need at least 3 point for a cost evaluation
            if (transformationPair.getMatchedContourRegions1().size() > 2) {
                trList.add(transformationPair);
            }
            
            deltaIdx++;
        }
        
        // chose best solution if there are more than 1
        
        /*
        set instance vars with best solution if any
        
        TransformationParameters solutionParameters = null;
        double solutionCost = Double.MAX_VALUE;
        List<BlobPerimeterRegion> solutionMatchedContourRegions1 = null;
        List<BlobPerimeterRegion> solutionMatchedContourRegions2 = null;
        List<FeatureComparisonStat> solutionMatchedCompStats = null;
        */
        
        return true;
    }

    private double distance(double x1, double y1, double x2, double y2) {

        double diffX = x1 - x2;
        double diffY = y1 - y2;

        double dist = Math.sqrt(diffX*diffX + diffY*diffY);

        return dist;
    }

    public TransformationParameters getSolvedParameters() {
        return solutionParameters;
    }

    /**
     * get the curvature scale space images shift between the
     * first set of contours and the second set.
     * @return
     */
    public double getSolvedCost() {
        return solutionCost;
    }

    public boolean scaleIsPossiblyAmbiguous() {
        return solutionHasSomeScalesSmallerThanOne;
    }

    public List<FeatureComparisonStat> getSolutionMatchedCompStats() {
        return solutionMatchedCompStats;
    }

}
