package algorithms.imageProcessing.features;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.imageProcessing.transform.EpipolarTransformer;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.ejml.simple.SimpleMatrix;

/**
 * class that uses color features and a lower resolution segmentation to
 * create groups of points and match those separately for their best.
 * The matching pattern should be better for objects which have moved with
 * respect to other objects in the image compared to the EpipolarSolver,
 * for example, which solves for a single projection of all matched points.
 * 
 * @author nichole
 */
public class EpipolarColorSegmentedSolver {

    private final ImageExt img1;
    private final ImageExt img2;

    private enum State {
        INITIALIZED, SOLVED, NO_SOLUTION
    }
    private State state = null;

    private List<List<PairInt>> solutionMatched1 = null;
    private List<List<PairInt>> solutionMatched2 = null;
    private List<List<FeatureComparisonStat>> solutionStats = null;
        
    private final FeatureMatcherSettings featureSettings;

    private Logger log = Logger.getLogger(this.getClass().getName());

    public EpipolarColorSegmentedSolver(ImageExt image1, ImageExt image2,
        FeatureMatcherSettings settings) {

        img1 = image1;
        img2 = image2;

        featureSettings = settings.copy();

        state = State.INITIALIZED;
    }

    public boolean solve() throws IOException,
        NoSuchAlgorithmException {

        if (!state.equals(State.INITIALIZED)) {
            throw new IllegalStateException("solve() has already been invoked");
        }

        TMPNonEuclideanSegmentFeatureMatcherColor matcher = 
            new TMPNonEuclideanSegmentFeatureMatcherColor(img1, img2, featureSettings);
        
        boolean solved = matcher.match();

        if (!solved) {
            state = State.NO_SOLUTION;
            return false;
        }
        
        solutionMatched1 = matcher.getSolutionMatched1();
        solutionMatched2 = matcher.getSolutionMatched2();
        solutionStats = matcher.getSolutionStats();
        
        assert(solutionMatched1.size() == solutionMatched2.size());
        assert(solutionStats.size() == solutionMatched1.size());
        
        state = State.SOLVED;

        if (featureSettings.debug()) {
            
            ImageExt img1Cp = img1.copyToImageExt();
            ImageExt img2Cp = img2.copyToImageExt();
            
            for (int i = 0; i < solutionMatched1.size(); ++i) {
                
                int[] clr = ImageIOHelper.getNextRGB(i);
                
                ImageIOHelper.<List<PairInt>>addToImage(solutionMatched1.get(i), 0, 0,
                    img1Cp, 3, clr[0], clr[1], clr[2]);
                
                ImageIOHelper.<List<PairInt>>addToImage(solutionMatched2.get(i), 0, 0,
                    img2Cp, 3, clr[0], clr[1], clr[2]);
            }
            
            MiscDebug.writeImage(img1Cp, featureSettings.getDebugTag() + "_final_1");
            MiscDebug.writeImage(img2Cp, featureSettings.getDebugTag() + "_final_2");
        }
        
        return true;
    }

}
