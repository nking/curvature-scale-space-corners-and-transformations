package algorithms.imageProcessing.features;

import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.SIGMA;
import algorithms.util.CornerArray;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author nichole
 */
public class FeatureHelper {
    
    public static List<CornerArray> filterByLocalizability(
        GreyscaleImage img, GreyscaleImage gradientImage, GreyscaleImage theta,
        List<CornerArray> cornersList) {
        
        RotatedOffsets rotatedOffsets = RotatedOffsets.getInstance();
        
        IntensityFeatures features = new IntensityFeatures(5, true, 
            rotatedOffsets);
        
        features.gXY = gradientImage;
        features.theta = theta;
        
        List<CornerArray> output = new ArrayList<CornerArray>();
        for (CornerArray corners : cornersList) {
         
            if (corners.getN() == 0) {
                output.add(corners);
                continue;
            }
            
            SIGMA sigma = corners.getSIGMA();
            CornerArray corners2 = new CornerArray(sigma);
            if (corners.isFromAClosedCurve()) {
                corners2.setIsClosedCurve();
            }
            
            for (int i = 0; i < corners.getN(); ++i) {
                
                int x = corners.getX(i);
                
                int y = corners.getY(i);
                
                if (!features.removeDueToLocalization(img, x, y,
                    features.calculateOrientation(x, y))) {
                    
                    corners2.add(
                        corners.getX(i), corners.getY(i),
                        corners.getCurvature(i),
                        corners.getXFirstDeriv(i),
                        corners.getXSecondDeriv(i),
                        corners.getYFirstDeriv(i),
                        corners.getYSecondDeriv(i), corners.getInt(i));
                }
            }
            
            output.add(corners2);
        }
        
        return output;
    }
}
