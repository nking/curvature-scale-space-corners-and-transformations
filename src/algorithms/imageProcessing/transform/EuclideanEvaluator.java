package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.features.FeatureComparisonStat;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author nichole
 */
public class EuclideanEvaluator {

    public EuclideanTransformationFit evaluate(PairIntArray xy1, 
        PairIntArray xy2, TransformationParameters params, int tolerance) {
        
        if (xy1 == null) {
            throw new IllegalArgumentException("xy1 cannot be null");
        }
        if (xy2 == null) {
            throw new IllegalArgumentException("xy2 cannot be null");
        }
        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
        if (xy1.getN() != xy2.getN()) {
            throw new IllegalArgumentException("xy1 and xy2 must be same length");
        }
        
        int n = xy1.getN();
        
        Transformer transformer = new Transformer();
            
        PairIntArray xy1Tr = transformer.applyTransformation(params, xy1);
        
        List<Integer> inlierIndexes = new ArrayList<Integer>();
        List<Double> distances = new ArrayList<Double>();
        
        for (int i = 0; i < n; ++i) {
            int x1 = xy1Tr.getX(i);
            int y1 = xy1Tr.getY(i);
            
            int x2 = xy2.getX(i);
            int y2 = xy2.getY(i);
            
            int diffX = x2 - x1;
            int diffY = y2 - y1;
            
            double dist = Math.sqrt(diffX*diffX + diffY*diffY);
            if (dist < tolerance) {
                inlierIndexes.add(Integer.valueOf(i));
                distances.add(Double.valueOf(dist));
            }
        }
        
        EuclideanTransformationFit fit = new EuclideanTransformationFit(
            params.copy(), inlierIndexes, distances, tolerance);
        
        return fit;
    }
    
    public EuclideanTransformationFit evaluate(List<FeatureComparisonStat> stats, 
        TransformationParameters params, int tolerance) {
        
        if (stats == null) {
            throw new IllegalArgumentException("xy1 cannot be null");
        }
        
        PairIntArray xy1 = new PairIntArray(stats.size());
        PairIntArray xy2 = new PairIntArray(stats.size());
        
        for (FeatureComparisonStat stat : stats) {
            
            PairInt p1 = stat.getImg1Point().copy();
            PairInt p2 = stat.getImg2Point().copy();
            
            xy1.add(p1.getX(), p1.getY());
            xy2.add(p2.getX(), p2.getY());
        }
        
        return evaluate(xy1, xy2, params, tolerance);
    }
   
}
