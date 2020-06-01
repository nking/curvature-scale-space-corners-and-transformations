package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.features.FeatureComparisonStat;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class EuclideanEvaluator {

    public EuclideanTransformationFit evaluate(PairIntArray xy1, 
        PairIntArray xy2, TransformationParameters params, double tolerance) {
        
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
        
        fit.calculateErrorStatistics();
        
        return fit;
    }
    
    public EuclideanTransformationFit evaluate(List<FeatureComparisonStat> stats, 
        TransformationParameters params, double tolerance) {
        
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
   
    /**
     * transform the set templateSetToTransform by given parameters, then
     * calculate the F1 score using precision and recall.
     * 
     * @param templateSetToTransform
     * @param set2
     * @param params
     * @param tolerance
     * @return 
     */
    public float transformAndCalculateF1Score(Set<PairInt> templateSetToTransform, 
        Set<PairInt> set2, TransformationParameters params,
        double tolerance) {
        
        Transformer transformer = new Transformer();
        
        Set<PairInt> templateSet = transformer.applyTransformation2(
            params, templateSetToTransform);
        
        return calculateF1Score(templateSet, set2, tolerance);
    }
    
    public float calculateF1Score(Set<PairInt> templateSet, Set<PairInt> set2, 
        double tolerance) {
        
        /* matching the aggregated adaptive means points to
           the expected template points which have been transformed to the 
           same reference frame.
        
        for metrics, the scores are 0 to 1 where 1 is best possible.
            accuracy = (T_p + T_n)/(T_p + T_n + F_p + F_n)
            precision = (T_p)/(T_p + F_p)
            recall = (T_p)/(T_p + F_n)
            F_1 = 2.* precision * recall/(precision + recall)

        where T_p = expected matches and found them 
              F_p = expected matches, but did not find them
              F_n = expected no matches, but did find matches
              T_n = expected no matches, and did not find them
       
        T_p, F_p : loop over trTemplateInner to find closest match within tolerance
                   in aggInner.
                   -- needs NearestNeighbor2D for aggInner points
        F_n:     : all the points remaining in aggInner that were not matched
        */
        
        int[] minMaxXY = MiscMath.findMinMaxXY(set2);
        int[] minMaxXY2 = MiscMath.findMinMaxXY(templateSet);
        int maxX = Math.max(minMaxXY[1], minMaxXY2[1]);
        int maxY = Math.max(minMaxXY[3], minMaxXY2[3]);
        
        //TODO: find the intersection of points who are each other's
        //     best nearest neighbor matched when transformed and
        //     reverse transformed.
        //     currently not evaluating for intersection with reverse best...
        
        NearestNeighbor2D nn = new NearestNeighbor2D(set2, maxX, maxY);
        
        Set<PairInt> matched = new HashSet<PairInt>();
        
        int d = (int)Math.ceil(tolerance);
        
        int tPos = 0;
        int fPos = 0;
        int fNeg = 0;
        for (PairInt trP : templateSet) {
            Set<PairInt> closest = nn.findClosest(
                trP.getX(), trP.getY(), d);
            if (closest == null || closest.isEmpty()) {
                fPos++;
                continue;
            }
            boolean found = false;
            for (PairInt p : closest) {
                if (!matched.contains(p)) {
                    found = true;
                    matched.add(p);
                    tPos++;
                    break;
                }
            }
            if (!found) {
                fPos++;
            }
        }
        
        if (set2.size() > tPos) {
            fNeg = set2.size() - tPos;
        }
     
        return fMeasure(tPos, fPos, fNeg, 1.0f);
    }
    
    protected float fMeasure(int tPos, int fPos, int fNeg, float beta) {
        
        float betaSq = beta * beta;
        
        float precision = (float)tPos/(float)(tPos + fPos);
        float recall = (float)tPos/(float)(tPos + fNeg);
        float f = (1.f + betaSq) * precision * recall/((betaSq * precision) + recall);

        return f;
    }
       
}
