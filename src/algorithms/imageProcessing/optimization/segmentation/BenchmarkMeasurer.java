package algorithms.imageProcessing.optimization.segmentation;

import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Misc;
import algorithms.search.KNearestNeighbors;
import algorithms.util.PairFloat;
import algorithms.util.PairInt;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import thirdparty.HungarianAlgorithm;

/**
 Class to estimate the difference between 
 regions in images segmented by humans and by a segmentation
 algorithm following the practices of the Berkeley
 Segmentation Group.  Details are on their web site
 and published:
 http://www.cs.berkeley.edu/projects/vision/grouping/papers/mfm-pami-boundary.pdf

 The class follows the benchmark map comparison section
 of "Learning to Detect Natural Image Boundaries
 Using Local Brightness, Color, and Texture Cues"
 by Martin, Fowlkes, and Malik.
 
 -- test data and human segmentation data (model) are 
    extracted to get the boundaries of the segmented regions.
 - any human drawn boundaries are valid ground truth
    - output to segmentation code in format
      of "soft boundary map".
        - 1 pixel wide boundary valued from 0 to 1 with 1 being
          greater confidence in the boundary.
        - treated as a min-cost bipartite problem
          w/ edges = distance between segm pixel and human boundary pixel.
          and a distance d_max beyond which distance is a non-match.
          any isolated node, that is node without neighbors, can
          be removed from the graph.
          For any node in the graph, only 6 of the nearest points are
          kept before matching.
          - traditionally, one would make the boundary 0 or 1
            by choosing a threshold, but that has 2 problems:
            - the optimal threshold is dependent upon app use
            - removing low level thresholds removes too much information
          - the "threshold" is calculated by trying many "levels"
            (e.g. 30) and treating those above the threshold as '1's
            - at each level calculate "precision" and "recall"
              (precision and recall being similar, but diff from ROC
              curve axes)
              precision: probability that a machine=generated boundary
                  pixel is a true boundary pixel.
                  it's the fraction of detections which are true positives
                  using a distance tolerance of 2 pixels.
                  (= number of correct positive results/number of all positive results)
              recall: probability that a true boundary pixel is detected.
                  (= number of correct positive results/number of all results)
                  for example, 7 objects total, and 4 are correctly identified
                  is recall 4/7
              - need a single number result.
                can do this where the precision-recall curves do not intersect
                and are roughly parallel.
                the curve furthest from the origin dominate and
                the summary statistic is that distance, called F-measure.
              F-measure: harmonic mean of precision and recall.
                defined at all points on the precision-recall curve.
                the maximum is the summary statistic.
                defined as F = P * R / (alpha * R + (1 - alpha)*P).
                   where P is precision and r is recall. authors use alpha=0.5.
                  The location of the maximum F-measure along
                  the curve provides the optimal threshold given,
                  which they set to 0.5 in our experiments.
              *note, they choose "precision" instead of the false positives
                in recall because "precision" scales linearly w/ resolution
                while false positives scale by square power w/ resolution.
              *note, also that "precision" "recall" curves are also
               useful for other comp vis applications in determining s/n thresholds.
  
 @author nichole
 */
public class BenchmarkMeasurer {
    
    public float evaluate(SegmentationResults data,
        SegmentationResults model, int dMax) {
        
        /*
        for each point in perimeters, find the 6 nearest
        neighbors in expected.perimeters.

        discard those more distant than dMax=2
        
        after have the set of candidates,
            discard those with no candidate neighbors        
        */
        
        int nDataPerimeterPoints = data.sumNPerimeters();
        
        if (nDataPerimeterPoints == 0) {
            return 0;
        }
        
        int dMaxSq = dMax * dMax;
        
        Set<PairInt> allExpectedPoints = model.getAllPoints();

        int k = 6;
        
        // searching the data for expected points
        KNearestNeighbors kNN = data.createKNN();
        
        //key = expected point, value = data point and distance
        Map<PairInt, Map<PairInt, Float>> nearestNeighbors
             = new HashMap<PairInt, Map<PairInt, Float>>();
         
        for (PairInt p : allExpectedPoints) {
                
            int x = p.getX();
            int y = p.getY();
     
            // search data for model (x, y)
            List<PairFloat> nearestMatches = 
                kNN.findNearest(k, x, y, dMax);

            if (nearestMatches == null) {
                continue;
            }
            
            for (PairFloat p2 : nearestMatches) {

                float dist = distance(p2, x, y);

                assert(dist <= dMax);
                    
                // Map<PairInt, Map<PairInt, Float>> nearestNeighbors
                Map<PairInt, Float> map = nearestNeighbors.get(p);

                if (map == null) {
                    map = new HashMap<PairInt, Float>();
                    nearestNeighbors.put(p, map);
                }

                PairInt p3 = new PairInt(Math.round(p2.getX()), 
                    Math.round(p2.getY()));

                map.put(p3, Float.valueOf(dist));
            }
        }
        
        // discard the isolated points
        removeIsolatedPoints(nearestNeighbors); 
         
        /*
        number these:
            Map<PairInt, Map<PairInt, Float>> nearestNeighbors
            Set<PairInt> nnValueKeys = new HashSet<PairInt>();
        
        and put in cost array
        */
        int n1 = nearestNeighbors.size();
        int n2 = countUniqueValues(nearestNeighbors);
        
        if (n1 == 0 || n2 == 0) {
            return 0;
        }
        
        boolean transposed = false;
        float[][] cost;
        if (n1 > n2) {
            transposed = true;
            cost = new float[n2][n1];
        } else {
            cost = new float[n1][n2];
        }
        for (int i = 0; i < cost.length; ++i) {
            cost[i] = new float[cost[0].length];
            Arrays.fill(cost[i], -1);
        }
        
        Map<PairInt, Integer> p1Map = new HashMap<PairInt,Integer>();
        Map<PairInt, Integer> p2Map = new HashMap<PairInt,Integer>();
        n1 = 0;
        n2 = 0;
        for (Entry<PairInt, Map<PairInt, Float>> entry :
            nearestNeighbors.entrySet()) {
            
            PairInt p1 = entry.getKey();
            
            Integer i1 = p1Map.get(p1);
            if (i1 == null) {
                i1 = Integer.valueOf(n1);
                p1Map.put(p1, i1);
                n1++;
            }
            
            for (Entry<PairInt, Float> entry2 : 
                entry.getValue().entrySet()) {
                
                PairInt p2 = entry2.getKey();
                Float c = entry2.getValue();
                
                Integer i2 = p2Map.get(p2);
                if (i2 == null) {
                    i2 = Integer.valueOf(n2);
                    p2Map.put(p2, i2);
                    n2++;
                }
                if (transposed) {
                    cost[i2.intValue()][i1.intValue()] 
                        = c.floatValue();
                } else {
                    cost[i1.intValue()][i2.intValue()] 
                        = c.floatValue();
                }
            }
        }
        
        // find min cost matches
        HungarianAlgorithm ha = new HungarianAlgorithm();
        int[][] match = ha.computeAssignments(cost);

        int nMatched = 0;
        for (int i = 0; i < match.length; i++) {
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }
            PairInt pI = new PairInt(idx1, idx2);
            
            nMatched++;
        }
        /*
        precision: probability that a machine=generated boundary
          pixel is a true boundary pixel.
          it's the fraction of detections which are true positives
          using a distance tolerance of 2 pixels.
          (= number of correct positive results/number of all positive results)
        recall: probability that a true boundary pixel is detected.
          (= number of correct positive results/number of all results)
          for example, 7 objects total, and 4 are correctly identified
          is recall 4/7
        recipe is to make a precision vs recall curve and look at the
          value furthest from the origin.
        F-measure: harmonic mean of precision and recall.
          defined at all points on the precision-recall curve.
          the maximum is the summary statistic.
          defined as F = P * R / (alpha * R + (1 - alpha)*P).
          where P is precision and r is recall. authors use alpha=0.5.
          The location of the maximum F-measure along
          the curve provides the optimal threshold given,
          which they set to 0.5 in their experiments.
        Here, not using several levels of thresholding.
        */
        
        float recall = (float)nMatched/
            (float)nDataPerimeterPoints;
            
        float precision = (float)nMatched/
            (float)allExpectedPoints.size();
        
        float alpha = 0.5f;
        
        float fMeasure = (precision * recall)/ 
            (alpha * recall + (1.f - alpha)*precision);
        
        System.out.println("fMeasure=" + fMeasure);
        
        return fMeasure;
    }
    
    private float distance(PairFloat p2, int x, int y) {

        float diffX = p2.getX() - x;
        float diffY = p2.getY() - y;
        
        double dist = Math.sqrt(diffX * diffX + diffY * diffY);
        
        return (float)dist;
    }

    private void removeIsolatedPoints(
        Map<PairInt, Map<PairInt, Float>> nearestNeighbors) {

        Set<PairInt> rm = new HashSet<PairInt>();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        for (Entry<PairInt, Map<PairInt, Float>> entry :
            nearestNeighbors.entrySet()) {
            
            PairInt p = entry.getKey();
            int n = 0;
            for (int i = 0; i < dxs.length; ++i) {
                PairInt p2 = new PairInt
                    (p.getX() + dxs[i], p.getY() + dys[i]);
                if (nearestNeighbors.containsKey(p2)) {
                    n++;
                    break;
                }
            }
            if (n == 0) {
                rm.add(p);
            }
        }
        
        for (PairInt p : rm) {
            nearestNeighbors.remove(p);
        }
    }

    private int countUniqueValues(
        Map<PairInt, Map<PairInt, Float>> nearestNeighbors) {
        
        Set<PairInt> v = new HashSet<PairInt>();
        
        for (Entry<PairInt, Map<PairInt, Float>> entry :
            nearestNeighbors.entrySet()) {
            
            v.addAll(entry.getValue().keySet());
        }
        
        return v.size();
    }
}
