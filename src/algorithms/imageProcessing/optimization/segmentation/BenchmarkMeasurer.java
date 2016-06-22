package algorithms.imageProcessing.optimization.segmentation;

/**
 Class to estimate the difference between 
 regions in images segmented by humans and by a segmentation
 algorithm.
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
          - traditionally, one would make the boundary 9 or 1
            by choosing a threshhold, but that has 2 problems:
            - the optimal threshhold is dependent upon app use
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
    
}
