package algorithms.imageProcessing.features;

/**
 * NOT READY FOR USE.  in design stage...
 *
 * An implementation of Histograms of Oriented Gradients
 * constructed from reading the papers
 * <pre>
 * "Histograms of Oriented Gradients for Human Detection" 
 * by Dalal and Triggs, 2010 
 * and
 * "Distinctive Image Features from Scale-Invariant Keypoints"
 * by Lowe, 2004
 * </pre>
 * 
 * The histograms of oriented gradients are constructed to compare
 * the 2 part data of gradient angle and gradient magnitude
 * between regions within different images as a feature descriptor.
 * The feature is defined over several pixels surrounding the central
 * keypoint in a pattern to improve the measurement.
   
   The IntegralHistogram is used to store histograms of values that
   that are counts defined by gradient magnitudes in bins of
   orientation.
   
   As recommended by Dalal & Triggs, 9 bins are used for 180 degrees of
   gradient angle range.
   
   The building of the integral image has a runtime complexity of
   O(N_pixels) for the gradient and O(N_pixels) for the histogram
   image.
   
   Extraction of histogram data is 4 steps, just as in summed area
   tables, but there is additionally the copy which is nBins steps.
  
   The extraction of data is at the "cell" level, which is recommended
   to be 6 X 6 pixels^2 by Dalal and Triggs.
   
   a block of cells is gathered for a point and that is an addition of
   the N_cells X N_cells histograms.
   
   Before the addition, block level normalization is calculated for each
   cell.
   
   The block level normalization uses L2NormHys.
   They found best results using normalization of each cell histogram
   by the total over the block.
   The total number of values is summed over all cells within the
   block and a normalization factor for each cell is then computed using 
   that block total.  The individually normalized cells are then added
   over the block to create the block histogram.
       for each cell calculate cell_total_count.
       total_block_count = sum over cells ( cell_total_count )
       for each cell, normalization is 1/total_block_count
   then the block histogram is added over the same bins in each cell.
   To keep the block histogram as integer but normalized to same 
   max value could apply a further factor of 
   max possible value for a block being, 
   for example (2X2)*(6X6)*(255) = 36720.

   Note that a shift is needed for identifying the bin that is equivalent
   the angle 0 for example, that is a shift specific to a dominant angle
   correction for the histogram.
   That shift will be applied during the addition stage to produce a
   canonicalized descriptor in reference frame w.r.t. rotation.
   (Note that for the use case here, the reference frame orientation will
   be supplied to the method. it's learned from the mser ellipse in one
   use case for example. so the application of a dominant orientation
   will be the same, but the calculation will not be performed for it.
   though that could be added as a method later...not necessary for
   current requirements).

   Comparison of block descriptors is then a histogram intersection,
   where normally 0 is no intersection, hence maximally different,
   and an intersection equal to the max value is maximally similar.
   (see the method ColorHistogram.intersection, but here, the normalization
   will already have been applied instead of determined in the method).
  
  Other details are in converting the intersection to a cost or score
  and specialized methods for that specific to this project will be
  present in this class.
   
  
 * @author nichole
 */
public class HOGs {
    
}
