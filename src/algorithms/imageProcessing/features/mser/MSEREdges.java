package algorithms.imageProcessing.features.mser;

/**
 * class to explore the boundaries of the accumulated points
 * in MSER regions.
 * 
 * The contents may in the future be replaced by an implementation
 * of 
 * http://www.vision.ee.ethz.ch/~rhayko/paper/aapr2009_boundary_detection_SBER_hayko.pdf
 * 
 * but for now is an exploration of the existing MSER and Region
 * class results.
 * 
 * @author nichole
 */
public class MSEREdges {
    
    /*
    given a color image
      -- makes greyscale image and a polar theta of cie luv
      -- creates positive and negative mser regions for both,
            filtering by variation for the specific mser type
      -- applies a spatial filter to keep one for an overlapping
          or close centered mser.
          -- keeps the mser w/ most accumulated points
             (this is where the component tree approach would be better).
      -- using PerimeterFinder2 extract the outer boundaries from each mser
         region points, and do the same for the contiguous embedded regions.
      
    */
}
