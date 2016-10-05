package algorithms.imageProcessing;

/**
 algorithm based upon paper
 "Integral Histogram: A Fast Way to Extract Histograms in
   Cartesian Spaces" by Porikli 2005

 Algorithm uses a pattern similar to a summed area table
 to cumulatively add the content from one lower pixel in
 x and y with a runtime complexity that is O(N_pixels),
 * but a space complexity that is O(N_pixels * nBins) where
 * nBins is the number of bins in a histogram (default is 16 for
 * default value range of 0 to 255, inclusive).
 
 * @author nichole
 */
public class IntegralHistograms {
    
}
