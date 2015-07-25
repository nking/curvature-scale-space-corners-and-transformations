package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;

public class ShapeMatcher {

    public ShapeMatcher() {
    }

    /**
    method to extract general shapes from the images and compare them in order to 
    match points.  It returns a rough Euclidean fit the transformation.
    NOTE that the images may need pre-processing steps before using this.  For example,
    the Brown & Lowe 200? panoramic images of a mountain need to have the sky masked
    out of the image first.
     * @param image1
     * @param image2
     * @param outputMatched1
     * @param outputMatched2
    */
    public TransformationPointFit findMatchingShapes(ImageExt image1, ImageExt image2,
    PairIntArray outputMatched1, PairIntArray outputMatched2) throws 
        IOException, NoSuchAlgorithmException {
        
        GreyscaleImage img1Grey = image1.copyToGreyscale();
        GreyscaleImage img2Grey = image2.copyToGreyscale();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        imageProcessor.applyImageSegmentation(img1Grey, 3);
        imageProcessor.applyImageSegmentation(img2Grey, 3);
        
        throw new UnsupportedOperationException("not yet implemented");
    
        /*
        -- segmentation by cieXY color into 3 bands
        -- dfs contiguous find of each of the 3 groups of sizes > <100?>
           for the largest groups:
              -- convex hull to make shapes.
              -- compare the hulls to the other image:  
                 -- by area and single pixel intensity?  (the cie xy histogram skipping
                    algorithm may assign different colors to same feature, so should
                    probably ignore intensity unless know the images are very similar.
                 -- by shape
                    -- furthest pairs?
                    -- a description of ellipticity?
                    -- points on the hull or centroid matching enough between the 2 images to be used like "corners"?
        */
    }

}
