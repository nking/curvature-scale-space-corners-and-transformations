package algorithms.imageProcessing.features.mser;

import java.util.List;
import algorithms.util.PairIntArray;
import algorithms.util.TrioInt;

/**
 * class to produce descriptors for MSER, usable for
 * matching same objects in other images.
 * 
 * @author nichole
 */
public class Canonicalizer {
   
    // NOTE: the arguments as lists of Regions may change to disjoint
    //    forests or other structure as design of these
    //    methods needed for matching is in progress
    
    public PairIntArray calculateCentroids(List<Region> regions,
        int[] greyscale, int imageWidth, int imageHeight) {
        
        return calculateIntensityCentroids(regions,
            greyscale, imageWidth, imageHeight);
    }
    
    PairIntArray extractRegionXYCenters(List<Region> regions,
        int[] greyscale, int imageWidth, int imageHeight) {
        
        int[] xyCen = new int[2];
        
        PairIntArray output = new PairIntArray(regions.size());
        
        for (int i = 0; i < regions.size(); ++i) {
            Region r = regions.get(i);
            r.calculateXYCentroid(greyscale, imageWidth, imageHeight, 
                xyCen);
     
            output.add(xyCen[0], xyCen[1]);
        }
        
        return output;
    }
    
    PairIntArray calculateIntensityCentroids(List<Region> regions,
        int[] greyscale, int imageWidth, int imageHeight) {
        
        int n = regions.size();
        
        PairIntArray output = new PairIntArray(n);
        
        int[] cenXY = new int[2];
        
        for (int i = 0; i < n; ++i) {
            
            Region region = regions.get(i);
            
            region.calculateIntensityCentroid(greyscale, imageWidth, 
                imageHeight, cenXY);
            
            output.add(cenXY[0], cenXY[1]);
        }
        
        return output;
    }
    
    PairIntArray calculateIntensityCentroids(List<Region> regions,
        int radius, int[] greyscale, int imageWidth, int imageHeight) {
        
        int n = regions.size();
        
        PairIntArray output = new PairIntArray(n);
        
        int[] cenXY = new int[2];
        
        for (int i = 0; i < n; ++i) {
            
            Region r = regions.get(i);
            
            r.calculateIntensityCentroid(greyscale, imageWidth, 
                imageHeight, cenXY, radius);
            
            output.add(cenXY[0], cenXY[1]);
        }
        
        return output;
    }
}
