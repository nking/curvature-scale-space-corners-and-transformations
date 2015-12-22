package algorithms.imageProcessing;

import algorithms.compGeometry.RotatedOffsets;

/**
 *
 * @author nichole
 */
public class IntensityFeatures2 extends IntensityFeatures {

    
    // temporary class to check whether gradient is needed for orientation
    //  if yes, need to refactor around that method
    
    public IntensityFeatures2(GreyscaleImage image, int blockHalfWidths, 
        boolean useNormalizedIntensities, RotatedOffsets rotatedOffsets) {
        super(image, blockHalfWidths, useNormalizedIntensities, rotatedOffsets);
    }
    
    public IntensityFeatures2(final Image image, 
        final int blockHalfWidths, final boolean useNormalizedIntensities,
        RotatedOffsets rotatedOffsets) {
        super(image, blockHalfWidths, useNormalizedIntensities, rotatedOffsets);
    }
    
    public IntensityFeatures2(final int blockHalfWidths, final boolean 
        useNormalizedIntensities, RotatedOffsets rotatedOffsets) {
        super(blockHalfWidths, useNormalizedIntensities, rotatedOffsets);
    }
    
    private GreyscaleImage gXY = null;
    public void calculateGradientWithGreyscale(GreyscaleImage gsImg) {
        
        // temporary check for whether gradient image for orientation is needed
        CannyEdgeFilter filter = new CannyEdgeFilter();
        filter.applyFilter(gsImg.copyImage());
        gXY = filter.getEdgeFilterProducts().getGradientXY();
    }

    @Override
    public int calculate45DegreeOrientation(GreyscaleImage img, 
        int xCenter, int yCenter) throws CornerRegion.CornerRegionDegneracyException {
        return super.calculate45DegreeOrientation(gXY, xCenter, yCenter); 
    }
    
}
