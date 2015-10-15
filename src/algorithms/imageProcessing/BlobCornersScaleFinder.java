package algorithms.imageProcessing;

import algorithms.imageProcessing.IBlobScaleFinder;
import algorithms.imageProcessing.ISegmentedImageHelper;
import algorithms.imageProcessing.SegmentationType;
import algorithms.imageProcessing.TransformationParameters;

/**
 * class to invoke methods needed to solve for euclidean scale between
 * image1 and image2 using methods specific to corners on closed curves.
 * 
 * @author nichole
 */
public class BlobCornersScaleFinder implements IBlobScaleFinder {

    @Override
    public void setToDebug() {
        throw new UnsupportedOperationException("Not supported yet."); 
    }

    @Override
    public TransformationParameters solveForScale(ISegmentedImageHelper sih, 
        SegmentationType st, boolean bln, 
        ISegmentedImageHelper sih1, SegmentationType st1, boolean bln1, float[] floats) {
        throw new UnsupportedOperationException("Not supported yet."); 
    }
    
}
