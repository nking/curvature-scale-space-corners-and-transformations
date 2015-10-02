package algorithms.imageProcessing;

import algorithms.compGeometry.clustering.KMeansPlusPlus;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;

/**
 * class holding several different image segmentation methods.  Note that
 * some other techniques involving contrast for example, are elsewhere.
 * 
 * @author nichole
 */
public class ImageSegmentation {
    
    /**
     * applies KMeansPlusPlus algorithm to the values in input 
     * (greyscale intensities) to create kBands of clustered pixels 
     * (operates on input).
     * @param input
     * @param kBands
     * @throws IOException
     * @throws NoSuchAlgorithmException
     */
    public void applyImageSegmentation(GreyscaleImage input, int kBands)
        throws IOException, NoSuchAlgorithmException {

        KMeansPlusPlus instance = new KMeansPlusPlus();
        instance.computeMeans(kBands, input);

        int[] binCenters = instance.getCenters();

        for (int col = 0; col < input.getWidth(); col++) {

            for (int row = 0; row < input.getHeight(); row++) {

                int v = input.getValue(col, row);

                for (int i = 0; i < binCenters.length; i++) {

                    int vc = binCenters[i];

                    int bisectorBelow = ((i - 1) > -1) ?
                        ((binCenters[i - 1] + vc) / 2) : 0;

                    int bisectorAbove = ((i + 1) > (binCenters.length - 1)) ?
                        255 : ((binCenters[i + 1] + vc) / 2);

                    if ((v >= bisectorBelow) && (v <= bisectorAbove)) {

                        input.setValue(col, row, vc);

                        break;
                    }
                }
            }
        }
    }

}
