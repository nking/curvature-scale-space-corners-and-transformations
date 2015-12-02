package algorithms.imageProcessing;

import algorithms.misc.MedianSmooth;
import java.util.List;

/**
 *
 * @author nichole
 */
public class MedianTransform {

    /**
     * A computationally expensive multiscale median transform.
     * see, the pyramidal mean transform version instead.
     * This method has a runtime complexity of n_iter * O(N_pixels * lg2(windowArea))
     * where windowArea grows from 1 to 2*2*lg2(imageDimension) + 1 
     * and nIter = lg2(imageDimension).
     * @param input
     * @param outputTransformed
     * @param outputCoeff
     */
    public void multiscaleMedianTransform(GreyscaleImage input,
        List<GreyscaleImage> outputTransformed, List<GreyscaleImage> outputCoeff) {

        int imgDimen = Math.min(input.getWidth(), input.getHeight());

        GreyscaleImage img0 = input.copyImage();

        int nr = (int)(Math.ceil(Math.log(imgDimen)/Math.log(2)));
        int s = 1;

        outputTransformed.add(img0.copyToSignedImage());
        
        outputCoeff.add(img0.createSignedWithDimensions());

        for (int j = 0; j < (nr - 1); ++j) {

            int winL = 2*s + 1;
            
            MedianSmooth med = new MedianSmooth();
            
            outputTransformed.add(med.calculate(outputTransformed.get(j), winL, winL));
            
            outputCoeff.add(outputTransformed.get(j).subtract(outputTransformed.get(j + 1)));

            s = 2*s;
        }
        
    }

    public GreyscaleImage reconstructMultiscaleMedianTransform(GreyscaleImage
        c0, List<GreyscaleImage> mmCoeff) {

        int nr = mmCoeff.size();

        GreyscaleImage output = c0.copyToSignedImage();

        for (int j = 0; j < nr; ++j) {
            output = output.add(mmCoeff.get(j));
        }

        return output;
    }

}
