package algorithms.imageProcessing;

/**
 * shows horizontal changes in an image
 * 
 * @author nichole
 */
public class SobelX implements IKernel {
     
    public float getNormalizationFactor() {
        return 1.0f;
    }
    
    public Kernel getKernel() {
        /*
                        | -1  0  1 |
           Sobel_x =    | -2  0  2 |
                        | -1  0  1 |
        
           this is the n=2 binomial filter for a Gaussian first derivative,
           that is sigma = sqrt(2)/2 = 0.707 = [1, 0, -1]
        */
        Kernel kernel = new Kernel(3, 3);
        kernel.setValue(0, 0, -1);
        kernel.setValue(0, 1, -2);
        kernel.setValue(0, 2, -1);
        kernel.setValue(2, 0, 1);
        kernel.setValue(2, 1, 2);
        kernel.setValue(2, 2, 1);
        
        return kernel;
    }
}
