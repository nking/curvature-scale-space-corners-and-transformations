package algorithms.imageProcessing;

/**
 * shows vertical changes in an image
 * 
 * @author nichole
 */
public class SobelY implements IKernel {
    
    public float getNormalizationFactor() {
        return 1.0f;
    }
    
    public Kernel getKernel() {
        /*
                        |  1  2  1 |
           Soloblev_y = |  0  0  0 |
                        | -1 -2 -1 |
        
        this is the n=2 binomial filter for a Gaussian first derivative,
           that is sigma = sqrt(2)/2 = 0.707
        */
        Kernel kernel = new Kernel(3, 3);
       
        kernel.setValue(0, 0, 1);
        kernel.setValue(1, 0, 2);
        kernel.setValue(2, 0, 1);
        kernel.setValue(0, 2, -1);
        kernel.setValue(1, 2, -2);
        kernel.setValue(2, 2, -1);
        
        return kernel;
    }

}
