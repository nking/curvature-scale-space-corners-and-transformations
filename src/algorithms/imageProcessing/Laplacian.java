package algorithms.imageProcessing;

import java.util.Arrays;

/** 
 * @author nichole
 */
public class Laplacian implements IKernel {
     
    public float getNormalizationFactor() {
        return 1.0f;
    }
    
    public Kernel getKernel() {
        /*
          [ -1 -1 -1 ]
          [ -1  8 -1 ]
          [ -1 -1 -1 ]
        
        */
        Kernel kernel = new Kernel(3, 3);
        Arrays.fill(kernel.a, -1);
        kernel.setValue(1, 1,  0);
        
        return kernel;
    }
}
