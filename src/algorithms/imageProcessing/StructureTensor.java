package algorithms.imageProcessing;

/**
 * create first derivative products and optionally second derivative
 * products of an image.
 * 
 * some amount of the code below has been adapted from
     https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/corner.py
     and replaced with existing local project functions.
   The scikit-image license is at
   https://github.com/scikit-image/scikit-image/blob/master/LICENSE.txt
   code is freely usable in any form with the copyright notice.
   -- begin scikit-image copyright --
   Unless otherwise specified by LICENSE.txt files in individual
    directories, all code is

    Copyright (C) 2011, the scikit-image team
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in
        the documentation and/or other materials provided with the
        distribution.
     3. Neither the name of skimage nor the names of its contributors may be
        used to endorse or promote products derived from this software without
        specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
    IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
    INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
    IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
   -- end scikit-image copyright -- 
   
 * @author nichole
 */
public class StructureTensor {
    
    private final float sigma;
    
    private final float[][] dXSq;
    
    private final float[][] dYSq;
    
    private final float[][] dXdY;
    
    // optional:
    private final float[][] dX;
    private final float[][] d2X;
    private final float[][] dY;
    private final float[][] d2Y;
    
    // created on demand
    private float[][] detA = null;
    private float[][] traceA = null;
    
    /**
     * create sobel x and y derivatives (== first deriv gaussians w/ sigma =
     * sqrt(2)/2) followed by smoothing with a gaussian of sigma=given sigma.
     * If create2ndDerivs is true, the second derivatives needed for curvature
     * are calculated too.
     * 
     * @param image
     * @param sigma
     * @param create2ndDerivs 
     */
    public StructureTensor(float[][] image, float sigma, boolean create2ndDerivs) {
        
        this.sigma = sigma;
        
        // --- create Sobel derivatives (gaussian 1st deriv sqrt(2)/2 = 0.707)----
        ImageProcessor imageProcessor = new ImageProcessor();

        // switch X and Y sobel operations to match scipy

        float[][] gX = imageProcessor.copy(image);
        imageProcessor.applySobelY(gX);

        float[][] gY = imageProcessor.copy(image);
        imageProcessor.applySobelX(gY);
        
        // --- create structure tensors ----
        float[] kernel = Gaussian1D.getKernel(sigma);
        
        dXSq = imageProcessor.multiply(gX, gX);
        imageProcessor.applyKernelTwo1Ds(dXSq, kernel);

        dYSq = imageProcessor.multiply(gY, gY);
        imageProcessor.applyKernelTwo1Ds(dYSq, kernel);

        dXdY = imageProcessor.multiply(gX, gY);
        imageProcessor.applyKernelTwo1Ds(dXdY, kernel);

        if (create2ndDerivs) {
            
            dX = gX;
            dY = gY;
            
            // for curvature, need d/dy(dy) and d/dx(dx)
            kernel = Gaussian1DFirstDeriv.getKernel(sigma);
            
            d2X = imageProcessor.copy(gX);
            d2Y = imageProcessor.copy(gY);
        
            //TODO: revisit this in detail:
            // row major, so need to use y operations for x and vice versa
            imageProcessor.applyKernel1D(d2X, kernel, false);
            imageProcessor.applyKernel1D(d2Y, kernel, true);
        
        } else {
            dX = null;
            d2X = null;
            dY = null;
            d2Y = null;
        }
    }
   
    public float[][] getDeterminant() {
        
        if (detA == null) {
            
            ImageProcessor imageProcessor = new ImageProcessor();
            float[][] axxyy = imageProcessor.multiply(dXSq, dYSq);

            float[][] axyxy = imageProcessor.multiply(dXdY, dXdY);

            detA = imageProcessor.subtract(axxyy, axyxy);
        }
        
        return detA;
    }
    
    public float[][] getTrace() {
        
        if (traceA == null) {
            ImageProcessor imageProcessor = new ImageProcessor();
            traceA = imageProcessor.add(dXSq, dXdY);
        }
        
        return traceA;
    }
    
    public float[][] getDXSquared() {
        return dXSq;
    }
    
    public float[][] getDYSquared() {
        return dYSq;
    }
    
    public float[][] getDXDY() {
        return dXdY;
    }
    
    public float[][] getDX() {
        return dX;
    }
    
    public float[][] getDY() {
        return dY;
    }
    
    public float[][] getDDX() {
        return d2X;
    }
    
    public float[][] getDDY() {
        return d2Y;
    }
}
