package algorithms.imageProcessing;

import algorithms.matrix.MatrixUtil;

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
     * @param image 2 dimensional image in row-major format
     * @param sigma
     * @param create2ndDerivs 
     */
    public StructureTensor(float[][] image, float sigma, boolean create2ndDerivs) {

        this.sigma = sigma;
        
        // --- create Sobel derivatives (gaussian 1st deriv sqrt(2)/2 = 0.707)----
        ImageProcessor imageProcessor = new ImageProcessor();

        // switch X and Y sobel operations to match scipy (column major to row major or vice versa)

        // NOTE: may need to revisit this
        float norm = 0.5f;//(0.707f * (float)Math.sqrt(2. * Math.PI));
        
        float[][] gX = imageProcessor.copy(image);
        imageProcessor.applySobelX(gX);
        MatrixUtil.multiply(gX, norm);

        float[][] gY = imageProcessor.copy(image);
        imageProcessor.applySobelY(gY);
        MatrixUtil.multiply(gY, norm);

        dX = gX;
        dY = gY;

        // --- create structure tensors ----
        float[] kernel = (sigma > 0) ? Gaussian1D.getKernel(sigma) : null;
        
        //Axx
        dXSq = algorithms.imageProcessing.util.MatrixUtil.multiplyPointwise(gX, gX);
        if (kernel != null) {
            imageProcessor.applyKernelTwo1Ds(dXSq, kernel);
        }
        
        //Ayy
        dYSq = algorithms.imageProcessing.util.MatrixUtil.multiplyPointwise(gY, gY);
        if (kernel != null) {
            imageProcessor.applyKernelTwo1Ds(dYSq, kernel);
        }
        
        //Axy
        dXdY = algorithms.imageProcessing.util.MatrixUtil.multiplyPointwise(gX, gY);
        if (kernel != null) {
            imageProcessor.applyKernelTwo1Ds(dXdY, kernel);
        }
        
        if (create2ndDerivs) {

            // for curvature, need d/dy(dy) and d/dx(dx)
            kernel = (sigma > 0) ?
                Gaussian1DFirstDeriv.getKernel(sigma) :
                Gaussian1DFirstDeriv.getKernel(SIGMA.ZEROPOINTSEVENONE);
            
            d2X = imageProcessor.copy(gX);
            d2Y = imageProcessor.copy(gY);
        
            //TODO: revisit this in detail:
            // row major, so need to use y operations for x and vice versa
            imageProcessor.applyKernel1D(d2X, kernel, true);
            imageProcessor.applyKernel1D(d2Y, kernel, false);
        
        } else {
            d2X = null;
            d2Y = null;
        }
    }
   
    public float[][] getDeterminant() {
        
        if (detA == null) {
                        
            // detA = Axx * Ayy - Axy ** 2
            
            float[][] axxyy = algorithms.imageProcessing.util.MatrixUtil.multiplyPointwise(dXSq, dYSq);

            float[][] axyxy = algorithms.imageProcessing.util.MatrixUtil.multiplyPointwise(dXdY, dXdY);

            detA = MatrixUtil.subtract(axxyy, axyxy);
        }
        
        return detA;
    }
    
    public float[][] getTrace() {
        
        if (traceA == null) {
            traceA = MatrixUtil.add(dXSq, dYSq);
        }
        
        return traceA;
    }
    
    /**
     * get Axx
     * 
     * @return 
     */
    public float[][] getDXSquared() {
        return dXSq;
    }
    
    /**
     * get Ayy
     * @return 
     */
    public float[][] getDYSquared() {
        return dYSq;
    }
    
    /**
     * get Axy
     * @return 
     */
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
