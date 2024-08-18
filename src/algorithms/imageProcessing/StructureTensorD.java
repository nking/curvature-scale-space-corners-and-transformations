package algorithms.imageProcessing;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath0;

/**
 * create first derivative products and optionally second derivative
 * products of an image.
 * tailored to row-major image format (i.e. image[row][col])
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
public class StructureTensorD {

    private final float sigma;

    private final double[][] dX;

    private final double[][] dXSq;

    private final double[][] dY;

    private final double[][] dYSq;

    private final double[][] dXdY;

    // optional:
    private final double[][] d2X;
    private final double[][] d2Y;

    // created on demand
    private double[][] detA = null;
    private double[][] traceA = null;

    /**
     * create sobel x and y derivatives (== first deriv gaussians w/ sigma =
     * sqrt(2)/2) followed by smoothing with a gaussian of sigma=given sigma.
     * If create2ndDerivs is true, the second derivatives needed for curvature
     * are calculated too.
     *
     * Note that column-major data is assumed, so if have row-major data and context, then
     * use the opposite variable when accessing results.  e.g. if image context is row major,
     * then use dY = getDX() and dX = getDY() etc.
     *
     * @param image assumed to be row major format
     * @param sigma
     * @param create2ndDerivs
     */
    public StructureTensorD(double[][] image, float sigma, boolean create2ndDerivs) {

        this.sigma = sigma;
        
        // --- create Sobel derivatives (gaussian 1st deriv sqrt(2)/2 = 0.707)----
        ImageProcessor imageProcessor = new ImageProcessor();

        // switch X and Y sobel operations to match scipy (column major to row major or vice versa)

        // NOTE: may need to revisit this
        double norm = 0.5;//(0.707f * Math.sqrt(2. * Math.PI));
        
        double[][] gX = imageProcessor.copy(image);
        imageProcessor.applySobelX(gX);
        MatrixUtil.multiply(gX, norm);

        double[][] gY = imageProcessor.copy(image);
        imageProcessor.applySobelY(gY);
        MatrixUtil.multiply(gY, norm);
        
        // --- create structure tensors ----
        double[] kernel = (sigma > 0) ? MiscMath0.convertFloatToDouble(Gaussian1D.getKernel(sigma)) : null;
        
        //Axx
        dXSq = MatrixUtil.pointwiseMultiplication(gX, gX);
        if (kernel != null) {
            imageProcessor.applyKernelTwo1Ds(dXSq, kernel);
        }
        
        //Ayy
        dYSq = MatrixUtil.pointwiseMultiplication(gY, gY);
        if (kernel != null) {
            imageProcessor.applyKernelTwo1Ds(dYSq, kernel);
        }
        
        //Axy
        dXdY = MatrixUtil.pointwiseMultiplication(gX, gY);
        if (kernel != null) {
            imageProcessor.applyKernelTwo1Ds(dXdY, kernel);
        }

        dX = gX;
        dY = gY;

        if (create2ndDerivs) {

            // for curvature, need d/dy(dy) and d/dx(dx)
            kernel = (sigma > 0) ?
                    MiscMath0.convertFloatToDouble(Gaussian1DFirstDeriv.getKernel(sigma)) :
                    MiscMath0.convertFloatToDouble(Gaussian1DFirstDeriv.getKernel(SIGMA.ZEROPOINTSEVENONE));
            
            d2X = imageProcessor.copy(gX);
            d2Y = imageProcessor.copy(gY);
        
            imageProcessor.applyKernel1D(d2X, kernel, true);
            imageProcessor.applyKernel1D(d2Y, kernel, false);
        
        } else {
            d2X = null;
            d2Y = null;
        }
    }
   
    public double[][] getDeterminant() {
        
        if (detA == null) {
                        
            // detA = Axx * Ayy - Axy ** 2
            
            double[][] axxyy = MatrixUtil.pointwiseMultiplication(dXSq, dYSq);

            double[][] axyxy = MatrixUtil.pointwiseMultiplication(dXdY, dXdY);

            detA = MatrixUtil.pointwiseSubtract(axxyy, axyxy);
        }
        
        return detA;
    }
    
    public double[][] getTrace() {
        
        if (traceA == null) {
            traceA = MatrixUtil.pointwiseAdd(dXSq, dYSq);
        }
        
        return traceA;
    }
    
    /**
     * get Axx
     * 
     * @return 
     */
    public double[][] getDXSquared() {
        return dXSq;
    }
    
    /**
     * get Ayy
     * @return 
     */
    public double[][] getDYSquared() {
        return dYSq;
    }
    
    /**
     * get Axy
     * @return 
     */
    public double[][] getDXDY() {
        return dXdY;
    }
    
    public double[][] getDX() {
        return dX;
    }
    
    public double[][] getDY() {
        return dY;
    }
    
    public double[][] getDDX() {
        return d2X;
    }
    
    public double[][] getDDY() {
        return d2Y;
    }
}
