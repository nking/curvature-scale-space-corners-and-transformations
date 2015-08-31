package algorithms.imageProcessing;

import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class Transformer {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * transform the given edges using the given parameters.
     * Note, the rotation is applied in a clockwise
      * direction (which is in the -z direction using right hand rule).
      * 
     * @param params
     * @param edges
     * @return 
     */
     public PairIntArray[] applyTransformation(TransformationParameters params,
        PairIntArray[] edges) {
        
        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
         
        double rotInRadians = params.getRotationInRadians();
        double scale = params.getScale();        
        double translationX = params.getTranslationX();
        double translationY = params.getTranslationY();
        double centroidX = params.getOriginX();
        double centroidY = params.getOriginY();
        
        return applyTransformation(rotInRadians, scale, translationX,
            translationY, centroidX, centroidY, edges);
     }
     
     /**
     * transform the given edges using the given parameters.
     * Note, the rotation is applied in a clockwise
      * direction (which is in the -z direction using right hand rule).
      * 
     * @param params
     * @param points
     * @return
     */
     public PairIntArray applyTransformation(TransformationParameters params,
        PairIntArray points) {
        
        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
         
        double rotInRadians = params.getRotationInRadians();
        double scale = params.getScale();        
        double translationX = params.getTranslationX();
        double translationY = params.getTranslationY();
        double centroidX = params.getOriginX();
        double centroidY = params.getOriginY();
        
        return applyTransformation(rotInRadians, scale, translationX,
            translationY, centroidX, centroidY, points);
     }
     
     /**
     * transform the given edges using the given parameters.
     * Note, the rotation is applied in a clockwise
      * direction (which is in the -z direction using right hand rule).
      * 
     * @param params
     * @param points
     * @return
     */
     public PairFloatArray applyTransformation2(TransformationParameters params,
        PairIntArray points) {
        
        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
        
        double rotInRadians = params.getRotationInRadians();
        double scale = params.getScale();        
        double translationX = params.getTranslationX();
        double translationY = params.getTranslationY();
        double centroidX = params.getOriginX();
        double centroidY = params.getOriginY();
        
        return applyTransformation2(rotInRadians, scale, translationX,
            translationY, centroidX, centroidY, points);
     }
    
     /**
      * transform the given edges using the given parameters. 
      * 
      * @param rotInRadians rotation in radians.  Note, the rotation is applied in a clockwise
      * direction (which is in the -z direction using right hand rule).
      * @param scale
      * @param translationX translation along x axis in pixels
      * @param translationY translation along y axis in pixels
      * @param centroidX the horizontal center of the reference frame for edges.  
      * this should be the center of the image if edges points are in the
      * original image reference frame.
      * @param centroidY the vertical center of the reference frame for edges.  
      * this should be the center of the image if edges points are in the
      * original image reference frame.
      * @param edges
      * @return 
      */
    public PairIntArray[] applyTransformation(double rotInRadians,
        double scale, double translationX, double translationY,
        double centroidX, double centroidY, PairIntArray[] edges) {
        
        if (edges == null) {
            throw new IllegalArgumentException("edges cannot be null");
        }
        
        PairIntArray[] transformedEdges = new PairIntArray[edges.length];

        for (int ii = 0; ii < edges.length; ii++) {

            PairIntArray edge = edges[ii];

            PairIntArray te = applyTransformation(rotInRadians, scale, 
                translationX, translationY, centroidX, centroidY, edge);

            transformedEdges[ii] = te;
        }
        
        return transformedEdges;
    }
    
    /**
      * transform the given edges using the given parameters. 
      * Note, the rotation is applied in a clockwise
      * direction (which is in the -z direction using right hand rule).
      * @param params
      * @param edges
      * @return 
      */
    public List<PairIntArray> applyTransformation(TransformationParameters
        params, List<PairIntArray> edges) {
        
        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
        if (edges == null) {
            throw new IllegalArgumentException("edges cannot be null");
        }
        
        List<PairIntArray> transformedEdges = new ArrayList<PairIntArray>();

        double scale = params.getScale();
        double rotInRadians = params.getRotationInRadians();
        double translationX = params.getTranslationX();
        double translationY = params.getTranslationY();
        double centroidX = params.getOriginX();
        double centroidY = params.getOriginY();
        
        for (PairIntArray edge : edges) {

            PairIntArray te = applyTransformation(rotInRadians, scale, 
                translationX, translationY, centroidX, centroidY, edge);

            transformedEdges.add(te);
        }
        
        return transformedEdges;
    }
    
    /**
      * transform the given edges using the given parameters. 
      * 
      * @param rotInRadians rotation in radians. Note, the rotation is applied in a clockwise
      * direction (which is in the -z direction using right hand rule).
      * @param scale
      * @param translationX translation along x axis in pixels
      * @param translationY translation along y axis in pixels
      * @param centroidX the horizontal center of the reference frame for edges.  
      * this should be the center of the image if edges points are in the
      * original image reference frame.
      * @param centroidY the vertical center of the reference frame for edges.  
      * this should be the center of the image if edges points are in the
      * original image reference frame.
      * @param edges
      * @return 
      */
    public List<PairIntArray> applyTransformation(double rotInRadians,
        double scale, double translationX, double translationY,
        double centroidX, double centroidY, List<PairIntArray> edges) {
        
        if (edges == null) {
            throw new IllegalArgumentException("edges cannot be null");
        }
        
        List<PairIntArray> transformedEdges = new ArrayList<PairIntArray>();

        for (PairIntArray edge : edges) {

            PairIntArray te = applyTransformation(rotInRadians, scale, 
                translationX, translationY, centroidX, centroidY, edge);

            transformedEdges.add(te);
        }
        
        return transformedEdges;
    }
    
     /**
      * transform the given edge using the given parameters. 
      * 
      * @param rotInRadians rotation in radians.  Note, the rotation is applied in a clockwise
      * direction (which is in the -z direction using right hand rule).
      * @param scale
      * @param translationX translation along x axis in pixels
      * @param translationY translation along y axis in pixels
      * @param centroidX the horizontal center of the reference frame for edges.  
      * this should be the center of the image if edges points are in the
      * original image reference frame.
      * @param centroidY the vertical center of the reference frame for edges.  
      * this should be the center of the image if edges points are in the
      * original image reference frame.
      * @param edge
      * @return 
      */
    public PairIntArray applyTransformation(double rotInRadians,
        double scale, double translationX, double translationY,
        double centroidX, double centroidY, PairIntArray edge) {
        
        if (edge == null) {
            throw new IllegalArgumentException("edge cannot be null");
        }
        
        double cos = Math.cos(rotInRadians);
        double sin = Math.sin(rotInRadians);
                
        /*
        scale, rotate, then translate.
        
        xr_0 = xc*scale + (((x0-xc)*scale*math.cos(theta)) + ((y0-yc)*scale*math.sin(theta)))

        xt_0 = xr_0 + transX = x1

        yr_0 = yc*scale + (-((x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))

        yt_0 = yr_0 + transY = y1
        */
        
        PairIntArray te = new PairIntArray();

        for (int i = 0; i < edge.getN(); i++) {

            double x = edge.getX(i);
            double y = edge.getY(i);

            double xr = centroidX * scale + ((x - centroidX) * scale * cos) 
                + ((y - centroidY) * scale * sin);

            double yr = centroidY * scale + (-(x - centroidX) * scale * sin) 
                + ((y - centroidY) * scale * cos);

            double xt = xr + translationX;
            double yt = yr + translationY;

            int xte = (int) Math.round(xt);
            int yte = (int) Math.round(yt);
            
            te.add(xte, yte);
        }
          
        return te;
    }
    
    /**
      * transform the given edge using the given parameters. 
      * 
      * @param rotInRadians rotation in radians. Note, the rotation is applied in a clockwise
      * direction (which is in the -z direction using right hand rule).
      * @param scale
      * @param translationX translation along x axis in pixels
      * @param translationY translation along y axis in pixels
      * @param centroidX the horizontal center of the reference frame for edges.  
      * this should be the center of the image if edges points are in the
      * original image reference frame.
      * @param centroidY the vertical center of the reference frame for edges.  
      * this should be the center of the image if edges points are in the
      * original image reference frame.
      * @param points
      * @return 
      */
    public PairFloatArray applyTransformation2(double rotInRadians,
        double scale, double translationX, double translationY,
        double centroidX, double centroidY, PairIntArray points) {
       
        if (points == null) {
            throw new IllegalArgumentException("points cannot be null");
        }
        
        double cos = Math.cos(rotInRadians);
        double sin = Math.sin(rotInRadians);
                
        /*
        scale, rotate, then translate.
        
        xr_0 = xc*scale + (((x0-xc)*scale*math.cos(theta)) + ((y0-yc)*scale*math.sin(theta)))

        xt_0 = xr_0 + transX = x1

        yr_0 = yc*scale + (-((x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))

        yt_0 = yr_0 + transY = y1
        */
        
        PairFloatArray te = new PairFloatArray();

        for (int i = 0; i < points.getN(); ++i) {

            double x = points.getX(i);
            double y = points.getY(i);

            double xr = centroidX * scale + ((x - centroidX) * scale * cos) 
                + ((y - centroidY) * scale * sin);

            double yr = centroidY * scale + (-(x - centroidX) * scale * sin) 
                + ((y - centroidY) * scale * cos);

            double xt = xr + translationX;
            double yt = yr + translationY;

            float xte = (float)(xt);
            float yte = (float)(yt);
            
            te.add(xte, yte);
        }
          
        return te;
    }
    
    /**
     * apply the parameters of transformation to point (xPt, yPt). Note, the rotation is applied in a clockwise
      * direction (which is in the -z direction using right hand rule).
     * 
     * @param params
     * @param xPt
     * @param yPt
     * @return 
     */
    public double[] applyTransformation(TransformationParameters params,
        double xPt, double yPt) {
        
        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
        
        double scale = params.getScale();
        double scaleTimesCosine = scale * Math.cos(params.getRotationInRadians());
        double scaleTimesSine = scale * Math.sin(params.getRotationInRadians());
                
        double originX = params.getOriginX();
        double originY = params.getOriginY();
        
        /*
        xr_0 = xc*scale + (((x0-xc)*scale*math.cos(theta)) + ((y0-yc)*scale*math.sin(theta)))

        xt_0 = xr_0 + transX = x1

        yr_0 = yc*scale + (-((x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))

        yt_0 = yr_0 + transY = y1
        */
        double xr = originX * scale + ((xPt - originX) * scaleTimesCosine) 
            + ((yPt - originY) * scaleTimesSine);

        double yr = originY * scale + (-(xPt - originX) * scaleTimesSine) 
            + ((yPt - originY) * scaleTimesCosine);

        double xt = xr + params.getTranslationX();
        double yt = yr + params.getTranslationY();

        return new double[]{xt, yt};
    }
    
    /**
     * apply the transformation params to set1 and populate outputTrX and
     * outputTrY.
     * Note, the rotation is applied in a clockwise
      * direction (which is in the -z direction using right hand rule).
     * @param params
     * @param set1
     * @param outputTrX
     * @param outputTrY 
     */
    public void applyTransformation(TransformationParameters params,
        PairIntArray set1, float[] outputTrX, float[] outputTrY) {
        
        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
        if (set1 == null) {
            throw new IllegalArgumentException("set1 cannot be null");
        }
        if (outputTrX == null) {
            throw new IllegalArgumentException("outputTrX cannot be null");
        }
        if (outputTrY == null) {
            throw new IllegalArgumentException("outputTrY cannot be null");
        }
        if (set1.getN() != outputTrX.length) {
            throw new IllegalArgumentException(
                "outputTrX must the same length as set1");
        }
        if (set1.getN() != outputTrY.length) {
            throw new IllegalArgumentException(
                "outputTrY must the same length as set1");
        }
        
        float scale = params.getScale();
        float scaleTimesCosine = (float)
            (scale * Math.cos(params.getRotationInRadians()));
        float scaleTimesSine = (float)
            (scale * Math.sin(params.getRotationInRadians()));
                
        float originX = params.getOriginX();
        float originY = params.getOriginY();
        
        float transX = params.getTranslationX();
        float transY = params.getTranslationY();
        
        /*
        xr_0 = xc*scale + (((x0-xc)*scale*math.cos(theta)) + ((y0-yc)*scale*math.sin(theta)))

        xt_0 = xr_0 + transX = x1

        yr_0 = yc*scale + (-((x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))

        yt_0 = yr_0 + transY = y1
        */
        
        for (int i = 0; i < set1.getN(); ++i) {
        
            int xPt = set1.getX(i);
            int yPt = set1.getY(i);
            
            float xr = originX * scale + ((xPt - originX) * scaleTimesCosine) 
                + ((yPt - originY) * scaleTimesSine);

            float yr = originY * scale + (-(xPt - originX) * scaleTimesSine) 
                + ((yPt - originY) * scaleTimesCosine);

            float xt = xr + transX;
            float yt = yr + transY;

            outputTrX[i] = xt;
            outputTrY[i] = yt;
        }
    }
    
    /**
     * apply transformation params to the image input.  Note, the rotation is applied in a clockwise
      * direction (which is in the -z direction using right hand rule).
     * @param input
     * @param params
     * @param outputWidth
     * @param outputHeight
     * @return 
     */
    public GreyscaleImage applyTransformation(final GreyscaleImage input, 
        TransformationParameters params, int outputWidth, int outputHeight) {
        
        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
        
        double rotInRadians = params.getRotationInRadians();
        double cos = Math.cos(rotInRadians);
        double sin = Math.sin(rotInRadians);
        
        float scale = params.getScale();
        
        float translationX = params.getTranslationX();
        float translationY = params.getTranslationY();
        
        double centroidX = params.getOriginX();
        double centroidY = params.getOriginY();
        
        /*
        xr_0 = xc*scale + (((x0-xc)*scale*math.cos(theta)) + ((y0-yc)*scale*math.sin(theta)))

        xt_0 = xr_0 + transX = x1

        yr_0 = yc*scale + (-((x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))

        yt_0 = yr_0 + transY = y1
        */
        
        GreyscaleImage output = new GreyscaleImage(outputWidth, 
            outputHeight);
        
        for (int x = 0; x < input.getWidth(); x++) {
            for (int y = 0; y < input.getHeight(); y++) {
                
                int pix = input.getValue(x, y);
                
                if (pix != 0) {
                    
                    double xr = centroidX * scale + ((x - centroidX) *scale * cos) 
                        + ((y - centroidY) * scale * sin);

                    double yr = centroidY * scale + (-(x - centroidX) * scale * sin)
                        + ((y - centroidY) * scale * cos);

                    double xt = xr + translationX;
                    double yt = yr + translationY;
                
                    int x2 = (int)Math.round(xt);
                    int y2 = (int)Math.round(yt);
                    
                    if ((x2 > -1) && (x2 < (output.getWidth() - 1) &&
                        (y2 > -1) && (y2 < (output.getHeight() - 1)))) {
                        
                        output.setValue(x2, y2, pix);
                    }
                }
            }
        }
        
        return output;
    }
    
    /**
     * apply transformation params to image input. Note, the rotation is applied in a clockwise
      * direction (which is in the -z direction using right hand rule).
     * @param input
     * @param params
     * @param outputWidth
     * @param outputHeight
     * @return 
     */
    public Image applyTransformation(Image input, 
        TransformationParameters params, int outputWidth, int outputHeight) {
        
        if (input == null) {
            throw new IllegalArgumentException("input cannot be null");
        }
        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
        
        double rotInRadians = params.getRotationInDegrees() * Math.PI/180.f;
        double cos = Math.cos(rotInRadians);
        double sin = Math.sin(rotInRadians);
        
        float scale = params.getScale();
        
        float translationX = params.getTranslationX();
        float translationY = params.getTranslationY();
        
        double centroidX = params.getOriginX();
        double centroidY = params.getOriginY();
        
        /*
        xr_0 = xc*scale + (((x0-xc)*scale*math.cos(theta)) + ((y0-yc)*scale*math.sin(theta)))

        xt_0 = xr_0 + transX = x1

        yr_0 = yc*scale + (-((x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))

        yt_0 = yr_0 + transY = y1
        */
        
        boolean isExt = (input instanceof ImageExt);
        
        Image output = isExt ? new ImageExt(outputWidth, outputHeight) :
            new Image(outputWidth, outputHeight);
        
        for (int x = 0; x < input.getWidth(); x++) {
            for (int y = 0; y < input.getHeight(); y++) {
                
                double xr = centroidX * scale + ((x - centroidX) * scale * cos)
                    + ((y - centroidY) * scale * sin);

                double yr = centroidY * scale + (-(x - centroidX) * scale * sin)
                    + ((y - centroidY) * scale * cos);

                double xt = xr + translationX;
                double yt = yr + translationY;
                
                int x2 = (int)Math.round(xt);
                int y2 = (int)Math.round(yt);

                if ((x2 < 0) || (x2 > (outputWidth - 1)) ||
                    (y2 < 0) || (y2 > (outputHeight - 1))) {
                    continue;
                }
                
                int r = input.getR(x, y);
                int g = input.getG(x, y);
                int b = input.getB(x, y);

                output.setRGB(x2, y2, r, g, b);

                if (isExt) {
                    // TODO: copy any already set arrays too
                }
            }
        }
        
        return output;
    }
    
    /**
     * apply the scale transformation to the existing transformation 
     * parameters params and return the result.
     * @param params
     * @param scale
     * @return 
     */
    public TransformationParameters applyScaleTransformation(
        TransformationParameters params, float scale) {
        
        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
        
        /*
        xr_0 = xc*scale + (((x0-xc)*scale*math.cos(theta)) + ((y0-yc)*scale*math.sin(theta)))

        xt_0 = xr_0 + transX = x1

        yr_0 = yc*scale + (-((x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))

        yt_0 = yr_0 + transY = y1
        */
        
        TransformationParameters tr2 = params.copy();
        
        float ox = params.getOriginX();
        float oy = params.getOriginY();
        float tx = params.getTranslationX();
        float ty = params.getTranslationY();
        double cosine = Math.cos(0);
        double sine = Math.sin(0);
        
        tr2.setScale(scale*params.getScale());
        tr2.setOriginX(ox * scale);
        tr2.setOriginY(oy * scale);
        
        double tx2 = (ox * scale) + ((tx - ox) * scale * cosine) 
            + ((ty - oy) * scale * sine);
        double ty2 = (oy * scale) + (-(tx - ox) * scale * sine)
            + ((ty - oy) * scale * cosine);
        
        tr2.setTranslationX((float)tx2);
        tr2.setTranslationY((float)ty2);
        
        return tr2;
    }
    
    /**
     * given a two-dimensional array of x and y, apply rotation to them
     * and return the result.  Note, the rotation is applied in a clockwise
      * direction (which is in the -z direction using right hand rule).
     * @param rotationInDegrees
     * @param xy
     * @return 
     */
    public float[][] transformXY(float rotationInDegrees, float[][] xy) {
        
        if (xy == null) {
            throw new IllegalArgumentException("xy cannot be null");
        }
        if (xy.length == 0) {
            return new float[0][];
        }
        if (xy[0].length != 2) {
            throw new IllegalArgumentException(
            "xy must be size 2 for the 2nd dimension");
        }
        
        double rotationInRadians = rotationInDegrees * Math.PI/180.;
        
        double cos = Math.cos(rotationInRadians);
        double sin = Math.sin(rotationInRadians);
        
        float[][] transformed = new float[xy.length][];
        /*
        xr_0 = xc*scale + (((x0-xc)*scale*math.cos(theta)) + ((y0-yc)*scale*math.sin(theta)))

        xt_0 = xr_0 + transX = x1

        yr_0 = yc*scale + (-((x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))

        yt_0 = yr_0 + transY = y1
        */
        for (int i = 0; i < xy.length; ++i) {
            
            float x = xy[i][0];
            float y = xy[i][1];
            
            double xr = (x * cos) + (y * sin);
            double yr = (-x * sin) + (y * cos);
            
            transformed[i] = new float[]{(float)xr, (float)yr};
        }
        
        return transformed;
    }
    
    /**
     * given a two-dimensional array of x and y, apply rotation to them
     * and return the result. Note, the rotation is applied in a clockwise
      * direction (which is in the -z direction using right hand rule).
     * @param rotationInDegrees
     * @param x
     * @param y
     * @return 
     */
    public float[] rotate(float rotationInDegrees, float x, float y) {
        
        double rotationInRadians = rotationInDegrees * Math.PI/180.;
        
        double cos = Math.cos(rotationInRadians);
        double sin = Math.sin(rotationInRadians);
                            
        double xr = (x * cos) + (y * sin);
        double yr = (-x * sin) + (y * cos);
            
        return new float[]{(float)xr, (float)yr};
    }

    /**
     * apply rotationInDegrees to (x, y) and return the result in the arrays
     * outputX and outputY.
     * Note, the rotation is applied in a clockwise
      * direction (which is in the -z direction using right hand rule).
      * 
     * @param rotationInDegrees
     * @param x
     * @param y
     * @param outputX
     * @param outputY
     * @param outputIndex 
     */
    public void rotate(int rotationInDegrees, int x, int y, 
        float[] outputX, float[] outputY, int outputIndex) {
        
        if (outputX == null) {
            throw new IllegalArgumentException("outputX cannot be null");
        }
        if (outputY == null) {
            throw new IllegalArgumentException("outputY cannot be null");
        }
        if (outputX.length != outputY.length) {
            throw new IllegalArgumentException(
            "outputX.length must be same length as outputY.length");
        }
        
        if (outputIndex < 0 || (outputIndex > (outputX.length - 1))) {
            throw new IllegalArgumentException(
            "outputIndex is out of bounds of outputX");
        }
        
        double rotationInRadians = rotationInDegrees * Math.PI/180.;
        
        double cos = Math.cos(rotationInRadians);
        double sin = Math.sin(rotationInRadians);
                            
        double xr = (x * cos) + (y * sin);
        double yr = (-x * sin) + (y * cos);
        
        outputX[outputIndex] = (float)xr;
        outputY[outputIndex] = (float)yr;
        
    }
    
}
