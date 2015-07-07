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
     * 
     * @param params
     * @param edges
     * @param centroidX the horizontal center of the reference frame for edges.  
     * this should be the center of the image if edges points are in the
     * original image reference frame.
     * @param centroidY the vertical center of the reference frame for edges.  
     * this should be the center of the image if edges points are in the
     * original image reference frame.
     * @return 
     */
     public PairIntArray[] applyTransformation(TransformationParameters params,
        PairIntArray[] edges, double centroidX, double centroidY) {
        
        double rotInRadians = params.getRotationInRadians();
        double scale = params.getScale();        
        double translationX = params.getTranslationX();
        double translationY = params.getTranslationY();
        
        return applyTransformation(rotInRadians, scale, translationX,
            translationY, centroidX, centroidY, edges);
     }
     
     /**
     * transform the given edges using the given parameters.
     * 
     * @param params
     * @param edge
     * @param centroidX the horizontal center of the reference frame for edges.  
     * this should be the center of the image if edges points are in the
     * original image reference frame.
     * @param centroidY the vertical center of the reference frame for edges.  
     * this should be the center of the image if edges points are in the
     * original image reference frame.
     * @return
     */
     public PairIntArray applyTransformation(TransformationParameters params,
        PairIntArray edge, double centroidX, double centroidY) {
        
        double rotInRadians = params.getRotationInRadians();
        double scale = params.getScale();        
        double translationX = params.getTranslationX();
        double translationY = params.getTranslationY();
        
        return applyTransformation(rotInRadians, scale, translationX,
            translationY, centroidX, centroidY, edge);
     }
     
     /**
     * transform the given edges using the given parameters.
     * 
     * @param params
     * @param edge
     * @param centroidX the horizontal center of the reference frame for edges.  
     * this should be the center of the image if edges points are in the
     * original image reference frame.
     * @param centroidY the vertical center of the reference frame for edges.  
     * this should be the center of the image if edges points are in the
     * original image reference frame.
     * @return
     */
     public PairFloatArray applyTransformation2(TransformationParameters params,
        PairIntArray edge, double centroidX, double centroidY) {
        
        double rotInRadians = params.getRotationInRadians();
        double scale = params.getScale();        
        double translationX = params.getTranslationX();
        double translationY = params.getTranslationY();
        
        return applyTransformation2(rotInRadians, scale, translationX,
            translationY, centroidX, centroidY, edge);
     }
    
     /**
      * transform the given edges using the given parameters. 
      * 
      * @param rotInRadians rotation in radians
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
      * 
      * @param params
      * @param centroidX the horizontal center of the reference frame for edges.  
      * this should be the center of the image if edges points are in the
      * original image reference frame.
      * @param centroidY the vertical center of the reference frame for edges.  
      * this should be the center of the image if edges points are in the
      * original image reference frame.
      * @param edges
      * @return 
      */
    public List<PairIntArray> applyTransformation(TransformationParameters
        params, List<PairIntArray> edges, double centroidX, double centroidY) {
        
        List<PairIntArray> transformedEdges = new ArrayList<PairIntArray>();

        double scale = params.getScale();
        double rotInRadians = params.getRotationInRadians();
        double translationX = params.getTranslationX();
        double translationY = params.getTranslationY();
        
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
      * @param rotInRadians rotation in radians
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
      * @param rotInRadians rotation in radians
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
        
        double cos = Math.cos(rotInRadians);
        double sin = Math.sin(rotInRadians);
                
        /*
        scale, rotate, then translate.
        */
        
        PairIntArray te = new PairIntArray();

        for (int i = 0; i < edge.getN(); i++) {

            double x = edge.getX(i);
            double y = edge.getY(i);

            double xr = centroidX * scale 
                + (((x - centroidX) *scale * cos) 
                + ((y - centroidY) * scale * sin));

            double yr = centroidY * scale 
                + ((-(x - centroidX) * scale * sin) 
                + ((y - centroidY) * scale * cos));

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
      * @param rotInRadians rotation in radians
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
       
        double cos = Math.cos(rotInRadians);
        double sin = Math.sin(rotInRadians);
                
        /*
        scale, rotate, then translate.
        */
        
        PairFloatArray te = new PairFloatArray();

        for (int i = 0; i < points.getN(); ++i) {

            double x = points.getX(i);
            double y = points.getY(i);

            double xr = centroidX * scale 
                + (((x - centroidX) * scale * cos) 
                + ((y - centroidY) * scale * sin));

            double yr = centroidY * scale 
                + ((-(x - centroidX) * scale * sin) 
                + ((y - centroidY) * scale * cos));

            double xt = xr + translationX;
            double yt = yr + translationY;

            float xte = (float)(xt);
            float yte = (float)(yt);
            
            te.add(xte, yte);
        }
          
        return te;
    }
    
    /**
      * transform the point using the given parameters. 
      * 
      * @param params euclidean transformation parameters
      * @param centroidX the horizontal center of the reference frame for edges.  
      * this should be the center of the image if edges points are in the
      * original image reference frame.
      * @param centroidY the vertical center of the reference frame for edges.  
      * this should be the center of the image if edges points are in the
      * original image reference frame.
      * @param xPt
      * @param yPt
      * @return 
      */
    public double[] applyTransformation(TransformationParameters params,
        double centroidX, double centroidY, double xPt, double yPt) {
        
        double scale = params.getScale();
        double scaleTimesCosine = scale * Math.cos(params.getRotationInRadians());
        double scaleTimesSine = scale * Math.sin(params.getRotationInRadians());
                
        double xr = centroidX * scale + (((xPt - centroidX) * scaleTimesCosine) 
            + ((yPt - centroidY) * scaleTimesSine));

        double yr = centroidY * scale 
            + ((-(xPt - centroidX) * scaleTimesSine) 
            + ((yPt - centroidY) * scaleTimesCosine));

        double xt = xr + params.getTranslationX();
        double yt = yr + params.getTranslationY();

        return new double[]{xt, yt};
    }
    
    public GreyscaleImage applyTransformation(final GreyscaleImage input, 
        TransformationParameters params, int outputWidth, int outputHeight) {
        
        double rotInRadians = params.getRotationInRadians();
        double cos = Math.cos(rotInRadians);
        double sin = Math.sin(rotInRadians);
        
        float scale = params.getScale();
        
        float translationX = params.getTranslationX();
        float translationY = params.getTranslationY();
        
        double centroidX = input.getWidth() >> 1;
        double centroidY = input.getHeight() >> 1;
        
        GreyscaleImage output = new GreyscaleImage(outputWidth, 
            outputHeight);
        
        for (int x = 0; x < input.getWidth(); x++) {
            for (int y = 0; y < input.getHeight(); y++) {
                
                int pix = input.getValue(x, y);
                
                if (pix != 0) {
                    
                    double xr = centroidX * scale 
                        + (((x - centroidX) *scale * cos) 
                        + ((y - centroidY) * scale * sin));

                    double yr = centroidY * scale
                        + ((-(x - centroidX) * scale * sin)
                        + ((y - centroidY) * scale * cos));

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
    
    public Image applyTransformation(Image input, 
        TransformationParameters params, int outputWidth, int outputHeight) {
        
        double rotInRadians = params.getRotationInDegrees() * Math.PI/180.f;
        double cos = Math.cos(rotInRadians);
        double sin = Math.sin(rotInRadians);
        
        float scale = params.getScale();
        
        float translationX = params.getTranslationX();
        float translationY = params.getTranslationY();
        
        double centroidX = input.getWidth() >> 1;
        double centroidY = input.getHeight() >> 1;
        
        Image output = new Image(outputWidth, outputHeight);
        
        for (int x = 0; x < input.getWidth(); x++) {
            for (int y = 0; y < input.getHeight(); y++) {
                
                int r = input.getR(x, y);
                int g = input.getG(x, y);
                int b = input.getB(x, y);
                
                double xr = centroidX * scale 
                    + (((x - centroidX) * scale * cos)
                    + ((y - centroidY) * scale * sin));

                double yr = centroidY * scale
                    + ((-(x - centroidX) * scale * sin)
                    + ((y - centroidY) * scale * cos));

                double xt = xr + translationX;
                double yt = yr + translationY;
                
                int x2 = (int)Math.round(xt);
                int y2 = (int)Math.round(yt);

                if ((x2 > -1) && (x2 < (output.getWidth() - 1) &&
                    (y2 > -1) && (y2 < (output.getHeight() - 1)))) {

                    output.setRGB(x2, y2, r, g, b);
                }
            }
        }
        
        return output;
    }
    
    public PairFloatArray applyTransformation(TransformationParameters params,
        int centroidX1, int centroidY1, PairIntArray input) {
        
        PairFloatArray output = new PairFloatArray();
        
        float rotation = params.getRotationInRadians();
        float scale = params.getScale();
        float scaleTimesCosine = (float)(scale * Math.cos(rotation));
        float scaleTimesSine = (float)(scale * Math.sin(rotation));
        float transX = params.getTranslationX();
        float transY = params.getTranslationY();
                        
        for (int i = 0; i < input.getN(); i++) {
            
            int x = input.getX(i);
            int y = input.getY(i);
            
            double transformedX = (centroidX1*scale + ( 
                ((x - centroidX1) * scaleTimesCosine) +
                ((y - centroidY1) * scaleTimesSine))) + transX;
            
            double transformedY = (centroidY1*scale + ( 
                (-(x - centroidX1) * scaleTimesSine) +
                ((y - centroidY1) * scaleTimesCosine))) + transY;
          
            int xt = (int)Math.round(transformedX);
            
            int yt = (int)Math.round(transformedY);
            
            output.add(xt, yt);
        }
        
        return output;
    }
}
