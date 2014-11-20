package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class Transformer {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
     public PairIntArray[] applyTransformation(TransformationParameters params,
        PairIntArray[] edges) {
        
        double rotInRadians = params.getRotationInDegrees() * Math.PI/180.f;
        
        double scale = params.getScale();
        
        double translationX = params.getTranslationX();
        double translationY = params.getTranslationY();
        
        return applyTransformation(rotInRadians, scale, translationX,
            translationY, edges);
    }
    
    public PairIntArray[] applyTransformation(double rotInRadians,
        double scale, double translationX, double translationY,
        PairIntArray[] edges) {
        
        double cos = Math.cos(rotInRadians);
        double sin = Math.sin(rotInRadians);
                
        PairIntArray[] transformedEdges = new PairIntArray[edges.length];
        
        for (int ii = 0; ii < edges.length; ii++) {
        
            PairIntArray edge = edges[ii];
            
            PairIntArray te = new PairIntArray();
            
            for (int i = 0; i < edge.getN(); i++) {
                
                int x = edge.getX(i);
                int y = edge.getY(i);
                
                double xs = x * scale;
                double ys = y * scale;
                                
                double xt = (xs * cos) + (ys * sin);
                double yt = (-1. * xs * sin) + (ys * cos);
                
                xt += translationX;
                yt += translationY;
                
                te.add((int)xt, (int)yt);
                               
            }
            
            transformedEdges[ii] = te;
        }
        
        return transformedEdges;
    }
    
    public GreyscaleImage applyTransformation(final GreyscaleImage input, 
        TransformationParameters params, int outputWidth, int outputHeight) {
        
        double rotInRadians = params.getRotationInDegrees() * Math.PI/180.f;
        double cos = Math.cos(rotInRadians);
        double sin = Math.sin(rotInRadians);
        
        float scale = params.getScale();
        
        float translationX = params.getTranslationX();
        float translationY = params.getTranslationY();
        
        GreyscaleImage output = new GreyscaleImage(outputWidth, 
            outputHeight);
        
        for (int x = 0; x < input.getWidth(); x++) {
            for (int y = 0; y < input.getHeight(); y++) {
                
                int pix = input.getValue(x, y);
                
                if (pix != 0) {
                    
                    double xs = x * scale;
                    double ys = y * scale;

                    double xt = (xs * cos) + (ys * sin);
                    double yt = (-1. * xs * sin) + (ys * cos);

                    xt += translationX;
                    yt += translationY;
                    
                    int x2 = (int)xt;
                    int y2 = (int)yt;
                    
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
        
        Image output = new Image(outputWidth, outputHeight);
        
        for (int x = 0; x < input.getWidth(); x++) {
            for (int y = 0; y < input.getHeight(); y++) {
                
                double xs = x * scale;
                double ys = y * scale;
                                
                double xt = (xs * cos) + (ys * sin);
                double yt = (-1. * xs * sin) + (ys * cos);
                
                xt += translationX;
                yt += translationY;
                
                int r = input.getR(x, y);
                int g = input.getG(x, y);
                int b = input.getB(x, y);
                
                int x2 = (int)xt;
                int y2 = (int)yt;

                if ((x2 > -1) && (x2 < (output.getWidth() - 1) &&
                    (y2 > -1) && (y2 < (output.getHeight() - 1)))) {

                    output.setRGB(x2, y2, r, g, b);
                }
            }
        }
        
        return output;
    }
}
