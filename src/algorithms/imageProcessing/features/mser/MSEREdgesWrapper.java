package algorithms.imageProcessing.features.mser;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageProcessor;

/**
 * a convenience class to use MSEREdges.  It reduces the image to a resolution
 * which the algorithm has been currently tailored for.
 * 
 * @author nichole
 */
public class MSEREdgesWrapper {
    
    private int maxDimension = 256;
    
    private int binFactor = 1;
    
    private boolean debug = false;
    
    public MSEREdges extractAndMergeEdges(ImageExt img) {
        
        img = binImage(img);
        
        MSEREdges mse = new MSEREdges(img);
        
        if (debug) {
            mse.setToDebug();
        }
        
        mse.extractAndMergeEdges();
        
        return mse;
    }

    public MSEREdges extractAndMergeEdges(ImageExt img, int maxDimension) {

        this.maxDimension = maxDimension;

        img = binImage(img);

        MSEREdges mse = new MSEREdges(img);

        if (debug) {
            mse.setToDebug();
        }

        mse.extractAndMergeEdges();

        return mse;
    }
    
    public MSEREdges extractEdges(ImageExt img) {
        
        img = binImage(img);
        
        MSEREdges mse = new MSEREdges(img);
        
        if (debug) {
            mse.setToDebug();
        }
        
        mse.extractEdges();
        
        return mse;
    }
    
    public void setToDebug() {
        debug = true;
    }
    
    public ImageExt binImage(ImageExt img) {
        
        int w = img.getWidth();        
        int h = img.getHeight();

        binFactor = (int) Math.ceil(Math.max((float) w / maxDimension,
            (float) h / maxDimension));

        ImageProcessor imageProcessor = new ImageProcessor();

        img = imageProcessor.binImage(img, binFactor);
        
        return img;
    }
    
}
