package algorithms.imageProcessing;

/**
 * class to hold an image instance and the original
 * size and offsets.
 * @author nichole
 */
public class TrimmedImage {
        
    private final int origWidth;
    private final int origHeight;
    
    private final int xOffset;
    private final int yOffset;
    
    private final Image imageTrimmed;
    
    private final boolean isImageExt;
    
    public TrimmedImage(Image img, int x0, int x1, int y0, int y1) {
        
        this.imageTrimmed = img.copySubImage(x0, x1, y0, y1);
        
        this.isImageExt = (img instanceof ImageExt);
        
        this.origWidth = img.getWidth();
        
        this.origHeight = img.getHeight();
        
        this.xOffset = x0;
        
        this.yOffset = y0;
    }
    
    public Image copyTrimmedToFullFrame() {
        
        Image img2;
        if (isImageExt) {
            img2 = new ImageExt(origWidth, origHeight);
        } else {
            img2 = new Image(origWidth, origHeight);
        }
        
        for (int i = 0; i < imageTrimmed.getWidth(); ++i) {
            for (int j = 0; j < imageTrimmed.getHeight(); ++j) {
                int x = i + xOffset;
                int y = j + yOffset;
                int rgb = imageTrimmed.getRGB(i, j);
                img2.setRGB(x, y, rgb);
            }
        }
        
        return img2;
    }

    /**
     * @return the width0rig
     */
    public int getOrigWidth() {
        return origWidth;
    }

    /**
     * @return the heightOrig
     */
    public int getOrigHeight() {
        return origHeight;
    }

    /**
     * @return the xOffset
     */
    public int getXOffset() {
        return xOffset;
    }

    /**
     * @return the yOffset
     */
    public int getYOffset() {
        return yOffset;
    }

    /**
     * @return the imageTrimmed
     */
    public Image getTrimmed() {
        return imageTrimmed;
    }
    
}
