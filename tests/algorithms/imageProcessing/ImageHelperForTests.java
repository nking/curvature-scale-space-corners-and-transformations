package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class ImageHelperForTests {
    
    private GreyscaleImage img;
    
    private EdgeFilterProducts edgeFilterProducts = null;
    
    private int trimmedXOffset = 0;
    
    private int trimmedYOffset = 0;
    
    private boolean useOutdoorMode = false;
    
    private boolean useLineDrawingMode = false;
    
    public ImageHelperForTests(ImageExt input) {
        
        img = input.copyToGreyscale();
        
        init();
    }
    
    public ImageHelperForTests(ImageExt input, boolean useOutdoorMode) {
        
        this.useOutdoorMode = useOutdoorMode;
        
        img = input.copyToGreyscale();
        
        init();
    }
    
    public ImageHelperForTests(ImageExt input, boolean useOutdoorMode, 
        boolean useLineDrawingMode) {
        
        this.useOutdoorMode = useOutdoorMode;
        
        this.useLineDrawingMode = useLineDrawingMode;
        
        img = input.copyToGreyscale();
        
        init();
    } 
        
    private void init() {
        
        ImageProcessor ImageProcessor = new ImageProcessor();
        
        int[] offsetXY = ImageProcessor.shrinkImageToFirstNonZeros(img);
        trimmedXOffset = offsetXY[0];
        trimmedYOffset = offsetXY[1];
        
        CannyEdgeFilterAdaptive filter = new CannyEdgeFilterAdaptive();
        
        CannyEdgeFilterSettings settings = getCannyEdgeFilterSettings();
        
        filter.setSetters(settings);
                
        filter.applyFilter(img);
        
        edgeFilterProducts = filter.getFilterProducts();
        
    }
    
    /**
     * @return the gradientXY
     */
    public GreyscaleImage getGradientXY() {
        if (edgeFilterProducts == null) {
            return null;
        }
        return edgeFilterProducts.getGradientXY();
    }

    /**
     * @return the gradientX
     */
    public GreyscaleImage getGradientX() {
        if (edgeFilterProducts == null) {
            return null;
        }
        return edgeFilterProducts.getGradientX();
    }

    /**
     * @return the gradientY
     */
    public GreyscaleImage getGradientY() {
        if (edgeFilterProducts == null) {
            return null;
        }
        return edgeFilterProducts.getGradientY();
    }
    
    /**
     * @return the theta
     */
    public GreyscaleImage getTheta() {
        if (edgeFilterProducts == null) {
            return null;
        }
        return edgeFilterProducts.getTheta();
    }

    /**
     * @return the xOffset
     */
    public int getXOffset() {
        return trimmedXOffset;
    }

    /**
     * @return the yOffset
     */
    public int getYOffset() {
        return trimmedYOffset;
    }
    
    public CannyEdgeFilterSettings getCannyEdgeFilterSettings() {
                
        CannyEdgeFilterSettings settings = new CannyEdgeFilterSettings();
        
        if (useLineDrawingMode) {
            settings.setUseLineDrawingMode();
        }
        
        return settings;
    }
    
}
