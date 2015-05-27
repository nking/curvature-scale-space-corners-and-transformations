package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class ImageHelperForTests {
    
    private GreyscaleImage img;
    
    private GreyscaleImage gradientXY = null;
    
    private GreyscaleImage gradientX = null;
    
    private GreyscaleImage gradientY = null;
    
    private GreyscaleImage theta = null;
    
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
        
        CannyEdgeFilter filter = new CannyEdgeFilter();
        
        CannyEdgeFilterSettings settings = getCannyEdgeFilterSettings();
        
        filter.setSetters(settings);
                
        filter.applyFilter(img);
        
        gradientXY = filter.getGradientXY();
        
        gradientX = filter.getGradientX();
        
        gradientY = filter.getGradientY();
        
        theta = filter.getTheta();        
    }
    
    /**
     * @return the gradientXY
     */
    public GreyscaleImage getGradientXY() {
        return gradientXY;
    }

    /**
     * @return the gradientX
     */
    public GreyscaleImage getGradientX() {
        return gradientX;
    }

    /**
     * @return the gradientY
     */
    public GreyscaleImage getGradientY() {
        return gradientY;
    }

    /**
     * @return the theta
     */
    public GreyscaleImage getTheta() {
        return theta;
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
        
        if (useOutdoorMode) {
            settings.setUseOutdoorMode();
        }
        
        if (useLineDrawingMode) {
            settings.setUseLineDrawingMode();
        }
        
        return settings;
    }
    
}
