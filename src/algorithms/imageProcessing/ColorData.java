package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class ColorData {

    protected final boolean skyIsRed;
    protected final double contrast;
    protected final double colorDiff;
    protected final double skyStDevContrast;
    protected final double skyStDevColorDiff;
    protected final double diffCIEX;
    protected final double diffCIEY;
    protected final double skyStdDevCIEX;
    protected final double skyStdDevCIEY;
    
    public ColorData(boolean skyIsRed,
        double pixContrast, double pixBlueOrRedDiff, 
        double pixCIEXDiff, double pixCIEYDiff, 
        double skyStDevContrast, double skyStDevBlueOrRedDiff,
        double skyStDevCIEX, double skyStDevCIEY) {
        
        this.skyIsRed = skyIsRed;
        this.contrast = pixContrast;
        this.colorDiff = pixBlueOrRedDiff;
        this.skyStDevColorDiff = skyStDevBlueOrRedDiff;
        this.skyStDevContrast = skyStDevContrast;
        this.diffCIEX = pixCIEXDiff;
        this.diffCIEY = pixCIEYDiff;
        this.skyStdDevCIEX = skyStDevCIEX;
        this.skyStdDevCIEY = skyStDevCIEY;
    }
    
    public double getParameter(PARAM parameter) {
        
        double param;
        
        switch (parameter) {
            case ABSOLUTE_CONTRAST:
                param = Math.abs(contrast);
                break;
            case ABSOLUTE_DIFF_BLUE_OR_RED:
                param = Math.abs(colorDiff);
                break;
            case STDEV_CONTRAST:
                param = skyStDevContrast;
                break;
            case STDEV_BLUE_OR_RED:
                param = skyStDevColorDiff;
                break;
            case DIFF_CIEX:
                param = diffCIEX;
                break;
            case DIFF_CIEY:
                param = diffCIEY;
                break;
            case STDEV_CIEX:
                param = skyStdDevCIEX;
                break;
            case STDEV_CIEY:
                param = skyStdDevCIEY;
                break;
            case INT_ONE:
                param = 1.0;
                break;
            default:
                param = 1.0;
        }
        
        return param;
    }
    
    public boolean skyIsRed() {
        return skyIsRed;
    }
}
