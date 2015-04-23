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
    protected final int red;
    protected final int green;
    protected final int blue;
    protected final double rgbTot;
    
    public ColorData(boolean skyIsRed,
        double pixContrast, double pixBlueOrRedDiff, 
        double pixCIEXDiff, double pixCIEYDiff, 
        double skyStDevContrast, double skyStDevBlueOrRedDiff,
        double skyStDevCIEX, double skyStDevCIEY, int r, int g, int b) {
        
        this.skyIsRed = skyIsRed;
        this.contrast = pixContrast;
        this.colorDiff = pixBlueOrRedDiff;
        this.skyStDevColorDiff = skyStDevBlueOrRedDiff;
        this.skyStDevContrast = skyStDevContrast;
        this.diffCIEX = pixCIEXDiff;
        this.diffCIEY = pixCIEYDiff;
        this.skyStdDevCIEX = skyStDevCIEX;
        this.skyStdDevCIEY = skyStDevCIEY;
        this.red = r;
        this.green = g;
        this.blue = b;
        this.rgbTot = (r + g + b);
    }
    
    public double getParameter(PARAM parameter) {
        
        double param;
        
        switch (parameter) {
            case CONTRAST:
                param = contrast;
                break;
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
            case RED:
                param = red;
                break;
            case GREEN:
                param = green;
                break;
            case BLUE:
                param = blue;
                break;
            case R_DIV_TOT:
                param = (double)red/rgbTot;
                break;
            case G_DIV_TOT:
                param = (double)green/rgbTot;
                break;
            case B_DIV_TOT:
                param = (double)blue/rgbTot;
                break;
            case DIFF_R_DIV_TOT_ONE_THIRD:
                param = Math.abs((1./3.) - ((double)red/rgbTot));
                break;
            case DIFF_G_DIV_TOT_ONE_THIRD:
                param = Math.abs((1./3.) - ((double)green/rgbTot));
                break;
            case DIFF_B_DIV_TOT_ONE_THIRD:
                param = Math.abs((1./3.) - ((double)blue/rgbTot));
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
