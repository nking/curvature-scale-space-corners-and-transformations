package algorithms.imageProcessing;

import algorithms.compGeometry.PointInPolygon;
import algorithms.util.ArrayPair;
import java.awt.Color;

/**
 *
 * @author nichole
 */
public class SunColors {

    private CIEChromaticity cieC = new CIEChromaticity();
        
    private PointInPolygon pInPoly = new PointInPolygon();
        
    private ArrayPair yellowBounds = cieC.getYellowPolynomial();
        
    private float totRGBLimit = 0.2f * 768.f;
    
    private boolean useDarkSkiesLogic = false;
   
    public boolean isSunCenterColor(int r, int g, int b) {
        
        if ((r < 15) && (g < 15) && (b < 15)) {
            return false;
        }

        float[] hsb = new float[3];
        Color.RGBtoHSB(r, g, b, hsb);

        if (useDarkSkiesLogic) {
            
            if ((hsb[1] > 0.1f) && (hsb[2] > 0.25)) {

                float[] pixCIEXY = cieC.rgbToXYChromaticity(r, g, b);
                
                if (pInPoly.isInSimpleCurve(pixCIEXY[0], pixCIEXY[1],
                    yellowBounds.getX(), yellowBounds.getY(),
                    yellowBounds.getX().length)) {

                    return true;
                }
            }
            
        } else {
            
            if ((r >= 240) && (hsb[2] >= 0.87) && (hsb[1] < .30)
                && ((r + g + b) > totRGBLimit)) {

                return true;
            }
        }
        
        return false;
    }
     
    public void useDarkSkiesLogic() {
        
        useDarkSkiesLogic = true;
    }
    
}
