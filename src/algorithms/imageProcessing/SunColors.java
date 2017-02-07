package algorithms.imageProcessing;

import java.awt.Color;

/**
 *
 * @author nichole
 */
public class SunColors {
    
    /*
    photosphere of sun, 5800K,
    
    average digital camera response from detector + filters...
    
    range of atmos conditions
        sunset/dawn
           - clear skies
           - with clouds
       mid day
          - clear skies
          - with clouds
    
    with sunset/dawn and no clouds, could write a simple radiative transfer
    model w/ source being sun, optically thin (tau < 1) atmos having air 
    and aerosols and no clouds
    then blue light scattered out of beam by an amount dependent upon
    the airmass which is dependent upon the altitude.
    
    but for the other cases, the conditions are not single scattering,
    and clouds may be optically thin or optically thick and in between
    contributing a range of absorption and sometimes reflection.
    
    */
       
    public boolean isSunCenterColor(int r, int g, int b) {
        
        if ((r < 15) && (g < 15) && (b < 15)) {
            return false;
        }

        float[] hsb = new float[3];
        Color.RGBtoHSB(r, g, b, hsb);

        return isSunCenterColor(hsb[0], hsb[1], hsb[2]);
    }
    
    public boolean isSunCenterColor(ImageExt colorImg, int pixelCol, int pixelRow) {
        
        int idx = colorImg.getInternalIndex(pixelCol, pixelRow);
        
        return isSunCenterColor(colorImg, idx);
    }
    
    public boolean isSunCenterColor(float h, float s, float v) {
        if (
            (s < 0.4f) && 
            (v > 0.25)
            && (h >= 0.0) && (h <= 0.18)) {
            return true;
        } else if ((s < 0.05) && (v > 0.95)) {
            return true;
        }
        return false;
    }
    
    public boolean isSunCenterColor(ImageExt colorImg, int pixelIndex) {
         
        if ((colorImg.getR(pixelIndex) < 15) && 
            (colorImg.getG(pixelIndex) < 15) && 
            (colorImg.getB(pixelIndex) < 15)) {
            return false;
        }
        
        float saturation = colorImg.getSaturation(pixelIndex);
        float brightness = colorImg.getBrightness(pixelIndex);
        float hue = colorImg.getHue(pixelIndex);
        
        boolean a = isSunCenterColor(hue, saturation, brightness);
        //System.out.format("(%d,%d) %.3f, %.3f, %.3f sunClr=%b\n", colorImg.getCol(pixelIndex), 
        //    colorImg.getRow(pixelIndex), hue, saturation, brightness, a);
        
        return a;
    }
     
}
