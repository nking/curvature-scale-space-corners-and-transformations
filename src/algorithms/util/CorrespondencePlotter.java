package algorithms.util;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.misc.MiscMath;
import java.io.IOException;

/**
 * a helper class to visualize correspondence
 * lists.
 * An instance contains the transformations
 * for the two difference datasets to place
 * them in the output image reference frame.
 * 
 * @author nichole
 */
public class CorrespondencePlotter {
   
    private int spacerWidth = 50;
    
    private final int width1;
    private final int height1;
    private final int width2;
    private final int height2;
    
    private final int xOffset1;
    private final int yOffset1;
    private final int xOffset2;
    private final int yOffset2;
    
    private final Image comb;
    
    private int prevClr = -1;
    
    public CorrespondencePlotter(Image img1, Image img2) {
    
        width1 = img1.getWidth();
        height1 = img1.getHeight();
        
        width2 = img2.getWidth();
        height2 = img2.getHeight();
        
        xOffset1 = spacerWidth;
        yOffset1 = spacerWidth;
        
        xOffset2 = xOffset1 + width1 + spacerWidth;
        yOffset2 = spacerWidth;
        
        int maxX = xOffset2 + width2 + spacerWidth;
        int maxY = Math.max(yOffset1 + height1, 
            yOffset2 + height2) + spacerWidth;
    
        comb = new Image(maxX, maxY);
        
        for (int i = 0; i < img1.getWidth(); ++i) {
            for (int j = 0; j < img1.getHeight(); ++j) {
                int x = i + xOffset1;
                int y = j + yOffset1;
                int rgb = img1.getRGB(i, j);
                comb.setRGB(x, y, rgb);
            }
        }
        
        for (int i = 0; i < img2.getWidth(); ++i) {
            for (int j = 0; j < img2.getHeight(); ++j) {
                int x = i + xOffset2;
                int y = j + yOffset2;
                int rgb = img2.getRGB(i, j);
                comb.setRGB(x, y, rgb);
            }
        }
    }
    
    public CorrespondencePlotter(PairIntArray 
        boundary1, PairIntArray boundary2) {
        
        width1 = 1 + MiscMath.findMax(boundary1.x, boundary1.getN());
        height1 = 1 + MiscMath.findMax(boundary1.y, boundary1.getN());
        
        width2 = 1 + MiscMath.findMax(boundary2.x, boundary2.getN());
        height2 = 1 + MiscMath.findMax(boundary2.y, boundary2.getN());
        
        xOffset1 = spacerWidth;
        yOffset1 = spacerWidth;
        
        xOffset2 = xOffset1 + width1 + spacerWidth;
        yOffset2 = spacerWidth;
        
        int maxX = xOffset2 + width2 + spacerWidth;
        int maxY = Math.max(yOffset1 + height1, 
            yOffset2 + height2) + spacerWidth;
    
        comb = new Image(maxX, maxY);

        for (int i = 0; i < boundary1.getN(); ++i) {        
            int x = boundary1.getX(i) + xOffset1;
            int y = boundary1.getY(i) + yOffset1;
            comb.setRGB(x, y, 255, 255, 255);
        }
        
        for (int i = 0; i < boundary2.getN(); ++i) {        
            int x = boundary2.getX(i) + xOffset2;
            int y = boundary2.getY(i) + yOffset2;
            comb.setRGB(x, y, 255, 255, 255);
        }
    }
    
    public void drawLineInAlternatingColors(
        int x1, int y1, int x2, int y2,
        int nExtraForDot) {
        
        x1 += xOffset1;
        y1 += yOffset1;
        
        x2 += xOffset2;
        y2 += yOffset2;
        
        prevClr++;
        
        int[] clr = ImageIOHelper.getNextRGB(prevClr);
        
        ImageIOHelper.drawLineInImage(x1, y1, x2, y2, 
            comb, nExtraForDot, clr[0], clr[1], clr[2]);
    }
    
    public void drawLine(int x1, int y1, int x2, int y2,
        int r, int g, int b, int nExtraForDot) {
        
        x1 += xOffset1;
        y1 += yOffset1;
        
        x2 += xOffset2;
        y2 += yOffset2;
        
        prevClr++;
                
        ImageIOHelper.drawLineInImage(x1, y1, x2, y2, 
            comb, nExtraForDot, r, g, b);
    }
    
    public String writeImage(String fileSuffix) throws IOException {
        
         String dirPath = ResourceFinder.findDirectory("bin");
            
         String filePath = dirPath + "/img" + fileSuffix 
             + ".png";
            
         ImageIOHelper.writeOutputImage(filePath, comb);
         
         return filePath;
    }
   
}
