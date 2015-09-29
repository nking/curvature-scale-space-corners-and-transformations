package algorithms.imageProcessing;

import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.awt.Color;
import java.awt.color.ColorSpace;
import java.awt.image.BufferedImage;
import java.awt.image.ColorConvertOp;
import java.awt.image.Raster;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;
import javax.imageio.ImageIO;

/**
 *
 * @author nichole
 */
public class ImageIOHelper {
    
    public static Image readImage(String filePath) throws Exception {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            BufferedImage img = ImageIO.read(file);
            
            //System.out.println("imageType=" + img.getType());
            
            int h = img.getHeight();
            int w = img.getWidth();
            
            Image image = new Image(w, h);
            
            convertImage(img, image);
            
            return image;
            
        } catch (IOException e) {
        }
        
        return null;
    }
    
    public static ImageExt readImageExt(String filePath) throws Exception {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            BufferedImage bufferedInput = ImageIO.read(file);
            
            //System.out.println("imageType=" + img.getType());
            
            int h = bufferedInput.getHeight();
            int w = bufferedInput.getWidth();
            
            ImageExt image = new ImageExt(w, h);
            
            convertImage(bufferedInput, image);
            
            return image;
            
        } catch (IOException e) {
        }
        
        return null;
    }
    
    private static void convertImage(BufferedImage fromImage, Image toImage) 
        throws Exception {
     
        if (fromImage == null) {
            return;
        }
                
        int h = fromImage.getHeight();
        int w = fromImage.getWidth();

        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {

                //TYPE_INT_ARGB) and default sRGB colorspace
                int rgb = fromImage.getRGB(i, j);

                int r = (rgb >> 16) & 0xFF;
                int g = (rgb >> 8) & 0xFF;
                int b = rgb & 0xFF;

                toImage.setRGB(i, j, r, g, b);
            }
        }
    }
    
    public static ImageExt convertImage(GreyscaleImage fromImage) {
     
        if (fromImage == null) {
            return null;
        }
        
        int h = fromImage.getHeight();
        int w = fromImage.getWidth();

        ImageExt toImage = new ImageExt(w, h);
                
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                
                int v = fromImage.getValue(i, j);

                toImage.setRGB(i, j, v, v, v);
            }
        }
        
        return toImage;
    }
    
    /**
     * 
     * @param filePath
     * @return
     * @throws Exception 
     */
    public static GreyscaleImage readImageAsGrayScaleG(String filePath) 
        throws Exception {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            BufferedImage img = ImageIO.read(file);
            
            ColorSpace cs = ColorSpace.getInstance(ColorSpace.CS_GRAY);  
            ColorConvertOp op = new ColorConvertOp(cs, null);
            try {
                BufferedImage image2 = op.filter(img, null);
                img = image2;
            } catch (NullPointerException e) {
                // if type is custom, the source color space destination is not defined.
            }       
            int h = img.getHeight();
            int w = img.getWidth();
            
            GreyscaleImage image = new GreyscaleImage(w, h);
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                    
                    int rgb = img.getRGB(i, j);
                                        
                    int g = (rgb >> 8) & 0xFF;
                    
                    image.setValue(i, j, g);
                }
            }
            
            return image;
            
        } catch (IOException e) {
        }
        
        return null;
    }
    
    /**
     * read the Red channel of the image at filePath into a single channel
     * GreyscaleImage image.
     * 
     * @param filePath
     * @return
     * @throws Exception 
     */
    public static GreyscaleImage readImageAsGrayScaleR(String filePath) 
        throws Exception {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            BufferedImage img = ImageIO.read(file);
                       
            ColorSpace cs = ColorSpace.getInstance(ColorSpace.CS_GRAY);  
            ColorConvertOp op = new ColorConvertOp(cs, null);
            try {
                BufferedImage image2 = op.filter(img, null);
                img = image2;
            } catch (NullPointerException e) {
                // if type is custom, the source color space destination is not defined.
            }       
            int h = img.getHeight();
            int w = img.getWidth();
            
            GreyscaleImage image = new GreyscaleImage(w, h);
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                    
                    int rgb = img.getRGB(i, j);
                                        
                    int r = (rgb >> 16) & 0xFF;
                    
                    image.setValue(i, j, r);
                }
            }
            
            return image;
            
        } catch (IOException e) {
        }
        
        return null;
    }
    
    /**
     * read the Blue channel of the image at filePath into a single channel
     * GreyscaleImage image.
     * 
     * @param filePath
     * @return
     * @throws Exception 
     */
    public static GreyscaleImage readImageAsGrayScaleB(String filePath) 
        throws Exception {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            BufferedImage img = ImageIO.read(file);
                       
            ColorSpace cs = ColorSpace.getInstance(ColorSpace.CS_GRAY);  
            ColorConvertOp op = new ColorConvertOp(cs, null);
            try {
                BufferedImage image2 = op.filter(img, null);
                img = image2;
            } catch (NullPointerException e) {
                // if type is custom, the source color space destination is not defined.
            }       
            int h = img.getHeight();
            int w = img.getWidth();
            
            GreyscaleImage image = new GreyscaleImage(w, h);
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                    
                    int rgb = img.getRGB(i, j);
                                        
                    int b = rgb & 0xFF;
                    
                    image.setValue(i, j, b);
                }
            }
            
            return image;
            
        } catch (IOException e) {
            Logger.getLogger(ImageIOHelper.class.getName()).severe(e.getMessage());
        }
        
        return null;
    }
    
    /**
     * read image at filePath into a greyscale image.
     * 
     * @param filePath
     * @return
     * @throws Exception 
     */
    public static Image readImageAsGrayScale(String filePath) 
        throws Exception {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            BufferedImage img = ImageIO.read(file);
                        
            ColorSpace cs = ColorSpace.getInstance(ColorSpace.CS_GRAY);  
            ColorConvertOp op = new ColorConvertOp(cs, null);
            try {
                BufferedImage image2 = op.filter(img, null);
                img = image2;
            } catch (NullPointerException e) {
                // if type is custom, the source color space destination is not defined.
            }       
            int h = img.getHeight();
            int w = img.getWidth();
            
            Image image = new Image(w, h);
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                    
                    int rgb = img.getRGB(i, j);
                                        
                    int r = (rgb >> 16) & 0xFF;
                    int g = (rgb >> 8) & 0xFF;
                    int b = rgb & 0xFF;         
                    
                    image.setRGB(i, j, r, g, b);
                }
            }
            
            return image;
            
        } catch (IOException e) {
        }
        
        return null;
    }
    
    /**
     * read image at filePath into a greyscale image with assumptions that pixel
     * value '0' is '0' and anything else gets stored in output as a '1'.
     * 
     * @param filePath
     * @return
     * @throws Exception 
     */
    public static GreyscaleImage readImageAsBinary(String filePath) 
        throws Exception {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            BufferedImage img = ImageIO.read(file);
                 
            int h = img.getHeight();
            int w = img.getWidth();
           
            Raster raster = img.getData();
            
            GreyscaleImage image = new GreyscaleImage(w, h);
            
            int[] data = null;
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                    data = raster.getPixel(i, j, data);
                    int value = data[0];
                    image.setValue(i, j, value);
                }
            }
            
            return image;
            
        } catch (IOException e) {
        }
        
        return null;
    }
    
    /**
     * 
     * @param filePath
     * @return
     * @throws Exception 
     */
    public static GreyscaleImage readImageAsGrayScaleAvgRGB(String filePath) throws Exception {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            BufferedImage img = ImageIO.read(file);
            
            ColorSpace cs = ColorSpace.getInstance(ColorSpace.CS_GRAY);  
            ColorConvertOp op = new ColorConvertOp(cs, null);  
            try {
                BufferedImage image2 = op.filter(img, null);
                img = image2;
            } catch (NullPointerException e) {
                // if type is custom, the source color space destination is not defined.
            }
                        
            int h = img.getHeight();
            int w = img.getWidth();
            
            GreyscaleImage image = new GreyscaleImage(w, h);
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                    
                    int rgb = img.getRGB(i, j);
                                        
                    int r = (rgb >> 16) & 0xFF;
                    int g = (rgb >> 8) & 0xFF;
                    int b = rgb & 0xFF;
                    float v = (r + g + b)/3.f;
                    
                    image.setValue(i, j, (int)v);
                }
            }
            
            return image;
            
        } catch (IOException e) {
        }
        
        return null;
    }
    
    public static void writeOutputImage(String filePath, Image data) 
        throws IOException {
        
        writeOutputImage(filePath, data,  BufferedImage.TYPE_INT_RGB);       
    }
    
    public static void writeOutputImage(String filePath, Image data, 
        final int imageType) throws IOException {
                
        if (data == null) {
            return;
        }
        if (filePath == null) {
            throw new IllegalArgumentException("filePath cannot be null");
        }
        
        int w = data.getWidth();
        
        int h = data.getHeight();
        
        BufferedImage outputImage = new BufferedImage(w, h, 
            BufferedImage.TYPE_INT_RGB);
        
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                int rgbValue = data.getRGB(i, j);
                outputImage.setRGB(i, j, rgbValue);
            }
        }
        
        ImageIO.write(outputImage, "PNG", new File(filePath));
    }
    
    public static void writeOutputImage(String filePath, GreyscaleImage data) 
        throws IOException {
                
        if (data == null) {
            return;
        }
        if (filePath == null) {
            throw new IllegalArgumentException("filePath cannot be null");
        }
        
        int w = data.getWidth();
        
        int h = data.getHeight();
        
        BufferedImage outputImage = new BufferedImage(w, h, 
            BufferedImage.TYPE_BYTE_GRAY);
        
        WritableRaster raster = outputImage.getRaster();
        
        for (int i = 0; i < w; i++) {
            
            for (int j = 0; j < h; j++) {
                
                int value = data.getValue(i, j);
                
                raster.setSample(i, j, 0, value);
            }
        }
        
        ImageIO.write(outputImage, "PNG", new File(filePath));
    }
    
    /**
     * draw the edges over input and mark the corners.  The edges are drawn
     * in rotating colors to make it easier to see closed curves.
     * The corners are drawn in red dots of size 2. For best results, the
     * input image should be greyscale.
     * @param input
     * @param edges
     * @param corners
     * @param fileName
     * @param width
     * @param height
     * @throws IOException 
     */
    public static void writeImage(Image input, PairIntArray[] edges, 
        PairIntArray corners, String fileName, int width, int height) throws 
        IOException {
                
        if (edges == null || input == null || corners == null) {
            return;
        }
        if (fileName == null) {
            throw new IllegalArgumentException("fileName cannot be null");
        }
        
        addAlternatingColorCurvesToImage(edges, fileName, false, input);
                
        addCurveToImage(corners, input, 1, 255, 0, 0);
        
        String dirPath = ResourceFinder.findDirectory("bin");
       
        ImageIOHelper.writeOutputImage(dirPath + "/" + fileName, input);
    }
    
    /**
     * draw the edge over the image using the given rgb colors and the size
     * of the dot beyond 1 pixel.
     * 
     * @param edge
     * @param input
     * @param nExtraForDot
     * @param rClr
     * @param gClr
     * @param bClr 
     */
    public static void addCurveToImage(PairIntArray edge, Image input, 
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        if (edge == null || input == null) {
            return;
        }
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        for (int i = 0; i < edge.getN(); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                
                int xx = x + dx;
                
                if ((xx < 0) || (xx > (w - 1))) {
                    continue;
                }
                for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                    int yy = y + dy;
                    if ((yy < 0) || (yy > (h - 1))) {
                        continue;
                    }
                    input.setRGB(xx, yy, rClr, gClr, bClr);
                }
            }
        }
    }
    
    /**
     * draw the edge over the image using the given rgb colors and the size
     * of the dot beyond 1 pixel.
     * 
     * @param points
     * @param input
     * @param nExtraForDot
     * @param rClr
     * @param gClr
     * @param bClr 
     */
    public static void addCurveToImage(Collection<PairInt> points, Image input, 
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        if (points == null || input == null) {
            return;
        }
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                
                int xx = x + dx;
                
                if ((xx < 0) || (xx > (w - 1))) {
                    continue;
                }
                for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                    int yy = y + dy;
                    if ((yy < 0) || (yy > (h - 1))) {
                        continue;
                    }
                    input.setRGB(xx, yy, rClr, gClr, bClr);
                }
            }
        }
    }
    
    /**
     * draw lines from points at index i to index i+1 in the requested color
     * onto the input.
     * 
     * @param xVertexes
     * @param yVertexes
     * @param input
     * @param nExtraForDot
     * @param rClr
     * @param gClr
     * @param bClr 
     */
    public static void drawLinesInImage(int[] xVertexes, int[] yVertexes, 
        Image input, int nExtraForDot, int rClr, int gClr, int bClr) {
        
        if (xVertexes == null || yVertexes == null || input == null) {
            return;
        }
        if (xVertexes.length != yVertexes.length) {
            return;
        }
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        for (int i = 0; i < (xVertexes.length - 1); i++) {
            
            int x1 = xVertexes[i];
            int y1 = yVertexes[i];
            
            int x2 = xVertexes[i + 1];
            int y2 = yVertexes[i + 1];
            
            if (x2 < x1) {
                int swap = x1;
                x1 = x2;
                x2 = swap;
                swap = y1;
                y1 = y2;
                y2 = swap;
            }
            
            int dx0 = x2 - x1;
            int dy0 = y2 - y1;
            int nLine = (int)Math.ceil(Math.sqrt(dx0*dx0 + dy0*dy0));
            if (nLine == 1 && (dx0 > 0 || dy0 > 0)) {
                nLine = 2;
            }
int z = 1;          
            for (int ii = 0; ii < nLine; ++ii) {
            
                int x, y;
                
                if (dx0 == 0) {
                    x = x1;
                    y = (y1 + Math.round(dy0*((float)ii/(float)(nLine - 1))));
                } else if (dy0 == 0) {
                    x = (x1 + Math.round(dx0*((float)ii/(float)(nLine - 1))));
                    y = y1;
                } else {
                    x = (x1 + Math.round(dx0*((float)ii/(float)(nLine - 1))));
                    y = (y1 + Math.round(dy0*((float)ii/(float)(nLine - 1))));
                }                
int z1 = 1;            
                for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {

                    int xx = x + dx;
                    
                    if ((xx < 0) || (xx > (w - 1))) {
                        continue;
                    }
                    for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                        int yy = y + dy;
                        if ((yy < 0) || (yy > (h - 1))) {
                            continue;
                        }
                        input.setRGB(xx, yy, rClr, gClr, bClr);
                    }
                }
            }
        }
    }
    
    /**
     * draw lines from points at index i to index i+1 in the requested color
     * onto the input.
     * 
     * @param xVertexes
     * @param yVertexes
     * @param input
     * @param nExtraForDot
     * @param rClr
     * @param gClr
     * @param bClr 
     */
    public static void addCurveToImage(int[] xVertexes, int[] yVertexes, 
        Image input, int nExtraForDot, int rClr, int gClr, int bClr) {
        
        if (xVertexes == null || yVertexes == null || input == null) {
            return;
        }
        if (xVertexes.length != yVertexes.length) {
            return;
        }
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        for (int i = 0; i < xVertexes.length; i++) {
            
            int x = xVertexes[i];
            int y = yVertexes[i];
            
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {

                int xx = x + dx;

                if ((xx < 0) || (xx > (w - 1))) {
                    continue;
                }
                for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                    int yy = y + dy;
                    if ((yy < 0) || (yy > (h - 1))) {
                        continue;
                    }
                    input.setRGB(xx, yy, rClr, gClr, bClr);
                }
            }
        }
    }
    
    /**
     * draw the edge over the image using the given rgb colors and the size
     * of the dot beyond 1 pixel.
     * 
     * @param edge
     * @param input
     * @param nExtraForDot
     * @param rClr
     * @param gClr
     * @param bClr 
     */
    public static void addCurveToImage(PairFloatArray edge, Image input, 
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        if (edge == null || input == null) {
            return;
        }
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        for (int i = 0; i < edge.getN(); i++) {
            int x = Math.round(edge.getX(i));
            int y = Math.round(edge.getY(i));
            
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                
                int xx = x + dx;
                
                if ((xx < 0) || (xx > (w - 1))) {
                    continue;
                }
                for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                    int yy = y + dy;
                    if ((yy < 0) || (yy > (h - 1))) {
                        continue;
                    }
                    input.setRGB(xx, yy, rClr, gClr, bClr);
                }
            }
        }
    }
    
    /**
     * draw the point over the image using the given rgb colors and the size
     * of the dot beyond 1 pixel.
     * 
     * @param x of point
     * @param y of point
     * @param input
     * @param nExtraForDot
     * @param rClr
     * @param gClr
     * @param bClr 
     */
    public static void addPointToImage(int x, int y, Image input, 
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {

            int xx = x + dx;

            if ((xx < 0) || (xx > (w - 1))) {
                continue;
            }
            for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                int yy = y + dy;
                if ((yy < 0) || (yy > (h - 1))) {
                    continue;
                }
                input.setRGB(xx, yy, rClr, gClr, bClr);
            }
        }
    }
    /**
     * draw the point over the image using the given rgb colors and the size
     * of the dot beyond 1 pixel.
     * 
     * @param x of point
     * @param y of point
     * @param input
     * @param nExtraForDot
     * @param rgbClr
     */
    public static void addPointToImage(int x, int y, Image input, 
        int nExtraForDot, int rgbClr) {
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {

            int xx = x + dx;

            if ((xx < 0) || (xx > (w - 1))) {
                continue;
            }
            for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                int yy = y + dy;
                if ((yy < 0) || (yy > (h - 1))) {
                    continue;
                }
                input.setRGB(xx, yy, rgbClr);
            }
        }
    }
    
    /**
     * draw the point over the image using the given rgb colors and the size
     * of the dot beyond 1 pixel.
     * 
     * @param x of point
     * @param y of point
     * @param input
     * @param nExtraForDot
     * @param rClr
     * @param gClr
     * @param bClr 
     */
    public static void addPointToImage(float x, float y, Image input, 
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {

            int xx = Math.round(x + dx);

            if ((xx < 0) || (xx > (w - 1))) {
                continue;
            }
            for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                int yy = Math.round(y + dy);
                if ((yy < 0) || (yy > (h - 1))) {
                    continue;
                }
                input.setRGB((int)xx, (int)yy, rClr, gClr, bClr);
            }
        }
    }
    
    /**
     * draw the edge over the image using the given rgb colors and the size
     * of the dot beyond 1 pixel.
     * 
     * @param edgeX x points of pair (edgeX, edgeY)
     * @param edgeY y points of pair (edgeX, edgeY)
     * @param nPointsInEdge number of points in edgeX 
     * @param input
     * @param nExtraForDot
     * @param rClr
     * @param gClr
     * @param bClr 
     */
    public static void addCurveToImage(int[] edgeX, int[] edgeY, 
        final int nPointsInEdge, Image input, int nExtraForDot, int rClr, 
        int gClr, int bClr) {
        
        if (edgeX == null || edgeY == null) {
            return;
        }
        if (edgeX.length != edgeY.length) {
            throw new IllegalArgumentException(
                "edgeX and edgeY should be same length");
        }
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        for (int i = 0; i < nPointsInEdge; i++) {
            
            int x = edgeX[i];
            int y = edgeY[i];
            
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                
                int xx = x + dx;
                
                if ((xx < 0) || (xx > (w - 1))) {
                    continue;
                }
                for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                    int yy = y + dy;
                    if ((yy < 0) || (yy > (h - 1))) {
                        continue;
                    }
                    input.setRGB(xx, yy, rClr, gClr, bClr);
                }
            }
        }
    }
    
    public static void addAlternatingColorCurvesToImage(
        PairIntArray[] curves, String fileName, boolean writeImage, 
        Image input) throws IOException {
        
        if (curves == null || input == null) {
            return;
        }
        if (fileName == null) {
            throw new IllegalArgumentException("fileName cannot be null");
        }
        
        List<PairIntArray> c = new ArrayList<PairIntArray>();
        
        for (PairIntArray a : curves) {
            c.add(a);
        }
        
        addAlternatingColorCurvesToImage(c, fileName, writeImage, input);
    }
    
    public static void addAlternatingColorCurvesToImage(
        List<PairIntArray> curves, Image input) throws IOException {
        
        if (curves == null || input == null) {
            return;
        }
        
        addAlternatingColorCurvesToImage(
            curves.toArray(new PairIntArray[curves.size()]),
            input);
    }
    
    public static int getNextColorRGB(int count) {
        
        int[] c = getNextRGB(count);
        
        int rgb = (((c[0] & 0x0ff) << 16) | ((c[1] & 0x0ff) << 8) | 
            (c[2] & 0x0ff));
        
        return rgb;
    }
    
    /*
    http://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
    {255,0,0}, {228,228,0}, {0,255,0}, 
    {0,255,255}, {176,176,255}, {255,0,255}, 
    {228,228,228}, {176,0,0}, {186,186,0}, 
    {0,176,0}, {0,176,176}, {132,132,255}, 
    {176,0,176}, {186,186,186}, {135,0,0}, 
    {135,135,0}, {0,135,0}, {0,135,135}, 
    {73,73,255}, {135,0,135}, {135,135,135}, 
    {85,0,0}, {84,84,0}, {0,85,0}, 
    {0,85,85}, {0,0,255}, {85,0,85}, {84,84,84}
    */
    public static int[] getNextRGB(int count) {
        
        if (count < 0) {
            count = 0;
        } else {
            count = count % 28;
        }
        
        int[][] c = new int[28][];
        c[0] = new int[]{255,0,0};
        c[1] = new int[]{228,228,0};
        c[2] = new int[]{0,255,0};
        c[3] = new int[]{0,255,255}; 
        c[4] = new int[]{176,176,255};
        c[5] = new int[]{255,0,255};
        c[6] = new int[]{228,228,228};
        c[7] = new int[]{176,0,0}; 
        c[8] = new int[]{186,186,0};
        c[9] = new int[]{0,176,0}; 
        c[10] = new int[]{0,176,176}; 
        c[11] = new int[]{132,132,255}; 
        c[12] = new int[]{176,0,176}; 
        c[13] = new int[]{186,186,186}; 
        c[14] = new int[]{135,0,0};
        c[15] = new int[]{135,135,0}; 
        c[16] = new int[]{0,135,0}; 
        c[17] = new int[]{0,135,135}; 
        c[18] = new int[]{73,73,255}; 
        c[19] = new int[]{135,0,135}; 
        c[20] = new int[]{135,135,135}; 
        c[21] = new int[]{85,0,0};
        c[22] = new int[]{84,84,0}; 
        c[23] = new int[]{0,85,0};
        c[24] = new int[]{0,85,85};
        c[25] = new int[]{0,0,255};
        c[26] = new int[]{85,0,85};
        c[27] = new int[]{84,84,84};
        
        return c[count];
    }
    
    public static void addAlternatingColorCurvesToImage(
        PairIntArray[] curves, Image input) throws IOException {
        
        if (curves == null || input == null) {
            return;
        }
        
        int clr = 0;
        
        for (int i = 0; i < curves.length; i++) {
            
            PairIntArray edge = curves[i];
            
            int c = getNextColorRGB(clr);
            
            for (int ii = 0; ii < edge.getN(); ii++) {
                
                int col = edge.getX(ii);
                int row = edge.getY(ii);
                
                if ((col > -1) && (col < input.getWidth()) &&
                    (row > -1) && (row < input.getHeight())) {
                
                    input.setRGB(col, row, c);
                }
            }
            clr++;
        }
    }
    
    public static void addAlternatingColorPointSetsToImage(
        List<Set<PairInt>> pointSets, int xOffset, int yOffset, Image input) 
        throws IOException {
        
        if (pointSets == null || input == null) {
            return;
        }
        
        int clr = 0;
        
        for (int i = 0; i < pointSets.size(); i++) {
            
            Set<PairInt> points = pointSets.get(i);
            
            int c = getNextColorRGB(clr);
            
            for (PairInt p : points) {
                
                int col = p.getX() + xOffset;
                int row = p.getY() + yOffset;
                
                if ((col > -1) && (col < input.getWidth()) &&
                    (row > -1) && (row < input.getHeight())) {
                
                    input.setRGB(col, row, c);
                }
            }
            clr++;
        }
    }
    
    public static void addToImage(Set<PairInt> points, int xOffsetToApply, 
        int yOffsetToApply, Image input) throws IOException {
        
        if (points == null || input == null) {
            return;
        }
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        int c = Color.BLUE.getRGB();
          
        for (PairInt p : points) {

            int col = p.getX() + xOffsetToApply;
            int row = p.getY() + yOffsetToApply;

            if ((col < 0) || (row < 0) || (col > (w - 1)) || (row > (h - 1))) {
                continue;
            }

            input.setRGB(col, row, c);
        }
    }
    
    public static void addToImage(Set<PairInt> points, int xOffsetToApply, 
        int yOffsetToApply, Image input, int rClr, int gClr, int bClr) throws 
        IOException {
        
        if (points == null || input == null) {
            return;
        }
        
        int w = input.getWidth();
        int h = input.getHeight();
                  
        for (PairInt p : points) {

            int col = p.getX() + xOffsetToApply;
            int row = p.getY() + yOffsetToApply;
            
            if ((col < 0) || (row < 0) || (col > (w - 1)) || (row > (h - 1))) {
                continue;
            }

            input.setRGB(col, row, rClr, gClr, bClr);
        }
    }
    
    public static void addToImage(Set<PairInt> points, int xOffsetToApply, 
        int yOffsetToApply, Image input, int nExtraForDot,
        int rClr, int gClr, int bClr) throws 
        IOException {
        
        if (points == null || input == null) {
            return;
        }
        
        int w = input.getWidth();
        int h = input.getHeight();
                  
        for (PairInt p : points) {

            int x = p.getX() + xOffsetToApply;
            int y = p.getY() + yOffsetToApply;
            
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); ++dx) {
                int xx = x + dx;
                if ((xx < 0) || (xx > (w - 1))) {
                    continue;
                }
                for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                    int yy = y + dy;
                    if ((yy < 0) || (yy > (h - 1))) {
                        continue;
                    }
                    input.setRGB(xx, yy, rClr, gClr, bClr);
                }
            }
        }
    }

    public static void addAlternatingColorCurvesToImage(
        List<PairIntArray> curves, String fileName, boolean writeImage, 
        Image input) throws IOException {
        
        if (curves == null || input == null) {
            return;
        }
        if (fileName == null) {
            throw new IllegalArgumentException("fileName cannot be null");
        }
        
        addAlternatingColorCurvesToImage(curves, input);
       
        if (!writeImage) {
            return;
        }
        
        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/" + fileName, input);
        } catch (IOException ex) {
            Logger.getLogger(ImageIOHelper.class.getName()).severe(
                "ERROR: l446" + ex.getMessage());
        }
    }

    public static void addCurveToImage(ScaleSpaceCurve curve, Image input, 
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        if (curve == null || input == null) {
            return;
        }
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        for (int i = 0; i < curve.getSize(); i++) {
            int x = Math.round(curve.getX(i));
            int y = Math.round(curve.getY(i));
            
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); ++dx) {
                int xx = x + dx;
                if ((xx < 0) || (xx > (w - 1))) {
                    continue;
                }
                for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                    int yy = y + dy;
                    if ((yy < 0) || (yy > (h - 1))) {
                        continue;
                    }
                    input.setRGB(xx, yy, rClr, gClr, bClr);
                }
            }
        }
    }

    public static void addToImage(float[] xP, float[] yP, 
        int xOffsetToApply, int yOffsetToApply, Image input,
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        if (xP == null || yP == null || input == null) {
            return;
        }
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        for (int i = 0; i < xP.length; i++) {
            
            int x = Math.round(xP[i]) + xOffsetToApply;
            
            int y = Math.round(yP[i]) + yOffsetToApply;
            
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); ++dx) {
                int xx = x + dx;
                if ((xx < 0) || (xx > (w - 1))) {
                    continue;
                }
                for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                    int yy = y + dy;
                    if ((yy < 0) || (yy > (h - 1))) {
                        continue;
                    }
                    input.setRGB(xx, yy, rClr, gClr, bClr);
                }
            }
        }
    }
    
    public static void addToImage(int[] xP, int[] yP, 
        int xOffsetToApply, int yOffsetToApply, Image input,
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        if (xP == null || yP == null || input == null) {
            return;
        }
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        for (int i = 0; i < xP.length; i++) {
            
            int x = xP[i] + xOffsetToApply;
            
            int y = yP[i] + yOffsetToApply;
            
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                int xx = x + dx;
                if ((xx < 0) || (xx > (w - 1))) {
                    continue;
                }
                for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                    int yy = y + dy;
                    if ((yy < 0) || (yy > (h - 1))) {
                        continue;
                    }
                    input.setRGB(xx, yy, rClr, gClr, bClr);
                }
            }
        }
    }
    
    /**
     * given data, apply a factor and offset to put it in range 0 to 255 and
     * return an image from it.
     * @param data
     * @return 
     */
    public static GreyscaleImage scaleToImgRange(double[][] data) {
        
        double maxValue = Double.MIN_VALUE;
        double minValue = Double.MAX_VALUE;
        for (int i = 0; i < data.length; ++i) {
            for (int j = 0; j < data[i].length; ++j) {
                double v = data[i][j];
                if (v < minValue) {
                    minValue = v;
                }
                if (v > maxValue) {
                    maxValue = v;
                }
            }
        }
        
        /* scale data between 0 and 255
               
        (max - min) * factor = 255 --> factor = 255/(max - min)
        then subtract min, min*factor
        for example 10  20, delta=10, factor = 25.5
        (10*25.5)-(10*25.5) = 0
        (20*25.5)-(10*25.5) = 255
        */
        double factor = 255./(maxValue - minValue);
        double zp = minValue * factor;
        
        GreyscaleImage img = new GreyscaleImage(data.length, data[0].length);
        
        for (int i = 0; i < data.length; ++i) {
            for (int j = 0; j < data[i].length; ++j) {
                double v = data[i][j];
                int vs = (int)Math.round((v * factor) - zp);
                img.setValue(i, j, vs);
            }
        }
        
        return img;
    }
    
    /**
     * given data, apply a factor and offset to put it in range 0 to 255 and
     * return an image from it.
     * @param data
     * @return 
     */
    public static GreyscaleImage scaleToImgRange(int[][] data) {
        
        int maxValue = Integer.MIN_VALUE;
        int minValue = Integer.MAX_VALUE;
        for (int i = 0; i < data.length; ++i) {
            for (int j = 0; j < data[i].length; ++j) {
                int v = data[i][j];
                if (v < minValue) {
                    minValue = v;
                }
                if (v > maxValue) {
                    maxValue = v;
                }
            }
        }
        
        /* scale data between 0 and 255
               
        (max - min) * factor = 255 --> factor = 255/(max - min)
        then subtract min, min*factor
        for example 10  20, delta=10, factor = 25.5
        (10*25.5)-(10*25.5) = 0
        (20*25.5)-(10*25.5) = 255
        */
        double factor = 255./(double)(maxValue - minValue);
        double zp = minValue * factor;
        
        GreyscaleImage img = new GreyscaleImage(data.length, data[0].length);
        
        for (int i = 0; i < data.length; ++i) {
            for (int j = 0; j < data[i].length; ++j) {
                int v = data[i][j];
                int vs = (int)Math.round((v * factor) - zp);
                img.setValue(i, j, vs);
            }
        }
        
        return img;
    }
}
