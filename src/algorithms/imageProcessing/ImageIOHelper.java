package algorithms.imageProcessing;

import algorithms.ResourceFinder;
import java.awt.Color;
import java.awt.color.ColorSpace;
import java.awt.image.BufferedImage;
import java.awt.image.ColorConvertOp;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
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
            
            int type = img.getType();
            
            /*
            System.out.println("BufferedImage.TYPE_INT_RGB=" + BufferedImage.TYPE_INT_RGB);
            System.out.println("BufferedImage.TYPE_INT_ARGB=" + BufferedImage.TYPE_INT_ARGB);
            System.out.println("BufferedImage.TYPE_INT_ARGB_PRE=" + BufferedImage.TYPE_INT_ARGB_PRE);
            System.out.println("BufferedImage.TYPE_BYTE_GRAY=" + BufferedImage.TYPE_BYTE_GRAY);
            System.out.println("BufferedImage.TYPE_BYTE_BINARY=" + BufferedImage.TYPE_BYTE_BINARY);
            System.out.println("BufferedImage.TYPE_BYTE_INDEXED=" + BufferedImage.TYPE_BYTE_INDEXED);
            System.out.println("BufferedImage.TYPE_USHORT_GRAY=" + BufferedImage.TYPE_USHORT_GRAY);
            System.out.println("BufferedImage.TYPE_CUSTOM=" + BufferedImage.TYPE_CUSTOM);
            System.out.println("IMAGE TYPE=" + type);
            */
            /*
            BufferedImage.TYPE_INT_RGB, BufferedImage.TYPE_INT_ARGB, 
            BufferedImage.TYPE_INT_ARGB_PRE, BufferedImage.TYPE_INT_BGR, 
            BufferedImage.TYPE_3BYTE_BGR, BufferedImage.TYPE_4BYTE_ABGR, 
            BufferedImage.TYPE_4BYTE_ABGR_PRE, 
            BufferedImage.TYPE_BYTE_GRAY, BufferedImage.TYPE_BYTE_BINARY, 
            BufferedImage.TYPE_BYTE_INDEXED, BufferedImage.TYPE_USHORT_GRAY, 
            BufferedImage.TYPE_USHORT_565_RGB, 
            BufferedImage.TYPE_USHORT_555_RGB, BufferedImage.TYPE_CUSTOM
            */
            
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
            
            int type = img.getType();
            
            /*
            BufferedImage.TYPE_INT_RGB, BufferedImage.TYPE_INT_ARGB, BufferedImage.TYPE_INT_ARGB_PRE, BufferedImage.TYPE_INT_BGR, BufferedImage.TYPE_3BYTE_BGR, BufferedImage.TYPE_4BYTE_ABGR, BufferedImage.TYPE_4BYTE_ABGR_PRE, 
            BufferedImage.TYPE_BYTE_GRAY, BufferedImage.TYPE_BYTE_BINARY, 
            BufferedImage.TYPE_BYTE_INDEXED, BufferedImage.TYPE_USHORT_GRAY, 
            BufferedImage.TYPE_USHORT_565_RGB, BufferedImage.TYPE_USHORT_555_RGB, BufferedImage.TYPE_CUSTOM
            */
            
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
        
        for (int i = 0; i < edge.getN(); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                
                float xx = x + dx;
                
                if ((xx > -1) && (xx < (input.getWidth() - 1))) {
                    for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); 
                        dy++) {
                        
                        float yy = y + dy;
                        if ((yy > -1) && (yy < (input.getHeight() - 1))) {
                            input.setRGB((int)xx, (int)yy, rClr, gClr, bClr);
                        }
                    }
                }
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
        
        for (int i = 0; i < nPointsInEdge; i++) {
            
            int x = edgeX[i];
            int y = edgeY[i];
            
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                
                float xx = x + dx;
                
                if ((xx > -1) && (xx < (input.getWidth() - 1))) {
                    for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); 
                        dy++) {
                        
                        float yy = y + dy;
                        if ((yy > -1) && (yy < (input.getHeight() - 1))) {
                            input.setRGB((int)xx, (int)yy, rClr, gClr, bClr);
                        }
                    }
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
        List<PairIntArray> curves, String fileName, boolean writeImage, 
        Image input) throws IOException {
        
        if (curves == null || input == null) {
            return;
        }
        if (fileName == null) {
            throw new IllegalArgumentException("fileName cannot be null");
        }
        
        int clr = 0;
        
        for (int i = 0; i < curves.size(); i++) {
            PairIntArray edge = curves.get(i);
            if (clr > 5) {
                clr = 0;
            }
            int c = Color.BLUE.getRGB();
            switch(clr) {
                case 1:
                    c = Color.PINK.getRGB();
                    break;
                case 2:
                    c = Color.GREEN.getRGB();
                    break;
                case 3:
                    c = Color.RED.getRGB();
                    break;
                case 4:
                    c = Color.CYAN.getRGB();
                    break;
                case 5:
                    c = Color.MAGENTA.getRGB();
                    break;
            }
            for (int ii = 0; ii < edge.getN(); ii++) {
                input.setRGB(edge.getX(ii), edge.getY(ii), c);
            }
            clr++;
        }
       
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

    /**
     * create an image with edges2 in green and edges1 in white and write
     * it to fileName.
     * @param edges1
     * @param edges2
     * @param img
     * @param fileName 
     */
    public static void writeImage(PairIntArray[] edges1,
        SearchableCurve[] edges2, Image img, String fileName) {
        
        if (edges1 == null || edges2 == null) {
            return;
        }
        if (edges1.length != edges2.length) {
            throw new IllegalArgumentException(
                "edgeX and edgeY should be same length");
        }
        if (fileName == null) {
            throw new IllegalArgumentException(
                "fileName cannot be null");
        }
        
        try {
            for (int i = 0; i < edges2.length; i++) {
                SearchableCurve edge = edges2[i];
                addCurveToImage(edge.getX(), edge.getY(), edge.getN(),
                    img, 0, 0, 255, 0);
            }

            for (int i = 0; i < edges1.length; i++) {
                PairIntArray edge = edges1[i];
                addCurveToImage(edge.getX(), edge.getY(), edge.getN(),
                    img, 0, 255, 255, 255);
            }

            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.writeOutputImage(
                dirPath + "/" + fileName, img);
            
        } catch(IOException e) {
            Logger.getLogger(ImageIOHelper.class.getName()).severe(
                e.getMessage());
        }
    }
    
    public static void addCurveToImage(ScaleSpaceCurve curve, Image input, 
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        if (curve == null || input == null) {
            return;
        }
        
        for (int i = 0; i < curve.getSize(); i++) {
            int x = (int) curve.getX(i);
            int y = (int) curve.getY(i);
            
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                float xx = x + dx;
                if ((xx > -1) && (xx < (input.getWidth() - 1))) {
                    for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); dy++) {
                        float yy = y + dy;
                        if ((yy > -1) && (yy < (input.getHeight() - 1))) {
                            input.setRGB((int)xx, (int)yy, rClr, gClr, bClr);
                        }
                    }
                }
            }
        }
    }
    
}
