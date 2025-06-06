package algorithms.imageProcessing;

import algorithms.compGeometry.BresenhamsLine;
import algorithms.imageProcessing.features.CornerRegion;
import algorithms.imageProcessing.features.FeatureComparisonStat;
import algorithms.imageProcessing.scaleSpace.CurvatureScaleSpaceContour;
import algorithms.imageProcessing.scaleSpace.CurvatureScaleSpaceImagePoint;
import algorithms.imageProcessing.scaleSpace.ScaleSpaceCurve;
import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import gnu.trove.TIntCollection;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.set.TIntSet;
import java.awt.Color;
import java.awt.color.ColorSpace;
import java.awt.image.BufferedImage;
import java.awt.image.ColorConvertOp;
import java.awt.image.Raster;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import javax.imageio.ImageReadParam;
import javax.imageio.ImageReader;
import javax.imageio.stream.FileImageInputStream;
import javax.swing.ImageIcon;

/**
 *
 * @author nichole
 */
public class ImageIOHelper {

    public static Image readImage(String filePath) throws IOException {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            String lstPth = null;
            Iterator<Path> iter0 = file.toPath().iterator();
            while (iter0.hasNext()) {
                lstPth = iter0.next().toString();
            }
            
            if (lstPth.endsWith("jpg") || lstPth.endsWith("jpeg") 
                || lstPth.endsWith("jPG") || lstPth.endsWith("JPEG")) {
                return readJPEGImage(filePath);
            }
            
            BufferedImage img;
            if (filePath.endsWith("png") || filePath.endsWith("PNG")) {
                ImageIcon imageIcon = new ImageIcon(filePath);
                java.awt.Image tmpImage = imageIcon.getImage();

                img = new BufferedImage(
                    imageIcon.getIconWidth(), imageIcon.getIconHeight(), BufferedImage.TYPE_INT_ARGB);
                img.getGraphics().drawImage(tmpImage, 0, 0, null);
                tmpImage.flush();
            } else {
                img = ImageIO.read(file);
            }
            
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
    
    public static Image readJPEGImage(String filePath) throws IOException {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        File file = null;
        
        FileImageInputStream in = null;
        
        try {
            file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            String lstPth = null;
            Iterator<Path> iter0 = file.toPath().iterator();
            while (iter0.hasNext()) {
                lstPth = iter0.next().toString();
            }
            
            assert(lstPth.endsWith("jpg") || lstPth.endsWith("jpeg") 
                || lstPth.endsWith("jPG") || lstPth.endsWith("JPEG"));
            
            in = new FileImageInputStream(file);
            
            Iterator<ImageReader> iter = 
                ImageIO.getImageReadersByMIMEType("image/jpeg");
                    
            Image output = null;
            
            while (iter.hasNext()) {
                
                ImageReader rdr = iter.next();
                System.out.println("jpegreader instance " + 
                    rdr.getClass().getName());
                
                rdr.setInput(in);
                
                ImageReadParam rdrParam = rdr.getDefaultReadParam();
               
                int nImages = rdr.getNumImages(true);
                
                System.out.println("nJPEG bands=" + nImages);
                
                assert(nImages >= 1);
                                
                int w = rdr.getWidth(0);
                int h = rdr.getHeight(0);
                                
                output = new Image(w, h);
                
                BufferedImage im = rdr.read(0, rdrParam);
                
                int[] col = new int[w];
                int offset = 0;
                for (int y = 0; y < h; ++y) {
                    im.getRGB(0, y, w, 1, col, offset, 1);
                    for (int x = 0; x < w; ++x) {
                        int rgb = col[x];
                        
                        int r = (rgb >> 16) & 0xFF;
                        int g = (rgb >> 8) & 0xFF;
                        int b = rgb & 0xFF;

                        output.setRGB(x, y, r, g, b);
                    }
                }
               
                return output;
            }
          
            return output;
            
        } catch (IOException e) {
      
            System.err.println(e.getMessage());
            
            if (in != null) {
                in.close();
            }
        }

        return null;        
    }
    
    public static GreyscaleImage readImageAsGreyscaleFullRange(String filePath) 
        throws IOException {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            BufferedImage img;
            if (filePath.endsWith("png") || filePath.endsWith("PNG")) {
                ImageIcon imageIcon = new ImageIcon(filePath);
                java.awt.Image tmpImage = imageIcon.getImage();

                img = new BufferedImage(
                    imageIcon.getIconWidth(), imageIcon.getIconHeight(), BufferedImage.TYPE_INT_ARGB);
                img.getGraphics().drawImage(tmpImage, 0, 0, null);
                tmpImage.flush();
            } else {
                img = ImageIO.read(file);
            }
            
            //System.out.println("imageType=" + img.getType());
            
            int h = img.getHeight();
            int w = img.getWidth();
            
            ColorSpace cs = ColorSpace.getInstance(ColorSpace.CS_GRAY);  
            ColorConvertOp op = new ColorConvertOp(cs, null);
            try {
                BufferedImage image2 = op.filter(img, null);
                img = image2;
            } catch (NullPointerException e) {
                // if type is custom, the source color space destination is not defined.
            } 
            
            GreyscaleImage image = new GreyscaleImage(w, h, 
                GreyscaleImage.Type.Bits32FullRangeInt);
                        
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
   
    public static ImageExt readImageExt(String filePath) throws IOException {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
            
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            BufferedImage img;
            if (filePath.endsWith("png") || filePath.endsWith("PNG")) {
                ImageIcon imageIcon = new ImageIcon(filePath);
                java.awt.Image tmpImage = imageIcon.getImage();

                img = new BufferedImage(
                    imageIcon.getIconWidth(), imageIcon.getIconHeight(), BufferedImage.TYPE_INT_ARGB);
                img.getGraphics().drawImage(tmpImage, 0, 0, null);
                tmpImage.flush();
            } else {
                img = ImageIO.read(file);
            }
            
            //System.out.println("imageType=" + img.getType());
            
            int h = img.getHeight();
            int w = img.getWidth();
            
            ImageExt image = new ImageExt(w, h);
            
            convertImage(img, image);
            
            return image;
            
        } catch (IOException e) {
        }
        
        return null;
    }
    
    private static void convertImage(BufferedImage fromImage, Image toImage) 
        throws IOException {
     
        if (fromImage == null) {
            return;
        }
               
        int h = fromImage.getHeight();
        int w = fromImage.getWidth();

        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {

                //TYPE_INT_RGB) and default sRGB colorspace
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
        throws IOException {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            BufferedImage img;
            if (filePath.endsWith("png") || filePath.endsWith("PNG")) {
                ImageIcon imageIcon = new ImageIcon(filePath);
                java.awt.Image tmpImage = imageIcon.getImage();

                img = new BufferedImage(
                    imageIcon.getIconWidth(), imageIcon.getIconHeight(), BufferedImage.TYPE_INT_ARGB);
                img.getGraphics().drawImage(tmpImage, 0, 0, null);
                tmpImage.flush();
            } else {
                img = ImageIO.read(file);
            }
            
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
        throws IOException {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            BufferedImage img;
            if (filePath.endsWith("png") || filePath.endsWith("PNG")) {
                ImageIcon imageIcon = new ImageIcon(filePath);
                java.awt.Image tmpImage = imageIcon.getImage();

                img = new BufferedImage(
                    imageIcon.getIconWidth(), imageIcon.getIconHeight(), BufferedImage.TYPE_INT_ARGB);
                img.getGraphics().drawImage(tmpImage, 0, 0, null);
                tmpImage.flush();
            } else {
                img = ImageIO.read(file);
            }
                       
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
        throws IOException {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            BufferedImage img;
            if (filePath.endsWith("png") || filePath.endsWith("PNG")) {
                ImageIcon imageIcon = new ImageIcon(filePath);
                java.awt.Image tmpImage = imageIcon.getImage();

                img = new BufferedImage(
                    imageIcon.getIconWidth(), imageIcon.getIconHeight(), BufferedImage.TYPE_INT_ARGB);
                img.getGraphics().drawImage(tmpImage, 0, 0, null);
                tmpImage.flush();
            } else {
                img = ImageIO.read(file);
            }
                       
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
     * @throws IOException 
     */
    public static Image readImageAsGrayScale(String filePath) 
        throws IOException {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            BufferedImage img;
            if (filePath.endsWith("png") || filePath.endsWith("PNG")) {
                ImageIcon imageIcon = new ImageIcon(filePath);
                java.awt.Image tmpImage = imageIcon.getImage();

                img = new BufferedImage(
                    imageIcon.getIconWidth(), imageIcon.getIconHeight(), BufferedImage.TYPE_INT_ARGB);
                img.getGraphics().drawImage(tmpImage, 0, 0, null);
                tmpImage.flush();
            } else {
                img = ImageIO.read(file);
            }

            //ImageDisplayer.displayDisposableImage("buffered", img);

            ColorSpace cs = ColorSpace.getInstance(ColorSpace.CS_GRAY);  
            ColorConvertOp op = new ColorConvertOp(cs, null);
            try {
                BufferedImage image2 = op.filter(img, null);
                img = image2;
 
                //ImageDisplayer.displayDisposableImage("buffered2", img);
 
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
       
            //ImageDisplayer.displayImageGrey("image", image);
       
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
        throws IOException {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            BufferedImage img;
            if (filePath.endsWith("png") || filePath.endsWith("PNG")) {
                ImageIcon imageIcon = new ImageIcon(filePath);
                java.awt.Image tmpImage = imageIcon.getImage();

                img = new BufferedImage(
                    imageIcon.getIconWidth(), imageIcon.getIconHeight(), BufferedImage.TYPE_INT_ARGB);
                img.getGraphics().drawImage(tmpImage, 0, 0, null);
                tmpImage.flush();
            } else {
                img = ImageIO.read(file);
            }
                 
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
    public static GreyscaleImage readImageAsGrayScaleAvgRGB(String filePath)
        throws IOException {
     
        if (filePath == null) {
            throw new IllegalStateException("filePath cannot be null");
        }
                
        try {
            File file = new File(filePath);
            if (!file.exists()) {
                throw new IllegalStateException(filePath + " does not exist");
            }
            
            BufferedImage img;
            if (filePath.endsWith("png") || filePath.endsWith("PNG")) {
                ImageIcon imageIcon = new ImageIcon(filePath);
                java.awt.Image tmpImage = imageIcon.getImage();

                img = new BufferedImage(
                    imageIcon.getIconWidth(), imageIcon.getIconHeight(), BufferedImage.TYPE_INT_ARGB);
                img.getGraphics().drawImage(tmpImage, 0, 0, null);
                tmpImage.flush();
            } else {
                img = ImageIO.read(file);
            }
            
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
    
    public static String writeOutputImage(String filePath, Image data) 
        throws IOException {
        
        return writeOutputImage(filePath, data,  BufferedImage.TYPE_INT_RGB);       
    }
    
    public static String writeOutputImage(String filePath, Image data, 
        final int imageType) throws IOException {
                
        if (data == null) {
            return "";
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
        
        return filePath.toString();
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
    
    public static void writeOutputGreyscaleImage(String filePath, Image data) 
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
                
                int rgb = data.getRGB(i, j);
                
                int r = (rgb >> 16) & 0xFF;
                int g = (rgb >> 8) & 0xFF;
                int b = rgb & 0xFF;  
                
                // from wikipedia, srgb conversion
                double yLinear = 0.2126 * r + 0.7152 * g + 0.0722 * b;
                /*if (yLinear > 0.0031308) {
                    yLinear = 1.033 * Math.pow(yLinear, 1./2.4) - 0.055f;
                } else {
                    yLinear *= 12.92;
                }*/
                int v = (int)Math.round(yLinear);
                
                raster.setSample(i, j, 0, v);
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
    
    public static void addContoursToImage(List<CurvatureScaleSpaceContour> result, 
        ImageExt img, int xOffset, int yOffset,
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        if (result.isEmpty()) {
            return;
        }
        
        for (int i = 0; i < result.size(); i++) {
            
            CurvatureScaleSpaceContour cssC = result.get(i);
            
            addContoursToImage(cssC, img, xOffset, yOffset, nExtraForDot, 
                rClr, gClr, bClr);
        }
    }
    
    public static void addContoursToImage(CurvatureScaleSpaceContour contour, 
        ImageExt img, int xOffset, int yOffset,
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        if (contour == null) {
            return;
        }
                    
        CurvatureScaleSpaceImagePoint[] peakDetails = contour.getPeakDetails();

        for (CurvatureScaleSpaceImagePoint peakDetail : peakDetails) {
            
            int x = peakDetail.getXCoord() + xOffset;
            int y = peakDetail.getYCoord() + yOffset;
            
            addPointToImage(x, y, img, xOffset, yOffset, nExtraForDot, 
                rClr, gClr, bClr);
        }            
    }
    
    public static void addPointToImage(final int x, final int y, 
        Image img, int xOffset, int yOffset,
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
            float xx = x + dx;
            if ((xx > -1) && (xx < (img.getWidth() - 1))) {
                for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); 
                    dy++) {
                    float yy = y + dy;
                    if ((yy > -1) && (yy < (img.getHeight() - 1))) {
                        img.setRGB((int)xx, (int)yy, rClr, gClr, bClr);
                    }
                }
            }
        }
    }
    
    private static PairFloatArray getXYOfContourPeaks(
        List<CurvatureScaleSpaceContour> contours, int xOffset, int yOffset) {
        
        PairFloatArray xy = new PairFloatArray();
        
        for (int i = 0; i < contours.size(); ++i) {
            
            CurvatureScaleSpaceContour contour = contours.get(i);
            
            CurvatureScaleSpaceImagePoint[] peakDetails = contour.getPeakDetails();
            
            for (CurvatureScaleSpaceImagePoint peakDetail : peakDetails) {
                int x = peakDetail.getXCoord() + xOffset;
                int y = peakDetail.getYCoord() + yOffset;
                xy.add(x, y);
            }
        }
        
        return xy;
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
     * @param pixIdxs
     * @param input
     * @param nExtraForDot
     * @param rClr
     * @param gClr
     * @param bClr 
     */
    public static void addCurveToImage(TIntCollection pixIdxs, Image input, 
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        if (pixIdxs == null || input == null) {
            return;
        }
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        TIntIterator iter = pixIdxs.iterator();
        while (iter.hasNext()) {
            
            int pixIdx = iter.next();
            
            int y = pixIdx/w;
            int x = pixIdx - (y * w);
            
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
     * draw the edge over the image using the given rgb colors and the size
     * of the dot beyond 1 pixel.
     * 
     * @param pointGroups
     * @param input
     * @param nExtraForDot
     */
    public static <T extends Collection<PairInt>> void addAlternatingColorCurvesToImage0(
        List<T> pointGroups, Image input, int nExtraForDot) {
        
        if (pointGroups == null || input == null) {
            return;
        }
        
        for (int i = 0; i < pointGroups.size(); ++i) {
            
            int[] clr = getNextRGB(i);
            
            addCurveToImage(pointGroups.get(i), input, nExtraForDot,
                clr[0], clr[1], clr[2]);
        }
        
    }
    
    /**
     * draw the edge over the image using the given rgb colors and the size
     * of the dot beyond 1 pixel.
     * 
     * @param pointGroups
     * @param input
     * @param nExtraForDot
     */
    public static void addAlternatingColorCurvesToImage3(
        List<TIntSet> pointGroups, Image input, int nExtraForDot) {
        
        if (pointGroups == null || input == null) {
            return;
        }
        
        for (int i = 0; i < pointGroups.size(); ++i) {
            
            int[] clr = getNextRGB(i);
            
            addCurveToImage(pointGroups.get(i), input, nExtraForDot,
                clr[0], clr[1], clr[2]);
        }
        
    }
    
    /**
     * draw lines from points at index i to index i+1 
     * in the requested color onto the input.
     * 
     * @param xVertexes
     * @param yVertexes
     * @param input
     * @param nExtraForDot
     * @param rClr
     * @param gClr
     * @param bClr 
     */
    public static void drawLinesInImage(
        int[] xVertexes, int[] yVertexes, 
        Image input, int nExtraForDot, 
        int rClr, int gClr, int bClr) {
        
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
            
            drawLineInImage(x1, y1, x2, y2, input, 
                nExtraForDot, rClr, gClr, bClr);
        }
    }
    
    /**
     * draw lines from points at index i to index i+1 in the requested color
     * onto the input.
     * 
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param input
     * @param nExtraForDot
     */
    public static void drawLineInImage(int x1, int y1, 
        int x2, int y2,
        Image input, int nExtraForDot, int rgbClr) {
        
        int r = (rgbClr >> 16) & 0xFF;
        int g = (rgbClr >> 8) & 0xFF;
        int b = rgbClr & 0xFF;
        
        drawLineInImage(x1, y1, x2, y2, input, nExtraForDot, r, g, b);
    }
    
    /**
     * draw lines from points at index i to index i+1 in the requested color
     * onto the input.
     * 
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param input
     * @param nExtraForDot
     * @param rClr
     * @param gClr
     * @param bClr 
     */
    public static void drawLineInImage(int x1, int y1, 
        int x2, int y2,
        Image input, int nExtraForDot, int rClr, int gClr, 
        int bClr) {
        
        int w = input.getWidth();
        int h = input.getHeight();

        Set<PairInt> output = new HashSet<PairInt>();
        
        BresenhamsLine.createLinePoints(x1, y1, x2, y2,
            output);
        
        for (PairInt p : output) {

            int x = p.getX();
            int y = p.getY();
            
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); 
                dx++) {

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
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param input
     * @param nExtraForDot
     * @param rClr
     * @param gClr
     * @param bClr 
     */
    public static void drawLineInImage(int x1, int y1, 
        int x2, int y2, Image input, int nExtraForDot, int gapLength,
        int rClr, int gClr, int bClr) {
        
        int w = input.getWidth();
        int h = input.getHeight();

        Set<PairInt> output = new HashSet<PairInt>();
        
        BresenhamsLine.createLinePoints(x1, y1, x2, y2, gapLength,
            output);
        
        for (PairInt p : output) {

            int x = p.getX();
            int y = p.getY();
            
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); 
                dx++) {

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
     * @param x of point (NOTE that if the points are from a row major reference frame,
     *          you will want to use 'y' here instead).
     * @param y of point (NOTE that if the points are from a row major reference frame,
     *      *          you will want to use 'x' here instead).
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
                input.setRGB(xx, yy, rClr, gClr, bClr);
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
    
    public static void addAlternatingColorLabelsToRegion(Image img, int[] labels) {
    
        int nPixels = img.getNPixels();
        
        if (labels.length != nPixels) {
            throw new IllegalArgumentException(
                "labels.length should == img.nPixels");
        }
        
        Map<Integer, Integer> labelColorMap = new HashMap<Integer, Integer>();
        
        int count = 0;
        
        for (int i = 0; i < nPixels; ++i) {
            
            int label = labels[i];
            
            Integer key = Integer.valueOf(label);            
            
            Integer clr = labelColorMap.get(key);
            
            if (clr == null) {
                                
                clr = Integer.valueOf(getNextColorRGB(count));
                
                labelColorMap.put(key, clr);
                count++;
            }
            
            img.setRGB(img.getCol(i), img.getRow(i), clr);
        }
    }
    
    public static void addAlternatingColorLabelsToRegion(
        Image img, TIntObjectMap<TIntSet> labelsMap) {
                    
        TIntObjectIterator<TIntSet> iter = labelsMap.iterator();
        for (int i = 0; i < labelsMap.size(); ++i) {

            iter.advance();
            
            int clr = getNextColorRGB(i);
            
            TIntSet indexes = iter.value();
            TIntIterator iter2 = indexes.iterator();
            while (iter2.hasNext()) {
                int pixIdx = iter2.next();
                img.setRGB(img.getCol(pixIdx), 
                    img.getRow(pixIdx), clr);
            }
        }
    }
    
    public static void addAlternatingColorPointsToImages(
        Map<PairInt, List<PairInt>> points, Image imgCp1, 
        Image imgCp2, int nExtraForDot) {
        
        int w1 = imgCp1.getWidth();
        int h1 = imgCp1.getHeight();
        int w2 = imgCp2.getWidth();
        int h2 = imgCp2.getHeight();
        
        int count = 0;
        for (Entry<PairInt, List<PairInt>> entry : points.entrySet()) {
            int clr = getNextColorRGB(count);
            PairInt p1 = entry.getKey();
            int x1 = p1.getX();
            int y1 = p1.getY();
            for (int dx = (-1 * nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                int xx1 = x1 + dx;
                if ((xx1 < 0) || (xx1 > (w1 - 1))) {
                    continue;
                }
                for (int dy = (-1 * nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                    int yy1 = y1 + dy;
                    if ((yy1 < 0) || (yy1 > (h1 - 1))) {
                        continue;
                    }
                    imgCp1.setRGB(xx1, yy1, clr);
                }
            }                    
            for (PairInt p2 : entry.getValue()) {
                int x2 = p2.getX();
                int y2 = p2.getY();
                for (int dx = (-1 * nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                    int xx2 = x2 + dx;
                    if ((xx2 < 0) || (xx2 > (w2 - 1))) {
                        continue;
                    }
                    for (int dy = (-1 * nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                        int yy2 = y2 + dy;
                        if ((yy2 < 0) || (yy2 > (h2 - 1))) {
                            continue;
                        }
                        imgCp2.setRGB(xx2, yy2, clr);
                    }
                }
            }
            count++;
        }
    }
    
    public static void addAlternatingColorPointsToImages(
        List<FeatureComparisonStat> stats, Image imgCp1, 
        Image imgCp2, int nExtraForDot) {
        
        int w1 = imgCp1.getWidth();
        int h1 = imgCp1.getHeight();
        int w2 = imgCp2.getWidth();
        int h2 = imgCp2.getHeight();
        
        int count = 0;
        for (FeatureComparisonStat stat : stats) {
            int clr = getNextColorRGB(count);
            PairInt p1 = stat.getImg1Point();
            PairInt p2 = stat.getImg2Point();
            int x1 = p1.getX();
            int y1 = p1.getY();
            for (int dx = (-1 * nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                int xx1 = x1 + dx;
                if ((xx1 < 0) || (xx1 > (w1 - 1))) {
                    continue;
                }
                for (int dy = (-1 * nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                    int yy1 = y1 + dy;
                    if ((yy1 < 0) || (yy1 > (h1 - 1))) {
                        continue;
                    }
                    imgCp1.setRGB(xx1, yy1, clr);
                }
            }                    
            int x2 = p2.getX();
            int y2 = p2.getY();
            for (int dx = (-1 * nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                int xx2 = x2 + dx;
                if ((xx2 < 0) || (xx2 > (w2 - 1))) {
                    continue;
                }
                for (int dy = (-1 * nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                    int yy2 = y2 + dy;
                    if ((yy2 < 0) || (yy2 > (h2 - 1))) {
                        continue;
                    }
                    imgCp2.setRGB(xx2, yy2, clr);
                }
            }
            count++;
        }
    }
    
    public static void addAlternatingColorCurvesToImage(
        List<PairIntArray> curves, Image input, int nExtraForDot) {
        
        if (curves == null || input == null) {
            return;
        }
        
        addAlternatingColorCurvesToImage(
            curves.toArray(new PairIntArray[curves.size()]),
            input, nExtraForDot);
    }
    
    public static void addAlternatingColorCurvesToImage(
        List<PairIntArray> curves, Image input, int xOffset, int yOffset,
        int nExtraForDot) {
        
        if (curves == null || input == null) {
            return;
        }
        
        if (curves == null || input == null || curves.isEmpty()) {
            return;
        }
                
        int clr = 0;
        int w = input.getWidth();
        int h = input.getHeight();
        
        for (int i = 0; i < curves.size(); i++) {
            
            PairIntArray edge = curves.get(i);
            
            int c = getNextColorRGB(clr);
            
            for (int ii = 0; ii < edge.getN(); ii++) {
                
                int col = edge.getX(ii) + xOffset;
                int row = edge.getY(ii) + yOffset;
                
                for (int dx = (-1 * nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                    int xx = col + dx;
                    if ((xx < 0) || (xx > (w - 1))) {
                        continue;
                    }
                    for (int dy = (-1 * nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                        int yy = row + dy;
                        if ((yy < 0) || (yy > (h - 1))) {
                            continue;
                        }
                        input.setRGB(xx, yy, c);
                    }
                }                
            }
            clr++;
        }
    }
    
    public static void addAlternatingColorCurvesToImage2(
        List<PairIntArray> curves, Image input, int nExtraForDot) {
        
        if (curves == null || input == null) {
            return;
        }
        
        addAlternatingColorCurvesToImage(
            curves.toArray(new PairIntArray[curves.size()]),
            input, nExtraForDot);
    }
    
    public static void addAlternatingColorCurvesToImage(
        List<PairIntArray> curves, Image input) {
        
        if (curves == null || input == null) {
            return;
        }
        
        addAlternatingColorCurvesToImage(
            curves.toArray(new PairIntArray[curves.size()]),
            input, 1);
    }
    
    public static int getNextColorRGB(int count) {
        
        int[] c = getNextRGB(count);
        
        int rgb = (((c[0] & 0x0ff) << 16) | ((c[1] & 0x0ff) << 8) | 
            (c[2] & 0x0ff));
        
        return rgb;
    }
    
    /*
    http://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
    
    */
    public static int[] getNextRGB(int count) {
        
        //http://godsnotwheregodsnot.blogspot.ru/2013/11/kmeans-color-quantization-seeding.html
        String[] indexcolors = new String[]{
            "FFFF00", "1CE6FF", "FF34FF", "FF4A46", "008941", "006FA6", "A30059",
            "FFDBE5", "7A4900", "0000A6", "63FFAC", "B79762", "004D43", "8FB0FF", "997D87",
            "5A0007", "809693", "FEFFE6", "1B4400", "4FC601", "3B5DFF", "4A3B53", "FF2F80",
            "61615A", "BA0900", "6B7900", "00C2A0", "FFAA92", "FF90C9", "B903AA", "D16100",
            "DDEFFF", "000035", "7B4F4B", "A1C299", "300018", "0AA6D8", "013349", "00846F",
            "372101", "FFB500", "C2FFED", "A079BF", "CC0744", "C0B9B2", "C2FF99", "001E09",
            "00489C", "6F0062", "0CBD66", "EEC3FF", "456D75", "B77B68", "7A87A1", "788D66",
            "885578", "FAD09F", "FF8A9A", "D157A0", "BEC459", "456648", "0086ED", "886F4C",
            "34362D", "B4A8BD", "00A6AA", "452C2C", "636375", "A3C8C9", "FF913F", "938A81",
            "575329", "00FECF", "B05B6F", "8CD0FF", "3B9700", "04F757", "C8A1A1", "1E6E00",
            "7900D7", "A77500", "6367A9", "A05837", "6B002C", "772600", "D790FF", "9B9700",
            "549E79", "FFF69F", "201625", "72418F", "BC23FF", "99ADC0", "3A2465", "922329",
            "5B4534", "FDE8DC", "404E55", "0089A3", "CB7E98", "A4E804", "324E72", "6A3A4C",
            "83AB58", "001C1E", "D1F7CE", "004B28", "C8D0F6", "A3A489", "806C66", "222800",
            "BF5650", "E83000", "66796D", "DA007C", "FF1A59", "8ADBB4", "1E0200", "5B4E51",
            "C895C5", "320033", "FF6832", "66E1D3", "CFCDAC", "D0AC94", "7ED379", "012C58",
            "7A7BFF", "D68E01", "353339", "78AFA1", "FEB2C6", "75797C", "837393", "943A4D",
            "B5F4FF", "D2DCD5", "9556BD", "6A714A", "001325", "02525F", "0AA3F7", "E98176",
            "DBD5DD", "5EBCD1", "3D4F44", "7E6405", "02684E", "962B75", "8D8546", "9695C5",
            "E773CE", "D86A78", "3E89BE", "CA834E", "518A87", "5B113C", "55813B", "E704C4",
            "00005F", "A97399", "4B8160", "59738A", "FF5DA7", "F7C9BF", "643127", "513A01",
            "6B94AA", "51A058", "A45B02", "1D1702", "E20027", "E7AB63", "4C6001", "9C6966",
            "64547B", "97979E", "006A66", "391406", "F4D749", "0045D2", "006C31", "DDB6D0",
            "7C6571", "9FB2A4", "00D891", "15A08A", "BC65E9", "FFFFFE", "C6DC99", "203B3C",
            "671190", "6B3A64", "F5E1FF", "FFA0F2", "CCAA35", "374527", "8BB400", "797868",
            "C6005A", "3B000A", "C86240", "29607C", "402334", "7D5A44", "CCB87C", "B88183",
            "AA5199", "B5D6C3", "A38469", "9F94F0", "A74571", "B894A6", "71BB8C", "00B433",
            "789EC9", "6D80BA", "953F00", "5EFF03", "E4FFFC", "1BE177", "BCB1E5", "76912F",
            "003109", "0060CD", "D20096", "895563", "29201D", "5B3213", "A76F42", "89412E",
            "1A3A2A", "494B5A", "A88C85", "F4ABAA", "A3F3AB", "00C6C8", "EA8B66", "958A9F",
            "BDC9D2", "9FA064", "BE4700", "658188", "83A485", "453C23", "47675D", "3A3F00",
            "061203", "DFFB71", "868E7E", "98D058", "6C8F7D", "D7BFC2", "3C3E6E", "D83D66",
            "2F5D9B", "6C5E46", "D25B88", "5B656C", "00B57F", "545C46", "866097", "365D25",
            "252F99", "00CCFF", "674E60", "FC009C", "92896B"
        };
        
        int n = indexcolors.length - 1;
        
        if (count < 0) {
            count = 0;
        } else {
            count = count % n;
        }
    
        String str = indexcolors[count];
        
        int rgb = Integer.parseInt(str, 16);
        int r = (rgb >> 16) & 0xFF;
        int g = (rgb >> 8) & 0xFF;
        int b = rgb & 0xFF;
                
        return new int[]{r, g, b};
    }
    
    public static void addAlternatingColorCurvesToImage(
        PairIntArray[] curves, Image input) {
        
        if (curves == null || input == null) {
            return;
        }
        
        int nExtraForDot = 1;
        
        addAlternatingColorCurvesToImage(curves, input, nExtraForDot);
    }
    
    public static void addAlternatingColorCurvesToImage(
        PairIntArray[] curves, Image input, int nExtraForDot) {
        
        if (curves == null || input == null) {
            return;
        }
                
        int clr = 0;
        int w = input.getWidth();
        int h = input.getHeight();
        
        for (int i = 0; i < curves.length; i++) {
            
            PairIntArray edge = curves[i];
            
            int c = getNextColorRGB(clr);
            
            for (int ii = 0; ii < edge.getN(); ii++) {
                
                int col = edge.getX(ii);
                int row = edge.getY(ii);
                
                for (int dx = (-1 * nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                    int xx = col + dx;
                    if ((xx < 0) || (xx > (w - 1))) {
                        continue;
                    }
                    for (int dy = (-1 * nExtraForDot); dy < (nExtraForDot + 1); ++dy) {
                        int yy = row + dy;
                        if ((yy < 0) || (yy > (h - 1))) {
                            continue;
                        }
                        input.setRGB(xx, yy, c);
                    }
                }                
            }
            clr++;
        }
    }
    
    public static void addAlternatingColorPointSetsToImage(
        List<Set<PairInt>> pointSets, int xOffset, int yOffset, 
        int nExtraForDot, Image input) {
        
        if (pointSets == null || input == null) {
            return;
        }
                
        for (int i = 0; i < pointSets.size(); i++) {
            
            Set<PairInt> points = pointSets.get(i);
            
            int[] c = getNextRGB(i);
            
            addToImage(points, 0, 0, input, nExtraForDot, c[0], c[1], c[2]);
        }
    }
    
    public static void addAlternatingColorPointSetsToImage2(
        List<TIntSet> pointSets, int xOffset, int yOffset, 
        int nExtraForDot, Image input) {
        
        if (pointSets == null || input == null) {
            return;
        }
                
        for (int i = 0; i < pointSets.size(); i++) {
            
            TIntSet pIdxs = pointSets.get(i);
            
            int[] c = getNextRGB(i);
            
            addToImage(pIdxs, 0, 0, input, nExtraForDot, c[0], c[1], c[2]);
        }
    }
    
    public static <T extends Collection<PairInt>> void addToImage(
        T points, int xOffsetToApply, int yOffsetToApply, Image input) {
        
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
    
    public static <T extends Collection<PairInt>> void addToImage(
        T points, int xOffsetToApply, 
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
    
    public static <T extends Collection<PairInt>> void addToImage(
        T points, int xOffsetToApply, 
        int yOffsetToApply, Image input, int nExtraForDot,
        int rClr, int gClr, int bClr) {
        
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

    public static void addToImage(
        TIntSet pIdxs, int xOffsetToApply, 
        int yOffsetToApply, Image input, int nExtraForDot,
        int rClr, int gClr, int bClr) {
        
        if (pIdxs == null || input == null) {
            return;
        }
        
        int w = input.getWidth();
        int h = input.getHeight();
             
        TIntIterator iter = pIdxs.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int y = pixIdx/w;
            int x = pixIdx - (y * w);
            x += xOffsetToApply;
            y += yOffsetToApply;
            
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

    public static void addAlternatingColorContoursToImage(
        List<CurvatureScaleSpaceContour> contours, ImageExt img, 
        int xOffset, int yOffset, int nExtraForDot) {
        
        //StringBuilder sb = new StringBuilder();
        
        for (int i = 0; i < contours.size(); i++) {
            
            CurvatureScaleSpaceContour cssC = contours.get(i);
            
            /*sb.append("[").append(Integer.toString(i))
                .append("] t=")
                .append(Float.toString(cssC.getPeakScaleFreeLength()))
                .append(", s=")
                .append(Float.toString(cssC.getPeakSigma()))
                .append(" ");*/
            
            int[] c = getNextRGB(i);
            
            addContoursToImage(cssC, img, xOffset, yOffset, nExtraForDot, 
                c[0], c[1], c[2]);
        }
        
        //Logger.getLogger(ImageIOHelper.class.getName()).info(sb.toString());
    }

    public static void addAlternatingColorCornerRegionsToImage(
        List<CornerRegion> regions, Image img, 
        int xOffset, int yOffset, int nExtraForDot) {
        
        for (int i = 0; i < regions.size(); i++) {
            
            CornerRegion br = regions.get(i);
            
            int[] c = getNextRGB(i);

            int x = br.getX()[br.getKMaxIdx()];
            int y = br.getY()[br.getKMaxIdx()];
            
            addPointToImage(x, y, img, xOffset, yOffset, nExtraForDot, 
                c[0], c[1], c[2]);
        }
    }
    
    public static <T extends CornerRegion> void 
        addAlternatingColorCornerRegionListsToImage(
        List<List<T>> regions, Image img, int xOffset, int yOffset, 
            int nExtraForDot) {
        
        for (int i = 0; i < regions.size(); i++) {
            
            List<T> corners = regions.get(i);
            
            int[] c = getNextRGB(i);
            
            addCornerRegionsToImage(corners, img, xOffset, yOffset, 
                nExtraForDot, c[0], c[1], c[2]);
        }
    }
    
    public static void addCornerRegionsToImage(Collection<CornerRegion> regions,
        Image img, int xOffset, int yOffset, int nExtraForDot) {

        for (CornerRegion br : regions) {

            int x = br.getX()[br.getKMaxIdx()];
            int y = br.getY()[br.getKMaxIdx()];

            addPointToImage(x, y, img, xOffset, yOffset, nExtraForDot, 255, 0, 0);
        }
    }
    
    public static <T extends CornerRegion> void addCornerRegionsToImage(
        Collection<T> regions, Image img, int xOffset, int yOffset, 
        int nExtraForDot, int rClr, int gClr, int bClr) {
        
        for (T br : regions ) {
            
            int x = br.getX()[br.getKMaxIdx()];
            int y = br.getY()[br.getKMaxIdx()];
            
            addPointToImage(x, y, img, xOffset, yOffset, nExtraForDot, 
                rClr, gClr, bClr);
        }
    }
}
