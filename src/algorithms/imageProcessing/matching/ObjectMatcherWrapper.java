package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.features.CorrespondenceList;
import algorithms.imageProcessing.features.ObjectMatcher;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * A class to pre-process images in order to use the ObjectMatcher
 * within a rough range of resolution.
 * The class currently expects color images.
 * 
 * @author nichole
 */
public class ObjectMatcherWrapper {
    
    private static int maxDimension = 256;
    private SIGMA sigma = SIGMA.ZEROPOINTFIVE;//SIGMA.ONE;

    private ImageExt[] templateImage = null;
    private ImageExt searchImage = null;
    
    private boolean debug = false;
    
    public void setToDebug() {
        debug = true;
    }
    
    public List<CorrespondenceList> find(String templateFilePath,
        String templateMaskFilePath,
        String searchFilePath, String debugLabel) throws Exception {
   
        Set<PairInt> shape0 = new HashSet<PairInt>();

        ImageExt[] imgs0 = maskAndBin2(templateFilePath, templateMaskFilePath,
            shape0);
            
        return find(imgs0, shape0, searchFilePath, debugLabel);
    }
    
    public List<CorrespondenceList> find(String templateFilePath,
        String searchFilePath, String debugLabel) throws Exception {
   
        Set<PairInt> shape0 = new HashSet<PairInt>();

        ImageExt img0 = maskAndBin2(templateFilePath, shape0);
            
        return find(new ImageExt[]{img0, img0}, shape0, searchFilePath, 
            debugLabel);
    }
    
    public List<CorrespondenceList> find(ImageExt[] binnedTemplateAndMask, 
        Set<PairInt> shape0,
        String searchFilePath, String debugLabel) throws IOException {
        
        ImageExt img = ImageIOHelper.readImageExt(searchFilePath);
        
        return find(binnedTemplateAndMask, shape0, img, debugLabel);
    }
    
    public List<CorrespondenceList> find(ImageExt[] binnedTemplateAndMask, 
        Set<PairInt> shape0, ImageExt searchImage, String debugLabel) throws IOException {
        
        int nShape0_0 = shape0.size();

        System.out.println("shape0 nPts=" + nShape0_0);

        long ts = MiscDebug.getCurrentTimeFormatted();

        searchImage = bin(searchImage);
        
        return _run_matcher(binnedTemplateAndMask, shape0, searchImage, debugLabel);
    }
   
    private List<CorrespondenceList> _run_matcher(ImageExt[] imgs0, Set<PairInt> shape0, 
        ImageExt img, String debugLabel) throws IOException {
        
        this.templateImage = imgs0;
        this.searchImage = img;

        ObjectMatcher.Settings settings = new ObjectMatcher.Settings();
        settings.setToUseLargerPyramid0();
        settings.setToUseLargerPyramid1();

        ObjectMatcher objMatcher = new ObjectMatcher();

        if (debug) {
            objMatcher.setToDebug();
        }
        if (debugLabel != null) {  
            settings.setDebugLabel(debugLabel);
        }

        //settings.setToExcludeColorFilter();

        long t0 = System.currentTimeMillis();
        
        List<CorrespondenceList> corresList 
            //= objMatcher.findObject11(
            = objMatcher.findObject12(imgs0[0], shape0, img, settings);

        long t1 = System.currentTimeMillis();
        System.out.println("matching took " + ((t1 - t0)/1000.) + " sec");

        return corresList;
    }
    
    public static ImageExt maskAndBin2(String templateFilePath, 
        Set<PairInt> outputShape) throws IOException {
        
        ImageProcessor imageProcessor = new ImageProcessor();

        ImageExt img0 = bin(ImageIOHelper.readImageExt(templateFilePath));
          
        for (int x = 0; x < img0.getWidth(); ++x) {
            for (int y = 0; y < img0.getHeight(); ++y) {
                if (img0.getRGB(x, y) != 0) {
                    outputShape.add(new PairInt(x, y));
                }
            }
        }
   
        return img0;
    }
    
    public static ImageExt bin(ImageExt img) throws IOException {
        
        ImageProcessor imageProcessor = new ImageProcessor();
    
        int w = img.getWidth();
        int h = img.getHeight();

        int binFactor = (int) Math.ceil(Math.max(
             (float) w / maxDimension,
             (float) h / maxDimension));
                
        if (binFactor != 1) {
            img = imageProcessor.binImage(img, binFactor);
        }
   
        return img;
    }

    public static ImageExt[] maskAndBin2(String templateFilePath,
        String templateMaskFilePath, Set<PairInt> outputShape) throws 
        IOException {
        
        ImageProcessor imageProcessor = new ImageProcessor();

        //String fileNameMask0 = fileNames[1];
        //String filePathMask0 = ResourceFinder
        //    .findFileInTestResources(fileNameMask0);
        ImageExt imgMask0 = bin(ImageIOHelper.readImageExt(templateMaskFilePath));

        //String fileName0 = fileNames[0];
        //String filePath0 = ResourceFinder
        //    .findFileInTestResources(fileName0);
        ImageExt img0 = bin(ImageIOHelper.readImageExt(templateFilePath));
    
        ImageExt img0Masked = img0.copyToImageExt();
                
        assert(imgMask0.getNPixels() == img0.getNPixels());

        for (int pixIdx = 0; pixIdx < imgMask0.getNPixels(); ++pixIdx) {
            if (imgMask0.getR(pixIdx) == 0) {
                img0Masked.setRGB(pixIdx, 0, 0, 0);
            } else {
                outputShape.add(new PairInt(imgMask0.getCol(pixIdx), 
                    imgMask0.getRow(pixIdx)));
            }
        }
   
        return new ImageExt[]{img0, img0Masked};
    }
    
    /**
     * @return the templateImage
     */
    public ImageExt[] getTemplateImage() {
        return templateImage;
    }

    /**
     * @return the searchImage
     */
    public ImageExt getSearchImage() {
        return searchImage;
    }
}
