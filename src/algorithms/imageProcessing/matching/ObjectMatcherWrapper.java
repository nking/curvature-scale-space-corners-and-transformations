package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.features.CorrespondenceList;
import algorithms.imageProcessing.features.ObjectMatcher;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
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
    
    private int maxDimension = 256;
    private SIGMA sigma = SIGMA.ZEROPOINTFIVE;//SIGMA.ONE;

    private ImageExt[] templateImage = null;
    private ImageExt searchImage = null;
    
    public List<CorrespondenceList> find(String templateFilePath,
        String templateMaskFilePath,
        String searchFilePath, String debugLabel) throws Exception {
   
        Set<PairInt> shape0 = new HashSet<PairInt>();

        ImageExt[] imgs0 = maskAndBin2(templateFilePath, templateMaskFilePath,
            shape0);
            
        return find(imgs0, shape0, searchFilePath, debugLabel);
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

        //String filePath1 = 
        //    ResourceFinder.findFileInTestResources(searchFilePath);
        
        // template img size is 218  163
        //img = (ImageExt) imageProcessor.bilinearDownSampling(
        //    img, 218, 163, 0, 255);

        long ts = MiscDebug.getCurrentTimeFormatted();

        int w1 = searchImage.getWidth();
        int h1 = searchImage.getHeight();

        /*
        int binFactor1 = (int) Math.ceil(Math.max(
            (float) w1 / maxDimension,
            (float) h1 / maxDimension));
        */
        int binFactor1 = (int) Math.ceil(Math.max(
            (float) w1 / maxDimension,
            (float) h1 / maxDimension));

        ImageProcessor imageProcessor = new ImageProcessor();

        searchImage = imageProcessor.binImage(searchImage, binFactor1);

        /*
        GreyscaleImage theta1 = imageProcessor.createCIELUVTheta(imgs0[0], 255);
        MiscDebug.writeImage(theta1, fileName1Root + "_theta_0");
        theta1 = imageProcessor.createCIELUVTheta(img, 255);
        MiscDebug.writeImage(theta1, fileName1Root + "_theta_1");
        */
        
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
        objMatcher.setToDebug();

        if (debugLabel != null) {
            objMatcher.setToDebug();  
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
    
    public ImageExt[] maskAndBin2(String templateFilePath,
        String templateMaskFilePath, Set<PairInt> outputShape) throws 
        IOException {
        
        ImageProcessor imageProcessor = new ImageProcessor();

        //String fileNameMask0 = fileNames[1];
        //String filePathMask0 = ResourceFinder
        //    .findFileInTestResources(fileNameMask0);
        ImageExt imgMask0 = ImageIOHelper.readImageExt(templateMaskFilePath);

        //String fileName0 = fileNames[0];
        //String filePath0 = ResourceFinder
        //    .findFileInTestResources(fileName0);
        ImageExt img0 = ImageIOHelper.readImageExt(templateFilePath);
    
        int w0 = img0.getWidth();
        int h0 = img0.getHeight();

        int binFactor0 = (int) Math.ceil(Math.max(
             (float) w0 / maxDimension,
             (float) h0 / maxDimension));
        
        if (binFactor0 != 1) {
            img0 = imageProcessor.binImage(img0, binFactor0);
            imgMask0 = imageProcessor.binImage(imgMask0, binFactor0);
        }
        
        ImageExt img0Masked = img0.copyToImageExt();
                
        assert(imgMask0.getNPixels() == img0.getNPixels());

        for (int i = 0; i < imgMask0.getNPixels(); ++i) {
            if (imgMask0.getR(i) == 0) {
                img0Masked.setRGB(i, 0, 0, 0);
            } else {
                outputShape.add(new PairInt(imgMask0.getCol(i), 
                    imgMask0.getRow(i)));
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
