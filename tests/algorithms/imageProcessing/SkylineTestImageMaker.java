package algorithms.imageProcessing;

import algorithms.imageProcessing.SkylineExtractor.RemovedSets;
import algorithms.util.ResourceFinder;
import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class SkylineTestImageMaker {

    public void makeThresholdedGradientXYImages() throws Exception {
        
        String[] fileNames = new String[] {
            //"brown_lowe_2003_image1.jpg",
            //"brown_lowe_2003_image1_rot.jpg",
            //"brown_lowe_2003_image2.jpg",
            "venturi_mountain_j6_0001.png",
            //"venturi_mountain_j6_0010.png",
            "seattle.jpg",
            "arches.jpg",
            "stinson_beach.jpg",
            "cloudy_san_jose.jpg",            
            "stonehenge.jpg",
            "norwegian_mtn_range.jpg",
            "halfdome.jpg",
            "costa_rica.jpg",
            "new-mexico-sunrise_w725_h490.jpg",
            "arizona-sunrise-1342919937GHz.jpg",
            "sky_with_rainbow.jpg",
            "sky_with_rainbow2.jpg",
            //"30.jpg",
            "arches_sun_01.jpg",
            "stlouis_arch.jpg", 
            "contrail.jpg"
        };
        
        for (String fileName : fileNames) {
                        
            // revisit infl points.  is there a threshold removing points?
            String filePath1 = ResourceFinder.findFileInTestResources(fileName);
            ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
           
            CurvatureScaleSpaceCornerDetector detector = new
                CurvatureScaleSpaceCornerDetector(img1);
            detector.useOutdoorModeAndExtractSkyline();
            detector.findCorners();


            SkylineExtractor skylineExtr = new SkylineExtractor();
            
            int binFactor = skylineExtr.determineBinFactorForSkyMask(
                detector.getTheta().getNPixels());
        
            Set<PairInt> points = new HashSet<PairInt>();
        
            RemovedSets removedSets = skylineExtr.new RemovedSets();
        
            GreyscaleImage threshholdedGXY = skylineExtr.filterAndExtractSkyFromGradient(
                (ImageExt)detector.getOriginalImage(), detector.getTheta(), 
                detector.getGradientXY(), binFactor, points,
                removedSets);

            int idx = fileName.lastIndexOf(".");
            String fileNameRoot = fileName.substring(0, idx);

            String dirPath = ResourceFinder.findDirectory("bin");
            String outFilePath = dirPath + "/tmp_threshholded_sky_" + 
                fileNameRoot + ".png";

            ImageIOHelper.writeOutputImage(outFilePath, threshholdedGXY);
        }
    }
    
    public static void main(String[] args) {
        
        try {
            SkylineTestImageMaker runner = new SkylineTestImageMaker();

            runner.makeThresholdedGradientXYImages();
        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }
    }
    
}
