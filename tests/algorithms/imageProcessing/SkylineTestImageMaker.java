package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.PerimeterFinder;
import algorithms.imageProcessing.SkylineExtractor.RemovedSets;
import algorithms.misc.MiscMath;
import algorithms.util.ResourceFinder;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
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
            //"contrail.jpg"
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
    
    public void makeLDAInputFiles() throws Exception {
        
        String[] fileNames = new String[] {
            "brown_lowe_2003_image1.jpg",
            /*//"brown_lowe_2003_image1_rot.jpg",
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
            "sky_with_rainbow2.jpg",*/
            //"30.jpg",
            //"arches_sun_01.jpg",
            //"stlouis_arch.jpg", 
            //"contrail.jpg"
        };
        
        for (String fileName : fileNames) {
                    
            writeLDADataFile(fileName);
            
        }
    }
    
    private void writeLDADataFile(String fileName) throws Exception {
                        
        // revisit infl points.  is there a threshold removing points?

        int idx = fileName.lastIndexOf(".");
        String fileNameRoot = fileName.substring(0, idx);

        String filePath1 = ResourceFinder.findFileInTestResources(
            fileNameRoot + "_sky.png");

        GreyscaleImage skyMask = ImageIOHelper.readImageAsBinary(filePath1);

        // read the contiguous points above 0:
        Set<PairInt> skyPoints = readSkyPixels(skyMask);

        int[] skyRowMinMax = new int[2];
            
        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();

        PerimeterFinder perimeterFinder = new PerimeterFinder();

        int width = skyMask.getWidth();
        int height = skyMask.getHeight();
        int imageMaxColumn = width - 1;
        int imageMaxRow = height - 1;

        Map<Integer, List<PairInt>> skyRowColRanges = perimeterFinder.find(
            skyPoints, skyRowMinMax, imageMaxColumn, 
            outputEmbeddedGapPoints);

        Set<PairInt> borderPixels = perimeterFinder.getBorderPixels(
            skyRowColRanges, skyRowMinMax, imageMaxColumn, imageMaxRow);

        filePath1 = ResourceFinder.findFileInTestResources(fileName);

        ImageExt img = ImageIOHelper.readImageExt(filePath1);
        
        if (img.getWidth() != width || img.getHeight() != height) {
            throw new IllegalStateException(
                "sky mask is not the same size a the test image");
        }
        
        GroupPixelColors allSkyColor = new GroupPixelColors(skyPoints, img, 
            0, 0);
            
        /*
        ((avg contrast of the 24 sky pixel neighbors + border pixel) - 
            (contrast of adjacent non border pixel)
            /(standard deviation of 24 neighbor contrast)
        */
        float borderDiff24Contrast;

        /*
        same, but for the 8 surrounding sky pixel neighbors)
        */
        float borderDiff8Contrast;

        /*
        same, but for the border sky pixel alone)
        */
        float borderDiff1Contrast;

        double rDivB = allSkyColor.getAvgRed() / allSkyColor.getAvgBlue();
        boolean skyIsRed = (rDivB > 1);

        float borderDiff24BlueOrRed;
        float borderDiff8BlueOrRed;
        float borderDiff1BlueOrRed;

        float borderDiff24Hue;
        float borderDiff8Hue;
        float borderDiff1Hue;

        float borderDiff24CIETheta;
        float borderDiff8CIETheta;
        float borderDiff1CIETheta;

        /*
            diffCIEXY = sqrt(diffCIEX*diffCIEX + diffCIEY*diffCIEY)
            theta = 2*arctan(diffCIEY/(diffCIEX + diffCIEXY))
        */
            
        /*
        for simplicity, choosing the adjacent non-sky border pixel to be
        the pixel directly below it if that pixel exists and is not a sky 
        pixel.
        */
        SkylineExtractor skylineExtractor = new SkylineExtractor();

        String dirPath = ResourceFinder.findDirectory("bin");
        String outFilePath = dirPath + "/tmp_lda_border_" + fileNameRoot 
            + ".csv";

        FileWriter fw = null;
        BufferedWriter writer = null;
        
Set<PairInt> inSkyPoints = new HashSet<PairInt>();
Set<PairInt> nonSkyBorderPoints = new HashSet<PairInt>();

        try {
            File file = new File(outFilePath);
            if (file.exists()) {
                file.delete();
            }
            file.createNewFile();

            fw = new FileWriter(file);
            writer = new BufferedWriter(fw);

            for (PairInt p : borderPixels) {
                
                int x = p.getX();
                int y = p.getY();

                int xNonSkyBorder = x;
                int yNonSkyBorder = y + 1;

                if (yNonSkyBorder > (imageMaxRow - 1)) {
                    continue;
                }

                if (skyPoints.contains(new PairInt(xNonSkyBorder, yNonSkyBorder))) {
                    continue;
                }
                
nonSkyBorderPoints.add(new PairInt(xNonSkyBorder, yNonSkyBorder));

                Set<PairInt> neighbors24 = skylineExtractor.getTheNeighborPixels(
                    p, skyPoints, imageMaxColumn, imageMaxRow, 2);
                neighbors24.add(p);

                Set<PairInt> neighbors8 = skylineExtractor.getTheNeighborPixels(
                    p, skyPoints, imageMaxColumn, imageMaxRow, 1);
                neighbors8.add(p);

                float nonSkyBorderLuma = img.getLuma(xNonSkyBorder, yNonSkyBorder);
                float nonSkyBorderRed = img.getR(xNonSkyBorder, yNonSkyBorder);
                float nonSkyBorderBlue = img.getB(xNonSkyBorder, yNonSkyBorder);
                float nonSkyBorderCIEX = img.getCIEX(xNonSkyBorder, yNonSkyBorder);
                float nonSkyBorderCIEY = img.getCIEY(xNonSkyBorder, yNonSkyBorder);
                float nonSkyBorderHue = img.getHue(xNonSkyBorder, yNonSkyBorder);

                GroupPixelColors colors24 = new GroupPixelColors(neighbors24, img, 0, 0);
          
                GroupPixelColors colors8 = new GroupPixelColors(neighbors8,
                    img, 0, 0);

                StringBuilder sb = new StringBuilder();

                borderDiff24Contrast = (colors24.getAvgContrast() - 
                    colors24.calcContrastToOther(nonSkyBorderLuma));
                borderDiff24Contrast /= colors24.getStdDevContrast();

                borderDiff8Contrast = (colors8.getAvgContrast() - 
                    colors8.calcContrastToOther(nonSkyBorderLuma));
                borderDiff8Contrast /= colors8.getStdDevContrast();

                borderDiff1Contrast =
                    (img.getLuma(x, y) - nonSkyBorderLuma)/nonSkyBorderLuma;                    
                borderDiff1Contrast /= allSkyColor.getStdDevContrast();

                if (skyIsRed) {

                    borderDiff24BlueOrRed = (colors24.getAvgRed()
                        - nonSkyBorderRed);
                    borderDiff24BlueOrRed /= colors24.getStdDevRed();

                    borderDiff8BlueOrRed = (colors8.getAvgRed()
                        - nonSkyBorderRed);
                    borderDiff8BlueOrRed /= colors8.getStdDevRed();

                    borderDiff1BlueOrRed = (img.getR(x, y)
                        - nonSkyBorderRed);
                    borderDiff1BlueOrRed /= allSkyColor.getStdDevRed();

                    if (colors24.getStdDevRed() == 0) {
                        continue; 
                    }
                    if (colors8.getStdDevRed() == 0) {
                        continue; 
                    }

                } else {

                    borderDiff24BlueOrRed = (colors24.getAvgBlue()
                        - nonSkyBorderBlue);
                    borderDiff24BlueOrRed /= colors24.getStdDevBlue();

                    borderDiff8BlueOrRed = (colors8.getAvgBlue()
                        - nonSkyBorderBlue);
                    borderDiff8BlueOrRed /= colors8.getStdDevBlue();

                    borderDiff1BlueOrRed = (img.getR(x, y)
                        - nonSkyBorderBlue);
                    borderDiff1BlueOrRed /= allSkyColor.getStdDevBlue();

                    if (colors24.getStdDevBlue() == 0) {
                        continue; 
                    }
                    if (colors8.getStdDevBlue() == 0) {
                        continue; 
                    }
                }
                
                float[] hsb24 = new float[3];
                Color.RGBtoHSB((int)colors24.getAvgRed(), (int)colors24.getAvgGreen(), 
                    (int)colors24.getAvgBlue(), hsb24);
                float[] hsb8 = new float[3];
                Color.RGBtoHSB((int)colors8.getAvgRed(), (int)colors8.getAvgGreen(), 
                    (int)colors8.getAvgBlue(), hsb8);
                float[] hsb1 = new float[3];
                Color.RGBtoHSB(img.getR(x, y), img.getG(x, y), img.getB(x, y), hsb1);
                
                borderDiff24Hue = hsb24[0] - nonSkyBorderHue;
                borderDiff8Hue = hsb8[0] - nonSkyBorderHue;
                borderDiff1Hue = hsb1[0] - nonSkyBorderHue;
                
                /*
                 diffCIEXY = sqrt(diffCIEX*diffCIEX + diffCIEY*diffCIEY)
                 theta = 2*arctan(diffCIEY/(diffCIEX + diffCIEXY))
                 */
                double diffCIEX24 = colors24.getAverageCIEX() - nonSkyBorderCIEX;
                double diffCIEX8 = colors8.getAverageCIEX() - nonSkyBorderCIEX;
                double diffCIEX1 = img.getCIEX(x, y) - nonSkyBorderCIEX;
                double diffCIEY24 = colors24.getAverageCIEY() - nonSkyBorderCIEY;
                double diffCIEY8 = colors8.getAverageCIEY() - nonSkyBorderCIEY;
                double diffCIEY1 = img.getCIEY(x, y) - nonSkyBorderCIEY;

                double diffTheta24 = 2. * Math.atan(diffCIEX24/
                    (diffCIEX24 + diffCIEY24));

                double diffTheta8 = 2. * Math.atan(diffCIEX8/
                    (diffCIEX8 + diffCIEY8));

                double diffTheta1 = 2. * Math.atan(diffCIEX1/
                    (diffCIEX1 + diffCIEY1));

                borderDiff24CIETheta = (float)diffTheta24;
                borderDiff8CIETheta = (float)diffTheta8;
                borderDiff1CIETheta = (float)diffTheta1;

                /*
                sb.append("border,");

                sb.append(Float.toString(borderDiff24Contrast)).append(",")
                    .append(Float.toString(borderDiff8Contrast)).append(",")
                    .append(Float.toString(borderDiff1Contrast)).append(",");
                
                sb.append(Float.toString(borderDiff24Hue)).append(",")
                    .append(Float.toString(borderDiff8Hue)).append(",")
                    .append(Float.toString(borderDiff1Hue)).append(",");

                sb.append(Float.toString(borderDiff24BlueOrRed)).append(",")
                    .append(Float.toString(borderDiff8BlueOrRed)).append(",")
                    .append(Float.toString(borderDiff1BlueOrRed)).append(",");

                sb.append(Float.toString(borderDiff24CIETheta)).append(",")
                    .append(Float.toString(borderDiff8CIETheta)).append(",")
                    .append(Float.toString(borderDiff1CIETheta)).append(",");
                */
sb.append("1,");
sb.append(Float.toString(borderDiff24Contrast));
sb.append(",");
sb.append(Float.toString(borderDiff24Hue));
sb.append(",");
sb.append(Float.toString(borderDiff24BlueOrRed));
//sb.append(",");
//sb.append(Float.toString(borderDiff24CIETheta));

                sb.append("\n");

                writer.write(sb.toString());

                // --------- do same for a comparison to a sky pixel ---
                // for simplicity, looking 3 pixels above current pixel
                // for nearby sky pixel
                int xSkyNearby = x;
                int ySkyNearby = y - 3;

                if (ySkyNearby < 0) {
                    continue;
                }
                if (!skyPoints.contains(new PairInt(xSkyNearby, ySkyNearby))) {
                    continue;
                }
                
inSkyPoints.add(new PairInt(xSkyNearby, ySkyNearby));

                float skyNearbyLuma = img.getLuma(xSkyNearby, ySkyNearby);
                float skyNearbyRed = img.getR(xSkyNearby, ySkyNearby);
                float skyNearbyBlue = img.getB(xSkyNearby, ySkyNearby);
                float skyNearbyCIEX = img.getCIEX(xSkyNearby, ySkyNearby);
                float skyNearbyCIEY = img.getCIEY(xSkyNearby, ySkyNearby);
                float skyNearbyHue = img.getHue(xSkyNearby, ySkyNearby);

                sb = new StringBuilder();

                borderDiff24Contrast = (colors24.getAvgContrast() - 
                    colors24.calcContrastToOther(skyNearbyLuma));
                borderDiff24Contrast /= colors24.getStdDevContrast();

                borderDiff8Contrast = (colors8.getAvgContrast() - 
                    colors8.calcContrastToOther(skyNearbyLuma));
                borderDiff8Contrast /= colors8.getStdDevContrast();

                borderDiff1Contrast =
                    (img.getLuma(x, y) - skyNearbyLuma)/skyNearbyLuma;                    
                borderDiff1Contrast /= allSkyColor.getStdDevContrast();

                if (skyIsRed) {

                    borderDiff24BlueOrRed = (colors24.getAvgRed()
                        - skyNearbyRed);
                    borderDiff24BlueOrRed /= colors24.getStdDevRed();

                    borderDiff8BlueOrRed = (colors8.getAvgRed()
                        - skyNearbyRed);
                    borderDiff8BlueOrRed /= colors8.getStdDevRed();

                    borderDiff1BlueOrRed = (img.getR(x, y) - skyNearbyRed);
                    borderDiff1BlueOrRed /= allSkyColor.getStdDevRed();

                } else {

                    borderDiff24BlueOrRed = (colors24.getAvgBlue()
                        - skyNearbyBlue);
                    borderDiff24BlueOrRed /= colors24.getStdDevBlue();

                    borderDiff8BlueOrRed = (colors8.getAvgBlue()
                        - skyNearbyBlue);
                    borderDiff8BlueOrRed /= colors8.getStdDevBlue();

                    borderDiff1BlueOrRed = (img.getR(x, y) - skyNearbyBlue);
                    borderDiff1BlueOrRed /= allSkyColor.getStdDevBlue();

                }

                borderDiff24Hue = hsb24[0] - skyNearbyHue;
                //borderDiff8Hue = hsb8[0] - skyNearbyHue;
                //borderDiff1Hue = hsb1[0] - skyNearbyHue;
                
                /*
                 diffCIEXY = sqrt(diffCIEX*diffCIEX + diffCIEY*diffCIEY)
                 theta = 2*arctan(diffCIEY/(diffCIEX + diffCIEXY))
                 */
                diffCIEX24 = colors24.getAverageCIEX() - skyNearbyCIEX;
                diffCIEX8 = colors8.getAverageCIEX() - skyNearbyCIEX;
                diffCIEX1 = img.getCIEX(x, y) - skyNearbyCIEX;
                diffCIEY24 = colors24.getAverageCIEY() - skyNearbyCIEY;
                diffCIEY8 = colors8.getAverageCIEY() - skyNearbyCIEY;
                diffCIEY1 = img.getCIEY(x, y) - skyNearbyCIEY;

                diffTheta24 = 2. * Math.atan(diffCIEX24/
                    (diffCIEX24 + diffCIEY24));

                diffTheta8 = 2. * Math.atan(diffCIEX8/
                    (diffCIEX8 + diffCIEY8));

                diffTheta1 = 2. * Math.atan(diffCIEX1/
                    (diffCIEX1 + diffCIEY1));

                borderDiff24CIETheta = (float)diffTheta24;
                borderDiff8CIETheta = (float)diffTheta8;
                borderDiff1CIETheta = (float)diffTheta1;

                /*
                sb.append("sky,");
                sb.append(Float.toString(borderDiff24Contrast)).append(",")
                    .append(Float.toString(borderDiff8Contrast)).append(",")
                    .append(Float.toString(borderDiff1Contrast)).append(",");

                sb.append(Float.toString(borderDiff24Hue)).append(",")
                    .append(Float.toString(borderDiff8Hue)).append(",")
                    .append(Float.toString(borderDiff1Hue)).append(",");

                sb.append(Float.toString(borderDiff24BlueOrRed)).append(",")
                    .append(Float.toString(borderDiff8BlueOrRed)).append(",")
                    .append(Float.toString(borderDiff1BlueOrRed)).append(",");
                sb.append(Float.toString(borderDiff24CIETheta)).append(",")
                    .append(Float.toString(borderDiff8CIETheta)).append(",")
                    .append(Float.toString(borderDiff1CIETheta)).append(",");
                */

sb.append("2,");
sb.append(Float.toString(borderDiff24Contrast));
sb.append(",");
sb.append(Float.toString(borderDiff24Hue));
sb.append(",");
sb.append(Float.toString(borderDiff24BlueOrRed));
//sb.append(",");
//sb.append(Float.toString(borderDiff14CIETheta));

                sb.append("\n");
                
                writer.write(sb.toString());
                
                writer.flush();
            }

            System.out.println("wrote to file " + outFilePath);
            
        } finally {

            if (writer != null) {
                writer.close();
            }
            if (fw != null) {
                fw.close();
            }
        }
       
try {
    ImageExt clr = (ImageExt)img.copyImage();
    ImageIOHelper.addToImage(inSkyPoints,        0, 0, clr, 0, 0, 255);
    ImageIOHelper.addToImage(nonSkyBorderPoints, 0, 0, clr, 255, 0, 0);
    ImageIOHelper.writeOutputImage(
        dirPath + "/tmp_border_" + fileNameRoot + ".png", clr);
} catch (IOException e) {
    System.err.println("ERROR: " + e.getMessage());
}       
        
    }
    
    public static void main(String[] args) {
        
        try {
            SkylineTestImageMaker runner = new SkylineTestImageMaker();

            //runner.makeThresholdedGradientXYImages();
            
            runner.makeLDAInputFiles();
        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }
    }

    private Set<PairInt> readSkyPixels(GreyscaleImage img) {
        
        DFSContiguousValueFinder zerosFinder = new DFSContiguousValueFinder(img);
        
        zerosFinder.findGroupsNotThisValue(0);
        
        int nGroups = zerosFinder.getNumberOfGroups();
        
        // ====== find the group(s) with the largest number of non-zero pixels
        int nMaxGroupN = Integer.MIN_VALUE;
        int[] groupIndexes = new int[nGroups];
        int[] groupN = new int[nGroups];
        for (int gId = 0; gId < nGroups; gId++) {
            int n = zerosFinder.getNumberofGroupMembers(gId);
            groupIndexes[gId] = gId;
            groupN[gId] = n;
            if (n > nMaxGroupN) {
                nMaxGroupN = n;
            }
        }
        
        if (nMaxGroupN > 10000000) {
            MultiArrayMergeSort.sortByDecr(groupN, groupIndexes);
        } else {
            int maxValue = MiscMath.findMax(groupN);
            CountingSort.sortByDecr(groupN, groupIndexes, maxValue);
        }
        
        Set<PairInt> skyPoints = new HashSet<PairInt>();
        
        PairIntArray points = zerosFinder.getXY(groupIndexes[0]);
        
        for (int i = 0; i < points.getN(); i++) {
            skyPoints.add(new PairInt(points.getX(i), points.getY(i)));
        }
        
        return skyPoints;
    }
    
}
