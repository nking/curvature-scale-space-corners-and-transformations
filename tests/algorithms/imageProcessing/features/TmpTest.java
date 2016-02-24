package algorithms.imageProcessing.features;

import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.CannyEdgeFilter;
import algorithms.imageProcessing.EdgeFilterProducts;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.HistogramEqualization;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.ImageSegmentation.BoundingRegions;
import algorithms.imageProcessing.SegmentedCellMerger;
import algorithms.imageProcessing.WaterShed;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.MedianSmooth;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntPair;
import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.TestCase;
import org.ejml.data.Complex64F;
import org.ejml.simple.SimpleEVD;
import org.ejml.simple.SimpleMatrix;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class TmpTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public TmpTest() {
    }

    public void test0() throws Exception {

        String fileName1;

        FeatureMatcherSettings settings = new FeatureMatcherSettings();
        settings.setDebug(true);
        settings.setStartWithBinnedImages(true);
        settings.setToUse2ndDerivCorners();
//trees and grass contributing too many 2nd deriv pts.  does assoc w/ blobs retain enough remaining pts?
        //for (int i = 0; i < 32; ++i) {
        for (int i = 0; i < 1; ++i) {
        //for (int i = 0; i < 20; ++i) {
        //for (int i = 7; i < 8; ++i) {
        //for (int i = 8; i < 9; ++i) {
        //for (int i = 9; i < 10; ++i) {
        //for (int i = 10; i < 11; ++i) {
            switch(i) {
                case 0: {
                    fileName1 = "android_statues_01.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 1: {
                    fileName1 = "android_statues_04.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 3: {
                    fileName1 = "patagonia_snowy_foreground.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 4: {
                    fileName1 = "patagonia_snowy_foreground.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 5: {
                    fileName1 = "mt_rainier_snowy_field.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 6: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 7: {
                    fileName1 = "brown_lowe_2003_image2.jpg";
                    break;
                }
                case 8: {
                    fileName1 = "venturi_mountain_j6_0001.png";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                /*
                case 1: {
                    fileName1 = "android_statues_02.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 2: {
                    fileName1 = "android_statues_03.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 4: {
                    fileName1 = "seattle.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 5: {
                    fileName1 = "stonehenge.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 6: {
                    fileName1 = "cloudy_san_jose.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                
                case 12: {
                    fileName1 = "campus_010.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 13: {
                    fileName1 = "merton_college_I_001.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 14: {
                    fileName1 = "arches.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 15: {
                    fileName1 = "stinson_beach.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }                
                case 16: {
                    fileName1 = "norwegian_mtn_range.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 17: {
                    fileName1 = "halfdome.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 18: {
                    fileName1 = "halfdome2.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 19: {
                    fileName1 = "halfdome3.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 20: {
                    fileName1 = "costa_rica.jpg";
                    fileName2 = fileName1;
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 21: {
                    fileName1 = "new-mexico-sunrise_w725_h490.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 22: {
                    fileName1 = "arizona-sunrise-1342919937GHz.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 23: {
                    fileName1 = "sky_with_rainbow.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 24: {
                    fileName1 = "sky_with_rainbow2.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 25: {
                    fileName1 = "books_illum3_v0_695x555.png";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 26: {
                    fileName1 = "klein_matterhorn_snowy_foreground.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 27: {
                    fileName1 = "30.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 28: {
                    fileName1 = "arches_sun_01.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 29: {
                    fileName1 = "stlouis_arch.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
                case 30: {
                    fileName1 = "contrail.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }*/
                default: {
                    fileName1 = "checkerboard_01.jpg";
                    settings.setUseNormalizedFeatures(true);
                    break;
                }
            }
            writeLabellingFile(fileName1, settings);
        }
    }

    private void writeLabellingFile(String fileName1,
        FeatureMatcherSettings settings) throws Exception {

        if (fileName1 == null) {
            return;
        }

        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);

        settings.setDebugTag(fileName1Root);

        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        ImageProcessor imageProcessor = new ImageProcessor();

        int w1 = img1.getWidth();
        int h1 = img1.getHeight();
        
        int maxDimension = 350;
        int binFactor1 = (int) Math.ceil(Math.max((float)w1/maxDimension,
            (float)h1/ maxDimension));

        ImageExt img1Binned = imageProcessor.binImage(img1, binFactor1);
                
        Set<PairIntPair> coords = new HashSet<PairIntPair>();
        populatePairs(coords, fileName1Root);
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        FileWriter writer = null;
        try {
            String outFileName1 = "label_pix_" + fileName1Root + ".csv";
            String outFilePath1 = ResourceFinder.findTmpDataDirectory() + "/" + outFileName1;
            File fl = new File(outFilePath1);
            writer = new FileWriter(fl);
        } catch (IOException ex) {
            Logger.getLogger(SegmentedCellMerger.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        writer.write("#r1,g1,b1,r2,g2,b2,ha1,o1_1,o2_1,o3_1,ha2,o1_2,o2_2,o3_2,deltaE1994,deltaE2000,deltaDelta");
            
        for (PairIntPair pp : coords) {
            
            int x1 = pp.getX1();
            int y1 = pp.getY1();
            
            int x2 = pp.getX2();
            int y2 = pp.getY2();
            
            int r1 = img1Binned.getR(x1, y1);
            int g1 = img1Binned.getG(x1, y1);
            int b1 = img1Binned.getB(x1, y1);
            float[] lab1 = img1Binned.getCIELAB(x1, y1);
            double o1_1 = (double)(r1 - g1)/Math.sqrt(2.);
            double o2_1 = (double)(r1 + g1 - 2.*b1)/Math.sqrt(6.);
            double o3_1 = (double)(r1 + g1 - b1)/Math.sqrt(2.);
            float ha1;
            if (lab1[1] == 0) {
                ha1 = 0;
            } else {
                ha1 = (float)(Math.atan(lab1[2]/lab1[1]) * 180. / Math.PI);
                if (ha1 < 0) {
                    ha1 += 360.;
                }
            }
            
            int r2 = img1Binned.getR(x2, y2);
            int g2 = img1Binned.getG(x2, y2);
            int b2 = img1Binned.getB(x2, y2);
            float[] lab2 = img1Binned.getCIELAB(x2, y2);
            double o1_2 = (double)(r1 - g1)/Math.sqrt(2.);
            double o2_2 = (double)(r1 + g1 - 2.*b1)/Math.sqrt(6.);
            double o3_2 = (double)(r1 + g1 - b1)/Math.sqrt(2.);
            float ha2;
            if (lab2[1] == 0) {
                ha2 = 0;
            } else {
                ha2 = (float)(Math.atan(lab2[2]/lab2[1]) * 180. / Math.PI);
                if (ha2 < 0) {
                    ha2 += 360.;
                }
            }
            
            double deltaE1994 = Math.abs(cieC.calcDeltaECIE94(lab1, lab2));
            double deltaE2000 = Math.abs(cieC.calcDeltaECIE2000(lab1, lab2));
            double deltaDelta = Math.abs(deltaE1994 - deltaE2000);
            
            writer.write(String.format("%d,%d,%d,%d,%d,%d",r1,g1,b1,r2,g2,b2));
            
            String str = String.format(
                "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                (float)ha1,(float)o1_1,(float)o2_1,(float)o3_1,
                (float)ha2,(float)o1_2,(float)o2_2,(float)o3_2,
                (float)deltaE1994, (float)deltaE2000, (float)deltaDelta
                );
            
            writer.write(str);
        }
        
        try {
            if (writer != null) {
                writer.close();
            }
        } catch (IOException ex) {
            Logger.getLogger(SegmentedCellMerger.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

     public static void main(String[] args) {

        try {
            TmpTest test = new TmpTest();
            test.test0();
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            fail(e.getMessage());
        }
    }

    private void populatePairs(Set<PairIntPair> pairCoords, 
        String fileNameRoot) throws IOException {
        
        BufferedReader bReader = null;
        FileReader reader = null;
        
        String fileName = "label_pix_" + fileNameRoot + "_coords.csv";
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        try {
            reader = new FileReader(new File(filePath));
            
            bReader = new BufferedReader(reader);
            
            //read comment line and discard
            String line = bReader.readLine();
            line = bReader.readLine();
                        
            while (line != null) {
                
                String[] items = line.split(",");
                if (items.length != 5) {
                    throw new IllegalStateException("Error while reading " +
                        fileName + " expecting 5 items in a line");
                }
 
                PairIntPair pp = new PairIntPair(
                    Integer.valueOf(items[0]).intValue(),
                    Integer.valueOf(items[1]).intValue(),
                    Integer.valueOf(items[2]).intValue(),
                    Integer.valueOf(items[3]).intValue());
                
                pairCoords.add(pp);
                
                line = bReader.readLine();                
            }
                        
        } catch (IOException e) {
            log.severe(e.getMessage());
        } finally {
            if (reader == null) {
                reader.close();
            }
            if (bReader == null) {
                bReader.close();
            }
        }        
    }

}
