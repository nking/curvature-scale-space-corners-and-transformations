package algorithms.imageProcessing.features;

import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.SegmentedCellMerger;
import algorithms.util.PairIntPair;
import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.TestCase;
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

        String fileName1 = "";
        
        List<PairIntPair> pairs = new ArrayList<PairIntPair>();

        for (int i = 0; i < 10; ++i) {
            switch(i) {
                case 0: {
                    fileName1 = "android_statues_01.jpg";
                    pairs.add(new PairIntPair(74,69,73,69));
                    pairs.add(new PairIntPair(74,70,73,70));
                    pairs.add(new PairIntPair(74,72,73,72));
                    pairs.add(new PairIntPair(73,75,72,75));
                    pairs.add(new PairIntPair(69,80,68,80));
                    break;
                }
                case 1: {
                    fileName1 = "android_statues_03.jpg";
                    pairs.add(new PairIntPair(33,65,32,65));
                    pairs.add(new PairIntPair(33,65,33,64));
                    pairs.add(new PairIntPair(51,63,52,63));
                    pairs.add(new PairIntPair(49,62,50,62));
                    pairs.add(new PairIntPair(49,62,29,61));
                    break;
                }
                case 2: {
                    fileName1 = "android_statues_04.jpg";
                    pairs.add(new PairIntPair(189,86,189,85));
                    pairs.add(new PairIntPair(193,85,194,85));
                    break;
                }
                case 3: {
                    fileName1 = "brown_lowe_2003_image1.jpg";
                    pairs.add(new PairIntPair(57,120,57,119));
                    pairs.add(new PairIntPair(79,112,79,111));
                    pairs.add(new PairIntPair(123,95,123,94));
                    pairs.add(new PairIntPair(198,76,199,76));
                    pairs.add(new PairIntPair(199,77,200,76));
                    break;
                }
                case 4: {
                    fileName1 = "brown_lowe_2003_image2.jpg";
                    pairs.add(new PairIntPair(144,113,144,111));
                    pairs.add(new PairIntPair(226,109,226,108));
                    break;
                }
                case 7: {
                    fileName1 = "patagonia_snowy_foreground.jpg";
                    pairs.add(new PairIntPair(176,95,176,94));
                    pairs.add(new PairIntPair(191,92,191,91));
                    pairs.add(new PairIntPair(223,83,223,82));
                    break;
                }
                case 8: {
                    fileName1 = "venturi_mountain_j6_0001.png";
                    pairs.add(new PairIntPair(216,92,216,91));
                    pairs.add(new PairIntPair(224,97,224,96));
                    pairs.add(new PairIntPair(137,93,137,92));
                    pairs.add(new PairIntPair(94,93,94,92));
                    break;
                }
                case 9: {
                    fileName1 = "norwegian_mtn_range.jpg";
                    pairs.add(new PairIntPair(115,107,115,106));
                    pairs.add(new PairIntPair(122,107,122,106));
                    pairs.add(new PairIntPair(123,107,123,106));
                    break;
                }
            }
            writeLabellingFile(fileName1, pairs);
        }
    }

    private void writeLabellingFile(String fileName1,
        List<PairIntPair> coords) throws Exception {

        if (fileName1 == null) {
            return;
        }

        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);

        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);

        ImageProcessor imageProcessor = new ImageProcessor();

        int w1 = img1.getWidth();
        int h1 = img1.getHeight();
        
        int maxDimension = 350;
        int binFactor1 = (int) Math.ceil(Math.max((float)w1/maxDimension,
            (float)h1/ maxDimension));

        ImageExt img1Binned = imageProcessor.binImage(img1, binFactor1);
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        FileWriter writer = null;
        try {
            String outFileName1 = "label_pix_" + fileName1Root + ".csv";
            String outFilePath1 = ResourceFinder.findTmpDataDirectory() + "/" + outFileName1;
            File fl = new File(outFilePath1);
            writer = new FileWriter(fl);
            log.info("writing to " + outFilePath1);
        } catch (IOException ex) {
            Logger.getLogger(SegmentedCellMerger.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        writer.write("#r1,g1,b1,r2,g2,b2,ha1,o1_1,o2_1,o3_1,ha2,o1_2,o2_2,o3_2,labL1,labL2,deltaE1994,deltaE2000,deltaDelta\n");
            
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
            double o2_1 = (r1 + g1 - 2.*b1)/Math.sqrt(6.);
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
            double o2_2 = (r1 + g1 - 2.*b1)/Math.sqrt(6.);
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
            
            writer.write(String.format("%d,%d,%d,%d,%d,%d,",r1,g1,b1,r2,g2,b2));
            
            String str = String.format(
                "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                ha1,(float)o1_1,(float)o2_1,(float)o3_1,
                ha2,(float)o1_2,(float)o2_2,(float)o3_2,
                lab1[0],lab2[0],
                (float)deltaE1994, (float)deltaE2000, (float)deltaDelta
                );
            
            writer.write(str);
        }
        
        writer.flush();
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
