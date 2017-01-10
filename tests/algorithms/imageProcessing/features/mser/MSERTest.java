package algorithms.imageProcessing.features.mser;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.MedianTransform;
import algorithms.imageProcessing.SummedAreaTable;
import algorithms.imageProcessing.features.ORB;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import algorithms.util.TwoDFloatArray;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;
import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class MSERTest extends TestCase {
    
    public void test0() throws IOException {
        
        String[] files = new String[] {
            "android_statues_01_sz1.jpg",
            "android_statues_02.jpg",
            "android_statues_03_sz3.jpg",
            "android_statues_04_sz1.jpg"
        };
        
        /*
        gbman height for scale reference
           28
           76 
           59
           48
        */
        
        int[][] specFeatureIdx = new int[4][2];
        
        List<List<CRegion>> cRegions1List = new ArrayList<List<CRegion>>();
        
        List<GreyscaleImage> mImgs = new ArrayList<GreyscaleImage>();
        
        Canonicalizer canonicalizer = new Canonicalizer();
        
        int fIdx = 0;
        
        System.out.println("no corrections for scale present yet");
        
        for (String file : files) {
        
            File fl = ResourceFinder.findFileInTmpData(file);
        
            Image img = 
                ImageIOHelper.readImage(fl.getPath());
            
            MSER mser = new MSER();
           
            List<List<Region>> regions = mser.findRegions(img.copyToGreyscale2());
            
            Image img1 = ImageIOHelper.readImage(fl.getPath());
            
            Image img2 = img1.copyImage();
            
            int nExtraDot = 0;
                        
            for (int i = 0; i < regions.get(0).size(); ++i) {
                regions.get(0).get(i).drawEllipse(img1, nExtraDot, 
                    255, 255, 255);
            }

            for (int i = 0; i < regions.get(1).size(); ++i) {
                regions.get(1).get(i).drawEllipse(img2, nExtraDot, 
                    255, 255, 255);
            }
            
            //MiscDebug.writeImage(img1, file + "_out_0_");
            
            MiscDebug.writeImage(img2, file + "_out_1_");
        
            System.out.println("num extracted regions=" +
                (regions.get(0).size() + regions.get(1).size()));
        
            //---------- making descriptors ----
            nExtraDot = 1;
            Canonicalizer cr = new Canonicalizer();

            GreyscaleImage gs1 = img1.copyToGreyscale2();
            int width = gs1.getWidth();
            int height = gs1.getHeight();
            int[] data = new int[width * height];
            for (int i = 0; i < gs1.getNPixels(); ++i) {
                data[i] = gs1.getValue(i);
            }
            
            Image img_1 = img2.copyImage();
            
            //PairIntArray xyCens = cr.extractRegionXYCenters(regions.get(1), 
            //    data, img2.getWidth(), img2.getHeight());

            //ImageIOHelper.addCurveToImage(xyCens, img_1, 
            //    nExtraDot, 255, 0, 0);

            //MiscDebug.writeImage(img_1, file + "_out_xy1_");

            // use a small window for averaging pixels,
            int halfDimension = 1;
            
            // create an image for use with the descriptors
            SummedAreaTable sumTable = new SummedAreaTable();
            GreyscaleImage imgM = sumTable.createAbsoluteSummedAreaTable(
                gs1);
            imgM = sumTable.applyMeanOfWindowFromSummedAreaTable(imgM,
                2 * halfDimension + 1);
       
            mImgs.add(imgM);
            
            //List<Canonicalizer.CRegion> descr0 = 
            //    canonicalizer.canonicalizeRegions(regions.get(0),
            //    img2.getWidth(), img2.getHeight());
    
            List<Canonicalizer.CRegion> descr1 = 
                canonicalizer.canonicalizeRegions(regions.get(1),
                    imgM);
            
            cRegions1List.add(descr1);
        
            
            specFeatureIdx[fIdx] = new int[2];
            
            for (int i = 0; i < descr1.size(); ++i) {
                CRegion r = descr1.get(i);
                
                switch (fIdx) {
                    case 0: {
                        if (95 == r.xC && 59 == r.yC) {
                            specFeatureIdx[fIdx][0] = i;
                        } else if (100 == r.xC && 59 == r.yC) {
                            specFeatureIdx[fIdx][1] = i;
                        }
                        break;
                    }
                    case 1: {
                        if (210 == r.xC && 58 == r.yC) {
                            specFeatureIdx[fIdx][0] = i;
                        } else if (220 == r.xC && 59 == r.yC) {
                            specFeatureIdx[fIdx][1] = i;
                        }
                        break;
                    }
                    case 2: {
                        if (37 == r.xC && 70 == r.yC) {
                            specFeatureIdx[fIdx][0] = i;
                        } else if (53 == r.xC && 70 == r.yC) {
                            specFeatureIdx[fIdx][1] = i;
                        }
                        break;
                    }
                    default: {
                        if (206 == r.xC && 89 == r.yC) {
                            specFeatureIdx[fIdx][0] = i;
                        } else if (220 == r.xC && 88 == r.yC) {
                            specFeatureIdx[fIdx][1] = i;
                        }
                        break;
                    }
                }
            }
            fIdx++;
        } 
        
        /*
        gbman height for scale reference
           28
           76 
           59
           48
        */
        
        // compare the regions in specFeatureIdx
        // before any consideration of scale differences
        for (int j = 0; j < specFeatureIdx[0].length; ++j) {
            for (int i = 0; i < specFeatureIdx.length; ++i) {
                GreyscaleImage mImg = mImgs.get(i);
                List<CRegion> list = cRegions1List.get(i);
                CRegion r = list.get(specFeatureIdx[i][j]);
                int n1 = r.offsetsToOrigCoords.size();
                
                for (int ii = 0; ii < specFeatureIdx.length; ++ii) {
                    GreyscaleImage mImg2 = mImgs.get(ii);
                    List<CRegion> list2 = cRegions1List.get(ii);
                    CRegion r2 = list2.get(specFeatureIdx[ii][j]);
                    
                    int maxMatchable = Math.min(n1, r2.offsetsToOrigCoords.size());
                    
                    double ssdSum = 0;
                    int ssdCount = 0;
                    
                    for (Entry<PairInt, PairInt> entry : 
                        r.offsetsToOrigCoords.entrySet()) {
                        
                        PairInt pOffsets = entry.getKey();
                        PairInt xy = entry.getValue();
                        
                        PairInt xy2 = r2.offsetsToOrigCoords.get(pOffsets);
                        if (xy2 == null) {
                            continue;
                        }
                        int v2 = mImg2.getValue(xy2);

                        int v1 = mImg.getValue(xy);

                        int diff = v1 - v2;

                        ssdSum += (diff * diff);
                        ssdCount++;
                    }
                    ssdSum /= (double)ssdCount;
                    ssdSum = Math.sqrt(ssdSum);
                    ssdSum /= 255.;
                    
                    double f = 1. - ((double)ssdCount/(double)maxMatchable);
                    
                    System.out.format(
  "(%d,%d) i=%d j=%d ssd=%.2f autoc=%.2f f=%.3f ssd.n=%d\n", 
                        r.xC, r.yC, i, ii, (float)ssdSum, 
                        (float)r.autocorrel, (float)f, ssdCount);
                }
            }
        }
    }
    
    public void est00() throws IOException {
        
        String[] files = new String[] {
            "android_statues_01_sz1.jpg",
            "android_statues_02.jpg",
            "android_statues_03_sz3.jpg",
            "android_statues_04_sz1.jpg"
        };
                
        /*
        gbman height for scale reference
           28
           76   2.71 X
           59   2.1  X
           48   1.7  X
        */
        
        // ---- build median windowed greyscale pyramid ----
        
        ImageProcessor imageProcessor = new ImageProcessor();
        SummedAreaTable sumTable = new SummedAreaTable();
        int halfDimension = 1;
        
        List<List<GreyscaleImage>> mImgs 
            = new ArrayList<List<GreyscaleImage>>();
        for (String file : files) {
            File fl = ResourceFinder.findFileInTmpData(file);
            Image img = 
                ImageIOHelper.readImage(fl.getPath());
            List<GreyscaleImage> pyr = imageProcessor.buildPyramid2(
                img.copyToGreyscale2(), 32);
             
            for (int i = 0; i < pyr.size(); ++i) {
                GreyscaleImage imgM = pyr.get(i);
                imgM = sumTable.createAbsoluteSummedAreaTable(imgM);
                imgM = sumTable.applyMeanOfWindowFromSummedAreaTable(imgM,
                    2 * halfDimension + 1);
                pyr.set(i, imgM);
            }            
        }
        
        Canonicalizer canonicalizer = new Canonicalizer();
        
        List<List<List<CRegion>>> cRegions1List 
            = new ArrayList<List<List<CRegion>>>();
                
        for (int i = 0; i < mImgs.size(); ++i) {
            
            List<GreyscaleImage> pyr = mImgs.get(i);
            List<List<CRegion>> list = new ArrayList<List<CRegion>>();
            cRegions1List.add(list);
            
            for (int j = 0; j < pyr.size(); ++j) {
                
                GreyscaleImage mImg = pyr.get(j);
              
                MSER mser = new MSER();
                List<List<Region>> regions = mser.findRegions(
                    mImg.copyImage());
                
                Image img1 = mImg.copyToColorGreyscale();
                Image img2 = img1.copyImage();

                int nExtraDot = 0;

                //for (int jj = 0; jj < regions.get(0).size(); ++jj) {
                //    regions.get(0).get(jj).drawEllipse(img1, nExtraDot, 
                //        255, 255, 255);
                //}

                for (int jj = 0; jj < regions.get(1).size(); ++jj) {
                    regions.get(1).get(jj).drawEllipse(img2, nExtraDot, 
                        255, 255, 255);
                }
            
                //MiscDebug.writeImage(img1, file + "_out_0_");

                MiscDebug.writeImage(img2, "_" + i + "_" + j + "_out_1_");

                System.out.println("num extracted regions=" +
                    (regions.get(0).size() + regions.get(1).size()));

                //---------- making descriptors ----
                nExtraDot = 1;
                Canonicalizer cr = new Canonicalizer();

                GreyscaleImage gs1 = img1.copyToGreyscale2();
                int width = gs1.getWidth();
                int height = gs1.getHeight();
                int[] data = new int[width * height];
                for (int jj = 0; jj < gs1.getNPixels(); ++jj) {
                    data[jj] = gs1.getValue(jj);
                }

                //List<Canonicalizer.CRegion> descr0 = 
                //    canonicalizer.canonicalizeRegions(regions.get(0),
                //    img2.getWidth(), img2.getHeight());

                List<Canonicalizer.CRegion> descr1 = 
                    canonicalizer.canonicalizeRegions(regions.get(1),
                        mImg);

                list.add(descr1);
            }
        } 

        // compare descriptors
        // may want to filter out size rations different from "1"
       
    }
    
    public void est000() throws IOException {
         
        // write a test similar to test00, but for HSV and normalizes the
        // btightness vector
    }
    
    public void est1() throws IOException {
        
        int width = 256;
        int height = 171;
        
        String file = "android_statues_02.jpg";
                    
        File fl = ResourceFinder.findFileInTmpData(file);
        
        int[] data = new int[width * height];
        File flTxt = ResourceFinder.findFileInTmpData("andr02.txt");
        BufferedReader in = null;
        int count = 0;
        try {
            in = new BufferedReader(new FileReader(flTxt));
            String line = in.readLine();
            while (line != null) {
                line = line.trim();
                data[count] = Integer.valueOf(line).intValue();
                count++;
                line = in.readLine();
            }
        } finally {
            if (in != null) {
                in.close();
            }
        }
        
        MSER mser = new MSER();
        
        List<List<Region>> regions = mser.findRegions(data, width, height);
            
        Image img2 = ImageIOHelper.readImage(fl.getPath());

        Region.drawEllipses(img2, regions, 0);

        MiscDebug.writeImage(img2, file + "_out_txt_");
    
        int nR = regions.get(0).size() + regions.get(1).size();
        
        System.out.println("num extracted regions=" + nR);
    
        assertEquals(317, nR);
        
    }
    
}
