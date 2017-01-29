package algorithms.imageProcessing.features.mser;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.MedianTransform;
import algorithms.imageProcessing.SummedAreaTable;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MSERTest extends TestCase {
    
    public void est0() throws IOException {
        
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
        
        List<TIntObjectMap<CRegion>> cRegions1List 
            = new ArrayList<TIntObjectMap<CRegion>>();
        
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
                TIntObjectMap<CRegion> map = cRegions1List.get(i);
                CRegion r = map.get(specFeatureIdx[i][j]);
                int n1 = r.offsetsToOrigCoords.size();
                
                for (int ii = 0; ii < specFeatureIdx.length; ++ii) {
                    GreyscaleImage mImg2 = mImgs.get(ii);
                    TIntObjectMap<CRegion> map2 = cRegions1List.get(ii);
                    CRegion r2 = map2.get(specFeatureIdx[ii][j]);
                    
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
                        r.ellipseParams.xC, r.ellipseParams.yC, i, ii, 
                        (float)ssdSum, (float)r.autocorrel, (float)f, ssdCount);
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
        
        SummedAreaTable sumTable = new SummedAreaTable();
        int halfDimension = 1;
        
        List<List<GreyscaleImage>> mImgs 
            = new ArrayList<List<GreyscaleImage>>();
        List<GreyscaleImage> mImgs0 = new ArrayList<GreyscaleImage>();
        
        for (String file : files) {
            
            File fl = ResourceFinder.findFileInTmpData(file);
            
            Image img = ImageIOHelper.readImage(fl.getPath());
            GreyscaleImage gsImg = img.copyToGreyscale2();
            mImgs0.add(gsImg);
            
            List<GreyscaleImage> pyr = new ArrayList<GreyscaleImage>();
            MedianTransform mt = new MedianTransform();
                mt.multiscalePyramidalMedianTransform2(gsImg, pyr, 32);
            mImgs.add(pyr);
            
            for (int i = 0; i < pyr.size(); ++i) {
                GreyscaleImage imgM = pyr.get(i);
                imgM = sumTable.createAbsoluteSummedAreaTable(imgM);
                imgM = sumTable.applyMeanOfWindowFromSummedAreaTable(imgM,
                    2 * halfDimension + 1);
                pyr.set(i, imgM);
            }            
        }
        
        // --- build coordinate maps needed for descriptors out of the regions[1] ----
            
        Canonicalizer canonicalizer = new Canonicalizer();
        
        // img, octave, regions
        List<List<TIntObjectMap<CRegion>>> cRegions1List 
            = new ArrayList<List<TIntObjectMap<CRegion>>>();
        
        for (int i = 0; i < mImgs.size(); ++i) {
            
            MSER mser = new MSER();
           
            List<List<Region>> regions = mser.findRegions(mImgs0.get(i));
        
            List<GreyscaleImage> pyr = mImgs.get(i);
            
            List<TIntObjectMap<CRegion>> cregions 
                = canonicalizer.canonicalizeRegions(regions.get(1), pyr);
            
            cRegions1List.add(cregions);
            
            // --- print out the ellipses and centers for images ----
            for (int j = 0; j < pyr.size(); ++j) {
            //for (int j = 0; j < 1; ++j) {
                
                Image img1 = pyr.get(j).copyToColorGreyscale();

                TIntObjectMap<CRegion> crMap = cregions.get(j);
                TIntObjectIterator<CRegion> iter = crMap.iterator();
                
                int nExtraDot = 0;

                for (int ii = 0; ii < crMap.size(); ++ii) {
                    iter.advance();
                    
                    int idx = iter.key();
                    
                    CRegion cr = iter.value();
                     
                    int[] clr = ImageIOHelper.getNextRGB(ii);
                    
                    cr.draw(img1, nExtraDot, clr[0], clr[1], clr[2]);
                }
                
                MiscDebug.writeImage(img1, "_" + i + "_" + j + "_crs_");
            }
        } 

        // compare descriptors
        // may want to filter out size ratios different from "1"
                      
        //for (int imgIdx = 0; imgIdx < cRegions1List.size(); ++imgIdx) {
        for (int imgIdx = 0; imgIdx < 1; ++imgIdx) {
            for (int pyrIdx = 0; pyrIdx < cRegions1List.get(imgIdx).size(); ++pyrIdx) {
                GreyscaleImage mImg = mImgs.get(imgIdx).get(pyrIdx);
                TIntObjectMap<CRegion> cMap = cRegions1List.get(imgIdx).get(pyrIdx);
                
                /* 
                for stereo image matching, would want to store the top few
                matches per octave and only keep the best if the 2nd best
                is larger by 10% or some other error value.
                
                for object finding rather than stereo imaging, would want to
                store the top k matches at each octave where k is the number
                of mser regions in the template object.
                storing the top k should help for objects such as the gbman 
                where there is summetry in the pattern and the best match
                for pose and lighting, etc might be another gbman point
                (false match, but correct object).
                */
                
                for (int imgIdx2 = 0; imgIdx2 < cRegions1List.size(); ++imgIdx2) {
                    for (int pyrIdx2 = 0; pyrIdx2 < cRegions1List.get(imgIdx2).size(); ++pyrIdx2) {
                        GreyscaleImage mImg2 = mImgs.get(imgIdx2).get(pyrIdx2);
                        TIntObjectMap<CRegion> cMap2 = cRegions1List.get(imgIdx2).get(pyrIdx2);
                    
                        TIntObjectIterator<CRegion> iteri1 = cMap.iterator();
                        for (int i1 = 0; i1 < cMap.size(); ++i1) {
                            iteri1.advance();
                            
                            int idx1 = iteri1.key();
                            CRegion cr1 = iteri1.value();
                            int n1 = cr1.offsetsToOrigCoords.size();
                            
                            TIntObjectIterator<CRegion> iteri2 = cMap2.iterator();
                            for (int i2 = 0; i2 < cMap2.size(); ++i2) {
                                iteri2.advance();
                               
                                int idx2 = iteri2.key();
                                CRegion cr2 = iteri2.value();
                                int n2 = cr2.offsetsToOrigCoords.size();
                                
                                int maxMatchable = Math.min(n1, n2);

                                double ssdSum = 0;
                                int ssdCount = 0;

                                for (Entry<PairInt, PairInt> entry
                                    : cr1.offsetsToOrigCoords.entrySet()) {

                                    PairInt pOffsets = entry.getKey();
                                    PairInt xy = entry.getValue();

                                    PairInt xy2 = cr2.offsetsToOrigCoords.get(pOffsets);
                                    if (xy2 == null) {
                                        continue;
                                    }
                                    int v2 = mImg2.getValue(xy2);

                                    int v1 = mImg.getValue(xy);

                                    int diff = v1 - v2;

                                    ssdSum += (diff * diff);
                                    ssdCount++;
                                }
                                if (ssdCount == 0) {
                                    continue;
                                }
                                
                                ssdSum /= (double) ssdCount;
                                ssdSum = Math.sqrt(ssdSum);
                                ssdSum /= 255.;

                                double f = 1. - ((double) ssdCount / (double) maxMatchable);

                                double err = Math.max(cr1.autocorrel, cr2.autocorrel);
                                
                                if (ssdSum <= 1.1 * err) {
                                    String str1 = String.format(
                                        "im1Idx=%d im2Idx=%d pyr1=%d pyr2=%d",
                                        imgIdx, imgIdx2, pyrIdx, pyrIdx2);
                                    String str2 = String.format(
                                        " (%d,%d) (%d,%d) ssd=%.2f autoc=%.2f,%.2f f=%.3f ssd.n=%d",
                                        cr1.ellipseParams.xC, cr1.ellipseParams.yC, 
                                        cr2.ellipseParams.xC, cr2.ellipseParams.yC, 
                                        (float) ssdSum,
                                        (float) cr1.autocorrel, (float) cr2.autocorrel,
                                        (float) f, ssdCount);
                                    System.out.println(str1 + str2);
                                }
                            }
                        }
                        System.out.println("");
                    }
                }
            }
        }  
    }
    
    public void est000() throws IOException {
         
        // write a test similar to test00, but for HSV and normalizes the
        // btightness vector
    }
    
    public void test1() throws IOException {
        
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
        } catch (IOException e) {
            throw new IOException(e);
        } finally {
            if (in != null) {
                in.close();
            }
        }
        
        /*
        MSER mser8(2, 0.0005, 0.1, 0.5, 0.5, true);
        MSER mser4(2, 0.0005, 0.1, 0.5, 0.5, false);
        */
        int delta = 2;
        double minArea = 0.0005; 
        double maxArea = 0.1;
        double maxVariation = 0.5;
        double minDiversity = 0.5;
        
        MSER mser = new MSER();
        
        List<List<Region>> regions = mser.findRegions(
            data, width, height, delta, minArea, maxArea, maxVariation, minDiversity);
            
        Image img2 = ImageIOHelper.readImage(fl.getPath());

        Region.drawEllipses(img2, regions, 0);

        MiscDebug.writeImage(img2, file + "_out_txt_");
    
        int nR = regions.get(0).size() + regions.get(1).size();
        
        System.out.println("num extracted regions=" + nR);
    
        assertEquals(317, nR);
        
    }
    
}
