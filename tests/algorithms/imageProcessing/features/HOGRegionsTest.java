package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import algorithms.imageProcessing.features.mser.Canonicalizer.RegionPoints;
import algorithms.imageProcessing.features.mser.MSER;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import algorithms.misc.MiscMath;
import algorithms.packing.Intersection2DPacking;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class HOGRegionsTest extends TestCase {
    
    public HOGRegionsTest() {
    }

    public void test1() throws IOException {
        
        int maxDimension = 256;
        
        int nCellsPerDim = 3;
        int nPixPerCellDim = 3; 
        int nAngleBins = 9;
        
        String templateFilePath = ResourceFinder.findFileInTestResources(
            "android_statues_01_honeycomb.png");
                
        ImageProcessor imageProcessor = new ImageProcessor();

        ImageExt img0 = ImageIOHelper.readImageExt(templateFilePath);
    
        int w0 = img0.getWidth();
        int h0 = img0.getHeight();

        int binFactor0 = (int) Math.ceil(Math.max(
             (float) w0 / maxDimension,
             (float) h0 / maxDimension));
        
        img0 = imageProcessor.binImage(img0, binFactor0);
        
        GreyscaleImage gsImg = img0.copyToGreyscale2();
        MSER.Threshold thrGs = MSER.Threshold.LEAST_SENSITIVE;
        MSER mser = new MSER();
        List<List<Region>> regions = mser.findRegions(gsImg, thrGs);
         
        Canonicalizer canonicalizer = new Canonicalizer();
        
        TIntObjectMap<Canonicalizer.RegionPoints> regionPoints0 =
            canonicalizer.canonicalizeRegions2(regions.get(1), 
            img0.getWidth(), img0.getHeight());
        
        ObjectMatcher.replaceWithAccumulatedPoints(regionPoints0);
        
        {
            int[] xyCen = new int[2];
            Image im0Cp;
            im0Cp = img0.copyImage();
            int n9 = regions.get(1).size();
            for (int i = 0; i < n9; ++i) {
                Region r = regions.get(1).get(i);
                int[] clr = ImageIOHelper.getNextRGB(i);
                r.drawEllipse(im0Cp, 0, clr[0], clr[1], clr[2]);
                r.calculateXYCentroid(xyCen, im0Cp.getWidth(), im0Cp.getHeight());
                ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], im0Cp,
                    1, 255, 0, 0);
            }
            MiscDebug.writeImage(im0Cp, "_hogs_regions_");
        }
                
        // instead of sobel, using 1st deriv
        GreyscaleImage[] gXgY = 
            imageProcessor.createCentralDifferenceGradients(gsImg);
        GreyscaleImage theta = imageProcessor.computeTheta180(gXgY[0], gXgY[1]);
        GreyscaleImage gXY = 
            imageProcessor.combineConvolvedImages(gXgY[0], gXgY[1]);

        TIntObjectMap<CRegion> regionMap = new TIntObjectHashMap<CRegion>();
        
        HOGRegions hr = new HOGRegions(regionMap, img0.getWidth(), 
            img0.getHeight(), nCellsPerDim, nPixPerCellDim, nAngleBins);
        
        TIntObjectIterator<Canonicalizer.RegionPoints> iter = regionPoints0.iterator();
        for (int i = 0; i < regionPoints0.size(); ++i) {
            iter.advance();
            RegionPoints regionPoints = iter.value();
            hr.addARegion(gXY, theta, regionPoints);
        }
        
        int[] h_0 = new int[hr.getNumberOfBins()];
        int[] h_1 = new int[hr.getNumberOfBins()];
        TIntObjectIterator<CRegion> iter2 = regionMap.iterator();
        for (int i = 0; i < regionMap.size(); ++i) {
            iter2.advance();
            
            CRegion cr0 = iter2.value();
            int rIndex0 = cr0.dataIdx;
            
            //TODO: use rotated image instead
            CRegion cr1 = regionMap.get(i);
            int rIndex1 = cr1.dataIdx;
            
            Intersection2DPacking ip = new Intersection2DPacking();
            Set<PairInt> intersectingKeys = ip.intersection(
                cr0.offsetsToOrigCoords.keySet(),
                cr1.offsetsToOrigCoords.keySet());
            Set<PairInt> offsets0 = ip.naiveStripPacking(
                intersectingKeys, nPixPerCellDim);
            
            int orientation0 = cr0.hogOrientation;
            int orientation1 = cr1.hogOrientation;
            Map<PairInt, PairInt> offsetMap1 = cr1.offsetsToOrigCoords;
            // key = transformed offsets, value = coords in image ref frame,
            // so, can compare dataset0 and dataset1 points with same
            //  keys
            int count = 0;
            for (PairInt pOffset0 : offsets0) {
                PairInt xy1 = offsetMap1.get(pOffset0);
                if (xy1 == null) {
                    continue;
                }
                PairInt xy0 = cr0.offsetsToOrigCoords.get(pOffset0);            
                hr.extractBlock(rIndex0, xy0.getX(), xy0.getY(), h_0);
                hr.extractBlock(rIndex1, xy1.getX(), xy1.getY(), h_1);
            
                float intersection = 
                    hr.intersection(h_0, orientation0, h_1, orientation1);
                
                System.out.println(count + " " + intersection);
                count++;
            }
        }
        
    }
}
