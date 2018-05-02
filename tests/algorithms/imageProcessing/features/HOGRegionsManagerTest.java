package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import algorithms.imageProcessing.features.mser.MSER;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import algorithms.packing.Intersection2DPacking;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class HOGRegionsManagerTest extends TestCase {
    
    public HOGRegionsManagerTest() {
    }

    public void test1() throws IOException {
        
        int maxDimension = 256;
        
        int N_CELLS_PER_BLOCK_DIM = 3;
        int N_PIX_PER_CELL_DIM = 3; 
        int nAngleBins = 9;
        int nHistBins = 12;
        
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
        
        w0 = img0.getWidth();
        h0 = img0.getHeight();
        
        GreyscaleImage gsImg = img0.copyToGreyscale2();
        MSER.Threshold thrGs = MSER.Threshold.LEAST_SENSITIVE;
        MSER mser = new MSER();
        List<List<Region>> regions = mser.findRegions(gsImg, thrGs);
         
        GreyscaleImage ptImg = imageProcessor
            .createCIELUVTheta_WideRangeLightness(img0, 255);
        
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
        
        TIntObjectMap<CRegion> regionMap = new TIntObjectHashMap<CRegion>();
        
        HOGRegionsManager hogMgr = new HOGRegionsManager(regionMap,
            w0, h0, N_CELLS_PER_BLOCK_DIM, N_PIX_PER_CELL_DIM, nAngleBins, 
            nHistBins);
       
        float scale = 1.0f;
        
        HOGRegionsManager.populateRegionsIfNeeded(
            regionPoints0, scale, hogMgr, gsImg, ptImg);
        
        System.out.println("region map.size=" + regionMap.size());
        
        int[] h_0 = new int[nAngleBins];
        int[] h_1 = new int[nAngleBins];
        int[] ha_0 = new int[nHistBins];
        int[] ha_1 = new int[nHistBins];
        TIntObjectIterator<CRegion> iter2 = regionMap.iterator();
        int n2 = regionMap.size();
        for (int i = 0; i < n2; ++i) {
            iter2.advance();
            
            CRegion cr0 = iter2.value();
            
            int rIndex0 = cr0.dataIdx;
            
            //TODO: use rotated image instead
            CRegion cr1 = regionMap.get(i);
            int rIndex1 = cr1.dataIdx;
            
            
            System.out.println("extract blocks for " + rIndex0 + ", " + rIndex1);
            
            
            Intersection2DPacking ip = new Intersection2DPacking();
            Set<PairInt> intersectingKeys = ip.intersection(
                cr0.offsetsToOrigCoords.keySet(),
                cr1.offsetsToOrigCoords.keySet());
            Set<PairInt> offsets0 = ip.naiveStripPacking(
                intersectingKeys, N_PIX_PER_CELL_DIM);
            
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
                
                if (!hogMgr.extractBlockHOG(rIndex0, xy0.getX(), xy0.getY(), h_0)) {
                    continue;
                }
                if (!hogMgr.extractBlockHOG(rIndex1, xy1.getX(), xy1.getY(), h_1)) {
                    continue;
                }
                float intersection0 = hogMgr.intersection(h_0, orientation0, h_1, orientation1);
                
                hogMgr.extractBlockHCPT(rIndex0, xy0.getX(), xy0.getY(), ha_0);
                hogMgr.extractBlockHCPT(rIndex1, xy1.getX(), xy1.getY(), ha_1);
                float intersection1 = hogMgr.intersection(ha_0, ha_1);
                
                hogMgr.extractBlockHGS(rIndex0, xy0.getX(), xy0.getY(), ha_0);
                hogMgr.extractBlockHGS(rIndex1, xy1.getX(), xy1.getY(), ha_1);
                float intersection2 = hogMgr.intersection(ha_0, ha_1);
                
                System.out.println(count + " " + intersection0 + ", " +
                    intersection1 + ", " + intersection2);
                
                count++;
            }
        }
        
    }
}
