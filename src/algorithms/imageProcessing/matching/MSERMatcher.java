package algorithms.imageProcessing.matching;

import algorithms.FixedSizeSortedVector;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.features.CorrespondenceList;
import algorithms.imageProcessing.features.HOGRegionsManager;
import algorithms.imageProcessing.features.HOGsManager;
import algorithms.imageProcessing.features.ObjectMatcher.Settings;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import algorithms.packing.Intersection2DPacking;
import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
import algorithms.util.QuadInt;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.iterator.TLongIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TLongSet;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class MSERMatcher {

    private boolean debug = false;
    
    private int N_PIX_PER_CELL_DIM = 3;
    private int N_CELLS_PER_BLOCK_DIM = 3;
    private int N_ANGLE_BINS = 9;
    private int N_HIST_BINS = 12;
    
    private static float eps = 0.000001f;

    public void setToDebug() {
        debug = true;
    }
    
    /**
     * This method is a work in progress.  
     * It uses Histogram of Oriented Gradients, histograms of 
     * images of cie luv converted to the polar angle, and
     * histograms of greyscale intensity to find the object in
     * regionPoints0 in the MSER regions of regionPoints1.
     * 
     * The method uses a cell size for the histograms and the results
     * are sensitive to that.
     * The input images have been pre-processed in several ways.
     * The images are binned down with preserved aspect ratios to an image 
     * size such that the largest of width or height is 256 or smaller.
     * 
     * Then ObjectMatcher.findObject12 is used.
     * ObjectMatcher.findObject12 creates the polar theta images and
     * then looks at the general black, white or other characteristics
     * of the template object in dataset0 to determine which MSER 
     * methods should be used (MSER has a positive and negative image 
     * search and several parameters that affect the threshold of the
     * results).
     * The dataset1 MSER regions are filtered to remove those very 
     * different in color than the template object.
     * Both MSER regions are then filtered to keep the strongest mser
     * when there are overlapping mser regions.
     * The results given to this method here are 3 or so mser for dataset0
     * and about 40 or less MSER for dataset1.
     * 
     * The sensitivity of ObjectMatcher.findObject12 and this method to image 
     * resolution and size mean that use of this method should probably be 
     * wrapped in a class that handles resolution and size logic in 
     * pre-processing steps.
     * Note that there may also be some color filter properties that would
     * need to change for extreme cases.
     * 
     * @param pyrRGB0
     * @param pyrPT0
     * @param regionPoints0
     * @param pyrRGB1
     * @param pyrPT1
     * @param regionPoints1
     * @param settings
     * @return 
     */
    public List<CorrespondenceList> matchObject0(
        CMODE clrMode0,
        List<List<GreyscaleImage>> pyrRGB0, List<GreyscaleImage> pyrPT0, 
        TIntObjectMap<Canonicalizer.RegionPoints> regionPoints0, 
        List<List<GreyscaleImage>> pyrRGB1, List<GreyscaleImage> pyrPT1, 
        TIntObjectMap<Canonicalizer.RegionPoints> regionPoints1, 
        Settings settings) {
            
        TIntObjectMap<HOGRegionsManager> hogsMap0 
            = new TIntObjectHashMap<HOGRegionsManager>();
        TIntObjectMap<HOGRegionsManager> hogsMap1 
            = new TIntObjectHashMap<HOGRegionsManager>();

        // a reference to the regions, is same as what the HOGRegions etc hold
        List<TIntObjectMap<CRegion>> cRegionsList0 = createList(pyrRGB0.size());
        List<TIntObjectMap<CRegion>> cRegionsList1 = createList(pyrRGB1.size());
               
        int n0 = pyrPT0.size();
        int n1 = pyrPT1.size();

        int w0 = pyrPT0.get(0).getWidth();
        int h0 = pyrPT0.get(0).getHeight();
        int w1 = pyrPT1.get(0).getWidth();
        int h1 = pyrPT1.get(0).getHeight();
              
        float sizeFactor = 1.2f;//2;//1.2f;

        FixedSizeSortedVector<Obj> bestOverallA =
            //new FixedSizeSortedVector<Obj>(100, Obj.class);
            new FixedSizeSortedVector<Obj>(n0, Obj.class);
           
        for (int pyrIdx0 = 0; pyrIdx0 < n0; ++pyrIdx0) {

            GreyscaleImage gsI0 = combineImages(pyrRGB0.get(pyrIdx0));
            GreyscaleImage ptI0 = pyrPT0.get(pyrIdx0);

            int w0_i = ptI0.getWidth();
            int h0_i = ptI0.getHeight();
            assert(gsI0.getWidth() == w0_i);
            assert(gsI0.getHeight() == h0_i);
            float scale0 = (
                ((float) w0 / (float) w0_i)
                + ((float) h0 / (float) h0_i)) / 2.f;

            TIntObjectMap<CRegion> regions0 = cRegionsList0.get(pyrIdx0);
            
            // instantiate and cache hogs0 if doesn't exist
            HOGRegionsManager hogsMgr0 = getOrCreate(hogsMap0, regions0, 
                w0_i, h0_i, pyrIdx0);
            
            assert(regionPoints0.size() > 0);
            
            hogsMgr0.populateRegionsIfNeeded(
                regionPoints0, scale0, hogsMgr0, gsI0, ptI0);
            
            int nr0 = regions0.size();
            assert(nr0 > 0);
            
            FixedSizeSortedVector<Obj> bestPerOctave =
                //new FixedSizeSortedVector<Obj>(2, Obj.class);
                new FixedSizeSortedVector<Obj>(1, Obj.class);
            
            int maxArea0 = getLargestArea(regions0);
            
            for (int pyrIdx1 = 0; pyrIdx1 < n1; ++pyrIdx1) {

                GreyscaleImage gsI1 = combineImages(pyrRGB1.get(pyrIdx1));
                GreyscaleImage ptI1 = pyrPT1.get(pyrIdx1);

                int w1_i = ptI1.getWidth();
                int h1_i = ptI1.getHeight();
                float scale1 = (((float) w1 / (float) w1_i)
                    + ((float) h1 / (float) h1_i)) / 2.f;

                TIntObjectMap<CRegion> regions1 = cRegionsList1.get(pyrIdx1);
                
                // instantiate and cache hogs0 if doesn't exist
                // instantiate and cache hogs0 if doesn't exist
                HOGRegionsManager hogsMgr1 = getOrCreate(hogsMap1, regions1, 
                    w1_i, h1_i, pyrIdx1);
            
                hogsMgr1.populateRegionsIfNeeded(
                    regionPoints1, scale1, hogsMgr1, gsI1, ptI1);
                
                int nr1 = regions1.size();
                assert(nr1 > 0);
                 
                TIntObjectIterator<CRegion> iter0 = regions0.iterator();
                for (int i0 = 0; i0 < nr0; ++i0) {
                    
                    iter0.advance();
                    
                    int rIdx0 = iter0.key();
                    
                    CRegion cr0 = iter0.value();
      
                    int sz0 = calculateObjectSizeByAvgDist(
                        cr0.ellipseParams.xC, cr0.ellipseParams.yC,
                        cr0.getPixelCoords(), w0_i );

                    if (sz0 == 0) {
                        continue;
                    }
                    
                    /*
                    CMODE cmode0 = CMODE.determineColorMode(
                        pyrRGB0.get(pyrIdx0).get(0),
                        pyrRGB0.get(pyrIdx0).get(1), pyrRGB0.get(pyrIdx0).get(2),
                        cr0.offsetsToOrigCoords.values());
                    CMODE cmodeLUV0 = CMODE.determinePolarThetaMode(ptI0, 
                        cr0.offsetsToOrigCoords.values());
                    System.out.println("cmode0=" + cmode0.name() + ", " + 
                        cmodeLUV0);
                    */
                    
                    //int area0_full = csr0.get(0).get(rIdx0).offsetsToOrigCoords.size();
                    TIntObjectIterator<CRegion> iter1 = regions1.iterator();
                    for (int i1 = 0; i1 < nr1; ++i1) {
                        
                        iter1.advance();
                        
                        // because these regions were made w/ hog orientations,
                        // there may be multiple regions with the same 
                        // original rIdx1 stored as cr.dataIdx, but having
                        // a different rIdx1 here.
                        // so cr.dataIdx is used below for the identity to keep
                        // the best match for the cr.dataIdx (== original rIdx)
                        int rIdx1 = iter1.key();
                        CRegion cr1 = iter1.value();

                        int sz1 = calculateObjectSizeByAvgDist(
                            cr1.ellipseParams.xC, cr1.ellipseParams.yC,
                            cr1.getPixelCoords(), w1_i);

                        if (sz1 == 0) {
                            continue;
                        }
                        
                        /*
                        int xp0 = -1; int yp0 = -1; int xp1 = -1; int yp1 = -1;
                        if (debug) {
                            float scale00 = (((float) w0 / (float) w0_i) + ((float) h0 / (float) h0_i)) / 2.f;
                            float scale01 = (((float) w1 / (float) w1_i) + ((float) h1 / (float) h1_i)) / 2.f;                
                            xp0 = (int)Math.round(scale00 * cr0.ellipseParams.xC);
                            yp0 = (int)Math.round(scale00 * cr0.ellipseParams.yC);
                            xp1 = (int)Math.round(scale01 * cr1.ellipseParams.xC);
                            yp1 = (int)Math.round(scale01 * cr1.ellipseParams.yC);
                        }*/
                        
                        // size filter
                        if ((sz1 > sz0 && ((sz1 / sz0) > sizeFactor))
                            || (sz0 > sz1 && ((sz0 / sz1) > sizeFactor))) {
                            
                            //System.out.format(
                            //    "  REMOVING (%d,%d) where sz0=%d sz1=%d\n",
                            //    xp1, yp1, sz0, sz1);
                            
                            continue;
                        }
                                                   
                        Intersection2DPacking ip = new Intersection2DPacking();
                        Set<PairInt> intersectingOffsetKeys = ip.intersection(
                            cr0.getOffsetKeys(), cr1.getOffsetKeys());
                        Set<PairInt> offsets0 = ip.naiveStripPacking(
                            intersectingOffsetKeys, N_PIX_PER_CELL_DIM);
                        
                        //DEBUG
                        /*{
                            Image tmp = gsI1.copyToColorGreyscale();
                            for (Entry<PairInt, PairInt> entry : 
                                cr1.getOffsetsToOrigCoords().entrySet()) {
                                ImageIOHelper.addPointToImage(entry.getValue().getX(),
                                    entry.getValue().getY(), tmp,
                                    1, 255, 0, 0);
                            };
                            MiscDebug.writeImage(tmp, "_DBG_" 
                                + pyrIdx0 + "_" + i0 + "__" + pyrIdx1 + "_" + i1);
                        }
                        {
                            System.out.println("intersectionKeys.size=" +
                                intersectingOffsetKeys.size() + 
                                " offsets0.size=" + offsets0.size());
                            
                            for (PairInt p : offsets0) {
                                assert(intersectingOffsetKeys.contains(p));
                            }
                        }*/

                        Obj obj = calculateHOGCosts(offsets0, intersectingOffsetKeys.size(), 
                            hogsMgr0, cr0, scale0, rIdx0, pyrIdx0,
                            hogsMgr1, cr1, scale1, rIdx1, pyrIdx1);
                        
                        if (obj == null) {
                            continue;
                        }
                        boolean added = bestPerOctave.add(obj);
                        
                        int x0 = Math.round(scale0 * obj.cr0.ellipseParams.xC);
                        int y0 = Math.round(scale0 * obj.cr0.ellipseParams.yC);
                        int x1 = Math.round(scale1 * obj.cr1.ellipseParams.xC);
                        int y1 = Math.round(scale1 * obj.cr1.ellipseParams.yC);

                        //NOTE: may need to consider the best match
                        //  for each rIdx, that is, consider multiple 
                        //  orientations for a region rather than keeping
                        //  the best orienation for a region only
                        
                        if (debug) {

                            String arrow;
                            if (added) {
                                arrow = "==>";
                            } else {
                                arrow = "   ";
                            }
                            System.out.format(
                            "%s octave %d %d] %s (%d,%d;%d) best: %.4f (%d,%d;%d) [%.3f,%.3f,%.3f,%.3f,%.3f] n=%d\n",
                            settings.getDebugLabel(), pyrIdx0, pyrIdx1,
                            arrow, x1, y1, obj.cr1.hogOrientation,
                            (float) obj.cost, x0, y0, 
                            obj.cr0.hogOrientation,
                            (float) obj.costs[0], (float) obj.costs[1], 
                            (float) obj.costs[2], (float) obj.costs[3], 
                            (float) obj.costs[4], 
                            obj.nMatched
                            );
                        }
                    }
                }
            } // end over dataset1 octaves
            
            // temporarily print the best of each octave0 to look at 
            //    scale biases
            for (int k = 0; k < bestPerOctave.getNumberOfItems(); ++k) {
                Obj obj0 = bestPerOctave.getArray()[k];                               
                bestOverallA.add(obj0);
            }
        }
       
        if (bestOverallA.getNumberOfItems() == 0) {
            return null;
        }
        // re-calculating fraction of whole:
        int maxN = Integer.MIN_VALUE;
        for (int i = 0; i < bestOverallA.getNumberOfItems(); ++i) {
            Obj objB = bestOverallA.getArray()[i];
            if (objB.nMatched > maxN) {
                maxN = objB.nMatched;
            }
        }
        
        Map<QuadInt, Obj> bestCombined = new HashMap<QuadInt, Obj>();        
        for (int i = 0; i < bestOverallA.getNumberOfItems(); ++i) {
            Obj objB = bestOverallA.getArray()[i];
            
            double f = 1. - (objB.nMatched/(double)maxN);
            objB.cost = (float) Math.sqrt(
                objB.costs[0] * objB.costs[0]
                //+ 2. * f * f
                + f * f
                + objB.costs[1] * objB.costs[1]
                + objB.costs[2] * objB.costs[2]
            );         
                        
            int imgIdx0 = objB.imgIdx0;
            int imgIdx1 = objB.imgIdx1;
            GreyscaleImage gsI0 = pyrRGB0.get(imgIdx0).get(1);
            GreyscaleImage gsI1 = pyrRGB1.get(imgIdx1).get(1);            
            int w0_i = gsI0.getWidth();
            int h0_i = gsI0.getHeight();
            float scale0 = (((float)w0/(float)w0_i) + ((float)h0/(float)h0_i))/2.f;
            int w1_i = gsI1.getWidth();
            int h1_i = gsI1.getHeight();
            float scale1 = (((float)w1/(float)w1_i) + ((float)h1/(float)h1_i))/2.f;
            int x0 = Math.round(scale0 * objB.cr0.ellipseParams.xC);
            int y0 = Math.round(scale0 * objB.cr0.ellipseParams.yC);
            int x1 = Math.round(scale1 * objB.cr1.ellipseParams.xC);
            int y1 = Math.round(scale1 * objB.cr1.ellipseParams.yC);
            QuadInt q = new QuadInt(x0, y0, x1, y1);
            Obj existing = bestCombined.get(q);
            if (existing != null) {
                if (existing.cost > objB.cost){
                    bestCombined.put(q, objB);
                }
            } else {
                bestCombined.put(q, objB);
            }
            
            if (debug) {
                //hogCost, f, hcptHgsCost, f0, f1, costHCPT, costHGS
                String lbl = "_" + objB.imgIdx0 + "_" + objB.imgIdx1 + "_"
                    + objB.r0Idx + "_" + objB.r1Idx;
                System.out.format(
                    "_final) %s %d (%d,%d) best: %.4f (%d,%d) %s [%.3f,%.3f,%.3f,%.3f,%.3f] n=%d\n",
                    settings.getDebugLabel(), i, x1, y1,
                    (float) objB.cost, x0, y0, lbl,
                    (float) objB.costs[0], (float) objB.costs[1],
                    (float) objB.costs[2], (float) objB.costs[3],
                    (float) objB.costs[4],
                    objB.nMatched
                );    
            }
        }
        
        List<Obj> bestOverall = new ArrayList<Obj>();
        bestOverall.addAll(bestCombined.values());
        Collections.sort(bestOverall, new CostComparator());
                
        Set<PairInt> pairs = new HashSet<PairInt>();
        List<CorrespondenceList> out = new ArrayList<CorrespondenceList>();
        
        for (int i = 0; i < bestOverall.size(); ++i) {
            
            List<QuadInt> qs = new ArrayList<QuadInt>();
            
            Obj obj = bestOverall.get(i);
            
            PairInt pair = new PairInt(obj.r0Idx, obj.r1Idx);
            if (!pairs.add(pair)) {
                continue;
            }
            
            int imgIdx0 = obj.imgIdx0;
            int imgIdx1 = obj.imgIdx1;
            
            GreyscaleImage gsI0 = pyrRGB0.get(imgIdx0).get(1);
            GreyscaleImage gsI1 = pyrRGB1.get(imgIdx1).get(1);
            
            float scale0, scale1;
            int w0_i = gsI0.getWidth();
            int h0_i = gsI0.getHeight();
            scale0 = (((float)w0/(float)w0_i) + ((float)h0/(float)h0_i))/2.f;

            int w1_i = gsI1.getWidth();
            int h1_i = gsI1.getHeight();
            scale1 = (((float)w1/(float)w1_i) + ((float)h1/(float)h1_i))/2.f;                
            
            // NOTE: for now, just mser centers,
            // but should fill out more than this, including centroid of points
            
            int x0 = Math.round(scale0 * obj.cr0.ellipseParams.xC);
            int y0 = Math.round(scale0 * obj.cr0.ellipseParams.yC);
            int x1 = Math.round(scale1 * obj.cr1.ellipseParams.xC);
            int y1 = Math.round(scale1 * obj.cr1.ellipseParams.yC);
            QuadInt q = new QuadInt(x0, y0, x1, y1);
            qs.add(q);
            
            CorrespondenceList cor = new CorrespondenceList(qs);
            out.add(cor);
            
            if (debug) {
                //hogCost, f, hcptHgsCost, f0, f1, costHCPT, costHGS
                String lbl = "_" + obj.imgIdx0 + "_" + obj.imgIdx1 + "_"
                    + obj.r0Idx + "_" + obj.r1Idx;
                System.out.format(
                    "final) %s %d (%d,%d) best: %.4f (%d,%d) %s [%.3f,%.3f,%.3f,%.3f,%.3f] n=%d\n",
                    settings.getDebugLabel(), i, x1, y1,
                    (float) obj.cost, x0, y0, lbl,
                    (float) obj.costs[0], (float) obj.costs[1],
                    (float) obj.costs[2], (float) obj.costs[3],
                    (float) obj.costs[4],
                    obj.nMatched
                );               
            }
        }
            
        return out;
    }
   
    
    private GreyscaleImage combineImages(List<GreyscaleImage> rgb) {

        GreyscaleImage r = rgb.get(0);
        GreyscaleImage g = rgb.get(1);
        GreyscaleImage b = rgb.get(2);

        if (r.getWidth() != g.getWidth() || r.getWidth() != b.getWidth() ||
            r.getHeight() != g.getHeight() || r.getHeight() != b.getHeight()) {
            throw new IllegalArgumentException("r, g, and b must have same"
                + " width and height");
        }

        int w = r.getWidth();
        int h = r.getHeight();

        GreyscaleImage comb = new GreyscaleImage(w, h, r.getType());

        for (int i = 0; i < r.getNPixels(); ++i) {
            float v0 = r.getValue(i);
            float v1 = g.getValue(i);
            float v2 = b.getValue(i);
            float avg = (v0 + v1 + v2)/3.f;

            comb.setValue(i, Math.round(avg));
        }

        return comb;
    }

    private int calculateObjectSizeByAvgDist(int x, int y, 
        TLongSet pixCoords, int imageWidth) {

        PixelHelper ph = new PixelHelper();
        int[] xy = new int[2];
        
        TLongIterator iter = pixCoords.iterator();
        
        int sumD = 0;
        while (iter.hasNext()) {
            long pix = iter.next();
            ph.toPixelCoords(pix, imageWidth, xy);
            int diffX = xy[0] - x;
            int diffY = xy[1] - y;
            sumD += (diffX * diffX + diffY * diffY);
        }
        sumD = (int)Math.ceil(Math.sqrt((double)sumD/(double)pixCoords.size()));
        
        return sumD;
    }
    
    //double[]{sumHOG, sumHCPT, sumHGS, f0, f1, intersectionCount, area0, area1};
    private double[] sumHOGCost(Set<PairInt> offsets0, int intersectionCount,
        HOGRegionsManager hogMgr0, CRegion cr0, 
        HOGRegionsManager hogMgr1, CRegion cr1) {
                
        if (offsets0.isEmpty()) {
            return null;
        }

        assert(cr0.dataIdx != -1);
        assert(cr1.dataIdx != -1);
        
        int orientation0 = cr0.hogOrientation;
        int orientation1 = cr1.hogOrientation;
        
        Map<PairInt, PairInt> offsetMap1 = cr1.getOffsetsToOrigCoords();

        double sumHOG = 0;
        double sumHCPT = 0;
        double sumHGS = 0;
        
        int[] h0 = new int[N_ANGLE_BINS];
        int[] h1 = new int[h0.length];
        
        int[] ha0 = new int[N_HIST_BINS];
        int[] ha1 = new int[ha0.length];
        float intersection;
        // key = transformed offsets, value = coords in image ref frame,
        // so, can compare dataset0 and dataset1 points with same
        //  keys
    
        int count = 0;
        
        for (PairInt pOffset0 : offsets0) {
            
            PairInt xy1 = offsetMap1.get(pOffset0);

            if (xy1 == null) {
                continue;
            }

            PairInt xy0 = cr0.getOffsetsToOrigCoords().get(pOffset0);
            assert(xy0 != null);
            
            if (!hogMgr0.extractBlockHOG(cr0.dataIdx, xy0.getX(), xy0.getY(), h0)) {
                continue;
            }
            if (!hogMgr1.extractBlockHOG(cr1.dataIdx, xy1.getX(), xy1.getY(), h1)) {
                continue;
            }
            // 1.0 is perfect similarity
            intersection = hogMgr0.intersection(h0, orientation0, 
                h1, orientation1);
            sumHOG += (intersection * intersection);
            
            hogMgr0.extractBlockHCPT(cr0.dataIdx, xy0.getX(), xy0.getY(), ha0);
            hogMgr1.extractBlockHCPT(cr1.dataIdx, xy1.getX(), xy1.getY(), ha1);
            // 1.0 is perfect similarity
            intersection = hogMgr0.intersection(ha0, ha1);
            sumHCPT += (intersection * intersection);
            
            hogMgr0.extractBlockHGS(cr0.dataIdx, xy0.getX(), xy0.getY(), ha0);
            hogMgr1.extractBlockHGS(cr1.dataIdx, xy1.getX(), xy1.getY(), ha1);
            // 1.0 is perfect similarity
            intersection = hogMgr0.intersection(ha0, ha1);
            sumHGS += (intersection * intersection);
            
            count++;
        }
        
        sumHOG /= (double)count;
        sumHOG = Math.sqrt(sumHOG);
        
        sumHCPT /= (double)count;
        sumHCPT = Math.sqrt(sumHCPT);
        
        sumHGS /= (double)count;
        sumHGS = Math.sqrt(sumHGS);
        
        double area1 = cr1.getOffsetsToOrigCoords().size() + eps;
        double f1 = 1. - ((double) intersectionCount / area1);
        
        double area0 = cr0.getOffsetsToOrigCoords().size() + eps;
        double f0 = 1. - ((double) intersectionCount / area0);
        
        return new double[]{sumHOG, sumHCPT, sumHGS, f0, f1, intersectionCount, area0, area1};
    }
    
    private double[] _sumHOGCost(Set<PairInt> offsets0, int intersectionCount,
        HOGsManager _hogMgr0, CRegion cr0, 
        HOGsManager _hogMgr1, CRegion cr1) {
                
        if (offsets0.isEmpty()) {
            return null;
        }

        assert(cr0.dataIdx != -1);
        assert(cr1.dataIdx != -1);
        
        int orientation0 = cr0.hogOrientation;
        int orientation1 = cr1.hogOrientation;
        
        Map<PairInt, PairInt> offsetMap1 = cr1.getOffsetsToOrigCoords();

        double sumHOG = 0;
        double sumHCPT = 0;
        double sumHGS = 0;
        
        int[] h0 = new int[N_ANGLE_BINS];
        int[] h1 = new int[h0.length];
        
        int[] ha0 = new int[N_HIST_BINS];
        int[] ha1 = new int[ha0.length];
        float intersection;
        // key = transformed offsets, value = coords in image ref frame,
        // so, can compare dataset0 and dataset1 points with same
        //  keys
    
        int count = 0;
        
        for (PairInt pOffset0 : offsets0) {
            
            PairInt xy1 = offsetMap1.get(pOffset0);

            if (xy1 == null) {
                continue;
            }

            PairInt xy0 = cr0.getOffsetsToOrigCoords().get(pOffset0);
            assert(xy0 != null);
            
            if (!_hogMgr0.extractBlockHOG(xy0.getX(), xy0.getY(), h0)) {
                continue;
            }
            if (!_hogMgr1.extractBlockHOG(xy1.getX(), xy1.getY(), h1)) {
                continue;
            }
            // 1.0 is perfect similarity
            intersection = _hogMgr0.intersection(h0, orientation0, 
                h1, orientation1);
            sumHOG += (intersection * intersection);
            
            _hogMgr0.extractBlockHCPT(xy0.getX(), xy0.getY(), ha0);
            _hogMgr1.extractBlockHCPT(xy1.getX(), xy1.getY(), ha1);
            // 1.0 is perfect similarity
            intersection = _hogMgr0.intersection(ha0, ha1);
            sumHCPT += (intersection * intersection);
            
            _hogMgr0.extractBlockHGS(xy0.getX(), xy0.getY(), ha0);
            _hogMgr1.extractBlockHGS(xy1.getX(), xy1.getY(), ha1);
            // 1.0 is perfect similarity
            intersection = _hogMgr0.intersection(ha0, ha1);
            sumHGS += (intersection * intersection);
            
            count++;
        }
        
        sumHOG /= (double)count;
        sumHOG = Math.sqrt(sumHOG);
        
        sumHCPT /= (double)count;
        sumHCPT = Math.sqrt(sumHCPT);
        
        sumHGS /= (double)count;
        sumHGS = Math.sqrt(sumHGS);
        
        double area1 = cr1.getOffsetsToOrigCoords().size() + eps;
        double f1 = 1. - ((double) intersectionCount / area1);
        
        double area0 = cr0.getOffsetsToOrigCoords().size() + eps;
        double f0 = 1. - ((double) intersectionCount / area0);
        
        return new double[]{sumHOG, sumHCPT, sumHGS, f0, f1, intersectionCount, 
            area0, area1};
    }
    
    private HOGRegionsManager getOrCreate(
        TIntObjectMap<HOGRegionsManager> hogsMap, 
        TIntObjectMap<CRegion> cRegions, 
        int imageWidth, int imageHeight, int idx) {
        
        HOGRegionsManager hogs = hogsMap.get(idx);
        if (hogs != null) {
            return hogs;
        }
        
        hogs = new HOGRegionsManager(cRegions, 
            imageWidth, imageHeight, 
            N_CELLS_PER_BLOCK_DIM, N_PIX_PER_CELL_DIM, N_ANGLE_BINS, N_HIST_BINS);
        
        hogsMap.put(idx, hogs);
        
        return hogs;
    }

    private int getLargestArea(TIntObjectMap<CRegion> regions) {
        int area = Integer.MIN_VALUE;    
    
        TIntObjectIterator<CRegion> iter = regions.iterator();
        for (int i = 0; i < regions.size(); ++i) {
            iter.advance();
            CRegion cr = iter.value();
            int n = cr.getOffsetsToOrigCoords().size();
            if (n > area) {
                area = n;
            }
        }
        return area;
    }

    private List<TIntObjectMap<CRegion>> createList(int n) {
        List<TIntObjectMap<CRegion>> out = new ArrayList<TIntObjectMap<CRegion>>();
        for (int i = 0; i < n; ++i) {
            out.add(new TIntObjectHashMap<CRegion>());
        }
        return out;
    }

    private Obj calculateHOGCosts(Set<PairInt> offsets0, int intersectionCount, 
        HOGRegionsManager hogsMgr0, CRegion cr0, float scale0, int rIdx0, int pyrIdx0, 
        HOGRegionsManager hogsMgr1, CRegion cr1, float scale1, int rIdx1, int pyrIdx1) {
        
        //double[]{sumHOG, sumHCPT, sumHGS, f0, f1, 
        //   intersectionCount, area0, area1};      
        double[] hogCosts = sumHOGCost(
            offsets0, intersectionCount, 
            hogsMgr0, cr0, hogsMgr1, cr1);
        if (hogCosts == null) {
            return null;
        }
        double hogCost = 1. - hogCosts[0];
        double hcptCost = 1. - hogCosts[1];
        double hgsCost = 1. - hogCosts[2];

        // 1 - fraction of whole (is coverage expressed as a cost)
        double f0 = Math.max(0, hogCosts[3]);
        double f1 = Math.max(0, hogCosts[4]);
        if (f0 > 0.85 || f1 > 0.85) {
            return null;
        }
        double f = (f0 + f1)/2;

        double cost = (float) Math.sqrt(
            hogCost * hogCost
            //+ 2. * f * f
            + f * f
            + hcptCost * hcptCost
            //+ hgsCost * hgsCost
        );

        Obj obj = new Obj();
        obj.cr0 = cr0;
        obj.cr1 = cr1;
        obj.r0Idx = rIdx0;
        obj.r1Idx = rIdx1;
        obj.imgIdx0 = pyrIdx0;
        obj.imgIdx1 = pyrIdx1;
        obj.nMatched = (int) hogCosts[5];
        obj.costs = new double[]{
            hogCost, hcptCost, hgsCost, f0, f1 
        };
        obj.cost = cost;
        obj.f = f;

        return obj;
    }

    private Obj _calculateHOGCosts(Set<PairInt> offsets0, int intersectionCount, 
        HOGsManager _hogsMgr0, CRegion cr0, float scale0, int rIdx0, int pyrIdx0, 
        HOGsManager _hogsMgr1, CRegion cr1, float scale1, int rIdx1, int pyrIdx1) {
        
        //double[]{sumHOG, sumHCPT, sumHGS, f0, f1, 
        //   intersectionCount, area0, area1};      
        double[] hogCosts = _sumHOGCost(
            offsets0, intersectionCount, 
            _hogsMgr0, cr0, _hogsMgr1, cr1);
        if (hogCosts == null) {
            return null;
        }
        double hogCost = 1. - hogCosts[0];
        double hcptCost = 1. - hogCosts[1];
        double hgsCost = 1. - hogCosts[2];

        // 1 - fraction of whole (is coverage expressed as a cost)
        double f0 = Math.max(0, hogCosts[3]);
        double f1 = Math.max(0, hogCosts[4]);
        if (f0 > 0.85 || f1 > 0.85) {
            return null;
        }
        double f = (f0 + f1)/2;

        double cost = (float) Math.sqrt(
            hogCost * hogCost
            //+ 2. * f * f
            + f * f
            + hcptCost * hcptCost
            //+ hgsCost * hgsCost
        );

        Obj obj = new Obj();
        obj.cr0 = cr0;
        obj.cr1 = cr1;
        obj.r0Idx = rIdx0;
        obj.r1Idx = rIdx1;
        obj.imgIdx0 = pyrIdx0;
        obj.imgIdx1 = pyrIdx1;
        obj.nMatched = (int) hogCosts[5];
        obj.costs = new double[]{
            hogCost, hcptCost, hgsCost, f0, f1 
        };
        obj.cost = cost;
        obj.f = f;

        return obj;
    }

    private static class CostComparator implements Comparator<Obj> {
        public CostComparator() {
        }
        @Override
        public int compare(Obj o1, Obj o2) {
            return o1.compareTo(o2);
        }
    }

    /* NOTE: recalculating this worsens the solutions.
    Using the MSER originally determined ellipse parameters has better
    results for the small number of tests here
    private void recalcOrientationAndTrans(HOGs hogs, 
        TIntObjectMap<CRegion> regions, GreyscaleImage gsImg, float scale) {
        
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        
        TIntObjectIterator<CRegion> iter = regions.iterator();
        for (int i = 0; i < regions.size(); ++i) {
            iter.advance();
            int rIdx = iter.key();
            CRegion cr = iter.value();
            
            int orientation = hogs.calculateDominantOrientation(
                cr.offsetsToOrigCoords.values());

            Collection<PairInt> xyp = cr.offsetsToOrigCoords.values();
            
            PairIntArray xy = Misc.convertWithoutOrder(xyp);
            
            PairInt xyCen = ch.calculateXYCentroids2(xyp);
            
            cr.ellipseParams.orientation = Math.PI * orientation/180.;
            cr.ellipseParams.xC = xyCen.getX();
            cr.ellipseParams.yC = xyCen.getY();
            
            Map<PairInt, PairInt> offsetToOrigMap = 
                Canonicalizer.createOffsetToOrigMap(
                xyCen.getX(), xyCen.getY(), 
                xy, gsImg.getWidth(), gsImg.getHeight(), orientation);
        
            cr.offsetsToOrigCoords = offsetToOrigMap;
        }
    }
    */  

    private class Obj implements Comparable<Obj>{
        CRegion cr0;
        CRegion cr1;
        int imgIdx0;
        int imgIdx1;
        int nMatched;
        double cost = Double.MAX_VALUE;
        double[] costs;
        double f;
                
        // might not be populatated:
        int r0Idx = -1;
        int r1Idx = -1;

        @Override
        public int compareTo(Obj other) {
            if (cost < other.cost) {
                return -1;
            } else if (cost > other.cost) {
                return 1;
            }
            if (nMatched > other.nMatched) {
                return -1;
            } else if (nMatched < other.nMatched) {
                return 1;
            }
            // NOTE: may revise this.  wanting to choose smallest scale
            //   or smaller fraction of whole
            if (imgIdx0 < other.imgIdx0 && imgIdx1 < other.imgIdx1) {
                return -1;
            } else if (imgIdx0 > other.imgIdx0 && imgIdx1 > other.imgIdx1) {
                return 1;
            }
            return 0;
        }
    }

    public static int distance(PairInt p1, PairInt p2) {
        int diffX = p1.getX() - p2.getX();
        int diffY = p1.getY() - p2.getY();
        return (int) Math.sqrt(diffX * diffX + diffY * diffY);
    }

    public static int distance(float x0, float y0, float x1, float y1) {
        float diffX = x0 - x1;
        float diffY = y0 - y1;
        return (int) Math.sqrt(diffX * diffX + diffY * diffY);
    }
}
