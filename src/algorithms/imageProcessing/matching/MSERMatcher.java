package algorithms.imageProcessing.matching;

import algorithms.compGeometry.PerimeterFinder2;
import algorithms.connected.ConnectedPointsFinder;
import algorithms.FixedSizeSortedVector;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.features.CorrespondenceList;
import algorithms.imageProcessing.features.HCPT;
import algorithms.imageProcessing.features.HGS;
import algorithms.imageProcessing.features.HOGs;
import algorithms.imageProcessing.features.ObjectMatcher.Settings;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import algorithms.imageProcessing.features.mser.Canonicalizer.RegionPoints;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PixelHelper;
import algorithms.util.QuadInt;
import com.climbwithyourfeet.clustering.ClusterFinder;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TLongIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TLongIntHashMap;
import gnu.trove.set.TIntSet;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class MSERMatcher {

    private boolean debug = false;
    
    private int N_PIX_PER_CELL_DIM = 3;
    private int N_CELLS_PER_BLOCK_DIM = 3;

    public void setToDebug() {
        debug = true;
    }

    private TIntObjectMap<CRegion> getOrCreate(
        TIntObjectMap<TIntObjectMap<CRegion>> csrs,
        int imgIdx, GreyscaleImage rgb, float scale) {

        TIntObjectMap<CRegion> csrMap = csrs.get(imgIdx);
        if (csrMap != null) {
            return csrMap;
        }
        csrMap = new TIntObjectHashMap<CRegion>();
        csrs.put(imgIdx, csrMap);

        TIntObjectMap<CRegion> csrMap0 = csrs.get(0);

        int w = rgb.getWidth();
        int h = rgb.getHeight();

        TIntObjectIterator<CRegion> iter = csrMap0.iterator();
        for (int i = 0; i < csrMap0.size(); ++i) {
            iter.advance();
            int idx = iter.key();
            CRegion csr = iter.value();

            if (csr.offsetsToOrigCoords.size() < 9) {
                continue;
            }

            // these are in scale of individual octave (not full reference frame)
            Set<PairInt> scaledSet = extractScaledPts(csr, w, h, scale);

            if (scaledSet.size() < 9) {
                continue;
            }

            PairIntArray xy = new PairIntArray();

            Region r = new Region();
            for (PairInt pl : scaledSet) {
                r.accumulate(pl.getX(), pl.getY());
                xy.add(pl.getX(), pl.getY());
            }

            int[] xyCen = new int[2];
            r.calculateXYCentroid(xyCen, w, h);
            int x = xyCen[0];
            int y = xyCen[1];
            assert(x >= 0 && x < w);
            assert(y >= 0 && y < h);
            double[] m = r.calcParamTransCoeff();

            double angle = Math.atan(m[0]/m[2]);
            if (angle < 0) {
                angle += Math.PI;
            }

            double major = 2. * m[4];
            double minor = 2. * m[5];

            double ecc = Math.sqrt(major * major - minor * minor)/major;
            assert(!Double.isNaN(ecc));

            Map<PairInt, PairInt> offsetMap = Canonicalizer.createOffsetToOrigMap(
                x, y, xy, w, h, angle);

            double autocorrel = Canonicalizer.calcAutoCorrel(rgb, x, y, offsetMap);

            CRegion csRegion = new CRegion();
            csRegion.ellipseParams.orientation = angle;
            csRegion.ellipseParams.eccentricity = ecc;
            csRegion.ellipseParams.major = major;
            csRegion.ellipseParams.minor = minor;
            csRegion.ellipseParams.xC = x;
            csRegion.ellipseParams.yC = y;
            csRegion.offsetsToOrigCoords = offsetMap;
            csRegion.autocorrel = Math.sqrt(autocorrel)/255.;
            csRegion.labels.addAll(csr.labels);
            csRegion.dataIdx = idx;
                
            csrMap.put(idx, csRegion);
        }

        return csrMap;
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

    private Set<PairInt> extractScaledPts(CRegion csr, int w, int h,
        float scale) {

        Set<PairInt> scaledSet = new HashSet<PairInt>();
        for (Entry<PairInt, PairInt> entry : csr.offsetsToOrigCoords.entrySet()) {

            PairInt p = entry.getValue();

            int xScaled = Math.round((float) p.getX() / scale);
            int yScaled = Math.round((float) p.getY() / scale);
            if (xScaled == -1) {
                xScaled = 0;
            }
            if (yScaled == -1) {
                yScaled = 0;
            }
            if (xScaled == w) {
                xScaled = w - 1;
            }
            if (yScaled == h) {
                yScaled = h - 1;
            }
            PairInt pOrigScaled = new PairInt(xScaled, yScaled);

            scaledSet.add(pOrigScaled);
        }

        return scaledSet;
    }

    private Set<PairInt> extractScaledPts(RegionPoints csr, int w, int h,
        float scale) {

        Set<PairInt> scaledSet = new HashSet<PairInt>();
        for (PairInt p : csr.points) {

            int xScaled = Math.round((float) p.getX() / scale);
            int yScaled = Math.round((float) p.getY() / scale);
            if (xScaled == -1) {
                xScaled = 0;
            }
            if (yScaled == -1) {
                yScaled = 0;
            }
            if (xScaled == w) {
                xScaled = w - 1;
            }
            if (yScaled == h) {
                yScaled = h - 1;
            }
            PairInt pOrigScaled = new PairInt(xScaled, yScaled);

            scaledSet.add(pOrigScaled);
        }

        return scaledSet;
    }
   
    private TIntIntMap calculateLabelSizes(List<Set<PairInt>> sets) {
    
        TIntIntMap map = new TIntIntHashMap();
        for (int i = 0; i < sets.size(); ++i) {
            int sz = MiscMath.calculateObjectSize(sets.get(i));
            map.put(i, sz);
        }
        return map;
    }

    
    private int calculateObjectSize(Collection<PairInt> values) {

        int[] minMaxXY = MiscMath.findMinMaxXY(values);
        int diffX = minMaxXY[1] - minMaxXY[0];
        int diffY = minMaxXY[3] - minMaxXY[2];
        double xy = Math.sqrt(diffX * diffX + diffY * diffY);
        
        return (int)Math.round(xy);
    }
    
    private int calculateObjectSizeByAvgDist(int x, int y, 
        Collection<PairInt> values) {

        int sumD = 0;
        for (PairInt p : values) {
            int diffX = p.getX() - x;
            int diffY = p.getY() - y;
            sumD += (diffX * diffX + diffY * diffY);
        }
        sumD = (int)Math.ceil(Math.sqrt((double)sumD/(double)values.size()));
        
        return sumD;
    }

    private void debugPrint3(List<List<Obj>> list, GreyscaleImage img0, 
        GreyscaleImage img1, float scale0, float scale1, String lbl) {
        
        for (int i = 0; i < list.size(); ++i) {
            
            List<Obj> objs = list.get(i);
            
            Image im0 = img0.copyToColorGreyscale();
            Image im1 = img1.copyToColorGreyscale();
            
            for (int j = 0; j < objs.size(); ++j) {
                int[] clr = ImageIOHelper.getNextRGB(j);
                Obj obj = objs.get(j);
                obj.cr0.drawEachPixel(im0, 0, clr[0], clr[1], clr[2]);
                obj.cr1.drawEachPixel(im1, 0, clr[0], clr[1], clr[2]);
            }
            
            MiscDebug.writeImage(im0, lbl + "__0__" + i);
            MiscDebug.writeImage(im1, lbl + "__1__" + i);
        }
    }

    //chordDiffSum, nMatches, p.n
    private double[] partialShapeCost(Obj obj, 
        float scale0, float scale1, 
        List<Set<PairInt>> labeledSets0, 
        List<Set<PairInt>> labeledSets1, 
        GreyscaleImage gsI0, GreyscaleImage gsI1) {
        
        // TODO: use caching if keep this method
        Set<PairInt> set0 = new HashSet<PairInt>();
        TIntIterator iter = obj.cr0.labels.iterator();
        while (iter.hasNext()) {
            int label = iter.next();
            Set<PairInt> s0 = labeledSets0.get(label);
            for (PairInt p : s0) {
                int x = Math.round((float)p.getX() / scale0);
                int y = Math.round((float)p.getY() / scale0);
                set0.add(new PairInt(x, y));
            }
        }
        
        set0 = reduceToContiguous(set0, gsI0.getWidth(), gsI0.getHeight());
        
        Set<PairInt> set1 = new HashSet<PairInt>();
        iter = obj.cr1.labels.iterator();
        while (iter.hasNext()) {
            int label = iter.next();
            Set<PairInt> s0 = labeledSets1.get(label);
            for (PairInt p : s0) {
                int x = Math.round((float)p.getX() / scale1);
                int y = Math.round((float)p.getY() / scale1);
                set1.add(new PairInt(x, y));
            }
        }
        
        set1 = reduceToContiguous(set1, gsI1.getWidth(), gsI1.getHeight());
        
        PerimeterFinder2 finder = new PerimeterFinder2();
        PairIntArray p = finder.extractOrderedBorder(set0);
        PairIntArray q = finder.extractOrderedBorder(set1);
    
        int dp = 1;
        if (p.getN() > 500 || q.getN() > 500) {
            int dn = Math.max(p.getN(), q.getN());
            dp += Math.ceil((float)dn/500.f);
        }
        
        PartialShapeMatcher2 matcher = new PartialShapeMatcher2();
        matcher.overrideSamplingDistance(dp);
        matcher.setToUseEuclidean();
        matcher.setToRemoveOutliers();
        PartialShapeMatcher2.Result result = matcher.match(p, q);

        double[] out = new double[] {
            result.chordDiffSum, result.getNumberOfMatches(), p.getN()
        };

        return out;        
    }

    private Set<PairInt> reduceToContiguous(Set<PairInt> set, int width,
        int height) {
        
        //TODO: convert user of this method to pixels indexes
        ImageProcessor imageProcessor = new ImageProcessor();
        TIntSet pixSet = imageProcessor.convertPointsToIndexes(set, width);
        
        ConnectedPointsFinder finder = new ConnectedPointsFinder(width, height);        
        finder.setMinimumNumberInCluster(1);
        finder.findConnectedPointGroups(pixSet);
        if (finder.getNumberOfGroups() == 1) {
            return new HashSet<PairInt>(set);
        }
        
        int maxIdx = -1;
        int nMax = Integer.MIN_VALUE;
        for (int i = 0; i < finder.getNumberOfGroups(); ++i) {
            int n = finder.getNumberofGroupMembers(i);
            if (n > nMax) {
                nMax = n;
                maxIdx = i;
            }
        }
        return imageProcessor.convertIndexesToPoints(
            finder.getXY(maxIdx), width);
    }
    
    //double[]{intersection, f0, f1, count};
    private double[] sumHOGCost2(HOGs hogs0, CRegion cr0, 
        HOGs hogs1, CRegion cr1) {
        
        int orientation0 = cr0.hogOrientation;
            
        int orientation1 = cr1.hogOrientation;
        
        Map<PairInt, PairInt> offsetMap1 = cr1.offsetsToOrigCoords;

        double sum = 0;
        int count = 0;
        
        int[] h0 = new int[hogs0.getNumberOfBins()];
        int[] h1 = new int[h0.length];
            
        // key = transformed offsets, value = coords in image ref frame,
        // so, can compare dataset0 and dataset1 points with same
        //  keys
        for (Entry<PairInt, PairInt> entry0 : cr0.offsetsToOrigCoords.entrySet()) {

            PairInt pOffset0 = entry0.getKey();

            PairInt xy1 = offsetMap1.get(pOffset0);

            if (xy1 == null) {
                continue;
            }

            PairInt xy0 = entry0.getValue();

            hogs0.extractBlock(xy0.getX(), xy0.getY(), h0);

            hogs1.extractBlock(xy1.getX(), xy1.getY(), h1);

            // 1.0 is perfect similarity
            float intersection = hogs0.intersection(h0, orientation0, 
                h1, orientation1);
            
            sum += (intersection * intersection);

            count++;
        }
        if (count == 0) {
            return null;
        }

        sum /= (double)count;
        sum = Math.sqrt(sum);
        
        double area1 = cr1.offsetsToOrigCoords.size();
        double f1 = 1. - ((double) count / area1);
        
        double area0 = cr0.offsetsToOrigCoords.size();
        double f0 = 1. - ((double) count / area0);
        
        return new double[]{sum, f0, f1, count};
    }

    //double[]{intersection, f0, f1, count};
    private double[] sumHOGCost3(HOGs hogs0, CRegion cr0, 
        HOGs hogs1, CRegion cr1) {
        
        int orientation0 = cr0.hogOrientation;
            
        int orientation1 = cr1.hogOrientation;
        
        Map<PairInt, PairInt> offsetMap1 = cr1.offsetsToOrigCoords;

        double sum = 0;
        double sumErrSq = 0;
        int count = 0;
        
        int[] h0 = new int[hogs0.getNumberOfBins()];
        int[] h1 = new int[h0.length];
            
        float[] diffAndErr;
        
        //NOTE: the coordinate set block could be changed to not visit
        //   every point, but to instead use a spacing of the HOG cell length
        //   which is what Dala & Trigg recommend for their detector windows.
        
        // key = transformed offsets, value = coords in image ref frame,
        // so, can compare dataset0 and dataset1 points with same
        //  keys
        for (Entry<PairInt, PairInt> entry0 : cr0.offsetsToOrigCoords.entrySet()) {

            PairInt pOffset0 = entry0.getKey();

            PairInt xy1 = offsetMap1.get(pOffset0);

            if (xy1 == null) {
                continue;
            }

            PairInt xy0 = entry0.getValue();

            hogs0.extractBlock(xy0.getX(), xy0.getY(), h0);

            hogs1.extractBlock(xy1.getX(), xy1.getY(), h1);

            // 1.0 is perfect similarity
            diffAndErr = hogs0.diff(h0, orientation0, h1, orientation1);
            
            sum += diffAndErr[0];
            sumErrSq += (diffAndErr[1] * diffAndErr[1]);
            
            count++;
        }
        if (count == 0) {
            return null;
        }

        sum /= (double)count;
        
        sumErrSq /= (double)count;
        sumErrSq = Math.sqrt(sumErrSq);
        
        double area1 = cr1.offsetsToOrigCoords.size();
        double f1 = 1. - ((double) count / area1);
        
        double area0 = cr0.offsetsToOrigCoords.size();
        double f0 = 1. - ((double) count / area0);
        
        return new double[]{sum, f0, f1, count, sumErrSq};
    }

    private HOGs getOrCreate(TIntObjectMap<HOGs> hogsMap, GreyscaleImage gs, 
        int idx) {
        
        HOGs hogs = hogsMap.get(idx);
        if (hogs != null) {
            return hogs;
        }
        hogs = new HOGs(gs, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM);
        
        hogsMap.put(idx, hogs);
        
        return hogs;
    }

    private HCPT getOrCreate2(TIntObjectMap<HCPT> hcptMap, GreyscaleImage pt, 
        int idx) {
        
        HCPT hcpt = hcptMap.get(idx);
        if (hcpt != null) {
            return hcpt;
        }
        hcpt = new HCPT(pt, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, 12);
        
        hcptMap.put(idx, hcpt);
        
        return hcpt;
    }

    private HGS getOrCreate3(TIntObjectMap<HGS> hgsMap, GreyscaleImage img, 
        int idx) {
        
        HGS hgs = hgsMap.get(idx);
        if (hgs != null) {
            return hgs;
        }
        hgs = new HGS(img, N_PIX_PER_CELL_DIM, N_CELLS_PER_BLOCK_DIM, 12);
        
        hgsMap.put(idx, hgs);
        
        return hgs;
    }
    
    private void calculateDominantOrientations(
        TIntObjectMap<RegionPoints> regionPoints, HOGs hogs) {

        TIntObjectIterator<RegionPoints> iter = regionPoints.iterator();
        for (int i = 0; i < regionPoints.size(); ++i) {
            iter.advance();
            
            RegionPoints r = iter.value();
            
            TIntSet orientations = hogs.calculateDominantOrientations(r.points);
        
            r.hogOrientations.addAll(orientations);
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
            double diffCost = Math.abs(other.cost - cost);
            if (diffCost < 0.01) {//0.001
                // NOTE: may revise this.  wanting to choose smallest scale
                //   or smaller fraction of whole
                if (imgIdx0 < other.imgIdx0 && imgIdx1 < other.imgIdx1) {
                    return -1;
                } else if (imgIdx0 > other.imgIdx0 && imgIdx1 > other.imgIdx1) {
                    return 1;
                }
                
                if (f < other.f) {
                    return -1;
                } else if (f > other.f) {
                    return 1;
                }
                return 0;
            } else if (cost < other.cost) {
                return -1;
            } else if (cost > other.cost) {
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

    private void debugPrint2(TIntObjectMap<CRegion> cRegions,
        List<GreyscaleImage> rgb, String label) {

        Image img1 = rgb.get(1).copyToColorGreyscale();

        Image img2 = rgb.get(1).copyToColorGreyscale();

        TIntObjectIterator<CRegion> iter = cRegions.iterator();

        int nExtraDot = 0;

        for (int ii = 0; ii < cRegions.size(); ++ii) {
            iter.advance();

            int idx = iter.key();

            CRegion cr = iter.value();

            int[] clr = ImageIOHelper.getNextRGB(ii);

            cr.draw(img1, nExtraDot, clr[0], clr[1], clr[2]);

            cr.drawEachPixel(img2, nExtraDot, clr[0], clr[1], clr[2]);
        }

        MiscDebug.writeImage(img1, label + "_" + "_csrs_");

        MiscDebug.writeImage(img2, label + "_" + "_csrs_pix_");

        System.out.println(cRegions.size() + " labeled regions for " + label);
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
        List<List<GreyscaleImage>> pyrRGB0, List<GreyscaleImage> pyrPT0, 
        TIntObjectMap<Canonicalizer.RegionPoints> regionPoints0, 
        List<List<GreyscaleImage>> pyrRGB1, List<GreyscaleImage> pyrPT1, 
        TIntObjectMap<Canonicalizer.RegionPoints> regionPoints1, 
        Settings settings) {
            
        TIntObjectMap<HOGs> hogsMap0 = new TIntObjectHashMap<HOGs>();
        TIntObjectMap<HOGs> hogsMap1 = new TIntObjectHashMap<HOGs>();
        
        TIntObjectMap<HCPT> hcptMap0 = new TIntObjectHashMap<HCPT>();
        TIntObjectMap<HGS> hgsMap0 = new TIntObjectHashMap<HGS>();
        TIntObjectMap<HCPT> hcptMap1 = new TIntObjectHashMap<HCPT>();
        TIntObjectMap<HGS> hgsMap1 = new TIntObjectHashMap<HGS>();

        // use hogs to calculate the dominant orientations
        calculateDominantOrientations(regionPoints0, 
            getOrCreate(hogsMap0, combineImages(pyrRGB0.get(0)), 0));
        
        calculateDominantOrientations(regionPoints1, 
            getOrCreate(hogsMap1, combineImages(pyrRGB1.get(0)), 0));
        
        Canonicalizer canonicalizer = new Canonicalizer();
        
        // create the CRegion objects which have the rotated points and 
        //    offsets in them
        TIntObjectMap<CRegion> cRegions0 = canonicalizer.canonicalizeRegions4(
            regionPoints0, pyrRGB0.get(0).get(1));
        
        TIntObjectMap<CRegion> cRegions1 = canonicalizer.canonicalizeRegions4(
            regionPoints1, pyrRGB1.get(0).get(1));
                
        // populated on demand, some are skipped for large size differences
        TIntObjectMap<TIntObjectMap<CRegion>> csr0
            = new TIntObjectHashMap<TIntObjectMap<CRegion>>();
        csr0.put(0, cRegions0);

        TIntObjectMap<TIntObjectMap<CRegion>> csr1
            = new TIntObjectHashMap<TIntObjectMap<CRegion>>();
        csr1.put(0, cRegions1);

        // key = region index, value = Obj w/ cost being hog intersection
        TIntObjectMap<FixedSizeSortedVector<Obj>> rIndexHOGMap1
            = new TIntObjectHashMap<FixedSizeSortedVector<Obj>>();
        
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
            float scale0 = (((float) w0 / (float) w0_i)
                + ((float) h0 / (float) h0_i)) / 2.f;

            HOGs hogs0 = getOrCreate(hogsMap0, gsI0, pyrIdx0);

            HCPT hcpt0 = getOrCreate2(hcptMap0, ptI0, pyrIdx0);
            HGS hgs0 = getOrCreate3(hgsMap0, gsI0, pyrIdx0);
            
            TIntObjectMap<CRegion> regions0 = getOrCreate(csr0, pyrIdx0, gsI0,
                scale0);
 
            FixedSizeSortedVector<Obj> bestPerOctave =
                //new FixedSizeSortedVector<Obj>(100, Obj.class);
                new FixedSizeSortedVector<Obj>(1, Obj.class);
            
            for (int pyrIdx1 = 0; pyrIdx1 < n1; ++pyrIdx1) {

                GreyscaleImage gsI1 = combineImages(pyrRGB1.get(pyrIdx1));
                GreyscaleImage ptI1 = pyrPT1.get(pyrIdx1);

                int w1_i = ptI1.getWidth();
                int h1_i = ptI1.getHeight();
                float scale1 = (((float) w1 / (float) w1_i)
                    + ((float) h1 / (float) h1_i)) / 2.f;

                TIntObjectMap<CRegion> regions1 = getOrCreate(csr1, pyrIdx1,
                    gsI1, scale1);

                HOGs hogs1 = getOrCreate(hogsMap1, gsI1, pyrIdx1);
                HCPT hcpt1 = getOrCreate2(hcptMap1, ptI1, pyrIdx1);
                HGS hgs1 = getOrCreate3(hgsMap1, gsI1, pyrIdx1);
                
                TIntObjectIterator<CRegion> iter0 = regions0.iterator();
                for (int i0 = 0; i0 < regions0.size(); ++i0) {
                    
                    iter0.advance();
                    
                    int rIdx0 = iter0.key();
                    
                    CRegion cr0 = iter0.value();
      
                    int sz0 = calculateObjectSizeByAvgDist(
                        cr0.ellipseParams.xC, cr0.ellipseParams.yC,
                        cr0.offsetsToOrigCoords.values());
                    
                    //int area0_full = csr0.get(0).get(rIdx0).offsetsToOrigCoords.size();
                    TIntObjectIterator<CRegion> iter1 = regions1.iterator();
                    for (int i1 = 0; i1 < regions1.size(); ++i1) {
                        
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
                            cr1.offsetsToOrigCoords.values());

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

                        double[] hogCosts;
                        float hogCost;
                        
                        if (false) {// use intersection    
                            //double[]{intersection, f0, f1, count};
                            hogCosts = sumHOGCost2(hogs0, cr0, hogs1, cr1);
                            if (hogCosts == null) {
                                continue;
                            }
                            hogCost = 1.f - (float) hogCosts[0];
                        } else {// use mean difference    
                            //double[]{diff, f0, f1, count, error};
                            hogCosts = sumHOGCost3(hogs0, cr0, hogs1, cr1);
                            if (hogCosts == null) {
                                continue;
                            }
                            hogCost = (float)hogCosts[0];
                        }
                        
                        // 1 - fraction of whole (is coverage expressed as a cost)
                        double f0 = Math.max(0, hogCosts[1]);
                        double f1 = Math.max(0, hogCosts[2]);
                        double f = (f0 + f1)/2;
                        
                        Obj obj = new Obj();
                        obj.cr0 = cr0;
                        obj.cr1 = cr1;
                        obj.r0Idx = rIdx0;
                        obj.r1Idx = rIdx1;
                        obj.imgIdx0 = pyrIdx0;
                        obj.imgIdx1 = pyrIdx1;
                        obj.nMatched = (int) hogCosts[3];
                        
                        double cost;
                        double[] costs2;
                        double hcptHgsCost;
                        double hcptCost;
                        double hgsCost;
                       
                        if (true) {
                            //double[]{combIntersection, f0, f1, 
                            //    intersectionHCPT, intersectionHGS, count}
                            costs2 = sumCost2(hcpt0, hgs0, cr0, hcpt1, hgs1, cr1);
                            hcptHgsCost = 1.f - costs2[3];//costs2[0];
                            hcptCost = 1.f - costs2[3];
                            hgsCost = 1.f - costs2[4];
                        } else {
                            costs2 = sumCost3(hcpt0, hgs0, cr0, hcpt1, hgs1, cr1);
                            hcptHgsCost = costs2[3];//costs2[0];
                            hcptCost = costs2[3];
                            hgsCost = costs2[4];
                        }
                        
                        cost = (float) Math.sqrt(
                            2. * hogCost * hogCost
                            + 2. * f * f
                            + hcptHgsCost * hcptHgsCost
                        );
                        obj.costs = new double[]{
                            hogCost, f, hcptHgsCost, f0, f1, hcptCost, hgsCost
                        };
                        
                        obj.cost = cost;
                        obj.f = f;

                        int x0 = Math.round(scale0 * obj.cr0.ellipseParams.xC);
                        int y0 = Math.round(scale0 * obj.cr0.ellipseParams.yC);
                        int x1 = Math.round(scale1 * obj.cr1.ellipseParams.xC);
                        int y1 = Math.round(scale1 * obj.cr1.ellipseParams.yC);
                        
                        //NOTE: may need to consider the best match
                        //  for each rIdx, that is, consider multiple 
                        //  orientations for a region rather than keeping
                        //  the best orienation for a region only
                        
                        // add to r1 map (which is actually cr.dataIdx)
                        FixedSizeSortedVector<Obj> objVec = rIndexHOGMap1.get(
                            cr1.dataIdx);
                        if (objVec == null) {
                            objVec = new FixedSizeSortedVector<Obj>(
                                n0, Obj.class);
                            rIndexHOGMap1.put(cr1.dataIdx, objVec);
                        }
                        boolean added = objVec.add(obj);
                        
                        added = bestPerOctave.add(obj);
                        
                        /*if (debug) {
                            double cost2 = (float) Math.sqrt(
                                obj.costs[0]*obj.costs[0] +
                                obj.costs[1]*obj.costs[1] +
                                obj.costs[2]*obj.costs[2]
                            );
                            System.out.format("%s octave %d %d] (%d,%d) best: %.4f (%d,%d) [%.3f,%.3f,%.3f] n=%d c2=%.3f\n",
                            settings.getDebugLabel(), pyrIdx0, pyrIdx1,
                            x1, y1, (float) obj.cost,
                            Math.round(scale0 * obj.cr0.ellipseParams.xC),
                            Math.round(scale0 * obj.cr0.ellipseParams.yC), 
                            (float) obj.costs[0], (float) obj.costs[1], (float) obj.costs[2], 
                            obj.cr0.offsetsToOrigCoords.size(),
                            (float)cost2
                            );
                        }*/
                    }
                }
            } // end over dataset1 octaves
            
            // temporarily print the best of each octave0 to look at 
            //    scale biases
            for (int k = 0; k < bestPerOctave.getNumberOfItems(); ++k) {
                Obj obj0 = bestPerOctave.getArray()[k];
                int imgIdx1 = obj0.imgIdx1;
                GreyscaleImage gsI1 = pyrRGB1.get(imgIdx1).get(1);
                int w1_i = gsI1.getWidth();
                int h1_i = gsI1.getHeight();
                float scale1 = (((float) w1 / (float) w1_i) + ((float) h1 / (float) h1_i)) / 2.f;
                int or0 = obj0.cr0.hogOrientation;
                int or1 = obj0.cr1.hogOrientation;
                //String str1 = String.format("angles=(%d,%d) s=(%.1f,%.1f)", 
                //    or0, or1, scale0, scale1);

                double cost2 = (float) Math.sqrt(
                    obj0.costs[0]*obj0.costs[0] +
                    obj0.costs[1]*obj0.costs[1] +
                    obj0.costs[2]*obj0.costs[2]
                );
                                
                bestOverallA.add(obj0);
            }
        }
         
        System.out.println("r1 points size = " + regionPoints1.size()
            + " r1 map size filtered = " + rIndexHOGMap1.size());

        if (debug) {
        //if (true) {
            
            Map<QuadInt, Obj> bestMap = new HashMap<QuadInt, Obj>();
            
            // re-ordering the best for each rIdx1:
            FixedSizeSortedVector<Obj> tmp1 = new FixedSizeSortedVector<Obj>(
                //rIndexHOGMap.size(), 
                5, Obj.class);

            // --- print out phog based rankings -----
            // printing range of hog values for a region1
            TIntObjectIterator<FixedSizeSortedVector<Obj>> iter2 = rIndexHOGMap1.iterator();

            StringBuilder sb = new StringBuilder();

            for (int i3 = 0; i3 < rIndexHOGMap1.size(); ++i3) {

                iter2.advance();

                //int rIdx = iter2.key();
                FixedSizeSortedVector<Obj> vec = iter2.value();

                int n = vec.getNumberOfItems();
                if (n == 0) {
                    continue;
                }

                Obj obj0 = vec.getArray()[0];

                tmp1.add(obj0);
            }

            System.out.println(sb.toString());

            StringBuilder sb2 = new StringBuilder();

            for (int i = 0; i < tmp1.getNumberOfItems(); ++i) {

                Obj obj0 = tmp1.getArray()[i];

                int imgIdx0 = obj0.imgIdx0;
                int imgIdx1 = obj0.imgIdx1;

                GreyscaleImage gsI0 = pyrRGB0.get(imgIdx0).get(1);
                GreyscaleImage gsI1 = pyrRGB1.get(imgIdx1).get(1);

                float scale00, scale01;
                {
                    int w0_i = gsI0.getWidth();
                    int h0_i = gsI0.getHeight();
                    scale00 = (((float) w0 / (float) w0_i) + ((float) h0 / (float) h0_i)) / 2.f;

                    int w1_i = gsI1.getWidth();
                    int h1_i = gsI1.getHeight();
                    scale01 = (((float) w1 / (float) w1_i) + ((float) h1 / (float) h1_i)) / 2.f;
                }

                String lbl = "_" + obj0.imgIdx0 + "_" + obj0.imgIdx1 + "_"
                    + obj0.r0Idx + "_" + obj0.r1Idx;

                int or0 = (int) Math.round(
                    obj0.cr0.ellipseParams.orientation * 180. / Math.PI);

                int or1 = (int) Math.round(
                    obj0.cr1.ellipseParams.orientation * 180. / Math.PI);

                String str1 = String.format("angles=(%d,%d ; %d,%d)",
                    or0, or1, obj0.cr0.hogOrientation, 
                    obj0.cr1.hogOrientation);

                int x1 = Math.round(scale01 * obj0.cr1.ellipseParams.xC);
                int y1 = Math.round(scale01 * obj0.cr1.ellipseParams.yC);
                int x0 = Math.round(scale00 * obj0.cr0.ellipseParams.xC);
                int y0 = Math.round(scale00 * obj0.cr0.ellipseParams.yC);
                
                //hogCost, f, hcptHgsCost, f0, f1, costHCPT, costHGS
                
                sb2.append(String.format(
 "2] r1 %s %d (%d,%d) best: %.4f (%d,%d) %s [%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f] %s n=%d\n",
                    settings.getDebugLabel(), i, 
                    x1, y1,
                    (float) obj0.cost, x0, y0, lbl,
                    (float) obj0.costs[0], (float) obj0.costs[1], 
                    (float) obj0.costs[2], (float) obj0.costs[3], 
                    (float) obj0.costs[4], (float) obj0.costs[5], (float) obj0.costs[6],
                    str1, obj0.cr0.offsetsToOrigCoords.size()
                ));
                
                QuadInt q = new QuadInt(x0, y0, x1, y1);
                Obj obj = bestMap.get(q);
                if (obj != null) {
                    if (obj0.compareTo(obj) < 0) {
                        bestMap.put(q, obj0);
                    }
                } else {
                    bestMap.put(q, obj0);
                }
                
                /*
                Image im0 = gsI0.copyToColorGreyscale();
                Image im1 = gsI1.copyToColorGreyscale();
                int[] clr = new int[]{255,0,0};
                obj0.cr0.drawEachPixel(im0, 0, clr[0], clr[1], clr[2]);
                obj0.cr1.drawEachPixel(im1, 0, clr[0], clr[1], clr[2]);
                //obj0.cr1.draw(im1, 1, 0, 0, 0);
                MiscDebug.writeImage(im0, debugLabel + "_" + lbl);
                MiscDebug.writeImage(im1, debugLabel + "_" + lbl);
                */
            }
            System.out.println(sb2.toString());
            
        //    bestOverallA = new FixedSizeSortedVector<Obj>(5, Obj.class);
        //    for (Entry<QuadInt, Obj> entry : bestMap.entrySet()) {
        //        bestOverallA.add(entry.getValue());
        //    }
        }
        
        if (bestOverallA.getNumberOfItems() == 0) {
            return null;
        }
        
        /*
        TODO: revisit with more tests.
        a few tests suggest that the correct answer is to order by cost,
        then walk down the array if hogs cost is lower for 2nd best
        */
        if (debug) {
            System.out.println("looking for smallest hogs cost:");
        }
        
        List<Obj> bestOverall = new ArrayList<Obj>(bestOverallA.getNumberOfItems());        
        for (int i = 0; i < bestOverallA.getNumberOfItems(); ++i) {
            Obj objB = bestOverallA.getArray()[i];
            bestOverall.add(objB);
        }
        
        Set<PairInt> pairs = new HashSet<PairInt>();
        
        // storing top 5 of r1 matches
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
            {
                int w0_i = gsI0.getWidth();
                int h0_i = gsI0.getHeight();
                scale0 = (((float)w0/(float)w0_i) + ((float)h0/(float)h0_i))/2.f;
                
                int w1_i = gsI1.getWidth();
                int h1_i = gsI1.getHeight();
                scale1 = (((float)w1/(float)w1_i) + ((float)h1/(float)h1_i))/2.f;                
            }
            
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
                    "final) %s %d (%d,%d) best: %.4f (%d,%d) %s [%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f] n=%d\n",
                    settings.getDebugLabel(), i, x1, y1,
                    (float) obj.cost, x0, y0, lbl,
                    (float) obj.costs[0], (float) obj.costs[1],
                    (float) obj.costs[2], (float) obj.costs[3],
                    (float) obj.costs[4], (float) obj.costs[5], (float) obj.costs[6],
                    obj.cr0.offsetsToOrigCoords.size()
                );               
            }
        }

        return out;
    }
    
    private void _debugPrint(Obj obj, GreyscaleImage gsI0, GreyscaleImage gsI1,
        int w0, int h0, int w1, int h1, String debugLabel) {
        
        if (debug) {
                
            StringBuilder sb2 = new StringBuilder();

            float scale00, scale01;
            {
                int w0_i = gsI0.getWidth();
                int h0_i = gsI0.getHeight();
                scale00 = (((float) w0 / (float) w0_i) + ((float) h0 / (float) h0_i)) / 2.f;

                int w1_i = gsI1.getWidth();
                int h1_i = gsI1.getHeight();
                scale01 = (((float) w1 / (float) w1_i) + ((float) h1 / (float) h1_i)) / 2.f;
            }

            String lbl = "_" + obj.imgIdx0 + "_" + obj.imgIdx1 + "_"
                    + obj.r0Idx + "_" + obj.r1Idx;

            int or0 = (int) Math.round(
                    obj.cr0.ellipseParams.orientation * 180. / Math.PI);

            int or1 = (int) Math.round(
                    obj.cr1.ellipseParams.orientation * 180. / Math.PI);

            String str1 = String.format("angles=(%d,%d ; %d,%d)",
                    or0, or1, obj.cr0.hogOrientation, obj.cr1.hogOrientation);

            //hogCost, f, hcptHgsCost, f0, f1, costHCPT, costHGS
            sb2.append(String.format(
                "  %s (%d,%d) best: %.4f (%d,%d) %s [%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f] %s n=%d",
                debugLabel,
                Math.round(scale01 * obj.cr1.ellipseParams.xC),
                Math.round(scale01 * obj.cr1.ellipseParams.yC),
                (float) obj.cost,
                Math.round(scale00 * obj.cr0.ellipseParams.xC),
                Math.round(scale00 * obj.cr0.ellipseParams.yC), lbl,
                (float) obj.costs[0], (float) obj.costs[1],
                (float) obj.costs[2], (float) obj.costs[3],
                (float) obj.costs[4], (float) obj.costs[5], (float) obj.costs[6],
                str1, obj.cr0.offsetsToOrigCoords.size()
            ));
            /*
                    Image im0 = gsI0.copyToColorGreyscale();
                    Image im1 = gsI1.copyToColorGreyscale();
                    int[] clr = new int[]{255,0,0};
                    obj0.cr0.drawEachPixel(im0, 0, clr[0], clr[1], clr[2]);
                    obj0.cr1.drawEachPixel(im1, 0, clr[0], clr[1], clr[2]);
                    //obj0.cr1.draw(im1, 1, 0, 0, 0);
                    MiscDebug.writeImage(im0, debugLabel + "_" + lbl);
                    MiscDebug.writeImage(im1, debugLabel + "_" + lbl);
             */
            System.out.println(sb2.toString());
        }
                
    }
    
    //double[]{combIntersection, f0, f1, intersectionHCPT, intersectionHGS, count}
    private double[] sumCost2(HCPT hcpt0, HGS hgs0, CRegion cr0, 
        HCPT hcpt1, HGS hgs1, CRegion cr1) {
        
        Map<PairInt, PairInt> offsetMap1 = cr1.offsetsToOrigCoords;

        double sumCombined = 0;
        double sumHCPT = 0;
        double sumHGS = 0;
        int count = 0;
        
        int[] h0 = new int[hcpt0.getNumberOfBins()];
        int[] h1 = new int[h0.length];
        
        // key = transformed offsets, value = coords in image ref frame,
        // so, can compare dataset0 and dataset1 points with same
        //  keys
        for (Entry<PairInt, PairInt> entry0 : cr0.offsetsToOrigCoords.entrySet()) {

            PairInt pOffset0 = entry0.getKey();

            PairInt xy1 = offsetMap1.get(pOffset0);

            if (xy1 == null) {
                continue;
            }

            PairInt xy0 = entry0.getValue();

            hcpt0.extractBlock(xy0.getX(), xy0.getY(), h0);

            hcpt1.extractBlock(xy1.getX(), xy1.getY(), h1);

            float intersection = hcpt0.intersection(h0, h1);
            
            sumCombined += (intersection * intersection);
            sumHCPT += (intersection * intersection);
            
            hgs0.extractBlock(xy0.getX(), xy0.getY(), h0);

            hgs1.extractBlock(xy1.getX(), xy1.getY(), h1);

            intersection = hgs0.intersection(h0, h1);
            
            sumCombined += (intersection * intersection);
            sumHGS += (intersection * intersection);

            count++;
        }
        if (count == 0) {
            return null;
        }

        sumCombined /= (2.*count);
        sumCombined = Math.sqrt(sumCombined);
        
        sumHCPT /= (double)count;
        sumHCPT = Math.sqrt(sumHCPT);
        
        sumHGS /= (double)count;
        sumHGS = Math.sqrt(sumHGS);
                     
        double area1 = cr1.offsetsToOrigCoords.size();
        double f1 = 1. - ((double) count / area1);
        
        double area0 = cr0.offsetsToOrigCoords.size();
        double f0 = 1. - ((double) count / area0);
        
        return new double[]{sumCombined, f0, f1, sumHCPT, sumHGS, count};            
    }
    
    //double[]{combIntersection, f0, f1, intersectionHCPT, intersectionHGS, count};
    private double[] sumCost3(HCPT hcpt0, HGS hgs0, CRegion cr0, 
        HCPT hcpt1, HGS hgs1, CRegion cr1) {
        
        int orientation0 = cr0.hogOrientation;
            
        int orientation1 = cr1.hogOrientation;
        
        Map<PairInt, PairInt> offsetMap1 = cr1.offsetsToOrigCoords;

        double sum = 0;
        double sumHCPT = 0;
        double sumHGS = 0;
        double sumErrSq = 0;
        int count = 0;
        
        int[] h0 = new int[hcpt0.getNumberOfBins()];
        int[] h1 = new int[h0.length];
            
        float[] diffAndErr;
        
        //NOTE: the coordinate set block could be changed to not visit
        //   every point, but to instead use a spacing of the HOG cell length
        //   which is what Dala & Trigg recommend for their detector windows.
        
        // key = transformed offsets, value = coords in image ref frame,
        // so, can compare dataset0 and dataset1 points with same
        //  keys
        for (Entry<PairInt, PairInt> entry0 : cr0.offsetsToOrigCoords.entrySet()) {

            PairInt pOffset0 = entry0.getKey();

            PairInt xy1 = offsetMap1.get(pOffset0);

            if (xy1 == null) {
                continue;
            }

            PairInt xy0 = entry0.getValue();

            hcpt0.extractBlock(xy0.getX(), xy0.getY(), h0);

            hcpt1.extractBlock(xy1.getX(), xy1.getY(), h1);

            // 1.0 is perfect similarity
            diffAndErr = hcpt0.diff(h0, h1);
            
            sum += diffAndErr[0];
            sumHCPT += diffAndErr[0];
            sumErrSq += (diffAndErr[1] * diffAndErr[1]);
            
            
            hgs0.extractBlock(xy0.getX(), xy0.getY(), h0);

            hgs1.extractBlock(xy1.getX(), xy1.getY(), h1);

            // 1.0 is perfect similarity
            diffAndErr = hgs0.diff(h0, h1);
            
            sum += diffAndErr[0];
            sumHGS += diffAndErr[0];
            sumErrSq += (diffAndErr[1] * diffAndErr[1]);
            
            count++;
        }
        if (count == 0) {
            return null;
        }

        sum /= (2.*count);
        sumHCPT /= (double)count;
        sumHGS /= (double)count;
        
        sumErrSq /= (2.*count);
        sumErrSq = Math.sqrt(sumErrSq);
        
        double area1 = cr1.offsetsToOrigCoords.size();
        double f1 = 1. - ((double) count / area1);
        
        double area0 = cr0.offsetsToOrigCoords.size();
        double f0 = 1. - ((double) count / area0);
        
        return new double[]{sum, f0, f1, count, sumErrSq};
    }
}
