package algorithms.imageProcessing.features.mser;

import algorithms.QuickSort;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.util.List;
import algorithms.util.PairIntArray;
import algorithms.util.TrioInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * class to produce descriptors for MSER, usable for
 * matching same objects in other images.
 *
 * @author nichole
 */
public class Canonicalizer {

    /*
    affine transformation
       see
          http://math.stackexchange.com/questions/310776/finding-the-affine-transformation-that-will-change-a-given-ellipse-into-the-unit

    x(t) = xCenter + aParam*cos(alpha)*cos(t) − bParam*sin(alpha)*sin(t)
    y(t) = yCenter + aParam*sin(alpha)*cos(t) + bParam*cos(alpha)*sin(t)

    v0x = aParam * cos(alpha);
    v1x = bParam * sin(alpha);
    v0y = aParam * sin(alpha);
    v1y = bParam * cos(alpha);

    so the semi-major axis length = 2 * Math.max(v0x, v1x)
    and the semi-minor axis length = 2 * Math.min(v0y, v1y)

    NOTE, the canonicalization currently preserves the
    ellipse shape, but transforms the coordinates to
    an orientation of 90 degrees for comparison in
    matching algorithms.

    The matching algorithm is still in progress, but can see
    -- might need canonicalization to be circular areas
       after all.  the radius would be the minor axis length.
       -- for the regions where minor axis length < 16
          should use a radius of 16.
    -- if the canoicalization into circular areas is successful
       at finding best true matches,
       will try to reduce them all to radius of 16 and if that
       is successful, then should be able to use the ORB descriptors
       in plae of these.
       the ORB descriptors are binary...comparisons
       are very fast.
    */

    public static class CRegion {
        public int xC;
        public int yC;
        public double orientation;// determined from ellipse alpha
        public double eccentricity;
        public double minor;
        public double major;
        public double autocorrel;

        // number of pixels in the transformed ellipse
        public int nTrEllipsePixels;

        // NOTE: this could probably be stored more efficiently
        /**
         * key = transformed xOffset, yOffset,
         * value = coordinate in the original untransformed reference frame.
         */
        public Map<PairInt, PairInt> offsetsToOrigCoords;

        public void draw(Image img, int nExtraDot, int rClr, int gClr, int bClr) {

            double mc = Math.cos(orientation - Math.PI/2.);
            double ms = Math.sin(orientation - Math.PI/2.);
            int x1 = (int)Math.round(xC - major * mc);
            int y1 = (int)Math.round(yC + major * ms);
            int x2 = (int)Math.round(xC + major * mc);
            int y2 = (int)Math.round(yC - major * ms);
            if (x1 < 0) { x1 = 0;}
            if (y1 < 0) { y1 = 0;}
            if (x2 < 0) { x2 = 0;}
            if (y2 < 0) { y2 = 0;}
            if (x1 >= img.getWidth()) { x1 = img.getWidth() - 1;}
            if (y1 >= img.getHeight()) { y1 = img.getHeight() - 1;}
            if (x2 >= img.getWidth()) { x2 = img.getWidth() - 1;}
            if (y2 >= img.getHeight()) { y2 = img.getHeight() - 1;}

            ImageIOHelper.drawLineInImage(x1, y1, x2, y2, img, nExtraDot,
                rClr, gClr, bClr);

            mc = Math.cos(orientation);
            ms = Math.sin(orientation);
            x1 = (int)Math.round(xC + minor * mc);
            y1 = (int)Math.round(yC - minor * ms);
            x2 = (int)Math.round(xC - minor * mc);
            y2 = (int)Math.round(yC + minor * ms);
            if (x1 < 0) { x1 = 0;}
            if (y1 < 0) { y1 = 0;}
            if (x2 < 0) { x2 = 0;}
            if (y2 < 0) { y2 = 0;}
            if (x1 >= img.getWidth()) { x1 = img.getWidth() - 1;}
            if (y1 >= img.getHeight()) { y1 = img.getHeight() - 1;}
            if (x2 >= img.getWidth()) { x2 = img.getWidth() - 1;}
            if (y2 >= img.getHeight()) { y2 = img.getHeight() - 1;}

            ImageIOHelper.drawLineInImage(x1, y1, x2, y2, img, nExtraDot,
                rClr, gClr, bClr);

        }

        @Override
        public String toString() {

            String str = String.format(
                "(%d, %d) ecc=%.3f angle=%.3f major=%d minor=%d area=%d",
                xC, yC, (float)eccentricity, (float)orientation,
                (int)Math.round(major), (int)Math.round(minor),
                nTrEllipsePixels);

            return str;
        }
    }

    /**
     * class composed of labeled segmented regions reformed to fit
     * the use of canonicalized mser region
     */
    public static class CSRegion extends CRegion {
        // the segmented cell label index
        public int label;
    }

    /**
     * NOTE, for best use, invoker should use this descriptor with
     * an image processed to create a window average for each pixel.
       For example:
       SummedAreaTable sumTable = new SummedAreaTable();
       GreyscaleImage imgM = sumTable.createAbsoluteSummedAreaTable(img);
       imgM = sumTable.applyMeanOfWindowFromSummedAreaTable(imgM,
            2*halfDimension + 1);

     * @return
     */
    public TIntObjectMap<CRegion> canonicalizeRegions(List<Region> regions,
        GreyscaleImage meanWindowedImg) {

        int imageWidth = meanWindowedImg.getWidth();
        int imageHeight = meanWindowedImg.getHeight();

        TIntObjectMap<CRegion> output = new TIntObjectHashMap<CRegion>();

        Transformer transformer = new Transformer();

        int[] xyCen = new int[2];
    
        for (int i = 0; i < regions.size(); ++i) {

            Region r = regions.get(i);

            r.calculateXYCentroid(xyCen, imageWidth, imageHeight);
            int x = xyCen[0];
            int y = xyCen[1];
            assert(x >= 0 && x < imageWidth);
            assert(y >= 0 && y < imageHeight);

            //v0x, v1x, v0y, v1y
            double[] m = r.calcParamTransCoeff();

            double angle = Math.atan(m[0]/m[2]);
            if (angle < 0) {
                angle += Math.PI;
            }

            double major = 2. * m[4];
            double minor = 2. * m[5];

            double ecc = Math.sqrt(major * major - minor * minor)/major;
            assert(!Double.isNaN(ecc));

            TIntList xs = new TIntArrayList();
            TIntList ys = new TIntArrayList();

            boolean createEllipse = true;
            double radius = minor;
            if (radius < 4) {
                radius = 4;
                createEllipse = false;
            }

            if (createEllipse) {
                // elliptical bounds
                // find the ranges of the untransformed ellipse first
                for (double t = 0.0; t < 2.0 * Math.PI; t += 0.001) {
                    int xE = (int)Math.round(x +
                        (Math.cos(t) * m[0] + Math.sin(t) * m[1]) * 2.0 + 0.5);
                    int yE = (int)Math.round(y + (Math.cos(t) * m[2]
                        + Math.sin(t) * m[3]) * 2.0 + 0.5);
                    if ((xE >= 0) && (xE < imageWidth) &&
                        (yE >= 0) && (yE < imageHeight)) {
                        xs.add(xE);
                        ys.add(yE);
                    }
                }
            } else {
                for (double t = 0.0; t < 2.0 * Math.PI; t += 0.001) {
                    double mc = Math.cos(t);
                    double ms = Math.sin(t);
                    int xE = (int)Math.round(x + (mc * radius));
                    int yE = (int)Math.round(y + (ms * radius));
                    if ((xE >= 0) && (xE < imageWidth) &&
                        (yE >= 0) && (yE < imageHeight)) {
                        xs.add(xE);
                        ys.add(yE);
                    }
                }
            }

            Map<PairInt, PairInt> offsetToOrigMap = createOffsetToOrigMap(x, y,
                xs, ys, meanWindowedImg.getWidth(), 
                meanWindowedImg.getHeight(), angle);

            double autocorrel = calcAutoCorrel(meanWindowedImg, x, y, 
                offsetToOrigMap);
                        
            CRegion cRegion = new CRegion();
            cRegion.eccentricity = ecc;
            cRegion.offsetsToOrigCoords = offsetToOrigMap;
            cRegion.major = major;
            cRegion.minor = minor;
            cRegion.orientation = angle;
            cRegion.xC = x;
            cRegion.yC = y;
            cRegion.nTrEllipsePixels = offsetToOrigMap.size();
            cRegion.autocorrel = Math.sqrt(autocorrel)/255.;

            output.put(i, cRegion);
        }

        return output;
    }

    /**
     * create canonicalized regions containing coordinate maps that can be used
     * to make descriptors.  Note that for best results, the pyramidImgs
     * should have been pre-processed so that a pixel contains the mean of
     * itself and its neighboring pixels.
     * @param regions
     * @param pyramidImgs
     * @return
     */
    public List<TIntObjectMap<CRegion>> canonicalizeRegions(
        List<Region> regions, List<GreyscaleImage> pyramidImgs) {

        List<TIntObjectMap<CRegion>> output = new ArrayList<TIntObjectMap<CRegion>>();

        TIntObjectMap<CRegion> crMap0 = canonicalizeRegions(regions, pyramidImgs.get(0));

        output.add(crMap0);

        Transformer transformer = new Transformer();

        TIntObjectIterator<CRegion> iter0;
        GreyscaleImage mImg0 = pyramidImgs.get(0);

        for (int imgIdx = 1; imgIdx < pyramidImgs.size(); ++imgIdx) {

            TIntObjectMap<CRegion> crMap = new TIntObjectHashMap<CRegion>();
            output.add(crMap);

            GreyscaleImage mImg = pyramidImgs.get(imgIdx);
            float scale = ((float)mImg0.getWidth()/(float)mImg.getWidth()) +
                ((float)mImg0.getHeight()/(float)mImg.getHeight());
            scale /= 2.f;

            iter0 = crMap0.iterator();
            for (int i = 0; i < crMap0.size(); ++i) {
                iter0.advance();

                int idx = iter0.key();
                CRegion cr = iter0.value();

                int xc2 = Math.round((float)cr.xC/scale);
                int yc2 = Math.round((float)cr.yC/scale);

                if (xc2 < 0 || yc2 < 0 || xc2 >= mImg.getWidth() ||
                    yc2 >= mImg.getHeight()) {
                    continue;
                }

                // transform the elliptical range so that angle2 is pointing
                // up in the image, that is angle2 = Math.PI/2
                // pi/2 = deltaR + angle --> deltaR = (pi/2) - angle
                TransformationParameters params = new TransformationParameters();
                params.setOriginX(xc2);
                params.setOriginY(yc2);
                params.setScale(1.0f);
                params.setTranslationX(0);
                params.setTranslationY(0);
                params.setRotationInDegrees(AngleUtil.getAngleDifference(90.f,
                    (float)(cr.orientation*180./Math.PI)));

                // copy and reduce structure in size by scale factor
                Map<PairInt, PairInt> offsetMap = new HashMap<PairInt, PairInt>();

                int vc = mImg.getValue(xc2, yc2);
                int nc = 0;
                double autocorSum = 0;

                for (Map.Entry<PairInt, PairInt> entry : cr.offsetsToOrigCoords.entrySet()) {

                    PairInt pOrig = entry.getValue();

                    int xScaled = Math.round((float)pOrig.getX()/scale);
                    int yScaled = Math.round((float)pOrig.getY()/scale);
                    if (xScaled == -1) {
                        xScaled = 0;
                    }
                    if (yScaled == -1) {
                        yScaled = 0;
                    }
                    if (xScaled == mImg.getWidth()) {
                        xScaled = mImg.getWidth() - 1;
                    }
                    if (yScaled == mImg.getHeight()) {
                        yScaled = mImg.getHeight() - 1;
                    }
                    PairInt pOrigScaled = new PairInt(xScaled, yScaled);

                    // TODO: review this...should be the same as entry.getKey/scale
                    // but slightly better integer rounding results
                    double[] xyETr = transformer.applyTransformation(params,
                        xScaled, yScaled);

                    int xETr = (int)Math.round(xyETr[0]);
                    int yETr = (int)Math.round(xyETr[1]);

                    if ((xETr >= 0) && (xETr < mImg.getWidth())
                        && (yETr >= 0) && (yETr < mImg.getHeight())) {

                        PairInt pTrOffset = new PairInt(xETr - xc2, yETr - yc2);

                        if (!offsetMap.containsKey(pTrOffset)) {

                            offsetMap.put(pTrOffset, pOrigScaled);

                            int diff = mImg.getValue(xETr, yETr) - vc;

                            autocorSum += (diff * diff);
                            nc++;
                        }
                    }
                }
                autocorSum /= (double)nc;

                double major = cr.major/scale;
                double minor = cr.minor/scale;

                double ecc = Math.sqrt(major * major - minor * minor)/major;
                assert(!Double.isNaN(ecc));

                CRegion cRegion = new CRegion();
                cRegion.eccentricity = ecc;
                cRegion.offsetsToOrigCoords = offsetMap;
                cRegion.major = major;
                cRegion.minor = minor;
                cRegion.orientation = cr.orientation;
                cRegion.xC = xc2;
                cRegion.yC = yc2;
                cRegion.nTrEllipsePixels = offsetMap.size();
                cRegion.autocorrel = Math.sqrt(autocorSum)/255.;

                crMap.put(idx, cRegion);
            }
        }

        return output;
    }

    public PairIntArray calculateCentroids(List<Region> regions,
        int[] greyscale, int imageWidth, int imageHeight) {

        //return calculateIntensityCentroids(regions,
        //    greyscale, imageWidth, imageHeight);

        return extractRegionXYCenters(regions, imageWidth, imageHeight);
    }

    PairIntArray extractRegionXYCenters(List<Region> regions,
        int imageWidth, int imageHeight) {

        int[] xyCen = new int[2];

        PairIntArray output = new PairIntArray(regions.size());

        for (int i = 0; i < regions.size(); ++i) {
            Region r = regions.get(i);
            r.calculateXYCentroid(xyCen, imageWidth, imageHeight);

            output.add(xyCen[0], xyCen[1]);
        }

        return output;
    }

    PairIntArray calculateIntensityCentroids(List<Region> regions,
        int[] greyscale, int imageWidth, int imageHeight) {

        int n = regions.size();

        PairIntArray output = new PairIntArray(n);

        int[] cenXY = new int[2];

        for (int i = 0; i < n; ++i) {

            Region region = regions.get(i);

            region.calculateIntensityCentroid(greyscale, imageWidth,
                imageHeight, cenXY);

            output.add(cenXY[0], cenXY[1]);
        }

        return output;
    }

    private PairIntArray calculateIntensityCentroids(List<Region> regions,
        int radius, int[] greyscale, int imageWidth, int imageHeight) {

        int n = regions.size();

        PairIntArray output = new PairIntArray(n);

        int[] cenXY = new int[2];

        for (int i = 0; i < n; ++i) {

            Region r = regions.get(i);

            r.calculateIntensityCentroid(greyscale, imageWidth,
                imageHeight, cenXY, radius);

            output.add(cenXY[0], cenXY[1]);
        }

        return output;
    }

    public TIntObjectMap<CSRegion> selectSets(TIntObjectMap<CRegion> cRegions,
        List<Set<PairInt>> labeledSets, TObjectIntMap<PairInt> pointLabelMap,
        GreyscaleImage img) {

        int imgWidth = img.getWidth();
        int imgHeight = img.getHeight();
        
        TIntObjectMap<TIntSet> labelRegionsMap = new TIntObjectHashMap<TIntSet>();
        TIntObjectIterator<CRegion> iter = cRegions.iterator();
        for (int i = 0; i < cRegions.size(); ++i) {
            iter.advance();
            int rIdx = iter.key();
            CRegion cr = iter.value();
            PairInt p = new PairInt(cr.xC, cr.yC);
            assert(pointLabelMap.containsKey(p));
            int label = pointLabelMap.get(p);

            TIntSet rSet = labelRegionsMap.get(label);
            if (rSet == null) {
                rSet = new TIntHashSet();
                labelRegionsMap.put(label, rSet);
            }
            rSet.add(rIdx);
        }

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        Transformer transformer = new Transformer();

        // for each key in labelRegionsMap, make a combined cregion
        TIntObjectMap<CSRegion> csrMap = new TIntObjectHashMap<CSRegion>();

        TIntObjectIterator<TIntSet> iter2 = labelRegionsMap.iterator();
        for (int i = 0; i < labelRegionsMap.size(); ++i) {
            iter2.advance();
            int label = iter2.key();
            TIntSet rIdxs = iter2.value();

            TIntIterator iter3 = rIdxs.iterator();
            int nTotal = 0;
            while (iter3.hasNext()) {
                int rIdx = iter3.next();
                CRegion cr = cRegions.get(rIdx);
                nTotal += cr.offsetsToOrigCoords.size();
            }

            double avgOrientation = 0;

            iter3 = rIdxs.iterator();
            while (iter3.hasNext()) {
                int rIdx = iter3.next();

                CRegion cr = cRegions.get(rIdx);

                int n = cr.offsetsToOrigCoords.size();

                double w = (double)n/(double)nTotal;

                avgOrientation += (w * cr.orientation);
            }

            // create a Region to use the accumulate method
            // and compare results

            TIntList xs = new TIntArrayList();
            TIntList ys = new TIntArrayList();

            Set<PairInt> labeledSet = labeledSets.get(label);
            Region r = new Region();
            for (PairInt pl : labeledSet) {
                r.accumulate(pl.getX(), pl.getY());
                xs.add(pl.getX());
                ys.add(pl.getY());
            }

            int[] xy = new int[2];
             r.calculateXYCentroid(xy, imgWidth, imgWidth);
            int x = xy[0];
            int y = xy[1];
            assert(x >= 0 && x < imgWidth);
            assert(y >= 0 && y < imgWidth);
            double[] m = r.calcParamTransCoeff();

            double angle = Math.atan(m[0]/m[2]);
            if (angle < 0) {
                angle += Math.PI;
            }

            double major = 2. * m[4];
            double minor = 2. * m[5];

            double ecc = Math.sqrt(major * major - minor * minor)/major;
            assert(!Double.isNaN(ecc));

            System.out.println("avg or=" +
                avgOrientation + " tensor orn="
                + angle + " (" + x + "," + y + ")");

            Map<PairInt, PairInt> offsetMap = createOffsetToOrigMap(x, y,
                xs, ys, imgWidth, imgHeight, angle);
            
            double autocorrel = calcAutoCorrel(img, x, y, offsetMap);
            
            CSRegion csRegion = new CSRegion();
            csRegion.label = label;
            csRegion.orientation = angle;
            csRegion.eccentricity = ecc;
            csRegion.major = major;
            csRegion.minor = minor;
            csRegion.xC = x;
            csRegion.yC = y;
            csRegion.offsetsToOrigCoords = offsetMap;
            csRegion.autocorrel = Math.sqrt(autocorrel)/255.;
            
            csrMap.put(label, csRegion);
        }

        return csrMap;
    }

    public static Map<PairInt, PairInt> createOffsetToOrigMap(
        int x, int y, TIntList xs, TIntList ys, int imgWidth, int imgHeight,
        double orientation) {

        // calc content for offsetsToOrigCoords
        int[] x2s = xs.toArray(new int[xs.size()]);
        int[] y2s = ys.toArray(new int[xs.size()]);

        QuickSort.sortBy1stThen2nd(y2s, x2s);

        int n = (new TIntHashSet(ys)).size();

        TrioInt[] ranges = new TrioInt[n];
        int count = 0;
        int idx = 0;
        for (idx = 0; idx < x2s.length; ++idx) {
            int yC = y2s[idx];
            int minX = Integer.MAX_VALUE;
            int maxX = Integer.MIN_VALUE;
            int idx2 = idx;
            for (idx2 = idx; idx2 < x2s.length; ++idx2) {
                if (y2s[idx2] > yC) {
                    idx2--;
                    break;
                }
                if (minX == Integer.MAX_VALUE) {
                    minX = x2s[idx2];
                }
                maxX = x2s[idx2];
            }
            idx = idx2;

            // subtract (xc,yc) from ranges
            ranges[count] = new TrioInt(yC, minX, maxX);
            count++;
        }
        assert(count == n);

        Set<PairInt> visited = new HashSet<PairInt>();

        // transform the elliptical range so that angle2 is pointing
        // up in the image, that is angle2 = Math.PI/2
        // pi/2 = deltaR + angle --> deltaR = (pi/2) - angle
        TransformationParameters params = new TransformationParameters();
        params.setOriginX(x);
        params.setOriginY(y);
        params.setScale(1.0f);
        params.setTranslationX(0);
        params.setTranslationY(0);
        params.setRotationInDegrees(AngleUtil.getAngleDifference(90.f,
            (float)(orientation*180./Math.PI)));

        Map<PairInt, PairInt> offsetToOrigMap = new HashMap<PairInt, PairInt>();

        Transformer transformer = new Transformer();

        // visit all points in region and transform them
        // also determine the auto-correlation

        int nc = 0;
        for (int j = 0; j < ranges.length; ++j) {
            int yE = ranges[j].getX();
            int xMin = ranges[j].getY();
            int xMax = ranges[j].getZ();
            for (int xE = xMin; xE <= xMax; ++xE) {

                double[] xyETr = transformer.applyTransformation(params, xE, yE);

                int xETr = (int)Math.round(xyETr[0]);
                int yETr = (int)Math.round(xyETr[1]);

                if ((xETr >= 0) && (xETr < imgWidth)
                    && (yETr >= 0) && (yETr < imgHeight)) {

                    PairInt pt = new PairInt(xETr, yETr);
                    if (visited.contains(pt)) {
                        continue;
                    }
                    visited.add(pt);

                    PairInt pOffsets = new PairInt(xETr - x, yETr - y);
                    offsetToOrigMap.put(pOffsets, new PairInt(xE, yE));
                    nc++;
                }
            }
        }

        return offsetToOrigMap;
    }
    
    public static double calcAutoCorrel(GreyscaleImage img, int x, int y, 
        Map<PairInt, PairInt> offsetMap) {
        
        int vc = img.getValue(x, y);
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        int nc = 0;
        
        double autocorSum = 0;
        
        Set<PairInt> visited = new HashSet<PairInt>();
        
        for (Entry<PairInt, PairInt> entry : offsetMap.entrySet()) {
            
            PairInt pOffset = entry.getKey();
            
            int x2 = pOffset.getX() + x;
            
            int y2 = pOffset.getY() + y;
            
            if (x2 == w) {
                x2--;
            }
            if (y2 == h) {
                y2--;
            }
            
            PairInt p2 = new PairInt(x2, y2);
            
            if (visited.contains(p2)) {
                continue;
            }
            visited.add(p2);
            
            int v = img.getValue(x2, y2);
           
            int diff = v - vc;

            autocorSum += (diff * diff);
            nc++;
        }
        
        autocorSum /= (double)nc;
        
        return autocorSum;
    }

}
