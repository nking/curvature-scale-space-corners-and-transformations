package algorithms.imageProcessing.features.mser;

import algorithms.QuickSort;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.util.List;
import algorithms.util.PairIntArray;
import algorithms.util.TrioInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
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

    public static class RegionGeometry {
        public int xC;
        public int yC;
        public double orientation;// determined from ellipse alpha
        public double eccentricity;
        public double minor;
        public double major;
    }
    
    public static class RegionPoints {
        
        public RegionGeometry ellipseParams;
        
        // NOTE: this could probably be stored more efficiently
        /**
         * key = transformed xOffset, yOffset,
         * value = coordinate in the original untransformed reference frame.
         */
        public Set<PairInt> points;

        // orientations in degrees in range 0 to 180
        public TIntList hogOrientations = new TIntArrayList();
    }

    public static class CRegion {
        
        public RegionGeometry ellipseParams = new RegionGeometry();

        public int hogOrientation;
        
        public double autocorrel;

        // NOTE: this could probably be stored more efficiently
        /**
         * key = transformed xOffset, yOffset,
         * value = coordinate in the original untransformed reference frame.
         */
        public Map<PairInt, PairInt> offsetsToOrigCoords;

        /**
         * when not empty, this holds label of segmented regions
         */
        public TIntSet labels = new TIntHashSet();
        
        public int dataIdx = -1;
        
        public void draw(Image img, int nExtraDot, int rClr, int gClr, int bClr) {

            double mc = Math.cos(ellipseParams.orientation - Math.PI/2.);
            double ms = Math.sin(ellipseParams.orientation - Math.PI/2.);
            int x1 = (int)Math.round(ellipseParams.xC - ellipseParams.major * mc);
            int y1 = (int)Math.round(ellipseParams.yC + ellipseParams.major * ms);
            int x2 = (int)Math.round(ellipseParams.xC + ellipseParams.major * mc);
            int y2 = (int)Math.round(ellipseParams.yC - ellipseParams.major * ms);
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

            mc = Math.cos(ellipseParams.orientation);
            ms = Math.sin(ellipseParams.orientation);
            x1 = (int)Math.round(ellipseParams.xC + ellipseParams.minor * mc);
            y1 = (int)Math.round(ellipseParams.yC - ellipseParams.minor * ms);
            x2 = (int)Math.round(ellipseParams.xC - ellipseParams.minor * mc);
            y2 = (int)Math.round(ellipseParams.yC + ellipseParams.minor * ms);
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

        public void drawEachPixel(Image img, int nExtraDot, int rClr, int gClr, int bClr) {

            for (Entry<PairInt, PairInt> entry : offsetsToOrigCoords.entrySet()) {
            
                PairInt p = entry.getValue();
                
                ImageIOHelper.addPointToImage(p.getX(), p.getY(), img, nExtraDot,
                    rClr, gClr, bClr);
            }
        }
        
        @Override
        public String toString() {

            String str = String.format(
                "(%d, %d) ecc=%.3f angle=%.3f major=%d minor=%d area=%d",
                ellipseParams.xC, ellipseParams.yC, 
                (float)ellipseParams.eccentricity, (float)ellipseParams.orientation,
                (int)Math.round(ellipseParams.major), 
                (int)Math.round(ellipseParams.minor),
                offsetsToOrigCoords.size());

            return str;
        }
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

        TIntObjectMap<CRegion> output = new TIntObjectHashMap<CRegion>();

        int[] xyCen = new int[2];
    
        for (int i = 0; i < regions.size(); ++i) {

            Region r = regions.get(i);
            
            CRegion cRegion = canonicalizeRegion(r, meanWindowedImg);
            
            if (cRegion != null) {
                output.put(i, cRegion);
            }
        }

        return output;
    }
    
    /**
     * uses RegionPoints.hogOrientations to make multiple cRegions for 
     * each RegionPoints instance.
     * 
     * @param regions
     * @param img
     * @return 
     */
    public TIntObjectMap<CRegion> canonicalizeRegions4(
        TIntObjectMap<RegionPoints> regions, GreyscaleImage img) {
        
        int addIdx = regions.size();
        
        TIntObjectMap<CRegion> output = new TIntObjectHashMap<CRegion>();

        int[] xyCen = new int[2];
    
        TIntObjectIterator<RegionPoints> iter = regions.iterator();
        for (int i = 0; i < regions.size(); ++i) {
            iter.advance();
            int rIdx = iter.key();
            RegionPoints r = iter.value();
            
            TIntList orientations = r.hogOrientations;
            
            PairIntArray points = Misc.convertWithoutOrder(r.points);
            
            for (int j = 0; j < orientations.size(); ++j) {
            
                double angle = orientations.get(j) * (Math.PI/180.);
                
                Map<PairInt, PairInt> offsetToOrigMap = createOffsetToOrigMap(
                    r.ellipseParams.xC, r.ellipseParams.yC,
                    points, img.getWidth(), img.getHeight(),
                    angle);

                CRegion cRegion = new CRegion();
                cRegion.ellipseParams = r.ellipseParams;
                cRegion.offsetsToOrigCoords = offsetToOrigMap;
                cRegion.dataIdx = rIdx;
                
                if (j == 0) {
                    output.put(rIdx, cRegion);
                } else {
                    output.put(addIdx, cRegion);
                    addIdx++;
                }
            }
        }

        return output;        
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
    public TIntObjectMap<CRegion> canonicalizeRegions3(
        TIntObjectMap<RegionPoints> regions, GreyscaleImage img) {

        TIntObjectMap<CRegion> output = new TIntObjectHashMap<CRegion>();

        int[] xyCen = new int[2];
    
        TIntObjectIterator<RegionPoints> iter = regions.iterator();
        for (int i = 0; i < regions.size(); ++i) {
            iter.advance();
            int label = iter.key();
            RegionPoints r = iter.value();
            
            Map<PairInt, PairInt> offsetToOrigMap = createOffsetToOrigMap(
                r.ellipseParams.xC, r.ellipseParams.yC,
                Misc.convertWithoutOrder(r.points), img.getWidth(), img.getHeight(), 
                r.ellipseParams.orientation);

            CRegion cRegion = new CRegion();
            cRegion.ellipseParams = r.ellipseParams;
            cRegion.offsetsToOrigCoords = offsetToOrigMap;
            
            output.put(label, cRegion);
        }

        return output;
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
    public TIntObjectMap<RegionPoints> canonicalizeRegions2(List<Region> regions,
        GreyscaleImage meanWindowedImg) {

        TIntObjectMap<RegionPoints> output = new TIntObjectHashMap<RegionPoints>();

        int[] xyCen = new int[2];
    
        for (int i = 0; i < regions.size(); ++i) {

            Region r = regions.get(i);
            
            RegionPoints cRegion = canonicalizeRegion2(r, 
                meanWindowedImg.getWidth(), meanWindowedImg.getHeight());
            
            output.put(i, cRegion);
        }

        return output;
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
    public CRegion canonicalizeRegion(Region r, GreyscaleImage meanWindowedImg) {

        int imageWidth = meanWindowedImg.getWidth();
        int imageHeight = meanWindowedImg.getHeight();

        CRegion cRegion = canonicalizeRegion(r, imageWidth, imageHeight);

        if (cRegion == null) {
            return null;
        }
        
        double autocorrel = calcAutoCorrel(meanWindowedImg, 
            cRegion.ellipseParams.xC, cRegion.ellipseParams.yC, 
            cRegion.offsetsToOrigCoords);

        cRegion.autocorrel = Math.sqrt(autocorrel)/255.;

        return cRegion;
    }
    
    /**
     * NOTE, for best use, invoker should use this descriptor with
     * an image processed to create a window average for each pixel.
       For example:
       SummedAreaTable sumTable = new SummedAreaTable();
       GreyscaleImage imgM = sumTable.createAbsoluteSummedAreaTable(img);
       imgM = sumTable.applyMeanOfWindowFromSummedAreaTable(imgM,
            2*halfDimension + 1);
     * @param r
     * @param imageWidth
     * @param imageHeight
     * @return
     */
    public RegionPoints canonicalizeRegion2(Region r, int imageWidth, 
        int imageHeight) {

        int[] xyCen = new int[2];
    
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

        PairIntArray xy = new PairIntArray();

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
                    xy.add(xE, yE);
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
                    xy.add(xE, yE);
                }
            }
        }
        
        fillInEllipse(xy);
        
        RegionGeometry rg = new RegionGeometry();
        rg.eccentricity = ecc;
        rg.major = major;
        rg.minor = minor;
        rg.orientation = angle;
        rg.xC = x;
        rg.yC = y;

        RegionPoints regionPoints = new RegionPoints();
        regionPoints.ellipseParams = rg;
        regionPoints.points = Misc.convert(xy);
        
        return regionPoints;
    }
    
    /**
     * NOTE, for best use, invoker should use this descriptor with
     * an image processed to create a window average for each pixel.
       For example:
       SummedAreaTable sumTable = new SummedAreaTable();
       GreyscaleImage imgM = sumTable.createAbsoluteSummedAreaTable(img);
       imgM = sumTable.applyMeanOfWindowFromSummedAreaTable(imgM,
            2*halfDimension + 1);

       NOTE the field utocorrelation is not clculated here and my be removed.
    
     * @return
     */
    public CRegion canonicalizeRegion(Region r, int imageWidth, int imageHeight) {

        int[] xyCen = new int[2];
    
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
        if (Double.isNaN(ecc)) {
            return null;
        }

        PairIntArray xy = new PairIntArray();

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
                    xy.add(xE, yE);
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
                    xy.add(xE, yE);
                }
            }
        }

        Map<PairInt, PairInt> offsetToOrigMap = createOffsetToOrigMap(x, y,
            xy, imageWidth, imageHeight, angle);

        RegionGeometry rg = new RegionGeometry();
        rg.eccentricity = ecc;
        rg.major = major;
        rg.minor = minor;
        rg.orientation = angle;
        rg.xC = x;
        rg.yC = y;

        CRegion cRegion = new CRegion();
        cRegion.ellipseParams = rg;
        cRegion.offsetsToOrigCoords = offsetToOrigMap;

        return cRegion;
    }
    
    public static RegionGeometry calculateEllipseParams(Region r, int imageWidth, int imageHeight) {

        int[] xyCen = new int[2];
    
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
        if (Double.isNaN(ecc)) {
            return null;
        }

        RegionGeometry rg = new RegionGeometry();
        rg.eccentricity = ecc;
        rg.major = major;
        rg.minor = minor;
        rg.orientation = angle;
        rg.xC = x;
        rg.yC = y;

        return rg;
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

                int xc2 = Math.round((float)cr.ellipseParams.xC/scale);
                int yc2 = Math.round((float)cr.ellipseParams.yC/scale);

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
                    (float)(cr.ellipseParams.orientation*180./Math.PI)));

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

                double major = cr.ellipseParams.major/scale;
                double minor = cr.ellipseParams.minor/scale;

                double ecc = Math.sqrt(major * major - minor * minor)/major;
                assert(!Double.isNaN(ecc));

                RegionGeometry rg = new RegionGeometry();
                rg.eccentricity = ecc;
                rg.major = major;
                rg.minor = minor;
                rg.orientation = cr.ellipseParams.orientation;
                rg.xC = xc2;
                rg.yC = yc2;
        
                CRegion cRegion = new CRegion();
                cRegion.ellipseParams = rg;
                cRegion.offsetsToOrigCoords = offsetMap;
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

    public static Map<PairInt, PairInt> createOffsetToOrigMap(
        int x, int y, PairIntArray xy, int imgWidth, int imgHeight,
        double orientation) {
        
        fillInEllipse(xy);

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
    
        // ellipse, rotated by orientation to create
        //   x and y offsets from center that are comparable to
        //   the same for CRegions in other dataset.
        PairIntArray xyTr = transformer.applyTransformation(params, xy);

        // visit all points in region and transform them
        // also determine the auto-correlation

        int nc = 0;
        for (int j = 0; j < xy.getN(); ++j) {
            
            int xp = xy.getX(j);
            int yp = xy.getY(j);
            
            int xpTr = xyTr.getX(j);
            int ypTr = xyTr.getY(j);
            
            if (xpTr == -1) {
                xpTr = 0;
            }
            if (ypTr == -1) {
                ypTr = 0;
            }
            if (xpTr == imgWidth) {
                xpTr--;
            }
            if (ypTr == imgHeight) {
                ypTr--;
            }
            if (xpTr < 0 || ypTr < 0 || xpTr >= imgWidth || ypTr >= imgHeight) {
                continue;
            }
            
            PairInt pt = new PairInt(xpTr, ypTr);
            if (visited.contains(pt)) {
                continue;
            }
            visited.add(pt);

            PairInt pOffsets = new PairInt(xpTr - x, ypTr - y);
            offsetToOrigMap.put(pOffsets, new PairInt(xp, yp));
            nc++;
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
            
            PairInt p = entry.getValue();
            
            int x2 = p.getX();
            
            int y2 = p.getY();
            
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
    
    private static void fillInEllipse(PairIntArray xy) {

        // key = row number, value = start and stop x range
        TIntObjectMap<PairInt> rowColRange = new TIntObjectHashMap<PairInt>();
        
        int minRow = Integer.MAX_VALUE;
        int maxRow = Integer.MIN_VALUE;
        
        for (int i = 0; i < xy.getN(); ++i) {
            
            int row = xy.getY(i);
            int col = xy.getX(i);
            
            PairInt xMinMax = rowColRange.get(row);
            if (xMinMax == null) {
                xMinMax = new PairInt(col, col);
                rowColRange.put(row, xMinMax);
            } else {
                if (col < xMinMax.getX()) {
                    xMinMax.setX(col);
                } else if (col > xMinMax.getY()) {
                    xMinMax.setY(col);
                }
            }
            
            if (row < minRow) {
                minRow = row;
            }
            if (row > maxRow) {
                maxRow = row;
            }
        }
    
        for (int i = minRow; i <= maxRow; ++i) {
            PairInt xMinMax = rowColRange.get(i);
            int x0, x1;
            if (xMinMax == null) {
                PairInt prev = null;
                int off = 0;
                while (prev == null) {
                    off--;
                    prev = rowColRange.get(i + off);
                }
                PairInt next = null;
                off = 0;
                while (next == null) {
                    off++;
                    next = rowColRange.get(i + off);
                }
                assert(prev != null);
                assert(next != null);
                x0 = Math.round(((float)prev.getX() + (float)next.getX())/2.f);
                x1 = Math.round(((float)prev.getY() + (float)next.getY())/2.f);
            } else {
                x0 = xMinMax.getX();
                x1 = xMinMax.getY();
            }
            for (int x = (x0 + 1); x < x1; ++x) {
                xy.add(x, i);
            }
        }
    }

}
