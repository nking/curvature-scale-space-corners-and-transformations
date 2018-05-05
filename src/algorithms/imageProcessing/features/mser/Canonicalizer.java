package algorithms.imageProcessing.features.mser;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.util.PairInt;
import java.util.List;
import algorithms.util.PairIntArray;
import algorithms.util.PixelHelper;
import com.climbwithyourfeet.clustering.ClusterFinder;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.iterator.TLongIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TLongIntMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TLongIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TIntHashSet;
import gnu.trove.set.hash.TLongHashSet;
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

    public static class RegionGeometry {
        public int xC;
        public int yC;
        public double orientation;// determined from ellipse alpha
        public double eccentricity;
        public double minor;
        public double major;
        
        /**
         * using the region moments and area, calculates parameters needed
         * for a bounding ellipse.
         * 
       The equation for rotation of an ellipse by angle alpha:
        x(t) = xCenter + aParam*cos(alpha)*cos(t) − bParam*sin(alpha)*sin(t)
        y(t) = yCenter + aParam*sin(alpha)*cos(t) + bParam*cos(alpha)*sin(t)
        
        returns 
            v0x = aParam * cos(alpha);
            v1x = bParam * sin(alpha);
            v0y = aParam * sin(alpha);
            v1y = bParam * cos(alpha);
            e0 = first eigenvalue
            e1 = 2nd eigenvalue
        m is []{v0x, v1x, v0y, v1y} */
        public double[] m;
       
        public RegionGeometry createNewDividedByScale(float scale,
            int maxX, int maxY) {
            RegionGeometry rg = new RegionGeometry();
            rg.xC = Math.round((float) xC / scale);
            if (rg.xC > maxX) {
                rg.xC = maxX;
            }
            rg.yC = Math.round((float) yC / scale);
            if (rg.yC > maxY) {
                rg.yC = maxY;
            }
            rg.orientation = orientation;
            rg.eccentricity = eccentricity;
            rg.minor = minor/scale;
            rg.major = major/scale;
            return rg;
        }
    }
    
    public static class RegionPoints {
        
        public RegionGeometry ellipseParams;
        
        // NOTE: this could probably be stored more efficiently
        /**
         * key = transformed xOffset, yOffset,
         * value = coordinate in the original untransformed reference frame.
         */
        //public Set<PairInt> points;

        // orientations in degrees in range 0 to 180
        public TIntList hogOrientations = new TIntArrayList();
        
        // temporary storage of original accumulated points to explore partial edges
        public TIntList accX = new TIntArrayList();
        public TIntList accY = new TIntArrayList();
        
        private int[] minMaxXY = null;
        
        public RegionPoints() {
            
        }
        
        public int[] getMinMaxXY() {
            if (minMaxXY == null) {
                minMaxXY = new int[]{Integer.MAX_VALUE, Integer.MIN_VALUE, 
                    Integer.MAX_VALUE, Integer.MIN_VALUE};
                for (int i = 0; i < accX.size(); ++i) {
                    if (accX.get(i) < minMaxXY[0]) {
                        minMaxXY[0] = accX.get(i);
                    }
                    if (accX.get(i) > minMaxXY[1]) {
                        minMaxXY[1] = accX.get(i);
                    }
                    if (accY.get(i) < minMaxXY[2]) {
                        minMaxXY[2] = accY.get(i);
                    }
                    if (accY.get(i) > minMaxXY[3]) {
                        minMaxXY[3] = accY.get(i);
                    }
                }
            }
            return minMaxXY;
        }
        
        public TLongSet createAccPixelCoords(int imageWidth) {
            TLongSet pixs = new TLongHashSet();
            PixelHelper ph = new PixelHelper();
            for (int i = 0; i < accX.size(); ++i) {
                pixs.add(ph.toPixelIndex(accX.get(i), accY.get(i), imageWidth));
            }
            return pixs;
        }
        
        public TLongSet createAccPixelCoords(int[] minMaxXY2) {
            int xOffset = minMaxXY2[0];
            int yOffset = minMaxXY2[1];
            int w2 = minMaxXY2[1] - minMaxXY2[0] + 1;
            int h2 = minMaxXY2[3] - minMaxXY2[2] + 1;
            TLongSet pixs = new TLongHashSet();
            PixelHelper ph = new PixelHelper();
            for (int i = 0; i < accX.size(); ++i) {
                int x = accX.get(i) - xOffset;
                int y = accY.get(i) - yOffset;
                pixs.add(ph.toPixelIndex(x, y, w2));
            }
            return pixs;
        }
        
        public RegionPoints createNewDividedByScale(float scale,
            int maxX, int maxY) {
            RegionPoints rp = new RegionPoints();
            rp.ellipseParams = ellipseParams.createNewDividedByScale(scale,
                maxX, maxY);
            rp.hogOrientations.addAll(hogOrientations);
            for (int i = 0; i < accX.size(); ++i) {
                int x = Math.round(accX.get(i)/scale);
                int y = Math.round(accY.get(i)/scale);
                if (x > maxX) {
                    x = maxX;
                }
                if (y > maxY) {
                    y = maxY;
                }
                rp.accX.add(x);
                rp.accY.add(y);
            }
            return rp;
        }
        
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
        private final Map<PairInt, PairInt> offsetsToOrigCoords =
            new HashMap<PairInt, PairInt>();

        public final int imgWidth;
        public final int imgHeight;
        
        // these are the values in offsetsToOrigCoords 
        private final TLongSet origPixs = new TLongHashSet();
        /**
         * when not empty, this holds label of segmented regions
         */
        public TIntSet labels = new TIntHashSet();
        
        public int dataIdx = -1;
        
        private int[] minMaxXY = null;
        
        public CRegion(int imageWidth, int imageHeight) {
            imgWidth = imageWidth;
            imgHeight = imageHeight;
        }
        
        public void addAllOffsets(Map<PairInt, PairInt> offsetCoords) {
            PixelHelper ph = new PixelHelper();
            for (Entry<PairInt, PairInt> entry : offsetCoords.entrySet()) {
                PairInt p = entry.getValue();
                offsetsToOrigCoords.put(entry.getKey(), p);
                origPixs.add(ph.toPixelIndex(p, imgWidth));
            }
        }
        public void resetToTheseOffsets(Map<PairInt, PairInt> offsetCoords) {
            offsetCoords.clear();
            origPixs.clear();
            addAllOffsets(offsetCoords);
        }
        
        public TLongSet getPixelCoords() {
            return origPixs;
        }
        
        public Set<PairInt> getOffsetKeys() {
            return offsetsToOrigCoords.keySet();
        }
        
        public Map<PairInt, PairInt> getOffsetsToOrigCoords() {
            return offsetsToOrigCoords;
        }
        
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
        
        public Set<PairInt> extractCoords() {
            return new HashSet<PairInt>(offsetsToOrigCoords.values());            
        }
        
        public CRegion createNewDividedByScale(float scale, 
            int maxX, int maxY) {
            
            int w2 = maxX + 1;
            int h2 = maxY + 1;
            
            PixelHelper ph = new PixelHelper();
            
            CRegion r = new CRegion(w2, h2);
            r.ellipseParams = ellipseParams.createNewDividedByScale(scale,
                maxX, maxY);
            r.hogOrientation = hogOrientation;
            r.autocorrel = autocorrel;
            r.dataIdx = dataIdx;
            if (minMaxXY != null) {
                r.minMaxXY = Arrays.copyOf(minMaxXY, minMaxXY.length);
            }
            r.labels.addAll(labels);
            int x0, y0, x1, y1;
            if (offsetsToOrigCoords != null) {
                for (Entry<PairInt, PairInt> entry : offsetsToOrigCoords.entrySet()) {
                    PairInt pOffset = entry.getKey();
                    x0 = Math.round(pOffset.getX()/scale);
                    y0 = Math.round(pOffset.getY()/scale);
                    if (x0 > maxX) {
                        x0 = maxX;
                    }
                    if (y0 > maxY) {
                        y0 = maxY;
                    }
                    PairInt p = entry.getValue();
                    x1 = Math.round(p.getX()/scale);
                    y1 = Math.round(p.getY()/scale);
                    if (x1 > maxX) {
                        x1 = maxX;
                    }
                    if (y1 > maxY) {
                        y1 = maxY;
                    }
                    r.offsetsToOrigCoords.put(new PairInt(x0, y0),
                        new PairInt(x1, y1));
                    r.origPixs.add(ph.toPixelIndex(x1, y1, w2));
                }
            }
            
            return r;
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

        public int[] getMinMaxXY() {
            if (minMaxXY == null) {
                int minX = Integer.MAX_VALUE;
                int minY = Integer.MAX_VALUE;
                int maxX = Integer.MIN_VALUE;
                int maxY = Integer.MIN_VALUE;
                for (Entry<PairInt, PairInt> entry : offsetsToOrigCoords.entrySet()) {
                    PairInt xy = entry.getValue();
                    if (xy.getX() < minX) {
                        minX = xy.getX();
                    }
                    if (xy.getX() > maxX) {
                        maxX = xy.getX();
                    }
                    if (xy.getY() < minY) {
                        minY = xy.getY();
                    }
                    if (xy.getY() > maxY) {
                        maxY = xy.getY();
                    }
                }
                this.minMaxXY = new int[]{minX, maxX, minY, maxY};
            }
            return minMaxXY;
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
     * @param regions
     * @param meanWindowedImg
     * @return
     */
    public TIntObjectMap<CRegion> canonicalizeRegions(List<Region> regions,
        GreyscaleImage meanWindowedImg) {

        TIntObjectMap<CRegion> output = new TIntObjectHashMap<CRegion>();
    
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
     * each RegionPoints instance.  it also make a region for the mser
     * ellipse derived orientation.
     * 
     * @param regions
     * @return 
     */
    public TIntObjectMap<CRegion> canonicalizeRegions4(
        TIntObjectMap<RegionPoints> regions, 
        int imageWidth, int imageHeight) {
        
        int addIdx = regions.size();
        
        TIntObjectMap<CRegion> output = new TIntObjectHashMap<CRegion>();
    
        TIntObjectIterator<RegionPoints> iter = regions.iterator();
        for (int i = 0; i < regions.size(); ++i) {
            iter.advance();
            
            int rIdx = iter.key();
            RegionPoints r = iter.value();
            
            Set<CRegion> cRegions = canonicalizeRegions4(
                r, imageWidth, imageHeight);
        
            for (CRegion cRegion : cRegions) {
            
                cRegion.dataIdx = rIdx;
                
                if (output.containsKey(rIdx)) {
                    cRegion.dataIdx = addIdx;
                    output.put(addIdx, cRegion);
                    addIdx++;
                } else {
                    output.put(rIdx, cRegion);
                }
            }
        }

        return output;        
    }
    
    /**
     * NOTE: the returned CRegions need to have their .dataIdx fields set for
     * their specific context.
     * @param regionPoints
     * @return 
     */
    public List<CRegion> canonicalizeRegions(
        int w0, int h0,
        Canonicalizer.RegionPoints regionPoints) {
        
        List<CRegion> out = new ArrayList<CRegion>();
            
        TIntSet orientations = new TIntHashSet(regionPoints.hogOrientations);
        
        /*
        NOTE that hog orientations have 90 pointing up and that is the
        direction of the major axis of points, that is 90 degrees is
        the direction from x,y = (0,0) to (1,0).

        NOTE also that the regionpoint ellipse orientation is the angle of the
        minor axis of the ellipse, so 90 degrees must be subtracted from
        it to use with the dominant orientations.
        */

        int eAngle = (int)Math.round(regionPoints.ellipseParams.orientation * 180./Math.PI);
        // put into 0 to 180 ref frame
        if (eAngle > 179) {
            eAngle -= 180;
        }
        // put into ref frame of dominant orientations (major axis direction)
        eAngle -= 90;
        if (eAngle < 0) {
            eAngle += 180;
        }
        orientations.add(eAngle);
        
        TIntIterator iter2 = orientations.iterator();
        while (iter2.hasNext()) {

            int or = iter2.next();

            double angle = or * (Math.PI/180.);

            Map<PairInt, PairInt> offsetToOrigMap = 
                createOffsetToOrigMap(regionPoints.ellipseParams.xC, 
                    regionPoints.ellipseParams.yC,
                    regionPoints.accX, regionPoints.accY, w0, h0, angle);

            CRegion cRegion = new CRegion(w0, h0);
            cRegion.ellipseParams = regionPoints.ellipseParams;
            cRegion.addAllOffsets(offsetToOrigMap);
            cRegion.hogOrientation = or;
            cRegion.minMaxXY = Arrays.copyOf(regionPoints.getMinMaxXY(), 4);
            
            out.add(cRegion);
        }
        
        return out;
    }
    
    /**
     * uses RegionPoints.hogOrientations to make multiple cRegions for 
     * each RegionPoints instance.  it also make a region for the mser
     * ellipse derived orientation.
     * 
     * NOTE" remember to set cRegion.dataIdx = rIdx in the results 
     * afterwards.
     * 
     * @return 
     */
    public Set<CRegion> canonicalizeRegions4(RegionPoints r, int imageWidth, 
        int imageHeight) {
        
        Set<CRegion> out = new HashSet<CRegion>();
            
        TIntSet orientations = new TIntHashSet(r.hogOrientations);
        
        /*
        NOTE that hog orientations have 90 pointing up and that is the
        direction of the major axis of points, that is 90 degrees is
        the direction from x,y = (0,0) to (1,0).

        NOTE also that the regionpoint ellipse orientation is the angle of the
        minor axis of the ellipse, so 90 degrees must be subtracted from
        it to use with the dominant orientations.
        */

        int eAngle = (int)Math.round(r.ellipseParams.orientation * 180./Math.PI);
        // put into 0 to 180 ref frame
        if (eAngle > 179) {
            eAngle -= 180;
        }
        // put into ref frame of dominant orientations (major axis direction)
        eAngle -= 90;
        if (eAngle < 0) {
            eAngle += 180;
        }
        orientations.add(eAngle);

        TIntIterator iter2 = orientations.iterator();
        while (iter2.hasNext()) {

            int or = iter2.next();

            double angle = or * (Math.PI/180.);

            Map<PairInt, PairInt> offsetToOrigMap = createOffsetToOrigMap(
                r.ellipseParams.xC, r.ellipseParams.yC,
                r.accX, r.accY, imageWidth, imageHeight, angle);

            CRegion cRegion = new CRegion(imageWidth, imageHeight);
            cRegion.ellipseParams = r.ellipseParams;
            cRegion.addAllOffsets(offsetToOrigMap);
            cRegion.hogOrientation = or;
            
            out.add(cRegion);
        }

        return out;        
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
                r.accX, r.accY, img.getWidth(), img.getHeight(), 
                r.ellipseParams.orientation);

            CRegion cRegion = new CRegion(img.getWidth(), img.getHeight());
            cRegion.ellipseParams = r.ellipseParams;
            cRegion.addAllOffsets(offsetToOrigMap);
            
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
        int imageWidth, int imageHeight) {

        TIntObjectMap<RegionPoints> output = new TIntObjectHashMap<RegionPoints>();
    
        for (int i = 0; i < regions.size(); ++i) {

            Region r = regions.get(i);
            
            RegionPoints cRegion = canonicalizeRegion2(r, imageWidth, imageHeight);
            
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
     * 
     * @param x
     * @param y
     * @param m ellipse coefficients derived from x,y moments.
       double[]{v0x, v1x, v0y, v1y}
     * @param imageWidth
     * @param imageHeight
     * @return 
     */
    public static PairIntArray createEllipse(int x, int y, 
        double[] m, int imageWidth, int imageHeight) {

        PairIntArray xy = new PairIntArray();
        
        assert(x >= 0 && x < imageWidth);
        assert(y >= 0 && y < imageHeight);

        //v0x, v1x, v0y, v1y
        //double[] m = r.calcParamTransCoeff();

        double angle = Math.atan(m[0]/m[2]);
        if (angle < 0) {
            angle += Math.PI;
        }

        double major = 2. * m[4];
        double minor = 2. * m[5];

        double ecc = Math.sqrt(major * major - minor * minor)/major;
        assert(!Double.isNaN(ecc));

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
        
        return xy;
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

        PairIntArray xy;

        boolean createEllipse = true;
        double radius = minor;
        if (radius < 4) {
            radius = 4;
            createEllipse = false;
        }

        if (createEllipse) {
            xy = createEllipse(x, y, m, imageWidth, imageHeight);
        } else {
            xy = new PairIntArray();
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
        
        fillInEllipse(r.accX, r.accY);
        
        RegionGeometry rg = new RegionGeometry();
        rg.eccentricity = ecc;
        rg.major = major;
        rg.minor = minor;
        rg.orientation = angle;
        rg.m = m;
        rg.xC = x;
        rg.yC = y;

        RegionPoints regionPoints = new RegionPoints();
        regionPoints.ellipseParams = rg;    
        regionPoints.accX.addAll(r.accX);
        regionPoints.accY.addAll(r.accY);
        
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

        PairIntArray xy;

        boolean createEllipse = true;
        double radius = minor;
        if (radius < 4) {
            radius = 4;
            createEllipse = false;
        }

        if (createEllipse) {
            xy = createEllipse(x, y, m, imageWidth, imageHeight);
        } else {
            xy = new PairIntArray();
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
            r.accX, r.accY, imageWidth, imageHeight, angle);

        RegionGeometry rg = new RegionGeometry();
        rg.eccentricity = ecc;
        rg.major = major;
        rg.minor = minor;
        rg.orientation = angle;
        rg.xC = x;
        rg.yC = y;
        rg.m = m;

        CRegion cRegion = new CRegion(imageWidth, imageHeight);
        cRegion.ellipseParams = rg;
        cRegion.addAllOffsets(offsetToOrigMap);

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
        rg.m = m;
        
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
            
            int w2 = mImg.getWidth();
            int h2 = mImg.getHeight();

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
                    if (xScaled >= w2) {
                        xScaled = w2 - 1;
                    }
                    if (yScaled >= h2) {
                        yScaled = h2 - 1;
                    }
                    PairInt pOrigScaled = new PairInt(xScaled, yScaled);

                    // TODO: review this...should be the same as entry.getKey/scale
                    // but slightly better integer rounding results
                    double[] xyETr = transformer.applyTransformation(params,
                        xScaled, yScaled);

                    int xETr = (int)Math.round(xyETr[0]);
                    int yETr = (int)Math.round(xyETr[1]);

                    if ((xETr >= 0) && (xETr < w2)
                        && (yETr >= 0) && (yETr < h2)) {

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
                //NOTE: these are not scaled:
                rg.m = cr.ellipseParams.m;
                
                CRegion cRegion = new CRegion(w2, h2);
                cRegion.ellipseParams = rg;
                cRegion.addAllOffsets(offsetMap);
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

    public static Map<PairInt, PairInt> createOffsetToOrigMap(int x, int y, 
        TIntList xList, TIntList yList, int imgWidth, int imgHeight, 
        double orientation) {
        
        fillInEllipse(xList, yList);

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
    
        TIntList xListTr = new TIntArrayList(xList);
        TIntList yListTr = new TIntArrayList(yList);
        // ellipse, rotated by orientation to create
        //   x and y offsets from center that are comparable to
        //   the same for CRegions in other dataset.
        transformer.applyTransformation(params, xListTr, yListTr);

        // visit all points in region and transform them
        // also determine the auto-correlation

        int nc = 0;
        for (int j = 0; j < xListTr.size(); ++j) {
            
            int xp = xList.get(j);
            int yp = yList.get(j);
            
            int xpTr = xListTr.get(j);
            int ypTr = yListTr.get(j);
            
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
    
    private static void fillInEllipse(TIntList xList, TIntList yList) {

        // key = row number, value = start and stop x range
        TIntObjectMap<PairInt> rowColRange = new TIntObjectHashMap<PairInt>();
        
        int minRow = Integer.MAX_VALUE;
        int maxRow = Integer.MIN_VALUE;
        
        for (int i = 0; i < xList.size(); ++i) {
            
            int row = yList.get(i);
            int col = xList.get(i);
            
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
                xList.add(x);
                yList.add(i);
            }
        }
    }
     
    public static void filterBySpatialProximity(float critSep, 
        List<Region> regions, int width, int height) {
        
        System.out.println("before spatial filter regions.n=" + regions.size());

        // index = pixIdx, value = rIdx
        TLongIntMap pixRIdxMap = new TLongIntHashMap();
        PixelHelper ph = new PixelHelper();
       
        TIntObjectMap<Canonicalizer.RegionGeometry> rgMap 
            = new TIntObjectHashMap<Canonicalizer.RegionGeometry>();
        
        for (int rIdx = 0; rIdx < regions.size(); ++rIdx) {

            Region region = regions.get(rIdx);
         
            Canonicalizer.RegionGeometry rg = Canonicalizer.calculateEllipseParams(
                region, width, height);
            
            if (rg == null) {
                continue;
            }
            
            long pixIdx = ph.toPixelIndex(rg.xC, rg.yC, width);
            pixRIdxMap.put(pixIdx, rIdx);
            
            rgMap.put(rIdx, rg);
        }
        
        int sep = Math.round(critSep);
        
        ClusterFinder cFinder
            = new ClusterFinder(pixRIdxMap.keySet(), width, height);
        cFinder.setMinimumNumberInCluster(2);
        cFinder.setThreshholdFactor(1.f);
        cFinder.setBackgroundSeparation(sep, sep);
        cFinder.findClusters();
        List<TLongSet> groupList = cFinder.getGroups();
        
        TIntList rm = new TIntArrayList();
        
        //NOTE: may need to revise how to choose best region to keep.
        for (int i = 0; i < groupList.size(); ++i) {
            
            TLongSet groupPixs = groupList.get(i);
            
            int maxArea = Integer.MIN_VALUE;
            int maxAreaIdx = -1;
            
            TLongIterator iter3 = groupPixs.iterator();
            while (iter3.hasNext()) {
                long pixIdx = iter3.next();
                int rIdx = pixRIdxMap.get(pixIdx);
            
                Canonicalizer.RegionGeometry rg = rgMap.get(rIdx);
                if (rg == null) {
                    continue;
                }
                //double area = rg.major * rg.minor;
                int area = regions.get(rIdx).accX.size();
                
                //if (rg.xC > 73 && rg.xC < 85 && rg.yC > 49 && rg.yC < 65) {
                //    System.out.format("(%d, %d) minor=%.3f major=%.3f area=%\n", 
                //        rg.xC, rg.yC, (float)rg.minor, (float)rg.major, area);
                //}
                
                if (area > maxArea) {
                    maxArea = area;
                    maxAreaIdx = rIdx;
                }               
            }
            assert(maxAreaIdx > -1);
            iter3 = groupPixs.iterator();
            while (iter3.hasNext()) {
                long pixIdx = iter3.next();
                int rIdx = pixRIdxMap.get(pixIdx);
                if (rIdx == maxAreaIdx) {
                    continue;
                }
                rm.add(rIdx);
            }
        }
        rm.sort();
        
        for (int i = (rm.size() - 1); i > -1; --i) {
            int rmIdx = rm.get(i);
            regions.remove(rmIdx);
        }
        
        System.out.println("after spatial filter regions.n=" + regions.size());
    }

}
