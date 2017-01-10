package algorithms.imageProcessing.features.mser;

import algorithms.QuickSort;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.util.PairInt;
import java.util.List;
import algorithms.util.PairIntArray;
import algorithms.util.TrioInt;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
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

    then the radius of the circle, before scaling would
    
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
     * NOTE, for best use, invoker should use this descriptor with
     * an image processed to create a window average for each pixel.
       For example:
       SummedAreaTable sumTable = new SummedAreaTable();
       GreyscaleImage imgM = sumTable.createAbsoluteSummedAreaTable(img);
       imgM = sumTable.applyMeanOfWindowFromSummedAreaTable(imgM,
            2*halfDimension + 1);
       
     * @return 
     */
    public List<CRegion> canonicalizeRegions(List<Region> regions,
        GreyscaleImage meanWindowedImg) {
        
        int imageWidth = meanWindowedImg.getWidth(); 
        int imageHeight = meanWindowedImg.getHeight();
            
        List<CRegion> output = new ArrayList<CRegion>();
       
        Transformer transformer = new Transformer();
        
        int[] xyCen = new int[2];
        
        for (int i = 0; i < regions.size(); ++i) {
         
            Region r = regions.get(i);
            
            r.calculateXYCentroid(xyCen, imageWidth, imageHeight, xyCen);
            int x = xyCen[0];
            int y = xyCen[1];

            //v0x, v1x, v0y, v1y
            double[] m = r.calcParamTransCoeff();

            double angle = Math.atan(m[0]/m[2]);
            if (angle < 0) {
                angle += Math.PI;
            }

            double major = 2. * Math.max(Math.abs(m[0]), Math.abs(m[1]));
            double minor = 2. * Math.min(Math.abs(m[2]), Math.abs(m[3]));

            double ecc = Math.sqrt(major * major - minor * minor)/major;
        
            TIntList xs = new TIntArrayList();
            TIntList ys = new TIntArrayList();
    
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
                (float)(angle*180./Math.PI)));
            
            Map<PairInt, PairInt> offsetToOrigMap = new HashMap<PairInt, PairInt>();
            
            // visit all points in region and transform them
            // also determine the auto-correlation
            
            int vc = meanWindowedImg.getValue(x, y);
            int nc = 0;
            double autocorSum = 0;
            
            for (int j = 0; j < ranges.length; ++j) {
                int yE = ranges[j].getX();
                int xMin = ranges[j].getY();
                int xMax = ranges[j].getZ();
                for (int xE = xMin; xE <= xMax; ++xE) {
                    
                    double[] xyETr = transformer.applyTransformation(params, xE, yE);
                    
                    if ((xyETr[0] >= 0) && (xyETr[0] < imageWidth) 
                        && (xyETr[1] >= 0) && (xyETr[1] < imageHeight)) {

                        int xETr = (int)Math.round(xyETr[0]);
                        int yETr = (int)Math.round(xyETr[1]);
                    
                        PairInt pt = new PairInt(xETr, yETr);
                        if (visited.contains(pt)) {
                            continue;
                        }
                        visited.add(pt);
                                       
                        PairInt pOffsets = new PairInt(xETr - x, yETr - y);
                        offsetToOrigMap.put(pOffsets, new PairInt(xE, yE));
                
                        int diff = meanWindowedImg.getValue(xE, yE) - vc;
                        
                        autocorSum += (diff * diff);
                        nc++;
                    } 
                }
            } 
            autocorSum /= (double)nc;
            
            CRegion cRegion = new CRegion();
            cRegion.eccentricity = ecc;
            cRegion.offsetsToOrigCoords = offsetToOrigMap;
            cRegion.major = major;
            cRegion.minor = minor;
            cRegion.orientation = angle;
            cRegion.xC = x;
            cRegion.yC = y;
            cRegion.nTrEllipsePixels = visited.size();
            cRegion.autocorrel = Math.sqrt(autocorSum)/255.;
            
            output.add(cRegion);
        }
        
        return output;
    }
    
    public PairIntArray calculateCentroids(List<Region> regions,
        int[] greyscale, int imageWidth, int imageHeight) {
        
        //return calculateIntensityCentroids(regions,
        //    greyscale, imageWidth, imageHeight);
        
        return extractRegionXYCenters(regions,
            greyscale, imageWidth, imageHeight);
    }
    
    PairIntArray extractRegionXYCenters(List<Region> regions,
        int[] greyscale, int imageWidth, int imageHeight) {
        
        int[] xyCen = new int[2];
        
        PairIntArray output = new PairIntArray(regions.size());
        
        for (int i = 0; i < regions.size(); ++i) {
            Region r = regions.get(i);
            r.calculateXYCentroid(greyscale, imageWidth, imageHeight, 
                xyCen);
     
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
}
