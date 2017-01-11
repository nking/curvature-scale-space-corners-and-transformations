package algorithms.imageProcessing.features.mser;

import algorithms.QuickSort;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.util.PairInt;
import algorithms.util.TrioInt;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
MSER.java and Region.java are java ports of the C++ MSER
implementation of MSER by Charles Dubout <charles.dubout@idiap.ch>
downloaded from https://github.com/idiap/mser

The C++ code has copyright:
--------------------------
GNU GENERAL PUBLIC LICENSE, Version 3

Copyright (c) 2011 Idiap Research Institute, http://www.idiap.ch/.
Written by Charles Dubout <charles.dubout@idiap.ch>.
 
MSER is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License version 3 as published by the Free Software Foundation.
 MSER is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with MSER. If not, see
<http://www.gnu.org/licenses/>.
--------------------------

Linear time Maximally Stable Extremal Regions (MSER) implementation as described
in D. Nistér and H. Stewénius, "Linear Time Maximally Stable Extremal Regions",
ECCV 2008.
The functionality is similar to that of VLFeat MSER feature detector
<http://www.vlfeat.org/overview/mser.html> but the code is several time faster.
MSER is a blob detector, like the Laplacian of Gaussian used by the SIFT
algorithm. It extracts stable connected regions of some level sets from an
image, and optionally fits ellipses to them.
* 
* ------
* author nichole ported the C++ code quoted above to java and added
* methods used in canonicalization (in progress).
*/
public class Region {

    /**
     * Level at which the region is processed.
     */
    public int level_;

    /**
     * Index of the initial pixel (y * width + x).
     */
    int pixel_;

    /**
     * Area of the region (moment zero), that is, the number of
     * pixels in this region.
     */
    int area_;
    
    /**
     * First and second moments of the region (x, y, x^2, xy, y^2).
     */
    double[] moments_;

    /**
     * MSER variation.
     */
    double variation_;

    /**
     * Flag indicating if the region is stable
     */
    private boolean stable_;

    private Region parent_ = null;
    private Region child_ = null;
    private Region next_ = null;
    
    /**
     * constructor with default level = 256 and pixel=0
     */
    public Region() {
        int level = 256;
        int pixel = 0;
        init(level, pixel);
    }

    /**
     * @param level Level at which the region is processed (default level = 256)
     * @param pixel Index of the initial pixel (y * width + x). (int pixel = 0)
     */
    public Region(int level, int pixel) {
        init(level, pixel);
    }

    private void init(int level, int pixel) {
        this.level_ = level;
        this.pixel_ = pixel;
        this.area_ = 0;

        this.moments_ = new double[5];

        this.variation_ = Double.POSITIVE_INFINITY;
        this.stable_ = false;

        this.parent_ = null;
        this.child_ = null;
        this.next_ = null;
        
    }

    public void accumulate(int x, int y) {
        ++area_;
        moments_[0] += x;
        moments_[1] += y;
        moments_[2] += x * x;
        moments_[3] += x * y;
        moments_[4] += y * y;
        
    }

    public void merge(Region child) {
        assert(child.parent_ == null);
        assert(child.next_ == null);

        // Add the moments together
        area_ += child.area_;
        moments_[0] += child.moments_[0];
        moments_[1] += child.moments_[1];
        moments_[2] += child.moments_[2];
        moments_[3] += child.moments_[3];
        moments_[4] += child.moments_[4];
        
        child.next_ = child_;
        child_ = child;
        child.parent_ = this;
    }

    public void detect(int delta, int minArea, int maxArea,
        double maxVariation, double minDiversity,
        List<Region> regions) {

        process(delta, minArea, maxArea, maxVariation);

        save(minDiversity, regions);
    }

    public void process(int delta, int minArea, int maxArea,
        double maxVariation) {

        // Find the last parent with level not higher than level + delta
        Region parent = this;

        while (parent.parent_ != null
            && (parent.parent_.level_ <= (level_ + delta))) {

            parent = parent.parent_;
        }

        // Calculate variation
        variation_ = (double)(parent.area_ - area_) / (double)area_;

        // Whether or not the region *could* be stable
        boolean stable
            = ((parent_ == null) || (variation_ <= parent_.variation_))
            && (area_ >= minArea) && (area_ <= maxArea)
            && (variation_ <= maxVariation);

        // Process all the children
        for (Region child = child_; child != null; child = child.next_) {

            child.process(delta, minArea, maxArea, maxVariation);

            if (stable && (variation_ < child.variation_)) {
                stable_ = true;
            }
        }

        // The region can be stable even without any children
        if ((child_ == null) && stable) {
            stable_ = true;
        }
    }

    public boolean check(double variation, int area) {

        if (area_ <= area) {
            return true;
        }

        if (stable_ && (variation_ < variation)) {
            return false;
        }

        for (Region child = child_; child != null; child = child.next_) {
            if (!child.check(variation, area)) {
                return false;
            }
        }

        return true;
    }

    public void save(double minDiversity, List<Region> regions) {

        int minParentArea = 0;
        if (stable_) {
            minParentArea = (int)(area_ / (1.0 - minDiversity) + 0.5);
        }

        Region parent = this;

        while ((parent.parent_ != null)
            && (parent.parent_.area_ < minParentArea)) {

            parent = parent.parent_;

            if (parent.stable_
                && (parent.variation_ <= variation_)) {
                stable_ = false;
                break;
            }
        }

        if (stable_) {
            int maxChildArea = (int)(area_ * (1.0 - minDiversity) + 0.5);

            if (!check(variation_, maxChildArea)) {
                stable_ = false;
            }
        }

        if (stable_) {
            regions.add(this.copy());
            Region last = regions.get(regions.size() - 1);
            last.parent_ = null;
            last.child_ = null;
            last.next_ = null;
        }

        for (Region child = child_; child != null; child  = child.next_) {
            child.save(minDiversity, regions);
        }
    }
    
    public static void drawEllipses(Image img, List<List<Region>> regions, int nExtraDot) {
        
        if (regions.size() != 2) {
            throw new IllegalArgumentException("expecting the bright then"
                + " dark msers in a list of size 2.");
        }
        
        for (int i = 0; i < regions.get(0).size(); ++i) {
            regions.get(0).get(i).drawEllipse(img, nExtraDot, 
                //127, 127, 127);
                255, 255, 255);
        }
        
        for (int i = 0; i < regions.get(1).size(); ++i) {
            regions.get(1).get(i).drawEllipse(img, nExtraDot, 
                255, 255, 255);
        }
    }
    
    /**
     * calculate the x,y centroid of the pixels within
     * this region.
     * @param imageWidth
     * @param imageHeight
     * @param outputXY 
     */
    public void calculateXYCentroid(int[] outputXY, int imageWidth, int imageHeight) {
                    
        outputXY[0] = (int)Math.round(moments_[0]/(double)area_);
        outputXY[1] = (int)Math.round(moments_[1]/(double)area_);
        if (outputXY[0] == -1) {
            outputXY[0] = 0;
        }
        if (outputXY[1] == -1) {
            outputXY[1] = 0;
        }
        if (outputXY[0] == imageWidth) {
            outputXY[0] = imageWidth - 1;
        }
        if (outputXY[1] == imageHeight) {
            outputXY[1] = imageHeight - 1;
        }
    }
    
    /**
     * calculate the intensity weighted centroid of the pixels within
     * this region.
     * @param greyscale
     * @param imageWidth
     * @param imageHeight
     * @param outputXY 
     */
    void calculateIntensityCentroid(int[] greyscale, int imageWidth, 
        int imageHeight, int[] outputXY) {
        
        TrioInt[] yXRanges = getEllipseRange(imageWidth, imageHeight);
        
        double iSum = 0;
        
        Set<PairInt> visited = new HashSet<PairInt>();
        
        for (TrioInt yXminMax : yXRanges) {
            int y = yXminMax.getX();
            for (int x = yXminMax.getY(); x <= yXminMax.getZ(); ++x) {
                PairInt p = new PairInt(x, y);
                if (visited.contains(p)) {
                    continue;
                }
                visited.add(p);
                int idx = (y * imageWidth) + x;
                int intensity = greyscale[idx];
                iSum += intensity;
            }
        }
        
        double xSum = 0;
        double ySum = 0;
        
        // for debugging only, assert sum of weights
        double wSum = 0;
        
        visited.clear();
        
        for (TrioInt yXminMax : yXRanges) {
            
            int y = yXminMax.getX();
            
            for (int x = yXminMax.getY(); x <= yXminMax.getZ(); ++x) {
                
                PairInt p = new PairInt(x, y);
                if (visited.contains(p)) {
                    continue;
                }
                visited.add(p);
                
                int idx = (y * imageWidth) + x;
                
                double intensity = greyscale[idx];
                double w = (intensity/iSum);
                
                xSum += (x * w);
                ySum += (y * w);
                
                wSum += w;
            }
        }
        assert(Math.abs(wSum - 1) < 0.01);
                
        outputXY[0] = (int)Math.round(xSum);
        outputXY[1] = (int)Math.round(ySum);
    }
    
    /**
     * calculate the intensity weighted centroid of pixels within
     * radius distance of center of region.
     * 
     * @param greyscale
     * @param imageWidth
     * @param imageHeight
     * @param outputXY
     * @param radius 
     */
    void calculateIntensityCentroid(int[] greyscale, int imageWidth, 
        int imageHeight, int[] outputXY, int radius) {
      
        int xC = (int)Math.round(moments_[0]/(double)area_);
        int yC = (int)Math.round(moments_[1]/(double)area_);
        
        //v0x, v1x, v0y, v1y
        double[] m = calcParamTransCoeff();
        
        // semi-major and semi-minor axes:
        double major = 2. * m[4];
        double minor = 2. * m[5];
            
        if (radius > major) {
            radius = (int)major; 
        }
        if (radius > minor) {
            radius = (int)minor; 
        }
        
        int yStart = yC - radius;
        if (yStart < 0) {
            yStart = 0;
        }
        int yStop = yC + radius;
        if (yStop >= imageHeight) {
            yStop = imageHeight - 1;
        }
        
        /*
        x^2 + y^2 = radius^2
        x^2 =  radius^2 - y^2
        (x - xC) = sqrt(radius^2 - (y - yc)^2)
        */
        double rSq = radius * radius;
        
        double iSum = 0;
        
        // for uniqueness
        Set<PairInt> visited = new HashSet<PairInt>();
        
        for (int y = yStart; y <= yStop; ++y) {
            
            double ySq = y - yC;
            ySq *= ySq;
            
            double xH = Math.sqrt(rSq - ySq);
            
            int xStart = (int)Math.round(xC - xH);
            if (xStart < 0) {
                xStart = 0;
            }
            int xStop = (int)Math.round(xC + xH);
            if (xStop >= imageWidth) {
                xStop = imageWidth - 1;
            }
            
            for (int x = xStart; x <= xStop; ++x) {
                
                PairInt p = new PairInt(x, y);
                if (visited.contains(p)) {
                    continue;
                }
                visited.add(p);
                
                int idx = (y * imageWidth) + x;
                iSum += greyscale[idx];
            }
        }
        
        double xSum = 0;
        double ySum = 0;
        
        // for debugging only, assert sum of weights
        double wSum = 0;
        
        visited.clear();
        
        for (int y = yStart; y <= yStop; ++y) {
            
            double ySq = y - yC;
            ySq *= ySq;
            
            double xH = Math.sqrt(rSq - ySq);
            
            int xStart = (int)Math.round(xC - xH);
            if (xStart < 0) {
                xStart = 0;
            }
            int xStop = (int)Math.round(xC + xH);
            if (xStop >= imageWidth) {
                xStop = imageWidth - 1;
            }
            
            for (int x = xStart; x <= xStop; ++x) {
                
                PairInt p = new PairInt(x, y);
                if (visited.contains(p)) {
                    continue;
                }
                visited.add(p);
                
                int idx = (y * imageWidth) + x;
                double intensity = greyscale[idx];
                double w = (intensity/iSum);
                
                xSum += (x * w);
                ySum += (y * w);
                
                wSum += w;
            }
        }
        assert(Math.abs(wSum - 1) < 0.01);
                
        outputXY[0] = (int)Math.round(xSum);
        outputXY[1] = (int)Math.round(ySum);
    }
    
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
     * @return 
     */
    public double[] calcParamTransCoeff() {
        
        /*
        moments_[0]  x;
        moments_[1]  y;
        moments_[2]  x * x;
        moments_[3]  x * y;
        moments_[4]  y * y;
        */
        
        // Centroid (mean)
        double x = moments_[0] / (double)area_;
        double y = moments_[1] / (double)area_;
        
        // a*x^2 + b*x*y + c*y^2 + d*x + e*y + f = 0
        
        // Covariance matrix [a b; b c]
        double a = moments_[2] / (double)area_ - x * x;
        double b = moments_[3] / (double)area_ - x * y;
        double c = moments_[4] / (double)area_ - y * y;

        /*
        
        covariance matrix as a real 2X2 matrix 
        xx  xy
        yx  yy
        
        from Strang "Linear Algebra", chap 10, have for real 2X2 matrix:
        
        | xx-eigenv   xy        |
        | yx          yy-eigenv | = eigenv^2 - (xx + yy)*eigenv 
                                      + ((xx yy) - (yx xy)) = 0 
        
        can solve for the zeroes with quardratic equation to
        get 2 eigenvalues:
        
        eigenv = (xx + yy +- sqrt( (xx + yy)^2 -4*(xx yy - xy yx) )) / 2
               = (d +- sqrt(d*d - 4*b*b))/2
        
        Looks like Dubout below uses orthogonal axes, replacing
        (xx, yy) with (yy, -xx).
        */
        // Eigenvalues of the covariance matrix
        double d  = a + c; // xx + yy
        double e  = a - c; // xx - yy
        double f  = Math.sqrt(4.0 * b * b + e * e);
        double e0 = (d + f) / 2.0; // First eigenvalue
        double e1 = (d - f) / 2.0; // Second eigenvalue

        // Desired norm of the eigenvectors
        double e0sq = Math.sqrt(e0);
        double e1sq = Math.sqrt(e1);
 
        // Eigenvectors
        double v0x = e0sq;
        double v0y = 0.0;
        double v1x = 0.0;
        double v1y = e1sq;

        if (b != 0.) {
            v0x = e0 - c;
            v0y = b;
            v1x = e1 - c;
            v1y = b;

            // Normalize the eigenvectors
            double n0 = e0sq / Math.sqrt(v0x * v0x + v0y * v0y);
            v0x *= n0;
            v0y *= n0;

            double n1 = e1sq / Math.sqrt(v1x * v1x + v1y * v1y);
            v1x *= n1;
            v1y *= n1;
        }
        
        /*
        rotation transformed ellipse:
        
        x(t) = xCenter + aParam*cos(alpha)*cos(t) − bParam*sin(alpha)*sin(t)
        y(t) = yCenter + aParam*sin(alpha)*cos(t) + bParam*cos(alpha)*sin(t)
        
        v0x = aParam * cos(alpha);  alpha = atan(v0x/v0y)
        v1x = bParam * sin(alpha);
        v0y = aParam * sin(alpha);
        v1y = bParam * cos(alpha);
        
        so the semi-major axis length = 2 * Math.max(v0x, v1x)
        
        and the semi-minor axis length = 2 * Math.min(v0y, v1y)
        */
        
        return new double[]{v0x, v1x, v0y, v1y, e0sq, e1sq};
    }
    
    /**
     * following the creation of an ellipse from the moments and
     * area (see drawEllipse), this method orders the ellipse points by
     * y, then x then aggregates the x to form the range xmin, xmax
     * per row of y.  The return is the ranges ordered by increasing
     * y in formation [y, xmin for row y, xmax for row y]
     * @return 
     */
    public TrioInt[] getEllipseRange(int imageWidth, int imageHeight) {
        
        TIntList xs = new TIntArrayList();
        TIntList ys = new TIntArrayList();
    
        // Centroid (mean)
        double x = moments_[0] / (double)area_;
        double y = moments_[1] / (double)area_;
        
        //v0x, v1x, v0y, v1y
        double[] coeffs = calcParamTransCoeff();
        
        double v0x = coeffs[0];
        double v1x = coeffs[1];
        double v0y = coeffs[2];
        double v1y = coeffs[3];

        Set<PairInt> visited = new HashSet<PairInt>();

        for (double t = 0.0; t < 2.0 * Math.PI; t += 0.001) {
            
            double mc = Math.cos(t);
            double ms = Math.sin(t);
            
            int x2 = (int)Math.round(x + 
                (mc * v0x + ms * v1x) * 2.0 + 0.5);
            int y2 = (int)Math.round(y + (mc * v0y 
                + ms * v1y) * 2.0 + 0.5);

            if ((x2 >= 0) && (x2 < imageWidth) 
                && (y2 >= 0) && (y2 < imageHeight)) {
                
                PairInt pt = new PairInt(x2, y2);
                if (visited.contains(pt)) {
                    continue;
                }
                visited.add(pt);
                
                xs.add(x2);
                ys.add(y2);
            }
        }
        
        int[] x2s = xs.toArray(new int[xs.size()]);
        int[] y2s = ys.toArray(new int[xs.size()]);
        
        QuickSort.sortBy1stThen2nd(y2s, x2s);
        
        int n = (new TIntHashSet(ys)).size();
        
        TrioInt[] output = new TrioInt[n];
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
            output[count] = new TrioInt(yC, minX, maxX);
            count++;
        }
        assert(count == n);
        
        return output;
    }
    
    public void drawEllipse(Image img, int nExtraDot,
        int rClr, int gClr, int bClr) {

        // Centroid (mean)
        double x = moments_[0] / (double)area_;
        double y = moments_[1] / (double)area_;
        
        //v0x, v1x, v0y, v1y
        double[] coeffs = calcParamTransCoeff();
        
        double v0x = coeffs[0];
        double v1x = coeffs[1];
        double v0y = coeffs[2];
        double v1y = coeffs[3];
        
        for (double t = 0.0; t < 2.0 * Math.PI; t += 0.001) {
            
            double mc = Math.cos(t);
            double ms = Math.sin(t);
            
            int x2 = (int)Math.round(x + (mc * v0x + ms * v1x) * 2.0 + 0.5);
            int y2 = (int)Math.round(y + (mc * v0y + ms * v1y) * 2.0 + 0.5);

            if ((x2 >= 0) && (x2 < img.getWidth()) 
                && (y2 >= 0) && (y2 < img.getHeight())) {
                
                ImageIOHelper.addPointToImage(x2, y2, img, 
                    nExtraDot, rClr, gClr, bClr);
            }
        }
    }

    public void drawCircle(Image img, int nExtraDot,
        int rClr, int gClr, int bClr) {

        // Centroid (mean)
        double x = moments_[0] / (double)area_;
        double y = moments_[1] / (double)area_;
        
        //v0x, v1x, v0y, v1y
        double[] m = calcParamTransCoeff();
        
        double major = 2. * m[4];
        double minor = 2. * m[5];
        double radius = minor;
        if (radius < 4) {
            radius = 4;
        }
            
        for (double t = 0.0; t < 2.0 * Math.PI; t += 0.001) {
            
            double mc = Math.cos(t);
            double ms = Math.sin(t);
            
            int x2 = (int)Math.round(x + (mc * radius));
            int y2 = (int)Math.round(y + (ms * radius));

            if ((x2 >= 0) && (x2 < img.getWidth()) 
                && (y2 >= 0) && (y2 < img.getHeight())) {
                
                ImageIOHelper.addPointToImage(x2, y2, img, 
                    nExtraDot, rClr, gClr, bClr);
            }
        }
    }
    
    public Region copy() {
        Region region = new Region(this.level_, this.pixel_);
        region.area_ = this.area_;
        System.arraycopy(moments_, 0,region. moments_, 0, moments_.length);
        region.variation_ = this.variation_;
        region.stable_ = this.stable_;
        region.parent_ = this.parent_;
        region.child_ = this.child_;
        region.next_ = this.next_;
        return region;
    }
    
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("level_=").append(level_);
        sb.append(" pixel_=").append(pixel_);
        sb.append(" area_=").append(area_);
        sb.append(" moments_=").append(Arrays.toString(moments_));
        sb.append(" variation_=").append(variation_);
        return sb.toString();
    }
    
}
