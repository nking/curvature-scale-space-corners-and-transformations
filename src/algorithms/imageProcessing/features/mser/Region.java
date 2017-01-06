package algorithms.imageProcessing.features.mser;

import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import java.util.Arrays;
import java.util.List;

/**
this package, mser, and its contents are java ports of the C++ MSER
implementation by 
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
*/
class Region {

    /**
     * Level at which the region is processed.
     */
    public int level_;

    /**
     * Index of the initial pixel (y * width + x).
     */
    int pixel_;

    /**
     * Area of the region (moment zero).
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
     * @param level Level at which the region is processed (default level = 256)
     * @param pixel Index of the initial pixel (y * width + x). (int pixel = 0)
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
        variation_ = (parent.area_ - area_) / area_;

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
            minParentArea = (int) Math.round(
                area_ / (1.0 - minDiversity) + 0.5);
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
            int maxChildArea = (int) Math.round(
                area_ * (1.0 - minDiversity) + 0.5);

            if (!check(variation_, maxChildArea)) {
                stable_ = false;
            }
        }

        if (stable_) {
            regions.add(this);
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
            regions.get(0).get(i).drawEllipse(img, nExtraDot, 127, 127, 127);
        }
        
        for (int i = 0; i < regions.get(1).size(); ++i) {
            regions.get(1).get(i).drawEllipse(img, nExtraDot, 255, 255, 255);
        }
    }
    
    public void drawEllipse(Image img, int nExtraDot,
        int rClr, int gClr, int bClr) {

        // Centroid (mean)
        double x = moments_[0] / area_;
        double y = moments_[1] / area_;

        // Covariance matrix [a b; b c]
        double a = moments_[2] / area_ - x * x;
        double b = moments_[3] / area_ - x * y;
        double c = moments_[4] / area_ - y * y;

        // Eigenvalues of the covariance matrix
        double d  = a + c;
        double e  = a - c;
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

        for (double t = 0.0; t < 2.0 * Math.PI; t += 0.001) {
            
            int x2 = (int)Math.round(x + 
                (Math.cos(t) * v0x + Math.sin(t) * v1x) * 2.0 + 0.5);
            int y2 = (int)Math.round(y + (Math.cos(t) * v0y 
                + Math.sin(t) * v1y) * 2.0 + 0.5);

            if ((x2 >= 0) && (x2 < img.getWidth()) 
                && (y2 >= 0) && (y2 < img.getHeight())) {
                
                ImageIOHelper.addPointToImage(x2, y2, img, 
                    nExtraDot, rClr, gClr, bClr);
            }
        }
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
