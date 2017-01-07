package algorithms.imageProcessing.features.mser;

import algorithms.imageProcessing.GreyscaleImage;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import java.util.ArrayList;
import java.util.List;

/**
this package, mser, and its contents are java ports of the C++ MSER
implementation by Charles Dubout <charles.dubout@idiap.ch>
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

The MSER class extracts maximally stable extremal regions
from a grayscale (8 bits) image.
note The MSER class is not reentrant, so if you want to
extract regions in parallel, each
thread needs to have its own MSER class instance.
*/
public class MSER {

    private int delta_;
    private double minArea_;
    private double maxArea_;
    private double maxVariation_;
    private double minDiversity_;
    private boolean eight_;
    
    private List<Region> regionStack = new ArrayList<Region>();

    /**
    constructor w/ the following default values:
    @param[in] delta DELTA parameter of the MSER algorithm.
               Roughly speaking, the stability of a
	       region is the relative variation of the region
               area when the intensity is changed by delta.
               (default delta = 2.0)
    @param[in] minArea Minimum area of any stable region
               relative to the image domain area.
               (double minArea = 0.0001)
    @param[in] maxArea Maximum area of any stable region
               relative to the image domain area.
               (double maxArea = 0.5)
    @param[in] maxVariation Maximum variation (absolute
               stability score) of the regions.
               (double maxVariation = 0.5)
    @param[in] minDiversity Minimum diversity of the regions.
               When the relative area of two
     	       nested regions is below this threshold,
               then only the most stable one is selected.
               (double minDiversity = 0.33)
    @param[in] eight Use 8-connected pixels instead of 4-connected.
    */
    public MSER() {

        int delta = 2;
        double minArea = 0.0001;
        double maxArea = 0.5;
        double maxVariation = 0.5;
        double minDiversity = 0.33;
        boolean eight = false;

        init(delta, minArea, maxArea, maxVariation, minDiversity,
            eight);
    }

    /**
    @param[in] delta DELTA parameter of the MSER algorithm.
               Roughly speaking, the stability of a
	       region is the relative variation of the region
               area when the intensity is changed by delta.
               (default delta = 2.0)
    @param[in] minArea Minimum area of any stable region
               relative to the image domain area.
               (double minArea = 0.0001)
    @param[in] maxArea Maximum area of any stable region
               relative to the image domain area.
               (double maxArea = 0.5)
    @param[in] maxVariation Maximum variation (absolute
               stability score) of the regions.
               (double maxVariation = 0.5)
    @param[in] minDiversity Minimum diversity of the regions.
               When the relative area of two
	       nested regions is below this threshold,
               then only the most stable one is selected.
               (double minDiversity = 0.33)
    @param[in] eight Use 8-connected pixels instead of 4-connected.
    */
    public MSER(int delta, double minArea, double maxArea,
        double maxVariation, double minDiversity,
        boolean eight) {

        init(delta, minArea, maxArea, maxVariation, minDiversity,
            eight);
    }

    private void init(int delta, double minArea, double maxArea,
        double maxVariation, double minDiversity,
        boolean eight) {

        this.eight_ = eight;
        this.delta_ = delta;
        this.minArea_ = minArea;
        this.maxArea_ = maxArea;
        this.maxVariation_ = maxVariation;
        this.minDiversity_ = minDiversity;
       
        if (delta <= 0) {
            throw new IllegalArgumentException("delta must be > 0");
        }
        if (minArea < 0.0) {
            throw new IllegalArgumentException("minArea must be >= 0");
        }
        if (maxArea > 1.0) {
            throw new IllegalArgumentException("maxArea must be <= 1");
        }
        if (minArea >= maxArea) {
            throw new IllegalArgumentException("minArea must be < maxArea");
        }
        if (maxVariation <= 0.0) {
            throw new IllegalArgumentException("maxVariation must be > 0");
        }
        if (minDiversity < 0.0) {
            throw new IllegalArgumentException("minDiversity must be < 0");
        }
        if (minDiversity >= 1) {
            throw new IllegalArgumentException("minDiversity must be < 1");
        }
    }

    /**
      Extracts maximally stable extremal regions from a
      grayscale (8 bits) image.
      @param[in] bits Pointer to the first scanline of the image.
      @param[in] width Width of the image.
      @param[in] height Height of the image.
      @param[out] regions Detected MSER.
      //public void operator()(const uint8_t * bits,
      //int width, int height, std::vector<Region> & regions);
      */
    public void operator(int[] bits, int width, int height,
        List<Region> regions) {

        // 1. Clear the accessible pixel mask,
        //    the heap of boundary pixels and
        //    the component stack. Push
        //    a dummy-component onto the stack, with grey-level
        //    higher than any allowed in the image.
        TIntList accessible = new TIntArrayList(width * height);
        for (int i = 0; i < width * height; ++i) {
            accessible.add(0);
        }
        TIntList[] boundaryPixels = new TIntArrayList[256];
        for (int i = 0; i < 256; ++i) {
            boundaryPixels[i] = new TIntArrayList();
        }
        int priority = 256;
        List<Region> regionStack = new ArrayList<Region>();

        Region tmp = new Region();
        regionStack.add(tmp);

        // 2. Make the source pixel (with its first edge) the current pixel, mark it as accessible and
        // store the grey-level of it in the variable current level.
        int curPixel = 0;
        int curEdge = 0;
        int curLevel = bits[0];
        accessible.set(0, 1);

        // 3. Push an empty component with current level onto the component stack.
        //step_3:

        tmp = new Region(curLevel, curPixel);
        regionStack.add(tmp);

        // 4. Explore the remaining edges to the neighbors of the current pixel, in order, as follows:
        // For each neighbor, check if the neighbor is already accessible. If it is not, mark it as
        // accessible and retrieve its grey-level. If the grey-level is not lower than the current one,
        // push it onto the heap of boundary pixels. If on the other hand the grey-level is lower than
        // the current one, enter the current pixel back into the queue of boundary pixels for later
        // processing (with the next edge number), consider the new pixel and its grey-level and go to 3.
        while (true) {

            int x = curPixel % width;
            int y = curPixel / width;

            boolean s3 = false;

            for (; curEdge < (eight_ ? 8 : 4); ++curEdge) {

                int neighborPixel = curPixel;

                if (eight_) {
                    switch (curEdge) {
                        case 0: if (x < width - 1) neighborPixel = curPixel + 1; break;
                        case 1: if ((x < width - 1) && (y > 0)) neighborPixel = curPixel - width + 1; break;
                        case 2: if (y > 0) neighborPixel = curPixel - width; break;
                        case 3: if ((x > 0) && (y > 0)) neighborPixel = curPixel - width - 1; break;
                        case 4: if (x > 0) neighborPixel = curPixel - 1; break;
                        case 5: if ((x > 0) && (y < height - 1)) neighborPixel = curPixel + width - 1; break;
                        case 6: if (y < height - 1) neighborPixel = curPixel + width; break;
                        default: if ((x < width - 1) && (y < height - 1)) neighborPixel = curPixel + width + 1; break;
                    }
                } else {
                    switch (curEdge) {
                        case 0: if (x < width - 1) neighborPixel = curPixel + 1; break;
                        case 1: if (y < height - 1) neighborPixel = curPixel + width; break;
                        case 2: if (x > 0) neighborPixel = curPixel - 1; break;
                        default: if (y > 0) neighborPixel = curPixel - width; break;
                    }
                }

                if (neighborPixel != curPixel
                    && accessible.get(neighborPixel) == 0) {

                    int neighborLevel = bits[neighborPixel];
                    accessible.set(neighborPixel, 1);
                    
                    if (neighborLevel >= curLevel) {

                        boundaryPixels[neighborLevel].add(neighborPixel << 4);
                        if (neighborLevel < priority) {
                            priority = neighborLevel;
                        }
                    } else {

                        boundaryPixels[curLevel].add((curPixel << 4) | (curEdge + 1));
                        if (curLevel < priority) {
                            priority = curLevel;
                        }
                        curPixel = neighborPixel;
                        curEdge = 0;
                        curLevel = neighborLevel;

                        //continue step_3;
                        s3 = true;
                        break;
                    }
                }
            } // end for loop over currEdge

            if (s3) {

                tmp = new Region(curLevel, curPixel);
                regionStack.add(tmp);

                continue;
            }

            // 5. Accumulate the current pixel to the component at the top of the stack (water
            // saturates the current pixel).
            regionStack.get(regionStack.size() - 1).accumulate(x, y);

            // 6. Pop the heap of boundary pixels. If the heap is empty, we are done. If the returned
            // pixel is at the same grey-level as the previous, go to 4.
            if (priority == 256) {

                regionStack.get(regionStack.size() - 1)
                    .detect(delta_, (int)(minArea_ * width * height),
                        (int)(maxArea_ * width * height),
                        maxVariation_, minDiversity_, regions);

                return;
            }

            curPixel = boundaryPixels[priority]
                .get(boundaryPixels[priority].size() - 1) >> 4;
            curEdge = boundaryPixels[priority]
                .get(boundaryPixels[priority].size() - 1) & 15;

            int sz = boundaryPixels[priority].size();
            boundaryPixels[priority].removeAt(sz - 1);

            while ((priority < 256) && boundaryPixels[priority].isEmpty()) {
                ++priority;
            }

            int newPixelGreyLevel = bits[curPixel];

            if (newPixelGreyLevel != curLevel) {

                curLevel = newPixelGreyLevel;
                
                // 7. The returned pixel is at a higher grey-level, so we must now process
                // all components on the component stack until we reach the higher
                // grey-level. This is done with the processStack sub-routine, see below.
                // Then go to 4.
                processStack(newPixelGreyLevel, curPixel, regionStack);
            }
        }// end outer no arg for loop
    }

    private void processStack(int newPixelGreyLevel, int pixel,
        List<Region> regionStack) {

        // 1. Process component on the top of the stack. The next grey-level
        // is the minimum of newPixelGreyLevel and the grey-level for the
        // second component on the stack.
        do {
            
            Region top = regionStack.remove(regionStack.size() - 1);

            // 2. If newPixelGreyLevel is smaller than the grey-level on the second component on the
            // stack, set the top of stack grey-level to newPixelGreyLevel and return from sub-routine
            // (This occurs when the new pixel is at a grey-level for which there is not yet a component
            // instantiated, so we let the top of stack be that level by just changing its grey-level.
            if (newPixelGreyLevel < 
                regionStack.get(regionStack.size() - 1).level_) {

                Region tmp = new Region(newPixelGreyLevel, pixel);
                regionStack.add(tmp);
            
                regionStack.get(regionStack.size() - 1).merge(top);

                return;
            }

            // 3. Remove the top of stack and merge it into the second component
            // on stack as follows:
            // Add the first and second moment accumulators together and/or
            // join the pixel lists.
            // Either merge the histories of the components, or take the history
            // from the winner. Note
            // here that the top of stack should be considered one ’time-step’
            // back, so its current
            // size is part of the history. Therefore the top of stack would
            // be the winner if its
            // current size is larger than the previous size of second on stack.
            
            regionStack.get(regionStack.size() - 1).merge(top);
            
        } // 4. If(newPixelGreyLevel>top of stack grey-level) go to 1.
        while (newPixelGreyLevel > 
            regionStack.get(regionStack.size() - 1).level_);
    }

    /**
     * given 8 bit image, calculate the MSER regions.
     * @param img
     * @return 
     */
    public List<List<Region>> findRegions(GreyscaleImage img) {

        int width = img.getWidth();
        int height = img.getHeight();

        int[] greyscale = new int[width * height];
        for (int i = 0; i < img.getNPixels(); ++i) {
            greyscale[i] = img.getValue(i);
        }
    
        return findRegions(greyscale, width, height);
    }
    
    /**
     * given 8 bit image, calculate the MSER regions.
     * @param img
     * @return 
     */
    public List<List<Region>> findRegions(int[] greyscale, int width,
        int height) {

        int delta = 2;
        double minArea = 0.0005;
        double maxArea = 0.1;
        double maxVariation = 0.5;
        double minDiversity = 0.5;
        
        // Extract MSER
        long start = System.currentTimeMillis();

        MSER mser8 = new MSER(delta, minArea, maxArea, maxVariation, 
            minDiversity, true);
        MSER mser4 = new MSER(delta, minArea, maxArea, maxVariation, 
            minDiversity, false);

        List<List<Region>> regions = new ArrayList<List<Region>>(2);
        regions.add(new ArrayList<Region>());
        regions.add(new ArrayList<Region>());

        mser8.operator(greyscale, width, height, regions.get(0));

        // Invert the pixel values
        for (int i = 0; i < width * height; ++i) {
            greyscale[i] = ~greyscale[i];
            if (greyscale[i] < 0) {
                greyscale[i] += 256;
            }
        }

        mser4.operator(greyscale, width, height, regions.get(1));

        long stop = System.currentTimeMillis();

        System.out.println(
            "Extracted " + (regions.get(0).size() + regions.get(1).size())
            + " regions  (" + width + 'x' + height + ") in "
            + ((stop - start) / 1000) + "s.");

        return regions;
    }
}

//--------------------------------------------------------------------------------------------------
// Linear time Maximally Stable Extremal Regions implementation as described in D. Nistér and
// H. Stewénius. Linear Time Maximally Stable Extremal Regions. Proceedings of the European
// Conference on Computer Vision (ECCV), 2008.
//
// Copyright (c) 2011 Idiap Research Institute, http://www.idiap.ch/.
// Written by Charles Dubout <charles.dubout@idiap.ch>.
//
// MSER is free software: you can redistribute it and/or modify it under the terms of the GNU
// General Public License version 3 as published by the Free Software Foundation.
//
// MSER is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU General Public License along with MSER. If not, see
// <http://www.gnu.org/licenses/>.
//--------------------------------------------------------------------------------------------------

