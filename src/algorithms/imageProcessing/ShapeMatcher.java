package algorithms.imageProcessing;

import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.misc.Histogram;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayDescendingComparator;
import algorithms.util.ResourceFinder;
import java.awt.Polygon;
import java.awt.geom.PathIterator;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Logger;

public class ShapeMatcher {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public ShapeMatcher() {
    }

    /**
    method to extract general shapes from the images and compare them in order to 
    match points.  It returns a fit to a rough Euclidean transformation.
    NOTE that the images may need pre-processing steps before using this.  For example,
    the Brown & Lowe 200? panoramic images of a mountain need to have the sky masked
    out of the image first.
     * @param image1
     * @param image2
     * @param outputMatched1
     * @param outputMatched2
    */
    public TransformationPointFit findMatchingShapes(ImageExt image1, ImageExt image2,
    PairIntArray outputMatched1, PairIntArray outputMatched2) throws 
        IOException, NoSuchAlgorithmException {
        
        GreyscaleImage img1Grey = image1.copyToGreyscale();
        GreyscaleImage img2Grey = image2.copyToGreyscale();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        boolean performBinning = false;
        int binFactor1 = 1;
        
        /*
        one could start with essentially no limits here and then
        looks at the distribution of resulting contiguous group
        sizes to decide the range to keep.
        For now, choosing limits.
        */
        int smallestGroupLimit = 50;
        int largestGroupLimit = 2000;
        
        if (performBinning) {
            binFactor1 = (int) Math.ceil(
                Math.max((float)img1Grey.getWidth()/200.f,
                (float)img2Grey.getHeight()/200.));
            smallestGroupLimit /= binFactor1;
            largestGroupLimit /= binFactor1;
            
            // prevent from being smaller than needed for a convex hull
            if (smallestGroupLimit <= 3) {
                smallestGroupLimit = 4;
            }
            img1Grey = imageProcessor.binImage(img1Grey, binFactor1);
            img2Grey = imageProcessor.binImage(img2Grey, binFactor1);
        }
        imageProcessor.applyImageSegmentation(img1Grey, 3);
        imageProcessor.applyImageSegmentation(img2Grey, 3);
           
        Map<Integer, Integer> freqMap1 = Histogram.createAFrequencyMap(img1Grey);
        Map<Integer, Integer> freqMap2 = Histogram.createAFrequencyMap(img2Grey);
        
        Map<Integer, List<PairIntArray>> contigMap1 
            = new HashMap<Integer, List<PairIntArray>>();
        Map<Integer, List<PairIntArray>> contigMap2 
            = new HashMap<Integer, List<PairIntArray>>();
        
        Map<Integer, List<GrahamScan>> hulls1 = 
            new HashMap<Integer, List<GrahamScan>>();
        Map<Integer, List<GrahamScan>> hulls2 = 
            new HashMap<Integer, List<GrahamScan>>();
        
        for (int im = 0; im < 2; ++im) {
            
            Map<Integer, Integer> freqMap = freqMap1;
            Map<Integer, List<PairIntArray>> contigMap = contigMap1;
            Map<Integer, List<GrahamScan>> hulls = hulls1;
            GreyscaleImage imgGrey = img1Grey;
            
            if (im == 1) {
                freqMap = freqMap2;
                contigMap = contigMap2;
                hulls = hulls2;
                imgGrey = img2Grey;
            }
 Image imgW = ImageIOHelper.convertImage(imgGrey);
 int c = 0;
 
            for (Entry<Integer, Integer> entry : freqMap.entrySet()) {

                Integer pixValue = entry.getKey();

                DFSContiguousValueFinder finder = new DFSContiguousValueFinder(
                    imgGrey);
                finder.setMinimumNumberInCluster(smallestGroupLimit);
                finder.findGroups(pixValue.intValue());

                int nGroups = finder.getNumberOfGroups();
                List<PairIntArray> list = new ArrayList<PairIntArray>();
                for (int i = 0; i < nGroups; i++) {
                    PairIntArray xy = finder.getXY(i);
                    if (xy.getN() < largestGroupLimit) {
                        list.add(xy);
if (c == 0) {
ImageIOHelper.addCurveToImage(xy, imgW, 0, 255, 0, 0);
} else if (c == 1) {
ImageIOHelper.addCurveToImage(xy, imgW, 0, 200, 255, 150);
} else  {
ImageIOHelper.addCurveToImage(xy, imgW, 0, 0, 0, 255);
}
                    }
                }
                Collections.sort(list, new PairIntArrayDescendingComparator());

                contigMap.put(pixValue, list);
                
                List<GrahamScan> listHulls = new ArrayList<GrahamScan>();
                for (PairIntArray xy : list) {
                    GrahamScan scan = new GrahamScan();
                    try {
                        float[] x = new float[xy.getN()];
                        float[] y = new float[x.length];
                        for (int i = 0; i < x.length; ++i) {
                            x[i] = xy.getX(i);
                            y[i] = xy.getY(i);
                        }

                        scan.computeHull(x, y);

                        listHulls.add(scan);
                    } catch (GrahamScanTooFewPointsException e) {
                        log.severe(e.getMessage());
                    }                    
                }
                hulls.put(pixValue, listHulls);
c++;
            }
 try {
    String dirPath = ResourceFinder.findDirectory("bin");
    ImageIOHelper.writeOutputImage(dirPath + "/img" + im + ".png", imgW);
} catch (Exception e) {
     e.printStackTrace();
    log.severe("ERROR: " + e.getMessage());
}
 
        }
        
        Image img1W = ImageIOHelper.convertImage(img1Grey);
        Image img2W = ImageIOHelper.convertImage(img2Grey);
        
        int c = 0;
        for (Entry<Integer, List<GrahamScan>> entry : hulls1.entrySet()) {
            List<GrahamScan> hulls = entry.getValue();
            for (GrahamScan hull : hulls) {
                int[] x = new int[hull.getXHull().length];
                int[] y = new int[x.length];
                for (int i = 0; i < x.length; ++i) {
                    x[i] = Math.round(hull.getXHull()[i]);
                    y[i] = Math.round(hull.getYHull()[i]);                   
                }
                if (c == 0) {
                    //ImageIOHelper.addCurveToImage(x, y, img1W, 1, 255, 0, 0);
                    ImageIOHelper.drawLinesInImage(x, y, img1W, 1, 255, 0, 0);
                } else if (c == 1) {
                    //ImageIOHelper.addCurveToImage(x, y, img1W, 1, 0, 255, 0);
                    ImageIOHelper.drawLinesInImage(x, y, img1W, 1, 0, 255, 0);
                } else {
                    //ImageIOHelper.addCurveToImage(x, y, img1W, 1, 0, 0, 255);
                    ImageIOHelper.drawLinesInImage(x, y, img1W, 1, 0, 0, 255);
                }
            }
            c++;
        }
        c = 0;
        for (Entry<Integer, List<GrahamScan>> entry : hulls2.entrySet()) {
            List<GrahamScan> hulls = entry.getValue();
            for (GrahamScan hull : hulls) {
                int[] x = new int[hull.getXHull().length];
                int[] y = new int[x.length];
                for (int i = 0; i < x.length; ++i) {
                    x[i] = Math.round(hull.getXHull()[i]);
                    y[i] = Math.round(hull.getYHull()[i]);
                }
                if (c == 0) {
                    //ImageIOHelper.addCurveToImage(x, y, img2W, 1, 255, 0, 0);
                    ImageIOHelper.drawLinesInImage(x, y, img2W, 1, 255, 0, 0);
                } else if (c == 1) {
                    //ImageIOHelper.addCurveToImage(x, y, img2W, 1, 0, 255, 0);
                    ImageIOHelper.drawLinesInImage(x, y, img2W, 1, 0, 255, 0);
                } else {
                    //ImageIOHelper.addCurveToImage(x, y, img2W, 1, 0, 0, 255);
                    ImageIOHelper.drawLinesInImage(x, y, img2W, 1, 0, 0, 255);
                }
            }
            c++;
        }
       
        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/img1_binned.png", img1W);
            ImageIOHelper.writeOutputImage(dirPath + "/img2_binned.png", img2W);
        } catch (Exception e) {
             e.printStackTrace();
            log.severe("ERROR: " + e.getMessage());
        }
        return null;
        //throw new UnsupportedOperationException("not yet implemented");
    
        /*
        -- segmentation by cieXY color into 3 bands
        -- dfs contiguous find of each of the 3 groups of sizes > <100?>
           for the largest groups:
              -- convex hull to make shapes.
              -- compare the hulls to the other image:  
                 -- by area and single pixel intensity?  (the cie xy histogram skipping
                    algorithm may assign different colors to same feature, so should
                    probably ignore intensity unless know the images are very similar).
                 -- by shape
                    -- furthest pairs?
                    -- a description of ellipticity?
                    -- points on the hull or centroid matching enough between 
                       the 2 images to be used like "corners"?
        */
    }

}
