package algorithms.imageProcessing.features.mser;

import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.DFSConnectedGroupsFinder;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.features.mser.MSER.Threshold;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * class to explore the boundaries of the accumulated points in MSER regions.
 * 
 * The contents may in the future be replaced by an implementation of 
 * http://www.vision.ee.ethz.ch/~rhayko/paper/aapr2009_boundary_detection_SBER_hayko.pdf
 * 
 * but for now is an exploration of the existing MSER and Region class results.
 * 
 * @author nichole
 */
public class MSEREdges {
    
    /*
    given a color image
      -- makes greyscale image and a polar theta of cie luv
      -- creates positive and negative mser regions for both,
            filtering by variation for the specific mser type
      -- applies a spatial filter to keep one for an overlapping
         or close centered mser.
         -- keeps the mser w/ most accumulated points
            (this is where the component tree approach would be better).
      -- using PerimeterFinder2 extract the outer boundaries from each mser
         region points, and do the same for the contiguous embedded regions.
    */
    
    private final GreyscaleImage gsImg;
    private final GreyscaleImage ptImg;
    private List<List<Region>> gsRegions = null;
    private List<List<Region>> ptRegions = null;
    private List<Set<PairInt>> edgeList = null;
    
    private boolean debug = false;
    
    private long ts = 0;
    
    public MSEREdges(ImageExt img) {
        
        ImageProcessor imageProcesor = new ImageProcessor();
        
        this.gsImg = img.copyToGreyscale2();
        this.ptImg = imageProcesor.createCIELUVTheta(img, 255);
        //this.ptImg = imageProcesor.createCIELABTheta(img, 255);
        
        //enhanceContrast();
    }
    
    public MSEREdges(GreyscaleImage greyscaleImage, GreyscaleImage 
        polarThetaCIELUV) {
        
        this.gsImg = greyscaleImage;
        this.ptImg = polarThetaCIELUV;  
        
        //enhanceContrast();
    }
    
    public void setToDebug() {
        debug = true;
        ts = MiscDebug.getCurrentTimeFormatted();
    }
    
    public void enhanceContrast() {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        imageProcessor.enhanceContrast(gsImg, 4);
        
        //TODO: for the polar theta, need to use deltaE2000 to enhance the
        // skyline contrast for seattle test image, for example.  
        // can see that in the result from
        // using an adapted sobel on cielab with deltaE2000,
        imageProcessor.enhanceContrast(ptImg, 4);        
    }
    
    public void extractEdges() {
        
        if (debug) {
            MiscDebug.writeImage(ptImg, "_" + ts + "_polartheta");
        }
        
        extractRegions();
        
        //NOTE: may prefer a component tree approach to embedded regions instead
        //   of this filter
        //filterBySpatialProximity();
                
        extractBoundaries();
        
        if (debug) {
            printEdges();
        }
    }
    
    private void extractRegions() {
        
        gsRegions = extractMSERRegions(gsImg, Threshold.LESS_SENSITIVE);
    
        ptRegions = extractMSERRegions(ptImg, Threshold.DEFAULT);
    
        // for 0.001, a gsRegions limit for var of 0. is used
        //   but for 0.0001, limit should be near 0.001
        
        for (int type = 0; type < 2; ++type) {
            List<Region> list = gsRegions.get(type);
            for (int i = (list.size() - 1); i > -1; --i) {
                Region r = list.get(i);
                if ((type == 1) && r.getVariation() > 0.99) {
                    list.remove(i);
                } else if ((type == 0) 
                    //&& r.getVariation() == 0.0) {
                    && r.getVariation() < 0.99) {
                    list.remove(i);
                }
            }
        }
        
        for (int type = 0; type < 2; ++type) {
            List<Region> list = ptRegions.get(type);
            for (int i = (list.size() - 1); i > -1; --i) {
                Region r = list.get(i);
                if ((type == 1) && r.getVariation() > 0.001) {
                    list.remove(i);
                } else if ((type == 0) && r.getVariation() == 0.0) {
                    list.remove(i);
                }
            }
        }
     
        if (debug) {
            
            int[] xyCen = new int[2];
            Image imCp;
            for (int type = 0; type < 2; ++type) {
                imCp = gsImg.copyToColorGreyscale();
                int n = gsRegions.get(type).size();
                for (int i = 0; i < n; ++i) {
                    Region r = gsRegions.get(type).get(i);
                    int[] clr = ImageIOHelper.getNextRGB(i);
                    r.drawEllipse(imCp, 0, clr[0], clr[1], clr[2]);
                    r.calculateXYCentroid(xyCen, imCp.getWidth(), imCp.getHeight());
                    ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], imCp,
                        1, 255, 0, 0);
                }
                MiscDebug.writeImage(imCp, "_" + ts + "_regions_gs_"+ type);
            }
            
            for (int type = 0; type < 2; ++type) {
                imCp = ptImg.copyToColorGreyscale();
                int n = ptRegions.get(type).size();
                for (int i = 0; i < n; ++i) {
                    Region r = ptRegions.get(type).get(i);
                    int[] clr = ImageIOHelper.getNextRGB(i);
                    r.drawEllipse(imCp, 0, clr[0], clr[1], clr[2]);
                    r.calculateXYCentroid(xyCen, imCp.getWidth(), imCp.getHeight());
                    ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], imCp,
                        1, 255, 0, 0);
                    //System.out.println(type + " xy=" + xyCen[0] + "," + xyCen[1] 
                    //    + " variation=" + r.getVariation());
                }
                MiscDebug.writeImage(imCp, "_" + ts + "_regions_pt_"+ type);
            }
        }
    }

    private List<List<Region>> extractMSERRegions(GreyscaleImage img,
        Threshold threshold) {
       
        int delta = 2;
        double minArea = 0.001;//.0001
        double maxArea = 0.99;//0.1;
        double maxVariation = 0.9;//0.5;
        double minDiversity = 0.75;//0.1;//0.5
        if (threshold.equals(Threshold.LESS_SENSITIVE)) {
            maxVariation = 0.95;
            minArea = 0.01;
        }
        
        MSER mser = new MSER();
        
        int[] a = mser.readIntoArray(img);
        
        List<List<Region>> regions = mser.findRegions(a, img.getWidth(),
            img.getHeight(), delta, minArea, maxArea, maxVariation, 
            minDiversity);
        
        return regions;
    }

    private void filterBySpatialProximity() {
        
        float critDens = 2.f/10.f;
        Canonicalizer.filterBySpatialProximity(critDens, gsRegions.get(0), 
            gsImg.getWidth(), gsImg.getHeight());
        
        Canonicalizer.filterBySpatialProximity(critDens, gsRegions.get(1), 
            gsImg.getWidth(), gsImg.getHeight());
        
        Canonicalizer.filterBySpatialProximity(critDens, ptRegions.get(0), 
            ptImg.getWidth(), ptImg.getHeight());
        
        Canonicalizer.filterBySpatialProximity(critDens, ptRegions.get(1), 
            ptImg.getWidth(), ptImg.getHeight());
    }

    private void extractBoundaries() {
        
        // without a component tree for now
        edgeList = new ArrayList<Set<PairInt>>();
        
        PerimeterFinder2 finder = new PerimeterFinder2();
    
        for (int type = 0; type < 2; ++type) {
            int n = gsRegions.get(type).size();
            for (int i = 0; i < n; ++i) {
                Region r = gsRegions.get(type).get(i);
                //NOTE: accumulated points are a connected group of points
                Set<PairInt> points = new HashSet<PairInt>();
                for (int j = 0; j < r.accX.size(); ++j) {
                    int x = r.accX.get(j);
                    int y = r.accY.get(j);
                    PairInt p2 = new PairInt(x, y);
                    points.add(p2);
                    // NOTE: here is where could test membership of point in
                    //   a disjoint forrest for the region connections
                }
                Set<PairInt> embedded = new HashSet<PairInt>();
                Set<PairInt> outerBorder = new HashSet<PairInt>();
                finder.extractBorder2(points, embedded, embedded);
                
                edgeList.add(outerBorder);
                
                DFSConnectedGroupsFinder dfsFinder = new DFSConnectedGroupsFinder();
                dfsFinder.setMinimumNumberInCluster(24);
                dfsFinder.findConnectedPointGroups(embedded);
                for (int j = 0; j < dfsFinder.getNumberOfGroups(); ++j) {
                    
                    Set<PairInt> eSet = dfsFinder.getXY(j);
                    
                    Set<PairInt> embedded2 = new HashSet<PairInt>();
                    Set<PairInt> outerBorder2 = new HashSet<PairInt>();
                    finder.extractBorder2(eSet, embedded2, outerBorder2);
                
                    if (outerBorder2.size() > 16) {
                        edgeList.add(outerBorder2);
                    }
                }
            }
        }
        
        for (int type = 0; type < 2; ++type) {
            int n = ptRegions.get(type).size();
            for (int i = 0; i < n; ++i) {
                Region r = ptRegions.get(type).get(i);
                //NOTE: accumulated points are a connected group of points
                Set<PairInt> points = new HashSet<PairInt>();
                for (int j = 0; j < r.accX.size(); ++j) {
                    int x = r.accX.get(j);
                    int y = r.accY.get(j);
                    PairInt p2 = new PairInt(x, y);
                    points.add(p2);
                    // NOTE: here is where could test membership of point in
                    //   a disjoint forrest for the region connections
                }
                Set<PairInt> embedded = new HashSet<PairInt>();
                Set<PairInt> outerBorder = new HashSet<PairInt>();
                finder.extractBorder2(points, embedded, outerBorder);
                //TODO: consider extracting the outer bounds of the embedded
                //  also, but separately
                edgeList.add(outerBorder);
            }
        }
    }

    private void printEdges() {
        
        Image im = gsImg.copyToColorGreyscale();
        
        int[] clr = new int[]{255, 0, 0};
        
        for (int i = 0; i < edgeList.size(); ++i) {
            //int[] clr = ImageIOHelper.getNextRGB(i);
            ImageIOHelper.addCurveToImage(edgeList.get(i), im, 0, clr[0], 
                clr[1], clr[2]);
        }
        
        MiscDebug.writeImage(im, "_" + ts + "_edges_");
    }
    
}
