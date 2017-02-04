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
    extractEdges:
        -- uses PerimeterFinder2 extract the outer boundaries from each mser
           region points, and the same for the contiguous embedded regions.
    mergeAndExtractEdges
        -- merges adjacent regions by similarity of color (and maybe texture in future)
        -- extracts edges using PerimeterFinder2 just as above
        
    */
    
    private enum STATE {
        INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
    }
    
    private final GreyscaleImage gsImg;
    private final GreyscaleImage ptImg;
    private List<Region> regions = null;
    // NOTE: the edges indexes do not correspond to the regions indexes
    private List<Set<PairInt>> edgeList = null;
    
    private boolean debug = false;
    
    private boolean useLowerContrastLimits = false;
    
    private long ts = 0;
    
    private STATE state = null;
    
    public MSEREdges(ImageExt img) {
        
        ImageProcessor imageProcesor = new ImageProcessor();
        
        this.gsImg = img.copyToGreyscale2();
        this.ptImg = imageProcesor.createCIELUVTheta(img, 255);
        
        state = STATE.INITIALIZED;
    }
    
    public MSEREdges(GreyscaleImage greyscaleImage, GreyscaleImage 
        polarThetaCIELUV) {
        
        this.gsImg = greyscaleImage;
        this.ptImg = polarThetaCIELUV;  
        
        state = STATE.INITIALIZED;
    }
    
    public void setToDebug() {
        debug = true;
        ts = MiscDebug.getCurrentTimeFormatted();
    }
    
    /**
     * use this to change the limits to keep lower contrast edges.
     * NOTE that this increases the number of edges that are noise.
     */
    public void setToLowerContrast() {
        useLowerContrastLimits = true;
    }
    
    public void extractEdges() {
        
        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (!state.equals(STATE.INITIALIZED)) {
            throw new IllegalStateException("can only perform extraction of "
                + "edges once");
        }
        
        if (debug) {
            MiscDebug.writeImage(ptImg, "_" + ts + "_polartheta");
        }
        
        extractRegions();
           
        extractBoundaries();
        
        if (debug) {
            printEdges();
        }
    }
    
    /**
     * merge regions and then extract edges.  this method is most useful when
     * setToLowerContrast() has been invoked before. 
     */
    public void mergeAndExtractEdges() {
        
        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (!state.equals(STATE.INITIALIZED)) {
            throw new IllegalStateException("can only perform extraction of "
                + "edges once");
        }
        
        if (debug) {
            MiscDebug.writeImage(ptImg, "_" + ts + "_polartheta");
        }
        
        extractRegions();
        
        mergeRegions();
        
        extractBoundaries();
        
        if (debug) {
            printEdges();
        }
    }
    
    private void extractRegions() {
        
        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("can only perform extraction of "
                + "edges once");
        }
        
        List<List<Region>> gsRegions = extractMSERRegions(gsImg, Threshold.LESS_SENSITIVE);
    
        List<List<Region>> ptRegions = extractMSERRegions(ptImg, Threshold.DEFAULT);
    
        regions = new ArrayList<Region>();
        
        // for minArea 0.001, a gsRegions limit for var of 0. is used
        //   but for 0.0001, limit should be near 0.001
        
        for (int type = 0; type < 2; ++type) {
            List<Region> list = gsRegions.get(type);
            for (int i = (list.size() - 1); i > -1; --i) {
                Region r = list.get(i);
                if ((type == 1) && r.getVariation() > 2.) {
                    list.remove(i);
                } else if ((type == 0) 
                    //&& r.getVariation() == 0.0) {
                    && r.getVariation() < 2.) {
                    list.remove(i);
                } else {
                    regions.add(r);
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
                } else {
                    regions.add(r);
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
        
        state = STATE.REGIONS_EXTRACTED;
    }

    private List<List<Region>> extractMSERRegions(GreyscaleImage img,
        Threshold threshold) {
       
        int delta = 2;
        double minArea = 0.001;//.0001
        double maxArea = 0.99;//0.1;
        double maxVariation = 0.9;//0.5;
        double minDiversity = 0.75;//0.1;//0.5
        if (threshold.equals(Threshold.LESS_SENSITIVE)) {
            maxVariation = 2.0;
            minArea = 0.0075;
            if (useLowerContrastLimits) {
                minDiversity = 0.1;//.4
            }
        }
        
        MSER mser = new MSER();
        
        int[] a = mser.readIntoArray(img);
        
        List<List<Region>> regions = mser.findRegions(a, img.getWidth(),
            img.getHeight(), delta, minArea, maxArea, maxVariation, 
            minDiversity);
        
        return regions;
    }

    private void extractBoundaries() {
        
        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("can only perform extraction of "
                + "edges once");
        }
        
        // without a component tree for now
        edgeList = new ArrayList<Set<PairInt>>();
        
        PerimeterFinder2 finder = new PerimeterFinder2();
    
        for (int rListIdx = 0; rListIdx < regions.size(); ++rListIdx) {
            Region r = regions.get(rListIdx);
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
        
        state = STATE.EDGES_EXTRACTED;
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
    
    
    private void mergeRegions() {
        /*
        -- create a point index map
        -- put the region indexes in a disjoint forest
        -- calculate the properties to use for merging for each region in lists
           that have the same indexes as regions
        -- create an adjacency map of region indexes
        -- bfs of regions:
           slighlty different than the merge pattern I usually use.
           usually I update the colors of a region as new points are merged
           and update the adjacency indexes.
     
           for this version, instead, will only write the adjacency indexes once
           and the color calcs once and rely on the disjoint forest to update
           merged region indexes.
          
           needs a set<PairInt> of already compared region indexes with the 
           smaller being first
        
           the regions get placed in a queue or list and iterated over.
           for each region
               get adjacent indexes
               for each neighbor region index
                  if already compared, continue
                  compare color (and/or texture) properties
                  if within merge limit
                      merge involves updating the disjoint forest index with a untion
                  store neighbor and current region index in compared set
         if more than one iteration is needed, then need to update the adjacency 
             indexes and colors and populate the queue with only the
             top most forest indexes.
         
        iterate over the regions list:
        use the disjoint forest findSet to find if a region's representative index
           has changed from original, and if so, put all of its points in the
           other region and clear the current region set.
           (adding to the other region involves using the accumulate method)
        
        iterate over the regions list backwrds and remove the regions which have no
            points
       */
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public List<Set<PairInt>> getEdges() {
        
        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (!state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("must use one of the extraction"
                + " methods first");
        }
        
        return edgeList;
    }
}
