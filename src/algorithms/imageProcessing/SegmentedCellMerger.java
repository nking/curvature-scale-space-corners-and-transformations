package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.PerimeterFinder;
import algorithms.disjointSets.DisjointSet2Helper;
import algorithms.disjointSets.DisjointSet2Node;
import algorithms.imageProcessing.ImageSegmentation.BoundingRegions;
import algorithms.imageProcessing.features.BlobMedialAxes;
import algorithms.imageProcessing.features.BlobsAndPerimeters;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class SegmentedCellMerger {
 
    private final int boundaryValue;
    private final boolean hasBoundaryValue;
    private final ImageExt img;
    private final GreyscaleImage segImg;
    private final String debugTag;
    private final Logger log = Logger.getLogger(this.getClass().getName());
    
    public SegmentedCellMerger(ImageExt img, GreyscaleImage segImg, 
        int boundaryValue, String debugTag) {
        
        this.boundaryValue = boundaryValue;
        if (boundaryValue > -1) {
            hasBoundaryValue = true;
        } else {
            hasBoundaryValue = false;
        }        
        this.img = img;
        this.segImg = segImg;
        this.debugTag = debugTag;
    }
    
    public void merge() {
                
        BoundingRegions br = extractPerimetersAndBounds();
        
        Map<PairInt, DisjointSet2Node<PairInt>> cellMap = new HashMap<
            PairInt, DisjointSet2Node<PairInt>>();
        
        Map<PairInt, Integer> cellIndexMap = new HashMap<PairInt, Integer>();
            
        createMergeMap(br, cellMap, cellIndexMap);
        
        List<PairIntArray> boundaries = new ArrayList<PairIntArray>();
        
        // this order results in visiting the smallest cells first
        Stack<PairInt> stack = new Stack<PairInt>();
        for (int i = 0; i < br.getPerimeterList().size(); ++i) {
            
            PairInt xyCen = br.getBlobMedialAxes().getOriginalBlobXYCentroid(i);
            
            if (cellIndexMap.containsKey(xyCen)) {
                stack.add(xyCen);
            } 
            
            boundaries.add(br.getPerimeterList().get(i).copy());
        }
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        // boundaries gets merged too and altered to be union - intersection.
        
        Set<PairInt> visited = new HashSet<PairInt>();
        Map<PairInt, Set<PairInt>> visitedMap = new HashMap<PairInt, Set<PairInt>>();
        
        DisjointSet2Helper disjointSetHelper = new DisjointSet2Helper();
        
        while (!stack.isEmpty()) {
            
            PairInt p = stack.pop();
            
            if (visited.contains(p)) {
                continue;
            }
                        
            DisjointSet2Node<PairInt> pParentNode = disjointSetHelper.findSet(cellMap.get(p));
            
            PairInt pParent = pParentNode.getMember();
            
            Integer pIndex = cellIndexMap.get(pParent);
            
            float[] labP = br.getBlobMedialAxes().getLABColors(pIndex.intValue());
            
            double hueAngleP = Math.atan(labP[2]/labP[1]) * 180./Math.PI;
            if (hueAngleP < 0) {
                hueAngleP += 360.;
            }
            
            boolean didMerge = false;
            
            // perimeter points
            PairIntArray p1Bounds = boundaries.get(pIndex.intValue());
            for (int i = 0; i < p1Bounds.getN(); ++i) {
                
                int x = p1Bounds.getX(i);
                int y = p1Bounds.getY(i);
                
                for (int k = 0; k < dxs.length; ++k) {
                    
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    
                    if (x2 < 0 || (x2 > (w - 1)) || (y2 < 0) || (y2 > (h - 1))) {
                        continue;
                    }
                    
                    PairInt p2 = new PairInt(x2, y2);
                    
                    //find p2 in the original bounding regions data
                    Integer p2Index = br.getPointIndexMap().get(p2);
                    
                    if (p2Index == null) {
                        // might be in the large watersheds that do not get merged
                        continue;
                    }
                    
                    PairInt xyCen2 = br.getBlobMedialAxes()
                        .getOriginalBlobXYCentroid(p2Index.intValue());
                    
                    // find the "representative" parent not of the cell
                    DisjointSet2Node<PairInt> p2ParentNode = disjointSetHelper.findSet(cellMap.get(xyCen2));
                    
                    if (p2ParentNode.equals(pParentNode)) {
                        continue;
                    }
                    
                    PairInt p2Parent = p2ParentNode.getMember();
                    
                    if (hasBeenVisited(visitedMap, pParent, p2Parent)) {
                        continue;
                    }
                                           
                    p2Index = cellIndexMap.get(p2Parent);
                    
                    float[] labP2 = br.getBlobMedialAxes().getLABColors(p2Index.intValue());
                
                    double hueAngleP2 = Math.atan(labP2[2]/labP2[1]) * 180./Math.PI;
                    if (hueAngleP2 < 0) {
                        hueAngleP2 += 360.;
                    }                            
                
                    double deltaE = cieC.calcDeltaECIE94(labP[0], labP[1], labP[2], 
                        labP2[0], labP2[1], labP2[2]);

                    float deltaHA = AngleUtil.getAngleDifference((float)hueAngleP, 
                        (float)hueAngleP2);

                    log.info(String.format("%s %s  deltaE=%.2f  dHA=%d",
                        pParent.toString(), p2Parent.toString(), (float)deltaE,
                        Math.round(deltaHA)));
                    
                    //-- if similar 
                    //   -- alter the boundaries for index of pParent: union - intersection of points
                    //      (intersection is found using point in polygon tests against other)
                    //   -- merge the cells:
                    //      DisjointSet2Node<PairInt> parent = 
                    //          disjointSetHelper.union(parentMap.get(reprPoint),
                    //          parentMap.get(qPoint));
                    //   -- reset the boundaries to the edited for the new parent index
                    //   -- set didMerge to tru
                    
                    addToVisited(visitedMap, pParent, p2Parent);
                }
            }
           
            if (didMerge) {
                visited.add(p);
            }
        }
        
        /*
        will use a DFS style visitor pattern w/ stack ordered to process smallest
        cells first.
        visit each cell and merge any similar neighbors.
        
        bounding regions have been ordered by descending size already, so can
        populate the stack from 0 to n-1.
        
        for deltaE, range of similar is 2.3 through 5.5 or ?  (max diff is 28.8).
        might need to use adjacent points instead of set averages.
        
        The gingerbread man and the background building have deltaE=5.7 and deltaL=3.64
        Two cells of the gingerbread, similar in color, but different shade
        have deltaE = 5.13 and deltaL=1.6
        
        the hue angles are very strong indicators for the two examples just given.
        o2 is also and looks useful for the icecream and cupcake.
        but o2 for the gingerbread man's feet would merge with the grass while it's deltaEis 9.17
       
        so need to look at the color spaces in detail to use more than deltaE for similarity...
        */
       
        /*
        consider a step to merge completely embedded with the encapsulating cell
        */
        
    }
    
    private BoundingRegions extractPerimetersAndBounds() {
        
        List<Set<PairInt>> boundaryValueSets = new ArrayList<Set<PairInt>>();
        
        int smallestGroupLimit = 1;
        int largestGroupLimit = Integer.MAX_VALUE;
        boolean filterOutImageBoundaryBlobs = false;
        boolean filterOutZeroPixels = false;
        
        //TODO: this may need revision.  wanting to exclude processing for
        // contiguous regions which are a large fraction of image.
        // these are usually background.
        
        
        // runtime complexity is N_freq * O(N) where N_freq is at most 256 and
        // the O(N) term may be as high as O(N*8) if highly connected.
        List<Set<PairInt>> blobs =  BlobsAndPerimeters.extractBlobsFromSegmentedImage(
            segImg, smallestGroupLimit, largestGroupLimit,
            filterOutImageBoundaryBlobs, filterOutZeroPixels, debugTag);
        
        // find the sets which are the boundary values
        if (hasBoundaryValue) {
            List<Integer> mv = new ArrayList<Integer>();
            for (int i = 0; i < blobs.size(); ++i) {
                Set<PairInt> blob = blobs.get(i);
                PairInt p = blob.iterator().next();
                if (segImg.getValue(p) == boundaryValue) {
                    mv.add(Integer.valueOf(i));
                }
            }
            for (int i = (mv.size() - 1); i > -1; --i) {
                Integer mvIndex = mv.get(i);
                Set<PairInt> blob = blobs.remove(mvIndex.intValue());
                boundaryValueSets.add(blob);
            }
        }
                
        // --- sort by descending sizes the remaining blobs ---- 
        int[] sizes = new int[blobs.size()];
        int[] indexes = new int[sizes.length];
        for (int i = 0; i < blobs.size(); ++i) {
            sizes[i] = blobs.get(i).size();
            indexes[i] = i;
        }
        
        // removing any blobs which are larger than 0.2 percent of image size also
        float nPixels = img.getNPixels();
        
        MultiArrayMergeSort.sortByDecr(sizes, indexes);
        List<Set<PairInt>> tmp = new ArrayList<Set<PairInt>>();
        
        for (int i = 0; i < sizes.length; ++i) {
            int idx = indexes[i];
            Set<PairInt> blob = blobs.get(idx);
            float frac = (float)blob.size()/nPixels;
            if (frac < 0.2) {
                tmp.add(blob);
            }
        }
        blobs.clear();
        blobs.addAll(tmp);
        
        //---- begin section to log colors to look at selecting matchable bounds by color ------
        CIEChromaticity cieC = new CIEChromaticity();
        //MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        List<Double> lAvg = new ArrayList<Double>();
        List<Double> aAvg = new ArrayList<Double>();
        List<Double> bAvg = new ArrayList<Double>();
        
        for (int i = 0; i < blobs.size(); ++i) {
            double redSum = 0;
            double greenSum = 0;
            double blueSum = 0;
            for (PairInt p : blobs.get(i)) {
                int x = p.getX();
                int y = p.getY();
                int red = img.getR(x, y);
                int green = img.getG(x, y);
                int blue = img.getB(x, y);
                redSum += red;
                greenSum += green;
                blueSum += blue;
            }
            double n = (double)blobs.get(i).size();
            redSum /= n;
            greenSum /= n;
            blueSum /= n;
            float[] avgLAB = cieC.rgbToCIELAB((int)Math.round(redSum), 
                (int)Math.round(greenSum), (int)Math.round(blueSum));
            
            lAvg.add(Double.valueOf(avgLAB[0]));
            aAvg.add(Double.valueOf(avgLAB[1]));
            bAvg.add(Double.valueOf(avgLAB[2]));
        
            /*PairInt xyCen = curveHelper.calculateXYCentroids(blobs.get(i));
            String str = String.format(
                "[%d] cen=(%d,%d) avgL=%.3f avgA=%.3f  avgB=%.3f  nPts=%d",
                i, xyCen.getX(), xyCen.getY(),
                avgLAB[0], avgLAB[1], avgLAB[2], blobs.get(i).size());
            log.info(str);*/
        }
        
        // place points in boundaryValueSets, individually, into the 
        // adjacent sets they best match
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        for (int i = 0; i < boundaryValueSets.size(); ++i) {
            
            Set<PairInt> bvSet = boundaryValueSets.get(i);
            
            for (PairInt p : bvSet) {
                
                float[] labP = img.getCIELAB(p.getX(), p.getY());
                
                double minDeltaE = Double.MAX_VALUE;
                int minDeltaEIdx = -1;
                
                // compare each point to colors in adjacent sets and choose closest.
                // TODO: can improve this with another datastructure
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = p.getX() + dxs[k];
                    int y2 = p.getY() + dys[k];
                    PairInt p2 = new PairInt(x2, y2);
                    
                    for (int j = 0; j < blobs.size(); ++j) {
                        Set<PairInt> blob = blobs.get(j);
                        if (blob.contains(p2)) {
                            Double ell = lAvg.get(j);
                            Double ay = aAvg.get(j);
                            Double bee = aAvg.get(j);
                            double deltaE = cieC.calcDeltaECIE94(labP[0], 
                                labP[1], labP[2], 
                                ell.floatValue(), ay.floatValue(), bee.floatValue());
                            if (deltaE < 0) {
                                deltaE *= -1;
                            }
                            if (deltaE < minDeltaE) {
                                minDeltaE = deltaE;
                                minDeltaEIdx = j;
                            }
                            break;
                        }
                    }
                }
                
                // some of the boundary value pixels are connected to large regions
                // so are not always adjacent to a non-boundary value region.
                if (minDeltaEIdx > -1) {
                    blobs.get(minDeltaEIdx).add(p);
                }
            }
        }
        
        // create a map for reverse look-ups later
        Map<PairInt, Integer> pointIndexMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < blobs.size(); ++i) {
            Integer key = Integer.valueOf(i);
            Set<PairInt> blob = blobs.get(i);
            for (PairInt p : blob) {
                pointIndexMap.put(p, key);
            }
        }
        
        // less than O(N)
        List<Set<PairInt>> borderPixelSets = BlobsAndPerimeters.extractBlobPerimeterAsPoints(
            blobs, segImg.getWidth(), segImg.getHeight());
        
        assert(blobs.size() == borderPixelSets.size());
        
        List<PairIntArray> perimetersList = new ArrayList<PairIntArray>();
        
        float srchRadius = (float)Math.sqrt(2);
        
        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        BlobMedialAxes bma = new BlobMedialAxes(blobs, lAvg, aAvg, bAvg);
        
        for (int i = 0; i < borderPixelSets.size(); ++i) {
                                    
            Set<PairInt> blob = blobs.get(i);
            Set<PairInt> borderPixels = borderPixelSets.get(i);
            
            // approx O(N_perimeter), but has factors during searches that could be improved
            PairIntArray orderedPerimeter = perimeterFinder.orderThePerimeter(
                borderPixels, blob, srchRadius, bma, i);
  
            /*Image imgCp = img.copyImage();
            ImageIOHelper.addCurveToImage(orderedPerimeter, imgCp, 2, 255, 0, 0);                  
            MiscDebug.writeImage(imgCp, "_" + i + "_" + MiscDebug.getCurrentTimeFormatted()); 
            */
            
            // runtime complexity is O(N_perimeter_pts * lg_2(N_perimeter_pts)
            // remove straight line segments except their endpoints to make simpler
            // polynomial for "point in polygon" tests
            imageSegmentation.makeStraightLinesHollow(orderedPerimeter, img.getWidth(), 
                img.getHeight(), srchRadius);
        
            /*
            Image imgCp = img.copyImage();
            ImageIOHelper.addCurveToImage(orderedPerimeter, imgCp, 2, 255, 0, 0);                  
            MiscDebug.writeImage(imgCp, "_" + i + "_" + MiscDebug.getCurrentTimeFormatted()); 
            */
            
            perimetersList.add(orderedPerimeter);
        }
        
        BoundingRegions br = new BoundingRegions(perimetersList, bma, 
            pointIndexMap);
        
        return br;
    }

    private void createMergeMap(BoundingRegions br, 
        Map<PairInt, DisjointSet2Node<PairInt>> outputParentMap,
        Map<PairInt, Integer> outputParentIndexMap) {
        
        DisjointSet2Helper disjointSetHelper = new DisjointSet2Helper();
        
        // store the centroids in a disjoint set/forrest
        
        // init map
        BlobMedialAxes bma = br.getBlobMedialAxes();
        int n = bma.getNumberOfItems();
        
        for (int i = 0; i < n; ++i) {
            
            PairInt xyCen = bma.getOriginalBlobXYCentroid(i);
            
            PairInt p = xyCen.copy();
            
            DisjointSet2Node<PairInt> pNode =
                disjointSetHelper.makeSet(new DisjointSet2Node<PairInt>(p));
            
            outputParentMap.put(p, pNode);
            
            outputParentIndexMap.put(p, Integer.valueOf(i));
        }        
    }

    private boolean hasBeenVisited(Map<PairInt, Set<PairInt>> visitedMap, 
        PairInt pParent, PairInt p2Parent) {
        
        Set<PairInt> set = visitedMap.get(pParent);
        if (set == null) {
            return false;
        }
        
        return set.contains(p2Parent);
    }

    private void addToVisited(Map<PairInt, Set<PairInt>> visitedMap, 
        PairInt pParent, PairInt p2Parent) {
        
        Set<PairInt> set = visitedMap.get(pParent);
        if (set == null) {
            set = new HashSet<PairInt>();
            visitedMap.put(pParent, set);
        }
        set.add(p2Parent);
        
        set = visitedMap.get(p2Parent);
        if (set == null) {
            set = new HashSet<PairInt>();
            visitedMap.put(p2Parent, set);
        }
        set.add(pParent);
    }

}
