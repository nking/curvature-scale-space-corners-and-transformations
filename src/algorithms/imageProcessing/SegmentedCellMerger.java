package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.ClosestPairBetweenSets;
import algorithms.compGeometry.ClosestPairBetweenSets.ClosestPairInt;
import algorithms.compGeometry.PerimeterFinder;
import algorithms.compGeometry.PointInPolygon;
import algorithms.disjointSets.DisjointSet2Helper;
import algorithms.disjointSets.DisjointSet2Node;
import algorithms.imageProcessing.ImageSegmentation.BoundingRegions;
import algorithms.imageProcessing.features.BlobMedialAxes;
import algorithms.imageProcessing.features.BlobsAndPerimeters;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
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
        
        /*        
        current steps (and changes to be made):
        
        -- the centroid of each blob, that is bounded region, from the segmented
           image is used to create a DisjointSet2Node.
           each of these is placed in a map w/ key being the centroid and
           value being the node.
        -- a map with key being centroid and value being cell bounds is made
        -- a map with key being centroid and value being cell points is made
        
        a stack is populated with the centroids, so that they are popped in
        order of smallest cells first.
        
        as each item is popped from the stack, the parent cell for that centroid
        is found.
        
           Then the adjacent cells are visited to compare colors and merge if similar.
           -- the adjacency map is determined by the boundaries and the points
               immediately next to the boundaries.
        
           -- once the adjacent cells are determined, each is visited.
              -- colors are compared and the similar cells are merged.
                 -- merging currently is done as union of the disjoint set forest.
                    then updating the maps for the merge and new parent.
        
           connectivity stored as adjacency sets and merged sets:
              adjacencyMap:
                  map of key=centroid, value=set of centroids of adjacent regions
              mergedMap:
                  map of key=centroid, value=set of centroids of regions merged into key region

           the two maps are altered upon a merge.     
              for merge of key A with key B,
                  mergedMap for keyA gets keyB and its values added to it and keyB gets removed.
                  adjacencyMap for keyA gains the entries for keyB, but then the contents of 
                      mergedMap for key A are deleted from the values to make sure that the
                      adjacent cell set doesn't include internal members.
        */
                
        BoundingRegions br = extractPerimetersAndBounds();
        
        final Map<PairInt, DisjointSet2Node<PairInt>> cellMap = 
            new HashMap<PairInt, DisjointSet2Node<PairInt>>();
        
        final Map<PairInt, Integer> cellIndexMap = new HashMap<PairInt, Integer>();
            
        populateCellMaps(br, cellMap, cellIndexMap);
        
        // key = cell centroid, value = set of adjacent cell centroids
        Map<PairInt, Set<PairInt>> adjacencyMap = createAdjacencyMap(br,
            img.getWidth(), img.getHeight());
        
        // key = cell centroid, value = set of cell centroids merged with this one
        Map<PairInt, Set<PairInt>> mergedMap = createMergeMap(br);
        
        
        // this order results in visiting the smallest cells first
        Stack<PairInt> stack = new Stack<PairInt>();
        for (int i = 0; i < br.getPerimeterList().size(); ++i) {
            
            PairInt xyCen = br.getBlobMedialAxes().getOriginalBlobXYCentroid(i);
            
            if (cellIndexMap.containsKey(xyCen)) {
                stack.add(xyCen);
            }            
        }
        
        if (debugTag != null && !debugTag.equals("")) {
            List<PairIntArray> boundaries = br.getPerimeterList();
            long ts = MiscDebug.getCurrentTimeFormatted();
            ImageExt imgCp = img.copyToImageExt();
            ImageIOHelper.addAlternatingColorCurvesToImage(
                boundaries.toArray(new PairIntArray[boundaries.size()]),
                imgCp, 0);
            MiscDebug.writeImage(imgCp, debugTag + "_boundaries_" + ts);
        }
                
        CIEChromaticity cieC = new CIEChromaticity();
        
        Set<PairInt> visited = new HashSet<PairInt>();
        Map<PairInt, Set<PairInt>> visitedMap = new HashMap<PairInt, Set<PairInt>>();
        
        DisjointSet2Helper disjointSetHelper = new DisjointSet2Helper();
        
        while (!stack.isEmpty()) {
            
            PairInt p = stack.pop();
            
            if (visited.contains(p)) {
                continue;
            }
            
            DisjointSet2Node<PairInt> pNode = cellMap.get(p);
                        
            DisjointSet2Node<PairInt> pParentNode = disjointSetHelper.findSet(pNode);
            
            PairInt pParent = pParentNode.getMember();
            
            Integer pIndex = cellIndexMap.get(pParent);
            
            float[] labP = br.getBlobMedialAxes().getLABColors(pIndex.intValue());
            
            double hueAngleP = Math.atan(labP[2]/labP[1]) * 180./Math.PI;
            if (hueAngleP < 0) {
                hueAngleP += 360.;
            }
            
            boolean didMerge = false;
            
            List<PairInt> neighborKeys = new ArrayList<PairInt>(
                adjacencyMap.get(pParent));
            
            for (PairInt p2 : neighborKeys) {
                
                DisjointSet2Node<PairInt> p2Node = cellMap.get(p2);
                
                // find the "representative" parent node of the cell
                DisjointSet2Node<PairInt> p2ParentNode = disjointSetHelper.findSet(p2Node);
                    
                if (p2ParentNode.equals(pParentNode)) {
                    continue;
                }
                    
                PairInt p2Parent = p2ParentNode.getMember();
                
                if (hasBeenVisited(visitedMap, pParent, p2Parent)) {
                    continue;
                }
                                           
                Integer p2ParentIndex = cellIndexMap.get(p2Parent);
                    
                float[] labP2 = br.getBlobMedialAxes().getLABColors(p2ParentIndex.intValue());

                double hueAngleP2 = Math.atan(labP2[2]/labP2[1]) * 180./Math.PI;
                if (hueAngleP2 < 0) {
                    hueAngleP2 += 360.;
                }                            

                double deltaE = cieC.calcDeltaECIE94(labP[0], labP[1], labP[2], 
                    labP2[0], labP2[1], labP2[2]);

                float deltaHA = AngleUtil.getAngleDifference((float)hueAngleP, 
                    (float)hueAngleP2);

                addToVisited(visitedMap, pParent, p2Parent);
                
                log.info(String.format("%s %s  deltaE=%.2f  dHA=%d",
                    pParent.toString(), p2Parent.toString(), (float)deltaE,
                    Math.round(deltaHA)));
                
                /*
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

                if (Math.abs(deltaE) > 5 && Math.abs(deltaHA) > 5) {
                    continue;
                }
                
                stack.add(p2Parent);
                
                DisjointSet2Node<PairInt> parentOfMergeNode = 
                    disjointSetHelper.union(pParentNode, p2ParentNode);

                PairInt parentOfMerge = parentOfMergeNode.getMember();

                if (parentOfMerge.equals(pParent)) {
                    
                    //mergedMap for key pParent gets values of p2Parent 
                    //and the p2Parent key gets removed
                    
                    Set<PairInt> p2ParentMergedSet = mergedMap.get(p2Parent);
                    
                    Set<PairInt> pParentMergedSet = mergedMap.get(pParent);
                    pParentMergedSet.addAll(p2ParentMergedSet);
                    
                    mergedMap.remove(p2Parent);
                    
                    // adjacencyMap for key pParent gets values of p2Parent, 
                    // and then subtracts all internal members from the final set
                    Set<PairInt> p2ParentAdjacencySet = adjacencyMap.get(p2Parent);
                    
                    Set<PairInt> pParentAdjacencySet = adjacencyMap.get(pParent);
                    pParentAdjacencySet.addAll(p2ParentAdjacencySet);
                    pParentAdjacencySet.removeAll(pParentMergedSet);
                    
                } else {
                    
                    if (!p.equals(pParent)) {
                        assert(!adjacencyMap.containsKey(p));
                        assert(!mergedMap.containsKey(p));
                    }
                    
                    //mergedMap for key p2Parent gets values of pParent 
                    //and the pParent key gets removed
                    
                    Set<PairInt> pParentMergedSet = mergedMap.get(pParent);
                    
                    Set<PairInt> p2ParentMergedSet = mergedMap.get(p2Parent);
                    p2ParentMergedSet.addAll(pParentMergedSet);
                    
                    mergedMap.remove(pParent);
                    
                    // adjacencyMap for key p2Parent gets values of pParent, 
                    // and then subtracts all internal members from the final set
                    Set<PairInt> pParentAdjacencySet = adjacencyMap.get(pParent);
                    
                    Set<PairInt> p2ParentAdjacencySet = adjacencyMap.get(p2Parent);
                    p2ParentAdjacencySet.addAll(pParentAdjacencySet);
                    p2ParentAdjacencySet.removeAll(p2ParentMergedSet);
                }
             
                didMerge = true;                  
            }
           
            if (didMerge) {
                visited.add(p);
            }
        }
               
        /*
        here, have final results in 
            // key = cell centroid, value = set of adjacent cell centroids
            Map<PairInt, Set<PairInt>> adjacencyMap
        
            // key = cell centroid, value = set of cell centroids merged with this one
            Map<PairInt, Set<PairInt>> mergedMap
        
        the cell indexes can then be gathered and then new blobs made from
            the points in those sets.
        
        line 500 to end of its method are what is then needed here to rebuild the merged as
        bounding regions. 
        
TODO:  consider encapsulating the maps in an extended disjoint helper
        to make the union and findSet updates also handle updating the
        maps.
        
        */
        
        int z = 1;
    }
    
    private BoundingRegions extractPerimetersAndBounds() {
        
        List<Set<PairInt>> boundaryValueSets = new ArrayList<Set<PairInt>>();
        
        int smallestGroupLimit = 1;
        int largestGroupLimit = Integer.MAX_VALUE;
        boolean filterOutImageBoundaryBlobs = false;
        boolean filterOutZeroPixels = false;
        boolean use8Neighbors = true;
        
        //TODO: this may need revision.  wanting to exclude processing for
        // contiguous regions which are a large fraction of image.
        // these are usually background.
        
        
        // runtime complexity is N_freq * O(N) where N_freq is at most 256 and
        // the O(N) term may be as high as O(N*8) if highly connected.
        List<Set<PairInt>> blobs =  BlobsAndPerimeters.extractBlobsFromSegmentedImage(
            segImg, smallestGroupLimit, largestGroupLimit,
            filterOutImageBoundaryBlobs, filterOutZeroPixels, 
            use8Neighbors, debugTag);
        
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
        
        BoundingRegions br = new BoundingRegions(perimetersList, bma, pointIndexMap);
        
        return br;
    }

    private void populateCellMaps(BoundingRegions br, 
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
            
            assert(pNode.getMember().equals(p));
            assert(pNode.getParent().equals(pNode));
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

    private Map<PairInt, Set<PairInt>> createAdjacencyMap(BoundingRegions br,
        int imageWidth, int imageHeight) {

        Map<PairInt, Set<PairInt>> adjacencyMap = new HashMap<PairInt, Set<PairInt>>();
        
        List<PairIntArray> boundaries = br.getPerimeterList();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        for (int i = 0; i < boundaries.size(); ++i) {
            
            PairInt key = br.getBlobMedialAxes().getOriginalBlobXYCentroid(i);
            
            Set<PairInt> adjacencySet = new HashSet<PairInt>();
            
            adjacencyMap.put(key, adjacencySet);
            
            PairIntArray boundary = boundaries.get(i);
           
            for (int j = 0; j < boundary.getN(); ++j) {
                
                int x = boundary.getX(j);
                int y = boundary.getY(j);
        
                for (int k = 0; k < dxs.length; ++k) {
                    
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    
                    if (x2 < 0 || (x2 > (imageWidth - 1)) || (y2 < 0) || 
                        (y2 > (imageHeight - 1))) {
                        continue;
                    }
                    
                    PairInt p2 = new PairInt(x2, y2);
                    
                    //find p2 in the original bounding regions data
                    Integer p2Index = br.getPointIndexMap().get(p2);
                    
                    if (p2Index == null || (p2Index.intValue() == i)) {
                        continue;
                    }
                    
                    PairInt xyCen2 = br.getBlobMedialAxes()
                        .getOriginalBlobXYCentroid(p2Index.intValue());
                    
                    adjacencySet.add(xyCen2);
                }
            }
            
            if (adjacencySet.isEmpty()) {
                // look at these one by one as debugging
                int z = 1;
            }
        }
        
        return adjacencyMap;
    }

    private Map<PairInt, Set<PairInt>> createMergeMap(BoundingRegions br) {
        
        Map<PairInt, Set<PairInt>> mergeMap = new HashMap<PairInt, Set<PairInt>>();
                        
        for (int i = 0; i < br.getPerimeterList().size(); ++i) {
            
            PairInt key = br.getBlobMedialAxes().getOriginalBlobXYCentroid(i);
           
            Set<PairInt> mergeSet = new HashSet<PairInt>();
            mergeSet.add(key);
            
            mergeMap.put(key, mergeSet);
        }
        
        return mergeMap;
    }

}
