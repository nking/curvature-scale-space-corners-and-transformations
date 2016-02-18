package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.ClosestPairBetweenSets;
import algorithms.compGeometry.ClosestPairBetweenSets.ClosestPairInt;
import algorithms.compGeometry.NearestPointsInLists;
import algorithms.compGeometry.PerimeterFinder;
import algorithms.compGeometry.PointInPolygon;
import algorithms.disjointSets.DisjointSet2Helper;
import algorithms.disjointSets.DisjointSet2Node;
import algorithms.imageProcessing.ImageSegmentation.BoundingRegions;
import algorithms.imageProcessing.features.BlobMedialAxes;
import algorithms.imageProcessing.features.BlobsAndPerimeters;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import algorithms.util.PairIntPair;
import algorithms.util.ResourceFinder;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author nichole
 */
public class SegmentedCellMerger {
 
    private final int boundaryValue;
    private final boolean hasBoundaryValue;
    private final ImageExt img;
    private final List<Set<PairInt>> segmentedCellList;
    private final String debugTag;
    private final Logger log = Logger.getLogger(this.getClass().getName());
    
    public SegmentedCellMerger(ImageExt img, List<Set<PairInt>> theSegmentedCellList, 
        int boundaryValue, String debugTag) {
        
        this.boundaryValue = boundaryValue;
        if (boundaryValue > -1) {
            hasBoundaryValue = true;
        } else {
            hasBoundaryValue = false;
        }        
        
        this.img = img;
        
        this.segmentedCellList = new ArrayList<Set<PairInt>>(theSegmentedCellList.size());
        for (Set<PairInt> set : theSegmentedCellList) {
            Set<PairInt> set2 = new HashSet<PairInt>(set);
            segmentedCellList.add(set2);
        }
        
        this.debugTag = debugTag;
    }
    
    private Set<PairIntPair> simClass = null;
    private Set<PairIntPair> diffClass = null;
    private String simFilePath = null;
    private FileWriter simWriter = null;
    
    /**
     * set up pairs of centroids to write data to text file for classes 
     * "similar" and "different".
     * @param similarClass
     * @param differentClass
     * @param outfileSuffix 
     */
    public void setClassPairs(Set<PairIntPair> similarClass, 
        Set<PairIntPair> differentClass, String outfileSuffix) throws IOException {
        
        simClass = similarClass;
        diffClass = differentClass;
        this.simFilePath = ResourceFinder.getAFilePathInTmpData("seg_color_"+ outfileSuffix + ".csv");
    }
    
    public void merge() {
        
        if (simFilePath != null) {
            try {
                File fl = new File(simFilePath);
                simWriter = new FileWriter(fl);
            } catch (IOException ex) {
                Logger.getLogger(SegmentedCellMerger.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        /*                
        -- the centroid of each blob, that is bounded region, from the segmented
           image is used to create a DisjointSet2Node.
           each of these is placed in a map w/ key being the centroid and
           value being the node.
        
        a stack is populated with the centroids, so that they are popped in
        order of smallest cells to largest.
        
        as each item is popped from the stack, the parent cell for that centroid
        is found.
        
           Then the adjacent cells are visited to compare colors and merge if similar.
           (the adjacency map is determined  by the boundaries and the points
           immediately next to the boundaries).
        
           -- each cell in the adjacent cells is visited.
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
        
//writeLabelFile(br, cellIndexMap);

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
        
        int nBefore = mergedMap.size();
                
        CIEChromaticity cieC = new CIEChromaticity();
        
        Set<PairInt> visited = new HashSet<PairInt>();
        Map<PairInt, Set<PairInt>> visitedMap = new HashMap<PairInt, Set<PairInt>>();
        
        DisjointSet2Helper disjointSetHelper = new DisjointSet2Helper();
        
        //TODO: improving this transformation
        double[][] ldaMatrix = getLDASegmentationMatrix();
               
        while (!stack.isEmpty()) {
            
            PairInt p = stack.pop();
            
            if (visited.contains(p)) {
                continue;
            }
            
            PairInt originalP = p.copy();
            Integer originalPIndex = cellIndexMap.get(originalP);
            
            DisjointSet2Node<PairInt> pNode = cellMap.get(p);
                        
            DisjointSet2Node<PairInt> pParentNode = disjointSetHelper.findSet(pNode);
            
            PairInt pParent = pParentNode.getMember();
            
            float[] labP = br.getBlobMedialAxes().getLABColors(originalPIndex.intValue());
            
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
                
                Integer p2Index = cellIndexMap.get(p2);
                float[] labP2 = br.getBlobMedialAxes().getLABColors(p2Index.intValue());
                
                double deltaE = cieC.calcDeltaECIE94(labP[0], labP[1], labP[2], 
                    labP2[0], labP2[1], labP2[2]);
                
                addToVisited(visitedMap, pParent, p2Parent);
                
                /*double[][] data = new double[4][1];
                data[0][0] = Math.abs(deltaE);
                data[1][0] = Math.abs(deltaO1);
                data[2][0] = Math.abs(deltaO2);
                data[3][0] = Math.abs(deltaO3);
                
                double[][] transformed = MatrixUtil.dot(ldaMatrix, data);
                double ldaX = transformed[0][0];
                double ldaY = transformed[1][0];
                */
                
                /*
                 for deltaE, range of similar is 2.3 through 5.5 or ?  (max diff is 28.8).
                 might need to use adjacent points instead of set averages.
        
                 The gingerbread man and the background building have deltaE=5.7 and deltaL=3.64
                 Two cells of the gingerbread, similar in color, but different shade
                 have deltaE = 5.13 and deltaL=1.6
                */

                //TODO: revising these vectors
                
                if (Math.abs(deltaE) > 2.3) {
                //if (ldaY > 24.5 || ldaX < -50) {
                    continue;
                }
                
                stack.add(p2Parent);
                
                DisjointSet2Node<PairInt> parentOfMergeNode = 
                    disjointSetHelper.union(pParentNode, p2ParentNode);

                PairInt parentOfMerge = parentOfMergeNode.getMember();

                if (parentOfMerge.equals(pParent)) {
                    
                    // parent is still the main comparator
                    
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
                    
                    // the new parent key has become a neighbor, so after the
                    // merges, have to re-assign the local variables within
                    // the upper block
                    
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
                    
                    // --- reassign upper block variables to use new parent information ---
                    p = p2;
                    
                    // NOTE: there's some inefficiency in that the visited check
                    // is now bypassed, but the number of original neighbors is
                    // small so the possible extra work should be idempotent.
                    
                    pNode = p2Node;
                    pParentNode = p2ParentNode;
                    pParent = p2Parent;
                    //pIndex = p2ParentIndex;
                    //labP = labP2;
                    //hueAngleP = hueAngleP2;
                }
             
                didMerge = true;                  
            }
           
            if (didMerge) {
                visited.add(originalP);
            }
        }
        
        log.info("size before=" + nBefore + " after=" + mergedMap.size());
        
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
        */
                   
        // making a debug image before tuning color conditionals
        Map<PairInt, Set<PairInt>> mergedPoints = new HashMap<PairInt, Set<PairInt>>();
        for (Entry<PairInt, Integer> entry : br.getPointIndexMap().entrySet()) {
            
            PairInt p = entry.getKey();
            
            Integer cellIndex = entry.getValue();
            PairInt cellCentroid = br.getBlobMedialAxes().getOriginalBlobXYCentroid(cellIndex.intValue());
            DisjointSet2Node<PairInt> cellNode = cellMap.get(cellCentroid);            
            DisjointSet2Node<PairInt> parentCellNode = disjointSetHelper.findSet(cellNode);  
            
            PairInt parentCellCentroid = parentCellNode.getMember();
            
            Set<PairInt> points = mergedPoints.get(parentCellCentroid);
            if (points == null) {
                points = new HashSet<PairInt>();
                mergedPoints.put(parentCellCentroid, points);
            }
            points.add(p);
        }
        
        Image imgCp = img.copyToGreyscale().copyToColorGreyscale();
        
        int nColors = mergedPoints.size();
        log.info("mergedPoints.size() = " + nColors);
        int delta = (int)Math.floor(256.f/(float)nColors);
        Random sr = null;
        long seed = System.currentTimeMillis();
        try {
            sr = SecureRandom.getInstance("SHA1PRNG");
            sr.setSeed(seed);
        } catch (NoSuchAlgorithmException e) {
            sr = new Random(seed);
        }
       
        Set<String> clrs = new HashSet<String>();
        for (Entry<PairInt, Set<PairInt>> entry : mergedPoints.entrySet()) {
            boolean alreadyChosen = true;
            int rClr = -1;
            int gClr = -1;
            int bClr = -1;
            while (alreadyChosen) {
                rClr = sr.nextInt(nColors)*delta;
                gClr = sr.nextInt(nColors)*delta;
                bClr = sr.nextInt(nColors)*delta;
                String str = "";
                if (rClr < 10) {
                    str = str + "00";
                } else if (rClr < 100) {
                    str = str + "0";
                }
                str = str + Integer.toString(rClr);
                if (gClr < 10) {
                    str = str + "00";
                } else if (gClr < 100) {
                    str = str + "0";
                }
                str = str + Integer.toString(bClr);
                if (bClr < 10) {
                    str = str + "00";
                } else if (bClr < 100) {
                    str = str + "0";
                }
                str = str + Integer.toString(bClr);
                alreadyChosen = clrs.contains(str);
                clrs.add(str);
            }
            for (PairInt p : entry.getValue()) {
                imgCp.setRGB(p.getX(), p.getY(), rClr, gClr, bClr);
            }
        }
        MiscDebug.writeImage(imgCp, "_" + debugTag + "_merged_");
        
        if (simWriter != null) {
            try {
                if (simWriter != null) {
                    simWriter.close();
                }
            } catch (IOException ex) {
                Logger.getLogger(SegmentedCellMerger.class.getName()).log(Level.SEVERE, null, ex);
            }
        }        
    }
    
    private BoundingRegions extractPerimetersAndBounds() {
                                
        // --- sort by descending sizes the remaining blobs ---- 
        int[] sizes = new int[segmentedCellList.size()];
        int[] indexes = new int[sizes.length];
        for (int i = 0; i < segmentedCellList.size(); ++i) {
            sizes[i] = segmentedCellList.get(i).size();
            indexes[i] = i;
        }
        
        // removing any blobs which are larger than 0.2 percent of image size also
        //float nPixels = img.getNPixels();
        
        MultiArrayMergeSort.sortByDecr(sizes, indexes);
        List<Set<PairInt>> tmp = new ArrayList<Set<PairInt>>();
        
        for (int i = 0; i < sizes.length; ++i) {
            int idx = indexes[i];
            Set<PairInt> blob = segmentedCellList.get(idx);
            //float frac = (float)blob.size()/nPixels;
            //if (frac < 0.2) {
                tmp.add(blob);
            //}
        }
        segmentedCellList.clear();
        segmentedCellList.addAll(tmp);
        
        //---- begin section to log colors to look at selecting matchable bounds by color ------
        CIEChromaticity cieC = new CIEChromaticity();
        //MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        List<Double> labLAvg = new ArrayList<Double>();
        List<Double> labAAvg = new ArrayList<Double>();
        List<Double> labBAvg = new ArrayList<Double>();
        List<Double> rAvg = new ArrayList<Double>();
        List<Double> gAvg = new ArrayList<Double>();
        List<Double> bAvg = new ArrayList<Double>();
        
        for (int i = 0; i < segmentedCellList.size(); ++i) {
            double redSum = 0;
            double greenSum = 0;
            double blueSum = 0;
            for (PairInt p : segmentedCellList.get(i)) {
                int x = p.getX();
                int y = p.getY();
                int red = img.getR(x, y);
                int green = img.getG(x, y);
                int blue = img.getB(x, y);
                redSum += red;
                greenSum += green;
                blueSum += blue;
            }
            double n = (double)segmentedCellList.get(i).size();
            redSum /= n;
            greenSum /= n;
            blueSum /= n;
            float[] avgLAB = cieC.rgbToCIELAB((int)Math.round(redSum), 
                (int)Math.round(greenSum), (int)Math.round(blueSum));
            
            labLAvg.add(Double.valueOf(avgLAB[0]));
            labAAvg.add(Double.valueOf(avgLAB[1]));
            labBAvg.add(Double.valueOf(avgLAB[2]));
            
            rAvg.add(Double.valueOf(redSum));
            gAvg.add(Double.valueOf(greenSum));
            bAvg.add(Double.valueOf(blueSum));
           
            /*PairInt xyCen = curveHelper.calculateXYCentroids(blobs.get(i));
            String str = String.format(
                "[%d] cen=(%d,%d) avgL=%.3f avgA=%.3f  avgB=%.3f  nPts=%d",
                i, xyCen.getX(), xyCen.getY(),
                avgLAB[0], avgLAB[1], avgLAB[2], blobs.get(i).size());
            log.info(str);*/
        }
        
        Map<PairInt, Integer> blobPointToListIndex = createBlobPointToListIndex(
            segmentedCellList);
                                
        // less than O(N)
        List<Set<PairInt>> borderPixelSets = BlobsAndPerimeters
            .extractBlobPerimeterAsPoints(segmentedCellList, img.getWidth(), 
            img.getHeight());
        
        assert(segmentedCellList.size() == borderPixelSets.size());
        
        List<PairIntArray> perimetersList = new ArrayList<PairIntArray>();
        
        float srchRadius = (float)Math.sqrt(2);
        
        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        //ImageSegmentation imageSegmentation = new ImageSegmentation();
                
        BlobMedialAxes bma = new BlobMedialAxes(segmentedCellList, labLAvg, labAAvg, labBAvg, 
            rAvg, gAvg, bAvg);
        
        for (int i = 0; i < borderPixelSets.size(); ++i) {
                                    
            Set<PairInt> blob = segmentedCellList.get(i);
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
            //imageSegmentation.makeStraightLinesHollow(orderedPerimeter, img.getWidth(), 
            //    img.getHeight(), srchRadius);
        
            /*
            Image imgCp = img.copyImage();
            ImageIOHelper.addCurveToImage(orderedPerimeter, imgCp, 2, 255, 0, 0);                  
            MiscDebug.writeImage(imgCp, "_" + i + "_" + MiscDebug.getCurrentTimeFormatted()); 
            */
            
            perimetersList.add(orderedPerimeter);
        }
        
        BoundingRegions br = new BoundingRegions(perimetersList, bma, 
            blobPointToListIndex);
        
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
        
        //debug(br, adjacencyMap);
        
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

    private double[][] getLDASegmentationMatrix() {
     
        double[][] m = new double[2][4];
        for (int i = 0; i < 2; ++i) {
            m[i] = new double[4];
        }
        
        m[0][0] = 0.343;
        m[1][0] = -0.482;
        m[0][1] = -0.575;
        m[1][1] = -0.381;
        m[0][2] = -0.439;
        m[1][2] = 0.771;
        m[0][3] = -0.599;
        m[1][3] = 0.166;
        
        return m;
    }

    private void debug(BoundingRegions br, Map<PairInt, Set<PairInt>> adjacencyMap) {
        
        Image imgCp = img.copyImage();
        Image imgCp2 = imgCp.copyToGreyscale().copyToColorGreyscale();
        
        List<PairIntArray> boundaries = br.getPerimeterList();
        
        int n = br.getBlobMedialAxes().getNumberOfItems();
        
        int nExtraForDot = 0;
        
        for (int i = 0; i < n; ++i) {
            
            PairIntArray bounds = boundaries.get(i);
            
            int[] clr = ImageIOHelper.getNextRGB(i);
          
            //-- draw boundaries as a curve --
            int nb = bounds.getN();
            int[] xVertexes = new int[nb + 1];
            int[] yVertexes = new int[nb + 1];
            for (int j = 0; j < nb; ++j) {
                xVertexes[j] = bounds.getX(j);
                yVertexes[j] = bounds.getY(j);
            }
            xVertexes[nb] = bounds.getX(0);
            yVertexes[nb] = bounds.getY(0);
            ImageIOHelper.drawLinesInImage(xVertexes, yVertexes, imgCp,
                nExtraForDot, clr[0], clr[1], clr[2]);
        }
        
        nExtraForDot = 1;
        
        //for (int i = 0; i < n; ++i) {
        for (int i = 1; i < 2; ++i) {
            
            PairIntArray bounds = boundaries.get(i);
            
            int[] clr = ImageIOHelper.getNextRGB(i);
          
            PairInt cenXY = br.getBlobMedialAxes().getOriginalBlobXYCentroid(i);
            
            ImageIOHelper.addPointToImage(cenXY.getX(), cenXY.getY(), imgCp, 
                0, clr[0], clr[1], clr[2]);
            
            /*
            -- draw boundaries as a curve --
            */
            int nb = bounds.getN();
            int[] xVertexes = new int[nb + 1];
            int[] yVertexes = new int[nb + 1];
            for (int j = 0; j < nb; ++j) {
                xVertexes[j] = bounds.getX(j);
                yVertexes[j] = bounds.getY(j);
            }
            xVertexes[nb] = bounds.getX(0);
            yVertexes[nb] = bounds.getY(0);
            ImageIOHelper.drawLinesInImage(xVertexes, yVertexes, imgCp,
                nExtraForDot, clr[0], clr[1], clr[2]);
            
            ImageIOHelper.addPointToImage(cenXY.getX(), cenXY.getY(), imgCp2, 
                nExtraForDot, clr[0], clr[1], clr[2]);
            
            Set<PairInt> adjacentCenXY = adjacencyMap.get(cenXY);
            if (adjacentCenXY != null) {
                for (PairInt cenXY2 : adjacentCenXY) {
                    ImageIOHelper.addPointToImage(cenXY2.getX(), cenXY2.getY(), imgCp2, 
                        nExtraForDot, clr[0], clr[1], clr[2]);
                    ImageIOHelper.addPointToImage(cenXY2.getX(), cenXY2.getY(), imgCp, 
                        0, clr[0], clr[1], clr[2]);
                }
            }
        }
        
        MiscDebug.writeImage(imgCp, "adj_0_" + debugTag);
        MiscDebug.writeImage(imgCp2, "adj_1_" + debugTag);
    }

    private Map<PairInt, Integer> createBlobPointToListIndex(List<Set<PairInt>> 
        blobs) {
        
        Map<PairInt, Integer> map = new HashMap<PairInt, Integer>();
        
        for (int i = 0; i < blobs.size(); ++i) {
            
            Set<PairInt> blob = blobs.get(i);
            
            Integer key = Integer.valueOf(i);
            
            for (PairInt p : blob) {
                map.put(p, key);
            }
        }
        
        return map;
    }

    private void writeLabelFile(BoundingRegions br, Map<PairInt, Integer> 
        cellIndexMap) {
        
        List<Set<PairInt>> centroidList = new ArrayList<Set<PairInt>>();
        for (int i = 0; i < br.getBlobMedialAxes().getNumberOfItems(); ++i) {
            PairInt xyCen = br.getBlobMedialAxes().getOriginalBlobXYCentroid(i);
            Set<PairInt> set = new HashSet<PairInt>();
            set.add(xyCen);
            centroidList.add(set);
        }
        
        NearestPointsInLists np = new NearestPointsInLists(centroidList);
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        float radius = 15.f;
        
        for (PairIntPair pp : simClass) {
            
            PairInt p0 = np.findClosest(pp.getX1(), pp.getY1(), radius);
            
            if (p0 == null) {
                log.severe("revise coords.  could not find x=" + pp.getX1() + 
                    ", y=" + pp.getY1());
                continue;
            }
            
            PairInt p1 = np.findClosest(pp.getX2(), pp.getY2(), radius);
            
            if (p1 == null) {
                log.severe("revise coords.  could not find x=" + pp.getX2() + 
                    ", y=" + pp.getY2());
                continue;
            }
            
            int list0Idx = cellIndexMap.get(p0).intValue();
            
            int list1Idx = cellIndexMap.get(p1).intValue();
            
            double r0 = br.getBlobMedialAxes().getR(list0Idx);
            double g0 = br.getBlobMedialAxes().getG(list0Idx);
            double b0 = br.getBlobMedialAxes().getB(list0Idx);
            float[] lab0 = br.getBlobMedialAxes().getLABColors(list0Idx);
            
            double r1 = br.getBlobMedialAxes().getR(list1Idx);
            double g1 = br.getBlobMedialAxes().getG(list1Idx);
            double b1 = br.getBlobMedialAxes().getB(list1Idx);
            float[] lab1 = br.getBlobMedialAxes().getLABColors(list1Idx);
            
            double absDeltaE = Math.abs(cieC.calcDeltaECIE94(lab0, lab1));
            double absDeltaL = Math.abs(lab0[0] - lab1[0]);
            double absDeltaBR = Math.abs((b0 - r0) - (b1 - r1));
            double absDeltaBG = Math.abs((b0 - g0) - (b1 - g1));
            
            double absDeltaO1 = Math.abs((double)((r0 - g0) - (r1 - g1))/Math.sqrt(2.));
            double absDeltaO2 = Math.abs((double)((r0 + g0 - 2*b0) - (r1 + g1 - 2*b1))/Math.sqrt(6.));
            double absDeltaO3 = Math.abs((double)((r0 + g0 + b0) - (r1 + g1 + b1))/Math.sqrt(2.));
            
            String str = String.format(
                "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,sim\n",
                (float)absDeltaE,
                (float)absDeltaO1,
                (float)absDeltaBR,
                (float)absDeltaBG,
                (float)absDeltaO2,
                (float)absDeltaO3,
                (float)absDeltaL
                );
            
            try {
                simWriter.write(str);
            } catch (IOException ex) {
                Logger.getLogger(SegmentedCellMerger.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        try {
            simWriter.flush();
        } catch (IOException ex) {
            Logger.getLogger(SegmentedCellMerger.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        for (PairIntPair pp : diffClass) {
            
            PairInt p0 = np.findClosest(pp.getX1(), pp.getY1(), radius);
            
            if (p0 == null) {
                log.severe("revise coords.  could not find x=" + pp.getX1() + 
                    ", y=" + pp.getY1());
                continue;
            }
            
            PairInt p1 = np.findClosest(pp.getX2(), pp.getY2(), radius);
            
            if (p1 == null) {
                log.severe("revise coords.  could not find x=" + pp.getX2() + 
                    ", y=" + pp.getY2());
                continue;
            }
            
            int list0Idx = cellIndexMap.get(p0).intValue();
            
            int list1Idx = cellIndexMap.get(p1).intValue();
            
            double r0 = br.getBlobMedialAxes().getR(list0Idx);
            double g0 = br.getBlobMedialAxes().getG(list0Idx);
            double b0 = br.getBlobMedialAxes().getB(list0Idx);
            float[] lab0 = br.getBlobMedialAxes().getLABColors(list0Idx);
            
            double r1 = br.getBlobMedialAxes().getR(list1Idx);
            double g1 = br.getBlobMedialAxes().getG(list1Idx);
            double b1 = br.getBlobMedialAxes().getB(list1Idx);
            float[] lab1 = br.getBlobMedialAxes().getLABColors(list1Idx);
            
            double absDeltaE = Math.abs(cieC.calcDeltaECIE94(lab0, lab1));
            double absDeltaL = Math.abs(lab0[0] - lab1[0]);
            double absDeltaBR = Math.abs((b0 - r0) - (b1 - r1));
            double absDeltaBG = Math.abs((b0 - g0) - (b1 - g1));
            
            double absDeltaO1 = Math.abs((double)((r0 - g0) - (r1 - g1))/Math.sqrt(2.));
            double absDeltaO2 = Math.abs((double)((r0 + g0 - 2*b0) - (r1 + g1 - 2*b1))/Math.sqrt(6.));
            double absDeltaO3 = Math.abs((double)((r0 + g0 + b0) - (r1 + g1 + b1))/Math.sqrt(2.));
            
            String str = String.format(
                "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,diff\n",
                (float)absDeltaE,
                (float)absDeltaO1,
                (float)absDeltaBR,
                (float)absDeltaBG,
                (float)absDeltaO2,
                (float)absDeltaO3,
                (float)absDeltaL
                );
            
            try {
                simWriter.write(str);
            } catch (IOException ex) {
                Logger.getLogger(SegmentedCellMerger.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        try {
            simWriter.flush();
        } catch (IOException ex) {
            Logger.getLogger(SegmentedCellMerger.class.getName()).log(Level.SEVERE, null, ex);
        }       
    }
    
}
