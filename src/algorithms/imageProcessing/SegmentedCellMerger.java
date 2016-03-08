package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.NearestPointsInLists;
import algorithms.disjointSets.DisjointSet2Helper;
import algorithms.disjointSets.DisjointSet2Node;
import algorithms.imageProcessing.ImageProcessor.Colors;
import algorithms.imageProcessing.ImageSegmentation.BoundingRegions;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import algorithms.util.PairIntPair;
import algorithms.util.ResourceFinder;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class SegmentedCellMerger {
 
    private final ImageExt img;
    private final List<Set<PairInt>> segmentedCellList;
    private final String debugTag;
    private final boolean useDeltaE2000;
    private final float deltaELimit;
    
    private final Logger log = Logger.getLogger(this.getClass().getName());
    
    public SegmentedCellMerger(ImageExt img, List<Set<PairInt>> theSegmentedCellList, 
        boolean useDeltaE2000, float deltaELimit, String debugTag) {
        
        this.useDeltaE2000 = useDeltaE2000;
        this.deltaELimit = deltaELimit;
        
        this.img = img;
        
        // sort theSegmentedCellList by size
        int[] indexes = new int[theSegmentedCellList.size()];
        int[] sizes = new int[indexes.length];
        for (int i = 0; i < theSegmentedCellList.size(); ++i) {
            indexes[i] = i;
            sizes[i] = theSegmentedCellList.get(i).size();
        }
        MultiArrayMergeSort.sortByDecr(sizes, indexes);
        
        this.segmentedCellList = new ArrayList<Set<PairInt>>(theSegmentedCellList.size());
        for (int i = 0; i < theSegmentedCellList.size(); ++i) {
            int idx = indexes[i];
            Set<PairInt> set = theSegmentedCellList.get(idx);
            Set<PairInt> set2 = new HashSet<PairInt>(set);
            segmentedCellList.add(set2);
        }
        
        this.debugTag = (debugTag != null) ? debugTag : "";
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
        
        long t0 = System.currentTimeMillis();
        
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
                
        // key = pairintwithindex holding xy centroid and the original list index
        final Map<PairIntWithIndex, DisjointSet2Node<PairIntWithIndex>> cellMap = 
            new HashMap<PairIntWithIndex, DisjointSet2Node<PairIntWithIndex>>();
        
        final Map<PairInt, Integer> pointIndexMap = new HashMap<PairInt, Integer>();
        
        final Map<PairIntWithIndex, Integer> cellIndexMap = new HashMap<PairIntWithIndex, Integer>();
        
        final Map<Integer, PairIntWithIndex> indexCellMap = new HashMap<Integer, PairIntWithIndex>();
        
        // cellMap is only initialized w/ the key's own centroid here
        populateCellMaps(cellMap, indexCellMap, cellIndexMap, pointIndexMap);
        
//writeLabelFile(br, cellIndexMap);
        
        // key = cell centroid, value = set of adjacent cell centroids
        final Map<PairIntWithIndex, Set<PairIntWithIndex>> adjacencyMap 
            = new HashMap<PairIntWithIndex, Set<PairIntWithIndex>>();
        
        populateAdjacencyMap(indexCellMap, pointIndexMap, adjacencyMap);

        // key = cell centroid, value = set of cell centroids merged with this one
        Map<PairIntWithIndex, Set<PairIntWithIndex>> mergedMap = 
            initializeMergeMap(cellIndexMap);
        
        Map<Integer, Colors> segmentedCellAvgLabColors 
            = new HashMap<Integer, Colors>();
        
        // this order results in visiting the largest cells first
        Stack<PairIntWithIndex> stack = new Stack<PairIntWithIndex>();
        for (int i = (segmentedCellList.size() - 1); i > -1; --i) {            
            PairIntWithIndex key = indexCellMap.get(Integer.valueOf(i));
            if (cellIndexMap.containsKey(key)) {
                stack.add(key);
            }            
        }
        
        int nBefore = mergedMap.size();
                
        CIEChromaticity cieC = new CIEChromaticity();
        
        Set<PairIntWithIndex> visited = new HashSet<PairIntWithIndex>();
        Map<PairIntWithIndex, Set<PairIntWithIndex>> visitedMap 
            = new HashMap<PairIntWithIndex, Set<PairIntWithIndex>>();
        
        DisjointSet2Helper disjointSetHelper = new DisjointSet2Helper();
                
        ImageProcessor imageProcessor = new ImageProcessor();
        
        //TODO: improving this transformation
        //double[][] ldaMatrix = getLDASegmentationMatrix();
          
        while (!stack.isEmpty()) {
            
            PairIntWithIndex p = stack.pop();
            
            if (visited.contains(p)) {
                continue;
            }
            
            PairIntWithIndex originalP = (PairIntWithIndex) p.copy();
            //Integer originalPIndex = cellIndexMap.get(originalP);
            
            DisjointSet2Node<PairIntWithIndex> pNode = cellMap.get(p);
                        
            DisjointSet2Node<PairIntWithIndex> pParentNode = disjointSetHelper.findSet(pNode);
            
            PairIntWithIndex pParent = pParentNode.getMember();
            
            Integer parentIndex = cellIndexMap.get(pParent);
            
            Colors colorsP = segmentedCellAvgLabColors.get(parentIndex);
            if (colorsP == null) {
                colorsP = imageProcessor.calculateAverageLAB(img, 
                    segmentedCellList.get(parentIndex.intValue()));
                segmentedCellAvgLabColors.put(parentIndex, colorsP);
            }
            float[] labP = colorsP.getColors();
            
            boolean didMerge = false;

            List<PairIntWithIndex> neighborKeys = new ArrayList<PairIntWithIndex>(
                adjacencyMap.get(pParent));
            
            for (PairIntWithIndex p2 : neighborKeys) {
                
                DisjointSet2Node<PairIntWithIndex> p2Node = cellMap.get(p2);
                
                // find the "representative" parent node of the cell
                DisjointSet2Node<PairIntWithIndex> p2ParentNode = disjointSetHelper.findSet(p2Node);
                    
                if (p2ParentNode.equals(pParentNode)) {
                    continue;
                }
                    
                PairIntWithIndex p2Parent = p2ParentNode.getMember();

                if (hasBeenVisited(visitedMap, pParent, p2Parent)) {
                    continue;
                }
                
                //Integer p2Index = cellIndexMap.get(p2);
                
                Integer p2ParentIndex = cellIndexMap.get(p2Parent);
                
                Colors colors2 = segmentedCellAvgLabColors.get(p2ParentIndex);
                if (colors2 == null) {
                    colors2 = imageProcessor.calculateAverageLAB(img,
                        segmentedCellList.get(p2ParentIndex.intValue()));
                    segmentedCellAvgLabColors.put(p2ParentIndex, colors2);
                }
                float[] labP2 = colors2.getColors();
                           
                double deltaE;
                if (useDeltaE2000) {
                    deltaE = Math.abs(cieC.calcDeltaECIE2000(labP, labP2));
                } else {
                    deltaE = Math.abs(cieC.calcDeltaECIE94(labP, labP2));
                }
                               
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
                
                if (Math.abs(deltaE) > deltaELimit) {
                    continue;
                }
    
                stack.add(p2Parent);
                
                DisjointSet2Node<PairIntWithIndex> parentOfMergeNode = 
                    disjointSetHelper.union(pParentNode, p2ParentNode);

                PairIntWithIndex parentOfMerge = parentOfMergeNode.getMember();

                if (parentOfMerge.equals(pParent)) {
                    
                    // parent is still the main comparator
                    
                    //mergedMap for key pParent gets values of p2Parent 
                    //and the p2Parent key gets removed
                    
                    Set<PairIntWithIndex> p2ParentMergedSet = mergedMap.get(p2Parent);
                    
                    Set<PairIntWithIndex> pParentMergedSet = mergedMap.get(pParent);
                    pParentMergedSet.addAll(p2ParentMergedSet);
                    
                    mergedMap.remove(p2Parent);
                    
                    // adjacencyMap for key pParent gets values of p2Parent, 
                    // and then subtracts all internal members from the final set
                    Set<PairIntWithIndex> p2ParentAdjacencySet = adjacencyMap.get(p2Parent);
                    
                    Set<PairIntWithIndex> pParentAdjacencySet = adjacencyMap.get(pParent);
                    pParentAdjacencySet.addAll(p2ParentAdjacencySet);
                    pParentAdjacencySet.removeAll(pParentMergedSet);
                    
                } else {

                    // the new parent key has become a neighbor, so after the
                    // merges, have to re-assign the local variables within
                    // the upper block
                    
                    //mergedMap for key p2Parent gets values of pParent 
                    //and the pParent key gets removed
                    
                    Set<PairIntWithIndex> pParentMergedSet = mergedMap.get(pParent);
                    
                    Set<PairIntWithIndex> p2ParentMergedSet = mergedMap.get(p2Parent);
                    p2ParentMergedSet.addAll(pParentMergedSet);
                    
                    mergedMap.remove(pParent);
                    
                    // adjacencyMap for key p2Parent gets values of pParent, 
                    // and then subtracts all internal members from the final set
                    Set<PairIntWithIndex> pParentAdjacencySet = adjacencyMap.get(pParent);
                    
                    Set<PairIntWithIndex> p2ParentAdjacencySet = adjacencyMap.get(p2Parent);
                    p2ParentAdjacencySet.addAll(pParentAdjacencySet);
                    p2ParentAdjacencySet.removeAll(p2ParentMergedSet);
   
                    // --- reassign upper block variables to use new parent information ---
                    p = p2;
                    
                    // NOTE: there's some inefficiency in that the visited check
                    // is now bypassed, but the number of original neighbors is
                    // small.  the possible extra work should be idempotent.
                    
                    pNode = p2Node;
                    pParentNode = p2ParentNode;
                    pParent = p2Parent;
                    parentIndex = p2ParentIndex;
                    labP = labP2;
                    //pIndex = p2ParentIndex;
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
        
        the cell indexes can then be gathered and then new blobs cells made from
            the points in those sets.
        */
                   
        // making a debug image before tuning color conditionals
        Map<PairIntWithIndex, Set<PairIntWithIndex>> mergedPoints 
            = new HashMap<PairIntWithIndex, Set<PairIntWithIndex>>();
        
        for (Entry<PairInt, Integer> entry : pointIndexMap.entrySet()) {
            
            Integer cellIndex = entry.getValue();
            
            PairIntWithIndex cellCentroid = indexCellMap.get(cellIndex);
           
            DisjointSet2Node<PairIntWithIndex> cellNode = cellMap.get(cellCentroid);            
            DisjointSet2Node<PairIntWithIndex> parentCellNode = disjointSetHelper.findSet(cellNode);  
            
            PairIntWithIndex parentCellCentroid = parentCellNode.getMember();
            
            Set<PairIntWithIndex> points = mergedPoints.get(parentCellCentroid);
            if (points == null) {
                points = new HashSet<PairIntWithIndex>();
                mergedPoints.put(parentCellCentroid, points);
            }
            
            PairIntWithIndex p = new PairIntWithIndex(entry.getKey(), 
                cellIndex.intValue());
            points.add(p);
        }
        
        int count = 0;
        segmentedCellList.clear();
        for (Entry<PairIntWithIndex, Set<PairIntWithIndex>> entry 
            : mergedPoints.entrySet()) {
            Set<PairInt> set = new HashSet<PairInt>();
            for (PairIntWithIndex p : entry.getValue()) {
                set.add(new PairInt(p.getX(), p.getY()));
            }
            count += set.size();
            segmentedCellList.add(set);
        }
        
        long t1 = System.currentTimeMillis();
        long t1Sec = (t1 - t0)/1000;
        log.info(debugTag + " " + t1Sec + " sec to merge " + mergedPoints.size() 
            + " cells in scm.  list size=" + segmentedCellList.size());
        
        log.info("pixels in cell lists=" + count + ", pixels in image=" 
            + img.getNPixels());
            
        try {
            if (simWriter != null) {
                simWriter.close();
            }
        } catch (IOException ex) {
            Logger.getLogger(SegmentedCellMerger.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public List<Set<PairInt>> getSegmentedCellList() {
        return segmentedCellList;
    }
    
    private void populateCellMaps( 
        Map<PairIntWithIndex, DisjointSet2Node<PairIntWithIndex>> outputParentMap,
        Map<Integer, PairIntWithIndex> outputIndexCellMap,
        Map<PairIntWithIndex, Integer> outputParentIndexMap,
        Map<PairInt, Integer> pointIndexMap) {
        
        DisjointSet2Helper disjointSetHelper = new DisjointSet2Helper();
        
        // store the centroids in a disjoint set/forest
        
        int n = this.segmentedCellList.size();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        // init map
        
        for (int i = 0; i < n; ++i) {
            
            Set<PairInt> set = segmentedCellList.get(i);
                        
            PairInt xyCen = curveHelper.calculateXYCentroids2(set);
            
            PairIntWithIndex pai = new PairIntWithIndex(xyCen, i);
            
            DisjointSet2Node<PairIntWithIndex> pNode =
                disjointSetHelper.makeSet(new DisjointSet2Node<PairIntWithIndex>(pai));
            
            outputParentMap.put(pai, pNode);
            
            Integer index = Integer.valueOf(i);
            
            outputParentIndexMap.put(pai, index);
            
            outputIndexCellMap.put(index, pai);
            
            assert(pNode.getMember().equals(pai));
            assert(pNode.getParent().equals(pNode));
            
            for (PairInt p : set) {
                pointIndexMap.put(p, index);
            }
        }   
        
    }

    private boolean hasBeenVisited(Map<PairIntWithIndex, Set<PairIntWithIndex>> 
        visitedMap, PairIntWithIndex pParent, PairIntWithIndex p2Parent) {
        
        Set<PairIntWithIndex> set = visitedMap.get(pParent);
        if (set == null) {
            return false;
        }
        
        return set.contains(p2Parent);
    }

    private void addToVisited(Map<PairIntWithIndex, Set<PairIntWithIndex>> 
        visitedMap, PairIntWithIndex pParent, PairIntWithIndex p2Parent) {
        
        Set<PairIntWithIndex> set = visitedMap.get(pParent);
        if (set == null) {
            set = new HashSet<PairIntWithIndex>();
            visitedMap.put(pParent, set);
        }
        set.add(p2Parent);
        
        set = visitedMap.get(p2Parent);
        if (set == null) {
            set = new HashSet<PairIntWithIndex>();
            visitedMap.put(p2Parent, set);
        }
        set.add(pParent);
    }

    private Map<PairIntWithIndex, Set<PairIntWithIndex>> initializeMergeMap(
        Map<PairIntWithIndex, Integer> cellIndexMap) {
        
        Map<PairIntWithIndex, Set<PairIntWithIndex>> mergeMap 
            = new HashMap<PairIntWithIndex, Set<PairIntWithIndex>>(cellIndexMap.size());
                        
        for (Entry<PairIntWithIndex, Integer> entry : cellIndexMap.entrySet()) {
                        
            PairIntWithIndex key = entry.getKey();
           
            Set<PairIntWithIndex> mergeSet = new HashSet<PairIntWithIndex>();
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
            
            double absDeltaO1 = Math.abs(((r0 - g0) - (r1 - g1))/Math.sqrt(2.));
            double absDeltaO2 = Math.abs(((r0 + g0 - 2*b0) - (r1 + g1 - 2*b1))/Math.sqrt(6.));
            double absDeltaO3 = Math.abs(((r0 + g0 + b0) - (r1 + g1 + b1))/Math.sqrt(2.));
            
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
            
            double absDeltaO1 = Math.abs(((r0 - g0) - (r1 - g1))/Math.sqrt(2.));
            double absDeltaO2 = Math.abs(((r0 + g0 - 2*b0) - (r1 + g1 - 2*b1))/Math.sqrt(6.));
            double absDeltaO3 = Math.abs(((r0 + g0 + b0) - (r1 + g1 + b1))/Math.sqrt(2.));
            
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

    private void populateAdjacencyMap(Map<Integer, PairIntWithIndex> indexCellMap, 
        Map<PairInt, Integer> pointIndexMap, 
        Map<PairIntWithIndex, Set<PairIntWithIndex>> adjacencyMap) {
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        int n = this.segmentedCellList.size();
                        
        for (int i = 0; i < n; ++i) {
            
            Set<PairInt> set = segmentedCellList.get(i);
            
            Integer index = Integer.valueOf(i);
            
            PairIntWithIndex xyCen = indexCellMap.get(index);
            
            Set<PairIntWithIndex> adjSet = adjacencyMap.get(xyCen);
            if (adjSet == null) {
                adjSet = new HashSet<PairIntWithIndex>();
                adjacencyMap.put(xyCen, adjSet);
            }
            
            // look for neighbors that are in other lists
            for (PairInt p : set) {
                int x = p.getX();
                int y = p.getY();
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    PairInt p2 = new PairInt(x2, y2);
                    Integer index2 = pointIndexMap.get(p2);
                    if (index2 == null || index.equals(index2)) {
                        continue;
                    }
                    
                    PairIntWithIndex xyCen2 = indexCellMap.get(index2);
                    assert(!xyCen.equals(xyCen2));
            
                    Set<PairIntWithIndex> adjSet2 = adjacencyMap.get(xyCen2);
                    if (adjSet2 == null) {
                        adjSet2 = new HashSet<PairIntWithIndex>();
                        adjacencyMap.put(xyCen2, adjSet2);
                    }
                    
                    adjSet.add(xyCen2);
                    adjSet2.add(xyCen);
                }
            }
        }
    }
    
    private static class PairIntWithIndex extends com.climbwithyourfeet.clustering.util.PairInt {
        int pixIdx;
        public PairIntWithIndex(int xPoint, int yPoint, int thePixIndex) {
            super(xPoint, yPoint);
            pixIdx = thePixIndex;
        }
        public PairIntWithIndex(algorithms.util.PairInt p, int thePixIndex) {
            super(p.getX(), p.getY());
            pixIdx = thePixIndex;
        }
        public int getPixIndex() {
            return pixIdx;
        }
        @Override
        public boolean equals(Object obj) {
            if (!(obj instanceof PairIntWithIndex)) {
                return false;
            }
            PairIntWithIndex other = (PairIntWithIndex) obj;
            return (x == other.getX()) && (y == other.getY()) &&
                (pixIdx == other.pixIdx);
        }
        @Override
        public int hashCode() {
            int hash = fnvHashCode(this.x, this.y, this.pixIdx);
            return hash;
        }
        @Override
        public com.climbwithyourfeet.clustering.util.PairInt copy() {
            return new PairIntWithIndex(x, y, pixIdx);
        }        
        protected int fnvHashCode(int i0, int i1, int i2) {
            /*
             * hash = offset_basis
             * for each octet_of_data to be hashed
             *     hash = hash xor octet_of_data
             *     hash = hash * FNV_prime
             * return hash
             *
             * Public domain:  http://www.isthe.com/chongo/src/fnv/hash_32a.c
             */
            int hash = 0;
            int sum = fnv321aInit;
            // xor the bottom with the current octet.
            sum ^= i0;
            // multiply by the 32 bit FNV magic prime mod 2^32
            sum *= fnv32Prime;
            sum ^= i1;
            sum *= fnv32Prime;            
            sum ^= i2;
            sum *= fnv32Prime;
            hash = sum;
            return hash;
        }
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder(super.toString());
            sb.append(" pixIdx=").append(Integer.toString(pixIdx));
            return sb.toString();
        }
    }
}
