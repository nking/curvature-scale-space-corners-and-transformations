package algorithms.misc;

import algorithms.CountingSort;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PointPairInt;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class MiscDebug {
    
    private static Logger log = Logger.getLogger(MiscDebug.class.getName());
    
      public static void printJoinPoints(Map<PairInt, PairInt> joinPoints,
        List<PairIntArray> edges) {
        
        StringBuilder sb = new StringBuilder("join points:\n");
        
        for (Map.Entry<PairInt, PairInt> entry : joinPoints.entrySet()) {
            
            PairInt loc0 = entry.getKey();
            
            PairInt loc1 = entry.getValue();
            
            PairIntArray edge0 = edges.get(loc0.getX());
            int x0 = edge0.getX(loc0.getY());
            int y0 = edge0.getY(loc0.getY());
            
            PairIntArray edge1 = edges.get(loc1.getX());
            int x1 = edge1.getX(loc1.getY());
            int y1 = edge1.getY(loc1.getY());
            
            sb.append(String.format("  (%d,%d) to (%d,%d) in edges %d and %d  at positions=%d out of %d and %d out of %d\n",
                x0, y0, x1, y1, loc0.getX(), loc1.getX(), 
                loc0.getY(), edges.get(loc0.getX()).getN(),
                loc1.getY(), edges.get(loc1.getX()).getN()
            ));
        }
        
        log.info(sb.toString());
    }
    
    public static void printJunctions(Map<Integer, Set<Integer>> jMap, 
        List<PairIntArray> edges, GreyscaleImage img) {
        
        try {
            
            Image img2 = img.copyImageToGreen();

            ImageIOHelper.addAlternatingColorCurvesToImage(edges, img2);
            int nExtraForDot = 1;
            int rClr = 255;
            int gClr = 0;
            int bClr = 100;
            for (Map.Entry<Integer, Set<Integer>> entry : jMap.entrySet()) {
                int pixIdx = entry.getKey().intValue();
                int col = img2.getCol(pixIdx);
                int row = img2.getRow(pixIdx);
                ImageIOHelper.addPointToImage(col, row, img2, nExtraForDot,
                    rClr, gClr, bClr);
            }

            String dirPath = algorithms.util.ResourceFinder.findDirectory("bin");
            String sep = System.getProperty("file.separator");
            ImageIOHelper.writeOutputImage(dirPath + sep + "junctions.png", img2);
            
        } catch (IOException e) {
            
        }
    }
    
     public static void printJoinPoints(PairInt[][] edgeJoins, int idxLo, int idxHi,
        Map<Integer, PairIntArray> edges) {
        
        // print array indexes
        int[] aIndexes = new int[edges.size()];
        int count = 0;
        for (Map.Entry<Integer, PairIntArray> entry : edges.entrySet()) {
            Integer edgeIndex = entry.getKey();
            aIndexes[count] = edgeIndex.intValue();
            count++;
        }
        CountingSort.sort(aIndexes, 2*edges.size());
        
        StringBuilder sb = new StringBuilder("output indexes of size ");
        sb.append(Integer.toString(edges.size())).append("\n");
        for (int i = 0; i < aIndexes.length; i++) {
            sb.append(Integer.toString(aIndexes[i])).append(" ");
        }
        log.info(sb.toString());
        
        sb = new StringBuilder("join points:\n");
        
        for (int i = idxLo; i <= idxHi; i++) {
                        
            PairInt loc0 = edgeJoins[i][0];
            PairInt loc1 = edgeJoins[i][1];

            PairIntArray edge0 = edges.get(Integer.valueOf(loc0.getX()));
                  
            int x0 = edge0.getX(loc0.getY());
            int y0 = edge0.getY(loc0.getY());
            
            PairIntArray edge1 = edges.get(Integer.valueOf(loc1.getX()));
            int x1 = edge1.getX(loc1.getY());
            int y1 = edge1.getY(loc1.getY());
            
            sb.append(String.format("  (%d,%d) to (%d,%d) in edges %d and %d  at positions=%d out of %d and %d out of %d\n",
                x0, y0, x1, y1, loc0.getX(), loc1.getX(), 
                loc0.getY(), edge0.getN(),
                loc1.getY(), edge1.getN()
            ));
        }
        
        log.info(sb.toString());
    }
    
    public static void assertConsistentEdgeCapacity(
        Map<Integer, Set<Integer>> theEdgeToPixelIndexMap, 
        Map<Integer, PairInt> junctionLocationMap, 
        Map<Integer, Set<Integer>> junctionMap, 
        List<PairIntArray> edges) {
        
        for (Map.Entry<Integer, PairInt> entry : junctionLocationMap.entrySet()) {
            Integer pixelIndex = entry.getKey();
            assert(pixelIndex != null);
            PairInt loc = entry.getValue();
            int edgeIdx = loc.getX();
            int idx = loc.getY();
            PairIntArray edge = edges.get(edgeIdx);
            assert(edge != null);
            assert (idx < edge.getN()) :
                "idx=" + idx + " edgeN=" + edge.getN() + " edgeIndex=" + edgeIdx;
        }
        /*
        previous edge 50 had n=24
        current edge 49 has n=24.
        
        java.lang.AssertionError: idx=28 edgeN=24 edgeIndex=49
        
        
        processing junction w/ center pixel index=50427 and loc=3:15
        ...[junit] edge=49 idx=28 (out of 39) pixIdx=42810 (13,247)
        
        before splice edge 4 (137 points) to edge 3 (16 points)
        
        -----------
        splice edge 49 (29 points) to edge 50 (8 points), that is append 50 to end of 49
        50 edge size is 8
        49 edge size is 49.  spliced to 29 and 20
          'off by 1'?  spliced last point is idx=28 within edge 49 before...
        
        splice0_0: update Y for pixIdx=42810   loc 49:28 to 49:28 (edgeN=29)  (**make sure edgeIdx is same)
        
        */
        
        for (Map.Entry<Integer, Set<Integer>> entry : junctionMap.entrySet()) {
            
            Integer pixelIndex = entry.getKey();
            Set<Integer> adjPixelIndexes = entry.getValue();
            
            PairInt loc = junctionLocationMap.get(pixelIndex);
            assert(loc != null);
            
            int edgeIdx = loc.getX();
            int idx = loc.getY();
            PairIntArray edge = edges.get(edgeIdx);
            assert(edge != null);
            assert(idx < edge.getN());
            
            for (Integer adjPixelIndex : adjPixelIndexes) {
                PairInt adjLoc = junctionLocationMap.get(adjPixelIndex);
                assert(adjLoc != null);
                edgeIdx = adjLoc.getX();
                idx = adjLoc.getY();
                edge = edges.get(edgeIdx);
                assert(edge != null);
                assert(idx < edge.getN());
            }
        }
        
        for (Map.Entry<Integer, Set<Integer>> entry : theEdgeToPixelIndexMap.entrySet()) {
            Integer edgeIndex = entry.getKey();
            Set<Integer> pixelIndexes = entry.getValue();
            
            PairIntArray edge = edges.get(edgeIndex.intValue());
            assert(edge != null);
            
            for (Integer pixelIndex : pixelIndexes) {
                PairInt loc = junctionLocationMap.get(pixelIndex);
                assert(loc != null);

                int edgeIdx = loc.getX();
                int idx = loc.getY();
                
                assert(edgeIdx == edgeIndex.intValue());
                
                assert(idx < edge.getN());
            }
        }
    }

    public static void assertConsistentJoinPointStructures(
        List<PairIntArray> edges,
        Map<Integer, Set<PointPairInt>> edgeFirstEndPointMap, 
        Map<Integer, Set<PointPairInt>> edgeLastEndPointMap, 
        Set<PointPairInt> theJoinPoints, boolean skipForSize3) {
        
        Set<PairInt> joinPointsSet = new HashSet<PairInt>();
        for (PointPairInt entry : theJoinPoints) {
            PairInt loc0 = entry.getKey();
            PairInt loc1 = entry.getValue();
            assert(!loc0.equals(loc1));
            
            PairIntArray edge0 = edges.get(loc0.getX());
            PairIntArray edge1 = edges.get(loc1.getX());
            int n0 = edge0.getN();
            int n1 = edge1.getN();
            
            assert(loc0.getY() < n0);
            assert(loc1.getY() < n1);
            
            if (!(skipForSize3 && (n0 == 3))) {
                assert(!joinPointsSet.contains(loc0));
            }
            if (!(skipForSize3 && (n1 == 3))) {
                assert(!joinPointsSet.contains(loc1));
            }
            
            joinPointsSet.add(loc0);
            joinPointsSet.add(loc1);
        }
        joinPointsSet.clear();
        for (PointPairInt entry : theJoinPoints) {
            PairInt loc0 = entry.getKey();
            PairInt loc1 = entry.getValue();
            assert(!loc0.equals(loc1));
            assert(!joinPointsSet.contains(loc0));
            assert(!joinPointsSet.contains(loc1));
            joinPointsSet.add(loc0);
            joinPointsSet.add(loc1);
        }

        for (Map.Entry<Integer, Set<PointPairInt>> entry : edgeFirstEndPointMap.entrySet()) {
            assert(entry.getValue().size() < 2);
            for (PointPairInt joinPoint : entry.getValue()) {
                PairInt loc0 = joinPoint.getKey();
                PairInt loc1 = joinPoint.getValue();
                assert(loc0 != null);
                assert(loc1 != null);
                                
                boolean contains = theJoinPoints.contains(joinPoint);                   
                if (!contains) {
                    int hash = joinPoint.hashCode();
                    for (PointPairInt entry2 : theJoinPoints) {
                        //log.info("   hash=" + entry2.hashCode());
                        if (hash == entry2.hashCode()) {
                            contains = entry2.equals(joinPoint);
                            //log.info("Set's use of HashMap.getEntry() did not find this point.");
                        }
                    }
                }
                assert(contains);
            }
        }
        for (Map.Entry<Integer, Set<PointPairInt>> entry : edgeLastEndPointMap.entrySet()) {
            assert(entry.getValue().size() < 2);
            for (PointPairInt joinPoint : entry.getValue()) {
                PairInt loc0 = joinPoint.getKey();
                PairInt loc1 = joinPoint.getValue();
                assert(loc0 != null);
                assert(loc1 != null);
                boolean contains = theJoinPoints.contains(joinPoint);                   
                if (!contains) {
                    int hash = joinPoint.hashCode();
                    for (PointPairInt entry2 : theJoinPoints) {
                        //log.info("   hash=" + entry2.hashCode());
                        if (hash == entry2.hashCode()) {
                            contains = entry2.equals(joinPoint);
                            //log.info("Set's use of HashMap.getEntry() did not find this point.");
                        }
                    }
                }
                assert(contains);
            }
        }
    }

    public static void writeJoinPointsImage(Map<PairInt, PairInt> theJoinPointMap, 
        List<PairIntArray> edges, GreyscaleImage img) {
        
        try {
            
            Image img2 = img.copyImageToGreen();

            ImageIOHelper.addAlternatingColorCurvesToImage(edges, img2);
            int nExtraForDot = 1;
            int rClr = 255;
            int gClr = 0;
            int bClr = 100;
            for (Map.Entry<PairInt, PairInt> entry : theJoinPointMap.entrySet()) {
                PairInt loc0 = entry.getKey();
                PairInt loc1 = entry.getValue();
                assert(!loc0.equals(loc1));

                PairIntArray edge0 = edges.get(loc0.getX());
                PairIntArray edge1 = edges.get(loc1.getX());
            
                ImageIOHelper.addPointToImage(
                    edge0.getX(loc0.getY()), edge0.getY(loc0.getY()), 
                    img2, nExtraForDot,
                    rClr, gClr, bClr);
                
                ImageIOHelper.addPointToImage(
                    edge1.getX(loc1.getY()), edge1.getY(loc1.getY()), 
                    img2, nExtraForDot,
                    rClr, gClr, bClr);
            }

            String dirPath = algorithms.util.ResourceFinder.findDirectory("bin");
            String sep = System.getProperty("file.separator");
            ImageIOHelper.writeOutputImage(dirPath + sep + "joinpoints.png", img2);
            
        } catch (IOException e) {
            
        }
    }

    public static void writeImageCopy(ImageExt img, String outfileName) {
        ImageExt img2 = (ImageExt)img.copyImage();
        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/" + outfileName, img2);
        } catch (IOException ex) {
            System.err.println(ex.getMessage());
        }
    }
    
    public static void writeImageCopy(GreyscaleImage img, String outfileName) {
        GreyscaleImage img2 = img.copyImage();
        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/" + outfileName, img2);
        } catch (IOException ex) {
            System.err.println(ex.getMessage());
        }
    }
    
    public static int findEdgeContainingPoint(List<PairIntArray> edges, 
        int startX, int stopX, int startY, int stopY) {
        
        for (int i = 0; i < edges.size(); ++i) {
            
            PairIntArray edge = edges.get(i);
            
            for (int eIdx = 0; eIdx < edge.getN(); ++eIdx) {
                int x = edge.getX(eIdx);
                int y = edge.getY(eIdx);
                if ((x >= startX) && (x <= stopX) && (y >= startY) && (y <= stopY)) {
                    return i;
                }
            }
        }
        return -1;
    }
    
    public static String getInvokingMethodName() {
        
        StackTraceElement[] ste = Thread.currentThread().getStackTrace();
        
        if (ste == null || ste.length == 0) {
            // should not happen. if the thread is not started can return null,
            // but would have had to start already to reach here.
            return "";
        }
                
        return ste[2].getMethodName();
    }
    
    public static void debugPrint(GreyscaleImage input, int xStart, int xStop,
        int yStart, int yStop) {
        
        StringBuilder sb = new StringBuilder();
                    
        for (int row = yStart; row <= yStop; row++) {
            sb.append(String.format("%3d: ", row));
            for (int col = xStart; col <= xStop; col++) {
                sb.append(String.format(" %3d ", input.getValue(col, row)));
            }
            sb.append(String.format("\n"));
        }
        
        System.out.println(sb.toString());
    }
    
    public static void debugPrint(Set<PairInt> points, 
        Set<PairInt> addedPoints, Set<PairInt> removedPoints,
        int xStart, int xStop,
        int yStart, int yStop) {
        
        for (int row = yStop; row >= yStart; row--) {
            StringBuilder sb = new StringBuilder(String.format("row %4d:  ", row));
            for (int col = xStart; col <= xStop; col++) {
                
                PairInt p = new PairInt(col, row);
                
                int v = 0;
                if (!removedPoints.contains(p) 
                    && (addedPoints.contains(p) || points.contains(p))) {
                    v = 1;
                }
                String str = (v == 0) ? String.format("     ") : String.format("%4d ", v);
                sb.append(str);
            }
            System.out.println(sb.toString());
        }
        StringBuilder sb = new StringBuilder(String.format("        "));
        for (int col = xStart; col <= xStop; col++) {
            sb.append(String.format("%4d ", col));
        }
        System.out.println(sb.toString());
        System.out.println("\n");
    }

}
