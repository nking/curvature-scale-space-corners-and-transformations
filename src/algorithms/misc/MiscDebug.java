package algorithms.misc;

import algorithms.sort.CountingSort;
import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.imageProcessing.features.CorrespondenceList;
import algorithms.imageProcessing.scaleSpace.CurvatureScaleSpaceContour;
import algorithms.imageProcessing.features.FeatureComparisonStat;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import static algorithms.imageProcessing.ImageIOHelper.getNextColorRGB;
import algorithms.imageProcessing.scaleSpace.ScaleSpaceCurve;
import algorithms.imageProcessing.scaleSpace.ScaleSpaceCurveImage;
import algorithms.util.Errors;
import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PointPairInt;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import gnu.trove.list.TIntList;
import java.awt.Color;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import no.uib.cipr.matrix.DenseMatrix;

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
            
            sb.append(String.format(
                "  (%d,%d) to (%d,%d) in edges %d and %d  at positions=%d out of %d and %d out of %d\n",
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
            
            Image img2 = img.copyToColorGreyscale();

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
        aIndexes = CountingSort.sort(aIndexes);
        
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
            //assert(entry.getValue().size() < 2);
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
            //assert(entry.getValue().size() < 2);
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
                //TODO: follow where data becomes inconsistent:
                //assert(contains);
            }
        }
    }

    public static void writeEdges(List<PairIntArray> edges, GreyscaleImage img,
        String fileName) {
        try {
            
            Image img2 = ImageIOHelper.convertImage(img);

            ImageIOHelper.addAlternatingColorCurvesToImage(edges, img2);
            
            if (!fileName.contains("\\.")) {
                fileName = fileName + ".png";
            }
            String dirPath = algorithms.util.ResourceFinder.findDirectory("bin");
            String sep = System.getProperty("file.separator");
            ImageIOHelper.writeOutputImage(dirPath + sep + fileName, img2);
            
        } catch (IOException e) {
            
        }
    }
    
    public static void writeJoinPointsImage(Map<PairInt, PairInt> theJoinPointMap, 
        List<PairIntArray> edges, GreyscaleImage img) {
        
        try {
            
            Image img2 = img.copyToColorGreyscale();

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

    public static String printJunctionsToString(Map<Integer, Set<Integer>> jMap, 
        List<PairIntArray> edges, GreyscaleImage img) {

        StringBuilder sb = new StringBuilder("junctions:\n");
        
        for (Entry<Integer, Set<Integer>> entry : jMap.entrySet()) {
            int pixIdx = entry.getKey().intValue();
            int col = img.getCol(pixIdx);
            int row = img.getRow(pixIdx);
            sb.append(String.format("%d (%d,%d)\n", pixIdx, col, row));
        }
        
        return sb.toString();
    }
    
    public static String printJunctionsToString(
        Map<Integer, PairInt> jLocationMap, Map<Integer, Set<Integer>> jMap, 
        List<PairIntArray> edges, GreyscaleImage img) {

        StringBuilder sb = new StringBuilder("junctions:\n");
        
        for (Entry<Integer, Set<Integer>> entry : jMap.entrySet()) {
            int pixIdx = entry.getKey().intValue();
            int col = img.getCol(pixIdx);
            int row = img.getRow(pixIdx);
            sb.append(String.format("%d (%d,%d)\n", pixIdx, col, row));
        }
        
        sb.append("junction locations:\n");
        for (Entry<Integer, PairInt> entry : jLocationMap.entrySet()) {
            Integer pixelIndex = entry.getKey();
            PairInt loc = entry.getValue();
            int edgeIdx = loc.getX();
            int indexWithinEdge = loc.getY();
            int edgeN = edges.get(edgeIdx).getN();
            int x = edges.get(edgeIdx).getX(indexWithinEdge);
            int y = edges.get(edgeIdx).getY(indexWithinEdge);
            sb.append(String.format("edge=%d idx=%d (out of %d) pixIdx=%d (%d,%d)\n", 
                edgeIdx, indexWithinEdge, edgeN, pixelIndex.intValue(), x, y));
        }
        
        return sb.toString();
    }

    public static void assertAllRowsPopulated(
        Map<Integer, List<PairInt>> rowColRanges, int[] rowMinMax, 
        int imageMaxColumn) {
        
        for (int row = rowMinMax[0]; row <= rowMinMax[1]; row++) {
            
            List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(row));
            
            boolean empty = (colRanges == null) || colRanges.isEmpty();
           
            if (empty) {
                System.out.println("error: row " + row + 
                    " has no columns.  row min=" +
                    rowMinMax[0] + " row max=" + 
                    rowMinMax[1] + " col max=" +
                    imageMaxColumn);
            }
            
            assert(!empty);
        }        
    }

    public static int[] findEdgeLocationOfPoint(List<PairIntArray> edges, 
        int x, int y) {
        
        for (int eIdx = 0; eIdx < edges.size(); ++eIdx) {
            
            PairIntArray edge = edges.get(eIdx);
            
            for (int idx = 0; idx < edge.getN(); ++idx) {
                if (edge.getX(idx) == x && edge.getY(idx) == y) {
                    return new int[]{eIdx, idx};
                }
            }
            
        }
        
        return null;
    }

    public static void plotPoints(GreyscaleImage image, DenseMatrix points, 
        int nExtraForDot, String fileNameSuffix) {
        
        Image img = image.copyToColorGreyscale();
        int w = img.getWidth();
        int h = img.getHeight();
        
        int n = points.numColumns();
        
        for (int i = 0; i < n; ++i) {
            int x = (int)Math.round(points.get(0, i));
            int y = (int)Math.round(points.get(1, i));
            for (int dx = -1 * nExtraForDot; dx <= nExtraForDot; ++dx) {
                for (int dy = -1 * nExtraForDot; dy <= nExtraForDot; ++dy) {
                    int x1 = x + dx;
                    int y1 = y + dy;
                    if ((x1 > -1) && (x1 < w) && (y1 > -1) && (y1 < h)) {
                        img.setRGB(x1, y1, 255, 0, 0);
                    }
                }
            }
        }
        
        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/img" + fileNameSuffix 
                + ".png", img);
        } catch (Exception e) {
             e.printStackTrace();
            log.severe("ERROR: " + e.getMessage());
        }
    }
    
    public static void plotPoints(GreyscaleImage image, PairFloatArray points, 
        int nExtraForDot, String fileNameSuffix) {
        
        Image img = image.copyToColorGreyscale();
        int w = img.getWidth();
        int h = img.getHeight();
        
        int n = points.getN();
        
        for (int i = 0; i < n; ++i) {
            int x = Math.round(points.getX(i));
            int y = Math.round(points.getY(i));
            for (int dx = -1 * nExtraForDot; dx <= nExtraForDot; ++dx) {
                for (int dy = -1 * nExtraForDot; dy <= nExtraForDot; ++dy) {
                    int x1 = x + dx;
                    int y1 = y + dy;
                    if ((x1 > -1) && (x1 < w) && (y1 > -1) && (y1 < h)) {
                        img.setRGB(x1, y1, 255, 0, 0);
                    }
                }
            }
        }
        
        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/img" + fileNameSuffix 
                + ".png", img);
        } catch (Exception e) {
             e.printStackTrace();
            log.severe("ERROR: " + e.getMessage());
        }
    }
    
    public static int getCurrentTimeFormatted() {
        double t0 = System.currentTimeMillis();
        double t = t0 - ((int)(t0/1.E9)) * 1E9;
        return (int)t;
    }

    public static void writeHullImages(GreyscaleImage imgGrey, 
        Map<Integer, List<GrahamScan>> hulls, String fileNameSuffix) {
        
        Image imgW = ImageIOHelper.convertImage(imgGrey);
        
        int c = 0;
        for (Entry<Integer, List<GrahamScan>> entry : hulls.entrySet()) {
            List<GrahamScan> hullList = entry.getValue();
            for (GrahamScan hull : hullList) {
                int[] x = new int[hull.getXHull().length];
                int[] y = new int[x.length];
                for (int i = 0; i < x.length; ++i) {
                    x[i] = Math.round(hull.getXHull()[i]);
                    y[i] = Math.round(hull.getYHull()[i]);                   
                }
                if (c == 0) {
                    ImageIOHelper.drawLinesInImage(x, y, imgW, 1, 255, 0, 0);
                } else if (c == 1) {
                    ImageIOHelper.drawLinesInImage(x, y, imgW, 1, 0, 255, 0);
                } else {
                    ImageIOHelper.drawLinesInImage(x, y, imgW, 1, 0, 0, 255);
                }
            }
            c++;
        }
       
        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/img" + fileNameSuffix 
                + ".png", imgW);
           
        } catch (Exception e) {
             e.printStackTrace();
            log.severe(e.getMessage());
        }
    }

    public static void writeImage(int[][] img, String fileName) {
        
        GreyscaleImage imgL = new GreyscaleImage(img.length, img[0].length);
        for (int i = 0; i < img.length; ++i) {
            for (int j = 0; j < img[0].length; ++j) {
                int v = img[i][j];
                imgL.setValue(i, j, v);
            }
        }
        try {
            String bin = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(bin + "/img" + fileName 
                + ".png", imgL);
        } catch (Exception e) {
             e.printStackTrace();
            log.severe("ERROR: " + e.getMessage());
        }
    }
    
    public static void writeImage(float[][] img, String fileName) {
        
        GreyscaleImage imgL = new GreyscaleImage(img.length, img[0].length);
        for (int i = 0; i < img.length; ++i) {
            for (int j = 0; j < img[0].length; ++j) {
                int v = Math.round(img[i][j]);
                imgL.setValue(i, j, v);
            }
        }
        try {
            String bin = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(bin + "/img" + fileName 
                + ".png", imgL);
        } catch (Exception e) {
             e.printStackTrace();
            log.severe("ERROR: " + e.getMessage());
        }
    }
    
    public static void writeImage(double[][] img, String fileName) {
        
        double[][] a = new double[img.length][];
        for (int i = 0; i < img.length; ++i) {
            a[i] = Arrays.copyOf(img[i], img[i].length);
        }
        MiscMath.applyRescale(a, 0, 255);
        
        GreyscaleImage imgL = new GreyscaleImage(img.length, img[0].length);
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[0].length; ++j) {
                int v = (int)Math.round(a[i][j]);
                imgL.setValue(i, j, v);
            }
        }
        try {
            String bin = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(bin + "/img" + fileName 
                + ".png", imgL);
        } catch (Exception e) {
             e.printStackTrace();
            log.severe("ERROR: " + e.getMessage());
        }
    }
        
    public static void writeGreyscaleImage(Image img, String fileNameSuffix) {
        
        try {
            
            String dirPath = ResourceFinder.findDirectory("bin");
            
            String filePath = dirPath + "/img" + fileNameSuffix 
                + ".png";
            
            ImageIOHelper.writeOutputGreyscaleImage(filePath, img);
                        
        } catch (Exception e) {
             e.printStackTrace();
            log.severe("ERROR: " + e.getMessage());
        }
    }
        
    public static void writeImage(Image img, String fileNameSuffix) {
        
        try {
            
            String dirPath = ResourceFinder.findDirectory("bin");
            
            String filePath = dirPath + "/img" + fileNameSuffix 
                + ".png";
            
            ImageIOHelper.writeOutputImage(filePath, img);
                        
        } catch (Exception e) {
             e.printStackTrace();
            log.severe("ERROR: " + e.getMessage());
        }
    }
    
    public static String writeImage(GreyscaleImage img, String fileNameSuffix) {
        
        String str = "";
        
        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/img" + fileNameSuffix 
                + ".png", img);
            str = "wrote to fle " +  dirPath + "/img" + fileNameSuffix 
                + ".png";
        } catch (Exception e) {
             e.printStackTrace();
            log.severe("ERROR: " + e.getMessage());
        }
        
        return str;
    }

    public static void plotCorners(GreyscaleImage imgGrey, PairIntArray corners, 
        String fileNameSuffix) {
        
        Image imgW = ImageIOHelper.convertImage(imgGrey);
        ImageIOHelper.addCurveToImage(corners, imgW, 1, 255, 0, 0);
        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/img" + fileNameSuffix 
                + ".png", imgW);
        } catch (Exception e) {
             e.printStackTrace();
            log.severe("ERROR: " + e.getMessage());
        }
    }
    
    public static void plotCorners(GreyscaleImage imgGrey, 
        Collection<PairInt> points, String fileNameSuffix, int pointSize) {
        
        Image imgW = ImageIOHelper.convertImage(imgGrey);
        ImageIOHelper.addCurveToImage(points, imgW, pointSize, 255, 0, 0);
        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/img" + fileNameSuffix 
                + ".png", imgW);
        } catch (Exception e) {
             e.printStackTrace();
            log.severe("ERROR: " + e.getMessage());
        }
    }

    public static void writeImage(Set<PairInt> points, Image img, 
        String fileSuffix) throws IOException {
       
        int rClr = 255;
        int gClr = 0;
        int bClr = 0;
        
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            int pointSize = 1;//+ Math.round((k - 0.1f)/0.1f);
            if (pointSize < 0) {
                pointSize = 1;
            }
            
            for (int dx = (-1*pointSize); dx < (pointSize + 1); dx++) {
                float xx = x + dx;
                if ((xx > -1) && (xx < (img.getWidth() - 1))) {
                    for (int dy = (-1*pointSize); dy < (pointSize + 1); dy++) {
                        float yy = y + dy;
                        if ((yy > -1) && (yy < (img.getHeight() - 1))) {
                            img.setRGB((int)xx, (int)yy, rClr, gClr, bClr);
                        }
                    }
                }
            }
        }
        
        String dirPath = algorithms.util.ResourceFinder.findDirectory("bin");
        String sep = System.getProperty("file.separator");
        ImageIOHelper.writeOutputImage(dirPath + sep + fileSuffix + ".png", img);
    }

    public static void debugPlot(List<CurvatureScaleSpaceContour> result, ImageExt 
        img, int xOffset, int yOffset, String fileSuffix) {
        
        if (result.isEmpty()) {
            return;
        }
        
        int nExtraForDot = 1;
        int rClr = 255;
        int gClr = 0;
        int bClr = 0;
        
        ImageIOHelper.addContoursToImage(result, img, xOffset, yOffset, 
            nExtraForDot, rClr, gClr, bClr);
        
        try {
            
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.writeOutputImage(dirPath + "/contours_" 
                + fileSuffix + ".png", img);
        
        } catch (IOException e) {}
    }

    public static void printScaleSpaceCurve(ScaleSpaceCurveImage scaleSpaceImage,
        String fileSuffix) throws IOException {
               
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        
        plotScaleSpaceCurve(plotter, scaleSpaceImage, 
            "t vs. sigma for inflection points");
        
        plotter.writeFile(fileSuffix);
    }
    
    public static void plotScaleSpaceCurve(PolygonAndPointPlotter plotter, 
        ScaleSpaceCurveImage scaleSpaceImage, String label) throws IOException {
        
        float[] sigmas = scaleSpaceImage.getImageSigmas();
        
        float[][] tVsSigma = scaleSpaceImage.getScaleSpaceImage();
        
        int n = 0;
        for (int col = 0; col < tVsSigma.length; col++) {
            n += tVsSigma[col].length;
        }
        
        float[] x = new float[n];
        float[] y = new float[n];
        n = 0;
        for (int col = 0; col < tVsSigma.length; col++) {
            float sigma = sigmas[col];
            for (int row = 0; row < tVsSigma[col].length; row++) {
                y[n] = sigma;
                x[n] = tVsSigma[col][row];
                n++;
            }
        }
        
        float xMin = 0;
        float xMax = 1.f;
        float yMin = 0;
        float yMax = 1.1f * algorithms.misc.MiscMath.findMax(y);
                    
        plotter.addPlot(xMin, xMax, yMin, yMax,
            x, y, null, null, label);        
    }
    
    public static void writeImage(PairIntArray[] edges,
        Image img, String suffix) throws IOException {
        
        for (int i = 0; i < edges.length; i++) {
            PairIntArray edge = edges[i];
            ImageIOHelper.addCurveToImage(edge, img, 0, 255, 0, 0);
        }     
        
        String dirPath = ResourceFinder.findDirectory("bin");
        
        ImageIOHelper.writeOutputImage(
            dirPath + "/" + suffix + ".png", img);
    }
    
    public static void writeImage(PairIntArray edge,
        Image img, String suffix) throws IOException {
        
        ImageIOHelper.addCurveToImage(edge, img, 1, 255, 0, 0);
        
        String dirPath = ResourceFinder.findDirectory("bin");
        
        ImageIOHelper.writeOutputImage(
            dirPath + "/" + suffix + ".png", img);
    }
    
    public static void writeImage(PairIntArray edge,
        Image img, int nExtraDot, String suffix) {
        try {
            ImageIOHelper.addCurveToImage(edge, img, nExtraDot, 255, 0, 0);
            
            String dirPath = ResourceFinder.findDirectory("bin");
            
            ImageIOHelper.writeOutputImage(
                dirPath + "/" + suffix + ".png", img);
        } catch (IOException ex) {
            Logger.getLogger(MiscDebug.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public static void print(CorrespondenceList cl) {
        
        List<PairInt> m1 = cl.getPoints1();
        List<PairInt> m2 = cl.getPoints2();
        
        StringBuilder sb = new StringBuilder();
        sb.append("int[] x1 = new int[]{");
        for (PairInt p : m1) {
            sb.append(Integer.valueOf(p.getX())).append(",");
        }
        sb.append("};\n");
        
        sb.append("int[] y1 = new int[]{");
        for (PairInt p : m1) {
            sb.append(Integer.valueOf(p.getY())).append(",");
        }
        sb.append("};\n");
        
        sb.append("int[] x2 = new int[]{");
        for (PairInt p : m2) {
            sb.append(Integer.valueOf(p.getX())).append(",");
        }
        sb.append("};\n");
        
        sb.append("int[] y2 = new int[]{");
        for (PairInt p : m2) {
            sb.append(Integer.valueOf(p.getY())).append(",");
        }
        sb.append("};\n");
        
        log.info(sb.toString());
    }

    public static void writeImages(GreyscaleImage img1, GreyscaleImage img2, 
        List<FeatureComparisonStat> stats, String imgSuffix) {
                
        Image img1Cp = img1.copyToColorGreyscaleExt();
        
        Image img2Cp = img2.copyToColorGreyscaleExt();
        
        writeImages(img1Cp, img2Cp, stats, imgSuffix);
    }
    
    public static void writeImages(GreyscaleImage img1, GreyscaleImage img2, 
        List<FeatureComparisonStat> stats, String imgSuffix,
        int nExtraForDot) {
                
        Image img1Cp = img1.copyToColorGreyscaleExt();
        
        Image img2Cp = img2.copyToColorGreyscaleExt();
        
        writeImages(img1Cp, img2Cp, stats, imgSuffix, nExtraForDot);
    }
    
    public static void writeImages(Image img1, Image img2, 
        List<FeatureComparisonStat> stats, String imgSuffix) {
        
        int nExtraForDot = 1;
        
        writeImages(img1, img2, stats, imgSuffix, nExtraForDot);
    }
    
    public static void writeImagesInAlternatingColor(ImageExt img1, 
        ImageExt img2, List<FeatureComparisonStat> stats, String imgSuffix,
        int nExtraForDot) {
                
        Image img1Cp = img1.copyImage();
        
        Image img2Cp = img2.copyImage();
        
        int w1 = img1Cp.getWidth();
        int h1 = img1Cp.getHeight();
        
        int w2 = img2Cp.getWidth();
        int h2 = img2Cp.getHeight();
        
        int count = 0;
        
        for (FeatureComparisonStat stat : stats) {
            
            PairInt p1 = stat.getImg1Point();
            
            PairInt p2 = stat.getImg2Point();
            
            int clr = getNextColorRGB(count);
            count++;
            
            for (int dx = -1*nExtraForDot; dx <= nExtraForDot; ++dx) {
                for (int dy = -1*nExtraForDot; dy <= nExtraForDot; ++dy) {
            
                    int x1 = p1.getX() + dx;
                    int y1 = p1.getY() + dy;
                    
                    if ((x1 > -1) && (x1 < w1) && (y1 > -1) && (y1 < h1)) {
                        img1Cp.setRGB(x1, y1, clr);
                    }
            
                    int x2 = p2.getX() + dx;
                    int y2 = p2.getY() + dy;
                    
                    if ((x2 > -1) && (x2 < w2) && (y2 > -1) && (y2 < h2)) {
                        img2Cp.setRGB(x2, y2, clr);
                    }
                }
            }
        }
        
        MiscDebug.writeImage(img1Cp, imgSuffix + "_1");
        MiscDebug.writeImage(img2Cp, imgSuffix + "_2");
    }
    
    public static void writeImagesInAlternatingColor(Image img1, 
        Image img2, PairIntArray pai1, PairIntArray pai2, String imgSuffix,
        int nExtraForDot) {
        
        if (pai1.getN() != pai2.getN()) {
            throw new IllegalArgumentException("pai1 and pai2 should be same size");
        }
                
        Image img1Cp = img1.copyImage();
        
        Image img2Cp = img2.copyImage();
        
        int w1 = img1Cp.getWidth();
        int h1 = img1Cp.getHeight();
        
        int w2 = img2Cp.getWidth();
        int h2 = img2Cp.getHeight();
        
        int count = 0;
        
        for (int i = 0; i < pai1.getN(); ++i) {
            
            int x1 = pai1.getX(i);
            int y1 = pai1.getY(i);
            int x2 = pai2.getX(i);
            int y2 = pai2.getY(i);
            
            int clr = getNextColorRGB(count);
            count++;
            
            for (int dx = -1*nExtraForDot; dx <= nExtraForDot; ++dx) {
                for (int dy = -1*nExtraForDot; dy <= nExtraForDot; ++dy) {
            
                    int x = x1 + dx;
                    int y = y1 + dy;
                    
                    if ((x > -1) && (x < w1) && (y > -1) && (y < h1)) {
                        img1Cp.setRGB(x, y, clr);
                    }
            
                    x = x2 + dx;
                    y = y2 + dy;
                    
                    if ((x > -1) && (x < w2) && (y > -1) && (y < h2)) {
                        img2Cp.setRGB(x, y, clr);
                    }
                }
            }
        }
        
        MiscDebug.writeImage(img1Cp, imgSuffix + "_1");
        MiscDebug.writeImage(img2Cp, imgSuffix + "_2");
    }
    
    public static void writeImages(Image img1, Image img2, 
        List<FeatureComparisonStat> stats, String imgSuffix, 
        int nExtraForDot) {
                
        Image img1Cp = img1.copyImage();
        
        Image img2Cp = img2.copyImage();
        
        int w1 = img1Cp.getWidth();
        int h1 = img1Cp.getHeight();
        
        int w2 = img2Cp.getWidth();
        int h2 = img2Cp.getHeight();
        
        for (FeatureComparisonStat stat : stats) {
            
            PairInt p1 = stat.getImg1Point();
            
            PairInt p2 = stat.getImg2Point();
            
            for (int dx = -1*nExtraForDot; dx <= nExtraForDot; ++dx) {
                for (int dy = -1*nExtraForDot; dy <= nExtraForDot; ++dy) {
            
                    int x1 = p1.getX() + dx;
                    int y1 = p1.getY() + dy;
                    
                    if ((x1 > -1) && (x1 < w1) && (y1 > -1) && (y1 < h1)) {
                        img1Cp.setRGB(x1, y1, 255, 0, 0);
                    }
            
                    int x2 = p2.getX() + dx;
                    int y2 = p2.getY() + dy;
                    
                    if ((x2 > -1) && (x2 < w2) && (y2 > -1) && (y2 < h2)) {
                        img2Cp.setRGB(x2, y2, 255, 0, 0);
                    }
                }
            }
        }
        
        MiscDebug.writeImage(img1Cp, imgSuffix + "_1");
        MiscDebug.writeImage(img2Cp, imgSuffix + "_2");
    }

    public static void plotImages(Map<PairInt, List<PairInt>> matched1Matched2, 
        GreyscaleImage img1, GreyscaleImage img2, int nExtraForDot, String fileNameSuffix) {
        
        Image imgCp1 = img1.copyToColorGreyscale();
        Image imgCp2 = img2.copyToColorGreyscale();
        
        ImageIOHelper.addAlternatingColorPointsToImages(matched1Matched2, imgCp1, 
            imgCp2, nExtraForDot);

        MiscDebug.writeImage(imgCp1, fileNameSuffix + "_1");
        MiscDebug.writeImage(imgCp2, fileNameSuffix + "_2");
    }

    public static void plotImages(List<FeatureComparisonStat> stats, 
        GreyscaleImage img1, GreyscaleImage img2, int nExtraForDot, 
        String fileNameSuffix) {
        
        Image imgCp1 = img1.copyToColorGreyscale();
        Image imgCp2 = img2.copyToColorGreyscale();
        
        ImageIOHelper.addAlternatingColorPointsToImages(stats, imgCp1, 
            imgCp2, nExtraForDot);

        MiscDebug.writeImage(imgCp1, fileNameSuffix + "_1");
        MiscDebug.writeImage(imgCp2, fileNameSuffix + "_2");
    }
    
    public static void debugAssertContiguous(List<Set<PairInt>> segmentedCellList) {
        
        for (int i = 0; i < segmentedCellList.size(); ++i) {
            Set<PairInt> set = segmentedCellList.get(i);
            debugAssertContiguous(set);
        }
    }
    
    public static void debugAssertContiguous(Set<PairInt> set) {
        
        int[] minMaxXY = MiscMath.findMinMaxXY(set);
        Set<Integer> xs = new HashSet<Integer>();
        Set<Integer> ys = new HashSet<Integer>();
        for (int k = minMaxXY[0]; k <= minMaxXY[1]; ++k) {
            xs.add(Integer.valueOf(k));
        }
        for (int k = minMaxXY[2]; k <= minMaxXY[3]; ++k) {
            ys.add(Integer.valueOf(k));
        }
        for (PairInt p : set) {
            Integer x = Integer.valueOf(p.getX());
            Integer y = Integer.valueOf(p.getY());
            xs.remove(x);
            ys.remove(y);
        }
        assert(xs.isEmpty());
        assert(ys.isEmpty());
    }
    
    public static void writeAlternatingColor(Image img, 
        List<Set<PairInt>> setList, String fileNameSuffix) {
                
        int nColors = setList.size();
        log.info("begin debug plot");
        int delta = (int)Math.floor(256.f/(float)nColors);
        if (delta == 0) {
            delta = 1;
            nColors = 256;
        }
        Random sr = null;
        long seed = 1234567;//System.currentTimeMillis();
        try {
            sr = SecureRandom.getInstance("SHA1PRNG");
            sr.setSeed(seed);
        } catch (NoSuchAlgorithmException e) {
            sr = new Random(seed);
        }
       
        Set<String> clrs = new HashSet<String>();
        for (int i = 0; i < setList.size(); ++i) {
            
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
                if (delta == 1) {
                    alreadyChosen = false;
                }
                clrs.add(str);
            }
            Set<PairInt> set = setList.get(i);
            for (PairInt p : set) {
                img.setRGB(p.getX(), p.getY(), rClr, gClr, bClr);
            }
        }
        
        MiscDebug.writeImage(img, fileNameSuffix);
    }
    
    public static GreyscaleImage rescaleToImageWithSwapMajor(double[][] a) {
        return rescaleToImage(a, false);
    }
    
    public static GreyscaleImage rescaleToImage(double[][] a) {
        return rescaleToImage(a, true);
    }
    
    private static GreyscaleImage rescaleToImage(double[][] a, boolean useDefaultAxesOrder) {
        
        double minV = Double.MAX_VALUE;
        double maxV = Double.MIN_VALUE;
        
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                double v = a[i][j];
                if (v < minV) {
                    minV = v;
                }
                if (v > maxV) {
                    maxV = v;
                }
            }
        }
        double range = maxV - minV;
        
        double scale = 255./range;
        
        GreyscaleImage output;
        if (useDefaultAxesOrder) {
            output = new GreyscaleImage(a.length, a[0].length);
        } else {
            output = new GreyscaleImage(a[0].length, a.length);
        }
        
//System.out.println("value 0 is rescaled to value=" + ((int)(-minV*scale))
//+ " minV=" + minV + " scale=" + scale);
        
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                
                double v = (a[i][j] - minV) * scale;
                
                if (useDefaultAxesOrder) {
                    output.setValue(i, j, (int)v);
                } else {
                    output.setValue(j, i, (int)v);
                }
            }
        }
        
        return output;
    }

    /**
     * 
     * @param a
     * @param plot0 indexes from dimension 0
     * in array a[dimension0][dimension1]
     * @param plot1 indexes from dimension 1
     * in array a[dimension0][dimension1]
     * @param lbl
     * @return
     * @throws IOException
     */
    public static String plot(double[][] a, TIntList plot0, 
        TIntList plot1, String lbl) throws IOException {
        
        double min = MiscMath.findMin(a);
        double max = MiscMath.findMax(a);
        
        int n1 = a[0].length;
        int n0 = a.length;
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        //    0, n0, 0, n1);
        
        float[] x = new float[n1];
        for (int ii = 0; ii < n1; ++ii) {
            x[ii] = ii;
        }
        float[] y = new float[n1];
        float[] xPolygon = null;
        float[] yPolygon = null;

        for (int i = 0; i < plot0.size(); ++i) {
            int index = plot0.get(i);
            for (int ii = 0; ii < n1; ++ii) {
                y[ii] = (float) a[index][ii];
            }
            float minY = MiscMath.findMin(y);
            float maxY = MiscMath.findMax(y);
            plotter.addPlot(-1, n1 + 1, minY, maxY, x, y, xPolygon,
                yPolygon, lbl + " index0=" + index);
        }
        
        x = new float[n0];
        for (int ii = 0; ii < n0; ++ii) {
            x[ii] = ii;
        }
        y = new float[n0];
        xPolygon = null;
        yPolygon = null;

        for (int i = 0; i < plot1.size(); ++i) {
            int index = plot1.get(i);
            for (int ii = 0; ii < n0; ++ii) {
                y[ii] = (float) a[ii][index];
            }
            float minY = MiscMath.findMin(y);
            float maxY = MiscMath.findMax(y);
            plotter.addPlot(-1, n0 + 1, minY, maxY, x, y, xPolygon,
                yPolygon, lbl + " index1=" + index);
        }
        
        return plotter.writeFile(lbl);
    }

    /**
     * 
     * @param scaleSpaceMap
     * @param fileSuffix
     * @param sigmaIndexDelta set this to 1 to plot every curve, set it to 2
     * to plot every other curve, etc.
     * 
     * @throws IOException 
     */
    public static void printScaleSpaceMap(Map<Float, ScaleSpaceCurve> scaleSpaceMap, 
        String fileSuffix, int sigmaIndexDelta) throws IOException {
        
        PolygonAndPointPlotter plotterC = new PolygonAndPointPlotter();
        PolygonAndPointPlotter plotterX = new PolygonAndPointPlotter();
        PolygonAndPointPlotter plotterY = new PolygonAndPointPlotter();

        // would like to print in order of key, sigma, so putting all into
        //   a tree that uses sort by key
        TreeMap<Float, ScaleSpaceCurve> sortedMap 
            = new TreeMap<Float, ScaleSpaceCurve>(scaleSpaceMap);
        
        final int delta = sigmaIndexDelta - 1;
        int ds = 0;
        
        // print the contour scale space images
        for (Entry<Float, ScaleSpaceCurve> entry : sortedMap.entrySet()) {

            while (ds < delta) {
                ds++;
                continue;
            }
            ds = 0;
            
            Float sigma = entry.getKey();
            ScaleSpaceCurve scaleSpaceCurveSigma = entry.getValue();

            if (scaleSpaceCurveSigma == null || scaleSpaceCurveSigma.getSize() < 2) {
                continue;
            }
            
            String sLabel = String.format(" sigma=%.2f", sigma.floatValue());

            float[] xPoints = new float[scaleSpaceCurveSigma.getSize()];
            float[] yPoints = new float[xPoints.length];
            for (int ii = 0; ii < xPoints.length; ii++) {
                xPoints[ii] = ii;
            }
            float xMin = 0;
            float xMax = 1.1f * algorithms.misc.MiscMath.findMax(xPoints);
            float yMin, yMax;

            // ============ draw X(t,sigma) =============
            for (int ii = 0; ii < xPoints.length; ++ii) {
                yPoints[ii] = scaleSpaceCurveSigma.getX(ii);
            }
            yMin = 0.9f * algorithms.misc.MiscMath.findMin(yPoints);
            yMax = 1.1f * algorithms.misc.MiscMath.findMax(yPoints);

            plotterX.addPlot(
                0, xMax, yMin, yMax,
                null, null, xPoints, yPoints,
                "t vs. X(t, sigma) " + sLabel);
            
            for (int ii = 0; ii < xPoints.length; ii++) {
                yPoints[ii] = scaleSpaceCurveSigma.getK(ii);
            }
            yMin = algorithms.misc.MiscMath.findMin(yPoints);
            if (yMin < 0) {
                yMin *= 1.1;
            } else {
                yMin *= 0.9;
            }
            yMax = algorithms.misc.MiscMath.findMax(yPoints);
            // ==== k vs t
            plotterC.addPlot(
                xMin, xMax, yMin, yMax,
                null, null, xPoints, yPoints,
                "t vs. curvature " + sLabel);

            // ============ draw Y(t,sigma) =============
            Arrays.fill(yPoints, 0);
            for (int ii = 0; ii < xPoints.length; ii++) {
                yPoints[ii] = scaleSpaceCurveSigma.getY(ii);
            }

            yMin = algorithms.misc.MiscMath.findMin(yPoints);
            if (yMin < 0) {
                yMin *= 1.1;
            } else {
                yMin *= 0.9;
            }
            yMax = 1.1f * algorithms.misc.MiscMath.findMax(yPoints);

            plotterY.addPlot(
                xMin, xMax, yMin, yMax,
                xPoints, yPoints, xPoints, yPoints,
                "t vs. Y(t, sigma) " + sLabel);
        }
        
        String filePathC = plotterC.writeFile("C_" + fileSuffix);
        String filePathX = plotterX.writeFile("X_" + fileSuffix);
        String filePath3 = plotterY.writeFile("Y_" + fileSuffix);
    }

    public static String getPrintRowMajor(float[][] a, String label) {
                
        StringBuilder sb = new StringBuilder(label);
        sb.append("\n");
        for (int i = 0; i < a.length; ++i) {
            sb.append(Arrays.toString(a[i])).append("\n");
        }

        return sb.toString();
    }
    
}
