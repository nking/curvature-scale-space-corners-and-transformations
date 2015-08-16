package algorithms.misc;

import algorithms.CountingSort;
import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.CornerRegion;
import algorithms.imageProcessing.CurvatureScaleSpaceContour;
import algorithms.imageProcessing.CurvatureScaleSpaceImagePoint;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageDisplayer;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.PixelColors;
import algorithms.imageProcessing.ScaleSpaceCurveImage;
import algorithms.util.Errors;
import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PointPairInt;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import algorithms.util.ScatterPointPlotterPNG;
import java.awt.Color;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
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

    public static void writeImageCopy(Image img, String outfileName) {
        Image img2 = img.copyImage();
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
        int imageMaxColumn, int imageMaxRow) {
        
        for (int row = rowMinMax[0]; row <= rowMinMax[1]; row++) {
            
            List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(row));
            
            boolean empty = (colRanges == null) || colRanges.isEmpty();
            
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

    /**
     * a debug method to print scatter diagrams and histograms if needed
     * of the colors in search of best indicator that blue or red skies
     * have many clouds.
     * @param points
     * @param colorImg
     * @param xOffset
     * @param yOffset
     */
    public static void plotSkyColor(Set<PairInt> points, ImageExt colorImg, 
        int xOffset, int yOffset) {

        int n = points.size();
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        float[] cieX = new float[n];
        float[] cieY = new float[n];
        float[] hue = new float[n];
        float[] saturation = new float[n];
        float[] brightness = new float[n];
        float[] r = new float[n];
        float[] g = new float[n];
        float[] b = new float[n];
        
        int i = 0;
        for (PairInt p : points) {
            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            
            int idx = colorImg.getInternalIndex(x, y);

            int rr = colorImg.getR(idx);
            int gg = colorImg.getG(idx);
            int bb = colorImg.getB(idx);
                        
            r[i] = rr;
            g[i] = gg;
            b[i] = bb;
            
            float[] cieXY = cieC.rgbToXYChromaticity(rr, gg, bb);
            cieX[i] = cieXY[0];
            cieY[i] = cieXY[1];
            
            float[] hsb = new float[3];
            Color.RGBtoHSB(rr, gg, bb, hsb);
            hue[i] = hsb[0];
            saturation[i] = hsb[1];
            brightness[i] = hsb[2];
            
            i++;
        }
        
        int plotNumber = getCurrentTimeFormatted();
        
        try {
            ScatterPointPlotterPNG plotter = new ScatterPointPlotterPNG();
            
            plotter.plot(0.2f, 0.6f, 0.2f, 0.6f, 
                cieX, cieY, "CIE X vs Y", "CIEX", "CIEY");
            plotter.writeFile(Integer.valueOf(plotNumber));
            
            plotter = new ScatterPointPlotterPNG();
            plotter.plot(0.0f, 256f, 0.0f, 1.1f, 
                b, brightness, "B vs brightness", "B", "brightness");
            plotter.writeFile(Integer.valueOf(plotNumber + 1));
            
            plotter = new ScatterPointPlotterPNG();
            plotter.plot(0.0f, 256f, 0.0f, 1.1f, 
                r, brightness, "R vs brightness", "R", "brightness");
            plotter.writeFile(Integer.valueOf(plotNumber + 2));
            
            plotter = new ScatterPointPlotterPNG();
            plotter.plot(0.0f, 256f, 0.0f, 256f, 
                b, r, "B vs R", "B", "R");
            plotter.writeFile(Integer.valueOf(plotNumber + 3));
            
            plotter = new ScatterPointPlotterPNG();
            plotter.plot(0.0f, 1.1f, 0.0f, 1.1f, 
                hue, saturation, "hue vs saturation", "hue", "saturation");
            plotter.writeFile(Integer.valueOf(plotNumber + 4));
            
            plotter = new ScatterPointPlotterPNG();
            plotter.plot(0.0f, 0.8f, 0.0f, 1.1f, 
                cieX, brightness, "CIE X vs brightness", "CIEX", "brightness");
            plotter.writeFile(Integer.valueOf(plotNumber + 5));
            
            plotter = new ScatterPointPlotterPNG();
            plotter.plot(0.0f, 0.9f, 0.0f, 1.1f, 
                cieY, brightness, "CIE Y vs brightness", "CIEY", "brightness");
            plotter.writeFile(Integer.valueOf(plotNumber + 6));
                        
            // for hue histograms, because it's space is formed in 0 to 360 degrees,
            // need to wrap around the bins.
            // for that reason, will look for an empty bin and cut and merge
            // the wrap
            
            HistogramHolder hueHist = Histogram.createSimpleHistogram(
                0.f, 1.0f, 10, hue, Errors.populateYErrorsBySqrt(hue));
                        
            int[] yh = Arrays.copyOf(hueHist.getYHist(), hueHist.getYHist().length);
            int zeroBinIdx = -1;
            for (int ii = 0; ii < yh.length; ii++) {
                if (yh[ii] == 0) {
                    zeroBinIdx = ii;
                    break;
                }
            }
            if (zeroBinIdx > -1) {
                int[] yShifted = new int[yh.length];
                float[] yShiftedF = new float[yh.length];
                int count = 0;
                for (int ii = (zeroBinIdx + 1); ii < yh.length; ii++) {
                    yShifted[count] = yh[ii];
                    yShiftedF[count] = yh[ii];
                    count++;
                }
                for (int ii = 0; ii <= zeroBinIdx; ii++) {
                    yShifted[count] = yh[ii];
                    yShiftedF[count] = yh[ii];
                    count++;
                }
                hueHist.setYHist(yShifted);
                hueHist.setYHistFloat(yShiftedF);
            }
            
            float fwhmHue = Histogram.measureFWHMOfStrongestPeak(hueHist);
            
            hueHist.plotHistogram("hue shifted by " + (zeroBinIdx + 0) + " bins",
                plotNumber);
            
            log.info("fwhm hue=" + fwhmHue);

            HistogramHolder saturationHist = Histogram.createSimpleHistogram(
                0.f, 1.0f, 10, 
                saturation, Errors.populateYErrorsBySqrt(saturation));
            
            saturationHist.plotHistogram("saturation", plotNumber + 1);
            
            float[] fwhmSaturation = Histogram.measureFWHMOfAllPeaks(
                saturationHist, 0.1f);
            
            log.info("fwhm saturation=" + Arrays.toString(fwhmSaturation));
            
        } catch (IOException e) {
            
            log.severe(e.getMessage());
        }
    }
    
    public static void plotPointSets(PairFloatArray set1, PairIntArray set2, 
        int imageWidth1, int imageHeight1, int imageWidth2, int imageHeight2, 
        int fileNumber) {
        
        float[] x1 = new float[set1.getN()];
        float[] y1 = new float[x1.length];
        for (int i = 0; i < set1.getN(); ++i) {
            x1[i] = set1.getX(i);
            y1[i] = set1.getY(i);
        }
        
        float[] x2 = new float[set2.getN()];
        float[] y2 = new float[x2.length];
        for (int i = 0; i < set2.getN(); ++i) {
            x2[i] = set2.getX(i);
            y2[i] = set2.getY(i);
        }
        
        try {
            ScatterPointPlotterPNG plotter = new ScatterPointPlotterPNG();
            
            plotter.plot(0, imageWidth1, 0, imageHeight1, 
                x1, y1, "set1", "X", "Y");
            plotter.writeFile(Integer.valueOf(fileNumber));
            
            plotter.plot(0, imageWidth2, 0, imageHeight2, 
                x2, y2, "set2", "X", "Y");
            plotter.writeFile(Integer.valueOf(fileNumber + 1));
           
        } catch (IOException e) {
            
            log.severe(e.getMessage());
        }
    }
    
    public static void plotPoints(Set<PairInt> set,
        int imageWidth, int imageHeight, int fileNumber) {
        
        float[] x1 = new float[set.size()];
        float[] y1 = new float[x1.length];
        int i = 0;
        for (PairInt p : set) {
            x1[i] = p.getX();
            y1[i] = p.getY();
            i++;
        }
        
        try {
            ScatterPointPlotterPNG plotter = new ScatterPointPlotterPNG();
            
            plotter.plot(0, imageWidth, 0, imageHeight, 
                x1, y1, "points", "X", "Y");
            plotter.writeFile(Integer.valueOf(fileNumber));
            
        } catch (IOException e) {
            
            log.severe(e.getMessage());
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

    public static void writeImage(ImageExt img, String fileNameSuffix) {
        
        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/img" + fileNameSuffix 
                + ".png", img);
        } catch (Exception e) {
             e.printStackTrace();
            log.severe("ERROR: " + e.getMessage());
        }
    }
    
    public static void writeImage(GreyscaleImage img, String fileNameSuffix) {
        
        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/img" + fileNameSuffix 
                + ".png", img);
        } catch (Exception e) {
             e.printStackTrace();
            log.severe("ERROR: " + e.getMessage());
        }
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

    public static void printCornerRegion(CornerRegion[] cr) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < cr.length; ++i) {
            if (cr[i] == null) {
                continue;
            }
            sb.append(Integer.valueOf(i)).append(") ").append(cr[i].toString());
        }
        
        log.info(sb.toString());
    }

    public static void display(CornerRegion cornerRegion1, 
        CornerRegion cornerRegion2, GreyscaleImage img1Grey, 
        GreyscaleImage img2Grey, String index, int halfWidthOfDisplay) 
        throws CornerRegion.CornerRegionDegneracyException, IOException {
        
        int sWidth = 2*halfWidthOfDisplay;
        int sHeight = 2*halfWidthOfDisplay;
        
        int kMaxIdx = cornerRegion1.getKMaxIdx();
        int x1 = cornerRegion1.getX()[kMaxIdx];
        int y1 = cornerRegion1.getY()[kMaxIdx];
        float orientation1 = cornerRegion1.getRelativeOrientationInDegrees();
        String label1 = String.format("[%s] (%d,%d) %.1f degrees", index, x1, y1, 
            orientation1);
        
        GreyscaleImage sImg1 = img1Grey.subImage(x1, y1, sWidth, sHeight);
        
        kMaxIdx = cornerRegion2.getKMaxIdx();
        int x2 = cornerRegion2.getX()[kMaxIdx];
        int y2 = cornerRegion2.getY()[kMaxIdx];
        float orientation2 = cornerRegion2.getRelativeOrientationInDegrees();
        String label2 = String.format("[%s] (%d,%d) %.1f degrees", index, x2, y2, 
            orientation2);
        
        GreyscaleImage sImg2 = img2Grey.subImage(x2, y2, sWidth, sHeight);
        
        ImageDisplayer.displayImage(label1, sImg1);
        ImageDisplayer.displayImage(label2, sImg2);
        
        int z = 1;
    }
    
    public static void writeImage(Set<CornerRegion> cornerRegions, 
        Image img, String fileSuffix) throws IOException {
       
        int rClr = 255;
        int gClr = 0;
        int bClr = 0;
        
        for (CornerRegion cr : cornerRegions) {
            int kMaxIdx = cr.getKMaxIdx();
            int x = cr.getX()[kMaxIdx];
            int y = cr.getY()[kMaxIdx];
            float k = cr.getK()[kMaxIdx];
            int pointSize = 1 + Math.round((k - 0.1f)/0.1f);
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
        
        for (int i = 0; i < result.size(); i++) {
            
            CurvatureScaleSpaceContour cssC = result.get(i);
            
            CurvatureScaleSpaceImagePoint[] peakDetails = cssC.getPeakDetails();
            
            for (CurvatureScaleSpaceImagePoint peakDetail : peakDetails) {
                int x = peakDetail.getXCoord() + xOffset;
                int y = peakDetail.getYCoord() + yOffset;
                for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                    float xx = x + dx;
                    if ((xx > -1) && (xx < (img.getWidth() - 1))) {
                        for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); 
                            dy++) {
                            float yy = y + dy;
                            if ((yy > -1) && (yy < (img.getHeight() - 1))) {
                                img.setRGB((int)xx, (int)yy, rClr, gClr, bClr);
                            }
                        }
                    }
                }
            }            
        }
        
        try {
            
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.writeOutputImage(dirPath + "/contours_" 
                + fileSuffix + ".png", img);
        
        } catch (IOException e) {}
    }

    public static void printScaleSpaceCurve(ScaleSpaceCurveImage scaleSpaceImage,
        int fileNumber) throws IOException {
        
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
            
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        
        plotter.addPlot(xMin, xMax, yMin, yMax,
            x, y, null, null, 
            "t vs. sigma for inflection points");
        
        plotter.writeFile(fileNumber);
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
}
