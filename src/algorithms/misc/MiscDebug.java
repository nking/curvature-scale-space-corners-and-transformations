package algorithms.misc;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.List;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class MiscDebug {
    
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
