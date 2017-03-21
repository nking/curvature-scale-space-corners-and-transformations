package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class PostLineThinnerCorrections {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    protected static final int[] fourNeighborsX = new int[]{0,  1,  0,  -1};
    protected static final int[] fourNeighborsY = new int[]{-1, 0,  1,   0};
    
    protected static final int[] eightNeighborsX = 
        new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
    protected static final int[] eightNeighborsY = 
        new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
   
    protected void correctForExtCorner(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*        
           0  0  #      2
        0  0  #  0      1
        0  #*<#  0      0
        0  0  0  0     -1
        
       -1  0  1  2  3
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
   
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(0, 1)); zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 1)); zeroes.add(new PairInt(2, 0)); zeroes.add(new PairInt(2, -1));
        
        ones.add(new PairInt(1, 0)); ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, -2));
        
        changeToZeroes.add(new PairInt(1, 0));
                    
        int nCorrections = 0;
        
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    protected void reverseXs(
        final TIntSet zeroes, final TIntSet ones, TIntSet changeToZeroes, 
        final TIntSet changeToOnes, int w, int h) {

        int[] xy = new int[2];
        
        // ----- change the sign of x to handle other direction -----
        TIntSet tmp;
        TIntSet set;
        TIntIterator iter;
        for (int i = 0; i < 4; ++i) {
            tmp = new TIntHashSet();
            switch (i) {
                case 0:
                    set = zeroes;
                    break;
                case 1:
                    set = ones;
                    break;
                case 2:
                    set = changeToZeroes;
                    break;
                default:
                    set = changeToOnes;
                    break;
            }
            iter = set.iterator();
            while (iter.hasNext()) {
                int ePixIdx = iter.next();
                decodePixelOffsets(ePixIdx, w, xy);
                tmp.add(encodePixelOffsets(-xy[0], xy[1], w));
            }
            set.clear();
            set.addAll(tmp);
        }
    }
    
    protected void reverseYs(
        final TIntSet zeroes, final TIntSet ones, TIntSet changeToZeroes, 
        final TIntSet changeToOnes, int w, int h) {

        int[] xy = new int[2];
        
        // ----- change the sign of y to handle other direction -----
        TIntSet tmp;
        TIntSet set;
        TIntIterator iter;
        for (int i = 0; i < 4; ++i) {
            tmp = new TIntHashSet();
            switch (i) {
                case 0:
                    set = zeroes;
                    break;
                case 1:
                    set = ones;
                    break;
                case 2:
                    set = changeToZeroes;
                    break;
                default:
                    set = changeToOnes;
                    break;
            }
            iter = set.iterator();
            while (iter.hasNext()) {
                int ePixIdx = iter.next();
                decodePixelOffsets(ePixIdx, w, xy);
                tmp.add(encodePixelOffsets(xy[0], -xy[1], w));
            }
            set.clear();
            set.addAll(tmp);
        }
    }
    
    protected void reverseXs(
        final Set<PairInt> zeroes, final Set<PairInt> ones, 
        Set<PairInt> changeToZeroes, final Set<PairInt> changeToOnes) {
        
        // ----- change the sign of x to handle other direction -----
        for (PairInt p : zeroes) {
            p.setX(-1 * p.getX());
        }
        
        for (PairInt p : ones) {
            p.setX(-1 * p.getX());
        }
          
        for (PairInt p : changeToZeroes) {
            p.setX(-1 * p.getX());
        }
        
        for (PairInt p : changeToOnes) {
            p.setX(-1 * p.getX());
        }
    }
    
    protected void reverseXs(final Set<PairInt> zeroes, final Set<PairInt> ones) {
        
        // ----- change the sign of x to handle other direction -----
        for (PairInt p : zeroes) {
            p.setX(-1 * p.getX());
        }
        
        for (PairInt p : ones) {
            p.setX(-1 * p.getX());
        }
        
    }
    
    protected void reverseYs(final Set<PairInt> zeroes, final Set<PairInt> ones, 
        Set<PairInt> changeToZeroes, final Set<PairInt> changeToOnes) {
        
        // ----- change the sign of y  -----
        for (PairInt p : zeroes) {
            p.setY(-1 * p.getY());
        }
        
        for (PairInt p : ones) {
            p.setY(-1 * p.getY());
        }
        
        for (PairInt p : changeToZeroes) {
            p.setY(-1 * p.getY());
        }
        
        for (PairInt p : changeToOnes) {
            p.setY(-1 * p.getY());
        }
        
    }
    
    protected void reverseYs(final Set<PairInt> zeroes, final Set<PairInt> ones) {
        
        // ----- change the sign of y  -----
        for (PairInt p : zeroes) {
            p.setY(-1 * p.getY());
        }
        
        for (PairInt p : ones) {
            p.setY(-1 * p.getY());
        }
        
    }
    
    private int rotate90ThreeTimes(
        TIntSet pixIdxs, int imageWidth, int imageHeight,
        final TIntSet zeroes, final TIntSet ones, TIntSet changeToZeroes, 
        final TIntSet changeToOnes) {
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes, imageWidth, imageHeight);
        
        int nCorrections = 0;
        
        nCorrections += replacePattern(
            pixIdxs, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
             
        // ----- change the sign of y to handle other direction -----
        reverseYs(zeroes, ones, changeToZeroes, changeToOnes, imageWidth, imageHeight);
        
        nCorrections += replacePattern(
            pixIdxs, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of x to handle another direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes, imageWidth, imageHeight);
                   
        nCorrections += replacePattern(
            pixIdxs, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        return nCorrections;
    }
    
    
    // the points which would be changed to zeroes or changed to ones are returned.
    private int rotate90ThreeTimes(
        Set<PairInt> points, int imageWidth, int imageHeight,
        final LinkedHashSet<PairInt> zeroes, final LinkedHashSet<PairInt> ones, 
        LinkedHashSet<PairInt> changeToZeroes, 
        final LinkedHashSet<PairInt> changeToOnes) {
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);
        
        int nCorrections = 0;
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
             
        // ----- change the sign of y to handle other direction -----
        reverseYs(zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        // ----- change the sign of x to handle another direction -----
        reverseXs(zeroes, ones, changeToZeroes, changeToOnes);
                    
        nCorrections += replacePattern(
            points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        return nCorrections;
    }
    
    
    // the points which would be changed to zeroes or changed to ones are returned.
    private Set<PairInt> findPattern(
        Set<PairInt> points, int imageWidth, int imageHeight,
        final LinkedHashSet<PairInt> zeroes, final LinkedHashSet<PairInt> ones, 
        final LinkedHashSet<PairInt> changeToZeroes, 
        final LinkedHashSet<PairInt> changeToOnes) {
        
        int w = imageWidth;
        int h = imageHeight;
        
        Set<PairInt> output = new HashSet<PairInt>();

        Set<PairInt> tmpPointsRemoved = new HashSet<PairInt>();
        Set<PairInt> tmpPointsAdded = new HashSet<PairInt>();
        
        for (PairInt p : points) {

            boolean foundPattern = foundPattern(p, points, imageWidth,
                imageHeight, zeroes, ones, tmpPointsRemoved, tmpPointsAdded);
            
            if (foundPattern) {
                for (PairInt p2 : changeToZeroes) {
                    PairInt p3 = new PairInt(p.getX() + p2.getX(),
                        p.getY() + p2.getY());
                    output.add(p3);
                }
                for (PairInt p2 : changeToOnes) {
                    PairInt p3 = new PairInt(p.getX() + p2.getX(),
                        p.getY() + p2.getY());
                    output.add(p3);
                }
            }
            
            tmpPointsAdded.clear();              
            tmpPointsRemoved.clear();
        }
        
        return output;
    }
    
    private int replacePattern(
        Set<PairInt> points, int imageWidth, int imageHeight,
        final LinkedHashSet<PairInt> zeroes, final LinkedHashSet<PairInt> ones, 
        final LinkedHashSet<PairInt> changeToZeroes, 
        final LinkedHashSet<PairInt> changeToOnes) {
        
        int w = imageWidth;
        int h = imageHeight;

        Set<PairInt> tmpPointsRemoved = new HashSet<PairInt>();
        Set<PairInt> tmpPointsAdded = new HashSet<PairInt>();
        
        for (PairInt p : points) {

            boolean foundPattern = foundPattern(p, points, imageWidth,
                imageHeight, zeroes, ones, tmpPointsRemoved, tmpPointsAdded);
            
            if (!foundPattern) {
                continue;
            }
            
            int col = p.getX();
            int row = p.getY();
            
            for (PairInt p2 : changeToZeroes) {
                int x = col + p2.getX();
                int y = row + p2.getY();
                if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                    continue;
                }
                PairInt p3 = new PairInt(x, y);
                tmpPointsAdded.remove(p3);              
                tmpPointsRemoved.add(p3);
            }

            for (PairInt p2 : changeToOnes) {
                int x = col + p2.getX();
                int y = row + p2.getY();
                if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                    continue;
                }
                PairInt p3 = new PairInt(x, y);                 
                tmpPointsRemoved.remove(p3);
                tmpPointsAdded.add(p3);
            }
        }
        
        int nCorrections = tmpPointsRemoved.size() + tmpPointsAdded.size();
        
        for (PairInt p2 : tmpPointsRemoved) {
            points.remove(p2);
        }
        for (PairInt p2 : tmpPointsAdded) {
            points.add(p2);
        }
        
        return nCorrections;
    }
    
    private int replacePattern(
        TIntSet pixIdxs, int imageWidth, int imageHeight,
        final TIntSet zeroes, final TIntSet ones, 
        final TIntSet changeToZeroes, final TIntSet changeToOnes) {
        
        int w = imageWidth;
        int h = imageHeight;

        TIntSet tmpPointsRemoved = new TIntHashSet();
        TIntSet tmpPointsAdded = new TIntHashSet();
    
        int[] xyOff = new int[2];
        
        TIntIterator iter = pixIdxs.iterator();
        while (iter.hasNext()) {

            int pixIdx = iter.next();
        
            boolean foundPattern = foundPattern(pixIdx, pixIdxs, imageWidth,
                imageHeight, zeroes, ones, tmpPointsRemoved, tmpPointsAdded);
            
            if (!foundPattern) {
                continue;
            }
                        
            int row = pixIdx/imageWidth;
            int col = pixIdx - (row * imageWidth);
            
            TIntIterator iter2 = changeToZeroes.iterator();
            while (iter2.hasNext()) {
                int ePixIdx2 = iter2.next();
                decodePixelOffsets(ePixIdx2, imageWidth, xyOff);
                int x2 = col + xyOff[0];
                int y2 = row + xyOff[1];
                if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                    continue;
                }
                int pixIdx3 = (y2 * w) + x2; 
                tmpPointsAdded.remove(pixIdx3);              
                tmpPointsRemoved.add(pixIdx3);
            }

            iter2 = changeToOnes.iterator();
            while (iter2.hasNext()) {
                int ePixIdx2 = iter2.next();
                decodePixelOffsets(ePixIdx2, imageWidth, xyOff);
                int x2 = col + xyOff[0];
                int y2 = row + xyOff[1];
                if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                    continue;
                }
                int pixIdx3 = (y2 * w) + x2;
                tmpPointsRemoved.remove(pixIdx3);
                tmpPointsAdded.add(pixIdx3);
            }
        }
        
        int nCorrections = tmpPointsRemoved.size() + tmpPointsAdded.size();
        
        pixIdxs.removeAll(tmpPointsRemoved);
        pixIdxs.addAll(tmpPointsAdded);
        
        return nCorrections;
    }
   
    /**
     * remove single isolated pixels
     * 
     * @param points
     */
    public static void correctForIsolatedPixels(Set<PairInt> points) {
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        Set<PairInt> rm = new HashSet<PairInt>();
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            boolean found = false;
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                PairInt p2 = new PairInt(x2, y2);
                if (points.contains(p2) && !rm.contains(p2)) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                // no neighbors
                rm.add(p);
            }
        }
        for (PairInt p : rm) {
            points.remove(p);
        }
    }
    
    /**
     * for gaps of one in lines, collect the gaps and return them.
     * @param points
     * @param imageWidth
     * @param imageHeight
     * @return 
     */
    public Set<PairInt> findGapsOf1(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        int w = imageWidth;
        int h = imageHeight;

        /*
        0  1  2
        7     3
        6  5  4
        fill in !value if these pairs are filled in:
            0:3, 0:4, 0:5
            1:4, 1:5, 1:6
            2:5, 2:6, 2:7
            3:6, 3:7, 3:0
            4:7
        so a +1 and -1 in x or y and a +1 or -1 in y or x
        */
        int[] dxs0 = new int[]{-1, -1, -1,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1};
        int[] dys0 = new int[]{+1, +1, +1,  1,  1,  1,  1,  1,  1,  0,  0,  0, -1};
        int[] dxs1 = new int[]{1,  +1,  0,  1,  0, -1,  0, -1, -1, -1, -1, -1, -1};
        int[] dys1 = new int[]{0,  -1, -1, -1, -1, -1, -1, -1,  0, -1,  0,  1,  0};
        
        Set<PairInt> gaps = new HashSet<PairInt>();
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {

                PairInt p = new PairInt(i, j);
                if (points.contains(p)) {
                    continue;
                }

                for (int k = 0; k < dxs0.length; ++k) {
                    int x1 = i + dxs0[k];
                    int y1 = j + dys0[k];
                    if (x1 < 0 || (x1 > (w - 1)) || y1 < 0 || (y1 > (h - 1))) {
                        continue;
                    }
                    PairInt p1 = new PairInt(x1, y1);
                    if (!points.contains(p1)) {
                        continue;
                    }
                    int x2 = i + dxs1[k];
                    int y2 = j + dys1[k];
                    if (x2 < 0 || (x2 > (w - 1)) || y2 < 0 || (y2 > (h - 1))) {
                        continue;
                    }
                    PairInt p2 = new PairInt(x2, y2);
                    if (!points.contains(p2)) {
                        continue;
                    }
                    gaps.add(p);
                    break;
                }
            }
        }
        
        return gaps;
    }

    private boolean foundPattern(PairInt p, Set<PairInt> points, 
        int imageWidth, int imageHeight, LinkedHashSet<PairInt> zeroes, 
        LinkedHashSet<PairInt> ones, Set<PairInt> tmpPointsRemoved, 
        Set<PairInt> tmpPointsAdded) {
        
        boolean isRemoved = tmpPointsRemoved.contains(p);
        if (isRemoved) {
            return false;
        }
            
        int col = p.getX();
        int row = p.getY();     
        
        int w = imageWidth;
        int h = imageHeight;

        boolean foundPattern = true;

        for (PairInt p2 : zeroes) {
            int x = col + p2.getX();
            int y = row + p2.getY();
            if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                return false;
            }
            PairInt p3 = new PairInt(x, y);
            if (!tmpPointsRemoved.contains(p3)
                && (points.contains(p3) || tmpPointsAdded.contains(p3))) {
                return false;
            }
        }

        for (PairInt p2 : ones) {
            int x = col + p2.getX();
            int y = row + p2.getY();
            if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                return false;
            }
            PairInt p3 = new PairInt(x, y);
            if (tmpPointsRemoved.contains(p3) ||
                (!points.contains(p3) && !tmpPointsAdded.contains(p3))) {
                return false;
            }
        }
        return true;
    }
    
    private boolean foundPattern(int pixIdx, TIntSet pixIdxs, 
        int imageWidth, int imageHeight, TIntSet zeroes, TIntSet ones, 
        TIntSet tmpPointsRemoved, TIntSet tmpPointsAdded) {
        
        boolean isRemoved = tmpPointsRemoved.contains(pixIdx);
        if (isRemoved) {
            return false;
        }
        
        int w = imageWidth;
        int h = imageHeight;
                    
        int row = pixIdx/w;
        int col = pixIdx - (row * w);
        
        int[] xyOff = new int[2];
        
        TIntIterator iter = zeroes.iterator();
        while (iter.hasNext()) {
            int ePixIdx2 = iter.next();
            decodePixelOffsets(ePixIdx2, w, xyOff);
            int x2 = col + xyOff[0];
            int y2 = row + xyOff[1];
            if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                return false;
            }
            int pixIdx3 = (y2 * w) + x2;
            if (!tmpPointsRemoved.contains(pixIdx3)
                && (pixIdxs.contains(pixIdx3) || tmpPointsAdded.contains(pixIdx3))) {
                return false;
            }
        }

        iter = ones.iterator();
        while (iter.hasNext()) {
            int ePixIdx2 = iter.next();
            decodePixelOffsets(ePixIdx2, w, xyOff);
            int x2 = col + xyOff[0];
            int y2 = row + xyOff[1];
            if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                return false;
            }
            int pixIdx3 = (y2 * w) + x2;
            if (tmpPointsRemoved.contains(pixIdx3) ||
                (!pixIdxs.contains(pixIdx3) && !tmpPointsAdded.contains(pixIdx3))) {
                return false;
            }
        }
        
        return true;
    }
   
    /**
     * use with caution
     * @param img
     */
    public void extremeStaircaseRemover(GreyscaleImage img) {
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        for (int i = 0; i < img.getNPixels(); ++i) {
            if (img.getValue(i) > 0) {
                points.add(new PairInt(img.getCol(i), img.getRow(i)));
            }
        }
        
        Set<PairInt> cp = new HashSet<PairInt>(points);
        
        extremeStaircaseRemover(points, img.getWidth(), img.getHeight());
        
        extremeCornerRemover(points, img.getWidth(), img.getHeight());
        
        cp.removeAll(points);
        
        for (PairInt p : cp) {
            img.setValue(p.getX(), p.getY(), 0);
        }
        
    }

    /**
     * use with caution
     * @param img
     */
    public void extremeCornerRemover(GreyscaleImage img) {
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        for (int i = 0; i < img.getNPixels(); ++i) {
            if (img.getValue(i) > 0) {
                points.add(new PairInt(img.getCol(i), img.getRow(i)));
            }
        }
        
        Set<PairInt> cp = new HashSet<PairInt>(points);
                
        extremeCornerRemover(points, img.getWidth(), img.getHeight());
        
        cp.removeAll(points);
        
        for (PairInt p : cp) {
            img.setValue(p.getX(), p.getY(), 0);
        }
        
    }
    
    /**
     * use with caution
     * @param points
     * @param imageWidth
     * @param imageHeight 
     */
    public void extremeCornerRemover(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
            0  0            2
            0  #< #         1
               #* 0         0
                           -1
                           -2
           -1  0  1  2  3        
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(1, 0)); zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, -2)); zeroes.add(new PairInt(-1, -2));
        
        ones.add(new PairInt(0, -1)); 
        ones.add(new PairInt(1, -1));
        
        changeToZeroes.add(new PairInt(0, -1));
                    
        int nCorrections = 0;
       
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    /**
     * NOTE: only use on closed curves because it has
     * a spur remover that will remove an entire line if
     * it's not closed.
     * @param input 
     */
    public void extremeThinning(GreyscaleImage input) {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        final int w = input.getWidth();
        final int h = input.getHeight();
        
        Set<PairInt> points = imageProcessor.readNonZeroPixels(input);
    
        Set<PairInt> cp = new HashSet<PairInt>(points);
        
        extremeThinning(points, w, h);
    
        cp.removeAll(points);
        
        for (PairInt p : cp) {
            input.setValue(p.getX(), p.getY(), 0);
        }
    }
    
    /**
     * NOTE: only use on closed curves because it has
     * a spur remover that will remove an entire line if
     * it's not closed.
     * @param input 
     */
    public void extremeThinning(Set<PairInt> points, int w, int h) {
            
        extremeStaircaseRemover(points, w, h);
        extremeCornerRemover(points, w, h);
        SpurRemover spurRemover = new SpurRemover();
        spurRemover.remove(points, w, h);                
    }
    
    /**
     * use with caution.  it has a small window it uses to determine which
     * pixel to remove.  
     * if precision is needed, consider the zig-zag methods because they
     * use larger windows and patterns.
     * @param points
     * @param imageWidth
     * @param imageHeight 
     */
    public void extremeStaircaseRemover(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
            0  0  0         2
            0  #< #         1
            #  #*           0
                           -1
           -1  0  1  2  3        
        */
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        /*       
            0  0  0         2
            0  #< #         1
            #  #*           0
                           -1
           -1  0  1  2  3        
        */
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, -2));
        
        ones.add(new PairInt(-1, 0)); 
        ones.add(new PairInt(0, 0));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, -1));
        
        changeToZeroes.add(new PairInt(0, -1));
                    
        int nCorrections = 0;
       
        nCorrections += replacePattern(points, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);
        
        nCorrections += rotate90ThreeTimes(points, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    /**
     * use with caution.  it has a small window it uses to determine which
     * pixel to remove.  
     * if precision is needed, consider the zig-zag methods because they
     * use larger windows and patterns.
     * @param pixIdxs
     * @param imageWidth
     * @param imageHeight 
     */
    public void extremeStaircaseRemover(TIntSet pixIdxs, int imageWidth, 
        int imageHeight) {
       
        /*       
            0  0  0         2
            0  #< #         1
            #  #*           0
                           -1
           -1  0  1  2  3        
        */
        TIntSet ones = new TIntHashSet();
        TIntSet zeroes = new TIntHashSet();
        TIntSet changeToZeroes = new TIntHashSet();
        TIntSet changeToOnes = new TIntHashSet();
    
        int w = imageWidth;
        
        /*       
            0  0  0         2
            0  #< #         1
            #  #*           0
                           -1
           -1  0  1  2  3        
        */
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(encodePixelOffsets(-1, -2, w));//new PairInt(-1, -2));
        zeroes.add(encodePixelOffsets(0, -2, w));//new PairInt(0, -2));
        zeroes.add(encodePixelOffsets(1, -2, w));//new PairInt(1, -2));
        
        ones.add(encodePixelOffsets(-1, 0, w));//new PairInt(-1, 0)); 
        ones.add(encodePixelOffsets(0, 0, w));//new PairInt(0, 0));
        ones.add(encodePixelOffsets(0, -1, w));//new PairInt(0, -1));
        ones.add(encodePixelOffsets(1, -1, w));//new PairInt(1, -1));
        
        changeToZeroes.add(encodePixelOffsets(0, -1, w));//new PairInt(0, -1));
                    
        int nCorrections = 0;
       
        nCorrections += replacePattern(pixIdxs, imageWidth, imageHeight,
            zeroes, ones, changeToZeroes, changeToOnes);

        nCorrections += rotate90ThreeTimes(pixIdxs, imageWidth, imageHeight, 
            zeroes, ones, changeToZeroes, changeToOnes);
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(nCorrections));
    }
    
    /**
     * removes spurs, that is anything with only one neighbor.
     * Note that the last "stump" of a spur may need to be removed with 
     * another method, maybe even RemoveSpurs class.
     * @param points
     * @return 
     */
    public static Set<PairInt> removeStragglers(Set<PairInt> points) {
        
        // remove anything with only one neighbor
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        Stack<PairInt> stack = new Stack<PairInt>();
        
        Set<PairInt> rm = new HashSet<PairInt>();
        
        // place single neighbor points at top of stack
        for (PairInt p : points) {
            
            int x = p.getX();
            int y = p.getY();
            
            int nn = 0;
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                PairInt p2 = new PairInt(x2, y2);
                if (points.contains(p2)) {
                    nn++;
                }
                if (nn > 1) {
                    break;
                }
            }
            if (nn == 1) {
                stack.add(p);
            } else if (nn == 0) {
                rm.add(p);
            }
        }
        
        points.removeAll(rm);
                
        Set<PairInt> visited = new HashSet<PairInt>();
                
        while (!stack.isEmpty()) {
            
            PairInt p = stack.pop();
            if (visited.contains(p)) {
                continue;
            }
            visited.add(p);
            
            int x = p.getX();
            int y = p.getY();
            
            PairInt pn = null;
            int nn = 0;
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                PairInt p2 = new PairInt(x2, y2);
                if (points.contains(p2)) {
                    if (pn == null) {
                        pn = p2;
                    } else {
                        nn = 100;
                        break;
                    }
                }
            }
            if (nn < 2) {
                // removes if no neighbors or 1 neighbor
                points.remove(p);
                rm.add(p);
                if (pn != null) {
                    stack.add(pn);
                }
            }
        }
    
        return rm;
    }

    /**
      find boundary pattern:
     <pre>
        0  0  0   
        #  #  0  0 
        0  0  #  # 
           0  0  0  
     </pre>
     * and return the center 2 points for each such pattern
     * @param points
     * @param imageWidth
     * @param imageHeight 
     */
    public Set<PairInt> findBoundaryPattern(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
       
        /*       
           0  0  0        2
           #  #  0  0     1
           0  0  #  #     0
              0  0  0    -1
                         -2
       -3 -2 -1  0  1  2     
        */ 
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> store = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
               
        // y's are inverted here because sketch above is top left is (0,0)
        ones.add(new PairInt(-2, -1));
        ones.add(new PairInt(-1, -1));
        ones.add(new PairInt(1, 0));
        
        zeroes.add(new PairInt(-2, 0)); zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-1, 1)); zeroes.add(new PairInt(-1, 0)); 
        zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(0, -1)); zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(0, 1));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -1));

        store.add(new PairInt(-1, -1)); store.add(new PairInt(0, 0));
                
        Set<PairInt> output = findPattern(points, 
            imageWidth, imageHeight, zeroes, ones, 
            store, changeToOnes);
        
        // ----- change the sign of x to handle other direction -----
        reverseXs(zeroes, ones, store, changeToOnes);
        output.addAll(findPattern(points, 
            imageWidth, imageHeight, zeroes, ones, 
            store, changeToOnes));
        
        // ----- change the sign of y to handle other direction -----
        reverseYs(zeroes, ones, store, changeToOnes);
        output.addAll(findPattern(points, 
            imageWidth, imageHeight, zeroes, ones, 
            store, changeToOnes));
        
        // ----- change the sign of x to handle another direction -----
        reverseXs(zeroes, ones, store, changeToOnes);
        output.addAll(findPattern(points, 
            imageWidth, imageHeight, zeroes, ones, 
            store, changeToOnes));            
        
        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" + 
            Integer.toString(output.size()));
        
        return output;
    }
    
    protected int encodePixelOffsets(int x, int y, int imageWidth) {
        
        int xt = (x < 0) ? ( (-1*x << 1) | 1) : (x << 1);
        
        int yt = (y < 0) ? ( (-1*y << 1) | 1) : (y << 1);
        
        imageWidth <<= 1;
        
        int pixIdxT = (yt * imageWidth) + xt;
        
        return pixIdxT;
    }
    
    protected void decodePixelOffsets(int encodedPixIdx, int imageWidth,
        int[] outputXY) {
        
        imageWidth <<= 1;
        
        int yt = encodedPixIdx/imageWidth;
        int xt = encodedPixIdx - (yt * imageWidth);
        
        int x = ((xt & 1) == 1 ) ? -1*(xt >> 1) : xt >> 1;
        int y = ((yt & 1) == 1 ) ? -1*(yt >> 1) : yt >> 1;
        
        outputXY[0] = x;
        outputXY[1] = y;
    }
}
