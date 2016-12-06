package algorithms.imageProcessing.matching;

import algorithms.QuickSort;
import algorithms.compGeometry.LinesAndAngles;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class LinesFinder {
    
    /**
     * array of x coordinates of points found as lines.  note that some of the
     * short segments are parts of curves, indistinguishable unless larger 
     * minimum line lengths are used to filter them out.
     */
    private TIntList xs = new TIntArrayList();
    
    /**
     * array of y coordinates of points found as lines.  note that some of the
     * short segments are parts of curves, indistinguishable unless larger 
     * minimum line lengths are used to filter them out.
     */
    private TIntList ys = new TIntArrayList();
    
    /**
     * map with key = unique segment index number, value =
     * indexes of the xs and ys arrays which are a contiguous line segment.
     */
    private TIntObjectMap<TIntList> segmentIndexes = new TIntObjectHashMap<TIntList>();
    
    /**
     * map with key = unique index number, value = the polar theta in degrees
     * and distance from origin for the line segment.
     */
    private TIntObjectMap<PairInt> segmentTRMap = new TIntObjectHashMap<PairInt>();
   
    /**
     * map with key = the polar theta in degrees
     * and distance from origin for the line segment,
     * value = segment indexes of lines with this combination of
     * theta and radius.
     */
    private Map<PairInt, TIntList> trSegmentIndexesMap = new HashMap<PairInt, TIntList>();
    
    private List<PairInt> orderedTRList = null;
    private List<TIntList> orderedTRXYIndexes = null;
    
    private int lastSegIdx = -1;
    
    private int lastCol = -1;
    private int lastRow = -1;
    
    private int minLineLength = 20;
    
    private boolean debug = false;
    
    private float thresh = (float)(1.e-7);
    
    /**
     * if this is set, vertical lines found at polar radius 0 and width from
     * origin are removed and so are horizontal lines found at polar radius
     * and height from origin.
     */
    public void setToRemoveBorderLines(int lastColumn, int lastRow) {
        this.lastCol = lastColumn;
        this.lastRow = lastRow;
    }
    
    /**
     * override the default minimum line length of 30 to the given value
     * @param length 
     */
    public void overrideMinimumLineLength(int length) {
        this.minLineLength = length;
    }
    
    /**
     * override default minimum length of 20, to this value
     * @param length 
     */
    public void overrideMinimumLength(int length) {
        minLineLength = length;
    }
   
    /**
     * override the default threshold of 1.e7.  the absolute average differernce
     * of chords over a window size must be less than the threshold in order
     * for the segment to be considered a line.  Increasing this number
     * too much will possibly result in including curves.
     * @param threshold 
     */
    public void overrideThreshold(float threshold) {
        this.thresh = threshold;
    }
    
    public void setToDebug() {
        debug = true;
    }
    
    public void find(List<Set<PairInt>> listOfContigousLabels) {
        
        // -- extract the boundaries of the sets
        // -- find the lines around each boundary using shape fitting
        // -- store results as member variables
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        List<PairIntArray> listOfBounds = new ArrayList<PairIntArray>();
        
        // extract ordered bounds
        PerimeterFinder2 pFinder = new PerimeterFinder2();
        List<PairIntArray> extractedBounds = new ArrayList<PairIntArray>();
        for (int i = 0; i < listOfContigousLabels.size(); ++i) {
            Set<PairInt> set = listOfContigousLabels.get(i);
            if (set.size() < 3) {
                continue;
            }
         
            PairIntArray b = pFinder.extractOrderedBorder(new HashSet<PairInt>(
                set));
            if (b != null && b.getN() > 2) {
                listOfBounds.add(b);
            }
        }
        
        find1(listOfBounds);
    }
    
    public void find1(List<PairIntArray> listOfOrderedBounds) {
        
        // -- find the lines around each boundary using shape fitting
        // -- store results as lists of x, y 
        //    and map with key=segIdx, value=list of x,y indexes
        //    and map with key-segIdx, value = theta (float)
        //    and map with key=segIdx, value = polar radius (int)
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
                
        // extract ordered bounds
        PerimeterFinder2 pFinder = new PerimeterFinder2();
        List<PairIntArray> extractedBounds = new ArrayList<PairIntArray>();
        for (int i = 0; i < listOfOrderedBounds.size(); ++i) {
            
            PairIntArray b = listOfOrderedBounds.get(i);
            if (b == null || b.getN() < 3) {
                continue;
            }
        
            LineFinder matcher = new LineFinder();
            if (debug) {
                matcher.setToDebug();
            }
            matcher.overrideMinimumLineLength(minLineLength);
            matcher._overrideToThreshhold(thresh);
            LineFinder.LineResult r = matcher.match(b);
            List<PairInt> lr = r.getLineIndexRanges();

            for (int ii = 0; ii < lr.size(); ++ii) {
                int startIdx = lr.get(ii).getX(); 
                int stopIdx = lr.get(ii).getY(); 
                
                if (debug) {
                    System.out.println("indexes: " + startIdx + ":" + stopIdx 
                        + "   " + " segIdx=" + ii);
                }
                
                int lineX0 = b.getX(startIdx);
                int lineY0 = b.getY(startIdx);
                int lineX1 = b.getX(stopIdx);
                int lineY1 = b.getY(stopIdx);

                if (debug) {
                    System.out.println("coords: (" + lineX0 + "," + lineY0 + ") "
                        + " (" + lineX1 + "," + lineY1 + ") ");
                }
                
                double polarR = curveHelper.distanceFromPointToALine(
                    lineX0, lineY0, lineX1, lineY1, 0, 0);
                int radius = (int)Math.round(polarR);
                
                // don't store lines on image boundaries if this is set
                if (lastCol > -1) {
                    if (radius < 5) {
                        continue;
                    }
                }
                
                // -180 to 180 are the ranges.  also need an offset by 90 perpendicular
                double theta = Math.atan2(lineY1 - lineY0, lineX1 - lineX0);
                int thetaDeg = (int)Math.round(theta * 180./Math.PI);
                // perpendicular to slipe:
                thetaDeg -= 90;
                // correction to place within range 0 to 180:
                if (thetaDeg < 0) {
                    while (thetaDeg < 0) {
                        // reverse the line direction
                        thetaDeg += 180;
                    }
                } else if (thetaDeg == 180) {
                    thetaDeg = 0;
                }
                
                if (debug) {
                    System.out.println("theta, radius: (" + thetaDeg + "," 
                        + radius + ")");
                }

                // don't store lines on image bundaries if this is set               
                if (lastCol > -1) {
                    if ((thetaDeg == 0 || thetaDeg == 180) && 
                        (radius > (lastCol - 5))) {
                        continue;
                    } else if ((thetaDeg == 90 || thetaDeg == 270) 
                        && (radius > (lastRow - 5))) {
                        continue;
                    }
                }
                
                lastSegIdx++;
                
                PairInt tr = new PairInt(thetaDeg, radius);
                segmentTRMap.put(lastSegIdx, tr);
                
                int xsIdx = xs.size();
                TIntList idxs = new TIntArrayList();
                for (int j = startIdx; j <= stopIdx; ++j) {
                    int x = b.getX(j);
                    int y = b.getY(j);
                    xs.add(x);
                    ys.add(y);
                    idxs.add(xsIdx);
                    xsIdx++;
                }
                segmentIndexes.put(lastSegIdx, idxs);
                
                TIntList segIdxs = trSegmentIndexesMap.get(tr);
                if (segIdxs == null) {
                    segIdxs = new TIntArrayList();
                    trSegmentIndexesMap.put(tr, segIdxs);
                }
                segIdxs.add(lastSegIdx);
            }
        }
    }
    
    /**
     * combine the entries for a theta and radius within theta tolerance
     * and radius tolerance into lists.
     */
    public void groupWithinTolerance() {
        
        orderedTRList = new ArrayList<PairInt>();
        orderedTRXYIndexes = new ArrayList<TIntList>();
    
        Map<PairInt, TIntList> trXYIndexesMap = new
            HashMap<PairInt, TIntList>();
                
        int n = trSegmentIndexesMap.size();
        PairInt[] trs = new PairInt[n];
        int[] nPoints = new int[n];
        int[] lIdxs = new int[n];
        
        int count = 0;
        for (Entry<PairInt, TIntList> entry : trSegmentIndexesMap.entrySet()) {
            trs[count] = entry.getKey();
            TIntList segIdxs = entry.getValue();
        
            int np = 0;
            for (int j = 0; j < segIdxs.size(); ++j) {
                int segIdx = segIdxs.get(j);
                TIntList idxs = segmentIndexes.get(segIdx);
                np += idxs.size();
                
                PairInt tr = entry.getKey();
                TIntList a = trXYIndexesMap.get(tr);
                if (a == null) {
                    a = new TIntArrayList();
                    trXYIndexesMap.put(tr, a);
                }
                a.addAll(idxs);
            }
            nPoints[count] = np;
            lIdxs[count] = count;
            count++;
        }
        QuickSort.sortBy1stArg(nPoints, lIdxs);
        
        Set<PairInt> skip = new HashSet<PairInt>();
        
        for (int i = (count - 1); i > -1; --i) {
            int lIdx = lIdxs[i];
            int np = nPoints[i];
            
            PairInt tr = trs[lIdx];
            
            if (skip.contains(tr)) {
                continue;
            }

            TIntList xyList = new TIntArrayList();

            int sumT = 0;
            int sumR = 0;
            int sumN = 0;
           
            for (int t0 = tr.getX() - 2; t0 <= tr.getX() + 2; ++t0) {
                for (int r0 = tr.getY() - 2; r0 <= tr.getY() + 2; ++r0) {
                    PairInt tr0 = new PairInt(t0, r0);
                    TIntList a = trXYIndexesMap.get(tr0);
                    if (a == null) {
                        continue;
                    }
                    xyList.addAll(a);
                    trXYIndexesMap.remove(tr0);
                    skip.add(tr0);
                    sumT += tr0.getX();
                    sumR += tr0.getY();
                    sumN++;
                }
            }
            sumR /= sumN;
            sumT /= sumN;
            PairInt tr0 = new PairInt(sumT, sumR);
            orderedTRList.add(tr0);
            orderedTRXYIndexes.add(xyList);
        }
    }
    
    public void debugPrintTRStats() {
        
        int n = trSegmentIndexesMap.size();
        PairInt[] trs = new PairInt[n];
        int[] nLines = new int[n];
        int[] nPoints = new int[n];
        int[] lIdxs = new int[n];
        
        int count = 0;
        for (Entry<PairInt, TIntList> entry : trSegmentIndexesMap.entrySet()) {
            trs[count] = entry.getKey();
            TIntList segIdxs = entry.getValue();
            nLines[count] = segIdxs.size();
            
            int np = 0;
            for (int j = 0; j < segIdxs.size(); ++j) {
                int segIdx = segIdxs.get(j);
                np += segmentIndexes.get(segIdx).size();
            }
            nPoints[count] = np;
            lIdxs[count] = count;
            count++;
        }
        QuickSort.sortBy1stArg(nPoints, lIdxs);
        
        for (int i = (count - 1); i > -1; --i) {
            int lIdx = lIdxs[i];
            System.out.println(String.format("np=%d nL=%d tr=%s", 
                nPoints[i], nLines[lIdx], trs[lIdx].toString()));
        }
    }
    
    public void debugDraw(algorithms.imageProcessing.Image img) {
        
        if (orderedTRList == null) {
            throw new IllegalStateException("groupWithinTolerance must be "
                + " invoked first");
        }
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        int end = 10;
        //if (end > (orderedTRList.size() - 1)) {
            end = orderedTRList.size();
        //}
        //for (int i = 0; i < orderedTRList.size(); ++i) {
        for (int i = 0; i < end; ++i) {
            
            PairInt tr = orderedTRList.get(i);
          
            int[] clr = ImageIOHelper.getNextRGB(i);
            
            boolean drawLines = false;
          
            if (drawLines) {
                int[] eps = LinesAndAngles.calcPolarLineEndPoints(
                    tr.getX(), tr.getY(), img.getWidth(), img.getHeight());

                System.out.println("tr=" + tr.toString() + " eps=" +
                    Arrays.toString(eps) + " w=" + img.getWidth() + 
                    " h=" + img.getHeight());

                ImageIOHelper.drawLineInImage(
                    eps[0], eps[1], eps[2], eps[3], img, 1, 
                    clr[0],clr[1], clr[2]);            
            } else {
                TIntList idxs = orderedTRXYIndexes.get(i);
                for (int k = 0; k < idxs.size(); ++k) {
                    int x = xs.get(k);
                    int y = ys.get(k);
                    ImageIOHelper.addPointToImage(x, y, img, 1, 
                        clr[0], clr[1], clr[2]);
                }
                System.out.println("  tr=" + tr + " n=" + idxs.size());
            }
        }
    }
}
