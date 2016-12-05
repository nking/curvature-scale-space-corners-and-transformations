package algorithms.imageProcessing.matching;

import algorithms.QuickSort;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import java.util.ArrayList;
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
    
    private int lastSegIdx = -1;
    
    private int lastCol = -1;
    private int lastRow = -1;
    
    /**
     * if this is set, vertical lines found at polar radius 0 and width from
     * origin are removed and so are horizontal lines found at polar radius
     * and height from origin.
     */
    public void setToRemoveBorderLines(int lastColumn, int lastRow) {
        this.lastCol = lastColumn;
        this.lastRow = lastRow;
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
            //matcher.setToDebug();
            LineFinder.LineResult r = matcher.match(b);
            List<PairInt> lr = r.getLineIndexRanges();

            for (int ii = 0; ii < lr.size(); ++ii) {
                int startIdx = lr.get(ii).getX(); 
                int stopIdx = lr.get(ii).getY(); 
                System.out.println(startIdx + ":" + stopIdx + "   " + " segIdx=" + ii);

                int lineX0 = b.getX(startIdx);
                int lineY0 = b.getY(startIdx);
                int lineX1 = b.getX(stopIdx);
                int lineY1 = b.getY(stopIdx);

                double polarR = curveHelper.distanceFromPointToALine(
                    lineX0, lineY0, lineX1, lineY1, 0, 0);
                int radius = (int)Math.round(polarR);
                
                // don't store lines on image boundaries if this is set
                if (lastCol > -1) {
                    if (radius == 0) {
                        continue;
                    }
                }
                
                double theta = Math.atan2(lineY1 - lineY0, lineX1 - lineX0);
                int thetaDeg = (int)Math.round(theta * 180./Math.PI);
                if (thetaDeg < 0) {
                    thetaDeg += 360;
                }

                // don't store lines on image bundaries if this is set               
                if (lastCol > -1) {
                    if ((thetaDeg == 0 || thetaDeg == 180) && (lastRow == radius)) {
                        continue;
                    } else if ((thetaDeg == 90 || thetaDeg == 270) && (lastCol == radius)) {
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
}
