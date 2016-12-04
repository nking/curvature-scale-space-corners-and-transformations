package algorithms.imageProcessing.matching;

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
                
                double theta = Math.atan2(lineY1 - lineY0, lineX1 - lineX0);
                int thetaDeg = (int)Math.round(theta * 180./Math.PI);
                if (thetaDeg < 0) {
                    thetaDeg += 360;
                }
                
                lastSegIdx++;
                
                PairInt tr = new PairInt(thetaDeg, (int)Math.round(polarR));
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
}
