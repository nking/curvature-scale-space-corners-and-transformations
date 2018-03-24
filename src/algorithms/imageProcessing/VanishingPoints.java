package algorithms.imageProcessing;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.imageProcessing.matching.LinesFinder;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.QuadInt;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * 
 * NOT READY FOR USE. 
 * a class to hold various methods for determining vanishing lines
 * and to hold the resulting vanishing points.
 * 
 * NOTE: should add an MSER implementation in here.  see snapshots in 
 * docs/colorSegmentation3.pdf to see how well the MSER regions
 * find the major vanishing lines.
 * 
 * The points require a 3d model so methods may be added for them
 * at a later time. 
 * meanwhile, the user can retrieve the vanishing points 
 * if any for a segmented cell.
 * 
 * @author nichole
 */
public class VanishingPoints {
        
    private boolean debug = false;
    
    /*
    key = segment idx,
       value = line endpoints and vanishing point
    */
    private TIntObjectMap<QuadInt> vanishingLines =
        new TIntObjectHashMap<QuadInt>();
    
    public PairInt[] points;
    
    public void setToDebug() {
        debug = true;
    }
public Image dbgImg = null;    
    public void find(List<Set<PairInt>> listOfContigousLabels,
        int imageWidth, int imageHeight) {
        
        LinesFinder finder = new LinesFinder();
        if (debug) {
            finder.setToDebug();
        }
        
        finder.setToRemoveBorderLines(imageWidth - 1, imageHeight - 1);
        finder.find(listOfContigousLabels);
        finder.groupWithinTolerance();

        finder.sortOrderedLists();
        
        if (dbgImg != null) {
            finder.debugDraw(dbgImg);
            MiscDebug.writeImage(dbgImg, "_all_" + MiscDebug.getCurrentTimeFormatted());
        }
    
        List<PairInt> orderedTRList = finder.getOrderedTRList();
        List<TIntList> orderedTRXYIndexes = finder.getOrderedTRXYIndexes();
        TIntList xs = finder.getXs();
        TIntList ys = finder.getYs();

        // key = point, value = segmented cell index
        TObjectIntMap<PairInt> pointSegmentedIndexMap = 
            new TObjectIntHashMap<PairInt>();
        for (int i = 0; i < listOfContigousLabels.size(); ++i) {
            Set<PairInt> set = listOfContigousLabels.get(i);
            for (PairInt p : set) {
                pointSegmentedIndexMap.put(p, i);
            }
        }
        
        // key = point, value = ordered tr index
        TObjectIntMap<PairInt> pointTRIndexMap = 
            new TObjectIntHashMap<PairInt>();
        for (int i = 0; i < orderedTRXYIndexes.size(); ++i) {
            TIntList idxs = orderedTRXYIndexes.get(i);
            for (int j = 0; j < idxs.size(); ++j) {
                int x = xs.get(j);
                int y = ys.get(j);
                pointTRIndexMap.put(new PairInt(x, y), i);
            }
        }
        
        // key = segmented cell index, value = ordered TR index
        TIntObjectMap<TIntList> segmentedLineIndexes = 
            new TIntObjectHashMap<TIntList>();
        for (int i = 0; i < xs.size(); ++i) {
            PairInt p = new PairInt(xs.get(i), ys.get(i));
            assert(pointSegmentedIndexMap.containsKey(p));
            int segIdx = pointSegmentedIndexMap.get(p);
            
            TIntList idxs = segmentedLineIndexes.get(segIdx);
            if (idxs == null) {
                idxs = new TIntArrayList();
                segmentedLineIndexes.put(segIdx, idxs);
            }
            int trIdx = pointTRIndexMap.get(p);
            idxs.add(trIdx);            
        }
                
        TIntObjectIterator<TIntList> iter = segmentedLineIndexes.iterator();
        for (int i = 0; i < segmentedLineIndexes.size(); ++i) {
            iter.advance();
            int segIdx = iter.key();
            TIntList trIdxs = iter.value();
            
            TIntSet uniqueT = new TIntHashSet();
            List<PairInt> trs = new ArrayList<PairInt>();
            for (int j = 0; j < trIdxs.size(); ++j) {
                int trIdx = trIdxs.get(j);                
                PairInt tr = orderedTRList.get(trIdx);
                
                if (uniqueT.contains(tr.getX())) {
                    continue;
                }
                uniqueT.add(tr.getX());
                trs.add(tr);
            }
            if (!trs.isEmpty()) {
                for (int j = 0; j < trs.size(); ++j) {
                    PairInt tr = trs.get(j);
                    int[] endpoints = LinesAndAngles.calcPolarLineEndPoints(
                        tr.getX(), tr.getY(), imageWidth, imageWidth);
                    
                    vanishingLines.put(segIdx, 
                        new QuadInt(endpoints[0], endpoints[1],
                        endpoints[2], endpoints[3]));
                }
            }
        }        
    
    }

    /**
     * @return the vanishingLines
     */
    public TIntObjectMap<QuadInt> getVanishingLines() {
        return vanishingLines;
    }
}
