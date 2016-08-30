package algorithms.imageProcessing.segmentation;

import algorithms.imageProcessing.DFSConnectedGroupsFinder;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class LabelToColorHelper {
    
    /**
     * calculate the average r,g,b of pixels grouped by their labels and
     * reassigne those pixels the average colors.
     * @param img
     * @param labels 
     */
    public static void applyLabels(ImageExt img, int[] labels) {
        
        if (img.getNPixels() != labels.length) {
            throw new IllegalArgumentException("labels.length must equal img.nPixels");
        }
        
        int maxLabel = MiscMath.findMax(labels);
        
        long[] rSum = new long[maxLabel + 1];
        long[] gSum = new long[maxLabel + 1];
        long[] bSum = new long[maxLabel + 1];
        
        int[] count = new int[rSum.length];
        
        for (int i = 0; i < labels.length; ++i) {
            int label = labels[i];
            rSum[label] += img.getR(i);
            gSum[label] += img.getG(i);
            bSum[label] += img.getB(i);
            
            count[label]++;
        }
        
        for (int i = 0; i < rSum.length; ++i) {
            if (count[i] > 0) {
                rSum[i] /= count[i];
                gSum[i] /= count[i];
                bSum[i] /= count[i];
            }
        }
        
        img.fill(0, 0, 0);
        
        for (int i = 0; i < labels.length; ++i) {
            int label = labels[i];
            img.setRGB(i, (int)rSum[label], (int)gSum[label], (int)bSum[label]);
        }
    }
    
    public static List<Set<PairInt>> 
        extractContiguousLabelPoints(Image img, int[] labels) {
        
        List<Set<PairInt>> out = new ArrayList<Set<PairInt>>();
        
        TIntObjectMap<Set<PairInt>> lMap = extractLabelPoints(img, labels);
      
        TIntObjectIterator<Set<PairInt>> iter = lMap.iterator();
        for (int i = 0; i < lMap.size(); ++i) {
            iter.advance();
            Set<PairInt> set = iter.value();
            DFSConnectedGroupsFinder finder = 
                new DFSConnectedGroupsFinder();
            finder.setMinimumNumberInCluster(1);
            finder.findConnectedPointGroups(set);
            for (int j = 0; j < finder.getNumberOfGroups(); ++j) {
                Set<PairInt> group = finder.getXY(j);
                out.add(group);
            }
        }
        
        return out;
    }
        
    /**
     * separate each labeled group to only contain 
     * contiguous members, that is relabel to further
     * segment by connectedness.
     * @param img
     * @param labels 
     */
    public static void applyContiguousSegmentation(Image img, 
        int[] labels) {

        assert(labels.length == img.getNPixels());
        
        int[] out = new int[labels.length];        
        
        TIntObjectMap<Set<PairInt>> lMap = 
            extractLabelPoints(img, labels);
        
        int count = 0;
        
        TIntObjectIterator<Set<PairInt>> iter = lMap.iterator();
        for (int i = 0; i < lMap.size(); ++i) {
            iter.advance();
            Set<PairInt> set = iter.value();
            DFSConnectedGroupsFinder finder = 
                new DFSConnectedGroupsFinder();
            finder.setMinimumNumberInCluster(1);
            finder.findConnectedPointGroups(set);
            for (int j = 0; j < finder.getNumberOfGroups(); ++j) {
                Set<PairInt> group = finder.getXY(j);
                for (PairInt p : group) {
                    int pixIdx = img.getInternalIndex(p);
                    out[pixIdx] = count;
                }
                count++;
            }
        }
        
        System.arraycopy(out, 0, labels, 0, labels.length);
    }
    
    /**
     * create map w/ key = label, index = all indexes in labels
     * with that label.  Note that the label indexes are
     * the pixel indexes also.
     * 
     * @param img
     * @param labels
     * @return map w/ key = label, index = all indexes in labels
     * with that label.  Note that the label indexes are
     * the pixel indexes also.
     */
    public static TIntObjectMap<TIntSet> createLabelIndexMap(
        Image img, int[] labels) {
        
        TIntObjectMap<TIntSet> map = new TIntObjectHashMap<TIntSet>();
        
        for (int i = 0; i < labels.length; ++i) {
            int label = labels[i];
            TIntSet set = map.get(label);
            if (set == null) {
                set = new TIntHashSet();
                map.put(label, set);
            }
            set.add(i);
        }
        
        return map;
    }
    
    /**
     * create map with key = label, value = all points w/ label
     * @param img
     * @param labels
     * @return 
     */
    public static TIntObjectMap<Set<PairInt>> extractLabelPoints(
        Image img, int[] labels) {
        
        TIntObjectMap<Set<PairInt>> out 
            = new TIntObjectHashMap<Set<PairInt>>();
        
        for (int i = 0; i < labels.length; ++i) {
            int label = labels[i];
            Set<PairInt> set = out.get(label);
            if (set == null) {
                set = new HashSet<PairInt>();
                out.put(label, set);
            }
            set.add(new PairInt(img.getCol(i), img.getRow(i)));
        }
        
        return out;
    }
    
    public static int[] createLabelsFromContiguousSets(
        List<Set<PairInt>> sets, Image img) {
        
        int[] labels = new int[img.getNPixels()];
        int count = 0;
        for (int i = 0; i < sets.size(); ++i) {
            for (Set<PairInt> set : sets) {
                for (PairInt p : set) {
                    int pixIdx = img.getInternalIndex(p.getX(), p.getY());
                    labels[pixIdx] = count;
                }
            }
            count++;
        }
        assert(sets.size() == count);
        
        return labels;
    }

    /**
     * create a map with key = point index, value = 
     * a set of the 8 adjacent neighbors which have a
     * different label.
     * 
     * @param img
     * @param labels
     * @return 
     */
    public static TIntObjectMap<TIntSet> createAdjacencyPointMap(
        ImageExt img, int[] labels) {
        
        TIntObjectMap<TIntSet> adjacencyMap =
            new TIntObjectHashMap<TIntSet>();

        int h = img.getHeight();
        int w = img.getWidth();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        for (int idx1 = 0; idx1 < labels.length; ++idx1) {
            
            TIntSet setIdx1 = new TIntHashSet();
            adjacencyMap.put(idx1, setIdx1);
            
            int l1 = labels[idx1];
            
            int x = img.getCol(idx1);
            int y = img.getRow(idx1);
            for (int j = 0; j < dxs.length; ++j) {
                int x2 = x + dxs[j];
                int y2 = y + dys[j];
                if (x2 < 0 || y2 < 0 || (x2 > (w - 1) ||
                    (y2 > (h - 1)))) {
                    continue;
                }
                int idx2 = img.getInternalIndex(x2, y2);
                int l2 = labels[idx2];
                if (l1 == l2) {
                    continue;
                }
                setIdx1.add(idx2);
            }
        }
        
        return adjacencyMap;
    }
    
    public static TIntObjectMap<TIntSet> createAdjacencyLabelMap(
        ImageExt img, int[] labels) {
        
        TIntObjectMap<TIntSet> adjacencyMap =
            new TIntObjectHashMap<TIntSet>();

        int h = img.getHeight();
        int w = img.getWidth();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        for (int idx1 = 0; idx1 < labels.length; ++idx1) {
            
            int l1 = labels[idx1];
            
            TIntSet set1 = adjacencyMap.get(l1);
            if (set1 == null) {
                set1 = new TIntHashSet();
                adjacencyMap.put(l1, set1);
            }
                        
            int x = img.getCol(idx1);
            int y = img.getRow(idx1);
            for (int j = 0; j < dxs.length; ++j) {
                int x2 = x + dxs[j];
                int y2 = y + dys[j];
                if (x2 < 0 || y2 < 0 || (x2 > (w - 1) ||
                    (y2 > (h - 1)))) {
                    continue;
                }
                int idx2 = img.getInternalIndex(x2, y2);
                int l2 = labels[idx2];
                if (l1 == l2) {
                    continue;
                }
                set1.add(l2);
            }
        }
        
        return adjacencyMap;
    }

    public static void condenseLabels(int[] labels) {
        
        int count = 0;
        TIntIntMap map = new TIntIntHashMap();
        for (int i = 0; i < labels.length; ++i) {
            int label = labels[i];
            int label2;
            if (map.containsKey(label)) {
                label2 = map.get(label);
            } else {
                label2 = count;
                count++;
                map.put(label, label2);
            }
            labels[i] = label2;
        }
    }
    
}