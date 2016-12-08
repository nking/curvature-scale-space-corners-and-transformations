package algorithms.imageProcessing.segmentation;

import algorithms.imageProcessing.DFSConnectedGroupsFinder;
import algorithms.imageProcessing.GreyscaleImage;
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
     * reassign those pixels the average colors.
     * @param img
     * @param labels 
     */
    public static void applyLabels(Image img, int[] labels) {
        
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
    
    /**
     extract contiguous points from the labeled regions and relabel
    labels to coincide with the returned list indexes
    */
    public static List<Set<PairInt>> extractContiguousLabelPoints(Image img, 
        int[] labels) {
                
        TIntObjectMap<Set<PairInt>> lMap = extractLabelPoints(img, labels);
              
        List<Set<PairInt>> out = extractContiguousLabelPoints(img, lMap);
        
        for (int i = 0; i < out.size(); ++i) {
            for (PairInt p : out.get(i)) {
                int pixIdx = img.getInternalIndex(p);
                labels[pixIdx] = i;
                assert(img.getCol(pixIdx) < img.getWidth());
                assert(img.getRow(pixIdx) < img.getHeight());
            }
        }
        
        return out;
    }
    
    /**
     extract contiguous points from the labeled regions and relabel
    labels to coincide with the returned list indexes
    * @param labels 2D array in format [row][col]
    */
    public static List<Set<PairInt>> extractContiguousLabelPoints(Image img, 
        int[][] labels) {
                
        TIntObjectMap<Set<PairInt>> lMap = extractRowMajorLabelPoints(img, labels);
              
        List<Set<PairInt>> out = extractContiguousLabelPoints(img, lMap);
        
        for (int i = 0; i < out.size(); ++i) {
            for (PairInt p : out.get(i)) {
                labels[p.getX()][p.getY()] = i;
                assert(p.getX() < img.getWidth());
                assert(p.getY() < img.getHeight());
            }
        }
       
        return out;
    }
    
    /**
     extract contiguous points from the labeled regions and relabel
    labels to coincide with the returned list indexes
    */
    private static List<Set<PairInt>> extractContiguousLabelPoints(Image img, 
        TIntObjectMap<Set<PairInt>> lMap) {
        
        List<Set<PairInt>> out = new ArrayList<Set<PairInt>>();
              
        TIntObjectIterator<Set<PairInt>> iter = lMap.iterator();
        for (int i = 0; i < lMap.size(); ++i) {
            iter.advance();
            Set<PairInt> set = iter.value();
            DFSConnectedGroupsFinder finder = new DFSConnectedGroupsFinder();
            // setting is for 4 neighbors
            finder.setMinimumNumberInCluster(1);
            finder.findConnectedPointGroups(set);
            for (int j = 0; j < finder.getNumberOfGroups(); ++j) {
                Set<PairInt> group = finder.getXY(j);
                out.add(group);
            }
        }
 
        assert(assertAllPointsFound(out, img.getWidth(), img.getHeight()));
        
        return out;
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
        int[] labels) {
        
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
        
        assert(labels.length == img.getNPixels());
        
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
    
    /**
     * create map with key = label, value = all points w/ label
     * @param img
     * @param labels two dimensional array of labels in format [row][col]
     * @return 
     */
    public static TIntObjectMap<Set<PairInt>> extractRowMajorLabelPoints(
        Image img, int[][] labels) {
        
        assert(labels.length*labels[0].length == img.getNPixels());
        
        TIntObjectMap<Set<PairInt>> out 
            = new TIntObjectHashMap<Set<PairInt>>();
        
        for (int j = 0; j < labels.length; ++j) {
            for (int i = 0; i < labels[j].length; ++i) {
                int label = labels[j][i];
                Set<PairInt> set = out.get(label);
                if (set == null) {
                    set = new HashSet<PairInt>();
                    out.put(label, set);
                }
                set.add(new PairInt(img.getCol(i), img.getRow(i)));
            }
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
     * create a map with key = label, value = set of adjacent
     * labels.
     * Note, this uses a 4 nieghbor region, not 8.
     * @param img
     * @param labels
     * @param excludeNegativeLabels
     * @return 
     */
    public static TIntObjectMap<TIntSet> createAdjacencyLabelMap(
        Image img, int[] labels, boolean excludeNegativeLabels) {
        
        TIntObjectMap<TIntSet> adjacencyMap =
            new TIntObjectHashMap<TIntSet>();

        int h = img.getHeight();
        int w = img.getWidth();
        
        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;
        
        for (int idx1 = 0; idx1 < labels.length; ++idx1) {
            
            int l1 = labels[idx1];
            
            if (excludeNegativeLabels && (l1 < 0)) {
                continue;
            }
            
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
                if (excludeNegativeLabels && (l2 < 0)) {
                    continue;
                }
                set1.add(l2);
            }
        }
        
        return adjacencyMap;
    }

    /**
     * create a map with key = label, value = set of adjacent
     * labels.
     * Note that it uses a 4 neighbor region, not 8.
     * @param img
     * @param labels
     * @param excludeNegativeLabels
     * @return 
     */
    public static TIntObjectMap<TIntSet> createAdjacencyLabelMap(
        GreyscaleImage img, int[] labels, boolean excludeNegativeLabels) {
        
        TIntObjectMap<TIntSet> adjacencyMap =
            new TIntObjectHashMap<TIntSet>();

        int h = img.getHeight();
        int w = img.getWidth();
        
        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;
        
        for (int idx1 = 0; idx1 < labels.length; ++idx1) {
            
            int l1 = labels[idx1];
            
            if (excludeNegativeLabels && (l1 < 0)) {
                continue;
            }
            
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
                if (excludeNegativeLabels && (l2 < 0)) {
                    continue;
                }
                set1.add(l2);
            }
        }
        
        return adjacencyMap;
    }

    /**
     * create a map with key = label, value = set of adjacent
     * labels.
     * Note uses 4 neighbor region, not 8.
     * @param img
     * @param labels
     * @return 
     */
    public static TIntObjectMap<TIntSet> createAdjacencyLabelMap(
        GreyscaleImage img, int[] labels) {
        
        TIntObjectMap<TIntSet> adjacencyMap =
            new TIntObjectHashMap<TIntSet>();

        int h = img.getHeight();
        int w = img.getWidth();
        
        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;
        
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

    public static boolean assertAllPointsFound(
        List<Set<PairInt>> listOfSets, int width, int height) {

        Set<PairInt> allPoints = new HashSet<PairInt>();
        for (Set<PairInt> set : listOfSets) {
            for (PairInt p : set) {
                boolean exists = allPoints.contains(p);
                if (exists) {
                    return false;
                }
                allPoints.add(p);
            }
        }
        if (allPoints.size() != width * height) {
            return false;
        }
        
        for (int i = 0; i < width; ++i) {
            for (int j = 0; j < height; ++j) {
                PairInt p = new PairInt(i, j);
                boolean exists = allPoints.contains(p);
                if (!exists) {
                    return false;
                }
                boolean rmvd = allPoints.remove(p);
                assert(rmvd);
            }
        }
        
        return allPoints.isEmpty();        
    }
    
    /**
     * create a map w/ key = contiguousSets index, value =
     *   set of indexes of contiguousSets which are adjacent to key's set.
     * Note that int[] labels usually have a value that is an index
     * of contiguousSets.
     * Note also that it uses a 4 neighbor region, not 8.
     * @param contiguousSets
     * @return 
     */
    public static TIntObjectMap<TIntSet> createAdjacencySetMap(
        List<Set<PairInt>> contiguousSets) {

        TObjectIntMap<PairInt> pointIndexMap =
            new TObjectIntHashMap<PairInt>();
        for (int i = 0; i < contiguousSets.size(); ++i) {
            for (PairInt p : contiguousSets.get(i)) {
                pointIndexMap.put(p, i);
            }
        }

        TIntObjectMap<TIntSet> contigAdjacencyMap =
            new TIntObjectHashMap<TIntSet>();

        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;
        for (int lIdx = 0; lIdx < contiguousSets.size(); ++lIdx) {

            TIntSet setIdx1 = contigAdjacencyMap.get(lIdx);
            if (setIdx1 == null) {
                setIdx1 = new TIntHashSet();
                contigAdjacencyMap.put(lIdx, setIdx1);
            }

            Set<PairInt> set = contiguousSets.get(lIdx);
            for (PairInt p : set) {
                int x = p.getX();
                int y = p.getY();
                for (int j = 0; j < dxs.length; ++j) {
                    int x2 = x + dxs[j];
                    int y2 = y + dys[j];
                    PairInt p2 = new PairInt(x2, y2);
                    if (!pointIndexMap.containsKey(p2)) {
                        continue;
                    }
                    int lIdx2 = pointIndexMap.get(p2);
                    if (lIdx == lIdx2) {
                        continue;
                    }
                
                    setIdx1.add(lIdx2);
                }
            }
        }

        for (int lIdx1 = 0; lIdx1 < contiguousSets.size(); ++lIdx1) {
            TIntSet setIdx1 = contigAdjacencyMap.get(lIdx1);
            setIdx1.remove(lIdx1);
        }

        return contigAdjacencyMap;
    }

}
