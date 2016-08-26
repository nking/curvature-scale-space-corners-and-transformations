package algorithms.imageProcessing.segmentation;

import algorithms.imageProcessing.DFSConnectedGroupsFinder;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
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
}