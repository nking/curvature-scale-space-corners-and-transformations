package algorithms.imageProcessing.segmentation;

import algorithms.imageProcessing.Image;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class LabelToColorHelperTest extends TestCase {
    
    public LabelToColorHelperTest() {
    }

    public void testApplyLabels() {
        int w = 5;
        int h = 7;
        
        Image img = new Image(w, h);
        
        int[] labels = new int[img.getNPixels()];
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int pixIdx = img.getInternalIndex(i, j);
                img.setRGB(pixIdx, i, i, i);
                labels[pixIdx] = i;
            }
            img.setRGB(i, h/2, i + 1, i + 1, i + 1);
        }
        
        LabelToColorHelper.applyLabels(img, labels);
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int pixIdx = img.getInternalIndex(i, j);
                assertEquals(i, img.getR(pixIdx));
                assertEquals(i, img.getG(pixIdx));
                assertEquals(i, img.getB(pixIdx));
            }
        }
    }
    
    public void testExtractContiguousLabelPoints() {
        
        /*
          0 0  1  1  0  0 
          0 0  1  1  0  0
        */
        
        int w = 6;
        int h = 2;
        
        Image img = new Image(w, h);
        int[] labels = new int[img.getNPixels()];
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int label = 0;
                if (i == 2 || i == 3) {
                    label = 1;
                }
                int pixIdx = img.getInternalIndex(i, j);
                labels[pixIdx] = label;
                img.setRGB(i, j, label, label, label);
            }
        }
        
        List<Set<PairInt>> list = LabelToColorHelper
            .extractContiguousLabelPoints(img, labels);
    
        /*
          0 0  1  1  0  0 
          0 0  1  1  0  0
        */
        assertEquals(3, list.size());
        
        for (Set<PairInt> set : list) {
            assertEquals(4, set.size());
            PairInt p0 = set.iterator().next();
            TIntList ys = new TIntArrayList();
            ys.add(0);
            ys.add(1);
            TIntList xs = new TIntArrayList();
            if (p0.getX() == 2 || p0.getX() == 3) {
                xs.add(2);
                xs.add(3);                
            } else if (p0.getX() == 0 || p0.getX() == 1) {
                xs.add(0);
                xs.add(1);
            } else if (p0.getX() == 4 || p0.getX() == 5) {
                xs.add(4);
                xs.add(5);
            } else {
                fail("unexpected x in p0: " + p0.toString());
            }
            for (int i = 0; i < xs.size(); ++i) {
                for (int j = 0; j < ys.size(); ++j) {
                    PairInt p = new PairInt(xs.get(i), ys.get(j));
                    assertTrue(set.contains(p));
                }
            }
        }
    }
    
    public void testApplyContiguousSegmentation() {
        
         /*
          0 0  1  1  0  0 
          0 0  1  1  0  0
        */
        
        int w = 6;
        int h = 2;
        
        Image img = new Image(w, h);
        int[] labels = new int[img.getNPixels()];
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int label = 0;
                if (i == 2 || i == 3) {
                    label = 1;
                }
                int pixIdx = img.getInternalIndex(i, j);
                labels[pixIdx] = label;
                img.setRGB(i, j, label, label, label);
            }
        }
        
        List<Set<PairInt>> list = 
            LabelToColorHelper.extractContiguousLabelPoints(img, labels);
        
        // spot checks for now...
        assertTrue(img.getR(0, 0) == img.getR(4, 0));
        assertTrue(img.getB(1, 1) != img.getB(3, 1));
        assertTrue(img.getG(2, 0) == img.getG(3, 1));
    }
  
    public void testCreateLabelIndexMap() {
        
        /*
          0 0  1  1  0  0 
          0 0  1  1  0  0
        */
        
        int w = 6;
        int h = 2;
        
        Image img = new Image(w, h);
        int[] labels = new int[img.getNPixels()];
        
        TIntSet set0 = new TIntHashSet();
        TIntSet set1 = new TIntHashSet();
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int pixIdx = img.getInternalIndex(i, j);
                int label = 0;
                if (i == 2 || i == 3) {
                    label = 1;
                    set1.add(pixIdx);
                } else {
                    set0.add(pixIdx);
                }
                labels[pixIdx] = label;
                img.setRGB(i, j, label, label, label);
            }
        }
        
        
        TIntObjectMap<TIntSet> labelToIndexMap = 
            LabelToColorHelper.createLabelIndexMap(img, labels);
        
        assertEquals(8, labelToIndexMap.get(0).size());
            
        assertEquals(4, labelToIndexMap.get(1).size());
    
        TIntIterator iter = labelToIndexMap.get(0).iterator();
        while (iter.hasNext()) {
            int idx = iter.next();
            assertTrue(set0.remove(idx));
        }
        assertEquals(0, set0.size());
        
        iter = labelToIndexMap.get(1).iterator();
        while (iter.hasNext()) {
            int idx = iter.next();
            assertTrue(set1.remove(idx));
        }
        assertEquals(0, set1.size());
    }
    
    public void testExtractLabelPoints() {
        
        /*
          0 0  1  1  0  0 
          0 0  1  1  0  0
        */
        
        int w = 6;
        int h = 2;
        
        Image img = new Image(w, h);
        int[] labels = new int[img.getNPixels()];
        
        Set<PairInt> set0 = new HashSet<PairInt>();
        Set<PairInt> set1 = new HashSet<PairInt>();
       
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int pixIdx = img.getInternalIndex(i, j);
                int label = 0;
                if (i == 2 || i == 3) {
                    label = 1;
                    set1.add(new PairInt(i, j));
                } else {
                    set0.add(new PairInt(i, j));
                }
                labels[pixIdx] = label;
                img.setRGB(i, j, label, label, label);
            }
        }
        
        TIntObjectMap<Set<PairInt>> labelPointsMap =
            LabelToColorHelper.extractLabelPoints(img, labels);
        
        assertEquals(8, labelPointsMap.get(0).size());
            
        assertEquals(4, labelPointsMap.get(1).size());
    
        Iterator<PairInt> iter  = labelPointsMap.get(0).iterator();
        while (iter.hasNext()) {
            PairInt p = iter.next();
            assertTrue(set0.remove(p));
        }
        assertEquals(0, set0.size());
        
        iter  = labelPointsMap.get(1).iterator();
        while (iter.hasNext()) {
            PairInt p = iter.next();
            assertTrue(set1.remove(p));
        }
        assertEquals(0, set1.size());
    }
    
    public void testCreateLabelFromContiguousSets() {
      
        /*
          0 0  1  1  0  0 
          0 0  1  1  0  0
        */
        
        int w = 6;
        int h = 2;
        
        Image img = new Image(w, h);
        
        List<Set<PairInt>> list = new ArrayList<Set<PairInt>>();
        Set<PairInt> set0 = new HashSet<PairInt>();
        Set<PairInt> set1 = new HashSet<PairInt>();
        Set<PairInt> set2 = new HashSet<PairInt>();
        list.add(set0); list.add(set1); list.add(set2);
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int label = 0;
                if (i == 0 || i == 1) {
                    set0.add(new PairInt(i, j));
                } else if (i == 2 || i == 3) {
                    label = 1;
                    set1.add(new PairInt(i, j));
                } else {
                    set2.add(new PairInt(i, j));
                }
                img.setRGB(i, j, label, label, label);
            }
        }
        
        int[] labels = LabelToColorHelper.
            createLabelsFromContiguousSets(list, img);
    
        assertEquals(img.getNPixels(), labels.length);
        
        int label0 = -1;
        int label1 = -1;
        int label2 = -1;
        int n0 = 0;
        int n1 = 0;
        int n2 = 0;
        for (int pixIdx = 0; pixIdx < labels.length; ++pixIdx) {
            int i = img.getCol(pixIdx);
            int j = img.getRow(pixIdx);
            int label = labels[pixIdx];
            PairInt p = new PairInt(i, j);
            if (set0.contains(p)) {
                n0++;
                if (label0 == -1) {
                    label0 = label;
                } else {
                    assert(label0 == label);
                }
            } else if (set1.contains(p)) {
                n1++;
                if (label1 == -1) {
                    label1 = label;
                } else {
                    assert(label1 == label);
                }
            } else if (set2.contains(p)) {
                n2++;
                if (label2 == -1) {
                    label2 = label;
                } else {
                    assert(label2 == label);
                }
            } else {
                fail("p not in set0,1,or 2: " + p.toString());
            }
        }
        assertEquals(4, n0);
        assertEquals(4, n1);
        assertEquals(4, n2);
    }
    
    public void testCreateAdjacencyLabelMap() {
        
        /*
          0 0  1  1  0  0 
          0 0  1  1  0  0
        */
        
        int w = 6;
        int h = 2;
        
        Image img = new Image(w, h);
        
        int[] labels = new int[img.getNPixels()];
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int label = 0;
                if (i == 0 || i == 1) {
                    label = 0;
                } else if (i == 2 || i == 3) {
                    label = 1;
                } else {
                    label = 2;
                }
                int pixIdx = img.getInternalIndex(i, j);
                labels[pixIdx] = label;
                img.setRGB(i, j, label, label, label);
            }
        }
        
        TIntObjectMap<TIntSet> adjLabelMap = 
            LabelToColorHelper.createAdjacencyLabelMap(img, 
                labels, true);
        
        assertEquals(3, adjLabelMap.size());
        assertEquals(1, adjLabelMap.get(0).size());
        assertEquals(2, adjLabelMap.get(1).size());
        assertEquals(1, adjLabelMap.get(2).size());
        
        assertTrue(adjLabelMap.get(0).contains(1));
        assertTrue(adjLabelMap.get(1).contains(0));
        assertTrue(adjLabelMap.get(1).contains(2));
        assertTrue(adjLabelMap.get(2).contains(1));
    }
 
    public void testCondenseLabels() {
        
        int n = 5;
        
        int[] labels = new int[n];
        
        for (int i = 0; i < n; ++i) {
            labels[i] = i + n;
        }
        
        LabelToColorHelper.condenseLabels(labels);
        
        for (int i = 0; i < n; ++i) {
            assertEquals(i, labels[i]);
        }
    }
  
    public void testCreateAdjacencySetMap() {
        
        /*
          0 0  1  1  0  0 
          0 0  1  1  0  0
        */
        
        int w = 6;
        int h = 2;
        
        Image img = new Image(w, h);
        
        List<Set<PairInt>> list = new ArrayList<Set<PairInt>>();
        Set<PairInt> set0 = new HashSet<PairInt>();
        Set<PairInt> set1 = new HashSet<PairInt>();
        Set<PairInt> set2 = new HashSet<PairInt>();
        list.add(set0); list.add(set1); list.add(set2);
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int label = 0;
                if (i == 0 || i == 1) {
                    set0.add(new PairInt(i, j));
                } else if (i == 2 || i == 3) {
                    label = 1;
                    set1.add(new PairInt(i, j));
                } else {
                    label = 2;
                    set2.add(new PairInt(i, j));
                }
                img.setRGB(i, j, label, label, label);
            }
        }
        
        TIntObjectMap<TIntSet> labelToLabelAdjMap = 
            LabelToColorHelper.createAdjacencySetMap(list);
    
        assertEquals(3, labelToLabelAdjMap.size());
        assertTrue(labelToLabelAdjMap.get(0).contains(1));
        assertTrue(labelToLabelAdjMap.get(1).contains(0));
        assertTrue(labelToLabelAdjMap.get(1).contains(2));
        assertTrue(labelToLabelAdjMap.get(2).contains(1));
    }
}
