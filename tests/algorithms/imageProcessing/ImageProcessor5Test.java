package algorithms.imageProcessing;

import algorithms.util.VeryLongBitString;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ImageProcessor5Test extends TestCase {

    public ImageProcessor5Test(String testName) {
        super(testName);
    }
    
    public void test0() throws Exception {
        
        /*
        testing adjacency map and point index maps
        
        making a grid of color cells that are 8 by 8 pixels and
           8 cells on each side with increasing color
       
        */
        int w = 64;
        int h = w;
        int nCells = 8;
        int cellSz = w/nCells;
        
        ImageProcessor imageProcessor = new ImageProcessor();

        ImageExt img = new ImageExt(w, h);
        
        TIntObjectMap<TIntSet> mapOfSets = new TIntObjectHashMap<TIntSet>();
        
        TIntIntMap pointIndexMap = new TIntIntHashMap();
        
        /*
        
        64
        
        56
        
        48
        
        40
        
        32
            24  25   26   27   28   29   30  31
        24
            16  17   18   19   20   21   22  23
        16
            8   9    10   11   12   13   14  15
        8
            0   1    2    3    4    5    6   7
        0   
          0   8   16   24   32   40   48  56  64
        */
        
        for (int j = 0; j < h; ++j) {
            int g = (j/cellSz);
            for (int i = 0; i < w; ++i) {
                int r = (i/cellSz);
                int b = r;
                
                img.setRGB(i, j, r * 8, g * 8, b * 8);
                
                int label = (g * nCells) + r;
                int pixIdx = (j * w) + i;
                //System.out.format("(%d,%d) label=%d pixIdx=%d\n", i, j, label, 
                //    pixIdx);
            
                TIntSet set = mapOfSets.get(label);
                if (set == null) {
                    set = new TIntHashSet();
                    mapOfSets.put(label, set);
                }
                set.add(pixIdx);
                
                pointIndexMap.put(pixIdx, label);
            }
        }

        TIntObjectMap<VeryLongBitString> adjMap = 
            imageProcessor.createAdjacencyMap(pointIndexMap, 
                mapOfSets, w, h);
        
        TIntObjectIterator<VeryLongBitString> iter = adjMap.iterator();
        
        for (int i = 0; i < adjMap.size(); ++i) {
            iter.advance();
            
            int label = iter.key();
            VeryLongBitString adj = iter.value();
            int[] setBits = adj.getSetBits();
            
            //System.out.format("label=%d adj=%s\n", label, 
            //    Arrays.toString(setBits));
        }
        
        /*
        
        64
        
        56
        
        48
        
        40
        
        32
            24  25   26   27   28   29   30  31
        24
            16  17   18   19   20   21   22  23
        16
            8   9    10   11   12   13   14  15
        8
            0   1    2    3    4    5    6   7
        0   
          0   8   16   24   32   40   48  56  64
        */
        
        assertEquals(0, pointIndexMap.get((1 * w) + 1));
        assertEquals(3, pointIndexMap.get((4 * w) + 25));
        assertEquals(27, pointIndexMap.get((25 * w) + 25));
    
        int[] expected = new int[]{1, 8};
        int[] found = adjMap.get(0).getSetBits();
        //System.out.println("found=" + Arrays.toString(found));
        assertTrue(Arrays.equals(expected, found));
        
        assertTrue(Arrays.equals(new int[]{6, 15}, 
            adjMap.get(7).getSetBits()));
        
        assertTrue(Arrays.equals(new int[]{11, 18, 20, 27}, 
            adjMap.get(19).getSetBits()));
    }
    
}
