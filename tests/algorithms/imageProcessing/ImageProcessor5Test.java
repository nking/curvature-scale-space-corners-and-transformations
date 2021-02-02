package algorithms.imageProcessing;

import algorithms.VeryLongBitString;
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
    }
    
}
