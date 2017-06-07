package algorithms.misc;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class MiscTest extends TestCase {
    
    public MiscTest() {
    }

    public void testPersistToFile() throws Exception {
        
        Set<PairInt> points = new HashSet<PairInt>();
        for (int i = 0; i < 10; ++i) {
            points.add(new PairInt(i, i));
        }
        
        Set<PairInt> points2 = new HashSet<PairInt>(points);
        
        String fileName = "blob_100_150.dat";
        
        String outPath = Misc.persistToFile(fileName, points);
        
        assertTrue((new File(outPath)).exists());
        
        String filePath = ResourceFinder.findDirectory("bin") + "/" + fileName;
        
        Set<PairInt> pointsR = Misc.deserializeSetPairInt(filePath);

        assertTrue(pointsR.size() == points.size());
        
        for (PairInt p : pointsR) {
            assertTrue(points.remove(p));
        }
        
        assertTrue(points.isEmpty());
        
        PairIntArray pointsR2 = Misc.deserializePairIntArray(filePath);

        assertTrue(pointsR2.getN() == points2.size());
        
        for (int i = 0; i < pointsR2.getN(); ++i) {
            PairInt p2 = new PairInt(pointsR2.getX(i), pointsR2.getY(i));
            assertTrue(points2.remove(p2));
        }
        
        assertTrue(points2.isEmpty());
    }    
    
    public void testReverse() {
        
        List<Integer> list = new ArrayList<Integer>();
        for (int i = 0; i < 4; ++i) {
            Integer index = Integer.valueOf(i);
            list.add(index);
        }
        Misc.<Integer>reverse(list);
        
        for (int i = 0; i < 4; ++i) {
            Integer index = list.get(i);
            assertEquals(3 - i, index.intValue());
        }
    }
   
    public void testHashCode() {
        
        int length = 100000;//256;
        
        int n = 1000000;
        
        TIntObjectMap<TIntSet> idxVMap = new TIntObjectHashMap<TIntSet>(length);
                
        Random rng = Misc.getSecureRandom();
        
        int nCollisions = 0;
        
        for (int i = 0; i < n; ++i) {
            int v = rng.nextInt(length);
            int idx = MiscMath.hashToIndex(v, length);
            assertTrue(idx < length);
            assertTrue(idx >= 0);
            if (idxVMap.containsKey(idx) && idxVMap.get(idx).contains(v)) {
                nCollisions++;
            }
            TIntSet set = idxVMap.get(idx);
            if (set == null) {
                set = new TIntHashSet();
                idxVMap.put(idx, set);
            }
            set.add(v);
        }
                
        System.out.println("nCollisons=" + nCollisions 
            + " idxs.size=" + idxVMap.size() + " nC/N=" 
            + ((double)nCollisions/n));
        
        assertTrue(idxVMap.size() > 0.9*length);  
        
        int len = 10;
        int[] hash = new int[len];
        int[] list = new int[len];
        Arrays.fill(hash, -1);
        Arrays.fill(list, -1);
    
        for (int i = 0; i < len; ++i) {
            list[i] = i;
            MiscMath.setHashValue(list, i, hash);
            assertTrue(MiscMath.hashContains(list, i, hash));
        }
        
    }
    
    public void testPermuteTheSetBits() {
         
        // 1 1 0 1
        int setBits = (1 << 0) + (1 << 2) + (1 << 3);
        
        int[] expected = new int[8];
        expected[1] = (1 << 0) + (1 << 2) + (1 << 3);//1 1 0 1
        expected[2] = (1 << 2) + (1 << 3);//1 1 0 0
        expected[3] = (1 << 0) + (1 << 3);//1 0 0 1
        expected[4] = (1 << 3);//1 0 0 0
        expected[5] = (1 << 0) + (1 << 2);//0 1 0 1
        expected[6] = (1 << 2);//0 1 0 0
        expected[7] = (1 << 0);//0 0 0 1
        
        Arrays.sort(expected);
   
        int[] list = MiscMath.permuteTheSetBits(setBits);

        Arrays.sort(list);
        
        assertTrue(Arrays.equals(expected, list));
    }
    
}
