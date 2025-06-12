package algorithms.disjointSets;

import junit.framework.TestCase;

import java.util.Map;
import java.util.Set;

import static org.junit.Assert.*;

public class UnionFindExtTest extends TestCase {

    public void test0() {
        int n = 4;
        UnionFind uf = new UnionFind(n);
        uf.union(1, 3);

        UnionFindExt uf2 = new UnionFindExt(uf);

        Map<Integer, Set<Integer>> map1 = uf.getComponents();
        Map<Integer, Set<Integer>> map2 = uf2.getComponents();
        assertEquals(3, map1.size());
        assertEquals(map1.size(), map2.size());
        for (Map.Entry<Integer, Set<Integer>> entry : map1.entrySet()) {
            assertTrue(map2.containsKey(entry.getKey()));
            for (int v : entry.getValue()) {
                assertTrue(map2.get(entry.getKey()).contains(v));
            }
            assertEquals(entry.getValue().size(), map2.get(entry.getKey()).size());
        }
    }
}