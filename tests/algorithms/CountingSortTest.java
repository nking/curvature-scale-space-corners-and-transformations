package algorithms;

import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import junit.framework.TestCase;

public class CountingSortTest extends TestCase {

    public void testsort() throws Exception {
        
        int[] a = new int[]{2, 5, 3, 0, 2, 3, 0, 3};
        
        int[] b = CountingSort.sort(a, 5);
        
        int[] expected = new int[]{0, 0, 2, 2, 3, 3, 3, 5};
        
        assertTrue(Arrays.equals(expected, b));
    }
    
    public void testsort1() throws Exception {
        
        int[] a = new int[]{2, 5, 3, 0, 2, 3, 0, 3};
        int[] b = new int[]{2, 5, 3, 0, 2, 3, 0, 3};
        
        CountingSort.sort(a, b, 5);
        
        int[] expected = new int[]{0, 0, 2, 2, 3, 3, 3, 5};
        
        assertTrue(Arrays.equals(expected, a));
        assertTrue(Arrays.equals(expected, b));
    }
    
    public void testsort2() throws Exception {
        
        int[] a = new int[]{2, 5, 3, 0, 2, 3, 0, 3};
        int[] b = new int[]{2, 5, 3, 0, 2, 3, 0, 3};
        
        CountingSort.sortByDecr(a, b, 5);
        
        int[] expected = new int[]{5, 3, 3, 3, 2, 2, 0, 0};
        
        assertTrue(Arrays.equals(expected, a));
        assertTrue(Arrays.equals(expected, b));
    }
    
    public void testsort3() throws Exception {
        
        // use more than 46340 random numbers whose value is higher than
        // 46340 to show that the internal long summations are safely
        // reduced back to the integer values
        
        List<Integer> list = new ArrayList<Integer>();
        
        SecureRandom sr = new SecureRandom();
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        int n = 46340*2;
        
        int[] a = new int[n];
        
        int max = Integer.MIN_VALUE;
        
        for (int i = 0; i < n; i++) {
            int r = sr.nextInt(67_000);
            while (r < 0) {
                r = sr.nextInt();
            }
            list.add(Integer.valueOf(r));
            a[i] = r;
            
            if (r > max) {
                max = r;
            }
        }
        
        int[] b = CountingSort.sort(a, max);
        
        Collections.sort(list);
        
        for (int i = 0; i < n; i++) {
            assertTrue(list.get(i).intValue() == b[i]);
        }
    }

}
