package algorithms.imageProcessing;

import java.lang.reflect.Type;
import java.lang.reflect.TypeVariable;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class FixedSizeSortedVectorTest extends TestCase {
    
    public void testAdd() throws Exception {
        
        FixedSizeSortedVector<Integer> sortedVector = 
            new FixedSizeSortedVector(4);
        
        sortedVector.add(Integer.valueOf(7));
        sortedVector.add(Integer.valueOf(6));
        sortedVector.add(Integer.valueOf(5));
        sortedVector.add(Integer.valueOf(4));
        
        Integer[] values = sortedVector.getArray();
        
        assertNotNull(values);
        
        assertTrue(values[0].intValue() == 4);
        assertTrue(values[1].intValue() == 5);
        assertTrue(values[2].intValue() == 6);
        assertTrue(values[3].intValue() == 7);
        
        sortedVector.add(Integer.valueOf(10));

        values = sortedVector.getArray();
        
        assertNotNull(values);

        assertTrue(values[0].intValue() == 4);
        assertTrue(values[1].intValue() == 5);
        assertTrue(values[2].intValue() == 6);
        assertTrue(values[3].intValue() == 7);
        
        sortedVector.add(Integer.valueOf(3));
        
        values = sortedVector.getArray();
        
        assertNotNull(values);

        assertTrue(values[0].intValue() == 3);
        assertTrue(values[1].intValue() == 4);
        assertTrue(values[2].intValue() == 5);
        assertTrue(values[3].intValue() == 6);
        
    }
    
    public void testTypes() throws Exception {
        
        FixedSizeSortedVector<A> vector = new FixedSizeSortedVector<A>(2);
        
        
        B b = new B(1);
        A a = new A(1);
        C c = new C(0);
        
        vector.add(a);
        vector.add(c);
        vector.add(c);
        
        A[] values = vector.getArray();
        
        assertNotNull(values);        
    }
    
    private static class A implements Comparable {
        int v = 0;
        public A(int v) {
            this.v = v;
        }
        @Override
        public int compareTo(Object o) {
            if (o == null) {
                return -1;
            }
            if (!(o instanceof A)) {
                return -1;
            }
            return Integer.compare(v, ((A)o).v);
        }
    }
    
    private static class B extends A {
        public B(int v) {
            super(v);
        }
    }
    
    private static class C extends A {
        public C(int v) {
            super(v);
        }
    }
}
