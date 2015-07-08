package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import java.lang.reflect.Array;
import java.lang.reflect.ParameterizedType;
import java.lang.reflect.Type;
import java.lang.reflect.TypeVariable;
import java.util.Arrays;

/**
 * O(k) + O(N*lg_2(k)) 
 * @author nichole
 */
public class FixedSizeSortedVector<T extends Comparable> {
    
    protected T[] a = null;
        
    protected final int size;
        
    protected int n;
        
    protected int availSlot;
        
    protected boolean sorted;
    
    public FixedSizeSortedVector(int fixedCapacity) {
        
        size = fixedCapacity;
            
        n = 0;

        availSlot = -1;

        sorted = false;
    }
        
    private void initializeA(T value) {

        if (a == null) {

            Class<?> cls = value.getClass();

            Type[] types = cls.getGenericInterfaces();

            if (types != null && types.length == 1) {
                ParameterizedType pt = (ParameterizedType)types[0];
                Type[] ts = pt.getActualTypeArguments();
                if (ts != null && ts.length == 1) {
                    cls = (Class<?>)ts[0];
                }
            }
           
            a = (T[]) Array.newInstance(cls, size);
        }
    }
  
    public void add(T value) {

        if (value == null) {
            return;
        }

        if (a == null) {
            initializeA(value);
        }

        if (!sorted) {

            assert(n < a.length);

            append(value);

        } else {

            assert(n == a.length);
               
            // compare to last item in array and return if this is a worse fit
            int comp = value.compareTo(a[n - 1]);

            if (comp != -1) {
                return;
            }

            // free up the last slot
            availSlot = n - 1;

            n--;

            // insert value into array a using position found by binarySearch
            // and compareTo
            insertIntoOpenSlot(value);
        }
    }
        
    /**
     * append item value onto the end of the list.  Note that if the item
     * is added to the last slot, the list is immediately sorted into 
     * ascending numerical order
     * afterwards as a side effect to keep the logic in the other
     * methods consistent.
     * runtime is usually O(1), but if append is used for the last item,
     * there is a sort adding O(N*log_2(N)).
     * For best use, append(v) the first size-1 items and thereafter use
     * insertIntoOpenSlot(v).
     * 
     * @param value 
     */
    private void append(T value) {

        if (n == (a.length)) {
            throw new IllegalArgumentException(
                "there must be an empty slot in order to append." +
                " remove and item then try insert again or construct larger list.");
        }

        a[n] = value;

        n++;

        if (n == a.length) {

            Arrays.sort(a);

            sorted = true;
        }
    }

    /**
     * Insert the value into the list while maintaining the sorted state
     * of the list.  Note that if there is not exactly one available slot
     * in the list, an IllegalArgumentException will be thrown.
     * runtime is usually O(log_2(N)) + less than O(N), but once per class lifetime
     * the sort may occur here adding O(N*log_2(N)).
     * @param value 
     */
    private void insertIntoOpenSlot(T value) {

        if (n != (a.length - 1)) {
            throw new IllegalArgumentException(
            "the method is meant to be used when there is only one avail slot");
        }

        int insIdx = Arrays.binarySearch(a, value);
        if (insIdx < 0) {
            insIdx *= -1;
            insIdx--;
        }

        if (insIdx == availSlot) {

            a[availSlot] = value;

        } else if (insIdx < availSlot) {

            // move all items from insIdx to availSlot down by 1
            for (int i = (availSlot - 1); i >= insIdx; i--) {
                a[i + 1] = a[i];
            }

            a[insIdx] = value;

        } else {

            int end = insIdx - 1;

            // move items up from availSlot +1 to insIdx - 1
            // then insert value into insIdx - 1
            for (int i = availSlot; i < end; i++) {
                a[i] = a[i + 1];
            }

            a[insIdx - 1] = value;
        }

        n++;

        availSlot = -1;           
    }

    public T[] getArray() {

        return a;
    }
}
