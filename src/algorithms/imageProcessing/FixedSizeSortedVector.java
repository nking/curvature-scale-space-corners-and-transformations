package algorithms.imageProcessing;

import java.lang.reflect.Array;
import java.util.Arrays;

/**
 * a class to hold a fixed number of items of type T which are kept in a sorted
 * stated.
 *<pre>
 * runtime complexity is:
 *     O(N * (lg_2(k) + smaller than k))
 * where k is the fixed capacity and N is the number of times add is used.
 * 
 * worse case runtime complexity is O(N * (k + lg_2(k)))
 * best case runtime complexity is O(N * (1 + lg_2(k)))
 *</pre>
 * @author nichole
 * @param <T> class type to be held and sorted by this class.  It must implement
 * Comparable.
 */
@SuppressWarnings({"unchecked"})
public class FixedSizeSortedVector<T extends Comparable<T>> {

    protected T[] a = null;

    protected final int size;

    protected int n;

    protected int availSlot;

    public FixedSizeSortedVector(int fixedCapacity, Class<T> classTypeToHold) {

        size = fixedCapacity;

        n = 0;

        availSlot = -1;

        a = (T[]) Array.newInstance(classTypeToHold, size);

    }

    /**
     * add value to the fixed size sorted list using (T).compareTo to order
     * the items in the internal list.
     * 
     * runtime complexity is O(log_2(capacity) + less than capacity).
     *
     * @param value
     * @return true if added, else false
     */
    public boolean add(T value) {

        if (value == null) {
            return false;
        }

        if (n < a.length) {

            availSlot = n;
            
            insertIntoOpenSlot(value);
            
        } else {

            // compare to last item in array and return if this is a worse fit
            int comp = value.compareTo(a[n - 1]);

            if (comp != -1) {
                return false;
            }

            // free up the last slot
            availSlot = n - 1;

            n--;

            // insert value into array at position found by binarySearch
            insertIntoOpenSlot(value);
        }
        
        return true;
    }

    /**
     * Insert the value into the list while maintaining the sorted state
     * of the list.  
     * @param value
     */
    private void insertIntoOpenSlot(T value) {

        int insIdx = Arrays.binarySearch(a, 0, n, value);
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

    /**
     * get the internal array for the sorted list.  note this is not a copy in
     * order to keep the use small, so do not edit it and continue to use
     * the add method.
     *
     * runtime complexity is O(1)
     *
     * @return
     */
    public T[] getArray() {

        return a;
    }
    
    /**
     * return the number of items in the internal array.  if the array is not
     * yet filled, the return will be less than the capacity, else will
     * be the same as the capacity.
     * @return 
     */
    public int getNumberOfItems() {
        return n;
    }
}
