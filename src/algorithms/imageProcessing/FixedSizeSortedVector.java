package algorithms.imageProcessing;

import java.lang.reflect.Array;
import java.util.Arrays;

/**
 * a class to hold a fixed number of items of type T which are kept in a sorted
 * stated after the list is filled, though a getArray before it is filled
 * will still return sorted results.
 * 
 * runtime complexity is:
 *     O(k) + O(N*lg_2(k)) 
 * where k is the fixedCapacity and N is the number of times add is used.
 * 
 * @author nichole
 */
@SuppressWarnings({"unchecked"})
public class FixedSizeSortedVector<T extends Comparable<T>> {
    
    protected T[] a = null;
        
    protected final int size;
        
    protected int n;
        
    protected int availSlot;
        
    protected boolean sorted;
    
    public FixedSizeSortedVector(int fixedCapacity, Class<T> classTypeToHold) {
        
        size = fixedCapacity;
            
        n = 0;

        availSlot = -1;

        sorted = false;
        
        a = (T[]) Array.newInstance(classTypeToHold, size);
        
    }
        
    /**
     * add value to the fixed size sorted list if the list is not full
     * or if the last item in the list is valued as worse than given value
     * (where (T).compareTo determines the value for the descending sort
     * of this class instance).
     * 
     * @param value 
     */
    public void add(T value) {

        if (value == null) {
            return;
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

            // insert value into array at position found by binarySearch
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

    /**
     * get the internal array for the sorted list.  note this is not a copy in
     * order to keep the use small, so do not edit it and continue to use
     * the add method.
     * 
     * runtime complexity is O(1) unless the internal list is not yet
     * filled, then the sort makes the complexity O(n*lg2(n)).
     * 
     * @return 
     */
    public T[] getArray() {

        if (!sorted && (a != null)) {
            
            // the sort will fail for nulls, so need to make a copy, sort it
            // and copy it back.  defeats the whole purpose of late construction of a!
            T[] c = Arrays.copyOfRange(a, 0, n);
            
            Arrays.sort(c);
            
            System.arraycopy(c, 0, a, 0, n);
            
            // do not set the 'sorted' flag because the list is not full yet
        }
        
        return a;
    }
}
