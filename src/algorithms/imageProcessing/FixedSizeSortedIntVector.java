package algorithms.imageProcessing;

import java.util.Arrays;

/**
 * a class to hold a fixed number of items of type int which are kept in a sorted
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
 */
public class FixedSizeSortedIntVector {

    protected int[] a = null;

    protected final int size;

    protected int n;

    protected int availSlot;

    public FixedSizeSortedIntVector(int fixedCapacity) {
        
        if (fixedCapacity < 1) {
            throw new IllegalArgumentException(
            "fixedCapacity must be a positive non zero (arg was " 
            + fixedCapacity + ")");
        }

        size = fixedCapacity;

        n = 0;

        availSlot = -1;

        a = new int[size];

    }

    /**
     * add value to the fixed size sorted list (sorted by increasing value).
     * 
     * runtime complexity is O(log_2(capacity) + less than capacity).
     *
     * @param value
     * @return true if added, else false
     */
    public boolean add(int value) {

        if (n < size) {

            if (availSlot == -1) {
                availSlot = n;
            }
            
            insertIntoOpenSlot(value);
            
        } else {
            
            int compareIdx = n - 1;
            
            if ((n == 1) && (size == 1)) {
                compareIdx = 0;
            }
            
            int comp = (value < a[compareIdx]) ? -1 : 
                ((value >a[compareIdx]) ? 1 : 0);
            
            if (comp != -1) {
                return false;
            }
            
            // free up the last slot
            availSlot = compareIdx;

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
    private void insertIntoOpenSlot(int value) {

        int insIdx = Arrays.binarySearch(a, 0, n, value);
        if (insIdx < 0) {
            insIdx *= -1;
            insIdx--;
        }

        if ((availSlot > -1) && (insIdx > availSlot) && (a[availSlot] == value)) {
            // this depends upon logic of previous remove setting availSlot 
            // to next value.
            boolean b = true;
            for (int i = insIdx; i > availSlot; --i) {
                if (a[i] != value) {
                    b = false;
                    break;
                }
            }
            
            if (b) {
            
                // no need to set value again
                n++;

                availSlot = -1;

                return;
            }
        }
        
        if (insIdx == availSlot) {

            a[availSlot] = value;
            
        } else if ((insIdx == (a.length - 1)) && (availSlot == (insIdx - 1))) {

            a[insIdx] = value;
            
        } else if (insIdx < availSlot) {

            // move all items from insIdx to availSlot down by 1
            for (int i = (availSlot - 1); i >= insIdx; i--) {
                a[i + 1] = a[i];
            }

            a[insIdx] = value;

        } else {

            int end = insIdx - 1;
           
            if (availSlot > -1) {
                while ((a[insIdx] == value) && ((end + 1) <= n) && (a[end + 1] == value)) {
                    end++;
                }
            }

            // move items up from availSlot +1 to insIdx - 1
            // then insert value into insIdx - 1
            for (int i = availSlot; i < end; i++) {
                a[i] = a[i + 1];
            }

            a[insIdx] = value;
        }

        n++;

        availSlot = -1;
    }

    /**
     * remove the item from the full list of items.
     * runtime is O(log_2(N)).
     * NOTE: this could be made O(1) runtime complexity at the expense
     * of 3 * space complexity.
     * @param value
     */
    public void remove(int value) {

        if (n != a.length) {
            throw new IllegalArgumentException(
            "the method is meant to be used only on a full list");
        }

        int rmIdx = Arrays.binarySearch(a, value);

        if (rmIdx < 0) {
            throw new IllegalArgumentException("could not find item " + value + 
                " in list");
        } else {
            // search for last in list with that value
            for (int i = rmIdx; i < n; ++i) {
                if (a[i] == value) {
                    rmIdx = i;
                } else {
                    break;
                }
            }
        }

        availSlot = rmIdx;

        // to keep the list in a state where the next binary search works,
        // set the empty slot value to the proceeding value or max integer.
        if (availSlot == (a.length - 1)) {
            a[availSlot] = Integer.MAX_VALUE;
        } else {
            a[availSlot] = a[availSlot + 1];
        }

        n--;
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
    public int[] getArray() {

        return a;
    }
    
    /**
     * get the internal array for the sorted list.  note this is not a copy in
     * order to keep the use small, so do not edit it and continue to use
     * the add method.
     *
     * runtime complexity is O(1)
     *
     * @param index
     * @return
     */
    public int getValue(int index) {

        if ((index < 0) || (index > (a.length - 1))) {
            throw new IllegalArgumentException("index is out of bounds");
        }
        
        return a[index];
    }
    
    public int getMedianValue() {

        int midIdx = ((n & 1) == 1) ? n / 2 : (n - 1) / 2;

        return a[midIdx];
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

    @Override
    public String toString() {
        return Arrays.toString(a);
    }
    
}
