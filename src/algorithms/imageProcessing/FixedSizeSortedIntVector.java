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
    
    public int getCapacity() {
        return size;
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
    
    public int findIndexFor(int value) {
        int idx = Arrays.binarySearch(a, 0, n, value);
        if (idx < 0) {
            //(-(insertion point) - 1)
            idx *= -1;
            idx--;
        } else {
            // in case there are more than one of same value, search prev
            while ((idx - 1) > -1 && a[idx - 1] == value) {
                idx--;
            }
        }
        return idx;
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

    /**
     * Split this into a new instance containing all elements less than or
     * equal to x and leave all elements greater than x in this instance
     * runtime complexity is binary search + copy = O(n) where n is the
     * number of items in this instance.
     * @param x
     * @return
     */
    public FixedSizeSortedIntVector split(int x) {

        int idx = findIndexFor(x);

        // copy: 0 : 1dx-1;  keep: idx : n-1
        FixedSizeSortedIntVector vecL = new FixedSizeSortedIntVector(size);
        int[] aL = vecL.getArray();

        // copy 0:idx-1 to aL
        System.arraycopy(a, 0, aL, 0, idx);
        vecL.n = (idx + 1);
        vecL.availSlot = vecL.n;

        // move idx:n-1 up to 0 
        System.arraycopy(a, idx, a, 0, n - idx);
        n = n - idx;
        availSlot = n;
        
        return vecL;
    }
}
