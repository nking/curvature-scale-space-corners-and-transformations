package algorithms;

/**
 * Methods rotate the values in an array of numbers by a distance in terms of
 * indexes.
 * 
 * @author nichole
 */
public class Rotate {
    
    /**
     * The solution is O(N) and uses the pattern flip top sub-array, flip bottom
     * sub-array, then flip all of the array.
     * The problem and suggested solution are from "Programming Pearls", 
     * Chapter 2, Problem B.
     * 
     * @param a array of numbers treated as a circular array.
     * @param left the number of spaces for which to shift left the values 
     * within array a. 
     */
    public void rotate(int[] a, int left) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if ((left == 0) || (left == a.length)) {
            return;
        }
                
        int n = a.length;
        
        if (left > 0) {
        
            if (left > n) {
                left = left % n;
            }

            reverse(a, 0, left - 1);
            reverse(a, left, n - 1);
            reverse(a, 0, n - 1);
            
        } else {
            
            left *= -1;
            if (left > n) {
                left = left % n;
            }
            
            reverse(a, n - left, n - 1);
            reverse(a, 0, left);
            reverse(a, 0, n - 1);
        }

    }
    
    public void rotate2(int[] a, int left) {
        rotate2(a, a.length, left);
    }
    
    /**
     * The solution is O(N) and uses the pattern of moving the first element
     * out of the array and moving the subsequent shiftee's into forward
     * shifted positions.  The algorithm has many more lines than the 
     * rotate(int[], int) method, but it has fewer iterations.  
     * The runtime is at most O(N).
     * 
     * Note that if left is a negative value, a reverse array before and
     * a reverse array afterwards are added.  For the current implementation,
     * it is not better than rotate(int[], int) method when left is
     * negative.
     * TODO: the method could be improved
     * for the case where left is a negative by writing different code for it, 
     * that is, code edited for left boundary logic and right shifts. 
     * 
     * The problem and suggested solution are from "Programming Pearls", 
     * Chapter 2, Problem B.
     * 
     * @param a array of numbers treated as a circular array.
     * @param n the number of items in the array to sort (i.e. ignoring beyond index n-1).
     * @param left the number of spaces for which to shift left the values 
     * within array a.
     */
    public void rotate2(int[] a, int n, int left) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if ((left == 0) || (left == a.length)) {
            return;
        }
                
        boolean leftIsNegative = (left < 0);
        
        if (leftIsNegative) {
            reverse(a, 0, n - 1);
            left *= -1;
        }        
        
        if (left > n) {
            left = left % n;
        }
                
        int firstIdx = 0;
        int lastIdx = n - 1;
        
        int offset = 0;
        int count = 0;
        int tmp = a[offset];
        
        boolean tmpHoldsFirstValue = true;
        boolean tmpHoldsLastValue = false;
        
        int nIter = 0;
        int prevResetOffset = 0;
        int prevResetNIter = 0;
                
        int idx, idx0;
        while (true) {
            
            if (nIter > n) {
                // this shouldn't happen!
                throw new IllegalStateException("the algorithm has an error.");
            }
            
            count++;
            idx = (left*count) + offset;
            idx0 = idx - left;
            
            // check for conditions to exit loop or change loop parameters:
            if (nIter == (n - left)) {
                
                a[idx0] = tmp;
                
                if (tmpHoldsFirstValue) {
                    firstIdx = idx0;
                    tmpHoldsFirstValue = false;
                } else if (tmpHoldsLastValue) {
                    lastIdx = idx0;
                    tmpHoldsLastValue = false;
                }
                
                if ((nIter == (n - 1)) || (nIter == firstIdx)){
                    
                    break;
                    
                } else {
                    
                    // recalc offset
                    if (prevResetOffset == 0) {
                        offset = lastIdx + 1;
                    } else {
                        offset = prevResetOffset + (nIter - prevResetNIter);
                    }
                    prevResetOffset = offset;
                    prevResetNIter = nIter;

                    // recalc left
                    left = firstIdx - offset;
                    
                    firstIdx = offset;
                    
                }
                
                count = 0;
                tmp = a[offset];
                
                if (offset == firstIdx) {
                    tmpHoldsFirstValue = true;
                } else if (offset == lastIdx) {
                    tmpHoldsLastValue = true;
                }
                
                continue;
            }

            if ((idx < n) && (idx > -1)) {
                
                a[idx0] = a[idx];
                
                if (idx == firstIdx) {
                    firstIdx = idx0;
                } else if (idx == lastIdx) {
                    lastIdx = idx0;
                }
                
                nIter++;
                
            } else {
                // idx has overrun a bounds
                
                a[idx0] = tmp;
                
                if (tmpHoldsFirstValue) {
                    firstIdx = idx0;
                    tmpHoldsFirstValue = false;
                } else if (tmpHoldsLastValue) {
                    lastIdx = idx0;
                    tmpHoldsLastValue = false;
                }
                
                offset++;
                count = 0;
                tmp = a[offset];
                
                if (offset == firstIdx) {
                    tmpHoldsFirstValue = true;
                } else if (offset == lastIdx) {
                    tmpHoldsLastValue = true;
                }                
            }
        }
                
        //System.out.println("    nIter=" + nIter);
        
        if (leftIsNegative) {
            reverse(a, 0, n - 1);
            //left *= -1;
        }
    }
  
    /**
     * reverse the array between indexes idxLo and idxHi, inclusive.
     * 
     * @param a array of numbers
     * @param idxLo the smallest index of the range to reverse in array a
     * @param idxHi the largest index, inclusive, of the range to reverse in 
     * array a
     */
    public void reverse(int[] a, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        int n = a.length;
        if ((idxLo < 0) || (idxLo > (n - 1))) {
            throw new IllegalArgumentException("idxLo is out of bounds of array");
        }
        if ((idxHi < 0) || (idxHi > (n - 1))) {
            throw new IllegalArgumentException("idxHi is out of bounds of array");
        }
        
        n = idxHi - idxLo + 1;
        
        int end = idxLo + (n/2);
        
        int count = 0;
        for (int i = idxLo; i < end; i++) {
            int idx2 = idxHi - count;
            int swap = a[i];
            a[i] = a[idx2];
            a[idx2] = swap;
            count++;
        }
    }
}
