package algorithms.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
//import java.util.logging.Level;
//import java.util.logging.Logger;

/**
 * NOTE: still testing this class.  it needs more testing for operations that
 * change the instance value from positive to negative and for cases when
 * a divisor is still larger than 64 bits.
 * 
 * A class to hold numbers that can be larger than 64 bits after adds.
 * 
 * The class was created w/ Pascal's triangle as a use case in order to be able
 * to print out numbers larger than ((1<<63) -1) for level >= 64.
 * 
 * the core of the constructor and add methods were started with
 *     http://www.cs.utexas.edu/users/djimenez/utsa/cs3343/lecture20.html
   and then added to and adapted afterwards.
 
 * the core of the Euclidean division started with:
 *     https://en.wikipedia.org/wiki/Division_algorithm
 *  
 * @author nichole
 */
public class VeryLargeNumber implements Comparable<VeryLargeNumber>, Cloneable {
    
    // could consider using a larger base for compaction, 
    // and making use of the space not used by positive values,
    // but for now, the simplicity is more important for testing and first uses.
    // adding faster alternative methods for the large integer division
    // could improve the largest bottleneck.
 
    public final static int BASE = (1<<30)-1;
    
    private int[] a;
    
    private int nLen = 0;
    
    private boolean isPositive = true;
    
    //private Logger log = Logger.getLogger(this.getClass().getName());
    
    public VeryLargeNumber(int number) {
        
        a = new int[10];
        
        createInternalNumber(number);
    }
    
    /**
     * add addThis to this instance.
     * 
     * @param addThis number to add to this
     */
    public void add(VeryLargeNumber addThis) {
                
        //log.log(Level.FINEST, "add: {0} + {1}", new String[]{toString(), 
        //    addThis.toString()});
        
        boolean thisIsLarger = (addThis.nLen <= nLen);
        
        int[] longer;
        int[] shorter;
        int n, idxOffsetShorter;
        boolean longerIsPositive, shorterIsPositive;
        
        if (thisIsLarger) {
            longer = a;
            shorter = addThis.a;
            n = nLen;
            idxOffsetShorter = nLen - addThis.nLen;
            longerIsPositive = isPositive;
            shorterIsPositive = addThis.isPositive();
        } else {
            longer = addThis.a;
            shorter = a;
            n = addThis.nLen;
            idxOffsetShorter = addThis.nLen - nLen;
            longerIsPositive = addThis.isPositive();
            shorterIsPositive = isPositive;
        }
                
        int	carry = 0;
        int sum = 0;
        int idxShorter;
        boolean resetCarry = true;
        
        int[] s = new int[n];
        
        for (int i = (n - 1); i > -1; i--) {
                   
            resetCarry = true;
            
            idxShorter = i - idxOffsetShorter;
           
            int tl = longerIsPositive ? longer[i] : -1*longer[i];
            
            if (idxShorter < 0) {
                
                if ((carry == -1) && (longer[i] == 0)) {
                    
                    carry = -1;
                    
                    resetCarry = false;
                    
                    sum = BASE - 1;
                    
                } else {
                    
                    sum = tl + carry;
                }
                
            } else {
                    
                int ts = shorterIsPositive ? shorter[idxShorter] : 
                    -1*shorter[idxShorter];
                    
                if (!shorterIsPositive && longerIsPositive && 
                    (shorter[idxShorter] > longer[i]) && (i > 0)) {
                    
                    sum = tl + BASE + ts + carry;
                    
                    carry = -1;
                    
                    resetCarry = false;
                    
                } else if ((longer[i] == 0) && !longerIsPositive) {
                        
                    //e.g. 4 - 10 
                    sum = ts - BASE;
                        
                    carry = 1;
                    
                    resetCarry = false;
                    
                } else {
                        
                    sum = tl + ts + carry;
                }
            }

            if ((sum < 0) && (sum < -1*BASE)) {
                
                sum += BASE;
                
                carry = -1;
                
            } else if (sum >= BASE) {

                carry = 1;

                sum -= BASE;
                
            } else {
                
                if (resetCarry) {
                    carry = 0;
                }
            }
            
            if (longerIsPositive && (sum < 0)) {
                longerIsPositive = false;
            } else if (!longerIsPositive && (sum > 0)) {
                longerIsPositive = true;
            }

            s[i] = longerIsPositive ? sum : -1*sum;
        }
        
        a = s;
        
        nLen = n;
        
        isPositive = longerIsPositive;
        
        if (carry != 0) {
            // expand a by 1 and move elements down
            insertSpaceAtTopOfArray();
            
            if (carry < 0) {
                isPositive = false;
                carry *= -1;
            }
            
            a[0] = carry;
            
        } else {
        
            moveUpIfStartsWithZeros();
        }
    }
    
    /**
     * if there are leading 0's in the array a, move items below
     * them up to fill them and reduce the size of nLen. if there are
     * no numbers left, the result is nLen=1 to result in a value of 0
     * for this instance's number.
     */
    private void moveUpIfStartsWithZeros() {
        
        // move up if needed
        int firstNonZeroIdx = -1;
        for (int i = 0; i < nLen; i++) {
            if (a[i] == 0) {
                firstNonZeroIdx = i + 1;
            } else {
                break;
            }
        }
        
        if (firstNonZeroIdx == nLen) {
            nLen = 1;
            isPositive = true;
            Arrays.fill(a, 0);
            return;
        }
        
        if (firstNonZeroIdx > -1) {
            
            for (int i = 0; i < (nLen - firstNonZeroIdx); i++) {
                a[i] = a[i + firstNonZeroIdx];
            }
            
            Arrays.fill(a, (nLen - firstNonZeroIdx), nLen, 0);
            
            nLen -= firstNonZeroIdx;
        }
    }
    
    /**
     * subtract subtractThis from this instance.
     * @param subtractThis 
     */
    public void subtract(VeryLargeNumber subtractThis) {
        
        //log.log(Level.FINEST, "subtract: {0} - {1}", new String[]{toString(), 
        //    subtractThis.toString()});
        
        subtractThis.reversePolarity();
        
        add(subtractThis);
        
        subtractThis.reversePolarity();
    }
 
    /**
     * divide internal number by the divisor and return a string.  note that the
     * method currently uses the simplest implementation, Euclidean division.
     * a faster internal implementation can be made upon need.
     * 
     * @param divisor
     * @return 
     */
    public String divideByAndPrint(VeryLargeNumber divisor) {
        
        //log.log(Level.FINEST, "divide: {0} / {1}", new String[]{toString(), 
        //    divisor.toString()});
               
        return divideByAndPrintEuclidean(divisor);
    }
    
    /**
     * whether this number value is 0
     * @return 
     */
    public boolean isZero() {
        long sum = 0;
        for (int i = 0; i < nLen; i++) {
            sum += a[i];
        }
        return (sum == 0);
    }
    
    /**
     * whether this is a positive number
     * @return 
     */
    public boolean isPositive() {
        return isPositive;
    }
    
    /**
     * reverse the sign of this number.
     */
    public void reversePolarity() {
        isPositive = !isPositive;
    }
    
    /**
     * using Euclidean division, divide this by divisor and return the result
     * as a string.  The string output is because the currently using code
     * needs only that.
     * 
     * @param divisor
     * @return 
     */
    private String divideByAndPrintEuclidean(VeryLargeNumber divisor) {
        
        if (divisor.isZero()) {
            throw new IllegalArgumentException("Cannot divide by zero");
        }
        
        boolean divisorIsNegative = !divisor.isPositive();
        
        boolean thisIsNegative = !isPositive;
        
        if (divisorIsNegative) {
            divisor.reversePolarity();
        }
        
        if (thisIsNegative) {
            reversePolarity();
        }
        
        // both this number and divisor are positive or zero
        VeryLargeNumber q = new VeryLargeNumber(0);
        
        VeryLargeNumber r = null;
        
        try {
            r = clone();
        } catch (CloneNotSupportedException e) {
            // this will never happen...
            throw new IllegalStateException("problem w/ native support for " +
                " cloneable? ", e);
        }

        // while  R ≥ D
        while (r.compareTo(divisor) > -1) {
            
            q.increment();
            
            r.subtract(divisor);
        }
        
        if (!r.isPositive()) {
            // the while loop proceeds one step too far so reverse by 1 loop
            q.reversePolarity();
            q.increment();
            q.reversePolarity();
            r.add(divisor);
        }
        
        if (thisIsNegative) {
            reversePolarity();
        }
        
        if (divisorIsNegative) {
            divisor.reversePolarity();            
        }
        
        if (divisorIsNegative && !thisIsNegative) {
            q.reversePolarity();            
        } else if (thisIsNegative && !divisorIsNegative) {
            q.reversePolarity();
        }
      
        return printQR(q, r, divisor);
    }
    
    /**
     * convenience method to create an instance with the value of Long.MAX_VALUE,
     * ((1<<63)-1)
     * 
     * @return 
     */
    public static VeryLargeNumber createMaxLong() {
        
        VeryLargeNumber maxLong = new VeryLargeNumber(0);
        
        //max long = 9223372036854775807
        if (BASE == 10) {
            
            maxLong.setInternalArray(new int[]{
                9, 2, 2, 3, 3, 7, 2, 0, 3, 6, 8, 5, 4, 7, 7, 5, 8, 0, 7
            }, 19, true);
            
        } else if (BASE == ((1<<30)-1)) {
           
            maxLong.setInternalArray(new int[]{8, 16, 7}, 3, true);
            
        } else {
            
            throw new IllegalStateException("code needs to be adapted "
                + " for BASE=" + BASE);
        }
        
        return maxLong;
    }
    
    /**
     * prints the result of division's quotient, remainder and divisor as a 
     * double number string.
     * 
     * @param q
     * @param r
     * @param divisor
     * @return 
     */
    private String printQR(VeryLargeNumber q, VeryLargeNumber r, 
        VeryLargeNumber divisor) {
                        
        StringBuilder sb = new StringBuilder(q.toString());
       
        sb.append(".");
        
        VeryLargeNumber maxLong = VeryLargeNumber.createMaxLong();
        if (divisor.compareTo(maxLong) > 0) {
            // this is effectively zero
            return sb.append("0").toString();
            //throw new IllegalStateException("divisor is larger than 2^63 - 1");
        }
       
        long numerator = Long.valueOf(r.toString());
        
        long denominator = Long.valueOf(divisor.toString());
        
        double mantissa = (double)numerator/(double)denominator;
        
        // trim off 0. or -0.
        String mantissaStr = Double.toString(mantissa);
        int idx = mantissaStr.indexOf(".");
        mantissaStr = mantissaStr.substring(idx + 1);
        
        sb.append(mantissaStr);
        
        return sb.toString();
    }

    /**
     * compare the number within this instance to the number within other and
     * return -1 if this is smaller, 0 if this is equal to other, and +1 if this
     * is larger than other.
     * 
     * @param other
     * @return 
     */
    @Override
    public int compareTo(VeryLargeNumber other) {
        
        int nLenOther = other.getInternalArraySize();
       
        int nLenThis = nLen;
       
        if (nLen > nLenOther) {
                
            if (isPositive) {
               return 1;
            }
            
            // compare all under nLenO in both
            nLenThis = nLenOther;
            
        } else if (nLen < nLenOther) {
            
            if (other.isPositive()) {
               return -1;
            }
            
            // compare all under nLen in both
            nLenThis = nLen;
        }
        
        int[] b = other.a;
        
        for (int i = 0; i < nLenThis; i++) {
            int ta = isPositive ? a[i] : -1*a[i];
            int tb = other.isPositive() ? b[i] : -1*b[i];
            if (ta > tb) {
                return 1;
            } else if (ta < tb) {
                return -1;
            }
        }
        
        return 0;
    }
    
    /**
     * populate this instance with the effective value of number.
     * 
     * @param number 
     */
    private void createInternalNumber(int number) {
        
        if (number == 0) {
            
            nLen = 1;
            
            return;
        }
                        
        if (number < 0) {
            
            // hold sign and adapt code to handle...
            isPositive = false;
            
            if (number == Integer.MIN_VALUE) {
                // special handling because with sign change, it overflows 
                // an int
                //2147483648
                
                if (BASE == 10) {
                    
                    a = new int[] {2, 1, 4, 7, 4, 8, 3, 6, 4, 8};
                
                    nLen = 10;
                    
                } else if (BASE == ((1<<30)-1)) {
                    //BASE=1073741823
                    
                    a = new int[] {2, 2};
                    
                    nLen = 2;
                    
                } else {
                    throw new IllegalStateException("code needs to be adapted "
                        + " for BASE=" + BASE);
                }
              
                return;
            }
            
            number *= -1;
        }
                
        int i = 0;
        
        while (number > 0) {

            expandIfNeeded(i + 1);

            //TODO: caveats w/ java modulus operator.  might be able to impl this faster too.
            int tmp = number % BASE;

            //2147483647
            //1073741823
            a[i] = tmp;

            number /= BASE;

            i++;
            
            nLen = i;
        }
        
        // reverse the array so that number=1234 results in a=[1, 2, 3, 4]
        reverse();
    }
    
    /**
     * increment this number value
     */
    public void increment() {
        
        if (isPositive) {
            incrementPositive();
        } else {
            incrementNegative();
        }
    }
    
    /**
     * an increment specifically for use when this number is a negative number
     */
    public void incrementNegative() {
	
        int i = (nLen - 1);
        
        boolean checkReduce = false;
                         
        while (i > -1) {
            
            if (a[i] == 0) {
                
                a[i] = BASE - 1;
                
                i--;
                
                if (i == 0) {
                    checkReduce = true;
                }
               
            } else {
                
                a[i]--;
                
                if (checkReduce) {
                    moveUpIfStartsWithZeros();
                }
                
                break;
            }
        }
        
        if (isZero()) {
            isPositive = true;
        }
	}
    
    /**
     * an increment specifically for use when this number is a positive number.
     */
    public void incrementPositive() {
	
        int i = (nLen - 1);

        while (i > -1) {
            
            a[i]++;

            if (a[i] == BASE) {

                a[i] = 0;  

                if (i == 0) {
                    
                    expandIfNeeded(nLen + 1);
                    
                    //carry over to a new power
                    nLen++;
                    
                } else {
                    
                    i--;
                }
                
            } else {
                
                break;
            }
        }
	}
    
    /**
     * reverse the order of items in a
     */
    private void reverse() {
        
        if (nLen < 2) {
            return;
        }
                
        int end = nLen >> 1;
        
        for (int i = 0; i < end; i++) {
            
            int idx2 = nLen - i - 1;
            
            int swap = a[i];
            a[i] = a[idx2];
            a[idx2] = swap;
        }
    }
    
    /**
     * expand the backing array a if needed so that it can hold nTotal items
     * @param nTotal 
     */
    private void expandIfNeeded(int nTotal) {
        
        if (nTotal > a.length) {
            
            int n2 = a.length + 10;
            
            if (nTotal > n2) {
                n2 = nTotal;
            }
            
            a = Arrays.copyOf(a, n2);            
        }
    }
    
    /**
     * insert an empty item at the top of the array.  the method internally
     * moves down all items currently in a by 1 after expanding the a if needed.
     */
    private void insertSpaceAtTopOfArray() {
        
        if (a.length >= (nLen + 1)) {
            
            for (int i = (nLen - 1); i > -1; i--) {
                a[i + 1] = a[i];
            }
            a[0] = 0;
            
        } else {
            int[] xx = new int[nLen + 1];
            System.arraycopy(a, 0, xx, 1, nLen);
            a = xx;
        }
        
        nLen++;
    }
    
    /**
     * make a copy of this instance with a different identity but same
     * values in the member variables.  it's expected that compareTo
     * is used for comparison, but equals will also return true for comparison
     * of the clone with the original instance.
     * 
     * @return
     * @throws CloneNotSupportedException 
     */
    @Override
    public VeryLargeNumber clone() throws CloneNotSupportedException {
                
        VeryLargeNumber clone = (VeryLargeNumber) super.clone();
        
        int[] b = Arrays.copyOf(a, a.length);
        
        clone.setInternalArray(b, nLen, isPositive);
        
        return clone;
    }

    /**
     * compare the number value of this instance to another and return true
     * if they are the same.
     * 
     * @param other
     * @return 
     */
    @Override
    public boolean equals(Object other) {
        
        if (other == null) {
            return false;
        }
        
        if (!(other instanceof VeryLargeNumber)) {
            return false;
        }
        
        int comp = compareTo((VeryLargeNumber)other);
        
        return (comp == 0);
    }

    @Override
    public int hashCode() {
        return super.hashCode();
    }
    
    /**
     * reset the content of this instance to be the same as the copied content of
     * copyThis.
     * 
     * @param copyThis 
     */
    public void resetTo(VeryLargeNumber copyThis) {
        
        a = Arrays.copyOf(copyThis.a, copyThis.nLen);
        
        isPositive = copyThis.isPositive;        
    }
    
    /**
     * method purely for testing.  TODO: should be added to an aspect woven for tests
     * only.
     * 
     * @param b array of numbers composing the large number.  note that no
     * checks are done to assert that the numbers are positive as this is a
     * method meant to be used in testing only.
     * @param newNLen 
     * @param sign the number array is positive or negative
     */
    protected void setInternalArray(int[] b, int newNLen, boolean sign) {
        a = b;
        nLen = newNLen;
        isPositive = sign;
    }
    
    /**
     * return the size of the internal array a
     * @return 
     */
    protected int getInternalArraySize() {
        return nLen;
    }
    
    /**
     * get a copy of the internal array a
     * 
     * @return 
     */
    protected int[] getInternalArray() {
        return Arrays.copyOf(a, nLen);
    }
    
    /**
     * return the number value of this instance as a string.  note that this
     * will be a problem if the value is > 9223372036854775807
     * @return 
     */
    @Override
    public String toString() {
        
        if (isZero()) {
            return "0";
        }
                
        StringBuilder sb = new StringBuilder();
        if (!isPositive) {
            sb.append("-");
        }
        List<Long> prevSums = new ArrayList<Long>();
        long prevSum = 0;
        long sum = 0;
        int m;
        for (int i = 0; i < nLen; i++) {
            m = nLen - 1 - i;
            long factor = 1;
            for (int ii = 0; ii < m; ii++) {
                factor *= BASE;
            }
            int ai = a[i];
            long v = ai * factor;
            if (sum < 0) {
                prevSums.add(prevSum);
                sum *= -1;
            }
            prevSum = sum;
            sum += v;
        }
        if (!prevSums.isEmpty()) {
            for (Long ps : prevSums) {
                sb.append(ps).append(" + ");
            }
            sb.append(sum);
        } else {
            sb.append(sum);
        }
        return sb.toString();
    }
}
