package algorithms.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;
//import java.util.logging.Level;
//import java.util.logging.Logger;

/**
 * NOTE: still testing this class.  it needs more testing for operations that
 * change the instance value from positive to negative.
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
 
    //1073741823 
    public final static int BASE = (1<<30) - 1;
    
    private int[] a;
    
    private int nLen = 0;
    
    private boolean isPositive = true;
    
    /** a temporary work around holding non-integer values as a result
    of inverse operations.
    */ 
    private boolean hasADoubleValue = false;
    /**
    a temporary work around holding non-integer values as a result
    of inverse operations.  this should not be used unless 
    * hasADoubleValue == true.
    * Note, that no other logic is in place to use this except as a
    * result of inverse operation, that is pow(-number) and any subsequent
    * operation is lost.
    * needs to be refactored for any real use.  (in which case
    * the inverse should be replaced with Picarte's or Newton's).
     */
    private double doubleValue = 0;
    
    //private Logger log = Logger.getLogger(this.getClass().getName());
    
    public VeryLargeNumber(int number) {
        
        a = new int[10];
        
        createInternalNumber(number);
    }
    
    /**
     * a temporary work around for holding the results of inverse operations.
     * retrieve the result with getDoubleValueIfExists().
     * NOTE that this instance is not compatible with any other operations
     * at this time.  It's merely a way to return a value from the pow operation.
     * @param value 
     */
    VeryLargeNumber(double value) {
        
        a = new int[1];
        nLen = 1;
        
        hasADoubleValue = true;
        doubleValue = value;
    }
    
    public double getDoubleValueIfExists() {
        if (hasADoubleValue) {
            return doubleValue;
        }
        return 0;
    }
    
    /**
     * add addThis to this instance.
     * 
     * @param addThis number to add to this
     */
    public void add(VeryLargeNumber addThis) {
                
        VeryLargeNumber result = add(this, addThis);
        
        resetTo(result);
    }
    
    static VeryLargeNumber add(final VeryLargeNumber number1, 
        final VeryLargeNumber addThis) {
                
        //log.log(Level.FINEST, "add: {0} + {1}", new String[]{toString(), 
        //    addThis.toString()});
        
        boolean thisIsLarger = (addThis.nLen <= number1.nLen);
        
        int[] longer;
        int[] shorter;
        int n, idxOffsetShorter;
        boolean longerIsPositive, shorterIsPositive;
        
        if (thisIsLarger) {
            longer = number1.getInternalArray();
            shorter = addThis.getInternalArray();
            n = number1.nLen;
            idxOffsetShorter = n - addThis.nLen;
            longerIsPositive = number1.isPositive;
            shorterIsPositive = addThis.isPositive();
        } else {
            longer = addThis.getInternalArray();
            shorter = number1.getInternalArray();
            n = addThis.nLen;
            idxOffsetShorter = addThis.nLen - number1.nLen;
            longerIsPositive = addThis.isPositive();
            shorterIsPositive = number1.isPositive;
        }
        
        int[] output = new int[n + 1];
        int[] outputLength = new int[1];
        
        boolean resultIsPositive = add(longer, shorter, n, idxOffsetShorter, 
            longerIsPositive, shorterIsPositive, output, outputLength);
        if (output[0] == 0) {
            int nLength = VeryLargeNumber.moveUpIfStartsWithZeros(output, 
                outputLength[0]);
            outputLength[0] = nLength;
        }
        
        VeryLargeNumber result = new VeryLargeNumber(0);
        result.setInternalArray(output, outputLength[0], resultIsPositive);
        
        return result;
    }
    
    /**
     * add shorter to longer.  NOTE that the use of this method requires the
     * invoker to have already made room for a carry, that is output.length
     * should be one larger than longer.length.
     * ALSO it is the responsibility of the invoker to move up the numbers
     * in output if they start w/ zero.
     * 
     * @return returns true if the result in output is a positive number, else
     * returns false.
     */
    private static boolean add(int[] longer, int[] shorter, int n, int idxOffsetShorter, 
        boolean longerIsPositive, boolean shorterIsPositive, int[] output,
        int[] outputLength) {
    
        int	carry = 0;
        int sum = 0;
        int idxShorter;
        boolean resetCarry = true;
                
        for (int i = (n - 1); i > -1; i--) {
                   
            resetCarry = true;
            
            idxShorter = i - idxOffsetShorter;
           
            int tl = longerIsPositive ? longer[i] : -1*longer[i];
            
            if (idxShorter < 0) {
                
                if ((carry == -1) && (longer[i] == 0)) {
                                        
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

            output[i] = longerIsPositive ? sum : -1*sum;
        }
        
        boolean resultIsPositive = longerIsPositive;
        
        outputLength[0] = n;
        
        if (carry != 0) {
            
            if (output.length <= longer.length) {
                throw new IllegalArgumentException(
                "output must be given to this method with a length that is at " 
                + "least one larger than longer.length");
            }
            
            VeryLargeNumber.moveDown(output, 1);
            
            if (carry < 0) {
                resultIsPositive = false;
                carry *= -1;
            }
            
            output[0] = carry;
            
            outputLength[0]++;
            
        }
        
        return resultIsPositive;
    }
    
    int countNumberOfDigits(int number) {
        
        int count = 0;
        
        if (number == 0) {
            return 1;
        }
        
        if (number < 0) {
            number *= -1;
        }
        
        while (number > 0) {
            number /= BASE;
        }
        
        return count;
    }
    
    /**
     * if there are leading 0's in the array a, move items below
     * them up to fill them and reduce the size of nLen. if there are
     * no numbers left, the result is nLen=1 to result in a value of 0
     * for this instance's number.
     */
    private void moveUpIfStartsWithZeros() {
        nLen = moveUpIfStartsWithZeros(this.a, this.nLen);
    }
    
    /**
     * if the given array in starts with 0's, moves the elements below it
     * (at higher indexes) up to remove the zero's and returns the new total
     * length of the usable portion of array in (that is highIndex + 1).
     * @param in
     * @param length
     * @return 
     */
    static int moveUpIfStartsWithZeros(int[] in, int length) {
        
        // move up if needed
        int firstNonZeroIdx = -1;
        for (int i = 0; i < length; i++) {
            if (in[i] == 0) {
                firstNonZeroIdx = i + 1;
            } else {
                break;
            }
        }
        
        if (firstNonZeroIdx > -1) {
            
            if (length == firstNonZeroIdx) {
                // it's all zeros
                return 1;
            }

            moveUp(in, length, firstNonZeroIdx);
            
            length -= firstNonZeroIdx;
        }
        
        return length;
    }
    
    static void moveUp(int[] in, int length, int shift) {
        
        for (int i = 0; i < (length - shift); i++) {
            in[i] = in[i + shift];
        }

        Arrays.fill(in, (length - shift), length, 0);
    }
    
    static void moveDown(int[] in, int shift) {
        
        if (shift > (in.length - 1)) {
            throw new IllegalArgumentException("shift is larger than array");
        }
        
        if (shift == 0) {
            return;
        }
        
        for (int i = (in.length - 1); i >= shift; i--) {
            in[i] = in[i - shift];
        }
        
        Arrays.fill(in, 0, shift, 0);
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

    /*
    https://en.wikipedia.org/wiki/Multiplication_algorithm
    summarizes use of:
    https://en.wikipedia.org/wiki/Karatsuba_algorithm
    https://en.wikipedia.org/wiki/Toom%E2%80%93Cook_multiplication
    https://en.wikipedia.org/wiki/Sch%C3%B6nhage%E2%80%93Strassen_algorithm
    Schonhage-Strassen is fastest for numbers larger than 2^32768.       
    */
    public void multiply(VeryLargeNumber number) {
        
        boolean numbersAreSmall = (this.nLen < 2) || (number.nLen < 2);
        
        if (numbersAreSmall) {
            VeryLargeNumber result = multiplySmall(number);
            this.resetTo(result);
        } else {
            VeryLargeNumber result = VeryLargeNumber.karatsuba(
                this, number);
            this.resetTo(result);
        }
        
    }
    
    /**
     * the O(N^2) method for multiplication:
     * http://www.cs.utexas.edu/users/djimenez/utsa/cs3343/lecture20.html
     * 
     * @param num2
     */
    protected VeryLargeNumber multiplySmall(VeryLargeNumber num2) {
        
        // zero check
        VeryLargeNumber zero = new VeryLargeNumber(0);
        if ((num2.compareTo(zero) == 0) || (compareTo(zero) == 0)) {
            return zero;
        }
        
        // one check:
        VeryLargeNumber one = new VeryLargeNumber(1);
        if (num2.compareTo(one) == 0) {
            VeryLargeNumber result = new VeryLargeNumber(0);
            result.resetTo(this);
            return result;
        } else if (compareTo(one) == 0) {
            VeryLargeNumber result = new VeryLargeNumber(0);
            result.resetTo(num2);
            return result;
        }
        
        int length0 = (nLen >= num2.nLen) ? nLen : num2.nLen;
        
        int shiftSum = 0;
        for (int i = 0; i <= length0; i++) {
            shiftSum += i;
        }
        // making the array longer to handle multiplication, and add a
        // right shift by one to allow space for a carry
        int length = length0 + shiftSum + 1;
        
        int[] bb = num2.getInternalArray();
        int[] aa;
        int[] cc;
        int[] pp;
        
        /*
        1294 is: a[0] is 1, a[1] is 2, a[2] is 9, and a[3] is 4
        
        so if multiplying by a larger number, need to shift the numbers down.
        */

        int[] tmp = new int[length];
        System.arraycopy(bb, 0, tmp, 1, num2.nLen);
        bb = tmp;
        
        tmp = new int[length];
        System.arraycopy(a, 0, tmp, 1, nLen);
        aa = tmp;
        cc = new int[length + 1];
        pp = new int[length];
     
        int[] outputLength = new int[1];
        
	    //cc will accumulate the sum of partial products.  It's initially 0.
        for (int i = 1; i < (length0 + 1); i++) {
            
            //multiply bb by digit aa[i]
            multiplyOneDigit(bb, pp, aa[i]);
            
            //shift the partial product i spaces up to increase by BASE steps
            moveDown(pp, i - 1);
            
            //add result to the running sum
            add(cc, pp, length, 0, true, true, cc, outputLength);        
        }
        
        // remove leading zeros
        length = moveUpIfStartsWithZeros(cc, length);
        
        boolean resultIsPositive = this.isPositive;
        if (!num2.isPositive) {
            resultIsPositive = !resultIsPositive;
        }
        
        // find the last 0
        for (int i = (length - 1); i > 0; i--) {
            if (cc[i] == 0) {
                length--;
            } else {
                break;
            }
        }
        
        VeryLargeNumber result = new VeryLargeNumber(0);
        result.setInternalArray(cc, length, resultIsPositive);
        
        return result;
    }
    
    /**
     * http://www.cs.utexas.edu/users/djimenez/utsa/cs3343/lecture20.html
     * 
     * NOTE: in and out must have same lengths (and have been created
     * with the same BASE).  To account for overflow, out and in should be
     * one larger than the original length of in.
     * 
     * output = input * d
     * 
     * @param in
     * @param out
     * @param d
     */
    protected void multiplyOneDigit(int in[], int out[], int d) {
        
        if (in == null) {
            throw new IllegalArgumentException("in cannot be null");
        }
        if (out == null) {
            throw new IllegalArgumentException("out cannot be null");
        }
        if (in.length != out.length) {
            throw new IllegalArgumentException(
                "in and out must have the same lengths");
        }
        
        int i, carry;

        // no extra overflow to add yet
        carry = 0;

        // for each digit, starting with least significant...
        for (i = 0; i < in.length; i++) {

            // multiply by digit d
            out[i] = d * in[i];

            // add in any overflow from the last digit
            out[i] += carry;

            // if this product is too big to fit in a digit...
            if (out[i] >= BASE) {

                // handle the overflow
                carry = out[i] / BASE;
                
                out[i] %= BASE;
                
            } else {
                
                // no overflow
                carry = 0;
            }
        }
        
        if (carry > 0) {
            throw new IllegalArgumentException("overflow in multiplication!  " 
            + "Increase the lengths of in and out by one.\n");
        }
    }
    
    /**
       https://en.wikipedia.org/wiki/Karatsuba_algorithm

       runtime O(n_digits^(lg2(3)))

       Let x and y be represented as n-digit strings in some base B.
       x = x_1*B^m + x_0
       y = y_1*B^m + y_0,  where x_0 and y_0 are less than B^m
       
       x*y = (x_1*B^m + x_0)(y_1*B^m + y_0)
       x*y = z_2*B^{2m} + z_1*B^m + z_0
       
       z_2 = x_1*y_1
       z_1 = x_1*y_0 + x_0*y_1
       z_0 = x_0*y_0

       z_1 = x_1*y_0 + x_0*y_1
       z_1 = (x_1 + x_0)(y_1 + y_0) - x_1*y_1 - x_0*y_0
       z_1 = (x_1 + x_0)(y_1 + y_0) - z_2 - z_0

       z_2 = x_1*y_1
       z_0 = x_0*y_0
       
       x*y = (b^2 + b)*x_1*y_1 - b*(x_1 - x_0)(y_1 - y_0) + (b + 1)*x_0*y_0
       where b is the power where the split occurs of x_1.
        
    Below, is an iterative version of karatsuba.  Java doesn't use tail recursion
    at this time, so cannot use a recursive version of the algorithm for very
    large numbers.
    */
    static VeryLargeNumber karatsuba(VeryLargeNumber num1, VeryLargeNumber num2) {
                
        try {
            num1 = num1.clone();
            num2 = num2.clone();
        } catch (CloneNotSupportedException e) {
            System.err.println(e.getMessage());
        }
        
        boolean resultIsPositive = true;
        if (!num1.isPositive) {
            resultIsPositive = !resultIsPositive;
            num1.reversePolarity();
        }
        if (!num2.isPositive) {
            resultIsPositive = !resultIsPositive;
            num2.reversePolarity();
        }
                
        java.util.Stack<VeryLargeNumber> stack1 = new java.util.Stack<VeryLargeNumber>();
        java.util.Stack<VeryLargeNumber> stack2 = new java.util.Stack<VeryLargeNumber>();
        java.util.Stack<String> stackKey = new java.util.Stack<String>();
        
        // key=current nIter, value = nIter+var to place results in
        Map<Integer, String> prevMap = new HashMap<Integer, String>();
                
        Map<String, Integer> keyM2Map = new HashMap<String, Integer>();
        
        stack1.push(num1);
        stack2.push(num2);
        stackKey.push(Integer.toString(0));
        
        Map<String, VeryLargeNumber> resultMap = new HashMap<String, VeryLargeNumber>();
              
        int nIter = -1;
        
        while(!stack1.empty()) {
            
            nIter++;
            
            num1 = stack1.pop();
            num2 = stack2.pop();
            String currentKey = stackKey.pop();
            
            // zero check
            VeryLargeNumber zero = new VeryLargeNumber(0);
            if ((num2.compareTo(zero) == 0) || (num2.compareTo(zero) == 0)) {
                return zero;
            }
 
            // one check:
            VeryLargeNumber one = new VeryLargeNumber(1);
            if (num2.compareTo(one) == 0) {
                
                VeryLargeNumber result = new VeryLargeNumber(0);
                result.resetTo(num1);
                
                resultMap.put(currentKey, result);
                
                continue;
                
            } else if (num1.compareTo(one) == 0) {
                
                VeryLargeNumber result = new VeryLargeNumber(0);
                result.resetTo(num2);
                
                resultMap.put(currentKey, result);
                              
                continue;
            }

            // calculates the size of the numbers
            int sz1 = num1.nLen;
            int sz2 = num2.nLen;

            if ((sz1 < 2) || (sz2 < 2)) {

                VeryLargeNumber result = num1.multiplySmall(num2);

                resultMap.put(currentKey, result);
                                
                continue;
            }

            int m = (sz2 > sz1) ? sz2 : sz1;
            int m2 = (m >> 1);
            if (m2 > (num1.nLen - 1)) {
                m2 = num1.nLen - 1;
            } else if (m2 > (num2.nLen - 1)) {
                m2 = num2.nLen - 1;
            }
            // split the digit sequences about the middle
            VeryLargeNumber[] highLow1 = num1.splitAt(m2);
            VeryLargeNumber[] highLow2 = num2.splitAt(m2);
            
            VeryLargeNumber high1 = highLow1[0];
            VeryLargeNumber low1 = highLow1[1];
            
            assert(high1.nLen + low1.nLen == num1.nLen);

            VeryLargeNumber high2 = highLow2[0];
            VeryLargeNumber low2 = highLow2[1];

            assert(high2.nLen + low2.nLen == num2.nLen);
            
            VeryLargeNumber z1pt1 = VeryLargeNumber.add(low1, high1);
            VeryLargeNumber z1pt2 = VeryLargeNumber.add(low2, high2);
            
            /*
            VeryLargeNumber z0 = VeryLargeNumber.karatsuba(low1, low2);
            VeryLargeNumber z1 = VeryLargeNumber.karatsuba(z1pt1, z1pt2);            
            VeryLargeNumber z2 = VeryLargeNumber.karatsuba(high1, high2);
            z1.subtract(z2);
            z1.subtract(z0);
            
            // result = (z2 * BASE^(2*m2)) + (z1 * BASE^(m2)) + (z0)

            int z2Length = 2*m2 + z2.nLen;
            z2.setInternalArray(Arrays.copyOf(z2.a, z2Length), z2Length, z2.isPositive);

            int z1Length = m2 + z1.nLen;
            z1.setInternalArray(Arrays.copyOf(z1.a, z1Length), z1Length, z1.isPositive);

            VeryLargeNumber result = VeryLargeNumber.add(z2, z1);
            result.add(z0);
            */
            
            String nIterStr = Integer.toString(nIter);
            
            String key = nIterStr + "z0";
            stack1.push(low1);  
            stack2.push(low2); 
            stackKey.push(key);
            keyM2Map.put(key, Integer.valueOf(0));
            
            key = nIterStr + "z1";
            stack1.push(z1pt1);  
            stack2.push(z1pt2); 
            stackKey.push(key);
            keyM2Map.put(key, Integer.valueOf(0));
            
            key = nIterStr + "z2";
            stack1.push(high1);
            stack2.push(high2);
            stackKey.push(key);
            keyM2Map.put(key, Integer.valueOf(m2));
            
            prevMap.put(Integer.valueOf(nIter), currentKey);
        }
        
        SortedSet<Integer> keys = new TreeSet<Integer>(prevMap.keySet());
        
        Integer keyNIter = keys.isEmpty() ? Integer.valueOf(0) : keys.last();
                
        while (keyNIter != null) {
                        
            if (keyNIter.equals(Integer.valueOf(0))) {
                
                if (resultMap.containsKey(keyNIter.toString())) {
                    VeryLargeNumber result = resultMap.get(keyNIter.toString());
                    if (!resultIsPositive) {
                        result.reversePolarity();
                    }
                    
                    assert(keys.isEmpty());
                    
                    return result;
                }
            }
                                                
            if (resultMap.containsKey(keyNIter.toString() + "z0") &&
                resultMap.containsKey(keyNIter.toString() + "z1") &&
                resultMap.containsKey(keyNIter.toString() + "z2")) {

                String prevKey = prevMap.get(keyNIter);

                if (prevKey != null) {

                    keys.remove(keyNIter);
                    
                    // process the whole set of z0, z1, and z2
                    VeryLargeNumber z0 = resultMap.get(keyNIter.toString() + "z0");
                    VeryLargeNumber z1 = resultMap.get(keyNIter.toString() + "z1");
                    VeryLargeNumber z2 = resultMap.get(keyNIter.toString() + "z2");
                    
                    //int shiftz0 = tmpKeyM2Map.get(keyNIter.toString() + "z0");
                    //int shiftz1 = tmpKeyM2Map.get(keyNIter.toString() + "z1");
                    int shiftz2 = keyM2Map.get(keyNIter.toString() + "z2");
                    
                    z1.subtract(z0);
                    z1.subtract(z2);
                    
                    // m2 will always be <= ((1<<30)-1)=1073741823
                    // result = (z2 * BASE^(2*m2)) + (z1 * BASE^(m2)) + (z0)

                    int z2Length = (2 * shiftz2) + z2.nLen;
                    z2.setInternalArray(Arrays.copyOf(z2.a, z2Length), z2Length, 
                        z2.isPositive);

                    int z1Length = shiftz2 + z1.nLen;
                    z1.setInternalArray(Arrays.copyOf(z1.a, z1Length), z1Length, 
                        z1.isPositive);
                    
                    VeryLargeNumber result2 = VeryLargeNumber.add(z2, z1);
                    result2.add(z0);
                    
                    resultMap.put(prevKey, result2);
            
                    keyNIter = keys.isEmpty() ? Integer.valueOf(0) : keys.last();
  
                }
            }
        }
            
        return null;
    }

    /**
     * a split of the array at BASE^m is performed and the returned array is
     * new VeryLargeNumber[]{highDigits, lowDigits}.
     * @param index
     * @return 
     */
    VeryLargeNumber[] splitAt(int m) {
               
        if (m > (nLen - 1)) {
            throw new IllegalArgumentException("m is larger than array size");
        }
        if (m < 1) {
            throw new IllegalArgumentException("m must be larger than 0");
        }
       
        int hi0 = 0;
        int hi1Excl = nLen - m;
        int hiLen = hi1Excl - hi0;
        
        int lo0 = nLen - m;
        int lo1Excl = nLen;
        int loLen = m;
      
        VeryLargeNumber[] result = new VeryLargeNumber[2];
        result[0] = new VeryLargeNumber(0);
        result[0].setInternalArray(Arrays.copyOfRange(a, hi0, hi1Excl), hiLen, isPositive);
        
        result[1] = new VeryLargeNumber(0);
        result[1].setInternalArray(Arrays.copyOfRange(a, lo0, lo1Excl), loLen, isPositive);
        
        return result;
    }
    
    private static void printDebug(String label, VeryLargeNumber number) {
        int[] tmp = Arrays.copyOf(number.a, number.nLen);
        System.out.println(label + " " + Arrays.toString(tmp));
    }
            
    /**
    For exponentiation:
        https://en.wikipedia.org/wiki/Exponentiation_by_squaring
    */
    public VeryLargeNumber pow(int x) {
        
        return VeryLargeNumber.expBySquaring(this, x);
    }
    
    /**     
     http://en.wikipedia.org/wiki/Exponentiation_by_squaring
     
     could replace with a 2^k method with precomputed values.
     * @param number
     * @param x
     * @return 
     */
    private static VeryLargeNumber expBySquaring(VeryLargeNumber number, int x) {
        
        //changed the recursion to use iteration...
        
        VeryLargeNumber correctionForOddX = null;
        
        int nIter = 0;
        
        while (true) {
            
            if (x < 0) {
                double inverted = number.inverse();
                double result0 = Math.pow(inverted, -1*x);
                VeryLargeNumber result = new VeryLargeNumber(result0);
                return result;
            } else if (x == 0) {
                return new VeryLargeNumber(1);
            } else if (x == 1) {
                if (correctionForOddX != null) {
                    number.multiply(correctionForOddX);
                }
                return number;
            } else {
                VeryLargeNumber num1 = new VeryLargeNumber(0);
                num1.resetTo(number);
                num1.multiply(number);
                if ((x & 1) == 1) {
                    VeryLargeNumber f = new VeryLargeNumber(0);
                    f.resetTo(number);
                    correctionForOddX = f;
                    x--;
                }
                                
                number = num1;
                
                x >>= 1;
                
                nIter++;
            }
        }
    }
    
    /*
    http://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Modular_integers
    
    function inverse(a, p)
    t := 0;     newt := 1;    
    r := p;     newr := a;    
    while newr ≠ 0
        quotient := r div newr
        (r, newr) := (newr, r - quotient * newr)
        (t, newt) := (newt, t - quotient * newt) 
    if degree(r) > 0 then 
        return "Either p is not irreducible or a is a multiple of p"
    return (1/r) * t
    */
        
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
        return (a[nLen - 1] == 0);
    }
    
    public boolean isOdd() {
        return ((a[nLen - 1] & 1) == 1);
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
     * 
     * 
     * @param number
     * @return 
     */
    double inverse() {
     
        return inverseByEuclidean(this);        
    }
    
    /**
     * using Euclidean division, divide this by divisor and return the result.
     * 
     * @param number
     * @return 
     */
    double inverseByEuclidean(VeryLargeNumber number) {
        
        if (number.isZero()) {
            throw new IllegalArgumentException("Cannot divide by zero");
        }
        
        boolean divisorIsNegative = !number.isPositive();
                
        if (divisorIsNegative) {
            number.reversePolarity();
        }
        
        // both this number and divisor are positive or zero
        VeryLargeNumber q = new VeryLargeNumber(0);
        
        VeryLargeNumber r = new VeryLargeNumber(1);
        
        // while  R ≥ D
        while (r.compareTo(number) > -1) {
            
            q.increment();
            
            r.subtract(number);
        }
        
        if (!r.isPositive()) {
            // the while loop proceeds one step too far so reverse by 1 loop
            q.reversePolarity();
            q.increment();
            q.reversePolarity();
            r.add(number);
        }
        
        if (divisorIsNegative) {
            number.reversePolarity();            
        }
        
        if (divisorIsNegative) {
            q.reversePolarity();            
        }
      
        long numerator = Long.valueOf(r.toString());
        
        long denominator = Long.valueOf(number.toString());
        
        double mantissa = (double)numerator/(double)denominator;
                
        return mantissa;
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
        
        if ((numerator == 1) && (mantissa < 0)) {
            if (sb.charAt(0) != '-') {
                sb.insert(0, "-");
            }
        }
        
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

            int tmp = number % BASE;

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
        
        VeryLargeNumber clone = new VeryLargeNumber(0);
        
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
        
        nLen = copyThis.nLen;
        
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
