package algorithms;

import algorithms.util.IntIntDouble;
import algorithms.util.PairInt;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import java.util.List;

/**
 *
 * @author nichole
 */
public class QuickSort {
    
     /**
     * sort a from index idxLo to idxHi, inclusive.
     * @param a
     */
    public static void sort(float[] a) {
        sort(a, 0, a.length - 1);
    }
    
    public static void descendingSort(int[] a, int[] b) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        
        descendingSort(a, b, 0, a.length - 1);
    }
    
    public static void descendingSort(TIntList a, List<? extends Object> b) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.size() != b.size()) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        
        descendingSort(a, b, 0, a.size() - 1);
    }
    
    public static void descendingSort(List<? extends Number> a, List<? extends Object> b) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.size() != b.size()) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        
        descendingSort(a, b, 0, a.size() - 1);
    }
    
    public static void descendingSort(double[] a) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        
        descendingSort(a, 0, a.length - 1);
    }
    
    public static void sortBy1stArg(int[] a, Object[][] b) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        
        sortBy1stArg(a, b, 0, a.length - 1);
    }
    
    public static void sortBy1stArg(int[] a, int[] b) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        
        sortBy1stArg(a, b, 0, a.length - 1);
    }
    
    public static void sortByA(IntIntDouble[] abc) {
        
        if (abc == null) {
            throw new IllegalArgumentException("abc cannot be null");
        }
        
        sortByA(abc, 0, abc.length - 1);
    }
    
    /**
     * sort a by ascending values and perform the same swap operation on b.
     * @param a
     * @param b 
     */
    public static void sortBy1stArg(float[] a, int[] b) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        
        sortBy1stArg(a, b, 0, a.length - 1);
    }
    
    public static void sortBy1stArg(int[] a, Object[] b) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        
        sortBy1stArg(a, b, 0, a.length - 1);
    }
    
    public static void sortBy1stArg(float[] a, Object[] b) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        
        sortBy1stArg(a, b, 0, a.length - 1);
    }
    
    public static void sortBy1stArg(TIntList a, TDoubleList b,
        TIntList c) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (c == null) {
            throw new IllegalArgumentException("c cannot be null");
        }
        if (a.size() != b.size()) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        if (a.size() != c.size()) {
            throw new IllegalArgumentException("a and v must be the same length");
        }
        
        sortBy1stArg(a, b, c, 0, a.size() - 1);
    }
    
    public static void sortBy1stArg(TIntList a, TIntList b) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.size() != b.size()) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        
        sortBy1stArg(a, b, 0, a.size() - 1);
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive.  Uses the optimized
     * qsort3 from the book "Programming in Pearls" by Jon Bentley.
     * @param a
     * @param b
     * @param idxLo
     * @param idxHi 
     */
    public static void sortBy1stArg(TIntList a, TDoubleList b, 
        TIntList c, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (c == null) {
            throw new IllegalArgumentException("c cannot be null");
        }
        if (a.size() != b.size()) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        if (a.size() != c.size()) {
            throw new IllegalArgumentException("a and c must be the same length");
        }
        
        if (a.size() < 2) {
            return;
        }
        
        if (idxLo < idxHi) {

            int x = a.get(idxLo);
            int store = idxLo;
            int idxMid = idxHi + 1;

            while (true) {
                do {
                    store++;     
                } while ((store <= idxHi) && (a.get(store) < x));
                do {
                    idxMid--;
                } while (a.get(idxMid) > x);
                if (store > idxMid) {
                    break;
                }
                int swap = a.get(store);
                a.set(store, a.get(idxMid));
                a.set(idxMid, swap);
                swap = c.get(store);
                c.set(store, c.get(idxMid));
                c.set(idxMid, swap);
                
                double bSwap = b.get(store);
                b.set(store, b.get(idxMid));
                b.set(idxMid, bSwap);
            }
            int swap = a.get(idxLo);
            a.set(idxLo, a.get(idxMid));
            a.set(idxMid, swap);
            swap = c.get(idxLo);
            c.set(idxLo, c.get(idxMid));
            c.set(idxMid, swap);
            
            double bSwap = b.get(idxLo);
            b.set(idxLo, b.get(idxMid));
            b.set(idxMid, bSwap);
         
            sortBy1stArg(a, b, c, idxLo, idxMid - 1);

            sortBy1stArg(a, b, c, idxMid + 1, idxHi);
        }
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive.  Uses the optimized
     * qsort3 from the book "Programming in Pearls" by Jon Bentley.
     * @param a
     * @param b
     * @param idxLo
     * @param idxHi 
     */
    public static void sortBy1stArg(TIntList a, TIntList b, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("v cannot be null");
        }
        if (a.size() != b.size()) {
            throw new IllegalArgumentException("a and v must be the same length");
        }
        
        if (a.size() < 2) {
            return;
        }
        
        if (idxLo < idxHi) {

            int x = a.get(idxLo);
            int store = idxLo;
            int idxMid = idxHi + 1;

            while (true) {
                do {
                    store++;     
                } while ((store <= idxHi) && (a.get(store) < x));
                do {
                    idxMid--;
                } while (a.get(idxMid) > x);
                if (store > idxMid) {
                    break;
                }
                int swap = a.get(store);
                a.set(store, a.get(idxMid));
                a.set(idxMid, swap);
                swap = b.get(store);
                b.set(store, b.get(idxMid));
                b.set(idxMid, swap);
            }
            int swap = a.get(idxLo);
            a.set(idxLo, a.get(idxMid));
            a.set(idxMid, swap);
            swap = b.get(idxLo);
            b.set(idxLo, b.get(idxMid));
            b.set(idxMid, swap);
         
            sortBy1stArg(a, b, idxLo, idxMid - 1);

            sortBy1stArg(a, b, idxMid + 1, idxHi);
        }
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive.  Uses the optimized
     * qsort3 from the book "Programming in Pearls" by Jon Bentley.
     * @param a
     * @param b
     * @param idxLo
     * @param idxHi 
     */
    public static void sortByA(IntIntDouble[] abc, int idxLo, 
        int idxHi) {
        
        if (abc == null) {
            throw new IllegalArgumentException("abc cannot be null");
        }
        
        if (abc.length < 2) {
            return;
        }
        
        if (idxLo < idxHi) {

            int x = abc[idxLo].getA();
            int store = idxLo;
            int idxMid = idxHi + 1;

            while (true) {
                do {
                    store++;     
                } while ((store <= idxHi) && (abc[store].getA() < x));
                do {
                    idxMid--;
                } while (abc[idxMid].getA() > x);
                if (store > idxMid) {
                    break;
                }
                IntIntDouble swap = abc[store];
                abc[store] = abc[idxMid];
                abc[idxMid] = swap;
            }
            IntIntDouble swap = abc[idxLo];
            abc[idxLo] = abc[idxMid];
            abc[idxMid] = swap;
         
            sortByA(abc, idxLo, idxMid - 1);

            sortByA(abc, idxMid + 1, idxHi);
        }
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive.  Uses the optimized
     * qsort3 from the book "Programming in Pearls" by Jon Bentley.
     * @param a
     * @param idxLo
     * @param idxHi 
     */
    public static void sort(float[] a, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
        if (idxLo < idxHi) {

            float x = a[idxLo];
            int store = idxLo;
            int idxMid = idxHi + 1;

            while (true) {
                do {
                    store++;     
                } while ((store <= idxHi) && (a[store] < x));
                do {
                    idxMid--;
                } while (a[idxMid] > x);
                if (store > idxMid) {
                    break;
                }
                float swap = a[store];
                a[store] = a[idxMid];
                a[idxMid] = swap;
            }
            float swap = a[idxLo];
            a[idxLo] = a[idxMid];
            a[idxMid] = swap;
         
            sort(a, idxLo, idxMid - 1);

            sort(a, idxMid + 1, idxHi);
        }
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive and by descending values.  
     * The swap operations performed on a are performed on b.  Uses the optimized
     * qsort3 from the book "Programming in Pearls" by Jon Bentley.
     * @param a
     * @param idxLo
     * @param idxHi 
     */
    public static void descendingSort(int[] a, int[] b, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be same length");
        }
        if (a.length < 2) {
            return;
        }
        
        if (idxLo < idxHi) {

            int x = a[idxLo];
            int store = idxLo;
            int idxMid = idxHi + 1;

            while (true) {
                do {
                    store++;     
                } while ((store <= idxHi) && (a[store] > x));
                do {
                    idxMid--;
                } while (a[idxMid] < x);
                
                if (store > idxMid) {
                    break;
                }
                int swap = a[store];
                a[store] = a[idxMid];
                a[idxMid] = swap;
                swap = b[store];
                b[store] = b[idxMid];
                b[idxMid] = swap;
            }
            int swap = a[idxLo];
            a[idxLo] = a[idxMid];
            a[idxMid] = swap;
            swap = b[idxLo];
            b[idxLo] = b[idxMid];
            b[idxMid] = swap;
         
            descendingSort(a, b, idxLo, idxMid - 1);

            descendingSort(a, b, idxMid + 1, idxHi);
        }
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive and by descending values.  
     * The swap operations performed on a are performed on b.  Uses the optimized
     * qsort3 from the book "Programming in Pearls" by Jon Bentley.
     * @param <T>
     * @param a
     * @param b
     * @param idxLo
     * @param idxHi 
     */
    public static <T extends Object> void descendingSort(
        TIntList a, List<T> b, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.size() != b.size()) {
            throw new IllegalArgumentException("a and b must be same length");
        }
        if (a.size() < 2) {
            return;
        }
        
        if (idxLo < idxHi) {

            int x = a.get(idxLo);
            int store = idxLo;
            int idxMid = idxHi + 1;

            while (true) {
                do {
                    store++;     
                } while ((store <= idxHi) && (a.get(store) > x));
                do {
                    idxMid--;
                } while (a.get(idxMid) < x);
                
                if (store > idxMid) {
                    break;
                }
                int swap = a.get(store);
                a.set(store, a.get(idxMid));
                a.set(idxMid, swap);
                T swap2 = b.get(store);
                b.set(store, b.get(idxMid));
                b.set(idxMid, swap2);
            }
            int swap = a.get(idxLo);
            a.set(idxLo, a.get(idxMid));
            a.set(idxMid, swap);
            T swap2 = b.get(idxLo);
            b.set(idxLo, b.get(idxMid));
            b.set(idxMid, swap2);
         
            descendingSort(a, b, idxLo, idxMid - 1);

            descendingSort(a, b, idxMid + 1, idxHi);
        }
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive and by descending values.  
     * The swap operations performed on a are performed on b.  Uses the optimized
     * qsort3 from the book "Programming in Pearls" by Jon Bentley.
     * @param a
     * @param idxLo
     * @param idxHi 
     */
    public static <S extends Number, T extends Object> void descendingSort(
        List<S> a, List<T> b, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.size() != b.size()) {
            throw new IllegalArgumentException("a and b must be same length");
        }
        if (a.size() < 2) {
            return;
        }
        
        if (idxLo < idxHi) {

            float x = a.get(idxLo).floatValue();
            int store = idxLo;
            int idxMid = idxHi + 1;

            while (true) {
                do {
                    store++;     
                } while ((store <= idxHi) && (a.get(store).floatValue() > x));
                do {
                    idxMid--;
                } while (a.get(idxMid).floatValue() < x);
                
                if (store > idxMid) {
                    break;
                }
                S swap = a.get(store);
                a.set(store, a.get(idxMid));
                a.set(idxMid, swap);
                T swap2 = b.get(store);
                b.set(store, b.get(idxMid));
                b.set(idxMid, swap2);
            }
            S swap = a.get(idxLo);
            a.set(idxLo, a.get(idxMid));
            a.set(idxMid, swap);
            T swap2 = b.get(idxLo);
            b.set(idxLo, b.get(idxMid));
            b.set(idxMid, swap2);
         
            descendingSort(a, b, idxLo, idxMid - 1);

            descendingSort(a, b, idxMid + 1, idxHi);
        }
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive and by descending values.  
     * The swap operations performed on a are performed on b.  Uses the optimized
     * qsort3 from the book "Programming in Pearls" by Jon Bentley.
     * @param a
     * @param idxLo
     * @param idxHi 
     */
    public static void descendingSort(double[] a, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
        
        if (idxLo < idxHi) {

            double x = a[idxLo];
            int store = idxLo;
            int idxMid = idxHi + 1;

            while (true) {
                do {
                    store++;     
                } while ((store <= idxHi) && (a[store] > x));
                do {
                    idxMid--;
                } while (a[idxMid] < x);
                
                if (store > idxMid) {
                    break;
                }
                double swap = a[store];
                a[store] = a[idxMid];
                a[idxMid] = swap;
            }
            double swap = a[idxLo];
            a[idxLo] = a[idxMid];
            a[idxMid] = swap;
         
            descendingSort(a, idxLo, idxMid - 1);

            descendingSort(a, idxMid + 1, idxHi);
        }
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive, by ascending values and
     * perform the same operations on b.  Uses the optimized
     * qsort3 from the book "Programming in Pearls" by Jon Bentley.
     * @param a
     * @param b
     * @param idxLo
     * @param idxHi 
     */
    public static void sortBy1stArg(float[] a, int[] b, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
        if (idxLo < idxHi) {

            float x = a[idxLo];
            int store = idxLo;
            int idxMid = idxHi + 1;

            while (true) {
                do {
                    store++;     
                } while ((store <= idxHi) && (a[store] < x));
                do {
                    idxMid--;
                } while (a[idxMid] > x);
                if (store > idxMid) {
                    break;
                }
                float swap = a[store];
                a[store] = a[idxMid];
                a[idxMid] = swap;
                int swap2 = b[store];
                b[store] = b[idxMid];
                b[idxMid] = swap2;
            }
            float swap = a[idxLo];
            a[idxLo] = a[idxMid];
            a[idxMid] = swap;
            int swap2 = b[idxLo];
            b[idxLo] = b[idxMid];
            b[idxMid] = swap2;
         
            sortBy1stArg(a, b, idxLo, idxMid - 1);

            sortBy1stArg(a, b, idxMid + 1, idxHi);
        }
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive, by ascending values and
     * perform the same operations on b.  Uses the optimized
     * qsort3 from the book "Programming in Pearls" by Jon Bentley.
     * @param a
     * @param b
     * @param idxLo
     * @param idxHi 
     */
    public static void sortBy1stArg(int[] a, int[] b, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
        if (idxLo < idxHi) {

            int x = a[idxLo];
            int store = idxLo;
            int idxMid = idxHi + 1;

            while (true) {
                do {
                    store++;     
                } while ((store <= idxHi) && (a[store] < x));
                do {
                    idxMid--;
                } while (a[idxMid] > x);
                if (store > idxMid) {
                    break;
                }
                int swap = a[store];
                a[store] = a[idxMid];
                a[idxMid] = swap;
                int swap2 = b[store];
                b[store] = b[idxMid];
                b[idxMid] = swap2;
            }
            int swap = a[idxLo];
            a[idxLo] = a[idxMid];
            a[idxMid] = swap;
            int swap2 = b[idxLo];
            b[idxLo] = b[idxMid];
            b[idxMid] = swap2;
         
            sortBy1stArg(a, b, idxLo, idxMid - 1);

            sortBy1stArg(a, b, idxMid + 1, idxHi);
        }
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive.  Uses the optimized
     * qsort3 from the book "Programming in Pearls" by Jon Bentley.
     * @param a
     * @param b
     * @param idxLo
     * @param idxHi 
     */
    public static void sortBy1stArg(int[] a, Object[][] b, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        
        if (a.length < 2) {
            return;
        }
        
        if (idxLo < idxHi) {

            int x = a[idxLo];
            int store = idxLo;
            int idxMid = idxHi + 1;

            while (true) {
                do {
                    store++;     
                } while ((store <= idxHi) && (a[store] < x));
                do {
                    idxMid--;
                } while (a[idxMid] > x);
                if (store > idxMid) {
                    break;
                }
                int swap = a[store];
                a[store] = a[idxMid];
                a[idxMid] = swap;
                
                Object[] bSwap = b[store];
                b[store] = b[idxMid];
                b[idxMid] = bSwap;
            }
            int swap = a[idxLo];
            a[idxLo] = a[idxMid];
            a[idxMid] = swap;
            
            Object[] bSwap = b[idxLo];
            b[idxLo] = b[idxMid];
            b[idxMid] = bSwap;
         
            sortBy1stArg(a, b, idxLo, idxMid - 1);

            sortBy1stArg(a, b, idxMid + 1, idxHi);
        }
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive.  Uses the optimized
     * qsort3 from the book "Programming in Pearls" by Jon Bentley.
     * @param a
     * @param b
     * @param idxLo
     * @param idxHi 
     */
    public static void sortBy1stArg(int[] a, Object[] b, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        
        if (a.length < 2) {
            return;
        }
        
        if (idxLo < idxHi) {

            int x = a[idxLo];
            int store = idxLo;
            int idxMid = idxHi + 1;

            while (true) {
                do {
                    store++;     
                } while ((store <= idxHi) && (a[store] < x));
                do {
                    idxMid--;
                } while (a[idxMid] > x);
                if (store > idxMid) {
                    break;
                }
                int swap = a[store];
                a[store] = a[idxMid];
                a[idxMid] = swap;
                
                Object bSwap = b[store];
                b[store] = b[idxMid];
                b[idxMid] = bSwap;
            }
            int swap = a[idxLo];
            a[idxLo] = a[idxMid];
            a[idxMid] = swap;
            
            Object bSwap = b[idxLo];
            b[idxLo] = b[idxMid];
            b[idxMid] = bSwap;
         
            sortBy1stArg(a, b, idxLo, idxMid - 1);

            sortBy1stArg(a, b, idxMid + 1, idxHi);
        }
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive.  Uses the optimized
     * qsort3 from the book "Programming in Pearls" by Jon Bentley.
     * @param a
     * @param b
     * @param idxLo
     * @param idxHi 
     */
    public static void sortBy1stArg(float[] a, Object[] b, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        
        if (a.length < 2) {
            return;
        }
        
        if (idxLo < idxHi) {

            float x = a[idxLo];
            int store = idxLo;
            int idxMid = idxHi + 1;

            while (true) {
                do {
                    store++;     
                } while ((store <= idxHi) && (a[store] < x));
                do {
                    idxMid--;
                } while (a[idxMid] > x);
                if (store > idxMid) {
                    break;
                }
                float swap = a[store];
                a[store] = a[idxMid];
                a[idxMid] = swap;
                
                Object bSwap = b[store];
                b[store] = b[idxMid];
                b[idxMid] = bSwap;
            }
            float swap = a[idxLo];
            a[idxLo] = a[idxMid];
            a[idxMid] = swap;
            
            Object bSwap = b[idxLo];
            b[idxLo] = b[idxMid];
            b[idxMid] = bSwap;
         
            sortBy1stArg(a, b, idxLo, idxMid - 1);

            sortBy1stArg(a, b, idxMid + 1, idxHi);
        }
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive.
     * It's an adaption of the optimized
     * qsort3 from the book "Programming in Pearls" by Jon Bentley.
     * @param a
     * @param b an array that will receive the same swap operations as are 
     * performed on a
     * @param c an array that will receive the same swap operations as are 
     * performed on a
     * @param idxLo
     * @param idxHi 
     */
    public static void sort(float[] a, int[] b, int[] c, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (c == null) {
            throw new IllegalArgumentException("c cannot be null");
        }
        if ((a.length != b.length) || (a.length != c.length)) {
            throw new IllegalArgumentException("array lengths must be the same");
        }
        
        if (idxLo < idxHi) {

            float x = a[idxLo];
            int store = idxLo;
            int idxMid = idxHi + 1;

            while (true) {
                do {
                    store++;     
                } while ((store <= idxHi) && (a[store] < x));
                do {
                    idxMid--;
                } while (a[idxMid] > x);
                if (store > idxMid) {
                    break;
                }
                float swap = a[store];
                a[store] = a[idxMid];
                a[idxMid] = swap;
                int swap2 = b[store];
                b[store] = b[idxMid];
                b[idxMid] = swap2;
                swap2 = c[store];
                c[store] = c[idxMid];
                c[idxMid] = swap2;
            }
            float swap = a[idxLo];
            a[idxLo] = a[idxMid];
            a[idxMid] = swap;
            int swap2 = b[idxLo];
            b[idxLo] = b[idxMid];
            b[idxMid] = swap2;
            swap2 = c[idxLo];
            c[idxLo] = c[idxMid];
            c[idxMid] = swap2;
                     
            sort(a, b, c, idxLo, idxMid - 1);

            sort(a, b, c, idxMid + 1, idxHi);
        }
    }
    
    /**
     * sorts along [0][index] and if there is a tie the value [1][index] is
     * used to decide order.
     * @param a 
     */
    public static void sortByDimension1FirstSecond(int[][] a) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length != 2) {
            throw new IllegalArgumentException("a first dimension length must be 2");
        }
        if (a[0].length < 2) {
            return;
        }
        sortByDimension1FirstSecond(a, 0, a[0].length - 1);
    }
    
    /**
     * sorts along [0][index] and if there is a tie the value [1][index] is
     * used to decide order.
     * @param a 
     * @param idxLo the first index in [0][index] to sort
     * @param idxHi the last index in [0][index] to sort, inclusive
     */
    public static void sortByDimension1FirstSecond(int[][] a, int idxLo, int idxHi) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length != 2) {
            throw new IllegalArgumentException("a first dimension length must be 2");
        }
        if (a[0].length < 2) {
            return;
        }
        if (idxLo < idxHi) {
            int idxMid = partitionByDimension1FirstSecond(a, idxLo, idxHi);
            sortByDimension1FirstSecond(a, idxLo, idxMid - 1);
            sortByDimension1FirstSecond(a, idxMid + 1, idxHi);
        }
    }
    
    private static int partitionByDimension1FirstSecond(int[][] a, int idxLo, 
        int idxHi) {
        
        int x = a[0][idxHi];
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            boolean doSwap = false;
            if (a[0][i] < x) {
                doSwap = true;
            } else if (a[0][i] == x) {
                if (a[1][i] <= a[1][idxHi]) {
                    doSwap = true;
                }
            }
            if (doSwap) {
                store++;
                int swap = a[0][store];
                a[0][store] = a[0][i];
                a[0][i] = swap;
                
                swap = a[1][store];
                a[1][store] = a[1][i];
                a[1][i] = swap;
            }
        }
        store++;
        int swap = a[0][store];
        a[0][store] = a[0][idxHi];
        a[0][idxHi] = swap;
        
        swap = a[1][store];
        a[1][store] = a[1][idxHi];
        a[1][idxHi] = swap;
        return store;
    }
    
    /**
     * sorts along [0][index] and if there is a tie the value [1][index] is
     * used to decide order.
     * @param a 
     */
    public static void sortByDimension1FirstSecondThird(int[][] a) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length != 3) {
            throw new IllegalArgumentException("a first dimension length must be 3");
        }
        if (a[0].length < 2) {
            return;
        }
        sortByDimension1FirstSecondThird(a, 0, a[0].length - 1);
    }
    
    /**
     * sorts along [0][index] and if there is a tie the value [1][index] is
     * used to decide order.
     * @param a 
     * @param idxLo the first index in [0][index] to sort
     * @param idxHi the last index in [0][index] to sort, inclusive
     */
    public static void sortByDimension1FirstSecondThird(int[][] a, int idxLo, int idxHi) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length != 3) {
            throw new IllegalArgumentException("a first dimension length must be 3");
        }
        if (a[0].length < 2) {
            return;
        }
        if (idxLo < idxHi) {
            int idxMid = partitionByDimension1FirstSecondThird(a, idxLo, idxHi);
            sortByDimension1FirstSecondThird(a, idxLo, idxMid - 1);
            sortByDimension1FirstSecondThird(a, idxMid + 1, idxHi);
        }
    }
    
    private static int partitionByDimension1FirstSecondThird(int[][] a, int idxLo, 
        int idxHi) {
        
        int x = a[0][idxHi];
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            boolean doSwap = false;
            if (a[0][i] < x) {
                doSwap = true;
            } else if (a[0][i] == x) {
                if (a[1][i] < a[1][idxHi]) {
                    doSwap = true;
                } else if (a[1][i] == x) {
                    if (a[2][i] <= a[2][idxHi]) {
                        doSwap = true;
                    }
                }
            }
            if (doSwap) {
                store++;
                for (int k = 0; k < 3; ++k) {
                    int swap = a[k][store];
                    a[k][store] = a[k][i];
                    a[k][i] = swap;
                }
            }
        }
        store++;
        for (int k = 0; k < 3; ++k) {
            int swap = a[k][store];
            a[k][store] = a[k][idxHi];
            a[k][idxHi] = swap;
        }
        return store;
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive, with next sorting by b and c
     * and all swap operations performed on all 3 arrays.  
     * @param a
     * @param b
     * @param c
     * @param idxLo
     * @param idxHi 
     */
    public static void sortBy1stThen2ndThen3rd(float[] a, float[] b, float[] c, 
        int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (c == null) {
            throw new IllegalArgumentException("c cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be same length");
        }
        if (a.length != c.length) {
            throw new IllegalArgumentException("a and c must be same length");
        }
        if (idxLo < idxHi) {
            int idxMid = partitionBy1stThen2ndThen3rd(a, b, c, idxLo, idxHi);
            sortBy1stThen2ndThen3rd(a, b, c, idxLo, idxMid - 1);
            sortBy1stThen2ndThen3rd(a, b, c, idxMid + 1, idxHi);
        }
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive, with next sorting by b and c
     * and all swap operations performed on all 3 arrays.  
     * @param a
     * @param b
     * @param c
     * @param idxLo
     * @param idxHi 
     */
    public static void sortBy1stThen2ndThen3rd(
        TIntList a, TIntList b, TIntList c) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.size() < 2) {
            return;
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (c == null) {
            throw new IllegalArgumentException("c cannot be null");
        }
        if (a.size() != b.size()) {
            throw new IllegalArgumentException("a and b must be same length");
        }
        if (a.size() != c.size()) {
            throw new IllegalArgumentException("a and c must be same length");
        }
        sortBy1stThen2ndThen3rd(a, b, c, 0, a.size() - 1);
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive, with next sorting by b and c
     * and all swap operations performed on all 3 arrays.  
     * @param a
     * @param b
     * @param c
     * @param idxLo
     * @param idxHi 
     */
    public static void sortBy1stThen2ndThen3rd(
        TIntList a, TIntList b, TIntList c, 
        int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.size() < 2) {
            return;
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (c == null) {
            throw new IllegalArgumentException("c cannot be null");
        }
        if (a.size() != b.size()) {
            throw new IllegalArgumentException("a and b must be same length");
        }
        if (a.size() != c.size()) {
            throw new IllegalArgumentException("a and c must be same length");
        }
        
        if (idxLo < idxHi) {
            int idxMid = partitionBy1stThen2ndThen3rd(a, b, c, idxLo, idxHi);
            sortBy1stThen2ndThen3rd(a, b, c, idxLo, idxMid - 1);
            sortBy1stThen2ndThen3rd(a, b, c, idxMid + 1, idxHi);
        }
    }
    
    public static void sortBy1stThen2nd(float[] a, float[] b) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be same length");
        }
        sortBy1stThen2nd(a, b, 0, a.length - 1);
    }
    
    public static void sortBy1stThen2nd(int[] a, int[] b) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be same length");
        }
        sortBy1stThen2nd(a, b, 0, a.length - 1);
    }
    
    public static <T extends PairInt> void sortByYThenX(T[] a) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
        sortByYThenX(a, 0, a.length - 1);
    }
    
    public static <T extends PairInt> void 
    sortByDecrYThenIncrX(T[] a, int[] b) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be same length");
        }
        
        sortByDecrYThenIncrX(a, b, 0, a.length - 1);
    }
    
    public static <T extends PairInt> void 
    sortByDecrYThenIncrX(T[] a, int[] b, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be same length");
        }
        
        if (idxLo < idxHi) {
            int idxMid = partitionByDecrYThenIncrX(a, b, idxLo, idxHi);
            sortByDecrYThenIncrX(a, b, idxLo, idxMid - 1);
            sortByDecrYThenIncrX(a, b, idxMid + 1, idxHi);
        }
    }
    
    public static <T extends PairInt> void sortByYThenX(T[] a, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
        
        if (idxLo < idxHi) {
            int idxMid = partitionByYThenX(a, idxLo, idxHi);
            sortByYThenX(a, idxLo, idxMid - 1);
            sortByYThenX(a, idxMid + 1, idxHi);
        }
    }
    
    public static void sortBy1stThen2nd(float[] a, float[] b, int idxLo, int idxHi) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be same length");
        }
        if (idxLo < idxHi) {
            int idxMid = partitionBy1stThen2nd(a, b, idxLo, idxHi);
            sortBy1stThen2nd(a, b, idxLo, idxMid - 1);
            sortBy1stThen2nd(a, b, idxMid + 1, idxHi);
        }
    }
    
    public static void sortBy1stThen2nd(int[] a, int[] b, int idxLo, int idxHi) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be same length");
        }
        if (idxLo < idxHi) {
            int idxMid = partitionBy1stThen2nd(a, b, idxLo, idxHi);
            sortBy1stThen2nd(a, b, idxLo, idxMid - 1);
            sortBy1stThen2nd(a, b, idxMid + 1, idxHi);
        }
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive, with ties sorted by b
     * and all swap operations performed on all arrays. The sorts are
     * ascending.
     * @param a
     * @param b
     * @param c
     * @param d
     * @param idxLo
     * @param idxHi 
     */
    public static void sortBy1stThen2nd(double[] a, double[] b, int[] c, int[] d, 
        int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (c == null) {
            throw new IllegalArgumentException("c cannot be null");
        }
        if (d == null) {
            throw new IllegalArgumentException("d cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be same length");
        }
        if (a.length != c.length) {
            throw new IllegalArgumentException("a and c must be same length");
        }
        if (a.length != d.length) {
            throw new IllegalArgumentException("a and d must be same length");
        }
        if (idxLo < idxHi) {
            int idxMid = partitionBy1stThen2nd(a, b, c, d, idxLo, idxHi);
            sortBy1stThen2nd(a, b, c, d, idxLo, idxMid - 1);
            sortBy1stThen2nd(a, b, c, d, idxMid + 1, idxHi);
        }
    }
    
    public static void sort(float[] a, float[] b, float[] c, int idxLo, 
        int idxHi) {
        
        if (idxLo < idxHi) {
            int idxMid = partition(a, b, c, idxLo, idxHi);
            sort(a, b, c, idxLo, idxMid - 1);
            sort(a, b, c, idxMid + 1, idxHi);
        }
    }

    private static int partition(float[] a, float[] b, float[] c, int idxLo, 
        int idxHi) {
        
        float x = a[idxHi];
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            if (a[i] <= x) {
                store++;
                float swap = a[store];
                a[store] = a[i];
                a[i] = swap;
                float swap2 = b[store];
                b[store] = b[i];
                b[i] = swap2;
                swap2 = c[store];
                c[store] = c[i];
                c[i] = swap2;
            }
        }
        store++;
        float swap = a[store];
        a[store] = a[idxHi];
        a[idxHi] = swap;
        float swap2 = b[store];
        b[store] = b[idxHi];
        b[idxHi] = swap2;
        swap2 = c[store];
        c[store] = c[idxHi];
        c[idxHi] = swap2;
        return store;
    }

    private static int partitionBy1stThen2ndThen3rd(float[] a, float[] b, 
        float[] c, int idxLo, int idxHi) {
        
        float x = a[idxHi];
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            boolean doSwap = false;
            if (a[i] < x) {
                doSwap = true;
            } else if (a[i] == x) {
                if (b[i] < b[idxHi]) {
                    doSwap = true;
                } else if (b[i] == b[idxHi]) {
                    if (c[i] <= c[idxHi]) {
                        doSwap = true;
                    }
                }
            }
            if (doSwap) {
                store++;
                float swap = a[store];
                a[store] = a[i];
                a[i] = swap;
                float swap2 = b[store];
                b[store] = b[i];
                b[i] = swap2;
                swap2 = c[store];
                c[store] = c[i];
                c[i] = swap2;
            }
        }
        store++;
        float swap = a[store];
        a[store] = a[idxHi];
        a[idxHi] = swap;
        float swap2 = b[store];
        b[store] = b[idxHi];
        b[idxHi] = swap2;
        swap2 = c[store];
        c[store] = c[idxHi];
        c[idxHi] = swap2;
        return store;
    }
    
    private static int partitionBy1stThen2ndThen3rd(
        TIntList a, TIntList b, 
        TIntList c, int idxLo, int idxHi) {
        
        int x = a.get(idxHi);
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            boolean doSwap = false;
            if (a.get(i) < x) {
                doSwap = true;
            } else if (a.get(i) == x) {
                if (b.get(i) < b.get(idxHi)) {
                    doSwap = true;
                } else if (b.get(i) == b.get(idxHi)) {
                    if (c.get(i) <= c.get(idxHi)) {
                        doSwap = true;
                    }
                }
            }
            if (doSwap) {
                store++;
                int swap = a.get(store);
                a.set(store, a.get(i));
                a.set(i, swap);
                int swap2 = b.get(store);
                b.set(store, b.get(i));
                b.set(i, swap2);
                swap2 = c.get(store);
                c.set(store, c.get(i));
                c.set(i, swap2);
            }
        }
        store++;
        int swap = a.get(store);
        a.set(store, a.get(idxHi));
        a.set(idxHi, swap);
        int swap2 = b.get(store);
        b.set(store, b.get(idxHi));
        b.set(idxHi, swap2);
        swap2 = c.get(store);
        c.set(store, c.get(idxHi));
        c.set(idxHi, swap2);
        return store;
    }
    
    private static int partitionBy1stThen2nd(float[] a, float[] b,
        int idxLo, int idxHi) {
        
        float x = a[idxHi];
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            boolean doSwap = false;
            if (a[i] < x) {
                doSwap = true;
            } else if (a[i] == x) {
                if (b[i] < b[idxHi]) {
                    doSwap = true;
                }
            }
            if (doSwap) {
                store++;
                float swap = a[store];
                a[store] = a[i];
                a[i] = swap;
                float swap2 = b[store];
                b[store] = b[i];
                b[i] = swap2;
            }
        }
        store++;
        float swap = a[store];
        a[store] = a[idxHi];
        a[idxHi] = swap;
        float swap2 = b[store];
        b[store] = b[idxHi];
        b[idxHi] = swap2;
        return store;
    }
    
    private static int partitionBy1stThen2nd(int[] a, int[] b,
        int idxLo, int idxHi) {
        
        int x = a[idxHi];
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            boolean doSwap = false;
            if (a[i] < x) {
                doSwap = true;
            } else if (a[i] == x) {
                if (b[i] < b[idxHi]) {
                    doSwap = true;
                }
            }
            if (doSwap) {
                store++;
                int swap = a[store];
                a[store] = a[i];
                a[i] = swap;
                int swap2 = b[store];
                b[store] = b[i];
                b[i] = swap2;
            }
        }
        store++;
        int swap = a[store];
        a[store] = a[idxHi];
        a[idxHi] = swap;
        int swap2 = b[store];
        b[store] = b[idxHi];
        b[idxHi] = swap2;
        
        return store;
    }
    
    private static <T extends PairInt> int partitionByYThenX(T[] a, int idxLo, 
        int idxHi) {
     
        T x = a[idxHi];
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            boolean doSwap = false;
            if (a[i].getY() < x.getY()) {
                doSwap = true;
            } else if (a[i].getY() == x.getY()) {
                if (a[i].getX() < x.getX()) {
                    doSwap = true;
                }
            }
            if (doSwap) {
                store++;
                T swap = a[store];
                a[store] = a[i];
                a[i] = swap;
            }
        }
        store++;
        T swap = a[store];
        a[store] = a[idxHi];
        a[idxHi] = swap;
        return store;
    }
    
    private static <T extends PairInt> int 
    partitionByDecrYThenIncrX(T[] a, int[] b, int idxLo, int idxHi) {
     
        T x = a[idxHi];
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            boolean doSwap = false;
            if (a[i].getY() > x.getY()) {
                doSwap = true;
            } else if (a[i].getY() == x.getY()) {
                if (a[i].getX() < x.getX()) {
                    doSwap = true;
                }
            }
            if (doSwap) {
                store++;
                T swap = a[store];
                a[store] = a[i];
                a[i] = swap;
                int swap2 = b[store];
                b[store] = b[i];
                b[i] = swap2;
            }
        }
        store++;
        T swap = a[store];
        a[store] = a[idxHi];
        a[idxHi] = swap;
        int swap2 = b[store];
        b[store] = b[idxHi];
        b[idxHi] = swap2;
        return store;
    }
    
    private static int partitionBy1stThen2nd(double[] a, double[] b, int[] c,
        int[] d, int idxLo, int idxHi) {
        
        double x = a[idxHi];
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            boolean doSwap = false;
            if (a[i] < x) {
                doSwap = true;
            } else if (a[i] == x) {
                if (b[i] < b[idxHi]) {
                    doSwap = true;
                }
            }
            if (doSwap) {
                store++;
                double swap = a[store];
                a[store] = a[i];
                a[i] = swap;
                swap = b[store];
                b[store] = b[i];
                b[i] = swap;
                int swap2 = c[store];
                c[store] = c[i];
                c[i] = swap2;
                swap2 = d[store];
                d[store] = d[i];
                d[i] = swap2;
            }
        }
        store++;
        double swap = a[store];
        a[store] = a[idxHi];
        a[idxHi] = swap;
        swap = b[store];
        b[store] = b[idxHi];
        b[idxHi] = swap;
        int swap2 = c[store];
        c[store] = c[idxHi];
        c[idxHi] = swap2;
        swap2 = d[store];
        d[store] = d[idxHi];
        d[idxHi] = swap2;        
        return store;
    }

}
