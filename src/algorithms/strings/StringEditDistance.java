package algorithms.strings;

import algorithms.util.PairIntArray;
import gnu.trove.list.TCharList;
import gnu.trove.list.array.TCharArrayList;
import gnu.trove.list.linked.TCharLinkedList;
import java.util.Arrays;

/**
 * find the number of insert, delete, and substitution operations to change
 * one string into another.

   note, may add a method for normalized edit distance in future:
     "An Efficient Uniform-Cost Normalized Edit Distance Algorithm"
      by Arslan and Egecioglu, 1999
      https://www.researchgate.net/publication/3820090_An_Efficient_Uniform-Cost_Normalized_Edit_Distance_Algorithm
 * 
 * NOTE: the alignment algorithms could be sped up to a runtime of O(N^2/log^2(N))
 * where N is max of string lengths, bu using
 * Arlazarov, Dinic, Kronrod, Faradzev Four Russian algorithm.
 * 
 * NOTE also that the linear space algorithms can be improved by using
 * "check point" algorithms.  See Optimal alignments in linear space", Myers and Miller, 1988.
 * @author nichole
 */
public class StringEditDistance {
   
    /**
     * calculate the number of insert, delete, and substitution operations to change
     * one string into another following Wagner and Fischer 1954 
     * "The String-to-String Correction Problem".
     * This uses dynamic programming with time O(a.length * b.length)
     * and space O(a.length * b.length).
     * http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.367.5281&rep=rep1&type=pdf
     * 
     * Also, see Levenshtein Distance:
     * https://en.wikipedia.org/wiki/Levenshtein_distance
     * 
     * see also https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
     * and
     * Vintsyuk TK (1968). "Speech discrimination by dynamic programming". Kibernetika. 4: 81–88.
     * 
     * @param a string to edit
     * @param b target string to change a into.
     * 
     * @param outIndexes a and b indexes of the solution in pairs (a_i, b_j)
     * @return 
     */
    public int calculateWithWagnerFischer(String a, String b, PairIntArray outIndexes) {
        
        final int cDel = 1;
        final int cIns = 1;
         /*
        The positive costs and minimization are the common forms in this
        algorithm in context of computer science.
        In biology negative costs and maximization are often used instead
        to represent similarity.
        "Reduced space sequence alignment"
         Grice, Hughey, & Speck, Bioinformatics, Volume 13, Issue 1, February 1997, 
         Pages 45–53, https://doi.org/10.1093/bioinformatics/13.1.45
        */
        
        int cChangeUnequal = 1;
        int m = a.length();
        int n = b.length();
        // (m+1) X (n+1)
        int[][] d = new int[m + 1][];
        for (int i = 0; i <= m; ++i) {
            d[i] = new int[n + 1];
        }
        
        //Theorem 2
        //D (i, j) is the cost of the least cost trace from A (i) to B (j).
        //D(i, j) = min( 
        //               D(i - 1, j - 1) + cost(A(i) -> B(j)), 
        //               D(i- 1,j) + cost(A(i) deleted),
        //               D(i,j -- 1) + cost(insert a character that is  B(j}) )
        //          for all i, j, where 1 <= i <= [A|, 1 <= j <= [B|

        //Theorem 3
        //D(0, 0) = 0; D(i, 0) = summation_over_r=1_to_i(delete of A(r)); 
        //         and D(O, j) = summation_over_r=1_to_j(insert of a character that is B(r))
        // where 1 <= i <= [A|, 1 <= j <= [B|
        
        //Algorithm X
        //D[0,0] := 0;
        //for i := 1 to |A| do D[i, 0] := D[i - 1, 0] + cost_rm(A(i));
        //for j := 1 to |B| do D[O, j] := D[O,j - 1] + cost_insert(B(j));
        //for i := 1 to |A| do
        //  for j := 1 to |B| do begin
        //    m1 := D[i - 1, j - 1] + cost_change(A(i), B(j));
        //    m2 := D[i - 1, j] + cost_rm(A(i));
        //    m3 := D[i, j - 1] + cost_ins(B(j));
        //    D[i, j] := min(m1, m2, m3);
        //    end;
        
        /*
        NOTE: this algorithm is like the LongestCommonSubsequences except that
        for string edit: cIns=cDel=1 and cEqual=0, and the cost is minimized
        for LCS:         cIns=cDel=0 and cEqual=1, and the cost is maximized
        */
        
        for (int i = 1; i <= m; ++i) {
            d[i][0] = i; //i*cDel
        }
        for (int i = 1; i <= n; ++i) {
            //a conversion of empty string 'a' into 'b' is all inserts
            d[0][i] = i; //i*cIns
        }
        
        int i, j;
        int c = 0;
        
        for (j = 1; j <= n; ++j) {
            for (i = 1; i <= m; ++i) {
                if (a.charAt(i-1) == b.charAt(j-1)) {
                    c = 0;
                } else {
                    c = cChangeUnequal;
                }
                d[i][j] = minimum(d[i-1][j-1] + c, d[i-1][j] + cDel, d[i][j-1] + cIns);
            }
        }
        
        // Algorithm Y: print results
        // i := |A|;j := |B|;
        // while(i != 0 && j != 0) do
        //    if D[i, j] = D[i - 1, j] + cost_rm(A(i)) then i := i - 1;
        //    else if D[i, j] = D[i,j - 1] + cost_insert(B(j)) then j :=j- 1;
        //    else begin
        //       print((i, j));
        //       i:=i- 1;j:=j-1;
        //       end;
        i = m;
        j = n;
        while (i != 0 && j != 0){
            if (d[i][j] == d[i-1][j] + cDel) {
                i--;
            } else if (d[i][j] == d[i][j-1] + cIns) {
                j--;
            } else {
                outIndexes.add(i-1, j-1);// subtracting 1 for 0-based indexes
                i--;
                j--;
            }
        }
        //System.out.printf("m=%d, n=%d, d=\n%s\n", m, n, FormatArray.toString(d, "%d"));
        return d[m][n];
    }
    
    /**
     * calculate the number of insert, delete, and substitution operations to change
     * one string into another following Wagner and Fischer 1954 
     * "The String-to-String Correction Problem", but modified to reduce the space complexity.
     * This uses dynamic programming with time O(a.length * b.length)
     * and space O(a.length).
     * <pre>
     * http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.367.5281&rep=rep1&type=pdf
     * 
     * for space complexity reduction, lecture slides of
     * CS 473: Fundamental Algorithms, Spring 2011
     https://courses.engr.illinois.edu/cs473/sp2011/Lectures/09_lec.pdf
     Slide Acknowledgments: Lecture notes are based on notes of Jeff Erickson, 
     Sariel Har-Peled, Chandra Chekuri and the slides of Mahesh Viswanathan and 
     the books by Kleinberg and Tardos, and by Dasgupta, Papadimitrioiu and Vazirani.

     * Also, see Levenshtein Distance:
     * https://en.wikipedia.org/wiki/Levenshtein_distance
     * </pre>
     * @param a string to edit
     * @param b target string to change a into.
     * 
     * @return the edit distance which is the minimum number of deletes of 
     * characters in string a and inserts of characters into string b to turn
     * string a into string b.
     */
    public int calculateWithModifiedWagnerFischer(String a, String b) {
       
        final int cDel = 1;
        final int cIns = 1;
        /*
        The positive costs and minimization are the common forms in this
        algorithm in context of computer science.
        In biology negative costs and maximization are often used instead
        to represent similarity.
        "Reduced space sequence alignment"
         Grice, Hughey, & Speck, Bioinformatics, Volume 13, Issue 1, February 1997, 
         Pages 45–53, https://doi.org/10.1093/bioinformatics/13.1.45
        */
       
        /*
            Entries in jth column only depend on (j − 1)’s column and earlier 
               entries in jth column.
            Only store the current column and the previous column reusing space; 
               N(i, 0) stores M(i, j − 1) and N(i, 1) stores M(i, j)
        */
        
        int cChangeUnequal = 1;
        int m = a.length();
        int n = b.length();
        int[][] dN = new int[m + 1][];
        int i;
        //dN(i, 0) stores M(i, j − 1) and dN(i, 1) stores M(i, j)
        for (i = 0; i <= m; ++i) {
            dN[i] = new int[2];
        }
        // (m+1) X (n+1)
        /*int[][] d = new int[m + 1][];
        for (int i = 0; i <= m; ++i) {
            d[i] = new int[n + 1];
        }*/
        
        //Theorem 2
        //D (i, j) is the cost of the least cost trace from A (i) to B (j).
        //D(i, j) = min( 
        //               D(i - 1, j - 1) + cost(A(i) -> B(j)), 
        //               D(i- 1,j) + cost(A(i) deleted),
        //               D(i,j -- 1) + cost(insert a character that is  B(j}) )
        //          for all i, j, where 1 <= i <= [A|, 1 <= j <= [B|

        //Theorem 3
        //D(0, 0) = 0; D(i, 0) = summation_over_r=1_to_i(delete of A(r)); 
        //         and D(O, j) = summation_over_r=1_to_j(insert of a character that is B(r))
        // where 1 <= i <= [A|, 1 <= j <= [B|
        
        //Algorithm X
        //D[0,0] := 0;
        //for i := 1 to |A| do D[i, 0] := D[i - 1, 0] + cost_rm(A(i));
        //for j := 1 to |B| do D[O, j] := D[O,j - 1] + cost_insert(B(j));
        //for i := 1 to |A| do
        //  for j := 1 to |B| do begin
        //    m1 := D[i - 1, j - 1] + cost_change(A(i), B(j));
        //    m2 := D[i - 1, j] + cost_rm(A(i));
        //    m3 := D[i, j - 1] + cost_ins(B(j));
        //    D[i, j] := min(m1, m2, m3);
        //    end;
        
        /*
        NOTE: this algorithm is like the LongestCommonSubsequences except that
        for string edit: cIns=cDel=1 and cEqual=0, and the cost is minimized
        for LCS:         cIns=cDel=0 and cEqual=1, and the cost is maximized
        */
        
        //dN(i, 0) stores M(i, j − 1) and dN(i, 1) stores M(i, j)
        for (i = 1; i <= m; ++i) {
            //d[i][0] = i; //i*cDel
            dN[i][0] = i; //i*cDel
        }
        //for (i = 1; i <= n; ++i) {
            //d[0][i] = i; //i*cIns
        //}
        
        int j;
        int c = 0;
        
        int minPos;
                
        for (j = 1; j <= n; ++j) {
            dN[0][1] = j; //j*cIns
            //System.out.printf("j=%d: ", j-1);
            for (i = 1; i <= m; ++i) {
                if (a.charAt(i-1) == b.charAt(j-1)) {
                    c = 0;
                    //System.out.printf("(%d,%d),  ", i-1, j-1); //subtracting 1 for 0-based indexes
                } else {
                    c = cChangeUnequal;
                }
                //d[i][j] = minimum(d[i-1][j-1] + c, d[i-1][j] + cDel, d[i][j-1] + cIns);
                minPos = minimumPos(dN[i-1][0] + c, dN[i-1][1] + cDel, dN[i][0] + cIns);
                switch(minPos) {
                    case 1: {
                        dN[i][1] = dN[i-1][1] + cDel;
                        break;
                    } case 2: {
                        dN[i][1] = dN[i][0] + cIns;
                        break;
                    } default: {
                        dN[i][1] = dN[i-1][0] + c; 
                        break;
                    }
                }
            }
            //System.out.printf("\n, ");
            
            // for i = 1 to m do Copy N[i, 0] = N[i, 1]
            for (i = 1; i <= m; ++i) {
                dN[i][0] = dN[i][1];
            }
        }
         
        /*
        can use checkpoint methods to recover the alignment in O(n) space where
        n is the length of the longest string:
        
        "Reduced space sequence alignment"
         Grice, Hughey, & Speck, Bioinformatics, Volume 13, Issue 1, February 1997, 
         Pages 45–53, https://doi.org/10.1093/bioinformatics/13.1.45
        
        The Ukkonen algorithm with the Powell checkpoint method can be used
        to recover alignment in addition to the edit distance using
        O(N) space:
        "A versatile divide and conquer technique for optimal string alignment"
         Powell, Allison, and Dix, Information Processing Letters 70 (1999) 127–139 
        
        alignment can be global, semi-lobal, local, optimal, approximate, 
        pairwise, mutlple sequence, etc.
        
        A 2004 review of use of it in biological context can be found the book 
        "Bioinformatics Sequence and Genome
        Analysis" by David Mount.
        */
        
        return dN[m][1];
    }
    
        /**
     * calculate the number of insert, delete, and substitution operations to change
     * one string into another as an implementation of the Needleman-Wunsch
     * scoring. (The scoring is a modification of the modified Wagner-Fischer algorithm.)
     * <pre>
     * The pseudocode is at https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm
     * see also
     * LongestCommonSubsequence.java
     * </pre>
     * @param a string to edit
     * @param b target string to change a into.
     * 
     * @return the edit distance which is the minimum number of deletes of 
     * characters in string a and inserts of characters into string b to turn
     * string a into string b.
     */
    public int calculateNeedlemanWunschScoring(char[] a, char[] b) {
       
        final int cDel = -2;
        final int cIns = -2;
        final int subUnequal = -1;
        final int subEqual = 2;
       
        /*
            Entries in jth column only depend on (j − 1)’s column and earlier 
               entries in jth column.
            Only store the current column and the previous column reusing space; 
               N(i, 0) stores M(i, j − 1) and N(i, 1) stores M(i, j)
        */
        
        int m = a.length;
        int n = b.length;
        int[][] dN = new int[m + 1][];
        int i;
        //dN(i, 0) stores M(i, j − 1) and dN(i, 1) stores M(i, j)
        for (i = 0; i <= m; ++i) {
            dN[i] = new int[2];
        }
        for (i = 1; i <= m; ++i) {
            dN[i][0] = i*cDel;
        }
        
        int j;
        int c = 0;
                        
        for (j = 1; j <= n; ++j) {
            dN[0][1] = j*cIns;
            for (i = 1; i <= m; ++i) {
                if (a[i-1] == b[j-1]) {
                    c = subEqual;
                    //System.out.printf("(%d,%d),  ", i-1, j-1); //subtracting 1 for 0-based indexes
                } else {
                    c = subUnequal;
                }
                dN[i][1] = maximum(dN[i-1][0] + c, dN[i-1][1] + cDel, dN[i][0] + cIns);
            }
            //System.out.printf("\n, ");
            
            // for i = 1 to m do Copy N[i, 0] = N[i, 1]
            for (i = 1; i <= m; ++i) {
                dN[i][0] = dN[i][1];
            }
        }
        
        return dN[m][1];
    }
    
    /**
     * find the optimal string metric distance between 2 strings using
     * Hirschberg's algorithm.  The scoring uses Levenshtein distance, 
     * defined to be the sum of the costs of insertions, replacements, deletions, 
     * and null actions needed to change one string into the other. 
     * <pre>
     * https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm
     * </pre>
     * @param a
     * @param b
     * @return 
     */
    public TCharList[] hirschbergOptimal(char[] a, char[] b) {
                
        TCharList z = new TCharArrayList();
        TCharList w = new TCharLinkedList();
        
        int m = a.length;
        int n = b.length;
        int i;
        if (m == 0) {
            for (i = 0; i < n; ++i) {
                z.add('-');
                w.add(b[i]);
            }
            return new TCharList[]{z, w};
        } else if (n == 0) {
            for (i = 0; i < m; ++i) {
                z.add(a[i]);
                w.add('-');
            }
            return new TCharList[]{z, w};
        } else if (m==1 || n==1) {
            PairIntArray outIndexes = new PairIntArray();
            int nEdits = calculateWithWagnerFischer(new String(a), new String(b), outIndexes);
            char[] resultA = new char[outIndexes.getN()];
            char[] resultB = new char[outIndexes.getN()];
            for (int ii = 0; ii < resultA.length; ++ii) {
                resultA[ii] = a[outIndexes.getX(ii)];
                resultB[ii] = b[outIndexes.getY(ii)];
            }
            z.addAll(resultA);
            w.addAll(resultB);
            return new TCharList[]{z, w};
        }
        //m: xlen = length(X)
        int aMid = m/2;
        //n: ylen = length(Y)
        //ScoreL = NWScore(X1:xmid, Y)
        //int scoreL = calculateNeedlemanWunschScoring(Arrays.copyOfRange(a, 0, aMid), b);
        //ScoreR = NWScore(rev(Xxmid+1:xlen), rev(Y))
        //int scoreR = calculateNeedlemanWunschScoring(
        //    reverse2(Arrays.copyOfRange(a, aMid + 1, m)), reverse2(b));
        //ymid = arg max ScoreL + rev(ScoreR)
        /* "Optimal alignments in linear space", Myers and Miller, 1988
        Let rev(A) denote the reverse of A, i.e. a_{M-1}a_M...a_0
        and A^T denote the suffix a_{i+1}a_{i+2}...a_M of A. 
        http://www.cs.tau.ac.il/~rshamir/algmb/98/scribe/html/lec02/node9.html
        find k in 0<=k<=n that maximizes D(m,n) = {D(m/2, k) + D^T(m/2, n-k) }
            where 
            D(m/2, k) is D(string1[1:m/2], string2[1:k])
            = D( a.substring(0, aMid), b.substring(0, k) )
            and D^T(m/2, n-k) is 
              D( reversed string1[(m/2)+1:m], reversed string2[k+1:n])
            = D( reverse( a.substring(aMid + 1, m)),reverse( b.substring(k + 1, n)) )
        bMid = ymid = k;
        */
        int bMid = argmax(aMid, a, b);
                
        //(Z,W) = Hirschberg(X1:xmid, y1:ymid) + Hirschberg(Xxmid+1:xlen, Yymid+1:ylen)
        
        TCharList[] zW0 = hirschbergOptimal(
            Arrays.copyOfRange(a, 0, aMid), Arrays.copyOfRange(b, 0, bMid));
        TCharList[] zW1 = hirschbergOptimal( 
            Arrays.copyOfRange(a, aMid, m), Arrays.copyOfRange(b, bMid, n));
        
        TCharList[] zW = new TCharList[2];
        zW[0] = new TCharArrayList(zW0[0]);
        zW[0].addAll(zW1[0]);
        zW[1] = new TCharArrayList(zW0[1]);
        zW[1].addAll(zW1[1]);
        
        return zW;
    }
    
    private String reverse(String a) {
        char[] r = a.toCharArray();
        reverse(r);
        return new String(r);
    }
    private void reverse(char[] a) {
        int n = a.length;
        int i, i2, mid = n/2;
        char swap;
        for (i = 0, i2=n-1; i < mid; ++i, i2--) {
            swap = a[i];
            a[i] = a[i2];
            a[i2] = swap;
        }
    }
    private char[] reverse2(char[] a) {
        char[] out = Arrays.copyOf(a, a.length);
        reverse(out);
        return out;
    }
    
    private int minimum(int a, int b, int c) {
        if (a <= b) {
            if (a <= c) {
                return a;
            } else {
                return c;
            }
        } else {
            // b < a
            if (b <= c) {
                return b;
            } else {
                return c;
            }
        } 
    }
    private int minimumPos(int a, int b, int c) {
        if (a <= b) {
            if (a <= c) {
                return 0;
            } else {
                return 2;
            }
        } else {
            // b < a
            if (b <= c) {
                return 1;
            } else {
                return 2;
            }
        } 
    }
    
    private int maximum(int a, int b, int c) {
        if (a >= b) {
            if (a >= c) {
                return a;
            } else {
                return c;
            }
        } else {
            // b > a
            if (b >= c) {
                return b;
            } else {
                return c;
            }
        } 
    }

    /**
     * find the index in b which maximizes the total score
     *    D( a.substring(0, aMid), b.substring(0, k) )
             + D( reverse( a.substring(aMid + 1, m)), reverse(b.substring(k + 1, n)) ) );
       where D is the score from algorithm NeedlemanWunsch.
     * @param aMid
     * @param a
     * @param b
     * @return 
     */
    private int argmax(int aMid, char[] a, char[] b) {
        /*
        D( a.substring(0, aMid), b.substring(0, k) )
             + D( reverse( a.substring(aMid + 1, m) ), 
                  reverse( b.substring(k + 1, n) ) );
        */
        char[] a0Mid = Arrays.copyOfRange(a, 0, aMid);
        char[] aMidMRev = Arrays.copyOfRange(a, aMid, a.length);
        reverse(aMidMRev);
        
        int maxK = -1;
        int maxVal = Integer.MIN_VALUE;
        int k, sL, sR, tot;
        char[] b0k, bkMidRev;
        for (k = 0; k < b.length; ++k) {
            b0k = Arrays.copyOfRange(b, 0, k);
            bkMidRev = Arrays.copyOfRange(b, k, b.length);
            reverse(bkMidRev);
            
            sL = calculateNeedlemanWunschScoring(a0Mid, b0k);
            sR = calculateNeedlemanWunschScoring(aMidMRev, bkMidRev);
            
            tot = sL + sR;
            if (tot > maxVal) {
                maxVal = tot;
                maxK = k;
            }
        }
        return maxK;
    }
}
