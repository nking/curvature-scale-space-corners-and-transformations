package algorithms.strings;

import algorithms.util.PairIntArray;
import gnu.trove.list.TIntList;
import java.util.Arrays;

/**
 * find the number of insert, delete, and substitution operations to change
 * one string into another.

   note, may add a method for normalized edit distance in future:
     "An Efficient Uniform-Cost Normalized Edit Distance Algorithm"
      by Arslan and Egecioglu, 1999
      https://www.researchgate.net/publication/3820090_An_Efficient_Uniform-Cost_Normalized_Edit_Distance_Algorithm
 * 
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
     * @param a string to edit
     * @param b target string to change a into.
     * 
     * @param outIndexes a and b indexes of the solution in pairs (a_i, b_j)
     * @return 
     */
    public int calculateWithWagnerFischer(String a, String b,
        PairIntArray outIndexes) {
        
        final int cDel = 1;
        final int cIns = 1;

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
                outIndexes.add(i, j);
                i--;
                j--;
            }
        }
        return d[m][n];
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
   
}
