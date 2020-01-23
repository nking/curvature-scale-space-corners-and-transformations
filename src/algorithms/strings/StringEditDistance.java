package algorithms.strings;

/**
 * find the number of insert, delete, and substitution operations to change
 * one string into another.
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
     * @return 
     */
    public int calculateWithWagnerFischer(String a, String b, boolean printToStdOut) {
        int m = a.length();
        int n = b.length();
        int[][] d = new int[m + 1][];
        for (int i = 0; i <= m; ++i) {
            d[i] = new int[n + 1];
        }

        //Algorithm X
        //for i := 1 to |A| do D[i, 0] := D[i - 1, 0] + cost_rm(A(i));
        //for j := 1 to |B| do D[O, j] := D[O,j - 1] + cost_insert(B(j));
        //for i := 1 to |A| do
        //  for j := 1 to |B| do begin
        //    m1 := D[i - 1, j - 1] + cost_change(A(i), B(j));
        //    m2 := D[i - 1, j] + cost_rm(A(i));
        //    m3 := D[i, j - 1] + cost_ins(B(j));
        //    D[i, j] := min(m1, m2, m3);
        //    end;
        
        for (int i = 1; i <= m; ++i) {
            d[i][0] = i;
        }
        for (int i = 1; i <= n; ++i) {
            d[0][i] = i;
        }
        
        int i, j;
        int c = 0;
        int cDel = 1;
        int cIns = 1;
        for (j = 1; j <= n; ++j) {
            for (i = 1; i <= m; ++i) {
                if (a.charAt(i-1) == b.charAt(j-1)) {
                    c = 0;
                } else {
                    c = 1;
                }
                d[i][j] = minimum(d[i-1][j-1] + c, d[i-1][j] + cDel, d[i][j-1] + cIns);
            }
        }
        
        if (!printToStdOut) {
            return d[m][n];
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
                System.out.printf("(%d,%d)\n", i, j);
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
