package algorithms.strings;

import algorithms.util.PairIntArray;
import java.util.Arrays;

/**
 * find the longest common subsequence between strings a and b.

 * 
 * @author nichole
 */
public class LongestCommonSubsequence {
    
    protected final static byte UP = (byte)1;
    protected final static byte UPLEFT = (byte)2;
    protected final static byte LEFT = (byte)3;
    
    /**
     * find the longest overlapping sequence between 2 sequences, allowing gaps
     * following Cormen et al. "Introduction to Algorithms".
     *
     * The solution takes time O(m * n) and uses dynamic programming.
     
     * @param x one string
     * @param y another string
     * @param printToStdOut
     * @return 
     */
    public static char[] calculateWithCormenEtAl(char[] x, char[] y, boolean printToStdOut) {
        
        /*
         * let c[i,j] be the length of a LCS of X_i and Y_j.
         *
         * the optimal solution gives the recursive solution:
         *
         *            { 0                              if i=0 or j = 0
         *   c[i,j] = { c[i-1, j-1] + 1                if i, j > 0 and x_i == y_i
         *            { max( c[i, j-1], c[i-1, j] )    if i, j > 0 and x_i != y_i
         *
         * A table of 'directions' is stored in b[0..m,0...n] to simplify reconstruction
         * of the matching characters, where directions are 'up', 'up and left' or 'left'.
         *
         * Below, will use 'U' for up, 'UL' for 'up and left
         * and 'L' for left, which are represented in
         * the code as bytes 1, 2, and 3 respectively.
        */
        int[][] c;
        byte[][] b;

        int m = x.length;
        int n = y.length;

        c = new int[m+1][n+1];
        b = new byte[m+1][n+1];
        for (int i = 0; i < c.length; i++) {
            c[i] = new int[n+1];
            b[i] = new byte[n+1];
        }

        // these are already 0
        /*
        for (int i = 1; i < xlen; i++) {
            c[i][0] = 0;
        }
        for (int j = 0; j < ylen; j++) {
            c[0][j] = 0;
        }*/

        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                if (x[i-1] == y[j-1]) {
                    c[i][j] = c[i-1][j-1] + 1;
                    b[i][j] = UPLEFT;
                } else if (c[i-1][j] >= c[i][j-1]) {
                    c[i][j] = c[i-1][j];
                    b[i][j] = UP;
                } else {
                    c[i][j] = c[i][j-1];
                    b[i][j] = LEFT;
                }
            }
        }
        
        if (printToStdOut) {
            for (int i = 0; i < c.length; i++) {
                StringBuilder sb = new StringBuilder("");
                for (int j = 0; j < c[0].length; j++) {
                    String d = "   ";
                    if (b[i][j] == UPLEFT) {
                        d = " UL";
                    } else if (b[i][j] == UP) {
                        d = "  U";
                    } else if (b[i][j] == LEFT) {
                        d = "  L";
                    }
                    sb.append(d).append(c[i][j]);
                }
                System.out.println(sb.toString());
            }

        }
        
        int maxLen = (x.length > y.length) ? y.length : x.length;
        char[] lcs = new char[maxLen];

        int len = getLCS(b, x, b.length - 1, b[0].length - 1, lcs, 0);

        return Arrays.copyOf(lcs, len);
    }
    
    /**
     * back track through the matrix to get the letters.  runtime complexity O(m + n).
     * @param b
     * @param x
     * @param row
     * @param col
     * @param str
     * @param currentLength
     * @return 
     */
    private static int getLCS(final byte[][]b, char[]x, int row, int col, char[] str, 
        int currentLength) {

        if (row == 0 || col == 0) {
            return currentLength;
        }

        if (b[row][col] == UPLEFT) {
            currentLength = getLCS(b, x, row-1, col-1, str, currentLength);
            str[currentLength] = x[row-1];
            currentLength++;
        } else if (b[row][col] == UP) {
            currentLength = getLCS(b, x, row-1, col, str, currentLength);
        } else {
            currentLength = getLCS(b, x, row, col-1, str, currentLength);
        }

        return currentLength;
    }
    
}
