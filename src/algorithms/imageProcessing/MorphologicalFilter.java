package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * other morphological features are located in other classes and might be moved
 * here in the future.
 * 
 * The class holds an implementation  of the Octave bwmorph function for "thin"
 * and the copyright for it is in the method comments.
 * 
 * @author nichole
 */
public class MorphologicalFilter {
    
    // array length is 512
    private static final int[] lut1 = 
              new int[]{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,1,1,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,1,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,0,0,0,0,0,0,0,1,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,
                       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,1,0,0};
    
    private static final int[] lut2
           = new int[]{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,0,1,0,1,1,0,0,0,0,1,0,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,0,0,1,0,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,0,0,0,1,1,0,1,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,0,1,0,1,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,0,1,0,0,1,1,1,1,1,1,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,0,1,1,0,1,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,1,0,0,1,0,1,1,0,1,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,1,0,1,1,1,0,1,0,1,0,1,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,0,1,0,1,1,1,1,1,1,1,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,0,0,1,0,1,1,1,1,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,0,1,0,0,1,0,1,1,1,1,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,
                       1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    
    /**
     * adapted from Octave's bwmorph "thin" algorithm (please see method comments
     * for copyright by Monas and Draug, GNU GPL).
     
      For a binary image bImg, performs the morphological operation,
      nIter times.  By default, nIter is 1.  If nIter is ~ Inf, the
      operation is continually performed until it no longer changes the image.

      Note that the output will always be binary too.

     * @param bImg binary image
     * @param nIter the number of times the algorithm will be performed or
     * fewer if no more changes occur
     * @return 
     */
    public int[][] bwMorphThin(int[][] bImg, int nIter) {
        
        /*
        ## Copyright (C) 2004 Josep Mones i Teixidor <jmones@puntbarra.com>
        ## Copyright (C) 2013 Carnë Draug <carandraug@octave.org>
        ##
        ## This program is free software; you can redistribute it and/or modify it under
        ## the terms of the GNU General Public License as published by the Free Software
        ## Foundation; either version 3 of the License, or (at your option) any later
        ## version.
        ##
        ## This program is distributed in the hope that it will be useful, but WITHOUT
        ## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
        ## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
        ## details.
        ##
        ## You should have received a copy of the GNU General Public License along with
        ## this program; if not, see <http://www.gnu.org/licenses/>.
        */
        
        if (nIter < 1) {
            throw new IllegalArgumentException("nIter must be >= 1");
        }
        
        int[][] bImg2 = copy(bImg);
        int i = 1;
        while (i <= nIter) {
            // morph = @(x) x & applylut (applylut (x, lut1), lut2);
            // int[][] applyLut(int[][] bImg, int[] lut)
            int[][] inner = applyLut(bImg, lut1);
            int[][] outer = applyLut(inner, lut2);
            // bImg & outer
            for (int col = 0; col < bImg.length; ++col) {
                for (int row = 0; row < bImg[col].length; ++row) {
                    bImg2[col][row] = bImg[col][row] & outer[col][row];
                }
            }
            if (isEqual (bImg, bImg2)) {
                //if it doesn't change we don't need to process it further
                 break;
            }
            
            bImg = bImg2;
            i++;
        }
        
        return bImg2;
        
        /*
        applylut performs a 3x3 neighborhood operation on the input with a 512 vector.
        For 3-by-3 neighborhoods, length(lut) is 512. There are nine pixels in 
        each neighborhood, and two possible states for each pixel, so the total 
        number of permutations is 2^9 = 512.

        To produce the matrix of indices, applylut convolves the binary image 
        BW with this matrix.

         256    32     4
         128    16     2
         64     8     1
         The resulting convolution contains integer values in the range [0,511]. 
         applylut uses the central part of the convolution, of the same size as 
         BW, and adds 1 to each value to shift the range to [1,512]. 
         It then constructs A by replacing the values in the cells of the index 
         matrix with the values in lut that the indices point to.
        
      morph = @(x) x & applylut (applylut (x, lut1), lut2);
        
      bw2_tmp = bw; ## make sure bw2_tmp will exist later, even if n == 0
      i = 1;
      while (i <= n) ## a for loop wouldn't work because n can be Inf
        bw2_tmp = morph (bw);
        if (isequal (bw, bw2_tmp))
          ## if it doesn't change we don't need to process it further
          break
        endif
        bw = bw2_tmp;
        i++;
      endwhile
        */
        
    }
    
    /**
     * Uses lookup tables to perform a neighbour operation on binary images.
     * 
     * adapted from Octave
     * https://sourceforge.net/p/octave/image/ci/8fe38c1c25c5ae89ffe8cb5e471cd43e22354847/tree/inst/applylut.m?format=raw
     * which has copyright:
     * Copyright (C) 2004 Josep Mones i Teixidor
    ##
    ## This program is free software; you can redistribute it and/or modify
    ## it under the terms of the GNU General Public License as published by
    ## the Free Software Foundation; either version 2 of the License, or
    ## (at your option) any later version.
    ##
    ## This program is distributed in the hope that it will be useful,
    ## but WITHOUT ANY WARRANTY; without even the implied warranty of
    ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    ## GNU General Public License for more details.
    ##
    ## You should have received a copy of the GNU General Public License
    ## along with this program; if not, write to the Free Software
    ## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

     * @param bImg
     * @param lut 
     */
    private int[][] applyLut(int[][] bImg, int[] lut) {
        
        // for 512, this is 9
        double nq = Math.log(lut.length)/Math.log(2);
        
        // for 512, this is 3
        int n = (int)Math.sqrt(nq);
        
        if (Math.floor(n) != n) {
            throw new IllegalArgumentException("lut is expected to be a power of 2");
        }
        
        //nq-1:-1:0 ->  [8, 7, 6, 5, 4, 3, 2, 1, 0]
        List<Integer> range = new ArrayList<Integer>();
        for (int i = (int)Math.round(nq) - 1; i >= 0; --i) {
            range.add(Integer.valueOf(i));
        }
        // 256, 128, 64, 32, 16, 8, 4, 2, 1
        for (int i = 0; i < range.size(); ++i) {
            int v = range.get(i);
            v = (1 << v);
            range.set(i, v);
        }
        
        //w=reshape(2.^[nq-1:-1:0],n,n);
        //w = reshape([256, 128, 64, 32, 16, 8, 4, 2, 1], 3, 3)
        // [[256, 128, 64].
        //  [32, 16, 8],
        //  [4, 2, 1]]
        int count = 0;
        int[][] w = new int[n][];
        for (int i = 0; i < n; ++i) {
            w[i] = new int[n];
            for (int j = 0; j < n; ++j) {
                w[i][j] = range.get(count).intValue();
                count++;
            }
        }
        
        /*
        A=LUT(filter2(w, BW) + 1);
        */
        int[][] f = filter2(w, bImg);
        
        // replace the values in f with the values of the lut[index]
        for (int col = 0; col < f.length; ++col) {
            for (int row = 0; row < f[col].length; ++row) {
                int idx = f[col][row];
                int v = lut[idx];
                f[col][row] = v;
            }
        }
        
        return f;
    }
    
    /**
     * apply the filter w to image bImg
     * 
     * adapted from Octave
     * https://sourceforge.net/p/octave/signal/ci/b40e74b9814bfcbcd770b9a12a852fd12e611995/tree/filter2.m?format=raw
     * which has copyright:
     * Copyright (C) 2001 Paul Kienzle
        ##
        ## This program is free software; you can redistribute it and/or modify
        ## it under the terms of the GNU General Public License as published by
        ## the Free Software Foundation; either version 2 of the License, or
        ## (at your option) any later version.
        ##
        ## This program is distributed in the hope that it will be useful,
        ## but WITHOUT ANY WARRANTY; without even the implied warranty of
        ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        ## GNU General Public License for more details.
        ##
        ## You should have received a copy of the GNU General Public License
        ## along with this program; if not, write to the Free Software
        ## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
        
        and from
        * https://sourceforge.net/p/octave/image/ci/a99c3050162e76ce511c2f02daf042ae326fe6d3/tree/conv2.cc?format=raw
        * which has copyright:
        * Copyright (C) 1999 Andy Adler
        * This code has no warrany whatsoever.
        * Do what you like with this code as long as you
        *     leave this copyright in place.
        * 
     * @param w
     * @param bImg
     * @return 
     */
    private int[][] filter2(int[][] w, int[][] bImg) {
        
        int nCols = w.length;
        int nRows = w[0].length;
        
        /*
        [nr, nc] = size(w);
        
        where w is
        // [[256, 128, 64],
        //  [32, 16, 8],
        //  [4, 2, 1]]
        nr=3
        nc=3
        w(nr:-1:1, nc:-1:1) -->  w([2, 1, 0], [2, 1, 0])
                        
        Y = conv2(bImg, w(nr:-1:1, nc:-1:1), shape);
        */
        int[][] w180 = new int[nCols][nRows];
        for (int i = 0; i < nCols; ++i) {
            w180[i] = Arrays.copyOf(w[nCols - i - 1], nRows);
        }
        //  [4, 2, 1]]
        //  [32, 16, 8],
        //  [256, 128, 64],
      
        // swap columns within w180[i][j]
        int end = nRows >> 1;
        for (int i = 0; i < nCols; ++i) {
            for (int j = 0; j < end; ++j) {
                int j2 = nRows - j - 1;
                int swap = w180[i][j2];
                w180[i][j2] = w180[i][j];
                w180[i][j] = swap;
            }
        }
        
        //    w180
        // [[ 1,   2,   4],
        //  [ 8,  16,  32],
        //  [64. 128, 256]]        
        
        // the values in y are the indexes of lut
        int[][] y = conv2(bImg, w180);
        
        return y;
    }
    
    /**
     * adapted from 
        * https://sourceforge.net/p/octave/image/ci/a99c3050162e76ce511c2f02daf042ae326fe6d3/tree/conv2.cc?format=raw
        * which has copyright:
        * Copyright (C) 1999 Andy Adler
        * This code has no warrany whatsoever.
        * Do what you like with this code as long as you
        *     leave this copyright in place.
        * 
     * @param w
     * @param bImg
     * @return 
     */
    private int[][] conv2(int[][] a, int[][] b) {
        
        int Am = a.length;
        int An = a[0].length;
        int Bm = b.length;
        int Bn = b[0].length;

        int outM = Am;
        int outN = An;
        int edgM = (Bm - 1)/2;
        int edgN = (Bn - 1)/2;
        
        int[][] output = new int[outM][outN];
        for (int oi = 0; oi < outM; oi++) {
            
            output[oi] = new int[outN];
            
            for (int oj = 0; oj < outN; oj++) {
                
                int sum = 0;
                
                for(int bj = Bn - 1 - Math.max(0, edgN-oj), 
                    aj = Math.max(0, oj-edgN); bj >= 0 && aj < An; bj--, aj++) {
                    
                    int bi = Bm - 1 - Math.max(0, edgM - oi);
                    int ai = Math.max(0, oi - edgM); 
                    
                    for ( ; bi >= 0 && ai < Am; bi--, ai++) {
                        //sum += (*Ad) * (*Bd);    
                        // Comment: it seems to be 2.5 x faster than this:
                        //        sum+= A(ai,aj) * B(bi,bj);
                        sum += a[ai][aj] * b[bi][bj];
                    }
                }
                output[oi][oj] = sum;
            }
        }
        
        return output;
    }

    private boolean isEqual(int[][] a, int[][] b) {
        
        if (a.length != b.length || a[0].length != b[0].length) {
            return false;
        }
        
        for (int i = 0; i < a.length; ++i) {
            if (!Arrays.equals(a[i], b[i])) {
                return false;
            }
        }
        
        return true;
    }

    private int[][] copy(int[][] a) {
     
        int[][] b = new int[a.length][];
        for (int i = 0; i < a.length; ++i) {
            b[i] = Arrays.copyOf(a[i], a[i].length);
        }
        
        return b;
    }
    
}
