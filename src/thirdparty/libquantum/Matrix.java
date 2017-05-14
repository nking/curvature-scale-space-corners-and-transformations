package thirdparty.libquantum;

import algorithms.misc.ComplexModifiable;

/* matrix.c: Matrix operations

   Copyright 2003, 2005 Bjoern Butscher, Hendrik Weimer

   This file is part of libquantum

   libquantum is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation; either version 3 of the License,
   or (at your option) any later version.

   libquantum is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with libquantum; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
   MA 02110-1301, USA

 */
public class Matrix {

    /**
     * extract the complex number m.t[(x) + (y) * m.cols]
     * without copying.
     * @param m
     * @param x
     * @param y
     * @return 
     */
    public static ComplexModifiable M(QuantumMatrix m, int x, int y) {
        return m.t[(x) + (y) * m.cols];
    }

    /**
     * Create a new COLS x ROWS matrix. 
     * Note that m.t is initialized and filled.
     */
    QuantumMatrix quantum_new_matrix(int cols, int rows) {
        
        QuantumMatrix m = new QuantumMatrix();

        m.rows = rows;
        m.cols = cols;
        int n = cols * rows;
        m.t = new ComplexModifiable[n];
        for (int i = 0; i < n; ++i) {
            m.t[i] = new ComplexModifiable(0, 0);
        }
        
        return m;
    }

    /**
     * Print the contents of a matrix to stdout
     */
    void quantum_print_matrix(QuantumMatrix m) {
        
        int i, j, z = 0;

        while ((1 << z++) < m.rows);
        z--;

        for (i = 0; i < m.rows; i++) {
            for (j = 0; j < m.cols; j++) {
                System.out.format("%g %+gi ", M(m, j, i).re(), M(m, j, i).im());
            }
            System.out.format("\n");
        }
        System.out.format("\n");
    }

    /**
     * Matrix multiplication
     */
    QuantumMatrix quantum_mmult(QuantumMatrix A, QuantumMatrix B) {
        
        int i, j, k;

        if (A.cols != B.rows) {
            throw new IllegalArgumentException("A.cols must == B.rows");
        }

        QuantumMatrix C = quantum_new_matrix(B.cols, A.rows);

        for (i = 0; i < B.cols; i++) {
            for (j = 0; j < A.rows; j++) {
                for (k = 0; k < B.rows; k++) {
                    ComplexModifiable tmp = M(A, k, j).copy();
                    tmp.times(M(B, i, k));
                    M(C, i, j).plus(tmp);
                }
            }
        }

        return C;
    }
}