package thirdparty.libquantum;

/* classic.h: Declarations for classic.c

   Copyright 2003 Bjoern Butscher, Hendrik Weimer

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
public class Classic {

    /**
     * Calculate A^B with A and B as integers
     */
    int quantum_ipow(int a, int b) {
        int i;
        int r = 1;

        for (i = 0; i < b; i++) {
            r *= a;
        }

        return r;
    }

    /**
     * Calculate the greatest common divisor with Euclid's algorithm
     */
    public int quantum_gcd(int u, int v) {
        int r;

        while (v != 0) {
            r = u % v;
            u = v;
            v = r;
        }
        return u;
    }

    /**
     * Fractional approximation of a decimal value
    
     * @param aInOut input and output array of size 1
     * @param bInOut input and output array of size 1
     * @param width 
     */
    void quantum_frac_approx(int[] aInOut, int[] bInOut, int width) {

        if (aInOut.length != 1 || bInOut.length != 1) {
            throw new IllegalArgumentException("arrays must be lenth 1");
        }
        
        int a = aInOut[0];
        int b = bInOut[0];

        float f = (float) a / b;
        float g = f;
        int i, num2 = 0, den2 = 1, num1 = 1, den1 = 0, num = 0, den = 0;

        do {
            i = (int) (g + 0.000005);

            g -= i - 0.000005f;
            g = 1.0f / g;

            if (i * den1 + den2 > (1 << width)) {
                break;
            }

            num = i * num1 + num2;
            den = i * den1 + den2;

            num2 = num1;
            den2 = den1;
            num1 = num;
            den1 = den;

        } while (Math.abs(((double) num / den) - f) > 1.0 / (2 * (1 << width)));

        aInOut[0] = num;
        bInOut[0] = den;

        return;
    }

    /**
     * Calculates the number of qubits required to store N
     */
    static int quantum_getwidth(int n) {
        
        int i;

        for (i = 1; 1<<i<n; i++);

        return i;
    }

    /**
     * Calculate the inverse modulus of N and C
     */
    int quantum_inverse_mod(int n, int c) {
        
        if (c == 0) {
            throw new IllegalArgumentException("c cannot == 0");
        }
        
        int i;

        for (i = 1; (i * c) % n != 1; i++);

        return i;
    }
}
