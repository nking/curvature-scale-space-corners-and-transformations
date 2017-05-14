package thirdparty.libquantum;

import algorithms.misc.Misc;
import java.util.Random;

 /* shor.c: Implementation of Shor's factoring algorithm

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

public class Shor {
        
    private QuantumReg qr;
    
    /**
     * the number to factorize
     */
    private final int N;
    
    /**
     * random factor of N
     */
    private int x = 0;
    
    private final Random rng;
    
    private final long rSeed;
    
    public Shor(int number) {
        
        if (number < 15) {
            throw new IllegalArgumentException("Invalid number\n\n");
        }
        
        this.N = number;
        
        this.rSeed = System.nanoTime();
        
        this.rng = Misc.getSecureRandom();
        this.rng.setSeed(rSeed);
    }
    
    public Shor(int number, long rSeed) {
        
        if (number < 15) {
            throw new IllegalArgumentException("Invalid number\n\n");
        }
        
        this.N = number;
        
        this.rSeed = rSeed;
        
        this.rng = Misc.getSecureRandom();
        this.rng.setSeed(rSeed);
    }
    
    /**
     * 
     * @param number
     * @param rSeed
     * @param x a known factor of x to use in initializing the quantum order 
     * finder
     */
    public Shor(int number, int x) {
        
        if (number < 15) {
            throw new IllegalArgumentException("Invalid number\n\n");
        }
        
        this.N = number;
        
        this.rSeed = System.currentTimeMillis();
        
        this.rng = Misc.getSecureRandom();
        this.rng.setSeed(rSeed);
                
        this.x = x;
    }
    
    /**
     * 
     * @param number
     * @param rSeed
     * @param x a known factor of x to use in initializing the quantum order 
     * finder
     */
    public Shor(int number, long rSeed, int x) {
        
        if (number < 15) {
            throw new IllegalArgumentException("Invalid number\n\n");
        }
        
        this.N = number;
        
        this.rSeed = rSeed;
        
        this.rng = Misc.getSecureRandom();
        this.rng.setSeed(rSeed);
        
        this.x = x;
    }
    
    /**
     * 
     * @return returns 2 factors of number, else returns a single item error code. 
     */
    public int[] run() {
            
        int width = Classic.quantum_getwidth(N*N);
        int swidth = Classic.quantum_getwidth(N);
    
        System.out.println("SEED=" + rSeed);
        System.out.format("N = %d, %d qubits required\n", N, width+3*swidth+2);
        
        if (x == 0) {
            Classic classic = new Classic();
            while ((classic.quantum_gcd(N, x) > 1) || (x < 2)) {
                //x = rand() % N;
                x = rng.nextInt((1<<16)-1) % N;
            }
        }
        
        int i;
        int c,q,a,b, factor;

        System.out.format("Random factor: %d of %d\n", x, N);

        QuReg qureg = new QuReg();
        
        QuantumReg qr = qureg.quantum_new_qureg(0, width);
 
        assert(qr.hash.length == QuReg.shiftLeftTruncate(qr.hashw));
        
        Gates gates = new Gates(rng);
       
        for (i = 0; i < width; i++) {
            gates.quantum_hadamard(i, qr);
        }
        
        assert(qr.hash.length == QuReg.shiftLeftTruncate(qr.hashw));
         
        int nbits = 3 * swidth + 2;

        qureg.quantum_addscratch(nbits, qr);

        gates.quantum_exp_mod_n(N, x, width, swidth,  qr);
        
        assert(qr.hash.length == QuReg.shiftLeftTruncate(qr.hashw));
     
        Measure measure = new Measure();
        
        for (i = 0; i < nbits; i++) {
            measure.quantum_bmeasure(0, qr, rng);
        }
        
        assert(qr.hash.length == QuReg.shiftLeftTruncate(qr.hashw));
 
        gates.quantum_qft(width,  qr);

        for (i = 0; i < width / 2; i++) {
            gates.quantum_cnot(i, width - i - 1, qr);
            gates.quantum_cnot(width - i - 1, i, qr);
            gates.quantum_cnot(i, width - i - 1, qr);
        }
        
        assert(qr.hash.length == QuReg.shiftLeftTruncate(qr.hashw));

        c = measure.quantum_measure(qr, rng);

        System.out.println("c=" + c);

        if (c == -1) {
            System.out.format("Impossible Measurement!\n");
            return new int[]{-1};
        }

        if (c == 0) {
            System.out.format("Measured zero, try again.\n");
            return new int[]{0};
        }

        q = 1 << (width);

        System.out.format("Measured %d (%f), ", c, (float) c / q);

        Classic classic = new Classic();
        int[] cInOut = new int[]{c};
        int[] qInOut = new int[]{q};
        classic.quantum_frac_approx(cInOut, qInOut, width);
        c = cInOut[0];
        q = qInOut[0];
        
        System.out.format("fractional approximation is %d/%d.\n", c, q);

        if ((q % 2 == 1) && (2 * q < (1 << width))) {
            System.out.format("Odd denominator, trying to expand by 2.\n");
            q *= 2;
        }

        if (q % 2 == 1) {
            System.out.format("Odd period, try again.\n");
            return new int[]{-2};
        }

        System.out.format("%d\n", 11);

        System.out.format("Possible period is %d.\n", q);

        a = classic.quantum_ipow(x, q / 2) + 1 % N;
        b = classic.quantum_ipow(x, q / 2) - 1 % N;

        a = classic.quantum_gcd(N, a);
        b = classic.quantum_gcd(N, b);

        if (a > b) {
            factor = a;
        } else {
            factor = b;
        }

        if ((factor < N) && (factor > 1)) {
            System.out.format("%d = %d * %d\n", N, factor, N / factor);
        } else {
            System.out.format("Unable to determine factors, try again.\n");
            return new int[]{-2};
        }

        return new int[]{factor, N/factor};
    }
}
