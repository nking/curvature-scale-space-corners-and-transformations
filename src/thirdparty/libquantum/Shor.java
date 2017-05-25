package thirdparty.libquantum;

import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import java.util.Random;
import java.util.logging.Logger;

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

/**
 * find the integer factors of a number.
 * N is an odd composite.
 * Also, ensure that N is not the power of a prime.
 * can check by:
 *    for all k .lte. log_2(N) that math.pow(N, k) is not an integer.
 *
 * the libquantum runtime complexity is approx N^2 * log_2(N) * log_2(N),
     but here, have reduced the number of qubits at initialization to
     2^(log2(N)) instead of 2^(log2(N*N)),
        so the runtime complexity is now 
        approx N * log_2(N) * log_2(N)
        
   For larger numbers, you might want to feedback the largest cofactor in a 
   single result into another instance.
    
   The code is limited to signed integers.
    
   Note, that in contrast the general number field sieve integer factorizations,
   has runtime complexity O( exp( ( (64/9)*b*(log b * log b) )^(1/3) ) )
   where b is bit size of N (which is ~log2(N)).
   
 */
public class Shor {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
        
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
    
    private boolean useLargerInit = false;
    
    private boolean retryMeasure0 = false;
    
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
     * @param x a known factor of x to use in initializing the quantum order 
     * finder
     */
    public Shor(int number, int x) {
        
        if (number < 15) {
            throw new IllegalArgumentException("Invalid number\n\n");
        } else if (number > 32768) {
            throw new IllegalArgumentException("number must"
                + "be les than 32768");
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
        } else if (number > 32768) {
            throw new IllegalArgumentException("number must"
                + "be les than 32768");
        }
        
        this.N = number;
        
        this.rSeed = rSeed;
        
        this.rng = Misc.getSecureRandom();
        this.rng.setSeed(rSeed);
        
        this.x = x;
    }
    
    public void overrideToUseLargerInitialization() {
        useLargerInit = true;
    }
    
    public void overrideToRetryMeasured0() {
        retryMeasure0 = true;
    }
    
    /**
     * essentially, makes an array of bitstrings of
       size 2^(N), calculates factors of
       the moduli and applies them, then
       examines the bit spacings by applying
       conditional phase shifts and hadamard gates
       followed by several rounds of measurement
       of the qubit 0 and collapse of superposed waveforms
       (reducing the states) then toggling the state bits.
      
     NOTE: the libquantum runtime complexity is approx N^2 * log_2(N) * log_2(N),
     but here, have reduced the number of qubits at initialization to
     2^(log2(N)) instead of 2^(log2(N*N)),
        so the runtime complexity is now 
        approx N * log_2(N) * log_2(N)
     
     @return returns 2 factors of number, else returns a single item error code. 
     */
    public int[] run() {
       
        //TODO: consider adding ability to cache cofactors, that is x,
        //   to not retry same x on subsequent run.
        
        
        // max width = 30 ==> max N is 32768, constrained by array length
        int width = MiscMath.numberOfBits(N * N);
        int swidth = MiscMath.numberOfBits(N);
        if (!useLargerInit) {
            width = swidth;
        }
        
        log.info("SEED=" + rSeed);
        log.info(String.format("N = %d, width=%d, swidth=%d, %d qubits required\n", 
            N, width, swidth, width+3*swidth+2));
        
        if (x == 0) {
            Classic classic = new Classic();
            while ((classic.quantum_gcd(N, x) > 1) || (x < 2)) {
                //x = rand() % N;
                x = rng.nextInt((1<<16)-1) % N;
            }
        }
        
        int i;
        int q,a,b, factor;

        log.info(String.format("Random factor: %d of %d\n", x, N));

        QuReg qureg = new QuReg();
        
        QuantumReg qr = qureg.quantum_new_qureg(0, width);
        
        assert(qr.hash.length == (1 << qr.hashw));
        
        Gates gates = new Gates(rng);
        
        // initialize 1<<width bitstrings and normalize them
        //   the bitstring is the bitstring representation of a node.state in
        //   register qr, and each bit is the figurative qubit in qr.
        //   each bitstring is stored as a node in register qr.
        //   the sum of the node amplitudes squared is approx 1.
        
        // ~O(qr.size) where qr.size is 2^(log2(N))
        for (i = 0; i < width; i++) {
            gates.quantum_hadamard(i, qr);
        }
        
        /*
        //log.info(
        System.out.println(
            "after init and norm with quantum_hadamard: "
            + "reg.size=" + qr.size
            + " hash.length=" + qr.hash.length);
        qureg.quantum_print_qureg(qr);
        */
        assert(qr.hash.length == (1 << qr.hashw));
        
        
        // ---- shift existing node states left by nbits ----
        
        int nbits = 3 * swidth + 2;
        qureg.quantum_addscratch(nbits, qr);
        
        /*
        //log.info(
        System.out.println(
            "after addscratch: "
            + "reg.size=" + qr.size
            + " hash.length=" + qr.hash.length);
        qureg.quantum_print_qureg(qr);
        */
        
        
        // ---- apply exp_mod_n ----
        
        //runtime complexity is width * swidth * O(reg.size).
        //           ~ log_2(N) * log_2(N) * 2^(log_2(N))
        //           ~ N * log_2(N) * log_2(N)
        gates.quantum_exp_mod_n(N, x, width, swidth,  qr);
        
        /*
        //log.info(
        System.out.println(
            "after exp_mod_n: "
            + "reg.size=" + qr.size
            + " hash.length=" + qr.hash.length);
        qureg.quantum_print_qureg(qr);
        */
        assert(qr.hash.length == (1 << qr.hashw));
     
        Measure measure = new Measure();
        
        for (i = 0; i < nbits; i++) {
            measure.quantum_bmeasure(0, qr, rng);
        }
       
        assert(qr.hash.length == (1 << qr.hashw));
 
        /*
        //log.info(
        System.out.println(
            "after bmeasure: "
            + "reg.size=" + qr.size
            + " hash.length=" + qr.hash.length);
        qureg.quantum_print_qureg(qr);
        */
        
        //log2(N) * (.lt. log2(N)) * 2^(log2(N))
        gates.quantum_qft(width,  qr);

        /*
        //log.info(
        System.out.println(
            "after qft: "
            + "reg.size=" + qr.size
            + " hash.length=" + qr.hash.length);
        qureg.quantum_print_qureg(qr);
        //for (i = 0; i < qr.size; i++) {
        //    System.out.format("IV %d %d\n", i, qr.node[i].state);
        //}
        */
        
        //SWAP = CNOT[i, j]CNOT[j, i]CNOT[i, j]
        for (i = 0; i < width / 2; i++) {
            gates.quantum_cnot(i, width - i - 1, qr);
            gates.quantum_cnot(width - i - 1, i, qr);
            gates.quantum_cnot(i, width - i - 1, qr);
        }

        /*
        //log.info(
        System.out.println(
            "after last swap: "
            + "reg.size=" + qr.size
            + " hash.length=" + qr.hash.length);
        qureg.quantum_print_qureg(qr);
        */
        assert(qr.hash.length == (1 << qr.hashw));

        long c = measure.quantum_measure(qr, rng);

        log.info("c=" + c);

        q = 1 << (width);
        
        if (c == -1) {
            log.info(String.format("Impossible Measurement!\n"));
            return new int[]{-1};
        }

        if (c == 0) {
            log.info(String.format("Measured zero, try again.\n"));
            return new int[]{0};
        }

        log.info(String.format("Measured %d (%f), ", c, (float) c / q));

        Classic classic = new Classic();
        int[] cInOut = new int[]{(int)c};
        int[] qInOut = new int[]{q};
        classic.quantum_frac_approx(cInOut, qInOut, width);
        c = cInOut[0];
        q = qInOut[0];
        
        log.info(String.format("fractional approximation is %d/%d.\n", c, q));

        if ((q % 2 == 1) && (2 * q < (1 << width))) {
            log.info(String.format("Odd denominator, trying to expand by 2.\n"));
            q *= 2;
        }

        if (q % 2 == 1) {
            log.info(String.format("Odd period, try again.\n"));
            return new int[]{-2};
        }

        log.info(String.format("%d\n", 11));

        log.info(String.format("Possible period is %d.\n", q));

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
            log.info(String.format("%d = %d * %d\n", N, factor, N / factor));
        } else {
            log.info(String.format("Unable to determine factors, try again.\n"));
            return new int[]{-2};
        }

        return new int[]{factor, N/factor};
    }
}
