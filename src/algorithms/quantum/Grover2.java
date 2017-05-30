package algorithms.quantum;

import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import java.util.Random;
import thirdparty.libquantum.*;

/* 
ported from libquantum's grover.c to java then
changed to reduce the runtime complexity.
NOT YET FINISHED

original file was grover.c: Implementation of Grover's search algorithm

  Copyright 2003 Bjoern Butscher, Hendrik Weimer

  This file is a port to java from a c file in libquantum.

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

public class Grover2 {

    /**
     * NOT READY FOR USE
     * 
     * @param number
     * @return 
     */
    public int run(int number) {

        int width = MiscMath.numberOfBits(number + 1);

        return run(number, width);
    }

    /**
     * NOT READY FOR USE
     * 
     * @param number
     * @param width
     * @return 
     */
    public int run(int number, int width) {

        /*
        Quantum Mechanics Helps in Searching for a Needle in a Haystack
          by Grover
        
        - N number of items in a list
        - n = log2(N) is number
          (below, width = n.  length of state bitstring = n)
        
        -- each state S_i is represented by an n-bit bitstring.
        -- there are 1<<n states which is ~N
        
        Sν is a unique state, where condition C(Sν) = 1, 
           whereas for all other states S, C(S) = 0 
        
        (assume that condition C(S) can be evaluated in O(1))
        
        (1) init system to superposition such that each of the
            N states has amplitude (1/sqrt(N)).
        
            QuantumReg reg = qureg.quantum_new_qureg(0, width);
        
            NOTE: algorithm is already O(N) if create all possible states,
            gates.quantum_hadamard(i, reg); for i = 0 to N-1
               but he may have wanted to perform discrete operations
        
        (2) Repeat the following unitary operations
            times (the precise number of repetitions is important
            as discussed in [5]):
        
            NOTE that each loop of this step increases the amplitude
            in the desired state by O(1/sqrt(N)) ???

            (a) Let the system be in any state S:
                In case C(S) = 1, rotate the phase by pi radians;
                In case C(S) = 0, leave the system unaltered.
        
                This is a phase rotation transformation.
                see sect 1.1
                
             selective rotation of phase of the amplitude in certain states.
             transformation describing this for a 2 state system is of the
             form:
                   | e^(j*phi_1)   0           | 
                   | 0             e^(j*phi_2) |
        
                   where j = sqrt(-1) and phi_1, phi_2 are real numbers
                   
        
            (b) Apply the diffusion transform D which is
                defined by the matrix D as follows:
        
                D_i_j = (2/N) if i != j 
                and D_i_i = -1 + (2/N)
                
                (D can be implemented as a product of 3 elementary
                matrices as discussed in section 5).
        
                This diffusion transform, D, can be
                implemented as D R W, where R the
                rotation matrix and W the Walsh-Hadamard
                Transform Matrix are defined as follows:
        
                R_i_j = 0 if i!=j
                R_i_i = 1 if i==0
                R_i_i = -1 if i!=0
              
                W_i_j = (1<<(-n/2)) * (-1)^(i dot j)
                     where i is the binary presentation of i
                           i dot j is the is the bitwise dot product of
                              the bitstring i and bitstring j (both of size n)
                              (bitwise dot product is '&')
        
        (3) Measure the resulting state. 
            This will be the state Sν
            (i.e. the desired state that satisfies the condition) 
            with a probability of at least 0.5.
        */
        
        int i;

        final int N = number;

        Random rng = Misc.getSecureRandom();

        Gates gates = new Gates(rng);

        int tmp = MiscMath.numberOfBits(N + 1);
        if (width < tmp) {
            width = tmp;
        }
        
        /*
        from wikipedia:
        
                   
                       ------------------------
                      /   diffuser              \
              _____   _____  ____________   _____
|0> -[H⊗n]---|     |--|H⊗n|--|2|0^n> -I_n|--|H⊗n|---- ...measure
             | U_w |  -----  ------------   ----
|1> -[H]-----|     |---------------------------------
             ------|
    
             Repeat U_w + diffuser O(sqrt(N)) times
        
        also see tutorial at:
        http://people.cs.umass.edu/~strubell/doc/quantum_tutorial.pdf
        */

        System.out.format("N = %d, width=%d\n", N, width);

        double rotation = Math.PI;

        QuReg qureg = new QuReg();
        
        QuantumReg reg = qureg.quantum_new_qureg(0, width);

        //DEBUG
        System.out.format("AFTER construction  reg.size=%d\n", reg.size);
        qureg.quantum_print_qureg(reg);
        
        //runtime complexity is O(reg.size * reg.width)
        for (i = 0; i < reg.width; i++) {
            gates.quantum_hadamard(i, reg);
        }
                
        //DEBUG
        System.out.format("AFTER 1st hadamard gates  reg.size=%d\n", reg.size);
        qureg.quantum_print_qureg(reg);
        
        int end = (int) ((Math.PI / 4) * Math.sqrt(1 << reg.width));
        
        System.out.format("Iterating %d times\n", end);
        
        i = rng.nextInt(reg.size);
        
        long bestMBits = 0;
        
        int nIter = 0;
        while (nIter < end) {
                        
            System.out.println("i=" + i);
            
            // U_w oracle
            //(2a) phase rotation by pi if in correct state
            
            // the oracle in this case, is the matching bits
            long mbits = reg.node[i].state & N;
            if (mbits > 0) {
                System.out.println("matched bits =" + Long.toBinaryString(mbits));
                //cos(phi) + IMAGINARY * sin(phi)
                double re = reg.node[i].amplitude.re() * Math.cos(rotation);
                //double im = reg.node[i].amplitude.im() * Math.sin(rotation);
                reg.node[i].amplitude.setReal(re);
            }
            
            //|H⊗n|  |2|0^n> -I_n|  |H⊗n|
                        
            //(2b) diffusion transform
            // determine the average amplitude (or keep a running calculation)
            // and then alter this amplitude by the negative of the difference.
            
            // inversion about the average, 
            // transforming the amplitude of each state so that
            // it is as far above the average as it was below the average prior 
            // to the transformation, and vice versa. 
            // This diffusion transform consists of another application of 
            // the Hadamard transform H⊗n, 
            // followed by a conditional phase shift that shifts every state 
            // except |0> by −1, followed by yet another Hadamard transform. 
            // The conditional phase shift can be
            
            double avg = 0;
            for (int j = 0; j < reg.size; ++j) {
                avg += reg.node[j].amplitude.re();
            }
            avg /= (double)reg.size;
            System.out.println("avg=" + avg);
            
            // need to change amplitudes by their difference from avg,
            // but also need the sum of squared amplitudes to be 1
            for (int j = 0; j < reg.size; ++j) {
                double a = reg.node[j].amplitude.re();
                // a = 2*avg - a
                reg.node[j].amplitude.setReal(2. * avg - a);
            }
            
            // renorm
            double sumsq = 0;
            for (int j = 0; j < reg.size; ++j) {
                sumsq += reg.node[j].amplitude.squareSum();
            }
            double div = Math.sqrt(sumsq);
            for (int j = 0; j < reg.size; ++j) {
                reg.node[j].amplitude.setReal(
                    reg.node[j].amplitude.re() / div);
                reg.node[j].amplitude.setImag(
                    reg.node[j].amplitude.im() / div);
            }
            
            //DEBUG
            System.out.format("AFTER grover nITer=%d\n", nIter);
            qureg.quantum_print_qureg(reg);
        
            //TODO: combine mbits w/ bestMBits to create new
            //   bestMBits
            
            // TODO: for choice of next index, base it on best mbits so far
            // or if that is 0, flip all bits
            i = rng.nextInt(reg.size);
            
            /*
            (b) Apply the diffusion transform D which is
                defined by the matrix D as follows:
        
                D_i_j = (2/N) if i != j 
                and D_i_i = -1 + (2/N)
                
                (D can be implemented as a product of 3 elementary
                matrices as discussed in section 5).
        
                This diffusion transform, D, can be
                implemented as D R W, where R the
                rotation matrix and W the Walsh-Hadamard
                Transform Matrix are defined as follows:
        
                R_i_j = 0 if i!=j
                R_i_i = 1 if i==0
                R_i_i = -1 if i!=0
              
                W_i_j = (1<<(-n/2)) * (-1)^(i dot j)
                    where i is the binary presentation of i
                        i dot j is the is the bitwise dot product of
                            the bitstring i and bitstring j (both of size n)
                            (bitwise dot product is '&')
            */
            
            
            nIter++;
        }
        
        //NOTE: the above isn't correct yet.
        //     just reading the first grover paper to impl what it
        //     describes, and the 2nd and another that it references...

        return 0;
    }
}
