package thirdparty.libquantum;

import algorithms.misc.ComplexModifiable;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import java.util.Random;

/* 
An implementation of the grover search algorithm, 
ported here to java from the libquantum file grover.c.
The method calls have been adapted for re-use by
other algorithms and methods to accept a list of
numbers have been created and an alternate oracle
and diffuser have been implemented by this project.

The file grover.c has copyright:
Implementation of Grover's search algorithm

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

public class Grover {
    
    private boolean debug = false;
    
    private int width0 = 0;

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
     */
    
    /**
     * 
     * runtime complexity is O(reg.size * reg.width),
       (because decoherence lambda is 0.0).
     * 
     * @param state (f(x) == 1 when x == state, else f(x) == 0)
     * @param reg
     */
    private void oracle(int query, QuantumReg reg, Gates gates) {
        int i;
        
        //NOTE: none of the quantum api methods has to check for unique
        //   quantum numbers before making a change, that is they do not 
        //   necessarily check pauli exclusion principle, or check for Cooper pairs
        
        /*
         function f(x)
                == 1 when x satisifies search criteria, 
                   that is, x == w
                   |U_w|x> = -|x>
                == 0 else is 0, that is, x != w
                   |U_w|x> = |x>

        // |x>|q> ----> (-1)^(f(x)) * |x>        
        */
        
        /*
        -- for each query bit: 
              if query bit i is 0, flips that bit in all states
        -- for each node state,
              if bits 0 and 1 are set, 
                  it flips the bit reg->width + 1
        -- for each node state,
              if bit reg->width + i is set, 
                  it flips the bit reg->width
        -- for each node.state (in reversed order): 
               if bits i and reg->width + i are set, 
                   it flips the bit reg.width + 1 + i
        -- for each node state,
               if bits 0 and 1 are set, 
                  it flips the bit reg->width + 1
        -- for each query bit:
               if query bit i is 0, 
                   flip bit i in all node states
        */
             
        //DEBUG
        StringBuilder[] sbs = new StringBuilder[reg.size];
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = new StringBuilder();
                sbs[ii] = sb;
            }
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
        }
        
        //runtime complexity is O(reg.size * reg.width),
        // (because decoherence lambda is 0.0).
        for (i = 0; i < width0; i++) {
            //if query bit i is 0, flip bit i in all node states
            if ((query & (1 << i)) == 0) {
                gates.quantum_sigma_x(i, reg);
            }
        }
        
        
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
        }

        //for each node.state: 
        // if bits 0 and 1 are set, it flips the bit reg->width + 1
        gates.quantum_toffoli(0, 1, width0 + 1, reg);
        
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
        }
        
        for (i = 1; i < width0; i++) {
            //for each node.state: 
            // if bits i and reg->width + i are set, 
            // it flips the bit reg->width + 1 + i
           
            gates.quantum_toffoli(i, width0 + i, width0 + i + 1, reg);
        }
        
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
        }
        
        //for each node.state: 
        // if bit reg->width + i is set, 
        // it flips the bit reg->width
        gates.quantum_cnot(width0 + i, width0, reg);
        
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
        }
        
        for (i = width0 - 1; i > 0; i--) {
            //for each node.state: 
            // if bits i and reg->width + i are set, 
            // it flips the bit reg.width + 1 + i
            
            gates.quantum_toffoli(i, width0 + i, width0 + i + 1, reg);
        }
        
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
        }
        
        //for each node.state: 
        // if bits 0 and 1 are set, 
        // it flips the bit reg->width + 1
        gates.quantum_toffoli(0, 1, width0 + 1, reg);

        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
        }
        
        for (i = 0; i < width0; i++) {
            //if query bit i is 0, flip bit i in all node states
            if ((query & (1 << i)) == 0) {
                gates.quantum_sigma_x(i, reg);
            }
        }
        
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
            
            System.out.println("END STATES");
            for (int ii = 0; ii < reg.size; ii++) {
                System.out.println(sbs[ii]);
            }
        }
    }   
    
    /*
    private void oracle3(int query, QuantumReg reg, Gates gates) {
        int i;
        
        //NOTE: since none of the quantum api methods has to check for unique
        //   quantum numbers before making a change, that is they do not 
        //   necessarily check pauli exclusion principle, or check for Cooper pairs
        //   it should be just as valid to
        //   use test and set bit operations alone rather than
        //   purely test and toggle bit operations as nicely demonstrated
        //   in oracle(...)
        
        //function f(x)
        //       == 1 when x satisifies search criteria, 
        //         that is, x == w
        //         |U_w|x> = -|x>
        //      == 0 else is 0, that is, x != w
        //         |U_w|x> = |x>
        //
        // |x>|q> ----> (-1)^(f(x)) * |x>        
             
        //DEBUG
        StringBuilder[] sbs = new StringBuilder[reg.size];
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = new StringBuilder();
                sbs[ii] = sb;
            }
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
        }
        
        // if 0 bit is same in query and state, set high bit
        for (i = 0; i < reg.size; ++i) {
            int qs = query & (1 << 0);
            long vs = reg.node[i].state & (1 << 0);
            if (((qs == 0) && (vs == 0)) || ((qs != 0) && (vs != 0))
                ) {
                reg.node[i].state  |= (1L << (reg.width - 1));
            }
        }
        
        // if high bit is set and query bit is same in query and state, 
        // keep high bit else unset.
        // NOTE: all of this conditional logic might not be available in
        //  a gate
        for (i = 0; i < reg.size; ++i) {
            if ((reg.node[i].state  & (1L << (reg.width - 1))) != 0) {
                int qs = query & (1 << 0);
                long vs = reg.node[i].state & (1 << 0);
                if (!(((qs == 0) && (vs == 0)) || ((qs != 0) && (vs != 0))
                    )) {
                    reg.node[i].state  &= ~(1L << (reg.width - 1));
                }
            }
        }
        
        // rotate by pi if high bit is set
        gates.quantum_phase_kick(reg.width - 1, Math.PI, reg);
        
        // unset the highest bit.
        // not reversible.
        for (i = 0; i < reg.size; i++) {
            reg.node[i].state  &= ~(1L << (reg.width - 1));
        }
        
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
            
            System.out.println("END STATES");
            for (int ii = 0; ii < reg.size; ii++) {
                System.out.println(sbs[ii]);
            }
        }
    }
    */
    
    private void oracle2(int query, QuantumReg reg, Gates gates) {
        int i;
        
        //function f(x)
        //       == 1 when x satisifies search criteria, 
        //         that is, x == w
        //         |U_w|x> = -|x>
        //      == 0 else is 0, that is, x != w
        //         |U_w|x> = |x>
        //
        // |x>|q> ----> (-1)^(f(x)) * |x>        
             
        //DEBUG
        StringBuilder[] sbs = new StringBuilder[reg.size];
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = new StringBuilder();
                sbs[ii] = sb;
            }
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
        }
        
        //runtime complexity is O(reg.size * reg.width),
        // (because decoherence lambda is 0.0).
        for (i = 0; i < width0; i++) {
            //if query bit i is 0, flip bit i in all node states
            if ((query & (1 << i)) == 0) {
                gates.quantum_sigma_x(i, reg);
            }
        }
        
        
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
        }

        //for each node.state: 
        // if bits 0 and 1 are set, it flips the bit reg->width + 1
        gates.quantum_toffoli(0, 1, width0 + 1, reg);
        
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
        }
        
        for (i = 1; i < width0; i++) {
            //for each node.state: 
            // if bits i and reg->width + i are set, 
            // it flips the bit reg->width + 1 + i
           
            gates.quantum_toffoli(i, width0 + i, width0 + i + 1, reg);
        }
        
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
        }
        
        //for each node.state: 
        // if bit reg->width + i is set, 
        // it flips the bit reg->width
        gates.quantum_cnot(width0 + i, width0, reg);
        
             
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append("*").append(str).append("  ");
            }
        }
        
        
        // anything that is all 1's at this point is the found number.
        // (all 1's except the high bit)
        
        // not reversible
        // set high bit if all bits are set
        gates.quantum_ccand(0, 0, reg.width-1, reg);
        for (i = 1; i < reg.width-2; i++) {
            gates.quantum_ccand(i, reg.width-1, reg.width-1, reg);
        }
        
        
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append("*").append(str).append("  ");
            }
        }
        
        for (i = width0 - 1; i > 0; i--) {
            //for each node.state: 
            // if bits i and reg->width + i are set, 
            // it flips the bit reg.width + 1 + i
            
            gates.quantum_toffoli(i, width0 + i, width0 + i + 1, reg);
        }
        
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
        }
        
        //for each node.state: 
        // if bits 0 and 1 are set, 
        // it flips the bit reg->width + 1
        gates.quantum_toffoli(0, 1, width0 + 1, reg);

        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
        }
        
        for (i = 0; i < width0; i++) {
            //if query bit i is 0, flip bit i in all node states
            if ((query & (1 << i)) == 0) {
                gates.quantum_sigma_x(i, reg);
            }
        }
        
        
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
            
        }
        
        //NOTE: since none of the quantum api methods has to check for unique
        //   quantum numbers before making a change, that is they do not 
        //   necessarily check pauli exclusion principle, or check for Cooper pairs
        //   it should be just as valid to
        //   use test and set pr clear bit operations rather than
        //   purely test and toggle bit operations as nicely demonstrated
        //   in oracle(...)
        
        // not reversible
        // unset the width0 bit
        for (i = 0; i < reg.size; i++) {
            reg.node[i].state  &= ~(1L << width0);
        }
       
        // rotate by pi if high bit is set
        gates.quantum_phase_kick(reg.width - 1, Math.PI, reg);
        
        // unset the highest bit.
        // not reversible.
        for (i = 0; i < reg.size; i++) {
            reg.node[i].state  &= ~(1L << (reg.width - 1));
        }
        
        
        if (debug) {//DEBUG
            for (int ii = 0; ii < reg.size; ii++) {
                StringBuilder sb = sbs[ii];
                String str = Long.toBinaryString(reg.node[ii].state);
                while (str.length() < reg.width) {
                    str = "0" + str;
                }
                sb.append(str).append("  ");
            }
            
            System.out.println("END STATES");
            for (int ii = 0; ii < reg.size; ii++) {
                System.out.println(sbs[ii]);
            }
        }
    }
    
    /**
      runtime complexity is O(reg.size * reg.width),
      (because decoherence lambda is 0.0).
     * @param reg
     */
    private void inversion(QuantumReg reg, Gates gates) {
        int i;

        //|2|0^n> -I_n|
        
        //Flip the target bit of each basis state, i
        for (i = 0; i < width0; i++) {
            gates.quantum_sigma_x(i, reg);
        }
       
        gates.quantum_hadamard(width0 - 1, reg);

        if (width0 == 3) {
        
            gates.quantum_toffoli(0, 1, 2, reg);
        
        } else {
            
            //If bits 0 and 1 are set, it flips the target bit.
            gates.quantum_toffoli(0, 1, width0 + 1, reg);

            for (i = 1; i < width0 - 1; i++) {
                //If bits i and reg.width+i are set, it flips the target bit.
                gates.quantum_toffoli(i, width0 + i, width0 + i + 1, reg);
            }

            //for each reg.state, 
            //   Flip the target bit of a basis state if 
            //   the control bit is set
            gates.quantum_cnot(width0 + i, width0 - 1, reg);

            for (i = width0 - 2; i > 0; i--) {
                //If bits i and reg.width+i are set, it flips the target bit.
                gates.quantum_toffoli(i, width0 + i, width0 + i + 1, reg);
            }

            //If bits 0 and 1 are set, it flips the target bit.
            gates.quantum_toffoli(0, 1, width0 + 1, reg);
        }

        gates.quantum_hadamard(width0 - 1, reg);
        
        //Flip the target bit of each basis state, i
        for (i = 0; i < width0; i++) {
            gates.quantum_sigma_x(i, reg);
        }
        
    }

    /**
     * runtime complexity is O(reg.size * reg.width)  (because decoherence lambda is 0.0).
     * 
     * @param target 
     *     (f(x) == 1 when x == target, else f(x) == 0)
     * @param reg
     */
    private void grover(int target, QuantumReg reg, Gates gates, QuReg qureg) {
        
        int i;  
        
        //unitary operator operating on two qubits, target and each i
        // |x>|q> ----> (-1)^(f(x)) * |x> 
        // (gives the found solutions negative signs)
        oracle(target, reg, gates);

        if (debug) {//DEBUG
            System.out.format(
                "AFTER oracle target=%d  reg.size=%d  hash.length=%d\n", 
                target, reg.size, 1 << reg.hashw);
            qureg.quantum_print_qureg(reg);
        }

        //   H⊗n   |2|0^n> -I_n|  H⊗n

        for (i = 0; i < width0; i++) {
            gates.quantum_hadamard(i, reg);
        }

        if (debug) {//DEBUG
            System.out.format(
                "AFTER hadamard target=%d hadamard reg.size=%d\n", 
                target, reg.size);
            qureg.quantum_print_qureg(reg);
        }

        inversion(reg, gates);

        if (debug) {//DEBUG
            System.out.format(
                "AFTER target=%d inversion reg.size=%d\n", 
                target, reg.size);
            qureg.quantum_print_qureg(reg);
        }

        for (i = 0; i < width0; i++) {
            gates.quantum_hadamard(i, reg);
        }

        if (debug) {//DEBUG
            System.out.format("AFTER target=%d 2nd hadamard  reg.size=%d\n", 
                target, reg.size);
            qureg.quantum_print_qureg(reg);
        }
       
    }
    
    /**
     * runtime complexity is O(reg.size * reg.width)  (because decoherence lambda is 0.0).
     * 
     * @param target 
     *     (f(x) == 1 when x == target, else f(x) == 0)
     * @param reg states with the high bit set are the found items
     */
    private void grover2(int target, QuantumReg reg, Gates gates, QuReg qureg) {
             
        int i;
        
        // matching state is rotated by PI
        oracle2(target, reg, gates);

        if (debug) {//DEBUG
            System.out.format(
                "AFTER oracle2 target=%d  reg.size=%d  hash.length=%d\n", 
                target, reg.size, 1 << reg.hashw);
            qureg.quantum_print_qureg(reg);
        }

        // amplify the state with high bit set       
                
        /*
        double avg = 0;
        for (int j = 0; j < reg.size; ++j) {
            avg += reg.node[j].amplitude.re();
        }
        avg /= (double) reg.size;
        */
        
        // NOTE: not knowing the avg from quantum methods, assume it
        //       to be the normalization
        // NOTE: also, assuming that an amplifier is available and a
        //       beam splitter or the equivalent
        double avg = 1./Math.sqrt(reg.size);

        // change amplitudes by their difference from avg
        for (int j = 0; j < reg.size; ++j) {
            double a = reg.node[j].amplitude.re();
            // a = 2*avg - a
            reg.node[j].amplitude.setReal(2. * avg - a);
            reg.node[j].amplitude.times(1./Math.sqrt(2.));
        }
        
        if (debug) {//DEBUG
            double sum = 0;
            for (i = 0; i < reg.size; ++i) {
                sum += reg.node[i].amplitude.squareSum();
            }
            System.out.println("sum of squares=" + sum);
        }
        
        if (debug) {//DEBUG
            System.out.format(
                "AFTER grover target=%d reg.size=%d reg.width=%d\n", 
                target, reg.size, reg.width);
            qureg.quantum_print_qureg(reg);
        }
    }

    /** runtime complexity is O(reg.size * reg.width) * nLoop
       (the runtime complexity of the preparation of the register 
     * is ignored.  it is O(2^width)).
     * Note that nLoop is (Math.PI / 4) * Math.sqrt(2^width)
     * where width is (the bit length of number) + 1
     * 
     * @param number a number to search for in the enumeration of numbers
     * from 0 to 2^(number bit length + 1)
    */
    public int run(int number) {

        int width = MiscMath.numberOfBits(number + 1);

        return run(number, width);
    }

    /**
     * runtime complexity is O(reg.size * reg.width) * nLoop.
     * Note that nLoop is (Math.PI / 4) * Math.sqrt(2^width).
     * (the runtime complexity of the preparation of the register 
     * is ignored.  it is O(2^width)).
     * 
     * @param number a number to search for in the enumeration of numbers
     * from 0 to 2^width.
     * @param width largest bit length to use in enumeration.
     * NOTE that if it is less than (the bit length of number) + 1,
     * it will be increased to that.
     * @return 
     */
    public int run(int number, int width) {

        int i;

        final int N = number;

        Random rng = Misc.getSecureRandom();

        Gates gates = new Gates(rng);

        int tmp = MiscMath.numberOfBits(N + 1);
        if (width < tmp) {
            width = tmp;
        }
        if (width < 2) {
            width = 2;
        }
        width0 = width;

        QuReg qureg = new QuReg();

        QuantumReg reg = qureg.quantum_new_qureg(0, width0);

        if (debug) {//DEBUG
            System.out.format(
                "AFTER construction  reg.size=%d reg.width=%d hash.length=%d\n", 
                reg.size, reg.width, 1 << reg.hashw);
            qureg.quantum_print_qureg(reg);
        }
        
        //Flip the target bit of each basis state, reg.width
        //runtime complexity is O(reg.size) (because decoherence lambda is 0.0).
        qureg.quantum_addscratch(1, reg);
        reg.width *= 2;
        reg.width += 2;
        qureg.quantum_expand_and_reconstruct_hash(reg);
        gates.quantum_sigma_x(width0, reg);
        
        System.out.format("N = %d, width0=%d reg.width=%d\n", N, width0,
            reg.width);

        if (debug) {
            //DEBUG
            System.out.format("AFTER sigma_x  reg.size=%d reg.width=%d\n", 
                reg.size, reg.width);
            qureg.quantum_print_qureg(reg);
        }

        //runtime complexity is O(reg.size * reg.width)
        for (i = 0; i < width0 + 1; i++) {
            gates.quantum_hadamard(i, reg);
        }
        
        if (debug) {//DEBUG
            System.out.format(
                "AFTER 1st hadamard gates  reg.size=%d reg.width=%d hash.length=%d\n", 
                reg.size, reg.width, 1 << reg.hashw);
            qureg.quantum_print_qureg(reg);
        }

        if (debug) {//DEBUG
            System.out.format("AFTER 2 1st hadamard gates  reg.size=%d reg.width=%d hash.length=%d\n", 
                reg.size, reg.width, 1 << reg.hashw);
            qureg.quantum_print_qureg(reg);
        }

        return search(number, rng, gates, qureg, reg);
    }
    
    private int search(int number, Random rng, Gates gates, QuReg qureg, 
        QuantumReg reg) {
        
        int i;
        
        // upper limit to number of iterations from:
        //"Tight Bounds on Quantum Searching" by Boyer, Brassard, Hoyer, and Tapp 
        int end = (int) (Math.PI / 4 * Math.sqrt(1 << width0));

        System.out.format("Iterating %d times\n", end);

        //runtime complexity is O(reg.size * reg.width) * nLoop
        for (i = 1; i <= end; i++) {
            
            if (debug) {
                System.out.format("Iteration #%d\n", i);
            }
            
            grover(number, reg, gates, qureg);
        }

        if (debug) { //DEBUG
            System.out.format(
                "AFTER grover  reg.size=%d reg.width=%d\n", 
                reg.size, reg.width);
            qureg.quantum_print_qureg(reg);
        }
        

        gates.quantum_hadamard(width0, reg);


        if (debug) {//DEBUG//DEBUG
            System.out.format(
                "AFTER last hadamard  reg.size=%d reg.width=%d\n", 
                reg.size, reg.width);
            qureg.quantum_print_qureg(reg);
        }
        
        Measure measure = new Measure();

        // runtime complexity is O(reg.size)
        measure.quantum_bmeasure(width0, reg, rng);

        //DEBUG
        System.out.format(
            "AFTER bmeasure reg.size=%d reg.width=%d\n", 
            reg.size, reg.width);
        qureg.quantum_print_qureg(reg);


        for (i = 0; i < reg.size; i++) {
            if (reg.node[i].state == number) {
                System.out.format(
                    "\nFound %d with a probability of %f\n\n", number,
                    reg.node[i].amplitude.squareSum());
                return number;
            }
        }

        return -1;
    }
    
    // ---- adding ability to find number within a list of numbers for use 
    //      within the quantum min algorithm ----
    
    /**
     * runtime complexity for the search 
     * is O(reg.size * reg.width) * nLoop
     * (the runtime complexity of the preparation of the register for the list, 
     * O(N), 
     * is ignored just as in the enumerated run method).
     * NOTE that the width should be set to the most number of bits needed
     * for any number in list. 
     * NOTE also that the largest number in the list must be
     * .lte. integer.max_value - 2^width.
     * @param number a number to search for in the enumeration of numbers
     * from 0 to 2^width.
     * @param width largest bit length to use in enumeration.
     * NOTE that if it is less than (the bit length of number) + 1,
     * it will be increased to that.
     * @param list a list of unordered numbers to search for number within.
     * NOTE that the list must be valid input.
     * @return 
     */
    public int run(int number, int width, int[] list) {

        int N = number;
        int i;
        
        int tmp = MiscMath.numberOfBits(N);
        if (width < tmp) {
            width = tmp;
        }
        if (width < 2) {
            width = 2;
        }
        width0 = width;
        
        QuReg qureg = new QuReg();

        Random rng = Misc.getSecureRandom();
        
        int ret;
        
        if (false) {
            
            QuantumReg reg = initializeRegister(qureg, list);

            System.out.format("N = %d, list.length=%d, width0=%d reg.width=%d\n", N, 
                list.length, width0, reg.width);

            ret = processInitialized(number, reg, rng);

        } else {
            
            QuantumReg reg = initializeRegister2(qureg, list);

            System.out.format("N = %d, list.length=%d, width0=%d reg.width=%d\n", N, 
                list.length, width0, reg.width);

            ret = processInitialized2(number, reg, rng);
        }
        
        return ret;
    }
    
    public int run(int number, int width, int setQuBits) {
    
        width0 = width;
        
        QuReg qureg = new QuReg();
        QuantumReg reg = initializeRegister(qureg, setQuBits);
        
        System.out.format("N = %d, width0=%d reg.width=%d\n", number, 
            width0, reg.width);
        
        Random rng = Misc.getSecureRandom();
        
        int ret = processInitialized2(number, reg, rng);
        
        return ret;
    }
    
    /**
     * Initialize the register with a list of numbers.
     * Note, that the register, as the possible states of superposition of
     * qubits, will have all possible permutation of the qubits up to the
     * power of 2 or next higher power of 2 in the list.
     * A continuous sequence of numbers from 0 up to a power of 2 is valid 
     * input for the current logic (can be unordered).
     * A continuous sequence of numbers from a power of 2 up to a power of 2.
     * is valid input.
     */
    
    /**
     * Initialize the register with a list of numbers as the eigenstate,
     * superposition, and their amplitudes.
     *
     * @param qureg
     * @param stateList
     * @param amplList amplitudes associated with the eigenstate at same index
     * in stateList
     * @return 
     */
    public QuantumReg initializeRegister(QuReg qureg, 
        ComplexModifiable[] amplList, int[] stateList) {
      
        final int initSize = 2 * amplList.length;

        QuantumReg reg = qureg.quantum_new_qureg_size(initSize, width0);
        reg.width *= 2;
        reg.width += 2;
        qureg.quantum_expand_and_reconstruct_hash(reg);
        
        //need to initialize a register to have the given states from list
        //and a set of the same numbers but with negative amplitude and
        //the next highest bit set, that is width + 1
        //    rest of the algorithm should proceed in same manner.
        
        int offset = 1 << width0;
        
        int i;
        double invSqrt = 1./Math.sqrt(2.);
        
        for (i = 0; i < amplList.length; ++i) {
            reg.node[i].state = stateList[i];
            reg.node[i].state |= offset;
            reg.node[i].amplitude.resetTo(amplList[i]);
            reg.node[i].amplitude.times(-invSqrt);
        }
        
        int idx = amplList.length;
        for (i = 0; i < amplList.length; ++i) {
            reg.node[idx].state = stateList[i];
            reg.node[idx].amplitude.resetTo(amplList[i]);
            reg.node[idx].amplitude.times(invSqrt);
            idx++;
        }
        
        return reg;
    }
    
    /**
     * Initialize the register with a list of numbers as the eigenstate,
     * superposition, and their amplitudes.
     *
     * @param qureg
     * @param setBits
     * @return 
     */
    public QuantumReg initializeRegister(QuReg qureg, int setBits) {
       
        int nBits = MiscMath.numberOfBits(setBits);
        
        int i;
        
        int nSetBits = 0;
        for (i = 0; i < nBits; ++i) {
            if ((setBits & (1 << i)) != 0) {
                nSetBits++;
            }
        }
        
        QuantumReg reg = qureg.quantum_new_qureg_size(
            nSetBits, width0);

        reg.width *= 2;
        reg.width += 2;
        qureg.quantum_expand_and_reconstruct_hash(reg);
        
        double norm = 1./Math.sqrt(nSetBits);  
        int ii = 0;
        for (i = 0; i < nBits; ++i) {
            if ((setBits & (1 << i)) != 0) {
                //initializing with same state + highbit off of register
                reg.node[ii].state = (1 << i);
                reg.node[ii].amplitude.setReal(norm);
                ++ii;
            }
        }
                
        if (debug) {//DEBUG
            System.out.format("initialized  reg.size=%d reg.width=%d\n", 
                reg.size, reg.width);
            qureg.quantum_print_qureg(reg);
        }
        
        return reg;
    }
    
    /**
     * Initialize the register with a list of numbers as the eigenstate,
     * superposition, and their amplitudes.
     *
     * @param qureg
     * @param setBits
     * @param width
     * @return 
     */
    private QuantumReg initializeRegister(QuReg qureg, int[] list) {
        
        int listLen = list.length;
        
        final int initSize = 2 * listLen;

        QuantumReg reg = qureg.quantum_new_qureg_size(initSize, width0);

        reg.width *= 2;
        reg.width += 2;
        qureg.quantum_expand_and_reconstruct_hash(reg);
        
        //need to initialize a register to have the given states from list
        //and a set of the same numbers but with negative amplitude and
        //the next highest bit set, that is width + 1
        //    rest of the algorithm should proceed in same manner.
                
        int offset = 1 << width0;
        
        int i;
        
        double norm = 1./Math.sqrt(initSize);  
        for (i = 0; i < list.length; ++i) {
            reg.node[i].state = list[i] + offset;
            reg.node[i].amplitude.setReal(-norm);
        }
        int ii = list.length;
        for (i = 0; i < list.length; ++i) {
            reg.node[ii].state = list[i];
            reg.node[ii].amplitude.setReal(norm);
            ii++;
        }
        
        if (debug) {//DEBUG
            System.out.format("AFTER init reg.size=%d "
                + " width0=%d reg.width=%d reg.hash.length=%d\n", reg.size,
                width0, reg.width, (1 << reg.hashw));
            qureg.quantum_print_qureg(reg);
        }
       
        return reg;
    }
    
    /**
     * Initialize the register with a list of numbers as the eigenstate,
     * superposition, and their amplitudes.
     *
     * @param qureg
     * @param setBits
     * @param width
     * @return 
     */
    private QuantumReg initializeRegister2(QuReg qureg, int[] list) {
        
        int listLen = list.length;
        
        final int initSize = listLen;

        QuantumReg reg = qureg.quantum_new_qureg_size(initSize, width0);

        reg.width *= 2;
        reg.width += 2;
        qureg.quantum_expand_and_reconstruct_hash(reg);
        
        //need to initialize a register to have the given states from list
        //and a set of the same numbers but with negative amplitude and
        //the next highest bit set, that is width + 1
        //    rest of the algorithm should proceed in same manner.
                
        int i;
        
        double norm = 1./Math.sqrt(initSize);  
        for (i = 0; i < list.length; ++i) {
            reg.node[i].state = list[i];
            reg.node[i].amplitude.setReal(norm);
        }
        
        if (debug) {//DEBUG
            System.out.format("AFTER init reg.size=%d "
                + " width0=%d reg.width=%d reg.hash.length=%d\n", reg.size,
                width0, reg.width, (1 << reg.hashw));
            qureg.quantum_print_qureg(reg);
        }
       
        return reg;
    }
    
    /**
     * runtime complexity for the processing 
     * is O(reg.size * reg.width) * nLoop
     * (the runtime complexity of the preparation of the register for the list, 
     * O(N),
     * is ignored just as in the enumerated run method).
     * NOTE that the width should be set to the most number of bits needed
     * for any number in list. 
     * NOTE also that the largest number in the list must be
     * .lte. integer.max_value - 2^width.
     * NOTE that measurements of register reg are not taken.
     * @param number a number to search for within the initialized register reg
     * @param reg initialized register which holds nodes of state which are 
     * searched and have amplitudes which when squared and summed over register 
     * are equal to 1.
     * @param rng
     * @return
     */
    public int processInitialized(int number, QuantumReg reg, Random rng) {
        
        int i;

        final int N = number;

        QuReg qureg = new QuReg();
        
        //DEBUG
        //System.out.format("AFTER construction  reg.size=%d\n", reg.size);
        //qureg.quantum_print_qureg(reg);

        // upper limit to number of iterations from:
        //"Tight Bounds on Quantum Searching" by Boyer, Brassard, Hoyer, and Tapp 
        //  NOTE that if the number of times number will appear in list
        //     is known ahead of time,
        //     the term in the sqrt can be divided by that multiplicity.
        int end = (int) (Math.PI / 4 * Math.sqrt(1 << width0));

        System.out.format("Iterating %d times\n", end);

        Gates gates = new Gates(rng);

        //runtime complexity is O(reg.size * reg.width) * nLoop
        for (i = 1; i <= end; i++) {
            
            if (debug) {
                System.out.format("Iteration #%d\n", i);
            }
            
            grover(N, reg, gates, qureg);
        }

        if (debug) { //DEBUG
            System.out.format(
                "AFTER grover  reg.size=%d reg.width=%d\n", 
                reg.size, reg.width);
            qureg.quantum_print_qureg(reg);
        }
        

        gates.quantum_hadamard(width0, reg);


        if (debug) {//DEBUG//DEBUG
            System.out.format(
                "AFTER last hadamard  reg.size=%d reg.width=%d\n", 
                reg.size, reg.width);
            qureg.quantum_print_qureg(reg);
        }
        
        Measure measure = new Measure();

        // runtime complexity is O(reg.size)
        measure.quantum_bmeasure(width0, reg, rng);

        //DEBUG
        System.out.format(
            "AFTER bmeasure reg.size=%d reg.width=%d\n", 
            reg.size, reg.width);
        qureg.quantum_print_qureg(reg);


        for (i = 0; i < reg.size; i++) {
            if (reg.node[i].state == N) {
                System.out.format(
                    "\nFound %d with a probability of %f\n\n", N,
                    reg.node[i].amplitude.squareSum());
                return number;
            }
        }

        return -1;
    }

    /**
     * runtime complexity for the processing 
     * is O(reg.size * reg.width) * nLoop
     * (the runtime complexity of the preparation of the register for the list, 
     * O(N),
     * is ignored just as in the enumerated run method).
     * NOTE that the width should be set to the most number of bits needed
     * for any number in list. 
     * NOTE also that the largest number in the list must be
     * .lte. integer.max_value - 2^width.
     * NOTE that measurements of register reg are not taken.
     * @param number a number to search for within the initialized register reg
     * @param reg initialized register which holds nodes of state which are 
     * searched and have amplitudes which when squared and summed over register 
     * are equal to 1.
     * @param rng
     * @return
     */
    public int processInitialized2(int number, QuantumReg reg, Random rng) {
        
        int i;

        final int N = number;

        QuReg qureg = new QuReg();
        
        //DEBUG
        //System.out.format("AFTER construction  reg.size=%d\n", reg.size);
        //qureg.quantum_print_qureg(reg);

        // upper limit to number of iterations from:
        //"Tight Bounds on Quantum Searching" by Boyer, Brassard, Hoyer, and Tapp 
        //  NOTE that if the number of times number will appear in list
        //     is known ahead of time,
        //     the term in the sqrt can be divided by that multiplicity.
        int end = (int) (Math.PI / 4 * Math.sqrt(1 << width0));

        System.out.format("Iterating %d times\n", end);

        Gates gates = new Gates(rng);

        //runtime complexity is O(reg.size * reg.width) * nLoop
        for (i = 1; i <= end; i++) {
            
            if (debug) {
                System.out.format("Iteration #%d\n", i);
            }
            
            grover2(N, reg, gates, qureg);
        }

        if (debug) { //DEBUG
            System.out.format(
                "AFTER grover  reg.size=%d reg.width=%d\n", 
                reg.size, reg.width);
            qureg.quantum_print_qureg(reg);
        }
        
        Measure measure = new Measure();

        // runtime complexity is O(reg.size)
        measure.quantum_bmeasure(width0, reg, rng);

        //DEBUG
        System.out.format(
            "AFTER bmeasure reg.size=%d reg.width=%d\n", 
            reg.size, reg.width);
        qureg.quantum_print_qureg(reg);


        for (i = 0; i < reg.size; i++) {
            if (reg.node[i].state == N) {
                System.out.format(
                    "\nFound %d with a probability of %f\n\n", N,
                    reg.node[i].amplitude.squareSum());
                return number;
            }
        }

        return -1;
    }
}
