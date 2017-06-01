package thirdparty.libquantum;

import algorithms.misc.ComplexModifiable;
import algorithms.misc.MiscMath;
import java.util.Arrays;

/* qureg
Quantum register management

   Copyright 2003, 2004 Bjoern Butscher, Hendrik Weimer

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

public class QuReg {

    /**
     * Our 64-bit multiplicative hash function
     * that is actually 31 bit.
     * TODO: refactor name when refactor for camel case.
     */
    /*static inline unsigned int quantum_hash64(MAX_UNSIGNED key, int width) {
        unsigned int k32;
        k32 = (key & 0xFFFFFFFF) ^ (key >> 32);
        k32 *= 0x9e370001UL;
        k32 = k32 >> (32 - width);
        return k32;
    }*/
    public static int quantum_hash64(long key, int width) {

        //this will always be == key value unless allow key type to be long
        long k = (key & 0x7FFFFFFF) ^ (key >> 31);
        k *= 1327202304L; // factor is 0x9e370001UL >> 1 which is 30.3057
        k &= ((1 << 31) - 1);
        k = k >> (31L - width);
        
        return (int) k;
    }
    
    /**
     * Get the position of a given base state via the hash table
     */
    //static inline int quantum_get_state(MAX_UNSIGNED a, quantum_reg reg) {
    static int quantum_get_state(long a, QuantumReg reg) {
        int i;

        if (reg.hashw == 0) {
            return (int)a;
        }

        i = quantum_hash64(a, reg.hashw);

        while (reg.hash[i] != 0) {
            if (reg.node[reg.hash[i] - 1].state == a) {
                return reg.hash[i] - 1;
            }
            i++;
            if (i == (1 << reg.hashw)) {
                i = 0;
            }
        }

        return -1;

    }

    /**
     * Add an element to the hash table
     */
    //static inline void quantum_add_hash(MAX_UNSIGNED a, int pos, quantum_reg *reg) {
    static void quantum_add_hash(long a, int pos, QuantumReg reg) {

        int i, mark = 0;

        i = quantum_hash64(a, reg.hashw);
        
        int end = 1 << reg.hashw;
        while (reg.hash[i] != 0) {
            i++;
            // if i is > last index
            if (i == end) {
                if (mark == 0) {
                    i = 0;
                    mark = 1;
                } else {
                    StackTraceElement[] st = Thread.currentThread().getStackTrace();
                    for (StackTraceElement s : st) {
                        System.out.println(s);
                    }
                    throw new IllegalStateException("hash is full.  i=" + i + 
                         " hash.lrngth=" + end);
                }
            }
        }
        
        //NOTE: since this value is the index of reg.node,
        //      the last entry is an invalid value.
        //      should the value be pos?
        reg.hash[i] = pos + 1;
    }

    /**
     * Reconstruct hash table
     */
    static void quantum_reconstruct_hash(QuantumReg reg) {
  
        int i;

        /* Check whether register is sorted */
        if (reg.hashw == 0) {
            return;
        }

        int end = 1 << reg.hashw;
        Arrays.fill(reg.hash, 0, end, 0);
   
        assert(reg.node.length >= reg.size);
        
        for (i = 0; i < reg.size; i++) {
            quantum_add_hash(reg.node[i].state, i, reg);
        }
    }
    
    /**
     * Return the reduced bitmask of a basis state
     */
    //static int quantum_bitmask(MAX_UNSIGNED a, int width, int  *bits) {
    static int quantum_bitmask(long a, int width, int[] bits) {
        
        int i;
        int mask = 0;

        for (i = 0; i < width; i++) {
            //if (a & ((MAX_UNSIGNED) 1 << bits[i])) {
            if ((a & (1L << bits[i])) != 0) {
                mask += (1 << i);
            }
        }

        return mask;
    }

    /**
     * Convert a vector to a quantum register
     */
    QuantumReg quantum_matrix2qureg(QuantumMatrix m, int w) {
  
        if (m.cols != 1) {
            throw new IllegalArgumentException("m nCols must == 1");
        }
        
        QuantumReg reg = new QuantumReg();
        int i, j, size = 0;

        reg.width = w;

        // Determine the size of the quantum register 
        for (i = 0; i < m.rows; i++) {
            if (m.t[i].squareSum()!= 0.0) {
                size++;
            }
        }

        // Allocate the required memory 
        reg.size = size;
        reg.hashw = w + 2;

        reg.node = new QuantumRegNode[size];
        for (int ii = 0; ii < reg.node.length; ii++) {
            reg.node[ii] = new QuantumRegNode();
            reg.node[ii].state = 0;
            reg.node[ii].amplitude = new ComplexModifiable(0, 0);
        }
        
        // Allocate the hash table 
        int nHash = 1 << reg.hashw;
        reg.hash = new int[nHash];

        // Copy the nonzero amplitudes of the vector into the 
        //quantum register 
        for (i = 0, j = 0; i < m.rows; i++) {
            if (m.t[i].squareSum() != 0.0) {
                reg.node[j].state = i;
                reg.node[j].amplitude.resetTo(m.t[i]);
                j++;
            }
        }

        return reg;
    }

    /**
     * Create a new quantum register from scratch
     */
    //QuantumReg quantum_new_qureg(MAX_UNSIGNED initval, int width) {
    public QuantumReg quantum_new_qureg(long initval, int w) {
        
        QuantumReg reg = new QuantumReg();

        reg.width = w;
        reg.size = 1;
        reg.hashw = w + 2;
        
        // Allocate memory for 1 base state 
        reg.node = new QuantumRegNode[1];

        // Allocate the hash table 
        int nHash = 1 << reg.hashw;
        reg.hash = new int[nHash];

        // Initialize the quantum register 
        reg.node[0] = new QuantumRegNode();
        reg.node[0].state = initval;
        reg.node[0].amplitude =  new ComplexModifiable(1, 0);

        System.out.format(
            "init reg: %d qubits, 1 node, and %d hash table length, hashw=%d\n",
            reg.width, nHash, reg.hashw);
        
        return reg;
    }

    /**
     * Returns an empty quantum register of size N
     * with instantiated nodes and hash.
     */
    public QuantumReg quantum_new_qureg_size(int n, int w) {
        
        QuantumReg reg = new QuantumReg();

        reg.width = w;
        reg.size = n;

        // Allocate memory for n basis states 
        reg.node = new QuantumRegNode[n];
        for (int ii = 0; ii < reg.node.length; ii++) {
            reg.node[ii] = new QuantumRegNode();
            reg.node[ii].state = 0;
            reg.node[ii].amplitude = new ComplexModifiable(0, 0);
        }
        
        reg.hashw = w + 2;
        
        // Allocate the hash table 
        int nHash = 1 << reg.hashw;
        reg.hash = new int[nHash];
        
        return reg;
    }

    /**
     * Convert a quantum register to a vector
     */
    QuantumMatrix quantum_qureg2matrix(QuantumReg reg) {
        
        int i;

        Matrix matrix = new Matrix();
        int sz = 1 << reg.width;
        QuantumMatrix m = matrix.quantum_new_matrix(1, sz);

        for (i = 0; i < reg.size; i++) {
            m.t[(int)reg.node[i].state].resetTo(reg.node[i].amplitude);
        }

        return m;
    }

    /**
     * Copy the contents of src to dst
     */
   public void quantum_copy_qureg(QuantumReg src, QuantumReg dst) {
   
        dst.hash = Arrays.copyOf(src.hash, src.hash.length);
        dst.hashw = src.hashw;
        dst.node = Arrays.copyOf(src.node, src.node.length);
        dst.size = src.size;
        dst.width = src.width;
    }

    /**
     * Print the contents of a quantum register to stdout
     */
   public void quantum_print_qureg(QuantumReg reg) {
        
        int i, j;

        for (i = 0; i < reg.size; i++) {
            System.out.format("%f %fi|%d> (%f) (|", 
                reg.node[i].amplitude.re(),
                reg.node[i].amplitude.im(), 
                reg.node[i].state,
                reg.node[i].amplitude.squareSum());
            
            //write bitstring of node's state from msb to lsb
            for (j = reg.width - 1; j >= 0; j--) {
                if (j % 4 == 3) {
                    System.out.format(" ");
                }
                // test if bit j is set in node[i]'s state:
                long bitJ = 1L << j;
                int b = ((bitJ & reg.node[i].state) > 0) ? 1 : 0;
                System.out.format("%d", b);
            }

            System.out.format(">)\n");
        }

        System.out.format("\n");
    }

    /**
     * Print the output of the modular exponentiation algorithm
     */
   public void quantum_print_expn(QuantumReg reg) {
        int i;
        for (i = 0; i < reg.size; i++) {
            long tmp = 1L << (reg.width / 2);
            System.out.format("%d: %d\n", i, 
                reg.node[i].state - i * tmp);
        }
    }

    /**
     * Add additional space to a qureg. It is initialized to zero and 
     * can be used as scratch space. Note that the space gets added at the LSB
     makes current bitstrings in reg.node states larger by a 
     * left bitshift of size bits.
     */
   public void quantum_addscratch(int bits, QuantumReg reg) {
        
        int i, oldwidth;
        //MAX_UNSIGNED l;
        long l;

        oldwidth = reg.width;
        reg.width += bits;

        long shift = bits;
        
        for (i = 0; i < reg.size; i++) {
            
            l = reg.node[i].state << shift;
            
            //System.out.format("-->%d (%d)\n", (int)l, i);
            
            reg.node[i].state = l;
        }
    }

    /**
     * Print the hash table to stdout and test if the hash table is corrupted
     */
   public void quantum_print_hash(QuantumReg reg) {
        int i;
        long tmp = 1L << reg.hashw;
        assert(reg.hash.length == tmp);
        int hashMax = MiscMath.findMax(reg.hash);
        System.out.println("hash length=" + reg.hash.length + " max value=" 
            + hashMax + " node.length=" + reg.node.length);
        assert(reg.node.length >= reg.size);
        for (i = 0; i < tmp; i++) {
            if (i > 0 && reg.hash[i] > 0) {
                int idx = reg.hash[i] - 1;
                if (idx < reg.node.length) {
                    System.out.format("%d: %d %d\n", i, idx,
                        reg.node[idx].state);
                }
            }
        }
    }

    /**
     * Compute the Kronecker product of two quantum registers.
     * 
     * <pre>
     *   if reg1 is A which is a 2^n vector and 
     *   reg2 is B which is a 2^m vector,
     *  
     *      the kronecker product is
     * 
     *      |A> ⊗ |B> is  | A_1 * B_1         |
     *                    | A_1 * B_2         |
     *                       ...
     *                    | A_1 * B_(2^m)     |
     *                        ...
     *                    | A_2 * B_1         |
     *                        ...
     *                    | A_(2^n) * B_(2^m) |
     * </pre>
     *
     */
    QuantumReg quantum_kronecker(QuantumReg reg1, QuantumReg reg2) {
  
        int i, j;
        QuantumReg reg = new QuantumReg();

        reg.width = reg1.width + reg2.width;
        reg.size = reg1.size * reg2.size;
        reg.hashw = reg.width + 2;

        // allocate memory for the new basis states 
        reg.node = new QuantumRegNode[reg.size];
        for (int ii = 0; ii < reg.node.length; ii++) {
            reg.node[ii] = new QuantumRegNode();
            reg.node[ii].state = 0;
            reg.node[ii].amplitude = new ComplexModifiable(0, 0);
        }
        
        // Allocate the hash table 
        int nHash = 1 << reg.hashw;
        reg.hash = new int[nHash];
        
        for (i = 0; i < reg1.size; i++) {
            for (j = 0; j < reg2.size; j++) {
                // printf("processing |%lli> x |%lli>\n", reg1->node[i].state, 
                //reg2 -> node[j].state
                //);
                //printf("%lli\n", (reg1 -> node[i].state) << reg2 -> width);
                
                long tmp = reg1.node[i].state << reg2.width;
                reg.node[i * reg2.size + j].state = ( 
                    tmp | reg2.node[j].state);
                
                reg.node[i * reg2.size + j].amplitude.resetTo(reg1.node[i].amplitude);
                reg.node[i * reg2.size + j].amplitude.times(reg2.node[j].amplitude);
            }
        }

        return reg;
    }

    /**
     * Reduce the state vector after measurement or partial trace
     */
    QuantumReg quantum_state_collapse(int pos, int value, QuantumReg reg) {
        int i, j, k;
        int size = 0;
        double d = 0;
        
        QuantumReg out = new QuantumReg();
        
        //MAX_UNSIGNED lpat = 0, rpat = 0, pos2;
        //pos2 = (MAX_UNSIGNED) 1 << pos;
        long pos2 = 1L << pos;
        
        // Eradicate all amplitudes of base states which have been ruled out
        //   by the measurement and get the norm of the new register 
        for (i = 0; i < reg.size; i++) {
        
            long posBit = reg.node[i].state & pos2;
            
            if (((posBit != 0) && (value != 0))
                || ((posBit == 0) && (value == 0))) {
                
                d += reg.node[i].amplitude.squareSum();
                size++;
            }
        }

        // Build the new quantum register 
        out.width = reg.width - 1;
        out.size = size;
        out.node = new QuantumRegNode[size];
        out.hashw = reg.hashw;
        out.hash = reg.hash;
        for (int ii = 0; ii < out.node.length; ii++) {
            out.node[ii] = new QuantumRegNode();
            out.node[ii].state = 0;
            out.node[ii].amplitude = new ComplexModifiable(0, 0);
        }
        
        // ---- if value == 0, keep the states with set bit pos
        //      else keep the states where bit pos is not set.
        //      d becomes the new normalization for the out states.
        
        long lpat = 0;
        long rpat = 0;
        
        // Determine the numbers of the new base states and norm 
        // the quantum register 
        for (i = 0, j = 0; i < reg.size; i++) {
            
            long posBit = reg.node[i].state & pos2;
            
            if (((posBit != 0) && (value != 0))
                || ((posBit == 0) && (value == 0))) {
             
                for (k = 0, rpat = 0; k < pos; k++) {
                    rpat += 1L << k;
                }

                rpat &= reg.node[i].state;
                                
                for (k = 63 - 1, lpat = 0; k > pos; k--) {
                    lpat += 1L << k;
                }

                lpat &= reg.node[i].state;

                out.node[j].state = ((lpat >> 1) | rpat);
                
                out.node[j].amplitude.resetTo(reg.node[i].amplitude);
                out.node[j].amplitude.times(1./Math.sqrt(d));
                
                j++;
            }
        }

        return out;
    }

    /**
     * Compute the dot product of two quantum registers.
     * 
     * <pre>
     * |psi> = summation over j of alpha*_j * |j>
     * 
     * |phi> = summation over j of beta*_j * |j>
     * 
     * dot product
     * (psi|phi> = summation over j of alpha*_j * beta_j
     * 
     * </pre>
     *
     */
    ComplexModifiable quantum_dot_product(QuantumReg reg1, QuantumReg reg2) {
        
        int i, j;
        ComplexModifiable f = new ComplexModifiable(0, 0);

        // Check whether quantum registers are sorted 
        if (reg2.hashw > 0) {
            quantum_reconstruct_hash(reg2);
        }
        
        for (i = 0; i < reg1.size; i++) {
            
            j = quantum_get_state(reg1.node[i].state, reg2);

            // state exists in reg2 
            if (j > -1) {
                ComplexModifiable tmp2 = reg1.node[i].amplitude.copy();
                tmp2.conjugate();
                tmp2.times(reg2.node[j].amplitude);
                
                f.plus(tmp2);
            }
        }
        return f;
    }

    /**
     * Same as above, but without complex conjugation
     */
    ComplexModifiable quantum_dot_product_noconj(QuantumReg reg1, QuantumReg reg2) {
  
        int i, j;
        ComplexModifiable f = new ComplexModifiable(0, 0);

        // Check whether quantum registers are sorted 
        if (reg2.hashw != 0) {
            quantum_reconstruct_hash(reg2);
        }

        for (i = 0; i < reg1.size; i++) {
            j = quantum_get_state(reg1.node[i].state, reg2);

            // state exists in reg2 
            if (j > -1) {
                ComplexModifiable tmp2 = reg1.node[i].amplitude.copy();
                tmp2.times(reg2.node[j].amplitude);
                f.plus(tmp2);
            }
        }
        return f;
    }

    /**
     * Vector addition of two quantum registers. This is a purely mathematical
     * operation without any physical meaning, so only use it if you know what
     * you are doing.
     */
    QuantumReg quantum_vectoradd(QuantumReg reg1, QuantumReg reg2) {
  
        int i, j, k;
        int addsize = 0;
        QuantumReg reg = new QuantumReg();

        quantum_copy_qureg(reg1, reg);

        if ((reg1.hashw != 0) || (reg2.hashw != 0)) {
            quantum_reconstruct_hash(reg1);
            quantum_copy_qureg(reg1, reg);

            // Calculate the number of additional basis states 
            for (i = 0; i < reg2.size; i++) {
                if (quantum_get_state(reg2.node[i].state, reg1) == -1) {
                    addsize++;
                }
            }
        }

        reg.size += addsize;
        if (reg.node.length != reg.size) {
            int len1 = reg.node.length;
            reg.node = Arrays.copyOf(reg.node, reg.size);
            for (i = len1; i < reg.size; i++) {
                reg.node[i] = new QuantumRegNode();
                reg.node[i].state = 0;
                reg.node[i].amplitude = new ComplexModifiable(0, 0);
            }
        }
        
        k = reg1.size;

        for (i = 0; i < reg2.size; i++) {
            j = quantum_get_state(reg2.node[i].state, reg1);

            if (j >= 0) {
                reg.node[j].amplitude.plus(reg2.node[i].amplitude);
            } else {
                reg.node[k].state = reg2.node[i].state;
                reg.node[k].amplitude.resetTo(reg2.node[i].amplitude);
                k++;
            }
        }

        return reg;
    }

    /**
     * Same as above, but the result is stored in the first register
     */
   public void quantum_vectoradd_inplace(QuantumReg reg1, QuantumReg reg2) {
  
        int i, j, k;
        int addsize = 0;

        if ((reg1.hashw != 0) || (reg2.hashw != 0)) {
            quantum_reconstruct_hash(reg1);

            // Calculate the number of additional basis states 
            for (i = 0; i < reg2.size; i++) {
                if (quantum_get_state(reg2.node[i].state, reg1) == -1) {
                    addsize++;
                }
            }
        }
        // Allocate memory for basis states 

        if (reg1.node.length != (reg1.size + addsize)) {
            int len1 = reg1.node.length;
            reg1.node = Arrays.copyOf(reg1.node, reg1.size + addsize);
            for (int ii = len1; ii < (reg1.size + addsize); ii++) {
                reg1.node[ii] = new QuantumRegNode();
                reg1.node[ii].state = 0;
                reg1.node[ii].amplitude = new ComplexModifiable(0, 0);
            }
        }
        
        // Allocate the hash table 
        k = reg1.size;

        for (i = 0; i < reg2.size; i++) {
            j = quantum_get_state(reg2.node[i].state, reg1);

            if (j >= 0) {
                reg1.node[j].amplitude.plus(reg2.node[i].amplitude);
            } else {
                reg1.node[k].state = reg2.node[i].state;
                reg1.node[k].amplitude.resetTo(reg2.node[i].amplitude);
                k++;
            }
        }
        reg1.size += addsize;
    }

    /**
     * Matrix-vector multiplication for a quantum register. A is a function
     * returning a quantum register containing the row given in the first
     * parameter. An additional parameter (e.g. time) may be supplied as well.
     */
    /*
    QuantumReg quantum_matrix_qureg(
        QuantumReg A (MAX_UNSIGNED, double), double t, QuantumReg reg) {
  
        MAX_UNSIGNED i;
        QuantumReg reg2 = new QuantumReg();
        QuantumReg tmp = new QuantumReg();

        reg2.width = reg.width;
        reg2.size = 1 << reg2.width;
        reg2.hashw = 0;
        reg2.hash = null;

        reg2.node = calloc(reg2.size, sizeof(quantum_reg_node));
        if (!reg2.node) {
            quantum_error(QUANTUM_ENOMEM);
        }

        quantum_memman(reg2.size * sizeof(quantum_reg_node));

        for (i = 0; i < (1 << reg -> width); i++) {
            reg2.node[i].state = i;
            tmp = A(i, t);
            reg2.node[i].amplitude = quantum_dot_product_noconj( & tmp, reg);
            quantum_delete_qureg( & tmp);
        }

        return reg2;

    }
    */
    
    /**
     * Scalar multiplication of a quantum register. This is a purely
     * mathematical operation without any physical meaning, so only use it if
     * you know what you are doing.
     */
   public void quantum_scalar_qureg(ComplexModifiable r, QuantumReg reg) {
  
        int i;

        for (i = 0; i < reg.size; i++) {
            reg.node[i].amplitude.times(r);
        }
    }
}
