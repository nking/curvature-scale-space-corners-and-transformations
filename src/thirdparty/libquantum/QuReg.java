package thirdparty.libquantum;

import algorithms.misc.ComplexModifiable;
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
    static int quantum_hash64(long key, int width) {
        //unsigned int k32;
        int k32;
        
        //             32 bit mask
        //k32 = (key & 0xFFFFFFFF) ^ (key >> 32);
        //             31 bit mask
        k32 = (int)((key & 0x7FFFFFFF) ^ (key >> 32));
        
        //     31.3057 bits
        //k32 *= 0x9e370001UL;
        
        //     unsigned int = unsigned int * unsigned long
        //        32 bits        32 bits        64 bits
        // in c, the righthand side are both converted to unsigned long
        //       then result is narrowed to unsigned int.
        //       overflow from unsigned long to unsigned int is truncated
        //       so max value is 1<<32 bits - 1
        
        //TODO: check the reduction of 0x9e370001UL by a bit:
        //   is the result a decent magic number at least as good as the previous?
        //k32 *= 1327202304.5;
        long tmp = k32 * 1327202304L;
        if (tmp > Integer.MAX_VALUE) {
            k32 = Integer.MAX_VALUE;
        } else {
            k32 = (int)tmp;
        }
        
        k32 = k32 >> (32 - width);
        
        return k32;
    }

    /**
     * Get the position of a given base state via the hash table
     */
    //static inline int quantum_get_state(MAX_UNSIGNED a, quantum_reg reg) {
    static int quantum_get_state(int a, QuantumReg reg) {
        int i;

        if (reg.hashw == 0) {
            return a;
        }

        i = quantum_hash64(a, reg.hashw);

        while (reg.hash[i] > 0) {
            if (reg.node[reg.hash[i] - 1].getState() == a) {
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
    static void quantum_add_hash(int a, int pos, QuantumReg reg) {

        int i, mark = 0;

        i = quantum_hash64(a, reg.hashw);

        while (reg.hash[i] > 0) {
            i++;
            // if i is > last index
            if (i == (1 << reg.hashw)) {
                if (mark == 0) {
                    i = 0;
                    mark = 1;
                } else {
                    throw new IllegalStateException("hash is full");
                }
            }
        }
        reg.hash[i] = pos + 1;

        System.out.format("   added hash i=%d\n", i);
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

        int nzd = 0;
        for (i = 0; i < (1 << reg.hashw); i++) {
            reg.hash[i] = 0;
            ++nzd;
        }
        System.out.format(" quantum_reconstruct_hash: zeroed TOTAL n=%d\n", nzd);

        for (i = 0; i < reg.size; i++) {
            quantum_add_hash(reg.node[i].getState(), i, reg);
        }
    }
    
    /**
     * convenienc method to handle c truncation due to a left shift
     * that produces a number larger than integer.max_value.
     * @param shift
     * @return 
     */
    public static int shiftLeftTruncate(int shift) {
        int tmp;
        if (shift > 30) {
            tmp = Integer.MAX_VALUE;
        } else {
            tmp = 1 << shift;
        }
        return tmp;
    }

    /**
     * Return the reduced bitmask of a basis state
     */
    //static int quantum_bitmask(MAX_UNSIGNED a, int width, int  *bits) {
    static int quantum_bitmask(int a, int width, int[] bits) {
        
        int i;
        int mask = 0;

        for (i = 0; i < width; i++) {
            //if (a & ((MAX_UNSIGNED) 1 << bits[i])) {
            int tmp = shiftLeftTruncate(bits[i]);
            if ((a & tmp) > 0) {
                mask += 1 << i;
            }
        }

        return mask;
    }

    /**
     * Convert a vector to a quantum register
     */
    QuantumReg quantum_matrix2qureg(QuantumMatrix m, int width) {
  
        QuantumReg reg = new QuantumReg();
        int i, j, size = 0;

        if (m.cols != 1) {
            throw new IllegalArgumentException("m nCols must == 1");
        }

        reg.width = width;

        // Determine the size of the quantum register 
        for (i = 0; i < m.rows; i++) {
            if (m.t[i] != null) {
                size++;
            }
        }

        // Allocate the required memory 
        reg.size = size;
        reg.hashw = width + 2;

        reg.node = new QuantumRegNode[size];
        for (int ii = 0; ii < reg.node.length; ii++) {
            reg.node[ii] = new QuantumRegNode();
            reg.node[ii].setState(0);
            reg.node[ii].amplitude = new ComplexModifiable(0, 0);
        }
        
        // Allocate the hash table 
        int nHash = shiftLeftTruncate(reg.hashw);
        reg.hash = new int[nHash];

        // Copy the nonzero amplitudes of the vector into the 
        //quantum register 
        for (i = 0, j = 0; i < m.rows; i++) {
            if (m.t[i] != null) {
                reg.node[j].setState(i);
                reg.node[j].amplitude = m.t[i];
                j++;
            }
        }

        return reg;
    }

    /**
     * Create a new quantum register from scratch
     */
    //QuantumReg quantum_new_qureg(MAX_UNSIGNED initval, int width) {
    QuantumReg quantum_new_qureg(int initval, int width) {
        
        QuantumReg reg = new QuantumReg();

        reg.width = width;
        reg.size = 1;
        reg.hashw = width + 2;

        // Allocate memory for 1 base state 
        reg.node = new QuantumRegNode[1];

        // Allocate the hash table 
        int nHash = shiftLeftTruncate(reg.hashw);
        reg.hash = new int[nHash];

        // Initialize the quantum register 
        reg.node[0] = new QuantumRegNode();
        reg.node[0].setState(initval);
        reg.node[0].amplitude =  new ComplexModifiable(1, 0);

        System.err.format(
            "init reg: %d qubits, 1 node, and %d hash table length, hashw=%d\n",
            reg.width, nHash, reg.hashw);

        return reg;
    }

    /**
     * Returns an empty quantum register of size N
     */
    QuantumReg quantum_new_qureg_size(int n, int width) {
        
        QuantumReg reg = new QuantumReg();

        reg.width = width;
        reg.size = n;
        reg.hashw = 0;
        reg.hash = null;

        // Allocate memory for n basis states 
        reg.node = new QuantumRegNode[n];
        for (int ii = 0; ii < reg.node.length; ii++) {
            reg.node[ii] = new QuantumRegNode();
            reg.node[ii].setState(0);
            reg.node[ii].amplitude = new ComplexModifiable(0, 0);
        }
        return reg;
    }

    /**
     * Convert a quantum register to a vector
     */
    QuantumMatrix quantum_qureg2matrix(QuantumReg reg) {
        
        QuantumMatrix m = new QuantumMatrix();
        int i;

        
        Matrix matrix = new Matrix();
        int sz = shiftLeftTruncate(reg.width);
        m = matrix.quantum_new_matrix(1, sz);

        for (i = 0; i < reg.size; i++) {
            m.t[reg.node[i].getState()] = reg.node[i].amplitude;
        }

        return m;
    }

    /**
     * Destroys the entire hash table of a quantum register
     */
    void quantum_destroy_hash(QuantumReg reg) {
        reg.hash = null;
    }

    /**
     * Delete a quantum register
     */
    void quantum_delete_qureg(QuantumReg reg) {
        reg.size = 0;
        reg.hash = null;
        reg.width = 0;
        reg.hashw = 0;
        reg.node = null;
    }

    /**
     * Copy the contents of src to dst
     */
    void quantum_copy_qureg(QuantumReg src, QuantumReg dst) {
   
        dst.hash = Arrays.copyOf(src.hash, src.hash.length);
        dst.hashw = src.hashw;
        dst.node = Arrays.copyOf(src.node, src.node.length);
        dst.size = src.size;
        dst.width = src.width;
    }

    /**
     * Print the contents of a quantum register to stdout
     */
    void quantum_print_qureg(QuantumReg reg) {
        
        int i, j;

        for (i = 0; i < reg.size; i++) {
            System.out.format("% f %+fi|%lli> (%e) (|", 
                reg.node[i].amplitude.re(),
                reg.node[i].amplitude.im(), 
                reg.node[i].getState(),
                reg.node[i].amplitude.abs());
            for (j = reg.width - 1; j >= 0; j--) {
                if (j % 4 == 3) {
                    System.out.format(" ");
                }
                int tmp = shiftLeftTruncate(j);
                System.out.format("%d", ((tmp & reg.node[i].getState()) > 0));
            }

            System.out.format(">)\n");
        }

        System.out.format("\n");
    }

    /**
     * Print the output of the modular exponentation algorithm
     */
    void quantum_print_expn(QuantumReg reg) {
        int i;
        for (i = 0; i < reg.size; i++) {
            int tmp = shiftLeftTruncate(reg.width / 2);
            System.out.format("%d: %lli\n", i, 
                reg.node[i].getState() - i * tmp);
        }
    }

    /**
     * Add additional space to a qureg. It is initialized to zero and can be
     * used as scratch space. Note that the space gets added at the LSB
     */
    void quantum_addscratch(int bits, QuantumReg reg) {
        
        int i, oldwidth;
        //MAX_UNSIGNED l;
        int l;

        oldwidth = reg.width;
        reg.width += bits;

        for (i = 0; i < reg.size; i++) {
            l = shiftLeftTruncate(reg.node[i].getState() << bits);
            reg.node[i].setState(l);
        }
    }

    /**
     * Print the hash table to stdout and test if the hash table is corrupted
     */
    void quantum_print_hash(QuantumReg reg) {
        int i;
        int tmp = shiftLeftTruncate(reg.hashw);
        for (i = 0; i < tmp; i++) {
            if (i > 0) {
                System.out.format("%d: %d %llu\n", i, reg.hash[i] - 1,
                    reg.node[reg.hash[i] - 1].getState());
            }
        }
    }

    /**
     * Compute the Kronecker product of two quantum registers
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
            reg.node[ii].setState(0);
            reg.node[ii].amplitude = new ComplexModifiable(0, 0);
        }
        
        // Allocate the hash table 
        int nHash = shiftLeftTruncate(reg.hashw);
        reg.hash = new int[nHash];
        
        for (i = 0; i < reg1.size; i++) {
            for (j = 0; j < reg2.size; j++) {
                // printf("processing |%lli> x |%lli>\n", reg1->node[i].state, 
                //reg2 -> node[j].state
                //);
                //printf("%lli\n", (reg1 -> node[i].state) << reg2 -> width);
                
                int tmp = shiftLeftTruncate((reg1.node[i].getState()) << reg2.width);
                reg.node[i * reg2.size + j].setState( 
                    tmp | reg2.node[j].getState());
                
                ComplexModifiable tmp2 = reg1.node[i].amplitude.copy();
                tmp2.times(reg2.node[j].amplitude);
                reg.node[i * reg2.size + j].amplitude = tmp2;
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
        int lpat = 0;
        int rpat = 0;
        int pos2 = shiftLeftTruncate(pos);
        
        // Eradicate all amplitudes of base states which have been ruled out
        //   by the measurement and get the norm of the new register 
        for (i = 0; i < reg.size; i++) {
            if ((((reg.node[i].getState() & pos2) == 0) && (value == 0))
                || (((reg.node[i].getState() & pos2) == 0) && (value == 0))) {
                
                d += reg.node[i].amplitude.abs();
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
            out.node[ii].setState(0);
            out.node[ii].amplitude = new ComplexModifiable(0, 0);
        }
        
        // Determine the numbers of the new base states and norm 
        // the quantum register 
        for (i = 0, j = 0; i < reg.size; i++) {
            if ((((reg.node[i].getState() & pos2) == 0) && (value == 0))
                || (((reg.node[i].getState() & pos2) == 0) && (value == 0))) {
                for (k = 0, rpat = 0; k < pos; k++) {
                    rpat += shiftLeftTruncate(k);
                }

                rpat &= reg.node[i].getState();

                /* in the c code compiled on 84 bit platform,
                   MAX_UNSIGNED is an unsigned long long
                   but since the sizeof gives that type in bytes
                   then k = number of bytes * 8,
                   k is the number of bits
                for (k = sizeof(MAX_UNSIGNED) * 8 - 1, lpat = 0; k > pos; k--) {
                    lpat += (MAX_UNSIGNED) 1 << k;
                }
                */
                
                for (k = 31 - 1, lpat = 0; k > pos; k--) {
                    lpat += shiftLeftTruncate(k);
                }

                lpat &= reg.node[i].getState();

                out.node[j].setState((lpat >> 1) | rpat);
                
                ComplexModifiable tmp2 = reg.node[i].amplitude.copy();
                tmp2.times(1./Math.sqrt(d));
                out.node[j].amplitude = tmp2;

                j++;
            }
        }

        return out;
    }

    /**
     * Compute the dot product of two quantum registers
     */
    ComplexModifiable quantum_dot_product(QuantumReg reg1, QuantumReg reg2) {
        
        int i, j;
        ComplexModifiable f = new ComplexModifiable(0, 0);

        // Check whether quantum registers are sorted 
        if (reg2.hashw > 0) {
            quantum_reconstruct_hash(reg2);
        }
        
        for (i = 0; i < reg1.size; i++) {
            
            j = quantum_get_state(reg1.node[i].getState(), reg2);

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
        if (reg2.hashw > 0) {
            quantum_reconstruct_hash(reg2);
        }

        for (i = 0; i < reg1.size; i++) {
            j = quantum_get_state(reg1.node[i].getState(), reg2);

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

        if ((reg1.hashw > 0) || (reg2.hashw > 0)) {
            quantum_reconstruct_hash(reg1);
            quantum_copy_qureg(reg1, reg);

            // Calculate the number of additional basis states 
            for (i = 0; i < reg2.size; i++) {
                if (quantum_get_state(reg2.node[i].getState(), reg1) == -1) {
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
                reg.node[i].setState(0);
                reg.node[i].amplitude = new ComplexModifiable(0, 0);
            }
        }
        
        k = reg1.size;

        for (i = 0; i < reg2.size; i++) {
            j = quantum_get_state(reg2.node[i].getState(), reg1);

            if (j >= 0) {
                reg.node[j].amplitude.plus(reg2.node[i].amplitude);
            } else {
                reg.node[k].setState(reg2.node[i].getState());
                reg.node[k].amplitude = reg2.node[i].amplitude;
                k++;
            }
        }

        return reg;
    }

    /**
     * Same as above, but the result is stored in the first register
     */
    void quantum_vectoradd_inplace(QuantumReg reg1, QuantumReg reg2) {
  
        int i, j, k;
        int addsize = 0;

        if ((reg1.hashw > 0) || (reg2.hashw > 0)) {
            quantum_reconstruct_hash(reg1);

            // Calculate the number of additional basis states 
            for (i = 0; i < reg2.size; i++) {
                if (quantum_get_state(reg2.node[i].getState(), reg1) == -1) {
                    addsize++;
                }
            }
        }
        // Allocate memory for basis states 

        if (reg1.node.length != reg1.size) {
            int len1 = reg1.node.length;
            reg1.node = Arrays.copyOf(reg1.node, reg1.size);
            for (i = len1; i < reg1.size; i++) {
                reg1.node[i] = new QuantumRegNode();
                reg1.node[i].setState(0);
                reg1.node[i].amplitude = new ComplexModifiable(0, 0);
            }
        }
        
        // Allocate the hash table 
        k = reg1.size;

        for (i = 0; i < reg2.size; i++) {
            j = quantum_get_state(reg2.node[i].getState(), reg1);

            if (j >= 0) {
                reg1.node[j].amplitude.plus(reg2.node[i].amplitude);
            } else {
                reg1.node[k].setState(reg2.node[i].getState());
                reg1.node[k].amplitude = reg2.node[i].amplitude;
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
    void quantum_scalar_qureg(ComplexModifiable r, QuantumReg reg) {
  
        int i;

        for (i = 0; i < reg.size; i++) {
            reg.node[i].amplitude.times(r);
        }
    }
}
