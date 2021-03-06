package thirdparty.libquantum;

import algorithms.misc.ComplexModifiable;
import java.util.Arrays;
import java.util.Random;

/* gates.c and qec.c: 
Basic gates for quantum register manipulation and 
Quantum Error Correction.

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
public class Gates {

    private final float epsilon = 0.000001f;
    
    private int counter = 0;
    
    private int counterQEC = 0;
    private static int freq = (1 << 30);
    
    /**
     * Type of the QEC. Currently implemented versions are: 
     * 0: no QEC (default)
     * 1: Steane's 3-bit code
     */
    private int type = 0;

    /**
     * How many qubits are protected
     */
    private int width = 0;

    private final Random rng;
    
    private final Decoherence decoherence = new Decoherence();
    
    public Gates(Random r) {
        this.rng = r;
    }
    
    /**
     * Apply a controlled-not gate.
     * for each reg.state, 
          Flip the target bit of a basis state if 
          the control bit is set 
     */
   public void quantum_cnot(int control, int target, QuantumReg reg) {
  
        int i;
        int[] qec = new int[1];

        quantum_qec_get_status(qec, null);
        
        if (qec[0] != 0) {
            quantum_cnot_ft(control, target, reg);
        } else {
            
            for (i = 0; i < reg.size; i++) {
                // Flip the target bit of a basis state if 
                //the control bit is set 

                if ((reg.node[i].state & (1L << control)) != 0) {
                    reg.node[i].state ^= (1L << target);
                }
            }
            decoherence.quantum_decohere(reg, rng, this);
        }
    }

    /**
     * Apply a toffoli (or controlled-controlled-not) gate.
     * <pre>
     * this is a ternary gate with the property 
     * CCNOT(x, y,z) = (x, y, (x ∧ y) ⊕ z).
       "From Reversible Logic Gates to Universal Quantum Bases"
          by Bocharov and Svore
           
       From wikipedia
          https://en.wikipedia.org/wiki/Toffoli_gate
         
       If the first two bits are set, it flips the third bit.
            Truth table
        INPUT	    OUTPUT
        0 	0 	0   0   0   0 
        0	0	1	0	0	1
        0	1	0	0	1	0
        0	1	1	0	1	1
        1	0	0	1	0	0
        1	0	1	1	0	1
        1	1	0	1	1	1
        1	1	1	1	1	0
     * </pre>
     * runtime complexity is O(reg.size * reg.width),
       but decoherence could be improved for case where lambda==0
     */
   public void quantum_toffoli(int control1, int control2, int target, QuantumReg reg) {

        int i;
        int[] qec = new int[1];

        quantum_qec_get_status(qec, null);

        if (qec[0] != 0) {
            quantum_toffoli_ft(control1, control2, target, reg);
        } else {

            for (i = 0; i < reg.size; i++) {
                // Flip the target bit of a basis state if 
                // both control bits are set 

                long st = reg.node[i].state;
                if ((st & (1L << control1)) != 0) {
                    if ((st & (1L << control2)) != 0) {
                        reg.node[i].state ^= (1L << target);
                    }
                }
            }
            decoherence.quantum_decohere(reg, rng, this);
        }
    }

    /**
     * Apply an unbounded toffoli gate. This gate is not considered elementary
     * and is not available on all physical realizations of a quantum computer.
     * Be sure to pass the function the correct number of controlling qubits.
     * The target is given in the last argument.
     */
    /*
   public void quantum_unbounded_toffoli(int controlling, QuantumReg reg, ...) {
  
        va_list bits;
        int target;
        int *controls;
        int i, j;

        controls = malloc(controlling * sizeof(int));

        if (!controls) {
            quantum_error(QUANTUM_ENOMEM);
        }

        quantum_memman(controlling * sizeof(int));

        va_start(bits, reg);

        for (i = 0; i < controlling; i++) {
            controls[i] = va_arg(bits, int);
        }

        target = va_arg(bits, int);

        va_end(bits);

        for (i = 0; i < reg.size; i++) {
            for (j = 0; (j < controlling)
                && (reg.node[i].state & (MAX_UNSIGNED) 1 << controls[j]); j++);

            // all control bits are set 
            if (j == controlling) {
                reg.node[i].state ^= ((MAX_UNSIGNED) 1 << target);
            }
        }

        free(controls);
        quantum_memman(-controlling * sizeof(int));

        quantum_decohere(reg);
    }
    */

    /**
     * Apply a sigma_x (or not) gate which is a Pauli spin 
     * operation for x axis.
     * runtime complexity is O(reg.width * reg.size).
     */
   public void quantum_sigma_x(final int target, QuantumReg reg) {
  
        int i;
        int[] qec = new int[1];

        quantum_qec_get_status(qec, null);

        if (qec[0] != 0) {
            //System.out.println("NEXT sigma_x_ft");
            quantum_sigma_x_ft(target, reg);
        } else {
            for (i = 0; i < reg.size; i++) {
                // Flip the target bit of each basis state 
                reg.node[i].state ^= (1L << target);
            }
            decoherence.quantum_decohere(reg, rng, this);
        }
    }

    /**
     * Apply a sigma_y gate which is a Pauli spin operation for the y axis.
     */
   public void quantum_sigma_y(int target, QuantumReg reg) {
  
        int i;

        for (i = 0; i < reg.size; i++) {
            // Flip the target bit of each basis state and 
            // multiply with 
            // +/- i 

            long tmp1 = 1L << target;
            reg.node[i].state ^= tmp1;

            if ((reg.node[i].state & tmp1) != 0) {
                double re = reg.node[i].amplitude.re();
                double im = reg.node[i].amplitude.im();
                // new im = re * i (no sign change)
                // new re = im * i (i*i == -1, so sign change)
                reg.node[i].amplitude.setReal(-1 * im);
                reg.node[i].amplitude.setImag(re);
            } else {
                double re = reg.node[i].amplitude.re();
                double im = reg.node[i].amplitude.im();
                // new im = re * -i (sign change)
                // new re = im * -i (no sign change)
                reg.node[i].amplitude.setReal(im);
                reg.node[i].amplitude.setImag(-1 * re);
            }
        }

        decoherence.quantum_decohere(reg, rng, this);
    }

    /**
     * Apply a sigma_z gate which is a Pauli spin operation for the z axis.
     */
   public void quantum_sigma_z(int target, QuantumReg reg) {
        
        int i;

        for (i = 0; i < reg.size; i++) {
            // Multiply with -1 if the target bit is set 

            if ((reg.node[i].state & (1L << target)) != 0) {
                reg.node[i].amplitude.times(-1);
            }
        }
        decoherence.quantum_decohere(reg, rng, this);
    }

    /**
     * Swap the first WIDTH bits of the quantum register. This is done
     * classically by renaming the bits, unless QEC is enabled.
     */
   public void quantum_swaptheleads(int width, QuantumReg reg) {
        int i, j;
        int pat1, pat2;
        int[] qec = new int[1];
        //MAX_UNSIGNED l;
        long l;

        quantum_qec_get_status(qec, null);

        if (qec[0] != 0) {
            for (i = 0; i < width; i++) {
                quantum_cnot(i, width + i, reg);
                quantum_cnot(width + i, i, reg);
                quantum_cnot(i, width + i, reg);
            }
        } else {
            for (i = 0; i < reg.size; i++) {

                // calculate left bit pattern 
                
                // TODO: revisit c narrowing
                pat1 = (int)(reg.node[i].state % (1L << width));
                
                //calculate right but pattern 
                pat2 = 0;

                for (j = 0; j < width; j++) {
                    pat2 += (reg.node[i].state & (1L << (width + j)));
                }

                // construct the new basis state 
                l = (reg.node[i].state - (pat1 + pat2));
                l += (pat1 << width);
                l += (pat2 >> width);
                reg.node[i].state = l;
            }
        }
    }

    /**
     * Swap WIDTH bits starting at WIDTH and 2*WIDTH+2 controlled by CONTROL.
     * runtime complexity is width * O(reg.size)
     */
   public void quantum_swaptheleads_omuln_controlled(int control, int width, 
        QuantumReg reg) {
   
        int i;

        // brief summary of 2-bit swap gate at
        // https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/FromReversibleLogicGatestoUniversalQuantumBases.pdf
        // "From Reversible Logic Gates to Universal Quantum Bases"
        //    by Bocharov and Svore
        // SWAP = CNOT[1, 2]CNOT[2, 1]CNOT[1, 2]
        for (i = 0; i < width; i++) {
            quantum_toffoli(control, width + i, 2 * width + i + 2, reg);
            quantum_toffoli(control, 2 * width + i + 2, width + i, reg);
            quantum_toffoli(control, width + i, 2 * width + i + 2, reg);
        }
    }

    /**
     * Apply the 2x2 matrix M to the target bit. M should be unitary.
     
     runtime complexity is O(reg.size), though reg.size is possibly increased
     and for the last decohere operation is
     is O(reg.size * reg.width)
     * @param target
     * @param m
     * @param reg
     */
   public void quantum_gate1(int target, QuantumMatrix m, QuantumReg reg) {
            
        int i, k;
        int addsize = 0;
        //COMPLEX_FLOAT t, tnot = 0;
        ComplexModifiable t = null;
        ComplexModifiable tnot = new ComplexModifiable(0, 0);
        float limit;
        int[] done;
        int j = 0;

        if ((m.cols != 2) || (m.rows != 2)) {
            throw new IllegalArgumentException("n.cols must == 2 and "
                + "m.rows must == 2");
        }
        
        assert(reg.size <= reg.node.length);
        
        assert(reg.hashw >= 0);
              
        if (reg.hashw > 0) {
                        
            QuReg.quantum_reconstruct_hash(reg);

            long shifted = 1L << target;
            // calculate the number of basis states to be added 
            for (i = 0; i < reg.size; i++) {
                // determine whether XORed basis state already exists 
                long tst = reg.node[i].state ^ shifted;
                if (QuReg.quantum_get_state(tst, reg) == -1) {
                    addsize++;
                }
            }
            
            // allocate memory for the new basis states 
            int sz2 = reg.size + addsize;
            if (reg.node.length != sz2) {
                reg.node = Arrays.copyOf(reg.node, sz2);
                assert(sz2 == reg.node.length);
            }

            for (i = 0; i < addsize; i++) {
                int idx2 = i + reg.size;
                if (reg.node[idx2] == null) {
                    reg.node[idx2] = new QuantumRegNode();
                    reg.node[idx2].amplitude = new ComplexModifiable(0, 0);
                } else {
                    reg.node[idx2].state = 0;
                    reg.node[idx2].amplitude.setReal(0);
                    reg.node[idx2].amplitude.setImag(0);
                }
            }
        }
        
        done = new int[reg.size + addsize];

        k = reg.size;

        limit = (1.0f / (1L << reg.width)) * epsilon;
        
        // perform the actual matrix multiplication 
        for (i = 0; i < reg.size; i++) {
            if (done[i] == 0) {
                // determine if the target of the basis state is set
                
                int iset = (int)(reg.node[i].state & (1L << target));
                tnot.setReal(0);
                tnot.setImag(0);

                //find state's node index for state where basis state is the
                //   opposite of node[i]'s basis state
                j = QuReg.quantum_get_state(
                    reg.node[i].state ^ (1L << target), reg);
                
                if (t == null) {
                    t = reg.node[i].amplitude.copy();
                } else {
                    t.resetTo(reg.node[i].amplitude);
                }
                
                if (j >= 0) {
                    tnot.resetTo(reg.node[j].amplitude);
                }
                
                
                if (iset != 0) {
                    ComplexModifiable t2 = m.t[2].copy();
                    t2.times(tnot);
                    ComplexModifiable t3 = m.t[3].copy();
                    t3.times(t);
                    t2.plus(t3);
                    reg.node[i].amplitude.resetTo(t2);
                } else {
                    ComplexModifiable t0 = m.t[0].copy();
                    t0.times(t);
                    ComplexModifiable t1 = m.t[1].copy();
                    t1.times(tnot);
                    t0.plus(t1);
                    reg.node[i].amplitude.resetTo(t0);
                }
                
                
                if (j >= 0) {
                
                    if (iset != 0) {
                        ComplexModifiable t0 = m.t[0].copy();
                        t0.times(tnot);
                        ComplexModifiable t1 = m.t[1].copy();
                        t1.times(t);
                        t0.plus(t1);
                        reg.node[j].amplitude.resetTo(t0);
                    } else {
                        ComplexModifiable t2 = m.t[2].copy();
                        t2.times(t);
                        ComplexModifiable t3 = m.t[3].copy();
                        t3.times(tnot);
                        t2.plus(t3);
                        reg.node[j].amplitude.resetTo(t2);
                    }
                                       
                } else {
                    // new basis state will be created 

                    if (((float)m.t[1].squareSum() == 0) && (iset != 0)) {
                        break;
                    }
                    if (((float)m.t[2].squareSum() == 0) && (iset == 0)) {
                        break;
                    }
                    reg.node[k].state = (reg.node[i].state ^ (1L << target));
                    
                    if (iset != 0) {
                        reg.node[k].amplitude.resetTo(m.t[1]);
                        reg.node[k].amplitude.times(t);
                    } else {
                        reg.node[k].amplitude.resetTo(m.t[2]);
                        reg.node[k].amplitude.times(t);
                    }
                    
                    k++;
                }

                if (j >= 0) {
                    done[j] = 1;
                }
            }
        }
        reg.size += addsize;
        
        assert(reg.size <= reg.node.length);
        
        // remove basis states with extremely small amplitude 
        if (reg.hashw != 0) {
            int decsize = 0;
            for (i = 0, j = 0; i < reg.size; i++) {
                if ((float)reg.node[i].amplitude.squareSum() < limit) {
                    j++;
                    decsize++;
                } else if (j != 0) {
                    reg.node[i - j].state = reg.node[i].state;
                    reg.node[i - j].amplitude.resetTo(reg.node[i].amplitude);
                }
            }

            if (decsize != 0) {
                reg.size -= decsize;
                int len1 = reg.node.length;
                reg.node = Arrays.copyOf(reg.node, reg.size);
                for (int ii = len1; ii < reg.size; ii++) {
                    reg.node[ii] = new QuantumRegNode();
                    reg.node[ii].state = 0;
                    reg.node[ii].amplitude = new ComplexModifiable(0, 0);
                }
            }
        }
        
        //runtime complexity is O(reg.size * reg.width)
        decoherence.quantum_decohere(reg, rng, this);        
    }

    /**
     * Apply the 4x4 matrix M to the bits TARGET1 and TARGET2. M should be
     * unitary. * Warning: code is mostly untested.
     */
   public void quantum_gate2(int target1, int target2, QuantumMatrix m,
        QuantumReg reg) {
        
        if ((m.cols != 4) || (m.rows != 4)) {
            throw new IllegalArgumentException("m.cols must == 4 and "
                + " m.rows must == 4");
        }

        int i, j, k, l;
        int addsize = 0, decsize = 0;
        ComplexModifiable[] psi_sub = new ComplexModifiable[4];
        for (i = 0; i < 4; ++i) {
            psi_sub[i] = new ComplexModifiable(0, 0);
        }
        int[] base = new int[4];
        int[] bits = new int[2];
        float limit;
        int[] done;

        // Build hash table 
        int end = (1 << reg.hashw);
        Arrays.fill(reg.hash, 0, end, 0);
      
        assert(reg.node.length >= reg.size);
        
        for (i = 0; i < reg.size; i++) {
            QuReg.quantum_add_hash(reg.node[i].state, i, reg);
        }
        
        // calculate the number of basis states to be added 
        for (i = 0; i < reg.size; i++) {
            if (QuReg.quantum_get_state(
                reg.node[i].state ^ (1L << target1), reg) == -1) {
                addsize++;
            }
            if (QuReg.quantum_get_state(
                reg.node[i].state ^ (1L << target2), reg) == -1) {
                addsize++;
            }
        }
        // allocate memory for the new basis states 

        if (reg.node.length != (reg.size + addsize)) {
            reg.node = Arrays.copyOf(reg.node, reg.size + addsize);
        }

        for (i = 0; i < addsize; i++) {
            int idx2 = i + reg.size;
            reg.node[idx2] = new QuantumRegNode();
            reg.node[idx2].state = 0;
            reg.node[idx2].amplitude = new ComplexModifiable(0, 0);
        }

        done = new int[reg.size + addsize];

        l = reg.size;

        // amplitudes smaller than this are removed below
        limit = (1.0f / (1L << reg.width)) / 1000000;
        
        bits[0] = target1;
        bits[1] = target2;

        // perform the actual matrix multiplication 
        for (i = 0; i < reg.size; i++) {
            if (done[i] == 0) {
                j = QuReg.quantum_bitmask(
                    reg.node[i].state, 2, bits);
                base[j] = i;
                long state = reg.node[i].state;
                long tmp2 = 1L << target2;
                long tmp1 = 1L << target1;
                base[j ^ 1] = QuReg.quantum_get_state(state ^ tmp2, reg);
                base[j ^ 2] = QuReg.quantum_get_state(state ^ tmp1, reg);
                base[j ^ 3] = QuReg.quantum_get_state(state ^ tmp1 ^ tmp2, reg);

                for (j = 0; j < 4; j++) {
                    if (base[j] == -1) {
                        base[j] = l;
                        //		  reg->node[l].state = reg->node[i].state
                        l++;
                    }
                    psi_sub[j].resetTo(reg.node[base[j]].amplitude);
                }

                for (j = 0; j < 4; j++) {
                    reg.node[base[j]].amplitude.setReal(0);
                    reg.node[base[j]].amplitude.setImag(0);
                    for (k = 0; k < 4; k++) {
                        ComplexModifiable tmp = Matrix.M(m, k, j).copy();
                        tmp.times(psi_sub[k]);
                        reg.node[base[j]].amplitude.plus(tmp);                        
                    }

                    done[base[j]] = 1;
                }
            }
        }
        reg.size += addsize;

        // remove basis states with extremely small amplitude 
        for (i = 0, j = 0; i < reg.size; i++) {
            if (reg.node[i].amplitude.squareSum() < limit) {
                j++;
                decsize++;
            } else if (j != 0) {
                reg.node[i - j].state = reg.node[i].state;
                reg.node[i - j].amplitude.resetTo(reg.node[i].amplitude);
            }
        }

        if (decsize != 0) {
            reg.size -= decsize;
            if (reg.node.length != reg.size) {
                int len1 = reg.node.length;
                reg.node = Arrays.copyOf(reg.node, reg.size);
                for (int ii = len1; ii < reg.size; ii++) {
                    reg.node[ii] = new QuantumRegNode();
                    reg.node[ii].state = 0;
                    reg.node[ii].amplitude = new ComplexModifiable(0, 0);
                }
            }
        }

        decoherence.quantum_decohere(reg, rng, this);
    }

    /**
     * Apply a hadamard gate.
     * this is a generalized class of fourier transforms.
     
     runtime complexity is O(reg.size), though reg.size is possibly increased
     and for the last decohere operation is
     is O(reg.size * reg.width).
    
     */
   public void quantum_hadamard(int target, QuantumReg reg) {
  
        Matrix matrix = new Matrix();
        
        QuantumMatrix m = matrix.quantum_new_matrix(2, 2);

        double a = Math.sqrt(1.0 / 2.);
        m.t[0] = new ComplexModifiable(a, 0);
        m.t[1] = new ComplexModifiable(a, 0);
        m.t[2] = new ComplexModifiable(a, 0);
        m.t[3] = new ComplexModifiable(-a, 0);

        quantum_gate1(target, m, reg);

    }

    /**
     * Apply a walsh-hadamard transform
     */
   public void quantum_walsh(int width, QuantumReg reg) {
  
        int i;

        for (i = 0; i < width; i++) {
            quantum_hadamard(i, reg);
        }
    }

    /**
     * Apply a rotation about the x-axis by the angle GAMMA
     */
   public void quantum_r_x(int target, float gamma, QuantumReg reg) {
  
        Matrix matrix = new Matrix();
        
        QuantumMatrix m = matrix.quantum_new_matrix(2, 2);

        m.t[0] = new ComplexModifiable(Math.cos(gamma / 2), 0);
        m.t[1] = new ComplexModifiable(0, -Math.sin(gamma / 2));
        m.t[2] = new ComplexModifiable(0, -Math.sin(gamma / 2));
        m.t[3] = new ComplexModifiable(Math.cos(gamma / 2), 0);

        quantum_gate1(target, m, reg);

    }

    /**
     * Apply a rotation about the y-axis by the angle GAMMA
     */
   public void quantum_r_y(int target, float gamma, QuantumReg reg) {
  
        Matrix matrix = new Matrix();
        
        QuantumMatrix m = matrix.quantum_new_matrix(2, 2);
        m.t[0] = new ComplexModifiable(Math.cos(gamma / 2), 0);
        m.t[1] = new ComplexModifiable(-Math.sin(gamma / 2), 0);
        m.t[2] = new ComplexModifiable(Math.sin(gamma / 2), 0);
        m.t[3] = new ComplexModifiable(Math.cos(gamma / 2), 0);

        quantum_gate1(target, m, reg);
    }

    /**
     * Apply a rotation about the z-axis by the angle GAMMA
     */
   public void quantum_r_z(int target, float gamma, QuantumReg reg) {
  
        int i;
        
        //return cos(phi) + IMAGINARY * sin(phi);
        double angle = gamma / 2;
        ComplexModifiable z = new ComplexModifiable(Math.cos(angle), Math.sin(angle));

        for (i = 0; i < reg.size; i++) {
            if ((reg.node[i].state & (1L << target)) != 0) {
                reg.node[i].amplitude.times(z);
            } else {
                reg.node[i].amplitude.divided(z);
            }
        }

        decoherence.quantum_decohere(reg, rng, this);
    }

    /**
     * Scale the phase of qubit
     */
   public void quantum_phase_scale(int target, float gamma, QuantumReg reg) {
  
        int i;
        //return cos(phi) + IMAGINARY * sin(phi);
        double angle = gamma;
        ComplexModifiable z = new ComplexModifiable(Math.cos(angle), Math.sin(angle));

        for (i = 0; i < reg.size; i++) {
            reg.node[i].amplitude.times(z);
        }

        decoherence.quantum_decohere(reg, rng, this);
    }

    /**
     * Apply a phase kick (== shift) by the angle GAMMA if target bit is set
     */
   public void quantum_phase_kick(int target, double gamma, QuantumReg reg) {
  
        int i;
        //return cos(phi) + IMAGINARY * sin(phi);
        double angle = gamma;
        ComplexModifiable z = new ComplexModifiable(Math.cos(angle), Math.sin(angle));

        for (i = 0; i < reg.size; i++) {
            if ((reg.node[i].state & (1L << target)) != 0){
                reg.node[i].amplitude.times(z);
            }
        }

        decoherence.quantum_decohere(reg, rng, this);
    }

    /**
     * Apply a conditional phase shift of PI / 2^(CONTROL - TARGET) if 
     * control and target bits are set
     */
   public void quantum_cond_phase(int control, int target, QuantumReg reg) {
  
        int i;
        
        //return cos(phi) + IMAGINARY * sin(phi);
        double angle = Math.PI / (1L << (control - target));
        ComplexModifiable z 
            = new ComplexModifiable(Math.cos(angle), Math.sin(angle));
        
        for (i = 0; i < reg.size; i++) {
            if ((reg.node[i].state & (1L << control)) != 0) {
                if ((reg.node[i].state & (1L << target)) != 0) {
                    reg.node[i].amplitude.times(z);
                }
            }
        }

        decoherence.quantum_decohere(reg, rng, this);
    
    }

   /**
     * Apply a conditional phase shift of -PI / 2^(CONTROL - TARGET) if 
     * control and target bits are set
     */
   public void quantum_cond_phase_inv(int control, int target, QuantumReg reg) {
  
        int i;
        //return cos(phi) + IMAGINARY * sin(phi);
        double angle = -Math.PI / (1L << (control - target));
        ComplexModifiable z = new ComplexModifiable(Math.cos(angle), Math.sin(angle));

        for (i = 0; i < reg.size; i++) {
            if ((reg.node[i].state & (1L << control)) != 0) {
                if ((reg.node[i].state & (1L << target)) != 0) {
                    reg.node[i].amplitude.times(z);
                }
            }
        }

        decoherence.quantum_decohere(reg, rng, this);
    }

   /**
     * Apply a conditional phase shift of gamma if 
     * control and target bits are set
     */
   public void quantum_cond_phase_kick(int control, int target, float gamma, 
        QuantumReg reg) {
  
        int i;
        //return cos(phi) + IMAGINARY * sin(phi);
        double angle = gamma;
        ComplexModifiable z = new ComplexModifiable(Math.cos(angle), Math.sin(angle));

        for (i = 0; i < reg.size; i++) {
            if ((reg.node[i].state & (1L << control)) != 0) {
                if ((reg.node[i].state & (1L << target)) != 0) {
                    reg.node[i].amplitude.times(z);
                }
            }
        }
        decoherence.quantum_decohere(reg, rng, this);
    }

    /**
     * Increase the gate counter by INC steps or reset it if INC .lt. 0. The
     * current value of the counter is returned.
     */
    int quantum_gate_counter(int inc) {
        //static int counter = 0;

        if (inc > 0) {
            counter += inc;
        } else if (inc < 0) {
            counter = 0;
        }

        return counter;
    }
    
    // --------- QEC methods -------
    
    /**
     * Change the status of the QEC.
     */
   public void quantum_qec_set_status(int stype, int swidth) {
        type = stype;
        width = swidth;
    }

    /**
     * Get the current QEC status
     */
   public void quantum_qec_get_status(int[] ptypeInOut, int[] pwidthInOut) {

        if (ptypeInOut != null && ptypeInOut.length > 0) {
            ptypeInOut[0] = type;
        }
        if (pwidthInOut != null && pwidthInOut.length > 0) {
            pwidthInOut[0] = width;
        }
    }

    /**
     * Encode a quantum register. All qubits up to SWIDTH are protected, the
     * rest is expanded with a repetition code.
     */
   public void quantum_qec_encode(int type, int width, QuantumReg reg) {
        int i;
        float lambda;

        lambda = decoherence.quantum_get_decoherence();

        decoherence.quantum_set_decoherence(0);

        for (i = 0; i < reg.width; i++) {
            if (i == reg.width - 1) {
                decoherence.quantum_set_decoherence(lambda);
            }

            if (i < width) {
                quantum_hadamard(reg.width + i, reg);
                quantum_hadamard(2 * reg.width + i, reg);

                quantum_cnot(reg.width + i, i, reg);
                quantum_cnot(2 * reg.width + i, i, reg);
            } else {
                quantum_cnot(i, reg.width + i, reg);
                quantum_cnot(i, 2 * reg.width + i, reg);
            }
        }

        quantum_qec_set_status(1, reg.width);
        reg.width *= 3;
    }

    /**
     * Decode a quantum register and perform Quantum Error Correction on it
     */
   public void quantum_qec_decode(int type, int width, QuantumReg reg) {

        int i, a, b;
        int swidth;
        float lambda;

        lambda = decoherence.quantum_get_decoherence();

        decoherence.quantum_set_decoherence(0);

        swidth = reg.width / 3;

        quantum_qec_set_status(0, 0);

        for (i = reg.width / 3 - 1; i >= 0; i--) {
            if (i == 0) {
                decoherence.quantum_set_decoherence(lambda);
            }

            if (i < width) {
                quantum_cnot(2 * swidth + i, i, reg);
                quantum_cnot(swidth + i, i, reg);

                quantum_hadamard(2 * swidth + i, reg);
                quantum_hadamard(swidth + i, reg);
            } else {
                quantum_cnot(i, 2 * swidth + i, reg);
                quantum_cnot(i, swidth + i, reg);
            }
        }
        
        Measure measure = new Measure();

        for (i = 1; i <= swidth; i++) {
            a = measure.quantum_bmeasure(swidth, reg, rng);
            b = measure.quantum_bmeasure(2 * swidth - i, reg, rng);
            if (a == 1 && b == 1 && i - 1 < width) {
                quantum_sigma_z(i - 1, reg);
            }
            /**
             * Z = HXH
             */
        }
    }

    /**
     * Counter which can be used to apply QEC periodically
     */
    int quantum_qec_counter(int inc, int frequency, QuantumReg reg) {

        //static int counter = 0;
        //static int freq = (1 << 30);
        if (inc > 0) {
            counterQEC += inc;
        } else if (inc < 0) {
            counterQEC = 0;
        }

        if (frequency > 0) {
            freq = frequency;
        }

        if (counterQEC >= freq) {
            counterQEC = 0;
            quantum_qec_decode(type, width, reg);
            quantum_qec_encode(type, width, reg);
        }

        return counterQEC;
    }

    /**
     * Fault-tolerant version of the NOT gate
     */
   public void quantum_sigma_x_ft(int target, QuantumReg reg) {
        int tmp;
        float lambda;

        tmp = type;
        type = 0;

        lambda = decoherence.quantum_get_decoherence();
        decoherence.quantum_set_decoherence(0);

        // These operations can be performed simultaneously
        quantum_sigma_x(target, reg);
        quantum_sigma_x(target + width, reg);
        decoherence.quantum_set_decoherence(lambda);
        quantum_sigma_x(target + 2 * width, reg);

        quantum_qec_counter(1, 0, reg);

        type = tmp;
    }

    /**
     * Fault-tolerant version of the Controlled NOT gate
     */
   public void quantum_cnot_ft(int control, int target, QuantumReg reg) {
        
        int tmp;
        float lambda;

        tmp = type;
        type = 0;
        
        // These operations can be performed simultaneously
        lambda = decoherence.quantum_get_decoherence();
        decoherence.quantum_set_decoherence(0);

        quantum_cnot(control, target, reg);
        quantum_cnot(control + width, target + width, reg);
        decoherence.quantum_set_decoherence(lambda);
        quantum_cnot(control + 2 * width, target + 2 * width, reg);

        quantum_qec_counter(1, 0, reg);

        type = tmp;

    }

    /**
     * Fault-tolerant version of the Toffoli gate
     */
   public void quantum_toffoli_ft(int control1, int control2, int target, 
        QuantumReg reg) {

        int i;
        int c1, c2;
        //MAX_UNSIGNED mask;
        long mask;
        
        mask = (1L << target)
            + (1L << (target + width))
            + (1L << (target + 2 * width));
        
        System.out.println("toffoli_ft mask=" + mask);
        
        for (i = 0; i < reg.size; i++) {
            c1 = 0;
            c2 = 0;

            if ((reg.node[i].state & (1L << control1)) != 0) {
                c1 = 1;
            }
            if ((reg.node[i].state
                & (1L << (control1 + width))) != 0) {
                c1 ^= 1;
            }
            if ((reg.node[i].state
                & (1L << (control1 + 2 * width))) != 0) {
                c1 ^= 1;
            }

            if ((reg.node[i].state & (1L << (control2))) != 0) {
                c2 = 1;
            }
            if ((reg.node[i].state & (1L << (control2 + width))) != 0) {
                c2 ^= 1;
            }
            if ((reg.node[i].state & (1L << (control2 + 2 * width))) != 0) {
                c2 ^= 1;
            }

            if (c1 == 1 && c2 == 1) {
                reg.node[i].state = (reg.node[i].state ^ mask);
            }
        }

        decoherence.quantum_decohere(reg, rng, this);

        quantum_qec_counter(1, 0, reg);
    }
    
    // ----- from expn.c, omuln.c, Multiplication modulo an integer N ----
    
    /**
       calculates f(a) = x^a mod N in the width of the working space.
       
       runtime complexity is runtime complexity is 
         width_input * swidth * O(reg.size).
     
     * @param N
     * @param x
     * @param width_input
     * @param swidth
     * @param reg 
     */
   public void quantum_exp_mod_n(int N, int x, int width_input, int swidth, 
        QuantumReg reg) {

        if (x == 0) {
            throw new IllegalArgumentException("x cannot == 0");
        }
        
        int i, j, f;
        
        //toggle bit '2 * width + 2' in each node state
        quantum_sigma_x(2 * swidth + 2, reg);
 
        f = x % N;
        int lastI = 1;
        for (i = 1; i <= width_input; i++) {
            f = x % N;			//compute
            for (j = lastI; j < i; j++) {
                f *= f;	//x^2^(i-1)
                f = f % N;
            }
            lastI = i;
            //apply f mod N = f - floor(f / N) * N
            mul_mod_n(N, f, 3 * swidth + 1 + i, swidth, reg);
        }
    }

   public void emul(int a, int L, int width, QuantumReg reg){

	    int i;
        for (i = width - 1; i >= 0; i--) {
            if (((a >> i) & 1) != 0) {
                quantum_toffoli(2 * width + 2, L, width + i, reg);
            }
        }
    }

   public void muln(int N, int a, int ctl, int w, QuantumReg reg) {
        //ctl tells, which bit is the external enable bit
	
        int i;
        int L = 2 * w + 1;

        quantum_toffoli(ctl, 2 * w + 2, L, reg);

        emul(a % N, L, w, reg);
        
        quantum_toffoli(ctl, 2 * w + 2, L, reg);
        
        for (i = 1; i < w; i++) {
        
            quantum_toffoli(ctl, 2 * w + 2 + i, L, reg);
            
            add_mod_n(N, ((1 << i) * a) % N, 
                w, reg);
           
            quantum_toffoli(ctl, 2 * w + 2 + i, L, reg);
        
        }
    }

    /**
     * runtime complexity is w * O(reg.size)
     * @param N
     * @param a
     * @param ctl
     * @param w
     * @param reg 
     */
   public void muln_inv(int N, int a, int ctl, int w, QuantumReg reg){

        //ctl tells, which bit is the external enable bit
	
        int i;
        int L = 2 * w + 1;

        Classic classic = new Classic();
        
        a = classic.quantum_inverse_mod(N, a);
 
        for (i = w - 1; i > 0; i--) {
            quantum_toffoli(ctl, 2 * w + 2 + i, L, reg);
             
            add_mod_n(N, N - ((1 << i) * a) % N, 
                w, reg);           
            
            quantum_toffoli(ctl, 2 * w + 2 + i, L, reg);
 
        }
          
        quantum_toffoli(ctl, 2 * w + 2, L, reg);
      
        emul(a % N, L, w, reg);
       
        quantum_toffoli(ctl, 2 * w + 2, L, reg);    
        
    }

    /**
     * runtime complexity is w * O(reg.size)
     * @param N
     * @param a
     * @param ctl
     * @param w
     * @param reg 
     */
   public void mul_mod_n(int N, int a, int ctl, int w, QuantumReg reg) {
        
        //runtime complexity is O(reg.size)
        muln(N, a, ctl, w, reg);

        //runtime complexity is w * O(reg.size)
        quantum_swaptheleads_omuln_controlled(ctl, w, reg);
        
        //runtime complexity is w * O(reg.size)
        muln_inv(N, a, ctl, w, reg);
    
    }
  
    // ------oaddn.c: Addition modulo an integer N ----

    private int num_regs = 4;
    
    /**
     * if bit "compare" - the global enable bit - is set, test_sums checks, if
     * the sum of the c-number and the q-number in register add_sum is greater
     * than n and sets the next lower bit to "compare"
     */
   public void test_sum(int compare, int w, QuantumReg reg) {
        int i;

        if ((compare & (1L << (w - 1))) != 0) {
            quantum_cnot(2 * w - 1, w - 1, reg);
            quantum_sigma_x(2 * w - 1, reg);
            quantum_cnot(2 * w - 1, 0, reg);
        } else {
            quantum_sigma_x(2 * w - 1, reg);
            quantum_cnot(2 * w - 1, w - 1, reg);
        }
        for (i = (w - 2); i > 0; i--) {
            if ((compare & (1 << i)) != 0) {
                //is bit i set in compare?
                quantum_toffoli(i + 1, w + i, i, reg);
                quantum_sigma_x(w + i, reg);
                quantum_toffoli(i + 1, w + i, 0, reg);
            } else {
                quantum_sigma_x(w + i, reg);
                quantum_toffoli(i + 1, w + i, i, reg);
            }
        }
        if ((compare & 1) != 0) {
            quantum_sigma_x(w, reg);
            quantum_toffoli(w, 1, 0, reg);
        }
        quantum_toffoli(2 * w + 1, 0, 2 * w, reg);//set output to 1 if enabled and b < compare

        if ((compare & 1) != 0) {
            quantum_toffoli(w, 1, 0, reg);
            quantum_sigma_x(w, reg);
        }

        for (i = 1; i <= (w - 2); i++) {
            if ((compare & (1 << i)) != 0) {
                //is bit i set in compare?
                quantum_toffoli(i + 1, w + i, 0, reg);
                quantum_sigma_x(w + i, reg);
                quantum_toffoli(i + 1, w + i, i, reg);
            } else {
                quantum_toffoli(i + 1, w + i, i, reg);
                quantum_sigma_x(w + i, reg);
            }
        }
        if ((compare & (1 << (w - 1))) != 0) {
            quantum_cnot(2 * w - 1, 0, reg);
            quantum_sigma_x(2 * w - 1, reg);
            quantum_cnot(2 * w - 1, w - 1, reg);
        } else {
            quantum_cnot(2 * w - 1, w - 1, reg);
            quantum_sigma_x(2 * w - 1, reg);
        }

    }

    /**This is a semi-quantum fulladder. It adds to b_in
    a c-number. Carry-in bit is c_in and carry_out is
    c_out. xlt-l and L are enablebits. See documentation
    for further information
   */
   public void muxfa(int a, int b_in, int c_in, int c_out, int xlt_l, int L, 
        int total, QuantumReg reg){

      //a,

        if (a == 0) {//00
            quantum_toffoli(b_in, c_in, c_out, reg);
            quantum_cnot(b_in, c_in, reg);
        }        
         
        if (a == 3) {//11
            quantum_toffoli(L, c_in, c_out, reg);
            quantum_cnot(L, c_in, reg);
            quantum_toffoli(b_in, c_in, c_out, reg);
            quantum_cnot(b_in, c_in, reg);          
        }
        
        if (a == 1) {//01
            quantum_toffoli(L, xlt_l, b_in, reg);
            quantum_toffoli(b_in, c_in, c_out, reg);
            quantum_toffoli(L, xlt_l, b_in, reg);
            quantum_toffoli(b_in, c_in, c_out, reg);
            quantum_toffoli(L, xlt_l, c_in, reg);
            quantum_toffoli(b_in, c_in, c_out, reg);
            quantum_cnot(b_in, c_in, reg);
        }
        
        if (a == 2) {//10
            quantum_sigma_x(xlt_l, reg);
            quantum_toffoli(L, xlt_l, b_in, reg);
            quantum_toffoli(b_in, c_in, c_out, reg);
            quantum_toffoli(L, xlt_l, b_in, reg);
            quantum_toffoli(b_in, c_in, c_out, reg);
            quantum_toffoli(L, xlt_l, c_in, reg);
            quantum_toffoli(b_in, c_in, c_out, reg);
            quantum_cnot(b_in, c_in, reg);
            quantum_sigma_x(xlt_l, reg);
        }
    }

    /**
    * This is just the inverse operation of the semi-quantum fulladder
    */
   public void muxfa_inv(int a, int b_in, int c_in, int c_out, int xlt_l, int L, 
        int total, QuantumReg reg){

        //a,

        if (a == 0) {//00
            quantum_cnot(b_in, c_in, reg);
            quantum_toffoli(b_in, c_in, c_out, reg);
        }

        if (a == 3) {//11
            quantum_cnot(b_in, c_in, reg);
            quantum_toffoli(b_in, c_in, c_out, reg);
            quantum_cnot(L, c_in, reg);
            quantum_toffoli(L, c_in, c_out, reg);
        }

        if (a == 1) {//01
            quantum_cnot(b_in, c_in, reg);
            quantum_toffoli(b_in, c_in, c_out, reg);
            quantum_toffoli(L, xlt_l, c_in, reg);
            quantum_toffoli(b_in, c_in, c_out, reg);
            quantum_toffoli(L, xlt_l, b_in, reg);
            quantum_toffoli(b_in, c_in, c_out, reg);
            quantum_toffoli(L, xlt_l, b_in, reg);
        }

        if (a == 2) {//10
            quantum_sigma_x(xlt_l, reg);
            quantum_cnot(b_in, c_in, reg);
            quantum_toffoli(b_in, c_in, c_out, reg);
            quantum_toffoli(L, xlt_l, c_in, reg);
            quantum_toffoli(b_in, c_in, c_out, reg);
            quantum_toffoli(L, xlt_l, b_in, reg);
            quantum_toffoli(b_in, c_in, c_out, reg);
            quantum_toffoli(L, xlt_l, b_in, reg);
            quantum_sigma_x(xlt_l, reg);
        }
    }

    /**This is a semi-quantum halfadder. It adds to b_in
    a c-number. Carry-in bit is c_in and carry_out is
    not necessary. xlt-l and L are enablebits. See
    documentation for further information
    */
   public void muxha(int a, int b_in, int c_in, int xlt_l, int L, int total, QuantumReg reg) {

        //a,

        if (a == 0) {//00
            quantum_cnot(b_in, c_in, reg);
        }

        if (a == 3) {//11
            quantum_cnot(L, c_in, reg);
            quantum_cnot(b_in, c_in, reg);
        }

        if (a == 1) {//01
            quantum_toffoli(L, xlt_l, c_in, reg);
            quantum_cnot(b_in, c_in, reg);
        }

        if (a == 2) {//10
            quantum_sigma_x(xlt_l, reg);
            quantum_toffoli(L, xlt_l, c_in, reg);
            quantum_cnot(b_in, c_in, reg);
            quantum_sigma_x(xlt_l, reg);
        }
    }

    //just the inverse of the semi quantum-halfadder
   public void muxha_inv(int a, int b_in, int c_in, int xlt_l, int L, int total, 
        QuantumReg reg){

        //a,

        if (a == 0) {//00
            quantum_cnot(b_in, c_in, reg);
        }

        if (a == 3) {//11
            quantum_cnot(b_in, c_in, reg);
            quantum_cnot(L, c_in, reg);
        }

        if (a == 1) {//01
            quantum_cnot(b_in, c_in, reg);
            quantum_toffoli(L, xlt_l, c_in, reg);
        }

        if (a == 2) {//10
            quantum_sigma_x(xlt_l, reg);
            quantum_cnot(b_in, c_in, reg);
            quantum_toffoli(L, xlt_l, c_in, reg);
            quantum_sigma_x(xlt_l, reg);
        }
    }

   public void madd(int a, int a_inv, int w, QuantumReg reg){
    
	    int i, j;
        int total;
        total = num_regs * w + 2;

        for (i = 0; i < w - 1; i++) {
            if (((1 << i) & a) != 0) {
                j = 2;
            } else {
                j = 0;
            }
            if (((1 << i) & a_inv) != 0) {
                j++;
            }
            muxfa(j, w + i, i, i + 1, 2 * w, 2 * w + 1, total, reg);
        }
        j = 0;
        if (((1 << (w - 1)) & a)  != 0) {
            j = 2;
        }
        if (((1 << (w - 1)) & a_inv) != 0) {
            j++;
        }
        muxha(j, 2 * w - 1, w - 1, 2 * w, 2 * w + 1, total, reg);
    }

   public void madd_inv(int a, int a_inv, int w, QuantumReg reg){
        
	    int i, j;
        int total;
        total = num_regs * w + 2;
        j = 0;

        if (((1<<(w-1)) & a) != 0) {
            j = 2;
        }
        if (((1<<(w-1)) & a_inv) != 0) {
            j++;
        }
        muxha_inv(j, w - 1, 2 * w - 1, 2 * w, 2 * w + 1, total, reg);

        for (i = w - 2; i >= 0; i--) {
            if (((1<<i) & a) != 0) {
                j = 2;
            } else {
                j = 0;
            }
            if (((1<<i) & a_inv) != 0) {
                j++;
            }
            muxfa_inv(j, i, w + i, w + 1 + i, 2 * w, 2 * w + 1, total, reg);
        }
    }

   public void addn(int N, int a, int w, QuantumReg reg){

        //add a to register reg (mod N)

	    test_sum(N - a, w, reg); //xlt N-a      
        
        madd((1 << w) + a - N, a, w, reg);//madd 2^K+a-N
    }

   public void addn_inv(int N, int a, int w, QuantumReg reg){

        //inverse of add a to register reg (mod N)

        quantum_cnot(2 * w + 1, 2 * w, reg);//Attention! cnot gate instead of not, as in description
        madd_inv((1 << w) - a, N - a, w, reg);//madd 2^K+(N-a)-N = 2^K-a

        quantum_swaptheleads(w, reg);

        test_sum(a, w, reg);
    }

   public void add_mod_n(int N, int a, int w, QuantumReg reg){
        
        //add a to register reg (mod N) and clear the scratch bits

	    addn(N, a, w, reg);
 
        addn_inv(N, a, w, reg);
    }

    // ------- quantum fourier transform ------    
    /**
     * Perform a QFT on a quantum register. This is done by application of
     * conditional phase shifts and hadamard gates. At the end, the position of
     * the bits is reversed.
     * 
     * runtime complexity is log2(N) * (.lt. log2(N)) * 2^(log2(N))
     */
   public void quantum_qft(int w, QuantumReg reg) {
        
        int i, j;

        System.out.format(" quantum_qft nloop=%d\n", w);
        
        //log2(N) * (< log2(N)) * 2^(log2(N))
        
        //2^(log2(N))
        for (i = w - 1; i >= 0; i--) {
            for (j = w - 1; j > i; j--) {
                quantum_cond_phase(j, i, reg);
            }

            System.out.format(" ..%d", i);

            //O(2^(log2(N)))
            quantum_hadamard(i, reg);
        }
        System.out.format("\n");        
    }

   public void quantum_qft_inv(int w, QuantumReg reg) {
  
        int i, j;

        for (i = 0; i < w; i++) {
            quantum_hadamard(i, reg);
            for (j = i + 1; j < w; j++) {
                quantum_cond_phase_inv(j, i, reg);
            }
        }
    }
    
   /**
     * Apply an irreversible AND gate.
     * 
     *                 Truth table
     * 
                 control1 control2 target
                    0        1       0
                    1        1       1
                    0        0       0
                    1        0       0
     
     NOTE: added by Nichole. Presumably the universal AND gate should be
     present on the computer.
    
     * @param control1
     * @param control2
     * @param target
     * @param reg
   */
   public void quantum_ccand(int control1, int control2, int target, 
       QuantumReg reg) {
  
        int i;
        int[] qec = new int[1];
       
        for (i = 0; i < reg.size; i++) {
            // set the target bit of a basis state to 1
            // if both control bits are set, else 0

            long st = reg.node[i].state;
            if ((st & (1L << control1)) != 0) {
                if ((st & (1L << control2)) != 0) {
                    reg.node[i].state |= (1L << target);
                } else {
                    reg.node[i].state &= ~(1L << target);
                }
            } else {
                reg.node[i].state &= ~(1L << target);
            }
        }
        
        decoherence.quantum_decohere(reg, rng, this);
    }

}