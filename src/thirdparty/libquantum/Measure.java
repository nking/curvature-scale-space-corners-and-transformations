package thirdparty.libquantum;

import algorithms.misc.ComplexModifiable;
import java.util.Random;

/* measure.c: Quantum register measurement

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
public class Measure {

    /**
     * Measure the contents of a quantum register
     */
    //MAX_UNSIGNED quantum_measure(QuantumReg reg, Random rng) {
    public long quantum_measure(QuantumReg reg, Random rng) {
        
        double r;
        int i;

        // Get a random number between 0 and 1 
        r = rng.nextDouble();

        for (i = 0; i < reg.size; i++) {
            // If the random number is less than the probability of the
            // given base state - r, return the base state as the
            // result. Otherwise, continue with the next base state. 

            r -= reg.node[i].amplitude.squareSum();
            if (0 >= r) {
                return reg.node[i].state;
            }
        }

        // The sum of all probabilities is less than 1. Usually, 
        // the cause for this is the application of a non-normalized 
        // matrix, but there is a slim chance that rounding errors 
        // may lead to this as well. 
        return -1;
    }

    /**
     * Measure a single bit of a quantum register. The bit measured is indicated
     * by its position POS, starting with 0 as the least significant bit. The
     * new state of the quantum register depends on the result of the
     * measurement.
     * runtime complexity is O(reg.size)
     */
    public int quantum_bmeasure(int pos, QuantumReg reg, Random rng) {
  
        int i;
        int result = 0;
        double pa = 0, r;
        //MAX_UNSIGNED pos2;

        long pos2 = 1L << pos;

        //System.out.format("REG.size=%d  pos2=%d\n", reg.size, pos2);
                
        // Sum up the probability for 0 being the result 
        for (i = 0; i < reg.size; i++) {
            
            long st = reg.node[i].state; 
            
            if ((st & pos2) == 0) {
                pa += reg.node[i].amplitude.squareSum();
            }
        }

        // Compare the probability for 0 with a random number 
        // and determine the result of the measurement 

        r = rng.nextDouble();
                
        if (r > pa) {
            // more states have a set bit in position pos
            result = 1;
        }
        
        QuReg qureg = new QuReg();

        // if result==1, removes all states without bit pos set
        //   and renormalizes the amplitudes,
        //   else, does the same for statues with bit pos set.
        // runtime complexity is O(reg.size)
        QuantumReg out = qureg.quantum_state_collapse(
            pos, result, reg);

        qureg.quantum_copy_qureg(out, reg);
    
        return result;
    }

    /**
     * Measure a single bit, but do not remove it from the quantum register
     */
    int quantum_bmeasure_bitpreserve(int pos, QuantumReg reg, Random rng) {
        int i, j;
        int size = 0, result = 0;
        double d = 0, pa = 0, r;
        //MAX_UNSIGNED pos2;

        long pos2 = 1L << pos;

        // Sum up the probability for 0 being the result 
        for (i = 0; i < reg.size; i++) {
            if ((reg.node[i].state & pos2) == 0) {
                pa += reg.node[i].amplitude.squareSum();
            }
        }

        // Compare the probability for 0 with a random number 
        // and determine the result of the measurement 
        r = rng.nextDouble();

        if (r > pa) {
            result = 1;
        }

        // Eradicate all amplitudes of base states which have been 
        // ruled out by the measurement and get the absolute 
        // of the new register 
        for (i = 0; i < reg.size; i++) {
            if ((reg.node[i].state & pos2) != 0) {
                if (result == 0) {
                    reg.node[i].amplitude.setReal(0);
                    reg.node[i].amplitude.setImag(0);
                } else {
                    d += reg.node[i].amplitude.squareSum();
                    size++;
                }
            } else {
                if (result != 0) {
                    reg.node[i].amplitude.setReal(0);
                    reg.node[i].amplitude.setImag(0);
                } else {
                    d += reg.node[i].amplitude.squareSum();
                    size++;
                }
            }
        }
        
        QuantumReg out = new QuantumReg();
        // Build the new quantum register 
        out.size = size;
        out.node = new QuantumRegNode[size];
        for (int ii = 0; ii < reg.node.length; ii++) {
            reg.node[ii] = new QuantumRegNode();
            reg.node[ii].state = 0;
            reg.node[ii].amplitude = new ComplexModifiable(0, 0);
        }
        out.hashw = reg.hashw;
        out.hash = reg.hash;
        out.width = reg.width;

        // Determine the numbers of the new base states and 
        // norm the quantum register 
        for (i = 0, j = 0; i < reg.size; i++) {
            if (reg.node[i].amplitude.abs() != 0) {
                out.node[j].state = reg.node[i].state;
                out.node[j].amplitude.resetTo(reg.node[i].amplitude);
                out.node[j].amplitude.times(1/Math.sqrt(d));

                j++;
            }
        }
        
        QuReg qureg = new QuReg();

        qureg.quantum_copy_qureg(out, reg);
        
        return result;
    }
}