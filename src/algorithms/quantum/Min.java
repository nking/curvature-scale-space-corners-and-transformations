package algorithms.quantum;

import algorithms.misc.Misc;
import java.util.Random;
import thirdparty.libquantum.Grover;
import thirdparty.libquantum.Measure;
import thirdparty.libquantum.QuReg;
import thirdparty.libquantum.QuantumReg;

/**
 * NOT READY FOR USE.
 * 
 * implementation of
 * "A quantum algorithm for finding the minimum"
 * by Durr and Hoyer
 * 
 * @author nichole
 */
public class Min {
    
    /**
     NOT READY FOR USE.
     NOTE: editing is in progress regarding the input list handling.
     
      find the index y such that list[y] is minimum.
       
      runtime complexity is roughly
          (22.5 * Math.sqrt(list.length) + 
          1.4 * Math.log(Math.log(list.length)))
          times O(grover + list_filtering)
           
      where O(Grover) is O((2^width) * width).
     
     * NOTE that the largest number in the list must be
     * .lte. integer.max_value - 2^width.
     * 
     * Also note that the input list may need to be composed of valid input.
     * still thinking about valid input to Grover's.
     * The register, as the possible states of superposition of
     * qubits, will have all possible permutation up to a power of 2 of
     * numbers in the list.
     * A straight sequence of numbers from 0 up to a power of 2 is valid input for
     * the current logic (can be unordered).
     * 
     * @param width largest bit length to use the register.
     * It should be the bit length of the largest expected number then add 1.
     * NOTE that width must be .gte. 2.
     * @param list a list of unordered numbers to search for number within
     * @return 
     */
    public int run(int width, int[] list) {

        int i;
        
        if (width < 2) {
            width = 2;
        }

        System.out.format("list.length=%d, width=%d\n", list.length, width);
        
        double limit = 22.5 * Math.sqrt(list.length) + 
            1.4 * Math.log(Math.log(list.length));
        
        /*
        find the index y such that T[y] is minimum.
        
        let N = T.length;
        
        (1) Choose threshold index 0 <= y <= N−1 uniformly at random
        (2) Repeat the following and interrupt it when the
            total running time is more than the limit.
            Mark every item T[j] < T[y]
            Then go to stage 2(2c)
          (a) initialize the memory as
              summation over j of (1/sqrt(N)) * |j>|y>
              where j are the indexes where T is < T[y]
          (b) apply Grover's quantum algorithm
          (c) observe the first register:
                let y' be the outcome. If T[y'] < T[y], then set threshold
                index y to y'
        (3) return y
        
        Note that each stage of initializing the register requires a valid
        input list.
        a valid list is most easily 0 to a power of 2,
        but doing so requires reading the list again... 
        still thinking about the the input list and filtering and that there
        will always be a state 0...
        
        the paper suggestion of a marked list doesn't exactly make sense to me
        yet... would need to consider the eigenstates that contribute to the
        state that one would like to filter out and those individually affect
        other states...cannot reduce the amplitude of a specific state
        without affecting other states representing superposition of some
        of the same qubits...
        */

        QuReg qureg = new QuReg();

        Grover grover = new Grover();
        
        Random rng = Misc.getSecureRandom();
        
        Measure measure = new Measure();
        
        QuantumReg reg;
        
        //a reverse hash to lookup index when given value.
        int[] hash = createHash(list);
        
        int y = rng.nextInt(list.length);
        
        int[] filtered = filter(list, y);
        
        int srch, y2;
        int[] hash2;
        
        for (int nIter = 0; nIter < limit; ++nIter) {
           
            System.out.println("current minIdx=" + y);
            
            if (filtered.length == 0) {
                break;
            }
            
            //TODO: editing the handling of the input list to the register
            //
            
            reg = grover.initializeRegister(qureg, filtered, width);
        
            /*
            choose srch uniformly at random among the nonnegative
                integers smaller than list[y].
            */
            srch = filtered[rng.nextInt(filtered.length)];
            
            grover.processInitialized(srch, reg, rng);
       
            //DEBUG
            //System.out.format("AFTER process  reg.size=%d\n", reg.size);
            //qureg.quantum_print_qureg(reg);

            measure.quantum_bmeasure(0, reg, rng);
            
            //DEBUG
            //System.out.format("AFTER measure  reg.size=%d\n", reg.size);
            //qureg.quantum_print_qureg(reg);

            int v = (int)reg.node[0].state;
        
            if (v < list[y]) {
                y = getIndex(list, v, hash);
            }
            //TODO: check this.
            if (srch < list[y]) {
                y = getIndex(list, srch, hash);
            }
            
            hash2 = createHash(filtered);
            y2 = getIndex(filtered, list[y], hash2);
            filtered = filter(filtered, y2);
                   
        }
        
        return y;
    }
    
    /**
     * get index of list from value using the hash.
     * 
     * This method is based upon a method from the libquantum library 
     * file qureg.c which has copyright:
     
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
     * 
     * 
     * @param list
     * @param value
     * @param hash 
     * @return index of list where value is stored
     */
    private int getIndex(final int[] list, final int value, int[] hash) {
        
        int i, mark = 0;

        i = hashToIndex(value, list.length);
        
        while (hash[i] != -1) {
            if (list[hash[i] - 1] == value) {
                return hash[i] - 1;
            }
            i++;
            if (i == list.length) {
                i = 0;
            }
        }

        return -1;
    }
    
    /**
     * set list[idx] into hash.
     * 
     * This method is based upon a method from the libquantum library 
     * file qureg.c which has copyright:
     
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
     * 
     * @param list
     * @param idx
     * @param hash
     * @return index of hash where list[idx] was stored
     */
    private int setValue(final int[] list, final int idx, int[] hash) {
        
        int i, mark = 0;

        i = hashToIndex(list[idx], list.length);
        
        int end = list.length;
        while (hash[i] != -1) {
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
        
        hash[i] = idx + 1;     

        return i;
    }
    
    private int[] createHash(int[] list) {
        
        int[] hash = new int[list.length];
        for (int i = 0; i < list.length; ++i) {
            hash[i] = -1;
        }
        
        for (int i = 0; i < list.length; ++i) {
            int hIdx = setValue(list, i, hash);
        }

        return hash;        
    }

    private int[] filter(int[] list, int threshIdx) {
         
        int[] m = new int[list.length];
        
        int c = list[threshIdx];
        
        int i;
        int j = 0;
                
        for (i = 0; i < list.length; ++i) {
            if (list[i] < c) {
                m[i] = list[i];
                ++j;
            }
        }
        
        int[] out = new int[j];
        for (i = 0; i < j; ++i) {
            out[i] = m[i];
        }
        
        return out;
    }

    /**
     * return an index for the given value
     * 
     * @param value
     * @param length size of the container in which the result is an index in
     * @return 
     */
    protected int hashToIndex(int value, int length) {
       
        return fnvRetry(value, length);
    }
    
    protected int fnvRetry(int value, int length) {
        
        // http://www.isthe.com/chongo/tech/comp/fnv/index.html#lazy-mod
        
        //#define TRUE_HASH_SIZE ((u_int32_t)50000) /* range top plus 1
        int FNV_32_PRIME = 16777619;
        //int FNV1_32_INIT = 2166136261;
        int FNV1_32_INIT = 1083068130; // which is 2166136261 >> 1
        //#define MAX_32BIT ((u_int32_t)0xffffffff) /* largest 32 bit unsigned value
        //#define RETRY_LEVEL ((MAX_32BIT / TRUE_HASH_SIZE) * TRUE_HASH_SIZE)
        int RETRY_LEVEL = (Integer.MAX_VALUE/length) * length;
        
        int hash = fnv_31(value, FNV1_32_INIT);
        //System.out.println("hash=" + hash);
        while (hash >= RETRY_LEVEL) {
            hash = (hash * FNV_32_PRIME) + FNV1_32_INIT;
            //System.out.println("  hash=" + hash);
        }
        
        //NOTE: forcing positive for use case
        hash = (hash ^ (hash >> 31)) + (hash >>> 31);
        
        hash %= length;
        
        return hash;
    }

    private int fnv_31(int value, int FNV1_31_INIT) {

        //http://www.isthe.com/chongo/src/fnv/hash_32.c
       
        int hval = FNV1_31_INIT;
        
        //multiply by the 32 bit FNV magic prime mod 2^32
        //if no opt: hval *= 0x01000193;//FNV_32_PRIME;
        //else:
        hval += (hval<<1) + (hval<<4) + (hval<<7) 
            + (hval<<8) + (hval<<24);

	    // xor the bottom with the current octet
	    hval ^= value;
  
        return hval;
    }
}
