package thirdparty.libquantum;

/* 
   Copyright 2003-2007 Bjoern Butscher, Hendrik Weimer

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
 The quantum register
 */
public class QuantumReg {
    
    /** number of qubits in the qureg */
    public int width;

    /** number of non-zero vectors */    
    public int size;
  
    /** width of the hash array */
    public int hashw;
    
    /** vectors w/ non-zero probability amplitudes */
    public QuantumRegNode[] node;
    
    /** hash nolding indexes of node.  the hash length
     is usually instantiated as (1 << hashw).*/
    public int[] hash;

}
