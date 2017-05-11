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

public class QuantumDensityOp {

    /** total number of state vectors */
    public int num;   
    
    /** probabilities of the state vectors */
    public float[] prob;      
   
    /** state vectors */
    public QuantumReg reg;
  
    /* for OBJCODE?
    #define quantum_density_operation(function, rho, ...) \
    do{ \
      int quantum_int; \
      for(quantum_int=0; quantum_int < rho.num; quantum_int++) \
        function(__VA_ARGS__, &rho.reg[quantum_int]); \
    } while(0)
    */

}
