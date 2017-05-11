package thirdparty.libquantum;

import algorithms.misc.ComplexModifiable;
import java.util.logging.Logger;

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

public class QuantumRegNode {
    
    private static transient Logger log = Logger.getLogger(
        QuantumRegNode.class.getName());

    //COMPLEX_FLOAT amplitude; /* alpha_j */
    //MAX_UNSIGNED state;      /* j */

    public ComplexModifiable amplitude;
    private int state = 0;
    
    public void setState(int s) {
        if (s < 0) {
            log.warning("s < 0 so truncating to 0");
            s = 0;
        }
        this.state = s;
    }
    
    public int getState() {
        return state;
    }
}
