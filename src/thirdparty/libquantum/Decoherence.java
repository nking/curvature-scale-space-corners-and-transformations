package thirdparty.libquantum;

import java.util.Random;

/* decoherence.c: Simulation of decoherence effects

ibquantum provides an efficient model to simulate the effects of 
decoherence. The effects are simulated by a random rotation about 
the z axis, where the angle of the rotation is a normal distributed 
value with the variance $2 \lambda$. $\lambda$ is the decoherence 
parameter, which depends on the experimental realization of the 
quantum computer. 
A list of values for $\lambda$ is given in [DiVincenzo, 1995].

   Copyright 2003 Bjoern Butscher, Hendrik Weimer

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
public class Decoherence {

    /**
     * Status of the decoherence simulation. Non-zero means enabled and
     * decoherence effects will be simulated.
     */
    private int quantum_status = 0;

    /**
     * Decoherence parameter. The higher the value, the greater the decoherence
     * impact.
     */
    private float quantum_lambda = 0;

    float quantum_get_decoherence() {
        return quantum_lambda;
    }

    /**
     * Initialize the decoherence simulation and set the decoherence parameter.
     */
    void quantum_set_decoherence(float l) {
        if (l != 0) {
            quantum_status = 1;
            quantum_lambda = l;
        } else {
            quantum_status = 0;
        }
    }

    /**
     * Perform the actual decoherence of a quantum register for a single step of
     * time. This is done by applying a phase shift by a normal distributed
     * angle with the variance LAMBDA.
     */
    void quantum_decohere(QuantumReg reg, Random rng, Gates gates) {
        double u, v, s, x;
        double[] nrands;
        double angle;
        int i, j;

        /**
         * Increase the gate counter
         */
        gates.quantum_gate_counter(1);
        
        if (quantum_status != 0) {

            nrands = new double[reg.width];

            for (i = 0; i < reg.width; i++) {
                
                // Generate normal distributed random numbers
                do {
                    u = 2 * rng.nextDouble() - 1;
                    v = 2 * rng.nextDouble() - 1;
                    s = u * u + v * v;
                } while (s >= 1);

                x = u * Math.sqrt(-2 * Math.log(s) / s);

                x *= Math.sqrt(2 * quantum_lambda);

                nrands[i] = x / 2;
            }

            // Apply the phase shifts for decoherence simulation */
            for (i = 0; i < reg.size; i++) {
                
                angle = 0;

                for (j = 0; j < reg.width; j++) {
                                        
                    if ((reg.node[i].state & (1L << j)) != 0) {
                        angle += nrands[j];
                    } else {
                        angle -= nrands[j];
                    }
                }
                //cos(phi) + IMAGINARY * sin(phi)
                double re = reg.node[i].amplitude.re() * Math.cos(angle);
                double im = reg.node[i].amplitude.im() * Math.sin(angle);
                reg.node[i].amplitude.setReal(re);
                reg.node[i].amplitude.setImag(im);
            }
        }
    }
}
