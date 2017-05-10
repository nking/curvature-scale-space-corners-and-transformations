package thirdparty.ardeleanasm.gates;

import java.util.Arrays;

import thirdparty.ardeleanasm.complexnumbers.ComplexMath;
import thirdparty.ardeleanasm.complexnumbers.ComplexNumber;
import thirdparty.ardeleanasm.qubits.Qubit;

/*
code is from
 https://github.com/ardeleanasm/quantum_computing.git
 which has license
 ardeleanasm/quantum_computing is licensed under the
 GNU General Public License v3.0
        
 Permissions of this strong copyleft license are conditioned
 on making available complete source code of licensed works
 and modifications, which include larger works using a
 licensed work, under the same license. Copyright and
 license notices must be preserved. Contributors provide
 an express grant of patent rights.
*/

/**
 * 
 * Implements the Hadamard Gate
*/
public final class HGate implements IGate {

	@Override
	public Qubit applyGate(Qubit inputQubit, int[] targetPosition,
	    int[] conditions, int noOfEntangledQubits) {
		
        //TODO: could consider refactoring to use ComplexModifiable
        
        System.out.println("noOfEntangledQubits=" + noOfEntangledQubits 
            + " targetPosition=" + 
            Arrays.toString(targetPosition)
        );
        
        int mask = 0;
		int newPosition = 0;
		
        double hCoefficient = 1.0 / Math.sqrt(2.0);
		
        ComplexNumber[] states = inputQubit.getQubit();
		ComplexNumber bufferState;
		
        int[] markedStates = new int[states.length];
                
		for (int i : targetPosition) {
			
            Arrays.fill(markedStates, 0);

            int xpon = (noOfEntangledQubits - 1 - i);
            assert(xpon < 32);
            
			mask = (1 << xpon);
            
            System.out.println("left shift=" + xpon 
                + " states.len=" + states.length);
                        
			for (int j = 0; j < states.length; j++) {
			
                if (markedStates[j] == 0) {
					
                    newPosition = j ^ mask;
                    
					bufferState = states[j];
                    
                    // mask = 2^((nQubits - 1 - i)
                    // H(|j>) = hCoeff * (states[j] + states[j^mask])
					states[j] = ComplexMath.multiply(
							ComplexMath.add(bufferState, states[newPosition]),
							hCoefficient);

                    // H(|j^mask>) = hCoeff * (states[j] - states[j^mask])
					states[newPosition] = ComplexMath.multiply(ComplexMath
							.subtract(bufferState, states[newPosition]),
							hCoefficient);
                    
					markedStates[j] = 1;
                    
					markedStates[newPosition] = 1;
				}
			}
			mask = 0;
		}
        
        Qubit output = new Qubit(states);
        assert(output.isValid());
		
        return output;
	}
    
}
