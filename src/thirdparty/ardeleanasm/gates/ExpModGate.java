package thirdparty.ardeleanasm.gates;

import algorithms.misc.ComplexModifiable;
import java.util.Arrays;

import thirdparty.ardeleanasm.complexnumbers.ComplexMath;
import thirdparty.ardeleanasm.complexnumbers.ComplexNumber;
import thirdparty.ardeleanasm.qubits.Qubit;

/**
 * applies outer product to each state of inputQubit
 * the value of a^x mod N
 * where a and N are given in conditions.
 * 
 * The operation is the summation over x (|x>  |f(x)>).
 *   where f(x) = a^x mod N
 *   
 * 
 * Note that there are alternate ways to implement this operation and they 
 * may be better numerically or more efficient.
 * 
 * for example, for this quantum operation, 
 * libquantum uses decohere, squaring, modulus, toffoli gate, and more
                         
 * 
 * @author nichole
*/
public final class ExpModGate implements IGate {

	@Override
	public Qubit applyGate(Qubit inputQubit, int[] targetPosition,
	    int[] conditions, int noOfEntangledQubits) {
		
        //TODO: could consider refactoring to use ComplexModifiable
        
        if (conditions == null || conditions.length < 2) {
            throw new IllegalArgumentException("conditions must be length=2 "
                + " and must have a and N in it for the operation a^x mod N");
        }
        
        System.out.println("noOfEntangledQubits=" + noOfEntangledQubits 
            + " targetPosition=" + 
            Arrays.toString(targetPosition)
        );
        
        int a = conditions[0];
        int N = conditions[1];
        				
        ComplexNumber[] states = inputQubit.getQubit();
		ComplexNumber bufferState;
		            
        for (int j = 0; j < states.length; j++) {

            bufferState = states[j];
                        
            double x = ComplexMath.mod(bufferState);
            
            //a^x mod N
            
            double y = Math.pow(a, x);
            y = y % N;
            
            bufferState = ComplexMath.multiply(bufferState, y);
            
            states[j] = bufferState;
		}
        
        Qubit output = new Qubit(states);
        assert(output.isValid());
		
        return output;
	}
    
}
