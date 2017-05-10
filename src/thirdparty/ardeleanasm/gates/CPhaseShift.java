package thirdparty.ardeleanasm.gates;

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
 */
public final class CPhaseShift implements IGate {

	@Override
	public Qubit applyGate(Qubit inputQubit, int[] targetPosition,
			int[] conditions, int noOfEntangledQubits) {
		int mask = 0;
		int new_position = 0;

		ComplexNumber[] states = inputQubit.getQubit();
		for (int i : conditions) {
			mask |= (1 << (noOfEntangledQubits - 1 - i));
		}
		mask |= (1 << (noOfEntangledQubits - targetPosition[0] - 1));
		new_position = mask | 0x01;
		
        states[mask] = ComplexMath.multiply(states[mask], new ComplexNumber(
				-1.0, 0.0));

		states[new_position] = ComplexMath.multiply(states[new_position],
				new ComplexNumber(-1.0, 0.0));

		return new Qubit(states);
	}

}
