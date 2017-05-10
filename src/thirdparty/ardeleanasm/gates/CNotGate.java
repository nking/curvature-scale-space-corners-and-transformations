package thirdparty.ardeleanasm.gates;

import thirdparty.ardeleanasm.complexnumbers.ComplexNumber;
import thirdparty.ardeleanasm.qubits.Qubit;

/**
 * Implements the CNOT Gate

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
public final class CNotGate implements IGate {

	@Override
	public Qubit applyGate(Qubit inputQubit, int[] targetPosition,
			int[] conditions, int noOfEntangledQubits) {

		int mask = 0;
		int newPosition = 0;
		ComplexNumber[] states = inputQubit.getQubit();
		for (int i : conditions) {
			mask |= (1 << (noOfEntangledQubits - 1 - i));
		}
		mask |= (1 << (noOfEntangledQubits - targetPosition[0]));
		newPosition = mask | 0x01;
		ComplexNumber state = states[mask];
		states[mask] = states[newPosition];
		states[newPosition] = state;

		return new Qubit(states);
	}

}
