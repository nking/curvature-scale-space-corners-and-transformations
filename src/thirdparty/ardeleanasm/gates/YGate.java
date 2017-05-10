package thirdparty.ardeleanasm.gates;

import java.util.Arrays;

import thirdparty.ardeleanasm.complexnumbers.ComplexMath;
import thirdparty.ardeleanasm.complexnumbers.ComplexNumber;
import thirdparty.ardeleanasm.qubits.Qubit;

/**
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
public final class YGate implements IGate {

	@Override
	public Qubit applyGate(Qubit inputQubit, int[] targetPosition,
			int[] conditions, int noOfEntangledQubits) {
		int mask = 0;
		int newPosition = 0;
		ComplexNumber[] states = inputQubit.getQubit();
		ComplexNumber bufferState;
		int[] markedStates = new int[states.length];
		for (int i : targetPosition) {
			Arrays.fill(markedStates, 0);
			mask = (1 << (noOfEntangledQubits - 1 - i));
			for (int j = 0; j < states.length; j++) {
				if (markedStates[j] == 0) {
					newPosition = j ^ mask;
					bufferState = states[j];
					states[j] = ComplexMath.multiply(new ComplexNumber(0.0,
							-1.0), states[newPosition]);
					states[newPosition] = ComplexMath.multiply(
							new ComplexNumber(0.0, 1.0), bufferState);
					markedStates[j] = 1;
					markedStates[newPosition] = 1;
				}
				continue;
			}
			mask = 0;
		}
		return new Qubit(states);

	}
}
