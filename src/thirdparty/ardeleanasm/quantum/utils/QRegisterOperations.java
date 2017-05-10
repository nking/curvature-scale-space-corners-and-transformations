package thirdparty.ardeleanasm.quantum.utils;

import java.util.List;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import thirdparty.ardeleanasm.quantum.exception.RegisterOverflowException;
import thirdparty.ardeleanasm.qubits.QRegister;
import thirdparty.ardeleanasm.qubits.Qubit;
import thirdparty.ardeleanasm.qubits.QubitOne;
import thirdparty.ardeleanasm.qubits.QubitZero;

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
public class QRegisterOperations {
	private static final QRegisterOperations INSTANCE = new QRegisterOperations();

	private QRegisterOperations() {
	}

	public static QRegisterOperations getInstance() {
		return INSTANCE;
	}

	public List<Qubit> fillWith(List<Qubit> list,
			Supplier<Qubit> qubitSupplier, int size) {
		list = Stream.generate(qubitSupplier).limit(size)
				.collect(Collectors.toList());
		return list;
	}

	public QRegister fillWith(QRegister reg, Supplier<Qubit> qubitSupplier) {
		reg.setQubits(Stream.generate(qubitSupplier).limit(reg.size())
				.collect(Collectors.toList()));
		return reg;
	}

	public QRegister fillWithPattern(String pattern)
			throws RegisterOverflowException {
        
		QRegister qreg = new QRegister(pattern.length()).initialize();
		for (int i = 0; i < pattern.length(); i++) {
			if (pattern.charAt(i) == '1') {
				qreg.change(i, new QubitOne());
			} //else {
                // already initialized to '0'
				//qreg.change(i, new QubitZero());
			//}
		}
		return qreg;
	}

	/**
	 * Perform the tensor product between two or more qubits. Example, for three
	 * qubits |0>, |0> and |1>, the result will be |001>.
	 * 
	 * @param quantumRegister
	 * @return qubit the tensor product of the two qubits.
	 */
	public Qubit entangle(QRegister quantumRegister) {
		if (quantumRegister.size() < 2) {
			return null;
		}
		Qubit bufferQubit = quantumRegister.get(0);
		for (int i = 1; i < quantumRegister.size(); i++) {
			bufferQubit = QuantumOperations.entangle(bufferQubit,
					quantumRegister.get(i));

		}
		return bufferQubit;
	}

}
