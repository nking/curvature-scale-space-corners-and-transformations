package thirdparty.ardeleanasm.qubits;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import thirdparty.ardeleanasm.quantum.exception.RegisterOverflowException;
import thirdparty.ardeleanasm.quantum.utils.QRegisterOperations;

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
 * a register of qubits where register is emulating 
 * a processor register is a quickly accessible location available 
 * to a computer's central processing unit (CPU). 
 registers are used for 
* arithmetic operations and are manipulated or tested by machine 
* instructions. Manipulated data is then often stored back to main memory, 
* either by the same instruction or by a subsequent one. Modern processors use 
* either static or dynamic RAM as main memory, 
* with the latter usually accessed via one or more cache levels.
* (https://en.wikipedia.org/wiki/Processor_register)
 */
public class QRegister implements Iterable<Qubit> {
	private List<Qubit>	qubitRegister;
	private int			registerSize;

	public QRegister(int size) {
		this.registerSize = size;
		qubitRegister = new ArrayList<>(size);

	}

	public QRegister initialize() {
		qubitRegister = QRegisterOperations.getInstance().fillWith(qubitRegister, 
            QubitZero::new, registerSize);
		return this;
	}

	public int size() {
		return registerSize;
	}

	public void add(Qubit e) throws RegisterOverflowException {
		if (qubitRegister.size() >= registerSize) {
			throw new RegisterOverflowException();
		}
		qubitRegister.add(e);
	}

	public QRegister change(int index, Qubit e) {
		qubitRegister.set(index, e);
		return this;
	}

	public Qubit get(int index) {
		return qubitRegister.get(index);
	}

	@Override
	public Iterator<Qubit> iterator() {
		return qubitRegister.iterator();
	}

	public List<Qubit> getQubits() {
		return qubitRegister;
	}

	public void setQubits(List<Qubit> qubits) {
		qubitRegister = qubits;
	}

	@Override
	public String toString() {
		StringBuffer buffer = new StringBuffer();
		for (Qubit q : qubitRegister)
			buffer.append(q);
		return buffer.toString();
	}
}
