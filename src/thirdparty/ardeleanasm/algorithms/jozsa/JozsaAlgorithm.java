package thirdparty.ardeleanasm.algorithms.jozsa;

import thirdparty.ardeleanasm.algorithms.MeasurementPerformer;
import thirdparty.ardeleanasm.algorithms.QuantumAlgorithms;
import thirdparty.ardeleanasm.gates.EGateTypes;
import thirdparty.ardeleanasm.gates.IGate;
import thirdparty.ardeleanasm.quantum.exception.RegisterOverflowException;
import thirdparty.ardeleanasm.quantum.utils.QRegisterOperations;
import thirdparty.ardeleanasm.quantum.utils.QuantumOperations;
import thirdparty.ardeleanasm.qubits.QRegister;

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
public class JozsaAlgorithm extends QuantumAlgorithms {

	private static final int NO_OF_QUBITS = 3;
	private IGate gateHadamard;
	private IGate oracle;

	public JozsaAlgorithm(){
	}

	@Override
	public void init() {
		gateHadamard = gateFactory.getGate(EGateTypes.E_HadamardGate);
		oracle = gateFactory.getGate(EGateTypes.E_CPhaseShift);
		QRegister reg=new QRegister(NO_OF_QUBITS);
		try {
			reg=QRegisterOperations.getInstance().fillWithPattern("000");
		} catch (RegisterOverflowException e) {
			e.printStackTrace();
		}
		resultQubit=QuantumOperations.entangle(reg);

	}

	@Override
	public void measure() {
		MeasurementPerformer measureObject = new MeasurementPerformer();
		// Measure qubit |000>
		measureObject = measureObject.configure(resultQubit);
		resultQubit = measureObject.measure();
	}

	@Override
	public void run() {
                
		// First Step: Apply H Gate to |000>
		resultQubit = gateHadamard.applyGate(resultQubit, new int[]{0,1,2}, 
            null, NO_OF_QUBITS);
		// Second Step:Apply oracle
		resultQubit = oracle.applyGate(resultQubit, new int[]{2}, new int[]{0,1}, NO_OF_QUBITS);
		// Third Step: Apply H Gate
		resultQubit = gateHadamard.applyGate(resultQubit, new int[]{0,1,2}, null, NO_OF_QUBITS);
	}


}
