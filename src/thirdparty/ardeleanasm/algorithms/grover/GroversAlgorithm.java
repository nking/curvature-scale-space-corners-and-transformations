package thirdparty.ardeleanasm.algorithms.grover;

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
public class GroversAlgorithm extends QuantumAlgorithms {
	private static final int NO_OF_INPUT = 3;
	private IGate gateH;
	private IGate gateX;
	private IGate gateCNot;
	private IGate gateCPhaseShift;
	

	public GroversAlgorithm() {


	}


	@Override
	public void init() {
		gateH = gateFactory.getGate(EGateTypes.E_HadamardGate);
		gateX = gateFactory.getGate(EGateTypes.E_XGate);
		gateCNot = gateFactory.getGate(EGateTypes.E_CNotGate);
		gateCPhaseShift = gateFactory.getGate(EGateTypes.E_CPhaseShift);
		QRegister qRegister=new QRegister(NO_OF_INPUT+1);
		QRegisterOperations qRegOps=QRegisterOperations.getInstance();
		try {
			qRegister=qRegOps.fillWithPattern("0001");
		} catch (RegisterOverflowException e) {
			e.printStackTrace();
		}
		resultQubit=QuantumOperations.entangle(qRegister);
	}

	@Override
	public void run() {
		resultQubit=gateH.applyGate(resultQubit, new int[]{0,1,2,3}, null, NO_OF_INPUT+1);
		resultQubit=gateCNot.applyGate(resultQubit, new int[]{3}, new int[]{0,1,2}, NO_OF_INPUT+1);
		
		
		resultQubit=gateH.applyGate(resultQubit, new int[]{0,1,2}, null, NO_OF_INPUT+1);
		resultQubit=gateX.applyGate(resultQubit, new int[]{0,1,2}, null, NO_OF_INPUT+1);
		resultQubit=gateCPhaseShift.applyGate(resultQubit, new int[]{2}, new int[]{0,1}, NO_OF_INPUT+1);
		
		resultQubit=gateX.applyGate(resultQubit, new int[]{0,1,2}, null, NO_OF_INPUT+1);
		resultQubit=gateH.applyGate(resultQubit, new int[]{0,1,2,3}, null, NO_OF_INPUT+1);
		
		
		
		
	}

	@Override
	public void measure() {
		MeasurementPerformer measurementPerformer = new MeasurementPerformer().configure(resultQubit);
		resultQubit = measurementPerformer.measure();
		System.out.println(resultQubit);

	}



}
