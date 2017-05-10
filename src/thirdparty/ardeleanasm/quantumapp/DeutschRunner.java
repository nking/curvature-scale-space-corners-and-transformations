package thirdparty.ardeleanasm.quantumapp;

import thirdparty.ardeleanasm.algorithms.QuantumAlgorithms;
import thirdparty.ardeleanasm.algorithms.deutsch.DeutschsAlgorithm;

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
public class DeutschRunner implements Runnable {
	private static final QuantumAlgorithms	ALGORITHM		= new DeutschsAlgorithm();
	

	@Override
	public void run() {
		ALGORITHM.init();
		ALGORITHM.run();
		ALGORITHM.measure();
		
	}
}
