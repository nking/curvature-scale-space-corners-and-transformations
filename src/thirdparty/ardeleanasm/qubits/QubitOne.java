package thirdparty.ardeleanasm.qubits;

import thirdparty.ardeleanasm.complexnumbers.ComplexNumber;

/**
* <h1>QubitOne</h1>
* Representation of the qubit |1>=[0,1].
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
public class QubitOne extends Qubit {

	/**
	 * Construct a new <code> QubitOne</code> object.
	 */
	public QubitOne() {
        // the real values are 0, 1
		super(new ComplexNumber(0.0, 0.0), new ComplexNumber(1.0, 0.0));
	}

}
