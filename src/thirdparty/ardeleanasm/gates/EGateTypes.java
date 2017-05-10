package thirdparty.ardeleanasm.gates;

/**
 * Implemented Quantum Gates.

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
public enum EGateTypes {
	/**
	 * Hadamard Gate
	 */
	E_HadamardGate,
	/**
	 * Pauli-X Gate
	 */
	E_XGate,
	/**
	 * Pauli-Y Gate
	 */
	E_YGate,
	/**
	 * Pauli-Z Gate
	 */
	E_ZGate,
	/**
	 * CNOT Gate
	 */
	E_CNotGate,
	/**
	 * Controlled Phase Shift
	 */
	E_CPhaseShift,
	/**
	 * Identity gate
	 */
	E_IGate
	
}
