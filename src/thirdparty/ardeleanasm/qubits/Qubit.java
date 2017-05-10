package thirdparty.ardeleanasm.qubits;

/**
 * <h1>QuantumComputingLib</h1>
 * QuantumComputingLib is a Java library for performing Quantum Computations.
 * This library provides methods for performing operations on qubits ( apply gates,
 * tensor products) and for performing basic operations with 2D-arrays ( matrices). 
 * @author Mihai Seba
 * @version 1.0-SNAPSHOT
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
import java.util.Arrays;

import thirdparty.ardeleanasm.complexnumbers.ComplexMath;
import thirdparty.ardeleanasm.complexnumbers.ComplexNumber;

public class Qubit {

	protected ComplexNumber[] qubitVector;
	
	/**
	 * Constructs a new <code>Qubit</code> object. 
	 * @param no0 complex number
	 * @param no1 complex number
	 * 
	 */
	public Qubit(ComplexNumber no0, ComplexNumber no1) {
		qubitVector = new ComplexNumber[]{no0, no1};
	}

	/**
	 * Constructs a new <code>Qubit</code> object.
	 * @param qubitVector an array of 2 complex numbers
	 */
	public Qubit(ComplexNumber[] qubitVector) {
		this.qubitVector = Arrays.copyOf(qubitVector, qubitVector.length);
	}

	/**
	 * Return the qubit represented as an array of 2 complex numbers.  
	 * @return qubit
	 */
	public ComplexNumber[] getQubit() {
		return qubitVector;
	}

	/**
	 * Return a string representation of the qubit.
	 * @return string the representation of the qubit
	 */
	@Override
	public String toString() {
		StringBuffer output=new StringBuffer();
		output.append("[ ");
		for (ComplexNumber i : qubitVector) {
			output.append(i);
			output.append(" ");
		}
		output.append("]");
		return output.toString();
	}

	/**
	 * Check if passed <code>Qubit</code> is equal to the current.
	 * @param o the qubit to be checked
	 * @return true if the two qubits are equals, otherwise false
	 */
	@Override
	public boolean equals(Object o) {
		if (o instanceof Qubit) {
			if (this.qubitVector.length != ((Qubit) o).getQubit().length) {
				return false;
			}
			for (int i = 0; i < this.qubitVector.length; i++) {
				if (this.qubitVector[i].equals(((Qubit) o).getQubit()[i])==false) {
					return false;
				}
			}
			return true;
		}
		return false;
	}

	/**
	 * Calculate the hashcode of the qubit.
	 * @return int hashcode
	 */
	@Override
	public int hashCode() {
		int hash = 5;
		hash += (this.qubitVector != null ? Arrays.hashCode(qubitVector) : 0);
		return hash;
	}
	
	/**
	 * Check if qubit state is valid, that is, the sum of the squares == 1,
     * that is, L2 norm is 1
	 * @return true if the state is valid, otherwise false
	 */
	public boolean isValid(){
        double eps = 1.e-6;
		double sum=0.0;
		for(ComplexNumber c:this.qubitVector){
			double mod = ComplexMath.mod(c);
			sum += mod*mod;
		}
		return (Math.abs(sum - 1.0) < eps);
	}
}
