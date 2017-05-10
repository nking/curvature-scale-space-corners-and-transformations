package thirdparty.ardeleanasm.algorithms;

import java.util.Random;

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
public class MeasurementPerformer {
	private Qubit resultQubit = null;
	private int length;
	private static final double Q_TOLERANCE = 0.1e-5;

	public MeasurementPerformer() {

	}

	public MeasurementPerformer configure(Qubit resultQubit) {
		this.resultQubit = resultQubit;
		return this;
	}

	public Qubit measure() {
		// collapse state;
		int position = fakeMeasure();
		collapseState(position);
		return resultQubit;
	}

	private void collapseState(int position) {
		ComplexNumber[] collapsedQubit = resultQubit.getQubit();
		for (int i = 0; i < collapsedQubit.length; i++) {
			collapsedQubit[i] = new ComplexNumber(0.0, 0.0);
		}
		collapsedQubit[position] = new ComplexNumber(1.0, 0.0);
		resultQubit = new Qubit(collapsedQubit);
	}

	private int fakeMeasure() {
		length = resultQubit.getQubit().length;
		double[] probabilities = new double[length];
		double value = 0.0;
		double total = 0.0;// total probability
		double choice;
		int position = 0;
		int idx = 0;
		int lastPositionWithNonZeroValue = 0;
		// fill probabilities array and calculate total probability
		for (ComplexNumber c : resultQubit.getQubit()) {
			value = ComplexMath.multiply(c, ComplexMath.conjugate(c)).getReal();
			total += value;
			probabilities[position++] = value;
		}
		assert (Math.abs(total - 1.0) < Q_TOLERANCE);
		// normalize probabilities array
		for (position = 0; position < length; position++) {
			probabilities[position] /= total;
		}
		// use random value to measure state
		choice = new Random(System.currentTimeMillis()).nextDouble();

		for (idx = 0, position = 0; idx < length; ++idx, ++position) {
			value = probabilities[idx];
			if (Double.compare(value, 0.0) != 0) {
				lastPositionWithNonZeroValue = idx;
			}
			if (Double.compare((choice - value), 0.0) == 0 || Double.compare((choice - value), 0.0) < 0) {
				return position;// proper entry, just return it
			}
			choice -= value;// update and continue
		}
		if (Math.abs(choice) < Q_TOLERANCE) {
			return lastPositionWithNonZeroValue;
		}
		return 0;
	}
}
