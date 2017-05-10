package thirdparty.ardeleanasm.complexnumbers;

/**
* <h1>ComplexMath</h1>
* ComplexMath implements the basic operations that can be performed on complex numbers. 
*
* @author  Mihai Seba
* @version 1.0

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
public class ComplexMath {

	/**
	 * Performs the addition of two complex numbers.
	 * @param z1 first complex number
	 * @param z2 second complex number
	 * @return the sum of the two complex numbers, z1 and z2, z3=z1+z2
	 */
	public static ComplexNumber add(ComplexNumber z1, ComplexNumber z2) {
		return new ComplexNumber(z1.getReal() + z2.getReal(), z1.getImaginary() + z2.getImaginary());
	}

	/**
	 * Performs the multiplication of two complex numbers.
	 * @param z1 first complex number
	 * @param z2 second complex number
	 * @return the result of multiplication, z1 and z2, z3=z1*z2
	 */
	public static ComplexNumber multiply(ComplexNumber z1, ComplexNumber z2) {
		return new ComplexNumber(z1.getReal() * z2.getReal() - z1.getImaginary() * z2.getImaginary(),
				z1.getReal() * z2.getImaginary() + z1.getImaginary() * z2.getReal());
	}

	/**
	 * Performs the subtraction of two complex numbers.
	 * @param z1 first complex number
	 * @param z2 second complex number
	 * @return the result of subtraction, z1 and z2, z3=z1-z2
	 */
	public static ComplexNumber subtract(ComplexNumber z1, ComplexNumber z2) {

		return new ComplexNumber(z1.getReal() - z2.getReal(), z1.getImaginary() - z2.getImaginary());
	}

	/**
	 * Performs the division of two complex numbers.
	 * @param z1 first complex number
	 * @param z2 second complex number
	 * @return the result of division, z1 and z2, z3=z1/z2
	 */
	public static ComplexNumber divide(ComplexNumber z1, ComplexNumber z2) {
		ComplexNumber numerator = multiply(z1, conjugate(z2));
		double denominator = Math.pow(mod(z2), 2);
		return new ComplexNumber(numerator.getReal() / denominator, numerator.getImaginary() / denominator);
	}

	/**
	 * Calculates the complex conjugate of a <code>ComplexNumber</code>.
	 * @param z the input complex number
	 * @return the complex conjugate number of z
	 */
	public static ComplexNumber conjugate(ComplexNumber z) {

		return new ComplexNumber(z.getReal(), -z.getImaginary());
	}

	/**
	 * The modulus, magnitude or the absolute value of current complex number.
	 * @param z the input complex number
	 * @return the modulus of number z1
	 */
	public static double mod(ComplexNumber z) {

		return Math.sqrt(Math.pow(z.getReal(), 2) + Math.pow(z.getImaginary(), 2));
	}

	/**
	 * Calculates the square of the <code>ComplexNumber</code>.
	 * @param z the input complex number
	 * @return the square of complex number z
	 */
	public static ComplexNumber square(ComplexNumber z) {
		return new ComplexNumber(z.getReal() * z.getReal() - z.getImaginary() * z.getImaginary(),
				2 * z.getReal() * z.getImaginary());
	}

	/**
	 * Calculates the sine of the <code>ComplexNumber</code>
	 * @param z the input complex number
	 * @return return a <code>ComplexNumber</code> which is sine of z
	 */
	public static ComplexNumber sin(ComplexNumber z) {
		double exp = Math.exp(z.getImaginary());
		return new ComplexNumber(Math.sin(z.getReal()) * (exp + 1 / exp) / 2.0,
				Math.cos(z.getReal()) * (exp - 1 / exp) / 2.0);
	}

	/**
	 * Calculates the cosine of the <code>ComplexNumber</code>
	 * @param z the input complex number
	 * @return return a <code>ComplexNumber</code> which is cosine of z
	 */
	public static ComplexNumber cos(ComplexNumber z) {
		double exp = Math.exp(z.getImaginary());
		return new ComplexNumber(Math.cos(z.getReal()) * (exp + 1 / exp) / 2.0,
				-Math.sin(z.getReal()) * (exp - 1 / exp) / 2.0);
	}

	/**
	 * Calculates the tangent of the <code>ComplexNumber</code>
	 * @param z the input complex number
	 * @return return a <code>ComplexNumber</code> which is tangent of z
	 */
	public static ComplexNumber tan(ComplexNumber z) {
		return divide(sin(z), cos(z));
	}

	/**
	 * Calculates the exponential  of the <code>ComplexNumber</code>
	 * @param z the input complex number
	 * @return return a <code>ComplexNumber</code> which is e^
	 */
	public static ComplexNumber exp(ComplexNumber z) {
		double r = Math.exp(z.getReal());
		return new ComplexNumber(r * Math.cos(z.getImaginary()), r * Math.sin(z.getImaginary()));
	}

	/**
	 * Calculates the <code>ComplexNumber</code> to the passed integer power.
	 * @param z the input complex number
	 * @param power the power
	 * @return return a <code>ComplexNumber</code> which is z^power
	 */
	public static ComplexNumber pow(ComplexNumber z, int power) {
		double realValue = z.getReal();
		double imaginaryValue = z.getImaginary();
		for (int i = 0; i < power - 1; i++) {
			double newRealValue = realValue * z.getReal() - imaginaryValue * z.getImaginary();
			double newImaginaryValue = realValue * z.getImaginary() + imaginaryValue * z.getReal();
			realValue = newRealValue;
			imaginaryValue = newImaginaryValue;
		}
		return new ComplexNumber(realValue, imaginaryValue);

	}
	/**
	 * Calculates the inverse of a <code>ComplexNumber</code>.
	 * @param z the input complex number
	 * @return return a <code>ComplexNumber</code> which is 1/z
	 */
	public static ComplexNumber inverse(ComplexNumber z) {
		return divide(new ComplexNumber(1, 0), z);
	}

	/**
	 * Performs the multiplication of a complex number with a double.
	 * @param z complex number
	 * @param constant double number
	 * @return the result of multiplication.
	 */
	public static ComplexNumber multiply(ComplexNumber z,double constant) {
		return new ComplexNumber(z.getReal()*constant,z.getImaginary()*constant);
	}
	
	/**
	 * Performs the addition of a complex number with a double.
	 * @param z complex number
	 * @param constant double number
	 * @return the result of addition.
	 */
	public static ComplexNumber add(ComplexNumber z, double constant) {
		return new ComplexNumber(z.getReal() + constant, z.getImaginary() );
	}


	/**
	 * Performs the subtraction of a complex number with a double.
	 * @param z complex number
	 * @param constant double number
	 * @return the result of subtraction.
	 */
	public static ComplexNumber subtract(ComplexNumber z, double constant) {

		return new ComplexNumber(z.getReal() - constant, z.getImaginary());
	}



	/**
	 * Performs the division of a complex number with a double.
	 * @param z complex number
	 * @param constant double number
	 * @return the result of division.
	 */
	public static ComplexNumber divide(ComplexNumber z, double constant) {
		return new ComplexNumber(z.getReal() / constant, z.getImaginary() / constant);
	}
}
