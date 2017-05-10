package thirdparty.ardeleanasm.complexnumbers;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import thirdparty.ardeleanasm.complexnumbers.ComplexNumber;

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
public class ComplexMathTest {
	private ComplexNumber		firstNumber;
	private ComplexNumber		secondNumber;
	private static final Double	DELTA						= 5e-10;
	private static final Double	REAL_VALUE_FIRST_NO			= 2.0;
	private static final Double	IMAGINARY_VALUE_FIRST_NO	= 3.0;
	private static final Double	REAL_VALUE_SECOND_NO		= 4.0;
	private static final Double	IMAGINARY_VALUE_SECOND_NO	= 5.0;

	@Before
	public void setUp() throws Exception {
		firstNumber = new ComplexNumber(REAL_VALUE_FIRST_NO, IMAGINARY_VALUE_FIRST_NO);
		secondNumber = new ComplexNumber(REAL_VALUE_SECOND_NO, IMAGINARY_VALUE_SECOND_NO);
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testConjugate() {
		ComplexNumber expectedNumber = new ComplexNumber(REAL_VALUE_FIRST_NO, -IMAGINARY_VALUE_FIRST_NO);
		ComplexNumber realNumber = null;
		realNumber = ComplexMath.conjugate(firstNumber);
		assertEquals(expectedNumber, realNumber);
	}

	@Test
	public void testMod() {
		double expectedNumber = Math.sqrt(Math.pow(REAL_VALUE_FIRST_NO, 2) + Math.pow(IMAGINARY_VALUE_FIRST_NO, 2));
		double realNumber = 0.0;
		realNumber = ComplexMath.mod(firstNumber);
		assertEquals(expectedNumber, realNumber, DELTA);
	}

	@Test
	public void testAdd() {
		ComplexNumber expectedNumber = new ComplexNumber();
		ComplexNumber realNumber = null;

		expectedNumber.setReal(REAL_VALUE_FIRST_NO + REAL_VALUE_SECOND_NO);
		expectedNumber.setImaginary(IMAGINARY_VALUE_FIRST_NO + IMAGINARY_VALUE_SECOND_NO);

		realNumber = ComplexMath.add(firstNumber, secondNumber);

		assertEquals(expectedNumber, realNumber);

	}

	@Test
	public void testSubtract() {
		ComplexNumber expectedNumber = new ComplexNumber();
		ComplexNumber realNumber = null;

		expectedNumber.setReal(REAL_VALUE_FIRST_NO - REAL_VALUE_SECOND_NO);
		expectedNumber.setImaginary(IMAGINARY_VALUE_FIRST_NO - IMAGINARY_VALUE_SECOND_NO);

		realNumber = ComplexMath.subtract(firstNumber, secondNumber);

		assertEquals(expectedNumber, realNumber);

	}

	@Test
	public void testMultiply() {
		/*
		 * a=2+3i b=4+5i c=(2+3i)(4+5i)=8+10i+12i-15=-7+22i
		 * 
		 */
		ComplexNumber expectedNumber = new ComplexNumber(-7.0, 22.0);
		ComplexNumber realValue = null;

		realValue = ComplexMath.multiply(firstNumber, secondNumber);
		assertEquals(true, realValue.equals(expectedNumber));
	}

	@Test
	public void testDivide() {
		/*
		 * a=2+3i b=4+5i
		 * c=a/b=(2+3i)(4-5i)/(16+25)=(8+15-10i+12i)/(41)=(23+2i)/41
		 * 
		 */
		ComplexNumber expectedNumber = new ComplexNumber(23 / 41.0, 2 / 41.0);
		ComplexNumber realNumber = null;

		realNumber = ComplexMath.divide(firstNumber, secondNumber);
		assertEquals(expectedNumber, realNumber);

	}

	@Test
	public void testSquare() {
		/*
		 * a=2+3i a^2=(2+3i)(2+3i)=4-9i+6i+6i=-5+12i
		 * 
		 */
		ComplexNumber expectedNumber = new ComplexNumber(-5, 12);
		ComplexNumber realNumber = null;

		realNumber = ComplexMath.square(firstNumber);
		assertEquals(expectedNumber, realNumber);
	}

	@Test
	public void testSine() {
		/*
		 * a=a1+a2i; x=e^a2 r=sin(a1)*(x+1/x)/2; ip=cos(a1)*(x-1/x)/2;
		 * sin(a)=r+ip*i;
		 * 
		 * if a=2+3i r=sin(2)+(e^3+1/e^3)/2 ip=cos(2)*(e^2+1/e^2)/2
		 */
		double x = Math.exp(IMAGINARY_VALUE_FIRST_NO);
		double r = Math.sin(REAL_VALUE_FIRST_NO) * (x + 1 / x) / 2.0;
		double ip = Math.cos(REAL_VALUE_FIRST_NO) * (x - 1 / x) / 2;
		ComplexNumber expectedNumber = new ComplexNumber(r, ip);
		ComplexNumber realNumber = null;

		realNumber = ComplexMath.sin(firstNumber);
		assertEquals(expectedNumber, realNumber);
	}

	@Test
	public void testCosine() {
		/*
		 * a=a1+a2i; x=e^a2 r=cos(a1)*(x+1/x)/2; ip=-sin(a1)*(x-1/x)/2;
		 * sin(a)=r+ip*i;
		 * 
		 * if a=2+3i r=cos(2)+(e^3+1/e^3)/2 ip=-sin(2)*(e^2+1/e^2)/2
		 */
		double x = Math.exp(IMAGINARY_VALUE_FIRST_NO);
		double r = Math.cos(REAL_VALUE_FIRST_NO) * (x + 1 / x) / 2.0;
		double ip = -Math.sin(REAL_VALUE_FIRST_NO) * (x - 1 / x) / 2;
		ComplexNumber expectedNumber = new ComplexNumber(r, ip);
		ComplexNumber realNumber = null;

		realNumber = ComplexMath.cos(firstNumber);
		assertEquals(expectedNumber, realNumber);
	}

	@Test
	public void testTan() {
		ComplexNumber expectedNumber = ComplexMath.divide(ComplexMath.sin(firstNumber), ComplexMath.cos(firstNumber));
		ComplexNumber realNumber;
		realNumber = ComplexMath.tan(firstNumber);
		assertEquals(expectedNumber, realNumber);
	}

	@Test
	public void testExp() {
		/*
		 * a=a1+a2i; b=b1+b2i r=e^a1
		 * 
		 * b1=r*cos(a2) b2=r*sin(a2)
		 */

		double r = Math.exp(REAL_VALUE_FIRST_NO);

		ComplexNumber expectedNumber = new ComplexNumber(r * Math.cos(IMAGINARY_VALUE_FIRST_NO),
				r * Math.sin(IMAGINARY_VALUE_FIRST_NO));
		ComplexNumber realNumber = null;

		realNumber = ComplexMath.exp(firstNumber);
		assertEquals(expectedNumber, realNumber);
	}

	@Test
	public void testPow() {
		ComplexNumber expectedNumber = ComplexMath.multiply(ComplexMath.square(firstNumber), firstNumber);
		ComplexNumber realNumber = null;

		realNumber = ComplexMath.pow(firstNumber, 3);
		assertEquals(expectedNumber, realNumber);
	}

	@Test
	public void testGetArg() {
		double expectedArg = Math.atan2(firstNumber.getImaginary(), firstNumber.getReal());
		double realArg = firstNumber.getArg();
		assertEquals(expectedArg, realArg, DELTA);
	}

	@Test
	public void testInverse() {
		ComplexNumber expectedNumber = ComplexMath.divide(new ComplexNumber(1, 0), firstNumber);
		ComplexNumber realNumber = null;

		realNumber = ComplexMath.inverse(firstNumber);
		assertEquals(expectedNumber, realNumber);
	}
	
	@Test
	public void testMultiplyConstant(){
		double constant=10.0;
		ComplexNumber expectedNumber = new ComplexNumber(REAL_VALUE_FIRST_NO*constant,IMAGINARY_VALUE_FIRST_NO*constant);
		ComplexNumber realNumber = null;

		realNumber = ComplexMath.multiply(firstNumber,constant);
		assertEquals(expectedNumber, realNumber);

	}

	@Test
	public void testDivisionConstant(){
		double constant=10.0;
		ComplexNumber expectedNumber = new ComplexNumber(REAL_VALUE_FIRST_NO/constant,IMAGINARY_VALUE_FIRST_NO/constant);
		ComplexNumber realNumber = null;

		realNumber = ComplexMath.divide(firstNumber,constant);
		assertEquals(expectedNumber, realNumber);

	}
	
	@Test
	public void testAddConstant(){
		double constant=10.0;
		ComplexNumber expectedNumber = new ComplexNumber(REAL_VALUE_FIRST_NO+constant,IMAGINARY_VALUE_FIRST_NO);
		ComplexNumber realNumber = null;

		realNumber = ComplexMath.add(firstNumber,constant);
		assertEquals(expectedNumber, realNumber);

	}
	
	@Test
	public void testSubtractConstant(){
		double constant=10.0;
		ComplexNumber expectedNumber = new ComplexNumber(REAL_VALUE_FIRST_NO-constant,IMAGINARY_VALUE_FIRST_NO);
		ComplexNumber realNumber = null;

		realNumber = ComplexMath.subtract(firstNumber,constant);
		assertEquals(expectedNumber, realNumber);

	}
}
