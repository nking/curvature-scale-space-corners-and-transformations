package thirdparty.ardeleanasm.qubits;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import thirdparty.ardeleanasm.quantum.exception.RegisterOverflowException;

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
public class QRegisterTest {
	private QRegister			qRegister;
	private static final int	REGISTER_LENGTH	= 3;

	@Before
	public void setUp() throws Exception {

	}

	@After
	public void tearDown() throws Exception {
		qRegister = null;
	}

	@Test
	public void testSizeNoElementsWithoutInit() {
		qRegister = new QRegister(REGISTER_LENGTH);
		assertEquals(3, qRegister.size());
	}

	@Test
	public void testAddElementsWithoutInit() {
		qRegister = new QRegister(REGISTER_LENGTH);
		try {
			qRegister.add(new QubitOne());
		} catch (RegisterOverflowException e) {
			e.printStackTrace();
		}
		assertEquals(REGISTER_LENGTH , qRegister.size());
	}

	@Test
	public void testInitialization() {
		qRegister = new QRegister(REGISTER_LENGTH).initialize();
		assertEquals(REGISTER_LENGTH, qRegister.size());
		for (Qubit q : qRegister) {
			assertEquals(new QubitZero(), q);
		}
	}

	@Test(expected = RegisterOverflowException.class)
	public void testAddElementException() throws RegisterOverflowException {
		qRegister = new QRegister(REGISTER_LENGTH).initialize();
		qRegister.add(new QubitOne());
	}

	@Test
	public void testGetWithoutInit() throws RegisterOverflowException {
		qRegister = new QRegister(REGISTER_LENGTH);
		qRegister.add(new QubitOne());
		qRegister.add(new QubitOne());
		qRegister.add(new QubitZero());
		assertEquals(new QubitZero(), qRegister.get(REGISTER_LENGTH - 1));
	}

	@Test
	public void testChange() {
		qRegister = new QRegister(REGISTER_LENGTH).initialize();
		qRegister = qRegister.change(1, new QubitOne());
		assertEquals(new QubitOne(), qRegister.get(1));
	}

}
