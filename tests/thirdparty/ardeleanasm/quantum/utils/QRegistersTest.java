package thirdparty.ardeleanasm.quantum.utils;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import thirdparty.ardeleanasm.quantum.exception.RegisterOverflowException;
import thirdparty.ardeleanasm.qubits.QRegister;
import thirdparty.ardeleanasm.qubits.Qubit;
import thirdparty.ardeleanasm.qubits.QubitOne;
import thirdparty.ardeleanasm.qubits.QubitZero;

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
public class QRegistersTest {
	private QRegister			qRegister;
	private static final int	REGISTER_LENGTH	= 3;

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testInitialization() {
		qRegister = new QRegister(REGISTER_LENGTH);
		qRegister = QRegisterOperations.getInstance().fillWith(qRegister, QubitZero::new);
		assertEquals(REGISTER_LENGTH, qRegister.size());
		for (Qubit q : qRegister) {
			assertEquals(new QubitZero(), q);
		}
	}

	@Test
	public void testInit() {
		Qubit q1 = new QubitOne();
		Qubit q0 = new QubitZero();
		try {
			qRegister = QRegisterOperations.getInstance().fillWithPattern("1101");
		} catch (RegisterOverflowException e) {

			e.printStackTrace();
		}
		assertEquals(4, qRegister.size());
		assertEquals(q1, qRegister.get(0));
		assertEquals(q1, qRegister.get(1));
		assertEquals(q0, qRegister.get(2));
		assertEquals(q1, qRegister.get(3));
	}
	

	@Test
	public void testEntangleWithThreeQubitsAllCombinations() {
		QRegister qRegister = new QRegister(3).initialize();
				Qubit qubit = new QubitZero();
		qubit = QuantumOperations.entangle(qubit, new QubitZero());
		qubit = QuantumOperations.entangle(qubit, new QubitZero());
		assertEquals(qubit, QRegisterOperations.getInstance().entangle(qRegister));

	}
}
