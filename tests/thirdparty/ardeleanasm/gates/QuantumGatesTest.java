package thirdparty.ardeleanasm.gates;


import static org.junit.Assert.*;


import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import thirdparty.ardeleanasm.complexnumbers.ComplexNumber;
import thirdparty.ardeleanasm.quantum.exception.RegisterOverflowException;
import thirdparty.ardeleanasm.quantum.utils.QRegisterOperations;
import thirdparty.ardeleanasm.quantum.utils.QuantumOperations;
import thirdparty.ardeleanasm.qubits.QRegister;
import thirdparty.ardeleanasm.qubits.Qubit;
import thirdparty.ardeleanasm.qubits.QubitPlus;

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
public class QuantumGatesTest {
	private GatesAbstractFactory factory;
	private IGate gate;

	@Before
	public void setUp() throws Exception {
		factory = GateProducer.getGateFactory();
	}

	@After
	public void tearDown() throws Exception {
	}


	@Test
	public void testXGateQubitZero() {
		QRegister targetRegister=new QRegister(3);
		QRegister expectedRegister=new QRegister(3);
		QRegisterOperations regOps=QRegisterOperations.getInstance();
		gate=factory.getGate(EGateTypes.E_XGate);	
		
		try {
			targetRegister=regOps.fillWithPattern("000");
			expectedRegister=regOps.fillWithPattern("100");
		} catch (RegisterOverflowException e) {
			e.printStackTrace();
		}
		
		Qubit target = regOps.entangle(targetRegister); 
		target=gate.applyGate(target, new int[]{0}, null, 3);	
		assertEquals(target,regOps.entangle(expectedRegister));

	}
	
	@Test
	public void testXGateQubitOne() {
		QRegister targetRegister=new QRegister(3);
		QRegister expectedRegister=new QRegister(3);
		QRegisterOperations regOps=QRegisterOperations.getInstance();
		gate=factory.getGate(EGateTypes.E_XGate);	
		
		try {
			targetRegister=regOps.fillWithPattern("111");
			expectedRegister=regOps.fillWithPattern("110");
		} catch (RegisterOverflowException e) {
			e.printStackTrace();
		}
		
		Qubit target = regOps.entangle(targetRegister); 
		target=gate.applyGate(target, new int[]{2}, null, 3);	
		assertEquals(target ,regOps.entangle(expectedRegister));

	}
	

	
	@Test
	public void testXGateTwoQubits() {
		QRegister targetRegister=new QRegister(3);
		QRegister expectedRegister=new QRegister(3);
		QRegisterOperations regOps=QRegisterOperations.getInstance();
		gate=factory.getGate(EGateTypes.E_XGate);	
		
		try {
			targetRegister=regOps.fillWithPattern("111");
			expectedRegister=regOps.fillWithPattern("001");
		} catch (RegisterOverflowException e) {
			e.printStackTrace();
		}
		
		Qubit target = regOps.entangle(targetRegister); 
		target=gate.applyGate(target, new int[]{0,1}, null, 3);	
		assertEquals(target ,regOps.entangle(expectedRegister));

	}

	
	@Test
	public void testYGateQubitZero(){
		QRegister targetRegister=new QRegister(3);
		QRegister expectedRegister=new QRegister(3);
		QRegisterOperations regOps=QRegisterOperations.getInstance();
		Qubit q=new Qubit(new ComplexNumber(0.0,0.0),new ComplexNumber(0.0,1.0));
		gate=factory.getGate(EGateTypes.E_YGate);	
		
		try {
			targetRegister=regOps.fillWithPattern("000");
			expectedRegister=regOps.fillWithPattern("000");
		} catch (RegisterOverflowException e) {
			e.printStackTrace();
		}
		expectedRegister.change(0, q);
		Qubit target = regOps.entangle(targetRegister); 
		target=gate.applyGate(target, new int[]{0}, null, 3);	
		assertEquals(target , regOps.entangle(expectedRegister));
	}
	
	@Test
	public void testYGateQubitOne(){
		QRegister targetRegister=new QRegister(3);
		QRegister expectedRegister=new QRegister(3);
		QRegisterOperations regOps=QRegisterOperations.getInstance();
		Qubit q=new Qubit(new ComplexNumber(0.0,-1.0),new ComplexNumber(0.0,0.0));
		gate=factory.getGate(EGateTypes.E_YGate);	
		
		try {
			targetRegister=regOps.fillWithPattern("111");
			expectedRegister=regOps.fillWithPattern("111");
		} catch (RegisterOverflowException e) {
			e.printStackTrace();
		}
		expectedRegister.change(0, q);
		Qubit target = regOps.entangle(targetRegister); 
		target=gate.applyGate(target, new int[]{0}, null, 3);	
		assertEquals(target , regOps.entangle(expectedRegister));
	}

	@Test
	public void testCnotGateQubitOne() {
		QRegister targetRegister=new QRegister(3);
		QRegister expectedRegister=new QRegister(3);
		QRegisterOperations regOps=QRegisterOperations.getInstance();
		gate=factory.getGate(EGateTypes.E_CNotGate);	
		
		try {
			targetRegister=regOps.fillWithPattern("111");
			expectedRegister=regOps.fillWithPattern("110");
		} catch (RegisterOverflowException e) {
			e.printStackTrace();
		}
		
		Qubit target = regOps.entangle(targetRegister); 
		target=gate.applyGate(target, new int[]{2}, new int[]{0,1}, 3);	
		assertEquals(target , regOps.entangle(expectedRegister));

	}
	
	@Test
	public void testCnotGateQubitZero() {
		QRegister targetRegister=new QRegister(3);
		QRegister expectedRegister=new QRegister(3);
		QRegisterOperations regOps=QRegisterOperations.getInstance();
		gate=factory.getGate(EGateTypes.E_CNotGate);	
		
		try {
			targetRegister=regOps.fillWithPattern("110");
			expectedRegister=regOps.fillWithPattern("111");
		} catch (RegisterOverflowException e) {
			e.printStackTrace();
		}
		
		Qubit target = regOps.entangle(targetRegister); 
		target=gate.applyGate(target, new int[]{2}, new int[]{0,1}, 3);	
		assertEquals(target ,regOps.entangle(expectedRegister));

	}
	
	@Test
	public void testCnotGateAllQubitZero() {
		QRegister targetRegister=new QRegister(3);
		QRegister expectedRegister=new QRegister(3);
		QRegisterOperations regOps=QRegisterOperations.getInstance();
		gate=factory.getGate(EGateTypes.E_CNotGate);	
		
		try {
			targetRegister=regOps.fillWithPattern("000");
			expectedRegister=regOps.fillWithPattern("000");
		} catch (RegisterOverflowException e) {
			e.printStackTrace();
		}
		
		Qubit target = regOps.entangle(targetRegister); 
		target=gate.applyGate(target, new int[]{2}, new int[]{0,1}, 3);	
		assertEquals(target ,regOps.entangle(expectedRegister));

	}
	
	@Test
	public void testCnotGateFirstQubitZero() {
		QRegister targetRegister=new QRegister(3);
		QRegister expectedRegister=new QRegister(3);
		QRegisterOperations regOps=QRegisterOperations.getInstance();
		gate=factory.getGate(EGateTypes.E_CNotGate);	
		
		try {
			targetRegister=regOps.fillWithPattern("010");
			expectedRegister=regOps.fillWithPattern("010");
		} catch (RegisterOverflowException e) {
			e.printStackTrace();
		}
		
		Qubit target = regOps.entangle(targetRegister); 
		target=gate.applyGate(target, new int[]{2}, new int[]{0,1}, 3);	
		assertEquals(target , regOps.entangle(expectedRegister));

	}
	@Test
	public void testZGateQubitZero(){
		QRegister targetRegister=new QRegister(3);
		QRegister expectedRegister=new QRegister(3);
		QRegisterOperations regOps=QRegisterOperations.getInstance();
		Qubit q=new Qubit(new ComplexNumber(1.0,0.0),new ComplexNumber(0.0,0.0));
		gate=factory.getGate(EGateTypes.E_ZGate);	
		
		try {
			targetRegister=regOps.fillWithPattern("000");
			expectedRegister=regOps.fillWithPattern("000");
		} catch (RegisterOverflowException e) {
			e.printStackTrace();
		}
		expectedRegister.change(0, q);
		Qubit target = regOps.entangle(targetRegister); 
		target=gate.applyGate(target, new int[]{0}, null, 3);	
		
		
		//assertEquals(target , QuantumOperations.entangle(expectedRegister));
	}
	
	@Test
	public void testZGateQubitOne(){
		QRegister targetRegister=new QRegister(3);
		QRegister expectedRegister=new QRegister(3);
		QRegisterOperations regOps=QRegisterOperations.getInstance();
		Qubit q=new Qubit(new ComplexNumber(0.0,0.0),new ComplexNumber(-1.0,0.0));
		gate=factory.getGate(EGateTypes.E_ZGate);	
		
		try {
			targetRegister=regOps.fillWithPattern("111");
			expectedRegister=regOps.fillWithPattern("111");
		} catch (RegisterOverflowException e) {
			e.printStackTrace();
		}
		expectedRegister.change(0, q);
		Qubit target = regOps.entangle(targetRegister); 
		target=gate.applyGate(target, new int[]{0}, null, 3);	
		assertEquals(target ,regOps.entangle(expectedRegister));
	}
	
	@Test
	public void testApplyHadamardGate(){
		QRegister targetRegister=new QRegister(3);
		QRegister expectedRegister=new QRegister(3);
		QRegisterOperations regOps=QRegisterOperations.getInstance();
		gate=factory.getGate(EGateTypes.E_HadamardGate);	
		
		try {
			targetRegister=regOps.fillWithPattern("000");
			expectedRegister=regOps.fillWithPattern("000");
		} catch (RegisterOverflowException e) {
			e.printStackTrace();
		}
		expectedRegister.change(0, new QubitPlus());
		Qubit target = regOps.entangle(targetRegister); 
		target=gate.applyGate(target, new int[]{0}, null, 3);	
		assertEquals(target , regOps.entangle(expectedRegister));
	}
	
	@Test
	public void testApplyIdentityGate(){
		QRegister targetRegister=new QRegister(3);
		QRegister expectedRegister=new QRegister(3);
		QRegisterOperations regOps=QRegisterOperations.getInstance();
		gate=factory.getGate(EGateTypes.E_IGate);	
		
		try {
			targetRegister=regOps.fillWithPattern("000");
			expectedRegister=regOps.fillWithPattern("000");
		} catch (RegisterOverflowException e) {
			e.printStackTrace();
		}
		
		Qubit target = regOps.entangle(targetRegister); 
		target=gate.applyGate(target, null, null, 3);	
		assertEquals(target , regOps.entangle(expectedRegister));
	}
	
	
}
