package thirdparty.ardeleanasm.gates;

/*
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

/**
 * 
 * Implementation of a Gate Factory.
*/
public class GateFactory extends GatesAbstractFactory {

	/**
	 * getGate return a new <code>Gate</code> object
	 * @param id specify the type of gate.
	 * @return gate 
	 *  
	 */
	public IGate getGate(EGateTypes id) {
		IGate gate = null;
		switch (id) {
		case E_HadamardGate:
			gate = new HGate();
			break;
		case E_XGate:
			gate = new XGate();
			break;
		case E_YGate:
			gate=new YGate();
			break;
		case E_ZGate:
			gate = new ZGate();
			break;
		case E_CNotGate:
			gate = new CNotGate();
			break;
		case E_CPhaseShift:
			gate=new CPhaseShift();
			break;
		case E_IGate:
			gate=new IdentityGate();
			break;
		}
		return gate;
	}
}
