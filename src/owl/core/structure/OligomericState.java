package owl.core.structure;

import java.util.HashSet;
import java.util.regex.Pattern;

/** 
 * This class encapsulates the oligomeric state as predicted by PQS and PISA 
 */
public class OligomericState {

	/*--------------------------- member variables --------------------------*/
	
	private String pqsBiolistState;	// the number of molecules in mmol as defined by PQS BIOLIST (can be null or number or string)
	private String pqsAsalistState;	// the number of molecules in mmol as defined by PQS ASALIST (can be null or number or list of comma-separated numbers)
	private int pisaState;			// the number of molecules in mmol as defined by PISA (if 0 undefined)
	
	/*----------------------------- constructors ----------------------------*/
	
	public OligomericState() {
		this.pqsBiolistState = null;
		this.pqsAsalistState = null;
		this.pisaState = 0;		
	}
	
	public OligomericState(String pqsBiolistState, String pqsAsalistState, int pisaState) {
		this.pqsBiolistState = pqsBiolistState;
		this.pqsAsalistState = pqsAsalistState;
		this.pisaState = pisaState;
	}	
	
	/*---------------------------- public methods ---------------------------*/
	
	public void setPqs(String biolistState, String asalistState) {
		this.pqsBiolistState = biolistState;
		this.pqsAsalistState = asalistState;
	}
	
	public void setPisa(int state) {
		this.pisaState = state;
	}
	
	public String getPqsBiolistState() {
		return pqsBiolistState;
	}
	
	public String getPqsAsalistState() {
		return pqsAsalistState;
	}
	
	public int getPisaState() {
		return pisaState;
	}
	
	public boolean hasPisa() {
		return (pisaState <= 0)?false:true;
	}
	
	public boolean hasPQS() {
		if ((pqsAsalistState != null) || (pqsBiolistState != null && Pattern.matches("^\\d+$", pqsBiolistState))) {
			return true;
		} else {
			return false;
		}
	}
	
	public int getPqsState() {
		if ((pqsAsalistState != null) && (pqsBiolistState != null && Pattern.matches("^\\d+$", pqsBiolistState))) {
			String[] asaStates = pqsAsalistState.split(",");
			HashSet<Integer> asaStatesSet = new HashSet<Integer>();
			for(String asaState: asaStates) {
				asaStatesSet.add(Integer.valueOf(asaState));
			}
			int bioState = Integer.valueOf(pqsBiolistState).intValue();
			if (asaStatesSet.contains(bioState)) {
				return bioState;
			} else {
				return Integer.valueOf(asaStates[0]);
			}
		}
		else if (pqsAsalistState != null) {
			return Integer.valueOf((pqsAsalistState.split(","))[0]).intValue();
		} else if (pqsBiolistState != null && Pattern.matches("^\\d+$", pqsBiolistState)) {
			return Integer.valueOf(pqsBiolistState).intValue();
		} else {
			return 0;
		}
	}
	
	public int getState() {
		if (hasPisa()) {
			return getPisaState();
		} else if (hasPQS()) {
			return getPqsState();
		} else {
			return 0;
		}
	}
}
