package proteinstructure;

public class AIGNode {
	
	private int atomSerial;
	private String atomName;
	private RIGNode parentResidue;

	public AIGNode(int atomSerial, String atomName, RIGNode parentResidue) {
		this.atomSerial = atomSerial;
		this.atomName = atomName;
		this.parentResidue = parentResidue;
	}
	
	public int getAtomSerial() {
		return this.atomSerial;
	}
	
	public String getAtomName() {
		return this.atomName;
	}
	
	public RIGNode getParent() {
		return this.parentResidue;
	}
	
	/**
	 * Convenience method to get the residue serial for this atom node
	 * @return
	 */
	public int getResidueSerial() {
		return this.parentResidue.getResidueSerial();
	}
}
