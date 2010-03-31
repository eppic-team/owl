package owl.core.structure;

import java.util.HashMap;
import java.util.Set;

/**
 * A particular catalytic within a protein structure
 */
public class CatalyticSite {
	
	/*--------------------------- member variables --------------------------*/

	HashMap<Integer,String> resser2chemfunc;
	int id;
	String evidence;
	String littEntryPdbCode;
	String littEntryPdbChainCode;
	
	/*----------------------------- constructors ----------------------------*/
	
	public CatalyticSite(int id, String evidence, String littEntryPdbCode, String littEntryPdbChainCode) {
		this.id = id;
		this.evidence = evidence;
		this.littEntryPdbCode = littEntryPdbCode;
		this.littEntryPdbChainCode = littEntryPdbChainCode;
		resser2chemfunc = new HashMap<Integer,String>();
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	public CatalyticSite copy() {
		return new CatalyticSite(id, evidence, littEntryPdbCode, littEntryPdbChainCode);
	}
	
	// add a residue with specific chemical function as a catalytic one
	public void addRes(int resser, String chemFunc) {
		resser2chemfunc.put(resser, chemFunc);
	}
	
	// remove a residue
	public void remRes(int resser) {
		resser2chemfunc.remove(resser);
	}
	
	// get all residues of the catalytic site
	public Set<Integer> getRes() {
		return resser2chemfunc.keySet();
	}

	// get the chemical function from residue serial
	public String getChemFuncFromResSerial(int resser) {
		return resser2chemfunc.get(resser);
	}
	
	/** Returns the serial for this site */ 
	public int getSerial() {
		return id;
	}
	
	/** Returns the evidence for this site */ 
	public String getEvidence() {
		return evidence;
	}
	
	/** Returns the literature entry pdb code for this site */ 
	public String getLittEntryPdbCode() {
		return littEntryPdbCode;
	}
	
	/** Returns the literature entry pdb chain code for this site */ 
	public String getLittEntryPdbChainCode() {
		return littEntryPdbChainCode;
	}
	
	public void print() {
		for(int resSer:resser2chemfunc.keySet()) {
			System.out.println("site:"+id+"-res:"+resSer+"-func:"+getChemFuncFromResSerial(resSer));
		}
	}
	
	public String toString() {
		String msg = "";
		for(int resSer:resser2chemfunc.keySet()) {
			msg += String.format("Site:%3d Res:%4d Function:%2s \n", id, resSer, getChemFuncFromResSerial(resSer));
		}
		return msg;
	}
}

