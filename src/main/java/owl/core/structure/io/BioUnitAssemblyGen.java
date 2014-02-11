package owl.core.structure.io;

import java.util.ArrayList;

/**
 * * A class to temporarily store the read in data of a biounit while parsing the pdb file
 * 
 * NOTE:
 * To store the data of a pdb or cif file assembly three temporary classes are maintained
 * BioUnitAssembly    : Contains details on the size, type and links to its generator
 * BioUnitAssemblyGen : Contains details on chains and its linked operators
 * BioUnitOperation   : Contains details on the operators of type Matrix4d
 * 
 * @author biyani_n
 *
 */
public class BioUnitAssemblyGen {
	
	private int id;
	private ArrayList<String> pdbChainCodes;
	private ArrayList<Integer> operationIds;
	private boolean reading;
	
	
	//Default constructor
	public BioUnitAssemblyGen(){
		this.id = 0;
		this.pdbChainCodes = new ArrayList<String>();
		this.operationIds = new ArrayList<Integer>();
		this.reading = false;
	}
	
	public BioUnitAssemblyGen(int id){
		this.id = id;
		this.pdbChainCodes = new ArrayList<String>();
		this.operationIds = new ArrayList<Integer>();
		this.reading = false;
	}
	
	//public methods
	public int getId() {
		return id;
	}
	
	public void setId(int id) {
		this.id = id;
	}
	
	public ArrayList<String> getPdbChainCodes() {
		return pdbChainCodes;
	}
	
	public void setPdbChainCodes(ArrayList<String> pdbChainCodes) {
		this.pdbChainCodes = pdbChainCodes;
	}
	
	public void addPdbChainCodes(String[] codes){
		for(String code:codes) this.pdbChainCodes.add(code.trim());
	}
	
	public ArrayList<Integer> getOperationIds() {
		return operationIds;
	}
	
	public void setOperationIds(ArrayList<Integer> operationIds) {
		this.operationIds = operationIds;
	}
	
	public void addOperationId(int id){
		this.operationIds.add(id);
	}
	
	public boolean isReading() {
		return reading;
	}
	
	public void setReading(boolean reading) {
		this.reading = reading;
	}

	

}
