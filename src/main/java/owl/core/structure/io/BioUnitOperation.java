package owl.core.structure.io;

/**
 * A class to temporarily store the read in data of a biounit operators while parsing the pdb file
 * 
 * NOTE:
 * To store the data of a pdb or cif file assembly three temporary classes are maintained
 * BioUnitAssembly    : Contains details on the size, type and links to its generator
 * BioUnitAssemblyGen : Contains details on chains and its linked operators
 * BioUnitOperation   : Contains details on the operators of type Matrix4d
 * @author biyani_n
 *
 */
public class BioUnitOperation {

	/**
	 * @param args
	 */
	private int id;
	private double[] operator;
	
	//Default constructor
	public BioUnitOperation(){
		this.id = 0;
		this.operator = new double[16]; operator[15]=1;
	}
	
	public BioUnitOperation(int id){
		this.id = id;
		this.operator = new double[16]; operator[15]=1;
	}
	
	//Getters and setters
	public int getId() {
		return id;
	}
	public void setId(int id) {
		this.id = id;
	}
	public double[] getOperator() {
		return operator;
	}
	public void setOperatorValue(int index, double value) {
		this.operator[index] = value;
	}
	
}
