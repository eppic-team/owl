/**
 * 
 */
package owl.core.structure.io;

import java.util.ArrayList;

/**
 * A class to temporarily store the read in data of a biounit while parsing the pdb file
 * 
 * NOTE:
 * To store the data of a pdb or cif file assembly two temporary classes are maintained
 * BioUnitAssembly    : Contains details on the size, type and links to its generator
 * BioUnitAssemblyGen : Contains details on chains and its linked operators
 * 
 * @author biyani_n
 *
 */
public class BioUnitAssembly {

	//Members
	private int id;
	private int size;
	private ArrayList<String> types;
	private ArrayList<Integer> generatorIds;
	private boolean reading;
	
	//Methods
	/**
	 * Main Constructor
	 */
	public BioUnitAssembly(){
		this.id = 0;
		this.size=-1;
		this.types = new ArrayList<String>();
		this.generatorIds = new ArrayList<Integer>();
		this.reading = false;
	}
	


	public void setId(int id){
		this.id = id;
	}
	
	public int getId() {
		return id;
	}
	
	public void setSize(int size){
		this.size = size;
	}
	
	/**
	 * Sets the size of biounit from a string
	 * @param oligomer
	 */
	public void setSize(String oligomer){
		this.size = (getSizefromString(oligomer));
	}
	
	public int getSize() {
		return size;
	}
	
	public void addType(String method){
		this.types.add(method);
	}
	
	public ArrayList<String> getTypes() {
		return types;
	}
	
	public void addGeneratorId(Integer id){
		this.generatorIds.add(id);
	}

	public ArrayList<Integer> getGeneratorIds() {
		return generatorIds;
	}

	public boolean isReading() {
		return reading;
	}

	public void setReading(boolean reading) {
		this.reading = reading;
	}
	
	//methods
	private static int getSizefromString(String oligomer){
		int size=0;
		
		oligomer = oligomer.toLowerCase();
		
		if (oligomer.equals("monomeric")) {
		    size = 1;
		} else if (oligomer.equals("dimeric")) {
		    size = 2;
		} else if (oligomer.equals("trimeric")) {
		    size = 3;
		} else if (oligomer.equals("tetrameric")) {
		    size = 4;
		} else if (oligomer.equals("pentameric")) {
		    size = 5;
		} else if (oligomer.equals("hexameric")) {
		    size = 6;
		} else if (oligomer.equals("heptameric")) {
		    size = 7;
		} else if (oligomer.equals("octameric")) {
		    size = 8;
		} else if (oligomer.equals("nonameric")) {
		    size = 9;
		} else if (oligomer.equals("decameric")) {
		    size = 10;
		} else if (oligomer.equals("undecameric")) {
		    size = 11;
		} else if (oligomer.equals("dodecameric")) {
		    size = 12;
		} else if (oligomer.equals("tridecameric")) {
		    size = 13;
		} else if (oligomer.equals("tetradecameric")) {
		    size = 14;
		} else if (oligomer.equals("pentadecameric")) {
		    size = 15;
		} else if (oligomer.equals("hexadecameric")) {
		    size = 16;
		} else if (oligomer.equals("heptadecameric")) {
		    size = 17;
		} else if (oligomer.equals("octadecameric")) {
		    size = 18;
		} else if (oligomer.equals("nonadecameric")) {
		    size = 19;
		} else if (oligomer.equals("eicosameric")) {
		    size = 20;
		} else if( oligomer.matches("(\\d+).*")) {
		    size = Integer.parseInt((oligomer.replaceAll("(\\d+).*", "$1")));
		}
		return size;
	}

}
