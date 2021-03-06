/**
 * 
 */
package owl.core.structure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import javax.vecmath.Matrix4d;

import owl.core.util.FileFormatException;
import owl.core.util.Goodies;

/**
 * A biounit annotation from PDB.
 * 
 * Note that the class implements the comparable interface so that 
 * biounits can be sorted and compared in parsing tests.
 * 
 * @author biyani_n
 *
 */
public class PdbBioUnit implements Serializable, Comparable<PdbBioUnit> { 
	
	private static final long serialVersionUID = 1L;
	
	private int size;										//size of the biounit
	private BioUnitAssignmentType type;						//assignment type: authors/pisa/pqs/eppic/none
	private TreeMap<String, List<Matrix4d>> operators;		//Map of PDB chaincodes to the list of operators
	
	//Methods
	
	//Main constructor
	public PdbBioUnit(){
		this.size = 0;
		this.type = null;
		this.operators = new TreeMap<String, List<Matrix4d>>();
	}
	
	/**
	 * Adds an entry to operators
	 * @param pdbChainCode
	 * @param operator
	 */
	public void addOperator(String pdbChainCode, Matrix4d operator){
		List<Matrix4d> tempList = new ArrayList<Matrix4d>();
		if(this.operators.containsKey(pdbChainCode)) tempList = this.operators.get(pdbChainCode);

		tempList.add(operator);
		this.operators.put(pdbChainCode, tempList);
	}
	
	/**
	 * Sets the size of biounit
	 * @param size
	 */
	public void setSize(int size){
		this.size=size;
	}
	
	/**
	 * Assigns the type of biounit
	 * @param type
	 */
	public void setType(BioUnitAssignmentType type){
		this.type=type;
	}
	
	/**
	 * Gets the size of the biounit
	 * @return
	 */
	public int getSize(){
		return this.size;
	}
	
	/**
	 * Gets the assignment type
	 * @return
	 */
	public BioUnitAssignmentType getType(){
		return this.type;
	}
	
	private TreeMap<String, List<Matrix4d>> getOperators(){
		return this.operators;
	}
	
	public Set<String> getMemberPdbChainCodes() {
		return this.operators.keySet();
	}
	
	public List<Matrix4d> getOperators(String pdbChainCode) {
		return this.operators.get(pdbChainCode);
	}
	
	//comparing methods
	public boolean equals(Object ob){
		if(ob == null) return false;
		if(!(ob instanceof PdbBioUnit)) return false;
		PdbBioUnit o = (PdbBioUnit) ob;
		
		if(this.size != o.getSize()) return false;
		if (this.type != o.getType()) return false;
		if (this.operators.size() != o.getOperators().size()) return false;
		if (!this.operators.keySet().equals(o.operators.keySet())) return false;
		
		for (String thisKey:this.operators.keySet()) {
			List<Matrix4d> thisList = this.operators.get(thisKey);
			List<Matrix4d> otherList = o.operators.get(thisKey);
			if (thisList.size()!=otherList.size()) return false;
			
			for(int j=0; j<thisList.size(); j++){
				// this assumes same order in both lists
				if (!thisList.get(j).epsilonEquals(otherList.get(j), 0.001)) 
					return false;
			}
			
		}
		return true;
	}
	
	@Override
	public int compareTo(PdbBioUnit o) {
		int ret=0;
		if(this.size > o.size){
			ret = 1;
		} else if(this.size < o.size){
			ret = -1;
		} else if(this.type.getId() > o.type.getId()){
			ret = 1;
		} else if(this.type.getId() < o.type.getId()){
			ret = -1;
		} else if(this.operators.size() > o.operators.size()){
			ret = 1;
		} else if(this.operators.size() < o.operators.size()){
			ret = -1;
		} else{
			String[] thisKeys = new String[this.operators.size()]; 
			this.operators.keySet().toArray(thisKeys);
			String[] oKeys = new String[o.operators.size()];
			o.operators.keySet().toArray(oKeys);
			for(int i=0; i < this.operators.keySet().size(); i++){
				int cmpStr = thisKeys[i].compareTo(oKeys[i]);
				if(cmpStr!=0){
					ret = cmpStr; break;
				}else if(this.operators.get(thisKeys[i]).size() > o.operators.get(oKeys[i]).size()){
					ret = 1; break;
				}else if(this.operators.get(thisKeys[i]).size() < o.operators.get(oKeys[i]).size()){
					ret = -1; break;
				}else{
					for(int j=0; j<this.operators.get(thisKeys[i]).size(); j++){
						for(int ii=0; ii<4; ii++){
							for(int jj=0; jj<4; jj++){
								Double f = this.operators.get(thisKeys[i]).get(j).getElement(ii, jj);
								Double s = o.operators.get(oKeys[i]).get(j).getElement(ii, jj);
								int cmpDble = f.compareTo(s);
								if(cmpDble!=0) {ret = cmpDble; break;}
								else continue;
							}
						}
					}
				}
			}	
		}

		return ret;
		
	}
	
	public PdbBioUnit copy() {
		PdbBioUnit unit = new PdbBioUnit();
		unit.size = this.size;
		unit.type = this.type;
		
		for(String code:this.operators.keySet())
			for(Matrix4d oper:this.operators.get(code))
				unit.addOperator(code, new Matrix4d(oper));
		
		
		return unit;
	}
	
	/**
	 * Test method
	 * @param args
	 * @throws IOException 
	 * @throws PdbLoadException 
	 * @throws FileFormatException 
	 */
	public static void main(String[] args) throws IOException, FileFormatException, PdbLoadException {
		File inputFile = new File(args[0]);
		@SuppressWarnings("resource")
		BufferedReader br = new BufferedReader(new FileReader(inputFile));
		String line;
		ArrayList<String> pdbCodes = new ArrayList<String>();
		while ((line=br.readLine())!=null){
			if (line.startsWith("#")) { continue; }
			for(String code:line.split(" ")) pdbCodes.add(code);
		}
		
		System.out.println("-----CIF FILE PARSER-------");
		for (String pdbCode:pdbCodes){
			try{
				pdbCode = pdbCode.toLowerCase();
				System.out.println("\nPdb: "+pdbCode);
				
				File repoGzFile = new File("/nfs/data/dbs/pdb/data/structures/all/mmCIF/",pdbCode+".cif.gz");
				String prefix = repoGzFile.getName().substring(0,repoGzFile.getName().lastIndexOf(".gz"));
				File pdbFile = File.createTempFile(prefix,"");
				pdbFile.deleteOnExit();
				Goodies.gunzipFile(repoGzFile, pdbFile);
				PdbAsymUnit asymUnit = new PdbAsymUnit(pdbFile);

				//Collections.sort(asymUnit.getPdbBioUnits());

				for (PdbBioUnit biounit:asymUnit.getPdbBioUnitList().getPdbBioUnits() ){
					System.out.println("BioUnit: size="+biounit.getSize()+" type:"+biounit.getType().getType());
					for(String chain:biounit.getOperators().keySet())
						System.out.println("Chain: "+chain+" has operators:\n"+biounit.getOperators().get(chain));
				}
			}catch(Throwable t)
			{
				t.printStackTrace();
			}
		}

		System.out.println("\n-----PDB FILE PARSER-------");
		for (String pdbCode:pdbCodes){
			try{
				pdbCode = pdbCode.toLowerCase();
				System.out.println("\nPdb: "+pdbCode);
				File repoGzFile = new File("/nfs/data/dbs/pdb/data/structures/all/pdb/","pdb"+pdbCode+".ent.gz");
				String prefix = repoGzFile.getName().substring(0,repoGzFile.getName().lastIndexOf(".gz"));
				File pdbFile = File.createTempFile(prefix,"");
				pdbFile.deleteOnExit();
				Goodies.gunzipFile(repoGzFile, pdbFile);
				PdbAsymUnit asymUnit = new PdbAsymUnit(pdbFile);

				//Collections.sort(asymUnit.getPdbBioUnits());

				for (PdbBioUnit biounit:asymUnit.getPdbBioUnitList().getPdbBioUnits() ){
					System.out.println("BioUnit: size="+biounit.getSize()+" type:"+biounit.getType().getType());
					for(String chain:biounit.getOperators().keySet())
						System.out.println("Chain: "+chain+" has operators:\n"+biounit.getOperators().get(chain));
				}
			}catch(Throwable t){
				t.printStackTrace();
			}
		}

	}

	
}
