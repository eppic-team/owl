/**
 * 
 */
package owl.core.structure;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import javax.vecmath.Matrix4d;

import owl.core.structure.io.BioUnitAssembly;
import owl.core.structure.io.BioUnitAssemblyGen;

/**
 * Contains a list of PdbBioUnits
 * @author biyani_n
 *
 */
public class PdbBioUnitList implements Serializable, Iterable<PdbBioUnit>{

	private static final long serialVersionUID = 1L;
	

	
	private List<PdbBioUnit> pdbBioUnits;
	
	//Default Constructor
	public PdbBioUnitList(){
		this.pdbBioUnits = new ArrayList<PdbBioUnit>();
	}
	
	/**
	 * Construct a list of pdb-bio-units from the raw data of parsers
	 * @param parent
	 * @param assemblies
	 * @param generators
	 * @param operations
	 * @param parser
	 */
	public PdbBioUnitList(PdbAsymUnit parent, ArrayList<BioUnitAssembly> assemblies, ArrayList<BioUnitAssemblyGen> generators, Map<Integer,Matrix4d> operations, String parser){
		this.pdbBioUnits = new ArrayList<PdbBioUnit>();
		boolean hasNucleotides = false;

		
		Map<Integer,Matrix4d> crystOperations = new TreeMap<Integer,Matrix4d>();		
		//Transform the orthogonal matrix to crystal coordinates
		if(parent.getCrystalCell()!=null && parent.isCrystallographicExpMethod()){
			for (int operId:operations.keySet()) {				
				Matrix4d crystOper = parent.getCrystalCell().transfToCrystal(operations.get(operId));
				crystOperations.put(operId,crystOper);
			}
		} else {
			crystOperations.putAll(operations);
		}


		//for cif parser reduce the chainCodes to pdbChainCodes
		if(parser.equals("cif")){
			for(BioUnitAssemblyGen gen:generators){
				ArrayList<String> tempList = new ArrayList<String>();
				for(String code:gen.getPdbChainCodes()){
					if(parent.getChainForChainCode(code)!=null && !(parent.getChainForChainCode(code).isNonPolyChain()))
						tempList.add(parent.getChainForChainCode(code).getPdbChainCode());
				}
				gen.setPdbChainCodes(tempList);
			}

		}

		//check if the chain code really exist in the pdbChainCodes 
		for(BioUnitAssemblyGen gen:generators){
			ArrayList<String> tempList = new ArrayList<String>();
			for(String code:gen.getPdbChainCodes()){
				if(!parent.getPdbChainCodes().contains(code)){
					//System.err.println("Warning: Chain with pdbChainCode "+code+" present in biounit-assembly but is not present in structure!");
					continue;
				}
				if(parent.getChain(code).getSequence().isProtein())
					tempList.add(code);
				else{
					//System.err.println("Warning: Chain with pdbChainCode "+code+" is not protein chain, excluding it from biounit assembly and over-writing size!");
					hasNucleotides = true;
				}
			}
			gen.setPdbChainCodes(tempList);
		}

		//Check for nulls
		if(assemblies.isEmpty()){
			//for NMR the biounit is same as the asymmetric unit
			if(parent.getExpMethod()!=null &&
					parent.getExpMethod().equals("SOLUTION NMR")){
				PdbBioUnit localUnit = new PdbBioUnit();
				
				Matrix4d idOperator = new Matrix4d();
				idOperator.m00=idOperator.m11=idOperator.m22=idOperator.m33=1.0; 
				
				int numProteinChains = 0;
				for(PdbChain chain:parent.getAllChains()) 
					if(chain.getSequence()!=null && chain.getSequence().isProtein()) {
						numProteinChains++;
						localUnit.addOperator(chain.getPdbChainCode(), idOperator);
					}
				localUnit.setSize(numProteinChains);
				localUnit.setType(BioUnitAssignmentType.getByString("authors"));	

				if(numProteinChains != 0) this.pdbBioUnits.add(localUnit);
			}
		}
		else{
			if(generators.isEmpty()) {
				//System.err.println("Warning: No bio unit assembly generator found; bio-units will not be added!");
			}
			else if(crystOperations.isEmpty()) {
				//System.err.println("Warning: No bio unit operator found; bio-units will not be added! ");
			}
			else
				for(BioUnitAssembly assembly:assemblies)
					if(assembly.getId()==0){
						//System.err.println("Warning: Error in reading Id for the assembly; will not add this biounit.");
						continue;
					} else
						for(String method:assembly.getTypes()){
							PdbBioUnit localUnit = new PdbBioUnit();
							boolean toAdd = true;

							if(assembly.getSize()>0) localUnit.setSize(assembly.getSize());
							else {
								//System.err.println("Warning: The size of biounit assembly with id:"+assembly.getId()+" seems to be 0; will not add this biounit!");
								continue;
							}

							if(BioUnitAssignmentType.getByString(method)!=null)
								localUnit.setType(BioUnitAssignmentType.getByString(method));
							else{
								//System.err.println("Warning: The assignment type of biounit assembly with id:"+assembly.getId()+" not understood; will not add this biounit!");
								continue;
							}

							for(int iGen:assembly.getGeneratorIds()){
								//Count the number of chains and number of operations
								if(generators.size() < iGen-1){ 
									//System.err.println("Warning: Assembly generator record "+iGen+" not present for assembly "+assembly.getId());
									toAdd = false;}
								else{
									BioUnitAssemblyGen gen = generators.get(iGen-1);  //ids starts from 1,2,3... and index of list starts from 0,1,2...
									int numChains = gen.getPdbChainCodes().size();
									int numOperations = gen.getOperationIds().size();

									if(numChains!=0 && numOperations!=0)
										for(String code:gen.getPdbChainCodes())
											for(int operId:gen.getOperationIds()){
												if(crystOperations.containsKey(operId)) 
													localUnit.addOperator(code, crystOperations.get(operId));
												else {
													//System.err.println("Warning: Operator record "+operId+" not present for generator "+iGen);
													toAdd = false;
												}
											}
									else {
										if(!hasNucleotides){
											//System.err.println("Warning: Either no chains or no operations in bio-unit assembly generator record "+iGen);
										}
										toAdd = false;
									}
								}
							}
							//Check if the number of operators is equal to the size of the biounit
							int numOperators = 0;
							for(String code:localUnit.getMemberPdbChainCodes())
								numOperators += localUnit.getOperators(code).size();


							if(localUnit.getSize() != numOperators) {
								if(hasNucleotides) localUnit.setSize(numOperators);
								else {
									//System.err.println("Warning: The size of the assembly for biounit "+assembly.getId()+" is not equal to it's number of operators; will not add this biounit.");
									toAdd=false;
								}
							}

							if(toAdd) this.pdbBioUnits.add(localUnit);
						}


		}

	}

	//Getters and Setters
	public List<PdbBioUnit> getPdbBioUnits() {
		return pdbBioUnits;
	}

	public void setPdbBioUnits(List<PdbBioUnit> pdbBioUnits) {
		this.pdbBioUnits = pdbBioUnits;
	}
	
	private void addPdbBioUnit(PdbBioUnit unit) {
		this.pdbBioUnits.add(unit);
		
	}
	
	public int size() {
		return pdbBioUnits.size();
	}

	public PdbBioUnitList copy() {
		PdbBioUnitList newList = new PdbBioUnitList();

		for(PdbBioUnit unit:this.pdbBioUnits)
			newList.pdbBioUnits.add(unit.copy());

		return newList;
	}

	@Override
	public Iterator<PdbBioUnit> iterator() {
		return this.pdbBioUnits.iterator();
	}

	/**
	 * Returns the PdbBioUnit with index iUnit starting from 0,1,2,...
	 * @param iUnit
	 * @return
	 */
	public PdbBioUnit get(int iUnit) {
		return this.pdbBioUnits.get(iUnit);
	}
	
	/**
	 * Returns the maximum size of the biounits in this list
	 */
	public int getMaximumSize(){
		int maxSize = 0;
		for(PdbBioUnit unit:this.pdbBioUnits)
			if(unit.getSize() > maxSize) maxSize = unit.getSize();
		return maxSize;
	}
	
	
	/**
	 * Returns a sub-list of PdbBioUnits with particular assignment type
	 */
	public PdbBioUnitList getSubsetByType(BioUnitAssignmentType type){
		PdbBioUnitList newList = new PdbBioUnitList();
		for(PdbBioUnit unit:this.pdbBioUnits)
			if(unit.getType().equals(type)) newList.addPdbBioUnit(unit);
		return newList;
	}
	
	/**
	 * Returns a sub-list of PdbBioUnits with particular size
	 */
	public PdbBioUnitList getSubsetBySize(int size){
		PdbBioUnitList newList = new PdbBioUnitList();
		for(PdbBioUnit unit:this.pdbBioUnits)
			if(unit.getSize() == size) newList.addPdbBioUnit(unit);
		return newList;
	}


}
