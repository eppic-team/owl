/**
 * 
 */
package owl.core.structure;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

import javax.vecmath.Matrix4d;

import owl.core.structure.io.BioUnitAssembly;
import owl.core.structure.io.BioUnitAssemblyGen;
import owl.core.structure.io.BioUnitOperation;

/**
 * Contains a list of PdbBioUnits
 * @author biyani_n
 *
 */
public class PdbBioUnitList implements Serializable, Iterable<PdbBioUnit>{

	private static final long serialVersionUID = 1L;
	
	/**
	 * A private class to reduce a list of PdbBioUnits to a single one 
	 * with unique type, mmSize and list of cluster ids  
	 * @author duarte_j
	 *
	 */
	private class BioUnitView {
		private int index;
		private BioUnitAssignmentType type;
		private int mmSize;
		private List<Integer> clusterIds;
		
		public BioUnitView(int index, BioUnitAssignmentType type, int mmSize, List<Integer> clusterIds) {
			this.index = index;
			this.type = type;
			this.mmSize = mmSize;
			this.clusterIds = clusterIds;
		}
		
		public boolean equals(Object o) {
			if (!(o instanceof BioUnitView)) return false;
			BioUnitView other = (BioUnitView) o;
			
			if (other.type!=this.type) return false;
			if (other.mmSize!=this.mmSize) return false;
			if (other.clusterIds.size()!=this.clusterIds.size()) return false;
			for (int i=0;i<clusterIds.size();i++) {
				if (other.clusterIds.get(i)!=this.clusterIds.get(i)) {
					return false;
				}
			}
			
			return true;
		}
	}
	
	
	private PdbAsymUnit parent;
	private List<PdbBioUnit> pdbBioUnits;
	
	//Default Constructor
	public PdbBioUnitList(PdbAsymUnit parent){
		this.parent = parent;
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
	public PdbBioUnitList(PdbAsymUnit parent, ArrayList<BioUnitAssembly> assemblies, ArrayList<BioUnitAssemblyGen> generators, ArrayList<BioUnitOperation> operations, String parser){
		this.parent = parent;
		this.pdbBioUnits = new ArrayList<PdbBioUnit>();
		boolean hasNucelotides = false;

		//Transform the orthogonal matrix to crystal coordinates
		if(parent.getCrystalCell()!=null && parent.isCrystallographicExpMethod()){
			for(BioUnitOperation oper:operations){				
				Matrix4d oper4d = new Matrix4d(oper.getOperator());
				oper4d = parent.getCrystalCell().transfToCrystal(oper4d);
				for(int i=0; i<4; i++)
					for(int j=0; j<4;j++)
						oper.setOperatorValue(4*i+j, oper4d.getElement(i, j));
			}
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
					hasNucelotides = true;
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
				localUnit.assignType(BioUnitAssignmentType.getByString("authors"));	

				if(numProteinChains != 0) this.pdbBioUnits.add(localUnit);
			}
		}
		else{
			if(generators.isEmpty()) {
				//System.err.println("Warning: No bio unit assembly generator found; bio-units will not be added!");
			}
			else if(operations.isEmpty()) {
				//System.err.println("Warning: No bio unit operator found; bio-units will not be added! ");
			}
			else
				for(BioUnitAssembly assembly:assemblies)
					if(assembly.getId()==0){
						//System.err.println("Warning: Error in reading Id for the assembly; will not add this biounit.");
						continue;
					}else
						for(String method:assembly.getTypes()){
							PdbBioUnit localUnit = new PdbBioUnit();
							boolean toAdd = true;

							if(assembly.getSize()>0) localUnit.setSize(assembly.getSize());
							else {
								//System.err.println("Warning: The size of biounit assembly with id:"+assembly.getId()+" seems to be 0; will not add this biounit!");
								continue;
							}

							if(BioUnitAssignmentType.getByString(method)!=BioUnitAssignmentType.none)
								localUnit.assignType(BioUnitAssignmentType.getByString(method));
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
												if(getOperator(operId, operations)!=null) localUnit.addOperator(code, getOperator(operId, operations));
												else {
													//System.err.println("Warning: Operator record "+operId+" not present for generator "+iGen);
													toAdd = false;
												}
											}
									else {
										if(!hasNucelotides){
											//System.err.println("Warning: Either no chains or no operations in bio-unit assembly generator record "+iGen);
										}
										toAdd = false;
									}
								}
							}
							//Check if the number of operators is equal to the size of the biounit
							int numOperators = 0;
							for(String code:localUnit.getOperators().keySet())
								numOperators += localUnit.getOperators().get(code).size();


							if(localUnit.getSize() != numOperators){
								if(hasNucelotides) localUnit.setSize(numOperators);
								else{
									//System.err.println("Warning: The size of the assembly for biounit "+assembly.getId()+" is not equal to it's number of operators; will not add this biounit.");
									toAdd=false;
									}
							}

							if(toAdd) this.pdbBioUnits.add(localUnit);
						}


		}

	}

	//private method
	private static Matrix4d getOperator(int id, ArrayList<BioUnitOperation> operations){
		Matrix4d operation = new Matrix4d();
		boolean present = false;
		for(BioUnitOperation oper:operations){
			if(id == oper.getId()) {
				operation = new Matrix4d(oper.getOperator());
				present = true;
				break;
			}
		}

		if(present) return operation;
		else return null;
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

	public PdbAsymUnit getParent() {
		return parent;
	}

	public void setParent(PdbAsymUnit parent) {
		this.parent = parent;
	}

	public PdbBioUnitList copy(PdbAsymUnit parent) {
		PdbBioUnitList newList = new PdbBioUnitList(parent);

		for(PdbBioUnit unit:this.pdbBioUnits)
			newList.pdbBioUnits.add(unit.copy());

		return newList;
	}

	@Override
	public Iterator<PdbBioUnit> iterator() {
		return this.pdbBioUnits.iterator();
	}
	
	/**
	 * For each biounit in the list returns the list of ids of interfaces that match that biounit.
	 * In case of no matches found for a certain biounit, the biounit will not be present in the returned map
	 * Ids for BioUnits start from 0,1,2,...
	 * Ids for Interfaces start from 1,2,3,...
	 * @param interfaces
	 * @return
	 */
	public TreeMap<Integer, List<Integer>> getInterfaceMatches(ChainInterfaceList interfaces){
		TreeMap<Integer,List<Integer>> matches = new TreeMap<Integer, List<Integer>>();
		for(PdbBioUnit bioUnit:this.pdbBioUnits){
			List<Integer> matchList = bioUnit.getInterfaceMatches(interfaces);
			// if nothing matches for this bioUnit, don't add
			if(bioUnit.getSize() > 1 && matchList.size() < 1) continue;
			else matches.put(this.pdbBioUnits.indexOf(bioUnit), matchList);
		}
		return matches;
	}
	
	/**
	 * Retunrs a map of PdbBioUnit indices to list of interface cluster ids that match
	 * that biounit assembly. Removes duplicate PdbBioUnits by grouping those with same
	 * type, size and list of cluster ids and considering them identical. 
	 * @param interfaces
	 * @return
	 */
	public TreeMap<Integer, List<Integer>> getInterfaceClusterMatches(ChainInterfaceList interfaces) {
		
		List<BioUnitView> buList = new ArrayList<BioUnitView>();
		
		for(PdbBioUnit bioUnit:this.pdbBioUnits){
			
			List<Integer> matchList = bioUnit.getInterfaceClusterMatches(interfaces);
			// if nothing matches for this bioUnit, don't add
			if(bioUnit.getSize() > 1 && matchList.size() < 1) {
				continue;
			}
			else {
				BioUnitView bu = new BioUnitView(pdbBioUnits.indexOf(bioUnit), bioUnit.getType(),bioUnit.getSize(),matchList);
				if (!buList.contains(bu)) {
					buList.add(bu);
				}
			}
		}
		TreeMap<Integer, List<Integer>> map = new TreeMap<Integer,List<Integer>>();
		
		for (BioUnitView bu:buList) {
			map.put(bu.index, bu.clusterIds);
		}
		
		return map;		
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
		PdbBioUnitList newList = new PdbBioUnitList(this.parent);
		for(PdbBioUnit unit:this.pdbBioUnits)
			if(unit.getType().equals(type)) newList.addPdbBioUnit(unit);
		return newList;
	}
	
	/**
	 * Returns a sub-list of PdbBioUnits with particular size
	 */
	public PdbBioUnitList getSubsetBySize(int size){
		PdbBioUnitList newList = new PdbBioUnitList(this.parent);
		for(PdbBioUnit unit:this.pdbBioUnits)
			if(unit.getSize() == size) newList.addPdbBioUnit(unit);
		return newList;
	}

	/**
	 * Returns a map with Assignment Type as the keys and corresponding array of gathered results.
	 * The array has the size of the number of interfaces and returns 
	 * 0: Match
	 * 1: Non-match
	 * -1: No results
	 * Rules:
	 * If present, the highest possible assembly of an assignment type,
	 * else take the union of all the given highest assemblies of the assignment type
	 * @param interfaces
	 * @return
	 */
	public TreeMap<BioUnitAssignmentType, int[]> gatherBioUnits(ChainInterfaceList interfaces) {
		
		TreeMap<Integer, List<Integer>> matchIds = getInterfaceMatches(interfaces);
		TreeMap<BioUnitAssignmentType, int[]> summary = new TreeMap<BioUnitAssignmentType, int[]>();

		PdbBioUnitList authorsList = this.getSubsetByType(BioUnitAssignmentType.authors);
		PdbBioUnitList pisaList = this.getSubsetByType(BioUnitAssignmentType.pisa);
		PdbBioUnitList pqsList = this.getSubsetByType(BioUnitAssignmentType.pqs);
		
		int[] authorFinal = new int[interfaces.getNumInterfaces()];
		int[] pisaFinal = new int[interfaces.getNumInterfaces()];
		int[] pqsFinal = new int[interfaces.getNumInterfaces()];
		
		for(int i=1; i<=interfaces.getNumInterfaces(); i++){
			
			authorFinal[i-1]=pisaFinal[i-1]=pqsFinal[i-1]=-1;
			
			//Check for authors
			for(PdbBioUnit bioUnit:authorsList.getSubsetBySize(authorsList.getMaximumSize())){
				int iUnit = this.pdbBioUnits.indexOf(bioUnit);
				if(matchIds.containsKey(iUnit)){
					if(matchIds.get(iUnit).contains(i)){
						authorFinal[i-1] = 0;
						break;
					} else authorFinal[i-1] = 1;
				}
			}

			//check for Pisa
			for(PdbBioUnit bioUnit:pisaList.getSubsetBySize(pisaList.getMaximumSize())){
				int iUnit = this.pdbBioUnits.indexOf(bioUnit);
				if(matchIds.containsKey(iUnit)){
					if(matchIds.get(iUnit).contains(i)){
						pisaFinal[i-1] = 0;
						break;
					} else pisaFinal[i-1] = 1;
				}
			}

			//check for pqs
			for(PdbBioUnit bioUnit:pqsList.getSubsetBySize(pqsList.getMaximumSize())){
				int iUnit = this.pdbBioUnits.indexOf(bioUnit);
				if(matchIds.containsKey(iUnit)){
					if(matchIds.get(iUnit).contains(i)){
						pqsFinal[i-1] = 0;
						break;
					} else pqsFinal[i-1] = 1;
				}
			}
		}

		summary.put(BioUnitAssignmentType.authors, authorFinal);
		summary.put(BioUnitAssignmentType.pisa, pisaFinal);
		summary.put(BioUnitAssignmentType.pqs, pqsFinal);

		return summary;
	}

}
