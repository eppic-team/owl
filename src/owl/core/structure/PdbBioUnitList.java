/**
 * 
 */
package owl.core.structure;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

import owl.core.structure.io.BioUnitAssembly;
import owl.core.structure.io.BioUnitAssemblyGen;
import owl.core.structure.io.BioUnitOperation;

/**
 * Contains a list of PdbBioUnits
 * @author biyani_n
 *
 */
public class PdbBioUnitList implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
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

		//Transform the translation matrix to crystal coordinates
		if(parent.getCrystalCell()!=null){
			for(BioUnitOperation oper:operations){
				Vector3d v=new Vector3d(oper.getOperator()[3],oper.getOperator()[7],oper.getOperator()[11]);
				parent.getCrystalCell().transfToCrystal(v);

				oper.setOperatorValue(3, v.x);
				oper.setOperatorValue(7, v.y);
				oper.setOperatorValue(11, v.z);
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
			for(String code:gen.getPdbChainCodes())
				if(parent.getPdbChainCodes().contains(code))
					tempList.add(code);
				else System.err.println("Warning: Chain with pdbChainCode "+code+" present in biounit-assembly but is not present in structure!");
			gen.setPdbChainCodes(tempList);
		}

		//Check for nulls
		if(assemblies.isEmpty()){
			//for NMR the biounit is same as the asymmetric unit
			if(parent.getExpMethod()!=null &&
					parent.getExpMethod().equals("SOLUTION NMR")){
				PdbBioUnit localUnit = new PdbBioUnit();
				localUnit.setSize(parent.getNumPolyChains());
				localUnit.assignType(BioUnitAssignmentType.getByString("authors"));

				Matrix4d idOperator = new Matrix4d();
				idOperator.m00=idOperator.m11=idOperator.m22=idOperator.m33=1.0; 
				for(String pdbChainCode:parent.getPdbChainCodes()) localUnit.addOperator(pdbChainCode, idOperator);

				this.pdbBioUnits.add(localUnit);
			}
		}
		else{
			if(generators.isEmpty()) System.err.println("Warning: No bio unit assembly generator found; bio-units will not be added!");
			else if(operations.isEmpty()) System.err.println("Warning: No bio unit operator found; bio-units will not be added! ");
			else
				for(BioUnitAssembly assembly:assemblies)
					if(assembly.getId()==0){
						System.err.println("Warning: Error in reading Id for the assembly; will not add this biounit.");
						continue;
					}else
						for(String method:assembly.getTypes()){
							PdbBioUnit localUnit = new PdbBioUnit();
							boolean toAdd = true;

							if(assembly.getSize()>0) localUnit.setSize(assembly.getSize());
							else {
								System.err.println("Warning: The size of biounit assembly with id:"+assembly.getId()+" seems to be 0; will not add this biounit!");
								continue;
							}

							if(BioUnitAssignmentType.getByString(method)!=BioUnitAssignmentType.none)
								localUnit.assignType(BioUnitAssignmentType.getByString(method));
							else{
								System.err.println("Warning: The assignment type of biounit assembly with id:"+assembly.getId()+" not understood; will not add this biounit!");
								continue;
							}

							for(int iGen:assembly.getGeneratorIds()){
								//Count the number of chains and number of operations
								if(generators.size() < iGen-1){ 
									System.err.println("Warning: Assembly generator record "+iGen+" not present for assembly "+assembly.getId());
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
													System.err.println("Warning: Operator record "+operId+" not present for generator "+iGen);
													toAdd = false;
												}
											}
									else {
										System.err.println("Warning: Either no chains or no operations in bio-unit assembly generator record "+iGen);
										toAdd = false;
									}
								}
							}
							//Check if the number of operators is equal to the size of the biounit
							int numOperators = 0;
							for(String code:localUnit.getOperators().keySet())
								numOperators += localUnit.getOperators().get(code).size();


							if(localUnit.getSize() != numOperators){
								System.err.println("Warning: The size of the assembly for biounit "+assembly.getId()+" is not equal to it's number of operators; will not add this biounit.");
								toAdd=false;
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

}
