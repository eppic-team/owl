package owl.core.structure;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;

import owl.core.runners.NaccessRunner;
import owl.core.structure.graphs.AICGraph;

public class ChainInterface implements Comparable<ChainInterface>, Serializable {
	
	private static final long serialVersionUID = 1L;

	private static final int FIRST = 0;
	private static final int SECOND = 1;
	
	private int id;
	private String name;
	private double interfaceArea;
	private AICGraph graph;
	
	private PdbChain firstMolecule;
	private PdbChain secondMolecule;
	
	private InterfaceRimCore firstRimCore;  // cached first molecule's rim and core
	private InterfaceRimCore secondRimCore; // cached second molecule's rim and core
	
	private double bsaToAsaCutoff;
	
	private Matrix4d firstTransf; 		// the transformation applied to first molecule expressed in crystal axes coordinates
	private Matrix4d firstTransfOrth;	// the transformation applied to first molecule expressed in orthonormal axes coordinates
	private Matrix4d secondTransf; 		// the transformation applied to second molecule expressed in crystal axes coordinates
	private Matrix4d secondTransfOrth;	// the transformation applied to second molecule expressed in orthonormal axes coordinates
	
	private Point3i transl;		// the crystal translation (excluding unit cell symmetry translations) applied to second molecule expressed in crystal axes coordinates
								// note that the translation component of secondTransf above contains this + any possible unit cell symmetry translations
	
	private double score; 			// a score value assigned to the interface (if from PISA this is the solvation energy) 
	
	private Vector3d connection; 	// the vector that connects the center of mass of first to center of mass of second
	
	private String chain2forOutput;
	
	/**
	 * Constructs an empty ChainInterface. Use the setters to set the values.
	 */
	public ChainInterface() {
		
	}
	
	public ChainInterface(PdbChain firstMolecule, PdbChain secondMolecule, AICGraph graph, Matrix4d firstTransf, Matrix4d secondTransf) {
		this.firstMolecule = firstMolecule;
		this.secondMolecule = secondMolecule;
		this.graph = graph;
		this.firstTransf = firstTransf;
		this.secondTransf = secondTransf;
		this.name = this.firstMolecule.getPdbChainCode()+"+"+this.secondMolecule.getPdbChainCode();
	}
	
	public int getId() {
		return id;
	}
	
	public void setId(int id) {
		this.id = id;
	}
	
	public String getName() {
		return name;
	}
	
	public void setName(String name) {
		this.name = name;
	}
	
	public double getInterfaceArea() {
		return interfaceArea;
	}
	
	public void setInterfaceArea(double interfaceArea) {
		this.interfaceArea = interfaceArea;
	}
	
	public PdbChain getFirstMolecule() {
		return firstMolecule;
	}
	
	public void setFirstMolecule(PdbChain firstMolecule) {
		this.firstMolecule = firstMolecule;
	}
	
	public PdbChain getSecondMolecule() {
		return secondMolecule;
	}
	
	public void setSecondMolecule(PdbChain secondMolecule) {
		this.secondMolecule = secondMolecule;
	}
	
	public double getCutoff() {
		return graph.getDistCutoff();
	}
	
	/**
	 * Calculates the AICGraph for this interface storing it locally in a cached variable.
	 * Get it subsequently with {@link #getAICGraph()}
	 * @param interfDistCutoff
	 */
	public void calcAICGraph(double interfDistCutoff) {
		this.graph = firstMolecule.getAICGraph(secondMolecule, interfDistCutoff);
	}
	
	public AICGraph getAICGraph() {
		return graph;
	}
	
	public Matrix4d getFirstTransf() {
		return firstTransf;
	}
	
	public void setFirstTransf(Matrix4d firstTransf) {
		this.firstTransf = firstTransf;
	}
	
	/**
	 * Returns the transformation used to generate the first molecule in orthonormal axes coordinates.
	 * @return
	 */
	public Matrix4d getFirstTransfOrth(){
		if (firstTransfOrth==null) {
			firstTransfOrth = firstMolecule.getParent().getCrystalCell().transfToOrthonormal(firstTransf);
		}
		return firstTransfOrth;
	}
	
	public void setFirstTransfOrth(Matrix4d firstTransfOrth) {
		this.firstTransfOrth = firstTransfOrth;
	}
	
	public Matrix4d getSecondTransf() {
		return secondTransf;
	}
	
	public void setSecondTransf(Matrix4d secondTransf) {
		this.secondTransf = secondTransf;
	} 
	
	public Matrix4d getSecondTransfOrth() {
		if (secondTransfOrth==null) {
			secondTransfOrth = secondMolecule.getParent().getCrystalCell().transfToOrthonormal(secondTransf);
		}
		return secondTransfOrth;
	}
	
	public void setSecondTransfOrth(Matrix4d secondTransfOrth) {
		this.secondTransfOrth = secondTransfOrth;
	}
	
	public Point3i getSecondTranslation() {
		return transl;
	}
	
	public void setSecondTranslation(Point3i transl) {
		this.transl = transl;
	}
	
	public double getScore() {
		return score;
	}
	
	public void setScore(double score) {
		this.score = score;
	}
	
	/**
	 * Runs the NACCESS program to calculate accessible surface areas of separate molecules
	 * and 2 molecules-complex setting both ASAs and BSAs of the residues of the two molecules making 
	 * up this interface. 
	 * The total interface area is also calculated from the individual residue values (use {@link #getInterfaceArea()} 
	 * to get it)
	 * @param naccessExecutable
	 * @param hetAtoms if true HET residues are considered, if false they aren't (NACCESS' -h option)
	 * @throws IOException
	 */
	public void calcSurfAccessNaccess(File naccessExecutable, boolean hetAtoms) throws IOException {

		// NOTE in principle it is more efficient to run naccess only once per isolated chain
		// BUT! surprisingly naccess gives slightly different values for same molecule in different 
		// orientations! (can't really understand why!)
		// That's why we run naccess always for 2 separate member of interface and the complex, otherwise 
		// we get (not very big but annoying) discrepancies and also things like negative (small) bsa values
		
		String naccessParams = "";
		if (hetAtoms) naccessParams = "-h";
		NaccessRunner nar = new NaccessRunner(naccessExecutable, naccessParams);
		
		PdbAsymUnit complex = new PdbAsymUnit();
		complex.setPdbCode(firstMolecule.getPdbCode());
		complex.setPolyChain("A", firstMolecule);
		complex.setPolyChain("B", secondMolecule);
		
		nar.runNaccess(firstMolecule);
		nar.runNaccess(secondMolecule);
		
		nar.runNaccess(complex); // this will set the bsa members of the Residues of firstMolecule and secondMolecule
	
		double totBuried = 0.0;
		for (Residue residue:firstMolecule) {
			totBuried+=residue.getBsa();
		}
		for (Residue residue:secondMolecule) {
			totBuried+=residue.getBsa();
		}

		this.setInterfaceArea(totBuried/2.0); // to use the same convention as in PISA we halve it
	}
	
	/**
	 * Calculates the buried surface area for the interface. 
	 * The two molecules must have already uncomplexed ASAs calculated, here we then
	 * calculate the ASAs for the complex and make the subtraction. Then we set the 
	 * BSAs of the residues of the two molecules making up this interface.
	 * The total interface area is also calculated from the individual residue values (use {@link #getInterfaceArea()} 
	 * to get it)
	 * @param nSpherePoints
	 * @param nThreads 
	 * @param hetAtoms if true HET residues are considered, if false they aren't, equivalent to 
	 * NACCESS' -h option
	 */
	protected void calcSurfAccess(int nSpherePoints, int nThreads, boolean hetAtoms) {

		// the ASAs of the uncomplexed chains must be already calculated by the caller
		
		PdbAsymUnit complex = new PdbAsymUnit();
		complex.setPdbCode(firstMolecule.getPdbCode());
		complex.setPolyChain("A", firstMolecule);
		complex.setPolyChain("B", secondMolecule);

		complex.calcBSAs(nSpherePoints, nThreads, hetAtoms);
		double totBuried = 0.0;
		
		for (Residue residue:firstMolecule) {
			totBuried+=residue.getBsa();
		}
		for (Residue residue:secondMolecule) {
			totBuried+=residue.getBsa();
		}

		this.setInterfaceArea(totBuried/2.0); // to use the same convention as in PISA we halve it
		
	}
	
	public int getNumAtomsInContact() {
		return this.graph.getVertexCount();
	}
	
	public int getNumContacts() {
		return this.graph.getEdgeCount();
	}
	
	public boolean isProtein() {
		return (isFirstProtein() && isSecondProtein());
	}
	
	public boolean isFirstProtein() {
		if (this.firstMolecule.isNonPolyChain()) return false;
		return this.firstMolecule.getSequence().isProtein();
	}
	
	public boolean isSecondProtein() {
		if (this.secondMolecule.isNonPolyChain()) return false;
		return this.secondMolecule.getSequence().isProtein();
	}
	
	public boolean hasClashes() {
		return this.graph.hasClashes();
	}
	
	public int getNumClashes() {
		return this.graph.getNumClashes();
	}
	
	public int compareTo(ChainInterface o) {
		// this will sort descending on interface areas
		return (Double.compare(o.interfaceArea,this.interfaceArea));
	}
	
	public String toString() {
		return firstMolecule.getPdbChainCode()+"-"+secondMolecule.getPdbChainCode()+" "+
		getAICGraph().getEdgeCount()+" "+SpaceGroup.getAlgebraicFromMatrix(secondTransf);
	}
	
	public boolean equals(Object o) {
		if (!(o instanceof ChainInterface)) return false;
		ChainInterface other = (ChainInterface) o;
		String tcc1 = firstMolecule.getPdbChainCode();
		String tcc2 = secondMolecule.getPdbChainCode();
		String occ1 = other.firstMolecule.getPdbChainCode();
		String occ2 = other.secondMolecule.getPdbChainCode();
		
		if ( !(occ1.equals(tcc1) && occ2.equals(tcc2)) &&
			 !(occ2.equals(tcc1) && occ1.equals(tcc2)) ) {
			return false;
		}
		return this.graph.equals(other.graph);
	}
	
	public int hashCode() {
		int hash = graph.getEdgeCount();
	    hash = hash * 31 + firstMolecule.getPdbChainCode().hashCode()+secondMolecule.getPdbChainCode().hashCode();
	    return hash; 
	}
	
	public void printTabular(PrintStream ps) {
		ps.print("# ");
		ps.printf("%d\t%9.2f\t%5.2f\t%s\t%s\n",this.getId(),this.getInterfaceArea(),this.getScore(), 
				this.getName(),
				SpaceGroup.getAlgebraicFromMatrix(this.getSecondTransf()));
		if (isFirstProtein()) {
			ps.print("## ");
			this.printFirstMolInfoTabular(ps);
		}
		if (isSecondProtein()) {
			ps.print("## ");
			this.printSecondMolInfoTabular(ps);
		}
	}
	
	private void printFirstMolInfoTabular(PrintStream ps) {
		String molType = null;
		if (this.getFirstMolecule().isNonPolyChain()) molType = "non-polymer";
		else if (this.getFirstMolecule().getSequence().isProtein()) molType = "protein";
		else molType = "nucleic acid";
		ps.println("1\t"+firstMolecule.getPdbChainCode()+"\t"+molType);

		InterfaceRimCore rimCore = getFirstRimCore();
		ps.printf("## %4.2f\n",bsaToAsaCutoff);
		ps.println("## rim : "+rimCore.getRimResString());
		if (rimCore.getCoreResidues().size()>0) {
			ps.println("## core: "+rimCore.getCoreResString());
		}
		ps.println("## seqres pdb res asa bsa burial(percent)");

		for (Residue residue:firstMolecule) {			
			ps.printf("%d\t%s\t%s\t%6.2f\t%6.2f",residue.getSerial(),residue.getPdbSerial(),residue.getLongCode(),residue.getAsa(),residue.getBsa());
			double percentBurial = 100.0*residue.getBsa()/residue.getAsa();
			if (percentBurial>0.1) {
				ps.printf("\t%5.1f\n",percentBurial);
			} else {
				ps.println();
			}
		}
	}

	private void printSecondMolInfoTabular(PrintStream ps) {
		String molType = null;
		if (this.getSecondMolecule().isNonPolyChain()) molType = "non-polymer";
		else if (this.getSecondMolecule().getSequence().isProtein()) molType = "protein";
		else molType = "nucleic acid";

		ps.println("2\t"+secondMolecule.getPdbChainCode()+"\t"+molType);

		InterfaceRimCore rimCore = getSecondRimCore();
		ps.printf("## %4.2f\n",bsaToAsaCutoff);
		ps.println("## rim : "+rimCore.getRimResString());
		if (rimCore.getCoreResidues().size()>0) {
			ps.println("## core: "+rimCore.getCoreResString());
		}
		ps.println("## seqres pdb res asa bsa burial(percent)");
		
		for (Residue residue:secondMolecule) {
			ps.printf("%d\t%s\t%s\t%6.2f\t%6.2f",residue.getSerial(),residue.getPdbSerial(),residue.getLongCode(),residue.getAsa(),residue.getBsa());
			double percentBurial = 100.0*residue.getBsa()/residue.getAsa();
			if (percentBurial>0.1) {
				ps.printf("\t%5.1f\n",percentBurial);
			} else {
				ps.println();
			}
		}
	}
	
	public void printRimCoreInfo(PrintStream ps) {
		

		ps.printf("%15s\t%6.1f",
				getId()+"("+getFirstMolecule().getPdbChainCode()+"+"+getSecondMolecule().getPdbChainCode()+")",
				getInterfaceArea());
		boolean isProt1 = isFirstProtein();
		boolean isProt2 = isSecondProtein();
		ps.printf("%5d\t%5d\t%5.2f", (!isProt1)?0:getFirstRimCore().getCoreSize(),
									 (!isProt2)?0:getSecondRimCore().getCoreSize(),
									 getBsaToAsaCutoff());
		ps.print("\t");

	}

	/**
	 * Calculates residues in rim and core using zooming for each of the 2 members of the 
	 * interface storing result in cached variables. 
	 * The zooming procedure is: the sum of the residues of the 2 cores is required to be at least minNumResidues, 
	 * if the minimum is not reached with the bsaToAsaSoftCutoff, then the cutoff is 
	 * relaxed in relaxationStep steps until reaching the bsaToAsaHardCutoff.
	 * If either of the 2 molecules of this interface is not a protein, its cached rimCore 
	 * object will be null.
	 * Use {@link #getFirstRimCore()} and {@link #getSecondRimCore()} to retrieve the core/rim lists.  
	 * @param bsaToAsaSoftCutoff
	 * @param bsaToAsaHardCutoff
	 * @param relaxationStep
	 * @param minNumResidues
	 * @return
	 */
	public void calcRimAndCoreZooming(double bsaToAsaSoftCutoff, double bsaToAsaHardCutoff, double relaxationStep, int minNumResidues) {
		//zoomingUsed = true;
		
		bsaToAsaCutoff = bsaToAsaHardCutoff;
		//this.bsaToAsaSoftCutoff = bsaToAsaSoftCutoff;
		//this.bsaToAsaRelaxStep = relaxationStep;
				
		if (!isFirstProtein() && !isSecondProtein()) {
			firstRimCore = null;
			secondRimCore = null;
			return;
		}
		
		// we introduce a margin of relaxationStep*0.10 to be sure we do go all the way down to bsaToAsaHardCutoff (necessary because of rounding)
		for (double cutoff=bsaToAsaSoftCutoff;cutoff>=bsaToAsaHardCutoff-relaxationStep*0.10;cutoff-=relaxationStep) {
			InterfaceRimCore rimCore1 = null;
			InterfaceRimCore rimCore2 = null;
			if (isFirstProtein()) rimCore1 = this.firstMolecule.getRimAndCore(cutoff);
			if (isSecondProtein()) rimCore2 = this.secondMolecule.getRimAndCore(cutoff);
			firstRimCore = rimCore1;
			secondRimCore = rimCore2;
			
			int totalCoreResidues = 0;
			if (isFirstProtein()) totalCoreResidues+=rimCore1.getCoreSize();
			if (isSecondProtein()) totalCoreResidues+=rimCore2.getCoreSize();
			if (totalCoreResidues>=minNumResidues) {
				bsaToAsaCutoff = cutoff;
				break;
			}
		}
		
	}

	/**
	 * Calculates residues in rim and core for given bsaToAsaCutoff, storing the 
	 * lists in cached variables.
	 * If either of the 2 molecules of this interface is not a protein, the InterfaceRimCore 
	 * object will be null for it  
	 * Use {@link #getFirstRimCore()} and {@link #getSecondRimCore()} to retrieve the core/rim lists.  
	 * @param bsaToAsaCutoff
	 * @return
	 */
	public void calcRimAndCore(double bsaToAsaCutoff) {
		//zoomingUsed = false;
		
		this.bsaToAsaCutoff = bsaToAsaCutoff;
		
		InterfaceRimCore rimCore1 = null;
		InterfaceRimCore rimCore2 = null;
		if (isFirstProtein()) { 
			rimCore1 = this.firstMolecule.getRimAndCore(bsaToAsaCutoff);
			firstRimCore = rimCore1;
		}
		
		if (isSecondProtein()) {		
			rimCore2 = this.secondMolecule.getRimAndCore(bsaToAsaCutoff);
			secondRimCore = rimCore2;
		}
	}
	
	/**
	 * Calculates core/rim lists of residues and caches them locally in {@link InterfaceRimCore} objects.
	 * Following the Chakrabarti definition (see Chakrabarti, Janin Proteins 2002)
	 * Use {@link #getFirstRimCore()} and {@link #getSecondRimCore()} to retrieve the core/rim lists. 
	 */
	public void calcRimAndCoreChakrabarti() {
		InterfaceRimCore rimCore1 = null;
		InterfaceRimCore rimCore2 = null;
		if (isFirstProtein()) { 
			rimCore1 = this.firstMolecule.getRimAndCoreChakrabarti();
			firstRimCore = rimCore1;
		}
		
		if (isSecondProtein()) {		
			rimCore2 = this.secondMolecule.getRimAndCoreChakrabarti();
			secondRimCore = rimCore2;
		}
		
	}
	
	/**
	 * Calculates core/rim lists of residues and caches them locally in {@link InterfaceRimCore} objects. 
	 * Following the Levy definition (see Levy JMB 2010). 
	 * Use {@link #getFirstRimCore()} and {@link #getSecondRimCore()} to retrieve the core/rim lists. 
	 * @param rASAcutoff 
	 */
	public void calcRimAndCoreLevy(double rASAcutoff) {
		InterfaceRimCore rimCore1 = null;
		InterfaceRimCore rimCore2 = null;
		if (isFirstProtein()) {
			rimCore1 = this.firstMolecule.getRimAndCoreLevy(rASAcutoff);
			firstRimCore = rimCore1;
		}
		
		if (isSecondProtein()) {		
			rimCore2 = this.secondMolecule.getRimAndCoreLevy(rASAcutoff);
			secondRimCore = rimCore2;
		}
		
	}
	
	public double getBsaToAsaCutoff() {
		return bsaToAsaCutoff;
	}
	
	public InterfaceRimCore getFirstRimCore() {
		return firstRimCore;
	}
	
	public InterfaceRimCore getSecondRimCore() {
		return secondRimCore;
	}
	
	public InterfaceRimCore getRimCore(int molecId) {
		if (molecId==FIRST) {
			return firstRimCore;
		} else if (molecId==SECOND) {
			return secondRimCore;
		}
		return null;
	}
	
	/**
	 * Writes this interface to given PDB file with original chain names (PDB chain codes),
	 * unless the 2 chains are the same where the second one is renamed to next letter in 
	 * alphabet.
	 * Subsequently the chain code for the renamed chain can be retrieved with {@link #getSecondPdbChainCodeForOutput()}
	 * @param file
	 * @throws FileNotFoundException 
	 */
	public void writeToPdbFile(File file) throws FileNotFoundException {
		chain2forOutput = secondMolecule.getPdbChainCode();
		if (secondMolecule.getPdbChainCode().equals(firstMolecule.getPdbChainCode())) {
			// if both chains are named equally we want to still named them differently in the output pdb file
			// so that molecular viewers can handle properly the 2 chains as separate entities 
			char letter = firstMolecule.getPdbChainCode().charAt(0);
			if (letter!='Z' && letter!='z') {
				chain2forOutput = Character.toString((char)(letter+1)); // i.e. next letter in alphabet
			} else {
				chain2forOutput = Character.toString((char)(letter-25)); //i.e. 'A' or 'a'
			}
		}
		
		PrintStream ps = new PrintStream(file);
		ps.println("HEADER");
		if (!firstMolecule.isNonPolyChain()) firstMolecule.writeSeqresRecord(ps, firstMolecule.getPdbChainCode());
		if (!secondMolecule.isNonPolyChain()) secondMolecule.writeSeqresRecord(ps,chain2forOutput);
		firstMolecule.writeAtomLines(ps, firstMolecule.getPdbChainCode());
		secondMolecule.writeAtomLines(ps,chain2forOutput);
		ps.println("END");
		ps.close();
	}
	
	public String getSecondPdbChainCodeForOutput() {
		return chain2forOutput;
	}
	
	public SubunitId getFirstSubunitId() {
		return new SubunitId(firstMolecule.getPdbChainCode().charAt(0),firstMolecule.getParent().getTransformId(),new Point3i(0,0,0));
	}

	public SubunitId getSecondSubunitId() {
		return new SubunitId(secondMolecule.getPdbChainCode().charAt(0),secondMolecule.getParent().getTransformId(),transl);
	}

	public boolean isParallel() {
		return this.getFirstSubunitId().isSymRelatedEquivalent(this.getSecondSubunitId());
	}
	
	/**
	 * Returns the vector that connects the center of mass of first 
	 * molecule to the center of mass of second molecule expressed in crystal coordinates
	 * @return
	 */
	public Vector3d getConnectionVector() {
		if (connection==null) {
			CrystalCell cell = firstMolecule.getParent().getCrystalCell();
			connection = new Vector3d(secondMolecule.getCentroid());
			connection.sub(firstMolecule.getCentroid());
			cell.transfToCrystal(connection);
		}
		return connection;
	}
	
	/**
	 * Returns true if both chains in this interface have same PDB chain code 
	 * and thus it's an interface of symmetry-related chains. Returns false otherwise.
	 * @return
	 */
	public boolean isSymRelated() {
		return (getFirstMolecule().getPdbChainCode().equals(getSecondMolecule().getPdbChainCode()));
	}
}
