package owl.core.structure;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import javax.vecmath.Matrix4d;
import javax.vecmath.Tuple3d;

import owl.core.runners.NaccessRunner;
import owl.core.structure.graphs.AICGraph;
import owl.core.util.GeometryTools;
import owl.core.util.OptSuperposition;

public class ChainInterface implements Comparable<ChainInterface>, Serializable {
	
	private static final long serialVersionUID = 1L;

	public static final int FIRST = 0;
	public static final int SECOND = 1;
	
	private int id;
	private String name;
	private double interfaceArea;
	private AICGraph graph;
	
	private PdbChain firstMolecule;
	private PdbChain secondMolecule;
	
	private List<PdbChain> firstCofactors;
	private List<PdbChain> secondCofactors;
	
	private InterfaceRimCore firstRimCore;  // cached first molecule's rim and core
	private InterfaceRimCore secondRimCore; // cached second molecule's rim and core
	
	private double bsaToAsaCutoff;
	
	private CrystalTransform firstTransf; 		// the transformation applied to first molecule
	private Matrix4d firstTransfOrth;	// the transformation applied to first molecule expressed in orthonormal axes coordinates
	private CrystalTransform secondTransf; 		// the transformation applied to second molecule
	private Matrix4d secondTransfOrth;	// the transformation applied to second molecule expressed in orthonormal axes coordinates
	
	private String chain2forOutput;
	
	/**
	 * Constructs an empty ChainInterface. Use the setters to set the values.
	 */
	public ChainInterface() {
		
	}
	
	public ChainInterface(PdbChain firstMolecule, PdbChain secondMolecule, AICGraph graph, CrystalTransform firstTransf, CrystalTransform secondTransf) {
		this.firstMolecule = firstMolecule;
		this.secondMolecule = secondMolecule;
		this.graph = graph;
		this.firstTransf = firstTransf;
		this.secondTransf = secondTransf;
		this.name = this.firstMolecule.getPdbChainCode()+"+"+this.secondMolecule.getPdbChainCode();
	}
	
	public List<PdbChain> getCofactors(int molecId) {
		if (molecId==FIRST) return getFirstCofactors();
		if (molecId==SECOND) return getSecondCofactors();
		return null;
	}
	
	public void setFirstCofactors(List<PdbChain> cofactors) {
		this.firstCofactors = cofactors;
	}
	
	public List<PdbChain> getFirstCofactors() {
		return firstCofactors;
	}
	
	public void setSecondCofactors(List<PdbChain> cofactors) {
		this.secondCofactors = cofactors;
	}
	
	public List<PdbChain> getSecondCofactors() {
		return secondCofactors;
	}
	
	public boolean hasCofactors() {
		if (firstCofactors==null && secondCofactors==null) return false;
		if (firstCofactors!=null && firstCofactors.size()>0) return true;
		if (secondCofactors!=null && secondCofactors.size()>0) return true;
		return false;
	}
	
	public boolean hasFirstCofactors() {
		if (firstCofactors!=null && firstCofactors.size()>0) return true;
		return false;
	}
	
	public boolean hasSecondCofactors() {
		if (secondCofactors!=null && secondCofactors.size()>0) return true;
		return false;
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
	
	public PdbChain getMolecule(int molecId) {
		if (molecId==FIRST) return getFirstMolecule();
		if (molecId==SECOND) return getSecondMolecule();
		return null;
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
	
	public CrystalTransform getFirstTransf() {
		return firstTransf;
	}
	
	public void setFirstTransf(CrystalTransform firstTransf) {
		this.firstTransf = firstTransf;
	}
	
	/**
	 * Returns the transformation used to generate the first molecule in orthonormal axes coordinates.
	 * @return
	 */
	public Matrix4d getFirstTransfOrth(){
		if (firstTransfOrth==null) {
			firstTransfOrth = firstMolecule.getParent().getCrystalCell().transfToOrthonormal(firstTransf.getMatTransform());
		}
		return firstTransfOrth;
	}
	
	public void setFirstTransfOrth(Matrix4d firstTransfOrth) {
		this.firstTransfOrth = firstTransfOrth;
	}
	
	public CrystalTransform getSecondTransf() {
		return secondTransf;
	}
	
	public void setSecondTransf(CrystalTransform secondTransf) {
		this.secondTransf = secondTransf;
	} 
	
	public Matrix4d getSecondTransfOrth() {
		if (secondTransfOrth==null) {
			secondTransfOrth = secondMolecule.getParent().getCrystalCell().transfToOrthonormal(secondTransf.getMatTransform());
		}
		return secondTransfOrth;
	}
	
	public void setSecondTransfOrth(Matrix4d secondTransfOrth) {
		this.secondTransfOrth = secondTransfOrth;
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
		
		// 1) we find the total number of atoms we have (two molecules plus possible cofactors)
		int numAtoms = 0;
		for (int molecId=0;molecId<2;molecId++) {
			PdbChain chain = getMolecule(molecId);
			chain.setAtomRadii();
			if (hetAtoms) {
				numAtoms += chain.getNumAtoms();
			} else {
				numAtoms += chain.getNumNonHetAtoms();
			}
			for (PdbChain cofactor:getCofactors(molecId)) {
				cofactor.setAtomRadii();
				numAtoms += cofactor.getNumAtoms();
			}
		}

		// 2) then we fill the atoms array with atoms in 2 molecules and possible cofactors
		Atom[] atoms = new Atom[numAtoms];
		
		int i = 0;
		for (int molecId=0;molecId<2;molecId++) {
			PdbChain chain = getMolecule(molecId);
			for (Residue residue:chain) {
				if (!hetAtoms && (residue instanceof HetResidue)) continue;
				for (Atom atom:residue) {
					atoms[i] = atom;
					i++;
				}
			}
			for (PdbChain cofactor:getCofactors(molecId)) {
				for (Residue residue:cofactor) {
					for (Atom atom:residue) {
						atoms[i] = atom;
						i++;
					}
				}
			}
		}
		
		// 3) we calculate asas for the complex
		double[] asas = Asa.calculateAsa(atoms, Asa.DEFAULT_PROBE_SIZE, nSpherePoints, nThreads);
		
		// 4) by subtraction to the uncomplex values (that should be present in atoms from previously calling calcASAs) we get bsas 
		for (i=0;i<atoms.length;i++){
			atoms[i].setBsa(atoms[i].getAsa()-asas[i]);
		}		
		// and the sums per residue
		for (int molecId=0;molecId<2;molecId++) {
			PdbChain chain = getMolecule(molecId);
			for (Residue residue:chain) {
				double tot = 0;
				for (Atom atom:residue) {
					tot+=atom.getBsa();
				}
				residue.setBsa(tot);
			}
			// even though we won't use the cofactors' bsas it doesn't hurt to calculate them
			for (PdbChain cofactor:getCofactors(molecId)) {
				for (Residue residue:cofactor) {
					double tot = 0;
					for (Atom atom:residue) {
						tot+=atom.getBsa();
					}				
					residue.setBsa(tot);
				}
			}
		}

		// 5) we finally calculate the interface bsa from summing all residues' bsas
		//    note here we count only the polymer molecules and not the cofactors
		double totBuried = 0.0;
		for (int molecId=0;molecId<2;molecId++) {
			for (Residue residue:getMolecule(molecId)) {
				totBuried+=residue.getBsa();
			}
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
		getAICGraph().getEdgeCount()+" "+SpaceGroup.getAlgebraicFromMatrix(secondTransf.getMatTransform());
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
	
	public void printTabular(PrintStream ps, boolean usePdbResSer) {
		ps.print("# ");
		ps.printf("%d\t%9.2f\t%s\t%s\n",this.getId(),this.getInterfaceArea(), 
				this.getName(),
				SpaceGroup.getAlgebraicFromMatrix(this.getSecondTransf().getMatTransform()));
		if (isFirstProtein()) {
			ps.print("## ");
			this.printMolInfoTabular(ps, FIRST, usePdbResSer);
		}
		if (isSecondProtein()) {
			ps.print("## ");
			this.printMolInfoTabular(ps, SECOND, usePdbResSer);
		}
	}
	
	private void printMolInfoTabular(PrintStream ps, int molecId, boolean usePdbResSer) {
		PdbChain molecule = null;
		InterfaceRimCore rimCore = null;
		if (molecId==FIRST) {
			molecule = this.getFirstMolecule();
			rimCore = getFirstRimCore();
		}
		else if (molecId==SECOND) {
			molecule = this.getSecondMolecule();
			rimCore = getSecondRimCore();
		}
		else {
			throw new IllegalArgumentException("Molecule id "+molecId+" is not valid");
		}
		
		String molType = null;
		if (molecule.isNonPolyChain()) molType = "non-polymer";
		else if (molecule.getSequence().isProtein()) molType = "protein";		
		else molType = "nucleic acid";
		
		ps.println((molecId+1)+"\t"+molecule.getPdbChainCode()+"\t"+molType);

		int numSurfRes0 = molecule.getSurfaceResidues(0).size();
		int numSurfResCutoff = molecule.getSurfaceResidues(rimCore.getMinAsaForSurface()).size();
		ps.printf("## surface residues: %d (min ASA for surface 0), %d (min ASA for surface %2.0f)\n",
				numSurfRes0,numSurfResCutoff,rimCore.getMinAsaForSurface());
		
		ps.printf("## %4.2f\n",bsaToAsaCutoff);
		ps.println("## rim : "+rimCore.getRimResString(usePdbResSer));
		if (rimCore.getCoreResidues().size()>0) {
			ps.println("## core: "+rimCore.getCoreResString(usePdbResSer));
		}
		ps.println("## seqres pdb res asa bsa burial(percent)");

		for (Residue residue:molecule) {			
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
	 * @param minAsaForSurface
	 * @return
	 */
	public void calcRimAndCoreZooming(double bsaToAsaSoftCutoff, double bsaToAsaHardCutoff, double relaxationStep, int minNumResidues, double minAsaForSurface) {
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
			if (isFirstProtein()) rimCore1 = this.firstMolecule.getRimAndCore(cutoff,minAsaForSurface);
			if (isSecondProtein()) rimCore2 = this.secondMolecule.getRimAndCore(cutoff,minAsaForSurface);
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
	 * @param minAsaForSurface
	 * @return
	 */
	public void calcRimAndCore(double bsaToAsaCutoff, double minAsaForSurface) {
		//zoomingUsed = false;
		
		this.bsaToAsaCutoff = bsaToAsaCutoff;
		
		InterfaceRimCore rimCore1 = null;
		InterfaceRimCore rimCore2 = null;
		if (isFirstProtein()) { 
			rimCore1 = this.firstMolecule.getRimAndCore(bsaToAsaCutoff, minAsaForSurface);
			firstRimCore = rimCore1;
		}
		
		if (isSecondProtein()) {		
			rimCore2 = this.secondMolecule.getRimAndCore(bsaToAsaCutoff, minAsaForSurface);
			secondRimCore = rimCore2;
		}
	}
	
	/**
	 * Calculates core/rim lists of residues and caches them locally in {@link InterfaceRimCore} objects.
	 * Following the Chakrabarti definition (see Chakrabarti, Janin Proteins 2002)
	 * Use {@link #getFirstRimCore()} and {@link #getSecondRimCore()} to retrieve the core/rim lists.
	 * @param minAsaForSurface 
	 */
	public void calcRimAndCoreChakrabarti(double minAsaForSurface) {
		InterfaceRimCore rimCore1 = null;
		InterfaceRimCore rimCore2 = null;
		if (isFirstProtein()) { 
			rimCore1 = this.firstMolecule.getRimAndCoreChakrabarti(minAsaForSurface);
			firstRimCore = rimCore1;
		}
		
		if (isSecondProtein()) {		
			rimCore2 = this.secondMolecule.getRimAndCoreChakrabarti(minAsaForSurface);
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
	 * alphabet. Subsequently the chain code for the renamed chain can be retrieved 
	 * with {@link #getSecondPdbChainCodeForOutput()} 
	 * @param file
	 * @param usePdbResSer if true PDB residue serials are written, if false CIF residue 
	 * serials are written
	 * @param gzip if true file will be gzipped, false it will be plain text
	 * @throws IOException 
	 */
	public void writeToPdbFile(File file, boolean usePdbResSer, boolean gzip) throws IOException {
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
			
		PrintStream ps = null;
		if (gzip) {
			ps = new PrintStream(new GZIPOutputStream(new FileOutputStream(file)));
		} else {
			ps = new PrintStream(file);
		}
		ps.println("HEADER");
		if (!firstMolecule.isNonPolyChain()) firstMolecule.writeSeqresRecord(ps, firstMolecule.getPdbChainCode());
		if (!secondMolecule.isNonPolyChain()) secondMolecule.writeSeqresRecord(ps,chain2forOutput);
		firstMolecule.writeAtomLines(ps, firstMolecule.getPdbChainCode(), usePdbResSer);
		if (firstCofactors!=null) {
			for (PdbChain cofactor:firstCofactors) {
				cofactor.writeAtomLines(ps, firstMolecule.getPdbChainCode(), usePdbResSer);
			}
		}
		secondMolecule.writeAtomLines(ps,chain2forOutput, usePdbResSer);
		if (secondCofactors!=null) {
			for (PdbChain cofactor:secondCofactors) {
				cofactor.writeAtomLines(ps, chain2forOutput, usePdbResSer);
			}
		}
		ps.println("END");
		ps.close();
	}
	
	public String getSecondPdbChainCodeForOutput() {
		return chain2forOutput;
	}
	
	public SubunitId getFirstSubunitId() {
		return new SubunitId(firstMolecule.getPdbChainCode().charAt(0),firstTransf);
	}

	public SubunitId getSecondSubunitId() {
		return new SubunitId(secondMolecule.getPdbChainCode().charAt(0),secondTransf);
	}
	
	/**
	 * Returns true if the transformation applied to the second chain of this interface
	 * has an infinite character (pure translation or screw rotation)
	 * and both chains of the interface have the same PDB chain code: in such cases the 
	 * interface would lead to infinite fiber-like (linear or helical) assemblies
	 * @return
	 */
	public boolean isInfinite() {
		return ((isSymRelated() && getSecondTransf().getTransformType().isInfinite()));
	}

	/**
	 * Returns true if both chains in this interface have same PDB chain code 
	 * and thus it's an interface of symmetry-related chains. Returns false otherwise.
	 * @return
	 */
	public boolean isSymRelated() {
		return (getFirstMolecule().getPdbChainCode().equals(getSecondMolecule().getPdbChainCode()));
	}
	
	/**
	 * Calculates the optimal superposition between the 2 chains of this interface based on 
	 * the (common) CA atoms only. The common CA atoms are takend from matching residue serials, 
	 * thus the 2 chains must correspond to the same sequence.
	 * @return
	 */
	public OptSuperposition getOptimalSuperposition() {
		Tuple3d[][] conformations = getCommonCAConformations();
		Tuple3d[] conformation1 = conformations[0];
		Tuple3d[] conformation2 = conformations[1];

		return GeometryTools.calcOptimalSuperposition(conformation1, conformation2, false);
	}
	
	/**
	 * Returns an array of size 2xn with the 2 conformations of n vectors of CA coordinates
	 * from both chains that are observed in both chains.
	 * The residue serials are used to match the residues from both sides, thus the 2 chains must
	 * correspond to same sequence (even if with different observed residues, i.e. NCS related)
	 * @return
	 */
	private Tuple3d[][] getCommonCAConformations() {
		 
		PdbChain shorterChain = null;
		PdbChain longerChain = null;
		boolean isFirstShorter = false;
		if (firstMolecule.getObsLength() < secondMolecule.getObsLength()) {
			shorterChain = this.firstMolecule;
			longerChain = this.secondMolecule;
			isFirstShorter = true;
		} else { // both equal or second shortest
			shorterChain = this.secondMolecule;
			longerChain = this.firstMolecule;
			isFirstShorter = false;
		}
		
		// we initialise the AL to the shorter chain length, that's the maximum size it can have
		ArrayList<Tuple3d> conf1AL = new ArrayList<Tuple3d>(shorterChain.getObsLength()); 
		ArrayList<Tuple3d> conf2AL = new ArrayList<Tuple3d>(shorterChain.getObsLength());	

		
		for (int resser:shorterChain.getAllResSerials()) {
			Residue shorterChainRes = shorterChain.getResidue(resser);
			if (shorterChainRes.containsAtom("CA")) {
				if (longerChain.containsResidue(resser)) {
					Residue longerChainRes = longerChain.getResidue(resser);
					if (longerChainRes.containsAtom("CA")) {
						Tuple3d shorterChainVector = shorterChainRes.getAtom("CA").getCoords();
						Tuple3d longerChainVector = longerChainRes.getAtom("CA").getCoords();
						if (isFirstShorter) {
							conf1AL.add(shorterChainVector);
							conf2AL.add(longerChainVector);
						} else {
							conf2AL.add(shorterChainVector);
							conf1AL.add(longerChainVector);							
						}
					}
				}
			}
		}
		if (conf1AL.size()!=conf2AL.size()) throw new NullPointerException("Conformations of different size!");
		
		// converting the ArrayLists to arrays		
		Tuple3d[] conformation1 = new Tuple3d[conf1AL.size()]; 
		Tuple3d[] conformation2 = new Tuple3d[conf2AL.size()];
		conf1AL.toArray(conformation1);
		conf2AL.toArray(conformation2);
		Tuple3d[][] conformations = {conformation1, conformation2};
		return conformations;
	}
	
}
