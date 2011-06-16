package owl.core.structure;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;

import javax.vecmath.Matrix4d;

import owl.core.runners.NaccessRunner;
import owl.core.runners.PymolRunner;
import owl.core.structure.graphs.AICGraph;

public class ChainInterface implements Comparable<ChainInterface>, Serializable {
	
	private static final long serialVersionUID = 1L;

	public static final String TYPE_PROTEIN = "Protein";
	
	private static final String TN_STYLE = "cartoon";
	private static final String TN_BG_COLOR = "white";
	private static final int[] TN_HEIGHTS = {75};
	private static final int[] TN_WIDTHS = {75};
	

	private int id;
	private String name;
	private double interfaceArea;
	private AICGraph graph;
	
	private PdbChain firstMolecule;
	private PdbChain secondMolecule;
	
	private InterfaceRimCore firstRimCore;  // cached first molecule's rim and core
	private InterfaceRimCore secondRimCore; // cached second molecule's rim and core
	
	private double bsaToAsaCutoff;
	private boolean zoomingUsed;
	private double bsaToAsaSoftCutoff; // the hard cutoff is stored in the bsaToAsaCutoffs var
	private double bsaToAsaRelaxStep;
	
	private Matrix4d firstTransf; 		// the transformation applied to first molecule expressed in crystal axes coordinates
	private Matrix4d firstTransfOrth;	// the transformation applied to first molecule expressed in orthonormal axes coordinates
	private Matrix4d secondTransf; 		// the transformation applied to second molecule expressed in crystal axes coordinates
	private Matrix4d secondTransfOrth;	// the transformation applied to second molecule expressed in orthonormal axes coordinates
	
	private String firstMolType;	// interface descriptions taken from PISA can have ligands as members, 
									// our own interfaces can only be protein (TYPE_PROTEIN constant)
	private String secondMolType;
	
	private double score; 			// a score value assigned to the interface (if from PISA this is the solvation energy) 
	
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
		this.firstMolType = TYPE_PROTEIN;
		this.secondMolType = TYPE_PROTEIN;
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
	
	public double getScore() {
		return score;
	}
	
	public void setScore(double score) {
		this.score = score;
	}
	
	public String getFirstMolType() {
		return firstMolType;
	}
	
	public void setFirstMolType(String firstMolType) {
		this.firstMolType = firstMolType;
	}
	
	public String getSecondMolType() {
		return secondMolType;
	}
	
	public void setSecondMolType(String secondMolType) {
		this.secondMolType = secondMolType;
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
	 * Calculates the accessible surface area for separate molecules and 2 molecules-complex setting
	 * both ASAs and BSAs of the residues of the two molecules making up this interface.
	 * The total interface area is also calculated from the individual residue values (use {@link #getInterfaceArea()} 
	 * to get it)
	 * @param nSpherePoints
	 * @param nThreads 
	 * @param hetAtoms if true HET residues are considered, if false they aren't, equivalent to 
	 * NACCESS' -h option
	 */
	public void calcSurfAccess(int nSpherePoints, int nThreads, boolean hetAtoms) {
		// NOTE in principle it is more efficient to calculate the ASAs only once per isolated chain
		// BUT! surprisingly the rolling ball algorithm gives slightly different values for same molecule in different 
		// orientations! (can't really understand why!)
		// That's why we calculate ASAs always for 2 separate member of interface and the complex, otherwise 
		// we get (not very big but annoying) discrepancies and also things like negative (small) bsa values

		
		firstMolecule.calcASAs(nSpherePoints, nThreads, hetAtoms);
		secondMolecule.calcASAs(nSpherePoints, nThreads, hetAtoms);
		
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
		return this.firstMolType.equals(TYPE_PROTEIN);
	}
	
	public boolean isSecondProtein() {
		return this.secondMolType.equals(TYPE_PROTEIN);
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
		ps.println("1\t"+firstMolecule.getPdbChainCode()+"\t"+this.getFirstMolType());

		InterfaceRimCore rimCore = getFirstRimCore();
		ps.printf("## %4.2f\n",bsaToAsaCutoff);
		ps.println("## rim : "+rimCore.getRimResString());
		if (rimCore.getCoreResidues().size()>0) {
			ps.println("## core: "+rimCore.getCoreResString());
		}

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
		ps.println("2\t"+secondMolecule.getPdbChainCode()+"\t"+this.getSecondMolType());

		InterfaceRimCore rimCore = getSecondRimCore();
		ps.printf("## %4.2f\n",bsaToAsaCutoff);
		ps.println("## rim : "+rimCore.getRimResString());
		if (rimCore.getCoreResidues().size()>0) {
			ps.println("## core: "+rimCore.getCoreResString());
		}

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
	 * interface storing result in cached arrays. 
	 * Use {@link #getFirstRimCores()} and {@link #getSecondRimCores()} to retrieve them.
	 * The zooming procedure is: the sum of the residues of the 2 cores is required to be at least minNumResidues, 
	 * if the minimum is not reached with the bsaToAsaSoftCutoff, then the cutoff is 
	 * relaxed in relaxationStep steps until reaching the bsaToAsaHardCutoff.
	 * If either of the 2 molecules of this interface is not a protein, its cached rimCore 
	 * object will be null.
	 * @param bsaToAsaSoftCutoff
	 * @param bsaToAsaHardCutoff
	 * @param relaxationStep
	 * @param minNumResidues
	 * @return
	 */
	public void calcRimAndCore(double bsaToAsaSoftCutoff, double bsaToAsaHardCutoff, double relaxationStep, int minNumResidues) {
		zoomingUsed = true;
		
		bsaToAsaCutoff = bsaToAsaHardCutoff;
		this.bsaToAsaSoftCutoff = bsaToAsaSoftCutoff;
		this.bsaToAsaRelaxStep = relaxationStep;
				
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
	 * lists in a cached variable
	 * Use {@link #getFirstRimCore()} and {@link #getSecondRimCore()} to retrieve them.
	 * (see getRimAndCore in {@link PdbChain})
	 * If either of the 2 molecules of this interface is not a protein, the map will be empty for it. 
	 * @param bsaToAsaCutoff
	 * @return
	 */
	public void calcRimAndCore(double bsaToAsaCutoff) {
		zoomingUsed = false;
		
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
	
	public double getBsaToAsaCutoff() {
		return bsaToAsaCutoff;
	}
	
	public double getBsaToAsaSoftCutoff() {
		return bsaToAsaSoftCutoff;
	}
	
	public double getBsaToAsaRelaxStep() {
		return bsaToAsaRelaxStep;
	}
	
	public InterfaceRimCore getFirstRimCore() {
		return firstRimCore;
	}
	
	public InterfaceRimCore getSecondRimCore() {
		return secondRimCore;
	}
	
	/**
	 * Writes this interface to given PDB file with original chain names (PDB chain codes),
	 * unless the 2 chains are the same where the second one is renamed to next letter in 
	 * alphabet.
	 * @param file
	 * @throws FileNotFoundException 
	 */
	public void writeToPdbFile(File file) throws FileNotFoundException {
		String chain2forOutput = secondMolecule.getPdbChainCode();
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
	
	/**
	 * Generate thumbnail files for this interface with PyMol for given pdbFile. 
	 * Output files will be written to same dir as pdbFile using given base name.
	 * @param pymolExe
	 * @param pdbFile
	 * @param base
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws PdbLoadException
	 */
	public void generateThumbnails(File pymolExe, File pdbFile, String base) throws IOException, InterruptedException, PdbLoadException {
		PymolRunner pr = new PymolRunner(pymolExe);
		File[] pngFiles = new File[TN_HEIGHTS.length];
		for (int i=0;i<TN_HEIGHTS.length;i++) {
			pngFiles[i] = new File(pdbFile.getParent(),base+"."+TN_WIDTHS[i]+"x"+TN_HEIGHTS[i]+".png");
		}
		pr.generatePng(pdbFile, pngFiles, TN_STYLE, TN_BG_COLOR, TN_HEIGHTS, TN_WIDTHS);
	}
	
	/**
	 * Tells whether rim/core residues were calculated with zooming (true)
	 * or with fixed cutoff (false).
	 * @return
	 */
	public boolean isRimAndCoreZoomed() {
		return zoomingUsed;
	}
	
}
