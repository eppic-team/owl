package owl.core.structure;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.ArrayList;
import java.util.TreeSet;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;

import owl.core.sequence.Sequence;
import owl.core.sequence.alignment.MultipleSequenceAlignment;
import owl.core.structure.features.CatalSiteSet;
import owl.core.structure.features.EC;
import owl.core.structure.features.Scop;
import owl.core.structure.features.ScopRegion;
import owl.core.structure.features.SecStrucElement;
import owl.core.structure.features.SecondaryStructure;
import owl.core.structure.graphs.AICGEdge;
import owl.core.structure.graphs.AICGraph;
import owl.core.structure.graphs.AIGEdge;
import owl.core.structure.graphs.AIGNode;
import owl.core.structure.graphs.AIGraph;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.BoundingBox;
import owl.core.util.FileFormatException;
import owl.core.util.Grid;
import owl.core.util.Interval;
import owl.core.util.IntervalSet;
import owl.core.util.MySQLConnection;

import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

import Jama.Matrix;
import Jama.SingularValueDecomposition;


import java.sql.SQLException;
import java.sql.Statement;

/**
 * A single chain of a PDB protein structure
 * 
 */
public class PdbChain implements Serializable, Iterable<Residue> {

	private static final long serialVersionUID = 1L;
	
	public static final String DEFAULT_CHAIN   = "A";		// used when constructing models, this will be the chain assigned to them
	public static final String NO_PDB_CHAIN_CODE = "";		// to specify no pdb chain code
	public static final String NO_CHAIN_CODE     = "";		// to specify no internal chain code
	public static final String DEFAULT_CASP_TS_CHAINCODE = " "; // Casp TS format allows only empty chain codes
	
	private static final double MIN_AVRG_NUM_ATOMS_RES = 6.5;	// the cutoff to consider that the average number of atoms per residue 
																// corresponds to that of an all atoms protein. See isAllAtom()
		
	/*-------------------------------------  members ---------------------------------------------*/

	// the asym unit parent
	private PdbAsymUnit parent;
	
	// atom/residue data
	private TreeMap<Integer, Residue> residues;			// residue serials to residue object references (only observed residues)
	private TreeMap<Integer, Atom>    atomser2atom;		// atom serials to Atom object references, we keep this is as a separate map to speed up searches
	private TreeMap<Integer,String> resser2pdbresser;	// internal residue serials to pdb (author) residue serials (can include insertion codes so they are strings)
	private TreeMap<String,Integer> pdbresser2resser; 	// pdb (author) residue serials (can include insertion codes so they are strings) to internal residue serials
	
	private Sequence sequence; 							// full sequence as it appears in SEQRES field

	// identifiers 
	/**
	 * The author assigned chain code (PDB chain code). 
	 * Called auth_asym_id in CIF files.
	 * Not always consistently named, must be one character exactly. 
	 */
	private String pdbChainCode;	
	
	/**
	 * The CIF chain code.
	 * Called asym_id or label_asym_id in CIF files.
	 * More consistently named than PDB chain codes (always consecutive 
	 * A,B,C,...,Z,AA,BA,...), every distinct non-polymeric chain has a different 
	 * one assigned and it can have more than one character.
	 * When we read from PDB files we assign CIF chain codes trying to follow the
	 * assignment criteria of the PDB but still in many cases, especially for 
	 * non-polymer chains, they don't coincide.
	 */
	private String chainCode;			

	private String sid;					// the scop id if PdbChain has been restricted (restrictToScopDomain)	
	
	// sequence features (annotations)
	private transient SecondaryStructure secondaryStructure;	// the secondary structure annotation for this pdb object (should never be null)
	private transient Scop scop;								// the scop annotation for this pdb object
	private transient EC ec;									// the ec annotation for this pdb object
	private transient CatalSiteSet catalSiteSet;				// the catalytic site annotation for this pdb object
	
	
	// optional fields for structures based on casp predictions
	private int targetNum;
	private int caspModelNum;
	private String caspAuthorStr;
	private String caspMethodStr;
	private int groupNum;
	private String[] caspParents;		// optional list of parents used for modelling, may be null

	private boolean hasASA; 			// true if naccess has been run and ASA values assigned
	private boolean hasBfactors; 		// true if atom b-factors have been assigned
	
	private boolean hasAltCodes;		// true if while parsing this chain we found alt codes (we don't store them)
	
	private boolean isNonPolyChain;		// true if chain is a non-polymer chain (purely a HET residues one, ligands and other hets), i.e. not polypeptide or nucleotide chain 

	private BoundingBox bounds; 		// cached bounding box (calculated in getAllAtoms() from old atoms in chain) to speed up getAICGraph()
	
	/*----------------------------------  constructors -----------------------------------------------*/

	/**
	 * Constructs an empty PdbChain with no sequence, no residues and no atoms 
	 */
	public PdbChain() {
		this.chainCode = DEFAULT_CHAIN;
		this.pdbChainCode = DEFAULT_CHAIN;
		
		this.hasASA = false;
		this.hasBfactors = false;
		this.hasAltCodes = false;
				
		this.initialiseResidues();
	}
	
	/**
	 * Constructs a PdbChain for given protein sequence setting coordinates of the given atom type to 
	 * the given coordinates coords
	 * @param sequence
	 * @param coords the array of all coordinates, must be ordered as the sequence and be
	 * of the same size 
	 * @param atom the atom code for which to set the coordinates, e.g. "CA"
	 * @throws IllegalArgumentException if sequence contains invalid characters or if array
	 * of coordinates is of different length as sequence
	 */
	public PdbChain(Sequence sequence, Vector3d[] coords, String atom) {
		this.chainCode = DEFAULT_CHAIN;
		this.pdbChainCode = DEFAULT_CHAIN;
		
		this.hasASA = false;
		this.hasBfactors = false;
		this.hasAltCodes = false;
		
		if (coords.length!=sequence.getLength()) {
			throw new IllegalArgumentException("Array of coordinates is not of same length as given sequence");
		}
		
		this.sequence = sequence;
		this.initialiseResidues();
		this.resser2pdbresser = new TreeMap<Integer, String>();
		this.pdbresser2resser = new TreeMap<String, Integer>();

		for (int i=0;i<sequence.getLength();i++) {
			int resser = i+1;
			char one = sequence.getSeq().charAt(i);
			if (!AminoAcid.isStandardAA(one)) {
				throw new IllegalArgumentException("Given input sequence to construct a pdb model contains an invalid aminoacid "+one);
			}
			this.addResidue(new AaResidue(AminoAcid.getByOneLetterCode(one),resser,this));
			resser2pdbresser.put(resser, String.valueOf(resser));
			pdbresser2resser.put(String.valueOf(resser), resser);
		}
		
		int i = 0;
		for (Residue residue:this.residues.values()) {
			// the element is set to first letter of given atom code, that will work only for atoms in standard aas (C,H,S,N,O)
			residue.addAtom(new Atom(i+1, atom, String.valueOf(atom.charAt(0)), new Point3d(coords[i]), residue, Atom.DEFAULT_OCCUPANCY, Atom.DEFAULT_B_FACTOR));
			i++;
		}
		this.initialiseMaps();
	}
	
	/*---------------------------------  public methods ----------------------------------------------*/

	@Override
	public Iterator<Residue> iterator() {
		return residues.values().iterator();
	}
	
	/**
	 * Initialises the residues to an empty map. 
	 */
	protected void initialiseResidues() {
		residues = new TreeMap<Integer, Residue>();
	}
	
	/**
	 * Adds a residue to this PdbChain
	 * @param residue
	 */
	public void addResidue(Residue residue) {
		Residue returnRes = residues.put(residue.getSerial(), residue);
		if (returnRes!=null) {
			System.err.println("Warning: a residue with the same residue serial already exists in this chain. New residue: "+residue+", old residue: "+returnRes);
		}
	}
	
	/**
	 * Gets a Residue given its serial
	 * @param resSerial
	 * @return the reference to the Residue object or null if residue serial doesn't exist
	 * for instance because it's not observed
	 */
	public Residue getResidue(int resSerial) {
		
		return residues.get(resSerial);

	}
	
	/**
	 * Tells whether residue of given residue number is a (observed) residue (standard 
	 * amino acid, nucleotide or het residue) in this PdbChain instance. 
	 * See also {@link #hasCoordinates(int)}
	 * @param resSerial
	 * @return
	 */
	public boolean containsResidue(int resSerial) {
		return residues.containsKey(resSerial);
	}
	
	/**
	 * Tells whether residue of given residue number is a (observed) standard amino 
	 * acid residue in this PdbChain instance. 
	 * See also {@link #hasCoordinates(int)} 
	 * @param resSerial
	 * @return
	 */
	public boolean containsStdAaResidue(int resSerial) {
		if (containsResidue(resSerial)) {
			if (getResidue(resSerial) instanceof AaResidue) {
				return true;
			}
		}
		return false;
	}	
	
	/**
	 * Gets all residues of specified amino acid type in a List
	 * @param aa the amino acid type (see AminoAcid class)
	 * @return
	 */
	public ArrayList<AaResidue> getResiduesOfType(AminoAcid aa) {
		ArrayList<AaResidue> list = new ArrayList<AaResidue>();
		for (Residue residue:residues.values()) {
			if (residue instanceof AaResidue) {
				AaResidue aares = (AaResidue) residue;
				if (aares.getAaType().equals(aa)) {
					list.add(aares);
				}
			}
		}
		return list;
	}
	
	/**
	 * Gets a new Map with residue serials to Residues that contain only the atoms for the 
	 * given contact type and given interval set (standard aminoacids only). The Residue objects 
	 * are new, but the Atom objects to which they point to are the same old references.
	 * @param ct the contact type
	 * @param intervSet only residues of this intervals will be considered, if null then
	 * all residues taken
	 * @return
	 */
	private TreeMap<Integer, AaResidue> getReducedResidues(String ct, IntervalSet intervSet) {
		TreeMap<Integer,AaResidue> reducedResidues = new TreeMap<Integer, AaResidue>();
		
		if (intervSet!=null) {
			for (Interval interv:intervSet) {
				for (int resser=interv.beg;resser<=interv.end;resser++) {
					Residue residue = getResidue(resser);
					if (residue==null) 
						throw new IllegalArgumentException("Invalid interval specified, residue "+resser+" is not part of this PdbChain");
					if (!(residue instanceof AaResidue)) continue;
					reducedResidues.put(resser,((AaResidue)residue).getReducedResidue(ct));
				}
			}
		} else { // we take all observed residues
			for (Residue residue:residues.values()) {
				if (!(residue instanceof AaResidue)) continue;
				reducedResidues.put(residue.getSerial(),((AaResidue)residue).getReducedResidue(ct));
			}
		}
		return reducedResidues;
	}
	
	/**
	 * Gets an array of Atoms that contain only the atoms for the given 
	 * contact type and given interval set for standard aminoacids only.
	 *
	 * @param ct the contact type
	 * @param intervSet only residues of this intervals will be considered, if null then 
	 * all residues taken
	 * @return
	 */
	private Atom[] getAtomsForCt(String ct, IntervalSet intervSet) {
		TreeMap<Integer, AaResidue> reducedResidues = getReducedResidues(ct, intervSet);
		int totalAtoms = 0;
		for (AaResidue residue:reducedResidues.values()) {
			totalAtoms+=residue.getNumAtoms(); 
		}
		Atom[] atoms = new Atom[totalAtoms];
		int i = 0;
		for (AaResidue residue:reducedResidues.values()) {
			for (Atom atom:residue.getAtoms()) {
				atoms[i]=atom;
				i++;
			}
		}
		return atoms;
	}
	
	/**
	 * Gets an array with all atoms present in this chain. Atoms of standard aminoacids,
	 * hetatoms and nucleotides will be taken.
	 * Calculates also the bounds array so that it gets cached and can be used in {@link #getAICGraph(PdbChain, double)}
	 * @return
	 */
	protected Atom[] getAllAtoms() {
 		Atom[] atoms = new Atom[this.getNumAtoms()];
		int i = 0;
		for (Residue residue:this) {
			for (Atom atom:residue) {
				atoms[i]=atom;
				i++;
			}
		}
		if (this.bounds==null) {
			this.bounds = new BoundingBox(atoms);
		}
		return atoms;		
	}
	
	protected BoundingBox getBoundingBox() {
		if (bounds!=null) {
			return bounds;
		} else {
			getAllAtoms();
			return bounds;
		}
	}
	
	protected boolean isNotOverlapping(PdbChain other, double cutoff) {
		BoundingBox thisbb = this.getBoundingBox();
		BoundingBox otherbb = other.getBoundingBox();
		return !thisbb.overlaps(otherbb, cutoff);
	}
	
	/**
	 * Gets an Atom given its serial
	 * @param atomSerial
	 * @return
	 */
	public Atom getAtom(int atomSerial) {
		return atomser2atom.get(atomSerial);
	}
	
	/**
	 * Populates the atomser2atom map (from the Residues and Atoms objects)
	 */
	protected void initialiseMaps() {
		atomser2atom = new TreeMap<Integer, Atom>();
		for (Residue residue:residues.values()) {
			for (Atom atom:residue.getAtoms()) {
				atomser2atom.put(atom.getSerial(), atom);
			}
		}
		
	}
	
	/**
	 * Sets the secondary structure members of the Residue objects from the SecondaryStructure object
	 */
	protected void initialiseResiduesSecStruct() {
		Iterator<SecStrucElement> it = this.getSecondaryStructure().iterator();
		while (it.hasNext()) {
			SecStrucElement ssElem = it.next();
			for (int resser=ssElem.getInterval().beg;resser<=ssElem.getInterval().end;resser++) {
				if (this.containsResidue(resser)) {
					Residue residue = this.getResidue(resser);
					if (residue instanceof AaResidue) {
						((AaResidue)residue).setSsElem(ssElem);
					} else if (residue instanceof HetResidue) {
						((HetResidue)residue).setSsElem(ssElem);
					}
				} else {
					// we don't warn because this does happen really often!
					//System.err.println("Warning: the residue serial "+resser+" can't be assigned with the secondary structure element "+
					//		ssElem.getId()+" with interval "+ssElem.getInterval()+" because it's not present in this PdbChain (e.g. it's not observed)");
				}
			}			
		}
	}
	
	/**
	 * Calculates compactness coefficient as defined by:
	 * Micheal H. Zehfus and George D. Rose (1986) Compact Units in Proteins. Biochemistry 25: 5759-5765.
	 * DOI: 10.1021/bi00367a062
	 * @param surface
	 * @param volume
	 * @return
	 */
	public double calcCompactnessCoefficient(double surface, double volume) {
		return (surface/Math.pow(36*Math.PI*Math.pow(volume,2), 1.0/3.0));
	}
	
	/**
	 * Returns the per-residue relative solvent accessible surface areas (SASA) as 
	 * calculated by NACCESS. Only considers standard amino acids.
	 * Returns null if SASA has not previously been calculated with {@link runner.NaccessRunner}.
	 * @return a map from residue serial to SASA value (null if not calculated yet)
	 */
	public HashMap<Integer, Double> getSurfaceAccessibilities() {
		if (hasASA()) {
			HashMap<Integer, Double> resser2allrsa = new HashMap<Integer, Double>();
			for (Residue residue:residues.values()) {
				if (residue instanceof AaResidue) { 
					resser2allrsa.put(residue.getSerial(), residue.getRsa());
				}
			}
			return resser2allrsa;
		} else {
			return null;
		}
	}

	/**
	 * Returns the per-residue side-chain relative solvent accessible surface areas (SASA)
	 * as calculated by NACCESS. Only considers standard amino acids.
	 * Returns null if SASA has not previously been calculated with {@link runner.NaccessRunner}.
	 * @return a map from residue serial to SASA value (null if not calculated yet)
	 */
	public HashMap<Integer, Double> getSideChainSurfaceAccessibilities() {
		if (hasASA()) {
			HashMap<Integer, Double> resser2scrsa = new HashMap<Integer, Double>();
			for (Residue residue:residues.values()) {
				if (residue instanceof AaResidue) {
					resser2scrsa.put(residue.getSerial(), ((AaResidue)residue).getScRsa());
				}
			}
			return resser2scrsa;
		} else {
			return null;
		}
	}
	
	/**
	 * Returns the per-residue all-atoms absolute solvent accessible areas as calculated by NACCESS.
	 * Only considers standard amino acids.
	 * @return a map from residue serial to absolute ASA value or null if ASA has not previously been 
	 * calculated with {@link runner.NaccessRunner}.
	 */
	public HashMap<Integer, Double> getAbsSurfaceAccessibilities() {
		if (hasASA()) {
			HashMap<Integer, Double> resser2allrsa = new HashMap<Integer, Double>();
			for (Residue residue:residues.values()) {
				if (residue instanceof AaResidue) {
					resser2allrsa.put(residue.getSerial(), residue.getAsa());
				}
			}
			return resser2allrsa;
		} else {
			return null;
		}
	}

	/**
	 * Assigns b-factor values to the atoms of this structure. If structure is written to pdb file,
	 * these values will appear in the b-factor column. Currently, if no b-factors are assigned, the
	 * default value 0 will be written.
	 * @param bfactorsPerAtom a map of atom serials to b-factors
	 */
	public void setBFactorsPerAtom(HashMap<Integer, Double> bfactorsPerAtom) {
		for(int atomser:bfactorsPerAtom.keySet()) {
			getAtom(atomser).setBfactor(bfactorsPerAtom.get(atomser));
		}
		hasBfactors = true;
	}
	
	/**
	 * Gets the total accessible surface area for this chain by summing up the individual
	 * residue ASAs. Use {@link #calcASAs(int, int, boolean)} before.
	 * @return
	 */
	public double getASA() {
		double asa = 0;
		for (Residue res:this) {
			asa+=res.getAsa();
		}
		return asa;
	}
	
	/**
	 * Gets the total buried surface area for this chain by summing up the individual
	 * residue BSAs. 
	 * The values are filled when ASAs are calculated for this chain and this chain in 
	 * a complex with another one. See PdbAsymUnit.calcBSAs()
	 * @return
	 */
	public double getBSA() {
		double bsa = 0;
		for (Residue res:this) {
			bsa+=res.getBsa();
		}
		return bsa;
	}
	
	/**
	 * Sets the VdW radius values of all atoms of this PdbChain (from standard amino acids,
	 * het residues and nucleotides) 
	 * Use subsequently Atom.getRadius() to get the value.
	 * This uses the AtomRadii parser of the vdw.radii resource file.
	 */
	protected void setAtomRadii() {
		for (Residue residue:residues.values()) {
			residue.setAtomRadii();
		}
	}
	
	/**
	 * Calculate the Accessible Surface Area using our implementation of the 
	 * rolling ball algorithm. Sets both the Atoms' and Residues' asa members 
	 * See Shrake, A., and J. A. Rupley. "Environment and Exposure to Solvent of Protein Atoms. 
	 * Lysozyme and Insulin." JMB (1973) 79:351-371.
	 * @param nSpherePoints the number of points to be used in generating the spherical 
	 * dot-density, the more points the more accurate (and slower) calculation
	 * @param nThreads number of threads to be used for ASA calculation
	 * @param hetAtoms if true HET residues are considered, if false they aren't, equivalent to 
	 * NACCESS' -h option
	 */
	public void calcASAs(int nSpherePoints, int nThreads, boolean hetAtoms) {
		this.setAtomRadii();
		int numAtoms = getNumNonHetAtoms();
		if (hetAtoms) {
			numAtoms = getNumAtoms();
		}
		Atom[] atoms = new Atom[numAtoms];
		int i = 0;
		for (Residue residue: this) {
			if (!hetAtoms && (residue instanceof HetResidue)) continue;
			for (Atom atom: residue) {
				atoms[i] = atom;
				i++;
			}
		}
		double[] asas = Asa.calculateAsa(atoms,Asa.DEFAULT_PROBE_SIZE,nSpherePoints,nThreads);
		
		for (i=0;i<atoms.length;i++) {
			atoms[i].setAsa(asas[i]);
		}

		// and finally sums per residue
		for (Residue residue: this) {
			if (!hetAtoms && (residue instanceof HetResidue)) continue;
			double tot = 0;
			for (Atom atom:residue.getAtoms()) {
				tot+=atom.getAsa();
			}
			residue.setAsa(tot);
			if (residue instanceof AaResidue) {
				AaResidue aares = (AaResidue) residue;
				residue.setRsa(tot/aares.getAaType().getAsaInExtTripept());
			}
		}
		this.hasASA = true;
	}
	
	protected void setASAs(PdbChain other) {
		if (!other.hasASA()) {
			System.err.println("Error! given chain has no ASAs from which to set the ASAs of this chain");
			return;
		}
		
		for (Residue residue: this) {
			
			for (Atom atom: residue) {
				if (other.getAtom(atom.getSerial())==null) {
					System.err.println("Error! atom serial "+atom.getSerial()+" does not exist in given chain. Won't set ASAs of this chain");
					return;
				}
				atom.setAsa(other.getAtom(atom.getSerial()).getAsa()); 
			}
		}
		// and finally sums per residue
		for (Residue residue: this) {
			double tot = 0;
			for (Atom atom:residue.getAtoms()) {
				tot+=atom.getAsa();
			}
			residue.setAsa(tot);
			if (residue instanceof AaResidue) {
				AaResidue aares = (AaResidue) residue;
				residue.setRsa(tot/aares.getAaType().getAsaInExtTripept());
			}
		}
		this.hasASA = true;
	}
	
	/**
	 * Calculates for each atom (of standard amino acids) in this structure the deviation to the corresponding atom
	 * in the reference structure and returns a map from atom serials to differences in Angstrom.
	 * The corresponding atom is the atom with the same residue serial and the same PDB residue code
	 * (e.g. CA, CB, N), ignoring residues which are mutated between this and the reference structure.
	 * If no corresponding atom is found, the value Atom.DEFAULT_B_FACTOR will be used. To give useful
	 * results the two structures need to have the same residue numbering (save mutations) and have to
	 * be properly superimposed. This holds for example for two consecutive states in an MD trajectory
	 * or for different FoldX mutant models based in the same original structure.
	 * The returned values can be assigned to the b-factor column using this.setBFactorsPerAtom()
	 * and can then be visualized in PyMol with the command 'spectrum b'. 
	 * @param referencePdb the reference structures to compare to
	 * @return the per atom distances between this and the reference structure
	 */
	public HashMap<Integer, Double> getPerAtomDistances(PdbChain referencePdb) {
		HashMap<Integer, Double> distances = new HashMap<Integer, Double>();
		
		// first set all distances to the default value, so that the map contains the full set of atoms
		for (Residue res:this) {
			if (!(res instanceof AaResidue)) continue;
			for (Atom atom:res) {
				distances.put(atom.getSerial(),Atom.DEFAULT_B_FACTOR);
			}
		}
		
		// then set the real distance whereever possible
		for(int resSer: this.getAllResSerials()) {
			Residue r = this.getResidue(resSer);
			if (!(r instanceof AaResidue)) continue;
			AminoAcid thisType = ((AaResidue)r).getAaType();
			if(referencePdb.hasCoordinates(resSer) && (referencePdb.getResidue(resSer) instanceof AaResidue)) {
				AminoAcid refType = ((AaResidue)referencePdb.getResidue(resSer)).getAaType();
				if(thisType == refType) {
					for(Atom a1:r.getAtoms()) {
						int atomSer = a1.getSerial();
						if(referencePdb.hasCoordinates(resSer, a1.getCode())) {
							Atom a2 = referencePdb.getAtom(referencePdb.getAtomSerFromResSerAndAtom(resSer, a1.getCode()));
							// a2 should never be null
							double dist = Math.abs(a1.getCoords().distance(a2.getCoords()));
							distances.put(atomSer,dist);
						}
					}
				} else {
					System.err.println("Skipping mismatching residue at position " + resSer + " ("+ thisType + "!=" + refType + ")");
				}
			} else {
				System.err.println("Skipping " + thisType.getThreeLetterCode() + " at position " + resSer + ". Not found in reference structure.");
			}
		}
		return distances;
	}
	
	/**
	 * Assigns b-factor values to all atoms for the given residues.
	 * @param bfactorsPerResidue a map of residue serials to b-factors
	 * @throws IllegalArgumentException when residue serials given in input are not
	 * present in this PdbChain instance
	 */
	public void setBFactorsPerResidue(HashMap<Integer, Double> bfactorsPerResidue) {
		for(int resser:bfactorsPerResidue.keySet()) {
			if (this.containsResidue(resser)) {
				Residue residue = this.getResidue(resser);
				double bfactor = bfactorsPerResidue.get(resser);
				for(Atom atom:residue) {
					atom.setBfactor(bfactor);
				}
			} else {
				System.err.println("Warning! Can't assign bfactor for residue serial "+resser+", it is not present in this PdbChain instance, pdbCode "+getPdbCode()+", pdbChainCode "+pdbChainCode);
			}
		}
		hasBfactors = true;
	}
	
	/**
	 * Write CASP TS file headers complying with the CASP 8 specification.
	 * Note that the CASP target number, CASP model number, CASP author
	 * and CASP method will be written from the internally 
	 * set values (targetNum, caspModelNum, caspAuthorStr, caspMethodStr) 
	 * so they must be set before trying to write them out.
	 * If the method string has new line characters it will be splitted into
	 * several METHOD lines
	 * @param Out
	 * @param parents PDB entries in which this homology prediction is based on or
	 * null if not applicable i.e. if this is an ab-initio prediction
	 */
	private void writeCaspTSHeader(PrintStream Out, String[] parents) {
		Out.println("PFRMAT TS");
		Out.printf("TARGET T%04d\n",targetNum);
		if(caspAuthorStr != null) Out.printf("AUTHOR %s\n", caspAuthorStr);
		if(caspMethodStr != null) {
			String[] caspMethodLines = caspMethodStr.split("\n");
			for (String caspMethodLine: caspMethodLines) {
				Out.printf("METHOD %s\n", caspMethodLine);
			}
		}
		Out.println("MODEL "+caspModelNum);
		String parentStr = "";
		if (parents==null) {
			parentStr = "N/A";
		} else {
			for (int i=0;i<parents.length;i++) {
				Pattern p = Pattern.compile("(\\d\\w\\w\\w)(\\w)");
				Matcher m = p.matcher(parents[i]);
				if (m.matches()) {
					parentStr += m.group(1)+"_"+m.group(2) + " ";
				} else {
					System.err.println("Parent string "+parents[i]+" to be written to file CASP TS file not in the expected format 1abcA ");
					parentStr += parents[i] + " ";
				}
				
			}
		}
		Out.println("PARENT "+parentStr);
	}
	
	/** 
	 * Writes atom lines for this structure to the given output stream.
	 * The chain code written is the CIF chain code. If the code is longer than 1 character 
	 * only first character is used and a warning issued.
	 * @param out
	 */
	public void writeAtomLines(PrintStream out) {
		writeAtomLines(out,chainCode);
	}
	
	/**
	 * Writes atom lines for this structure to the given output stream using the given
	 * chain code instead of the internal CIF chain code.
	 * @param out
	 * @param chainCode
	 */
	public void writeAtomLines(PrintStream out, String chainCode) {
		String chainCodeStr = chainCode;
		if (chainCode.length()>1) {
			System.err.println("Warning! Chain code with more than 1 character ("+chainCode+"), only first character will be written to ATOM lines");
			chainCodeStr = chainCode.substring(0,1);
		}
		// we write atoms sorted by atom serial, 
		// if we were iterating over residues and then atoms we'd get them in order of residue and atom code alphabetical order
		for (int atomser:this.getAllAtomSerials()) {
			Atom atom = getAtom(atomser);
			int resser = atom.getParentResSerial();
			String atomCode = atom.getCode();
			String atomType = atom.getType().getSymbol();
			String res = atom.getParentResidue().getLongCode();
			Point3d coords = atom.getCoords();
			double occupancy = atom.getOccupancy();
			double bFactor = atom.getBfactor();
			String lineType = "ATOM  ";
			if (atom.getParentResidue() instanceof HetResidue) {
				lineType = "HETATM";
			}
			String printfStr = lineType+"%5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n";
			if (atomCode.length()==4) { // some hydrogens have a 4 letter code and it is not aligned to column 14 but to 13 instead 
				printfStr =    lineType+"%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n";
			}
			if (atomType.length()==2) { // for atoms with 2 letter codes (CL, NA, ....) the alignment is to the left again
				printfStr =    lineType+"%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n";
			}
			// Local.US is necessary, otherwise java prints the doubles locale-dependant (i.e. with ',' for some locales)
			out.printf(Locale.US, printfStr,
					atomser, atomCode, res, chainCodeStr, resser, coords.x, coords.y, coords.z, occupancy, bFactor, atomType);
		}
		if (!isNonPolyChain()) out.println("TER");
	}

	/**
	 * Writes to given PrintStream the SEQRES record for this chain using the given 
	 * chain code instead of the internal CIF chain code.
	 * @param out
	 * @param chainCode
	 */
	protected void writeSeqresRecord(PrintStream out, String chainCode) {
		String chainCodeStr = chainCode;
		if (chainCode.length()>1) {
			System.err.println("Warning! Chain code with more than 1 character ("+chainCode+"), only first character will be written to ATOM lines");
			chainCodeStr = chainCode.substring(0,1);
		}
		ArrayList<String> seqLines = new ArrayList<String>();
		String seqLine = null;
		for (int i=0;i<this.sequence.getSeq().length();i++){
			String longCode = AminoAcid.XXX.getThreeLetterCode();
			if (this.sequence.isProtein()) {
				AminoAcid aa = AminoAcid.getByOneLetterCode(this.sequence.getSeq().charAt(i));
				longCode = aa.getThreeLetterCode();
			} else {
				Nucleotide nuc = Nucleotide.getByOneLetterCode(this.sequence.getSeq().charAt(i));
				longCode = nuc.getTwoLetterCode();				
			}
			if (longCode.equals(AminoAcid.XXX.getThreeLetterCode())) {
				if (containsResidue(i+1)) {
					longCode = getResidue(i+1).getLongCode();				
				}
			}
			if (i%13==0) {
				if (seqLine!=null) seqLines.add(seqLine);
				seqLine = "";
			}
			seqLine += String.format("%3s ",longCode);
		}
		seqLines.add(seqLine); // the last line wasn't added by the loop
		for (int i=0;i<seqLines.size();i++){
			out.printf("SEQRES %3d %s %4d  %s\n",i+1,chainCodeStr,this.getFullLength(),seqLines.get(i));
		}
	}
	
	/**
	 * Writes coordinate data into a file in PDB format (ATOM lines only)
	 * The residue serials written correspond to the SEQRES sequence (as in CIF files).
	 * The chain code written is the CIF chain code, if it is longer than one character
	 * then only first character is written and a warning issued.
	 * @param outFile
	 * @throws FileNotFoundException
	 */
	public void writeToPDBFile(File outFile) throws FileNotFoundException {
		PrintStream out = new PrintStream(new FileOutputStream(outFile));
		if (!isNonPolyChain()) {
			writeSeqresRecord(out, chainCode);
		}
		writeAtomLines(out);
		if (!isNonPolyChain()) {
			out.println("TER");
		}
		out.println("END");
		out.close();
	}

	/**
	 * Writes coordinates to given File in CASP TS format.
     * Note that the CASP target number, CASP model number, CASP author, 
     * CASP method and CASP parents will be written from the internally set values
     * (targetNum, caspModelNum, caspAuthorStr, caspMethodStr, caspParents) so they
     * must be set before trying to write them out.
	 * @param outFile
	 * @throws FileNotFoundException
	 */
	public void writeToCaspTSFile(File outFile) throws FileNotFoundException {
		PrintStream out = new PrintStream(new FileOutputStream(outFile));
		writeCaspTSHeader(out, this.caspParents);
		writeAtomLines(out,DEFAULT_CASP_TS_CHAINCODE);
		out.println("END");
		out.close();
	}
	
	/**
	 * Dump the full sequence of this PdbChain object in FASTA file format
	 * The FASTA tag is written as the concatenation of pdbCode and pdbChainCode
	 * @param seqfile
	 * @throws IOException if file can't be written
	 */
	public void writeSeqToFasta(String seqfile) throws IOException {
		PrintStream Out = new PrintStream(new FileOutputStream(seqfile));
		int len = 80;
		Out.println(">"+getPdbCode()+pdbChainCode);
		for(int i=0; i<sequence.getLength(); i+=len) {
			Out.println(sequence.getSeq().substring(i, Math.min(i+len,sequence.getLength())));
		}		
		Out.close();
	}

	/** 
	 * Returns the number of observed residues (standard amino acids, nucleotides and hets)
	 * @return number of observed residues
	 */
	public int getObsLength(){
		return residues.size();
	}
	
	/**
	 * Returns the number of observed standard amino acid residues.
	 * @return
	 */
	public int getStdAaObsLength() {
		int count=0;
		for (Residue residue:this) {
			if (residue instanceof AaResidue) count++;
		}
		return count;
	}

	/** 
	 * Returns the number of residues in the SEQRES sequence of this protein.
	 * @return number of residues in the full sequence
	 */
	public int getFullLength() {
		return sequence.getLength();
	}

	/**
	 * Returns number of atoms in the protein, including Hydrogens if they are present
	 * Includes all residues: standard amino acids, peptide linked het residues and ligands 
	 * @return number of atoms
	 */
	public int getNumAtoms() {
		int numAtoms = 0;
		for (Residue residue:this) {
			numAtoms+=residue.getNumAtoms();
		}
		return numAtoms;
	}
	
	/**
	 * Returns number of atoms in the protein, including Hydrogens if they are present 
	 * only considering standard amino acid residues or nucleotides if the chains is a nucleic acid.
	 * @return number of atoms
	 */	
	public int getNumNonHetAtoms() {
		int numAtoms = 0;
		for (Residue res:this) {
			if ((res instanceof HetResidue)) continue;
			numAtoms+=res.getNumAtoms();
		}
		return numAtoms;
	}
	
	/**
	 * Returns number of heavy (non-Hydrogen) atoms in the protein
	 * Includes all residues: standard amino acids, peptide linked het residues and ligands  
	 * @return number of (non-Hydrogen) atoms
	 */
	public int getNumHeavyAtoms() {
		int numAtoms = 0;
		for (Residue residue: this) {
			numAtoms+=residue.getNumHeavyAtoms();
		}
		return numAtoms;
	}
	
	/**
	 * Returs the number of heavy (non-Hydrogen) atoms of standard amino acids in this protein chain
	 * @return
	 */
	public int getNumStdAaHeavyAtoms() {
		int numAtoms = 0;
		for (Residue res:this) {
			if (!(res instanceof AaResidue)) continue;
			numAtoms+=res.getNumHeavyAtoms();
		}
		return numAtoms;
	}

	/**
	 * Removes all Hydrogen atoms from this chain
	 */
	public void removeHatoms() {
		for (Residue residue:this) {
			residue.removeHatoms();
		}
		initialiseMaps();
	}
	
	/**
	 * Returns the number of chain breaks in this structure, i.e. how many unobserved 
	 * gaps are there in the chain (excluding unobserved regions in the terminals)
	 * @return
	 */
	public int getNumChainBreaks() {
		int numChainBreaks = 0;
		int lastResser = 0;
		for (int resser:getAllResSerials()) {
			if (lastResser!=0 && resser-lastResser>1) 
				numChainBreaks++;
			lastResser = resser;
		}
		return numChainBreaks;
	}
	
	/**
	 * Returns the distance matrix of standard amino acid residues as a Jama Matrix object. 
	 * For multi atom contact types the distance matrix has the minimum distance for each pair of
	 * residues. 
	 * The indices of the matrix will be 0 to get_length()-1 and will be in the same order 
	 * as the residue serials
	 * @param ct contact type for which distances are being calculated
	 * @return
	 */
	public Matrix calcDistMatrixJamaFormat(String ct) {
		HashMap<Integer,Integer> resser2ind = new HashMap<Integer, Integer>();
		int i = 0;
		for (Residue residue:this) {
			if (residue instanceof AaResidue) {
				resser2ind.put(residue.getSerial(),i);
				i++;
			}
		}
		HashMap<Pair<Integer>, Double> distHM = this.calcDistMatrix(ct);
		Matrix matrix = new Matrix(this.getStdAaObsLength(), this.getStdAaObsLength());
		for (Pair<Integer> pair: distHM.keySet()) {
			double currentElem = distHM.get(pair);
			matrix.set(resser2ind.get(pair.getFirst()), resser2ind.get(pair.getSecond()), currentElem);
			matrix.set(resser2ind.get(pair.getSecond()), resser2ind.get(pair.getFirst()), currentElem);
		}
		return matrix;
	}
	
	/**
	 * Returns the distance matrix of standard amino acid residues as a HashMap with residue serial pairs as keys
	 * For multi atom contact types the distance matrix has the minimum distance for each pair of
	 * residues 
	 * AAinfo.isValidSingleAtomCT(ct) can be used to check before calling.
	 * @param ct contact type for which distances are being calculated
	 * @return a map which assigns to each edge the corresponding distance 
	 */
	public HashMap<Pair<Integer>, Double> calcDistMatrix(String ct){
		HashMap<Pair<Integer>,Double> distMatrixAtoms = calcAtomDistMatrix(ct);

		 // mapping atom serials to residue serials
		 // TODO: we could integrate this with the code in calcAtomDistMatrix to avoid storing two distance maps in memory
		HashMap<Pair<Integer>,Double> distMatrixRes = new HashMap<Pair<Integer>, Double>();
		for (Pair<Integer> cont: distMatrixAtoms.keySet()){
			int i_resser = getResSerFromAtomSer(cont.getFirst());
			int j_resser = getResSerFromAtomSer(cont.getSecond());
			Pair<Integer> edge = new Pair<Integer>(i_resser,j_resser);
			if (distMatrixRes.containsKey(edge)) {
				distMatrixRes.put(edge, Math.min(distMatrixRes.get(edge), distMatrixAtoms.get(cont)));
			} else {
				distMatrixRes.put(edge, distMatrixAtoms.get(cont));
			}
		}

		return distMatrixRes;
	}
	
	/**
	 * Returns the distance matrix as a HashMap with atom serial pairs as keys
	 * Only standard aminoacids are considered.
	 * This method can be used for any contact type 
	 * AAinfo.isValidSingleAtomCT(ct) can be used to check before calling.
	 * @param ct contact type for which distances are being calculated
	 * @return a map which assings to each atom pair the corresponding distance
	 */
	public HashMap<Pair<Integer>, Double> calcAtomDistMatrix(String ct){
		HashMap<Pair<Integer>,Double> distMatrixAtoms = new HashMap<Pair<Integer>,Double>();
		if (!ct.contains("/")){
			Atom[] atoms = getAtomsForCt(ct, null);
			for (Atom iAtom:atoms){
				for (Atom jAtom:atoms){
					int iAtomSer = iAtom.getSerial();
					int jAtomSer = jAtom.getSerial();
					if (jAtomSer>iAtomSer) {
						Pair<Integer> pair = new Pair<Integer>(iAtomSer,jAtomSer);
						distMatrixAtoms.put(pair, iAtom.getCoords().distance(jAtom.getCoords()));
					}
				}
			}
		} else {
			String i_ct = ct.split("/")[0];
			String j_ct = ct.split("/")[1];
			Atom[] iAtoms = getAtomsForCt(i_ct,null);
			Atom[] jAtoms = getAtomsForCt(j_ct,null);
			for (Atom iAtom:iAtoms){
				for (Atom jAtom:jAtoms){
					int iAtomSer = iAtom.getSerial();
					int jAtomSer = jAtom.getSerial();
					if (jAtomSer!=iAtomSer){
						Pair<Integer> pair = new Pair<Integer>(iAtomSer,jAtomSer);
						distMatrixAtoms.put(pair, iAtom.getCoords().distance(jAtom.getCoords()));
					}
				}
			}
		}

		return distMatrixAtoms;
	}
	
	/**
	 * Calculates the radius of gyration of this PdbChain 
	 * (defined as half of the maximum distance between any 2 CA atoms)
	 * @return
	 */
	public double calcRadGyration() {
		//TODO this is a very raw implementation o(n^2): should optimise it if that's really needed
		return Collections.max(this.calcAtomDistMatrix("Ca").values())/2;
	}
	
	/**
	 * Get the atom interaction graph for given contact type and cutoff for this PdbChain object.
	 * Only standard aminoacids are considered.
	 * Returns a Graph object with the contacts
	 * A geometric hashing algorithm is used for fast contact computation (without needing 
	 * to calculate full distance matrix)
	 * @param ct
	 * @param cutoff
	 * @return
	 */
	private AIGraph getAIGraph(String ct, double cutoff){ 
		Atom[] iAtoms = null;
		Atom[] jAtoms = null;		// only relevant for asymetric edge types
		boolean crossed = false;
		if (!ct.contains("/")){
			iAtoms = getAtomsForCt(ct, null);
			crossed = false;
		} else {
			String i_ct = ct.split("/")[0];
			String j_ct = ct.split("/")[1];
			iAtoms = getAtomsForCt(i_ct,null);
			jAtoms = getAtomsForCt(j_ct,null);
			crossed = true;
		}

		Grid grid = new Grid(cutoff);
		
		if (crossed) {
			grid.addAtoms(iAtoms,jAtoms);
		} else {
			grid.addAtoms(iAtoms,iAtoms);
		}
		
		float[][] distMatrix = grid.getDistMatrix(crossed);

		// creating the AIGraph
		AIGraph graph = new AIGraph();
		TreeMap<Integer,RIGNode> rignodemap = new TreeMap<Integer,RIGNode>();
		TreeSet<Integer> atomSerials = new TreeSet<Integer>();
		for (Atom atom:iAtoms){
			atomSerials.add(atom.getSerial());
		}
		if (jAtoms!=null){
			for (Atom atom:jAtoms) {
				atomSerials.add(atom.getSerial());
			}
		}
		// adding the AIGNodes (including parent RIGNode references)
		SecondaryStructure secondaryStructureCopy = secondaryStructure.copy();
		for (int atomSer:atomSerials) {
			int resser = getResSerFromAtomSer(atomSer);
			SecStrucElement sselem = secondaryStructureCopy.getSecStrucElement(resser);
			if (!rignodemap.containsKey(resser)) {
				// NOTE!: we are passing references to the SecStrucElement objects! they point to the same objects as secondaryStructureCopy 
				RIGNode resNode = new RIGNode(resser, ((AaResidue)getResidue(resser)).getAaType().getThreeLetterCode(), sselem);
				rignodemap.put(resser,resNode);
			}
			AIGNode atomNode = new AIGNode(atomSer,getAtom(atomSer).getCode(),rignodemap.get(resser));
			graph.addVertex(atomNode); // this also adds the atomNode to the serials2nodes map
		}
		
		graph.setSecondaryStructure(secondaryStructureCopy);
		graph.setCutoff(cutoff);
		graph.setSequence(sequence.getSeq());
		graph.setPdbCode(getPdbCode());
		graph.setChainCode(chainCode);
		graph.setPdbChainCode(pdbChainCode);
		graph.setModel(parent.getModel());
		graph.setSid(sid);
		graph.setTargetNum(targetNum);
		graph.setGroupNum(groupNum);
		graph.setCaspModelNum(caspModelNum);
		graph.setMethodStr(caspMethodStr);
		graph.setAuthorStr(caspAuthorStr);
		graph.setParents(caspParents);
		graph.setCrossed(crossed);
		
		// populating the AIGraph with AIGEdges 
		for (int i=0;i<distMatrix.length;i++){
			for (int j=0;j<distMatrix[i].length;j++){
				// the condition distMatrix[i][j]!=0.0 takes care of skipping several things: 
				// - diagonal of the matrix in case of non-crossed
				// - lower half of matrix in case of non-crossed
				// - cells for which we didn't calculate a distance because the 2 points were not in same or neighbouring boxes (i.e. too far apart)
				if (distMatrix[i][j]!=0.0f && distMatrix[i][j]<=cutoff){
					if (!crossed) {
						graph.addEdge(new AIGEdge(distMatrix[i][j]), graph.getNodeFromSerial(iAtoms[i].getSerial()), graph.getNodeFromSerial(iAtoms[j].getSerial()), EdgeType.UNDIRECTED);
					}
					// This condition is to take care of crossed contact types that have overlapping sets of atoms: 
					//   the matrix would contain both i,j and j,i but that's only 1 edge in the AIGraph
					//TODO if our AIGraph didn't allow parallel edges, this extra check wouldn't be necessary
					else if (!graph.containsEdgeIJ(iAtoms[i].getSerial(), jAtoms[j].getSerial())) {
						graph.addEdge(new AIGEdge(distMatrix[i][j]), graph.getNodeFromSerial(iAtoms[i].getSerial()), graph.getNodeFromSerial(jAtoms[j].getSerial()), EdgeType.UNDIRECTED);
					}
				}

			}
		}

		return graph;
	}

	/**
	 * Returns a RIGraph for given contact type, cutoff and directionality
	 * Only standard aminoacids are considered.
	 * @param ct  the contact type
	 * @param cutoff  the distance cutoff
	 * @param directed  true if we want a directed graph, false for undirected
	 * @return
	 */
	public RIGraph getRIGraph(String ct, double cutoff, boolean directed) {
		//NOTE: At the moment we don't allow directed graphs for overlapping contact types e.g. directed ALL/BB
		//      because the code wouldn't be able to cope correctly with them
		// To lift this restriction, one possibility would be to make AIGraph directed and
		//- remove the check for adding parallel atomic edges if crossed
		//- add atomic edges in both directions if !crossed
		//- make sure in getRIGraph for undirected graphs not to double-count atomic edges
		if (directed && ContactType.isOverlapping(ct)) {
			throw new IllegalArgumentException("Contact type "+ct+" is overlapping. Generating directed graphs for it is unsupported");
		}
		
		String[] cts = ct.split("\\+");		
		AIGraph atomGraph = getAIGraph(cts[0], cutoff);
		for(int i=1; i< cts.length; i++) {
			atomGraph.addGraph(getAIGraph(cts[i], cutoff));
		}
		RIGraph graph = atomGraph.getRIGraph(directed);
		
		graph.setContactType(ct);
		graph.setCutoff(cutoff);

		return graph;
	}
	
	/**
	 * Returns a RIGraph for given contact type and cutoff
	 * Crossed contact types (i.e. those containing a "/") will be considered always as 
	 * directed. If one wants an undirected graph for a crossed contact type then method
	 * {@link #getRIGraph(String, double, boolean)} should be used instead.
	 * @param ct  the contact type
	 * @param cutoff  the distance cutoff
	 * @return
	 */
	public RIGraph getRIGraph(String ct, double cutoff) {
		// TODO eventually we should use the ContactType class as parameter so we could encapsulate all properties of a contact type in there
		boolean crossed = false;
		if (ct.contains("/")) {
			crossed = true;
		}		
		return getRIGraph(ct, cutoff, crossed);
	}
	
	/**
	 * Returns an all atom graph in a AIGraph object
	 * Only standard aminoacids are considered.
	 * @param cutoff  the distance cutoff
	 * @return
	 */
	public AIGraph getAllAtomGraph(double cutoff) {
		return this.getAIGraph("ALL", cutoff);
	}
	
	/**
	 * Calculates the grid density.
	 * Only standard aminoacids are considered.
	 * @param ct
	 * @param cutoff
	 * @param densityCount
	 */
	public void calcGridDensity(String ct, double cutoff, Map<Integer, Integer> densityCount) { 
		Atom[] atoms = getAtomsForCt(ct, null);

		Grid grid = new Grid(cutoff);
		grid.addAtoms(atoms,atoms);

		grid.countDensity(densityCount);
		
	}

	/**
	 * Calculates the difference distance map of this structure and 
	 * another pdb object given a sequence alignment of the structures. The 
	 * resulting difference distance map may contains non-defined distances. 
	 * This behavior is due to the alignment. If any residue in either 
	 * structures is aligned with a gap one cannot assign a "difference 
	 * distance" to another residue pair.   
	 * @param contactType1  contact type of this structure
	 * @param pdb2  the second structure
	 * @param contactType2  contact type of the second structure
	 * @param ali  sequence alignment of both structures
	 * @param name1  sequence tag of the this structure in the alignment
	 * @param name2  sequence tag og the second structure in the alignment
	 * @return the difference distance map
	 */
	public HashMap<Pair<Integer>,Double> getDiffDistMap(String contactType1, PdbChain pdb2, String contactType2, MultipleSequenceAlignment ali, String name1, String name2) {

		HashMap<Pair<Integer>,Double> otherDistMatrix = pdb2.calcDistMatrix(contactType2);
		HashMap<Pair<Integer>,Double> thisDistMatrix = this.calcDistMatrix(contactType1);
		HashMap<Pair<Integer>,Double> alignedDistMatrix = new HashMap<Pair<Integer>, Double>(Math.min(this.getFullLength(), pdb2.getFullLength()));
		int i1,i2,j1,j2;
		TreeSet<Integer> unobserved1 = new TreeSet<Integer>();
		TreeSet<Integer> unobserved2 = new TreeSet<Integer>();

		// detect all unobserved residues
		for(int i = 1; i <= ali.getAlignmentLength(); ++i) {
			i1 = ali.al2seq(name1, i);
			i2 = ali.al2seq(name2, i);
			if( i1 != -1 && !containsStdAaResidue(i1) ) {
				unobserved1.add(i1);
			}
			if( i2 != -1 && !pdb2.containsStdAaResidue(i2) ) {
				unobserved2.add(i2);
			}
		}

		// strategy: we always have to look through the alignment to say 
		// whether a difference distance can be assigned to a pair of 
		// corresponding residues. To put it differently, for any two 
		// alignment columns one always has to ensure that both columns 
		// only contain observed residues (no gaps!), otherwise the one 
		// cannot obtain a distance in at least one structure as a gap 
		// indicates "no coordinates available".  

		for(int i = 1; i <= ali.getAlignmentLength()-1; ++i) {

			i1 = ali.al2seq(name1, i);
			i2 = ali.al2seq(name2, i);

			// alignment columns must not contain gap characters and both 
			// residues in the current column have to be observed!
			if( i1 == -1 || i2 == -1 || unobserved1.contains(i1) || unobserved2.contains(i2) ) {
				continue;
			}

			for(int j = i + 1; j <= ali.getAlignmentLength(); ++j) {

				j1 = ali.al2seq(name1, j);
				j2 = ali.al2seq(name2, j);

				if( j1 == -1 || j2 == -1 || unobserved1.contains(j1) || unobserved2.contains(j2) ) {
					continue;
				}

				// make the edges
				Pair<Integer> e1 = new Pair<Integer>(i1,j1);
				Pair<Integer> e2 = new Pair<Integer>(i2,j2);

				alignedDistMatrix.put(new Pair<Integer>(i,j),Math.abs(thisDistMatrix.get(e1)-otherDistMatrix.get(e2)));
			}
		}
		return alignedDistMatrix;
	}
	// TODO: Version of this where already buffered distance matrices are passed as paremeters
	
	/**
	 * Computes the atom interaction graph between this and given protein chain for all
	 * atoms (standard aminoacids, hetatoms and nucleotides) and given cutoff. 
	 * A geometric hashing algorithm is used for fast contact computation (without needing 
	 * to calculate full distance matrix) 
	 * @param other
	 * @param cutoff
	 * @return a graph containing one edge per atom interaction i.e. per atom pair falling 
	 * under the distance cutoff
	 */
	public AICGraph getAICGraph(PdbChain other, double cutoff) {
		Atom[] thisAtoms = this.getAllAtoms();
		Atom[] otherAtoms = other.getAllAtoms();
		
		Grid grid = new Grid(cutoff);
		grid.addAtoms(thisAtoms,bounds,otherAtoms,other.bounds);
		
		AICGraph graph = new AICGraph();
		graph.setDistCutoff(cutoff);

		if (grid.isNoOverlap()) {
			return graph;
		}
		
		float[][] distMatrix = grid.getDistMatrix(true);
		
		for (int i=0;i<distMatrix.length;i++){ 
			for (int j=0;j<distMatrix[i].length;j++){
				// the condition distMatrix[i][j]!=0.0 takes care of skipping cells for which we 
				// didn't calculate a distance because the 2 points were not in same or neighbouring boxes (i.e. too far apart)
				if (distMatrix[i][j]!=0.0f && distMatrix[i][j]<=cutoff){
					graph.addEdge(new AICGEdge(distMatrix[i][j]), thisAtoms[i], otherAtoms[j], EdgeType.UNDIRECTED);
				}

			}
		}
		return graph;
	}

	/**
	 * Gets the all atoms rsa given an internal residue serial corresponding to a standard amino acid
	 * @param resser
	 * @return the rsa or null if rsa has not been calculated yet or the residue number cannot be found
	 * @throws NullPointerException if residue serial not present in this PdbChain instance 
	 */
	public Double getAllRsaFromResSerial(int resser){
		if (hasASA() && containsResidue(resser)) {
			return getResidue(resser).getRsa();
		} else {
			return null;
		}
	}

	/**
	 * Gets the sc rsa given an internal residue serial corresponding to a standard amino acid
	 * @param resser
	 * @return the sc rsa or null if rsa has not been calculated yet or the residue number cannot be found
	 * @throws NullPointerException if residue serial not present in this PdbChain instance 
	 */
	public Double getScRsaFromResSerial(int resser){
		if (hasASA() && containsResidue(resser)) {
			return getResidue(resser).getScRsa();
		} else {
			return null;
		}
	}

	/**
	 * Gets the internal residue serial (cif) given a pdb residue serial (author assignment)
	 * @param pdbresser
	 * @return the residue serial or -1 if no mapping exists
	 */
	public int getResSerFromPdbResSer (String pdbresser){ 
		if (pdbresser2resser.containsKey(pdbresser)) {
			return pdbresser2resser.get(pdbresser);
		} else {
			return -1;
		}
	}

	/**
	 * Gets the pdb residue serial (author assignment) given an internal residue serial (cif)
	 * @param resser
	 * @return the pdb residue serial or null if no mapping exists
	 */
	public String getPdbResSerFromResSer (int resser){
		if (resser2pdbresser.containsKey(resser)) {
			return this.resser2pdbresser.get(resser);
		} else {
			return null;
		}
	}

	/**
	 * Gets the residue serial given an atom serial
	 * @param atomser
	 * @return
	 */
	public int getResSerFromAtomSer(int atomser){
		return atomser2atom.get(atomser).getParentResSerial();
	}

	/**
	 * Gets the atom serial given the residue serial and atom code.
	 * The caller of this method needs to check whether the resser, atom and combination of the two exists
	 * using {@link #hasCoordinates(int, String)}. Otherwise, if this function
	 * is called and no atom serial exists, a null pointer exception will be thrown.
	 * @param resser
	 * @param atomCode
	 * @return the atom serial
	 * @throws NullPointerException if residue serial not present in this PdbChain 
	 * instance or if no atom of given type exists for given residue
	 */
	public int getAtomSerFromResSerAndAtom(int resser, String atomCode) {
		return getResidue(resser).getAtom(atomCode).getSerial();			
	}

	/**
	 * Returns the set of all atom serials for observed atoms of the given residue serial.
	 * @param resser the residue serial
	 * @return an ordered set of serials of the (observed) atoms in this residue
	 * @throws NullPointerException if given residue serial not present in this PdbChain
	 * instance
	 */
	public Set<Integer> getAtomSersFromResSer(int resser) {
		TreeSet<Integer> atomSers = new TreeSet<Integer>();
		for (Atom atom:this.getResidue(resser).getAtoms()) {
			atomSers.add(atom.getSerial());
		}
		return atomSers;
	}
	
	/**
	 * Checks whether this PdbChain is an all-atom one (disregarding Hydrogens), i.e. its standard amino acids
	 * have coordinates for most atoms (thus excludes CA-only or BB only or any chain with major group of atoms missing).
	 * Even PDB structures with all atoms can still have missing atoms for some residues, 
	 * here what we check is that the average number of (non-Hydrogen) atoms per standard amino acid is above the 
	 * threshold {@value #MIN_AVRG_NUM_ATOMS_RES} . This threshold has been obtained from statistics of a set of non-redundant 
	 * PDB structures.
	 * 
	 * @return true if above average atoms per residue threshold, false otherwise
	 */
	public boolean isAllAtom() {
		if (((double)this.getNumStdAaHeavyAtoms()/(double)this.getStdAaObsLength())<MIN_AVRG_NUM_ATOMS_RES) {
			return false;
		} else {
			return true;
		}
	}
	
	/**
	 * Checks whether the given residue serial has any associated coordinates
	 * for its atoms.
	 * @param resser  the residue serial
	 * @return true if there is at least one atom with coordinates, else false 
	 */
	public boolean hasCoordinates(int resser) {
		if (this.containsResidue(resser)) {
			if (this.getResidue(resser).getNumAtoms()>0) {
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}
	
	/**
	 * Checks whether the given atom type for the given residue serial has any associated coordinates
	 * @param resser the residue serial
	 * @param atomCode  
	 * @return true if there is coordinates for given resser-atom, else false 
	 */
	public boolean hasCoordinates(int resser, String atomCode) {
		if (!this.containsResidue(resser)) 
			return false;
		return getResidue(resser).containsAtom(atomCode);
	}
	
	/**
	 * Returns true if at least one Hydrogen atom is present in a standard amino-acid of this chain
	 * @see #removeHatoms()
	 * @return
	 */
	public boolean hasHydrogens() {
		for (Residue res:this) {
			if (res instanceof AaResidue) {
				for (Atom atom:res) {
					if (atom.getType()==AtomType.H) {
						return true;
					}
				}
			}
		}
		return false;
	}
	
	/**
	 * Returns true if this PdbChain has been restricted to a specific SCOP domain 
	 * @return
	 */
	public boolean isRestrictedToScopDomain() {
		return sid!=null;
	}
	
	/**
	 * Returns the sid of this PdbChain 
	 * It is set when restrictToScopDomain is run
	 */
	public String getSid() {
		return sid;
	}
	
	/**
	 * Gets the atom coordinates (Point3d object) given the atom serial
	 * @param atomser
	 * @return the atom coordinates or null if no atom exist for the atomser
	 */
	public Point3d getAtomCoord(int atomser) {
		Atom atom = getAtom(atomser);
		if (atom!=null){
			return atom.getCoords();
		}
		return null;
	}

	/**
	 * Gets the atom coordinates (a Point3d) given the residue serial and atom code 
	 * (standard PDB atom name, e.g. CA, N, C, O, CB, ...) 
	 * @param resser
	 * @param atomCode
	 * @return the coordinates or null if no such atomCode and resser combination exist
	 */
	public Point3d getAtomCoord(int resser, String atomCode) {
		if (containsResidue(resser)) {
			Residue residue = getResidue(resser);
			if (residue.containsAtom(atomCode)) {
				return residue.getAtom(atomCode).getCoords();
			}
		}
		return null;
	}
	
	/**
	 * Returns a Set of all ordered atom serials (only observed atoms)
	 * The order is according to atom serials (not necessarily ordered by residue
	 * although they almost always are)
	 * Atoms of all residues are considered: standard amino acids, het residues or 
	 * nucleotides (if chain is nucleic acid) 
	 * @return
	 */
	public Set<Integer> getAllAtomSerials() {
		return this.atomser2atom.keySet();
	}

	/**
	 * Returns a Set of all ordered residue serials (only observed residues)
	 * All residues are considered: standard amino acids, het residues or nucleotides (if chain
	 * is nucleic acid)
	 * @return
	 */
	public Set<Integer> getAllResSerials() {
		return residues.keySet();
	}
	
	/**
	 * Returns a Set of all ordered standard amino acid residue serials (only observed residues)
	 * @return
	 */
	public Set<Integer> getAllStdAaResSerials() {
		Set<Integer> set = new TreeSet<Integer>();
		for (Residue residue:this) {
			if (residue instanceof AaResidue) {
				set.add(residue.getSerial());
			}
		}
		return set;
	}
	
	/**
	 * Returns the first observed residue in the chain (can be standard aminoacid, nucleotide or het)
	 * @return
	 */
	public Residue getFirstResidue() {
		return residues.firstEntry().getValue();
	}
	
	/**
	 * Returns the last observed residue in the chain (can be standard aminoacid, nucleotide or het)
	 * @return
	 */
	public Residue getLastResidue() {
		return residues.lastEntry().getValue();
	}
	
	
	/**
	 * Returns the parent PdbAsymUnit to which this chain belongs
	 * @return
	 */
	public PdbAsymUnit getParent() {
		return this.parent;
	}
	
	public void setParent(PdbAsymUnit parent) {
		this.parent = parent;
	}
	
	/**
	 * Gets the 4 letter pdb code identifying this structure.
	 * If this chain doesn't have a parent PDB entry it will return PdbAsymUnit.NO_PDB_CODE
	 * @return
	 */
	public String getPdbCode() {
		if (parent==null) {
			return PdbAsymUnit.NO_PDB_CODE;
		}
		return this.parent.getPdbCode();
	}

	/**
	 * Gets the internal chain code (cif)
	 * @return
	 */
	public String getChainCode(){
		return this.chainCode;
	}
	
	/**
	 * Sets the internal chain code (cif)
	 * @param chainCode
	 */
	public void setChainCode(String chainCode) {
		this.chainCode = chainCode;
	}

	/**
	 * Gets the pdb chain code (author assignment code)
	 * @return
	 */
	public String getPdbChainCode(){
		return this.pdbChainCode;
	}
	
	/**
	 * Sets the pdb chain code (author assignment code) 
	 * @param pdbChainCode
	 */
	public void setPdbChainCode(String pdbChainCode) {
		this.pdbChainCode = pdbChainCode;
	}
	
	protected void setPdbresser2resserMap(TreeMap<String,Integer> pdbresser2resser){
		this.pdbresser2resser = pdbresser2resser;
	}
	
	protected void setResser2pdbresserMap(TreeMap<Integer,String> resser2pdbresser){
		this.resser2pdbresser = resser2pdbresser;
	}
	
	/**
	 * Tells whether this chain contains atoms with alt codes. 
	 * We don't store them but catch the case when we parse from PDB/mmCIF/pdbase
	 * @return
	 */
	public boolean hasAltCodes(){
		return hasAltCodes;
	}
	
	/**
	 * Sets the hasAltCodes field
	 * @param hasAltCodes
	 */
	protected void setHasAltCodes(boolean hasAltCodes) {
		this.hasAltCodes = hasAltCodes;
	}
	
	/**
	 * Gets the CASP target number
	 * @return
	 */
	public int getTargetNum() {
		return this.targetNum;
	}
	
	/**
	 * Sets the CASP target number
	 * @param targetNum
	 */
	public void setTargetNum(int targetNum) {
		this.targetNum = targetNum;
	}
	
	/**
	 * Gets the CASP model number
	 * @return
	 */
	public int getCaspModelNum() {
		return this.caspModelNum;
	}

	/**
	 * Sets the CASP model number
	 * @param caspModelNum
	 */
	public void setCaspModelNum(int caspModelNum) {
		this.caspModelNum = caspModelNum;
	}
	
	/**
	 * Sets the author string required for casp submissions, set to null to supress output.
	 * @param authorStr the author string to be set
	 */
	public void setCaspAuthorStr(String authorStr) {
		this.caspAuthorStr = authorStr;
	}

	/**
	 * Sets the method string required for casp submissions, set to null to supress output.
	 * @param methodStr the method string to be set
	 */
	public void setCaspMethodStr(String methodStr) {
		this.caspMethodStr = methodStr;
	}
	
	/**
	 * Sets the parents record which Casp models may optionally have to specify parents used in modelling the strucure.
	 * The value is used when writing to Casp TS files or by the getParents() method.
	 * @param parents an array of parent strings used in modelling this structure.
	 */
	public void setCaspParents(String[] parents) {
		this.caspParents = parents;
	}	
	
	/**
	 * Returns the list of parents which Casp models may optionally have.
	 * The value will be set when reading from Casp TS files or by the setParents method.
	 * @return an array of parent strings or null if no parents have been specified.
	 */
	public String[] getCaspParents() {
		return this.caspParents;
	}
	
	/**
	 * Gets the sequence
	 * @return
	 */
	public Sequence getSequence() {
		return sequence;
	}
	
	/**
	 * Gets the sequence as a String using Ms instead of Xs for MSEs (seleno-methionines)
	 * @return
	 */
	public String getSequenceMSEtoMET() {
		StringBuffer sb = new StringBuffer();
		for (int i=0;i<sequence.getLength();i++) {
			char letter = sequence.getSeq().charAt(i);
			Residue res = getResidue(i+1);
			if (letter=='X' && res!=null && res.getLongCode().equals("MSE")) {
				sb.append('M');
				continue; // we don't want to add the X in this case
			} 
			sb.append(letter);
		}
		
		return sb.toString();
	}
	
	/**
	 * Sets the sequence for the first time in this PdbChain. To use only from parser classes when there's no sequence set yet.
	 * @param sequence
	 * @param protein true if sequence is a protein, false if it is a nucleotide
	 */
	protected void setSequence(String seq, boolean protein) {
		this.sequence = new Sequence(getPdbCode()+this.pdbChainCode,seq,protein);
	}
	
	/**
	 * Sets the sequence for this PdbChain, overriding any current sequence information.
	 * If the observed residues (those having 3D coordinates) do not match the new sequence,
	 * an exception will be thrown.
	 * No coordinates or observed residues are modified.
	 * @param seq the new sequence
	 * @throws PdbLoadException if the given sequence does not match observed sequence from ATOM lines or if 
	 * the given sequence is of different type than existing one (nucleotide/protein)
	 */
	public void setSequence(Sequence seq) throws PdbLoadException {
		if ((seq.isNucleotide() && !this.sequence.isNucleotide()) || (seq.isProtein() && !this.sequence.isProtein())) {
			throw new PdbLoadException("Given sequence is of a different type ("+(seq.isProtein()?"protein":"nucleotide")+") than current sequence ("
					+(sequence.isProtein()?"protein":"nucleotide")+")");
		}
		// we check that the sequences from ATOM lines and the new sequence coincide (except for unobserved residues)
		for (int resser:getAllResSerials()) {
			Residue residue = this.getResidue(resser);
			char seqLetter = 0;
			if (residue instanceof AaResidue) {
				seqLetter = ((AaResidue)residue).getAaType().getOneLetterCode();
			} else if (residue instanceof NucResidue) {
				seqLetter = ((NucResidue)residue).getNucType().getOneLetterCode();
			} else if (residue instanceof HetResidue) {
				seqLetter = AminoAcid.XXX.getOneLetterCode();
			}
			if (seq.getSeq().charAt(resser-1)!=seqLetter) {
				throw new PdbLoadException("Given sequence does not match observed sequence from ATOM lines for position "+resser+".");
			}
		}
		this.sequence = seq;
	}

	
	/**
	 * Gets the observed sequence, i.e. the sequence as it appears in the ATOM 
	 * lines of the PDB file (observed non-standard amino acids will be shown as 'X')
	 * @return
	 */
	public String getObsSequence() {
		String obsSequence = "";
		for (Residue residue:this) {
			obsSequence+=residue.getShortCode();
		}
		return obsSequence;
	}
	
	/**
	 * Returns true if this chain is a non-polymer chain (purely HET residues), i.e. not a peptide or nucleotide chain
	 * @return
	 */
	public boolean isNonPolyChain() {
		return isNonPolyChain;
	}
	
	/**
	 * Sets the isNonPolyChain flag, true for non-polymers (purely HET residues chain) false for polymer (peptide/nucleotide) chains
	 * @param isNonPolyChain
	 */
	public void setIsNonPolyChain(boolean isNonPolyChain) {
		this.isNonPolyChain = isNonPolyChain;
	}
	
	/**
	 * Returns true if this PdbChain instance has been assigned atom b-factors 
	 * @return
	 */
	public boolean hasBfactors() {
		return hasBfactors;
	}
	
	/**
	 * Returns true if naccess has been run and thus ASA values are present
	 * in this PdbChain instance
	 * @return
	 */
	public boolean hasASA() {
		return hasASA;
	}

	/**
	 * Sets the hasASA flag. To be used when running Naccess and setting the ASA values 
	 * of residues.
	 * @param hasASA
	 */
	public void setHasASA(boolean hasASA) {
		this.hasASA = hasASA;
	}
	
	/** 
	 * Returns true if csa information is available, false otherwise. 
	 * @return
	 */
	public boolean hasCSA() {
		if (catalSiteSet == null) {
			return false;
		} else if (catalSiteSet.isEmpty()) {
			return false;
		} else {
			return true;
		}
	}
	
	/**
	 * Returns the csa annotation object of this PdbChain.
	 * @return
	 */
	public CatalSiteSet getCSA() {
		return catalSiteSet;
	}

	/** 
	 * Returns true if ec information is available, false otherwise. 
	 * @return
	 */
	public boolean hasEC() {
		if (ec == null) {
			return false;
		} else if (ec.isEmpty()) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 * Returns the ec annotation object of this graph.
	 * @return
	 */
	public EC getEC() {
		return ec;
	}

	/** 
	 * Returns true if scop information is available, false otherwise.
	 * @return 
	 */
	public boolean hasScop() {
		if (scop == null) {
			return false;
		} else if (scop.isEmpty()) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 * Returns the scop annotation object of this graph.
	 * @return
	 */
	public Scop getScop() {
		return scop;
	}

	// secondary structure related methods

	/** 
	 * Returns true if secondary structure information is available, false otherwise.
	 * @return 
	 */
	public boolean hasSecondaryStructure() {
		return !this.secondaryStructure.isEmpty();
	}

	/**
	 * Returns the secondary structure annotation object of this PdbChain.
	 * @return
	 */
	public SecondaryStructure getSecondaryStructure() {
		return this.secondaryStructure;
	}
	
	/**
	 * Sets the secondary structure annotation for this PdbChain, overwriting the existing one. 
	 * @param secondaryStructure
	 */
	public void setSecondaryStructure(SecondaryStructure secondaryStructure) {
		this.secondaryStructure = secondaryStructure;
		initialiseResiduesSecStruct();
	}

	// end of secondary structure related methods
	
	/**
	 * Sets the catalitic site set annotation object of this PdbChain
	 * @param catalSiteSet
	 */
	public void setCatalSiteSet(CatalSiteSet catalSiteSet) {
		this.catalSiteSet = catalSiteSet;
	}

	/**
	 * Sets the scop annotation object of this PdbChain
	 * @param scop
	 */
	public void setScop(Scop scop) {
		this.scop = scop;
	}
	
	/**
	 * Sets the EC annotation object of this PdbChain
	 * @param ec
	 */
	public void setEC(EC ec) {
		this.ec = ec;
	}
	
	/**
	 * Calculates rmsd (on atoms given by ct) of this PdbChain object to otherPdb object
	 * Both objects must represent structures with same sequence (save unobserved residues or missing atoms)
	 * 
	 * @param otherPdb
	 * @param ct the contact type (crossed contact types don't make sense here)
	 * @return
	 * @throws ConformationsNotSameSizeException
	 */
	public double rmsd(PdbChain otherPdb, String ct) throws ConformationsNotSameSizeException {
		return rmsd(otherPdb, ct, null);
	}
	
	/**
	 * Calculates rmsd (on atoms of given by contact type ct) of this PdbChain object to otherPdb object
	 * restricted to the given set of intervals
	 * Both objects must represent structures with same sequence (save unobserved residues or missing atoms)
	 * 
	 * @param otherPdb
	 * @param ct the contact type (crossed contact types don't make sense here)
	 * @param intervSet an interval set of residues for which the rmsd will be calculated, if 
	 * null then all residues will be considered
	 * @return
	 * @throws ConformationsNotSameSizeException
	 */
	public double rmsd(PdbChain otherPdb, String ct, IntervalSet intervSet) throws ConformationsNotSameSizeException {
		TreeMap<Integer, AaResidue> thisResidues = this.getReducedResidues(ct,intervSet);
		TreeMap<Integer, AaResidue> otherResidues = otherPdb.getReducedResidues(ct,intervSet);

		ArrayList<Vector3d> conf1AL = new ArrayList<Vector3d>();
		ArrayList<Vector3d> conf2AL = new ArrayList<Vector3d>();	
		// there might be unobserved residues or some missing atoms for a residue
		// here we get the ones that are in common
		for (int resser:thisResidues.keySet()) {
			AaResidue thisRes = thisResidues.get(resser);
			if (otherResidues.containsKey(resser)) {
				AaResidue otherRes = otherResidues.get(resser);
				for (Atom atom:thisRes.getAtoms()) {
					if (otherRes.containsAtom(atom.getCode())) {
						conf1AL.add(new Vector3d(atom.getCoords()));
						conf2AL.add(new Vector3d(otherRes.getAtom(atom.getCode()).getCoords()));
					}
				}
			}
		}

		// converting the ArrayLists to arrays
		Vector3d[] conformation1 = new Vector3d[conf1AL.size()]; 
		Vector3d[] conformation2 = new Vector3d[conf2AL.size()];
		conf1AL.toArray(conformation1);
		conf2AL.toArray(conformation2);
		
		// this as well as calculating the rmsd, changes conformation1 and conformation2 to be optimally superposed
		double rmsd = calculate_rmsd(conformation1, conformation2);

//		// printing out individual distances (conformation1 and conformation2 are now optimally superposed)
//		for (i=0;i<conformation1.length;i++){
//		Point3d point1 = new Point3d(conformation1[i].x,conformation1[i].y,conformation1[i].z);
//		Point3d point2 = new Point3d(conformation2[i].x,conformation2[i].y,conformation2[i].z);
//		System.out.println(point1.distance(point2));
//		}

		return rmsd;

	}

	/**
	 * Calculates the RMSD between two conformations.      
	 * conformation1: Vector3d array (matrix of dimensions [N,3])       
	 * conformation2: Vector3d array (matrix of dimensions [N,3]) 
	 * 
	 * Both conformation1 and conformation2 are modified to be optimally superposed
	 * 
	 * Implementation taken (python) from http://bosco.infogami.com/Root_Mean_Square_Deviation, 
	 * then ported to java using Jama matrix package 
	 * (page has moved to: http://boscoh.com/protein/rmsd-root-mean-square-deviation)                
	 * @param conformation1
	 * @param conformation2
	 * @return
	 * @throws ConformationsNotSameSizeException
	 */
	public static double calculate_rmsd(Vector3d[] conformation1, Vector3d[] conformation2) throws ConformationsNotSameSizeException{
		if (conformation1.length!=conformation2.length) {
			//System.err.println("Conformations not the same size");
			throw new ConformationsNotSameSizeException(
					"Given conformations have different size: conformation1: "+conformation1.length+", conformation2: "+conformation2.length);
		}
		int n_vec = conformation1.length;

		// 1st we bring both conformations to the same centre by subtracting their respectives centres
		Vector3d center1 = new Vector3d();
		Vector3d center2 = new Vector3d();
		for (int i=0;i<n_vec;i++){ // summing all vectors in each conformation
			center1.add(conformation1[i]);
			center2.add(conformation2[i]);
		}
		// dividing by n_vec (average)
		center1.scale((double)1/n_vec);
		center2.scale((double)1/n_vec);
		// translating our conformations to the same coordinate system by subtracting centers
		for (Vector3d vec:conformation1){
			vec.sub(center1);
		}
		for (Vector3d vec:conformation2){
			vec.sub(center2);
		}

		//E0: initial sum of squared lengths of both conformations
		double sum1 = 0.0;
		double sum2 = 0.0;
		for (int i=0;i<n_vec;i++){
			sum1 += conformation1[i].lengthSquared();
			sum2 += conformation2[i].lengthSquared();
		}
		double E0 = sum1 + sum2;

		// singular value decomposition
		Matrix vecs1 = vector3dAr2matrix(conformation1);
		Matrix vecs2 = vector3dAr2matrix(conformation2);

		Matrix correlation_matrix = vecs2.transpose().times(vecs1); //gives a 3x3 matrix

		SingularValueDecomposition svd = correlation_matrix.svd();
		Matrix U = svd.getU();
		Matrix V_trans = svd.getV().transpose(); 
		double[] singularValues = svd.getSingularValues();

		boolean is_reflection = false;
		if (U.det()*V_trans.det()<0.0){ 
			is_reflection = true;
		}
		if (is_reflection){
			// reflect along smallest principal axis:
			// we change sign of last coordinate (smallest singular value)
			singularValues[singularValues.length-1]=(-1)*singularValues[singularValues.length-1];  			
		}

		// getting sum of singular values
		double sumSV = 0.0;
		for (int i=0;i<singularValues.length;i++){
			sumSV += singularValues[i];
		}

		// rmsd square: Kabsch formula
		double rmsd_sq = (E0 - 2.0*sumSV)/((double) n_vec);
		rmsd_sq = Math.max(rmsd_sq, 0.0);

		// finally we modify conformation2 to be aligned to conformation1
		if (is_reflection) { // first we check if we are in is_reflection case: we need to change sign to last row of U
			for (int j=0;j<U.getColumnDimension();j++){
				// we change sign to last row of U
				int lastRow = U.getRowDimension()-1;
				U.set(lastRow, j, (-1)*U.get(lastRow,j));
			}
		}
		Matrix optimal_rotation = U.times(V_trans); 
		Matrix conf2 = vecs2.times(optimal_rotation);
		for (int i=0;i<n_vec;i++){
			conformation2[i].x = conf2.get(i,0);
			conformation2[i].y = conf2.get(i,1);
			conformation2[i].z = conf2.get(i,2);
		}

		return Math.sqrt(rmsd_sq);
	}

	/** Gets a Jama.Matrix object from a Vector3d[] (deep copies) */
	private static Matrix vector3dAr2matrix(Vector3d[] vecArray) {
		double[][] array = new double[vecArray.length][3];
		for (int i=0;i<vecArray.length;i++){
			vecArray[i].get(array[i]);
		}
		return new Matrix(array);
	}

	/**
	 * Write residue info to given db, using our db graph OWL format, 
	 * i.e. tables: residue_info
	 * Only considers standard amino acids
	 * @param conn
	 * @param db
	 * @throws SQLException
	 */
	public void writeToDb(MySQLConnection conn, String db) throws SQLException{

		Statement stmt;
		String sql = "";

		conn.setSqlMode("NO_UNSIGNED_SUBTRACTION,TRADITIONAL");

		for (Residue residue:this) {
			if (!(residue instanceof AaResidue)) continue;
			int resser = residue.getSerial();
			String resType = String.valueOf(((AaResidue)residue).getAaType().getOneLetterCode());
			String pdbresser = residue.getPdbSerial();

			String secStructType = null;
			String secStructId = null;
			if (secondaryStructure != null) {
				if (secondaryStructure.getSecStrucElement(resser) != null) {
					secStructType = quote(String.valueOf(secondaryStructure.getSecStrucElement(resser).getType()));
					secStructId = quote(secondaryStructure.getSecStrucElement(resser).getId());
				}
			}

			String scopId = null;
			String sccs = null;
			String sunid = null;
			String orderIn = null;
			String domainType = null;
			String domainNumReg = null;
			if (scop != null) {
				if (scop.getScopRegion(resser)!=null) {
					ScopRegion sr = scop.getScopRegion(resser);
					scopId = quote(sr.getSId());
					sccs = quote(sr.getSccs());
					sunid = String.valueOf(sr.getSunid());
					orderIn = String.valueOf(sr.getOrder());
					domainType = quote(String.valueOf(sr.getDomainType()));
					domainNumReg = String.valueOf(sr.getNumRegions());
				}
			}

			Double allRsa = getAllRsaFromResSerial(resser);
			Double scRsa = getScRsaFromResSerial(resser);

			// support for consurf parsing now discontinued: web site changed format, need a local executable or similar
			Double consurfhsspScore = null;//((AaResidue)residue).getConsurfScore();
			Integer consurfhsspColor = null;//((AaResidue)residue).getConsurfColor();

			String ecId = null;
			if (ec != null) {
				if (ec.getECRegion(resser) != null) {
					ecId = quote(ec.getECNum(resser));
				}
			}

			String csaNums = null;
			String csaChemFuncs = null;
			String csaEvids = null;
			if (catalSiteSet != null) {
				if (catalSiteSet.getCatalSite(resser) != null) {
					csaNums = quote(catalSiteSet.getCatalSiteNum(resser));
					csaChemFuncs = quote(catalSiteSet.getCatalSiteChemFunc(resser));
					csaEvids = quote(catalSiteSet.getCatalSiteEvid(resser));
				}
			}

			sql = "INSERT IGNORE INTO "+db+".pdb_residue_info (pdb_code, chain_code, pdb_chain_code, res_ser, pdb_res_ser, res_type, sstype, ssid, scop_id, sccs, sunid, order_in, domain_type, domain_num_reg, all_rsa, sc_rsa, consurfhssp_score, consurfhssp_color, ec, csa_site_nums, csa_chem_funcs, csa_evid) " +
			" VALUES ("+quote(getPdbCode())+", "+quote(chainCode)+", "+(pdbChainCode.equals(PdbAsymUnit.NULL_CHAIN_CODE)?quote("-"):quote(pdbChainCode))+","+resser+", "+quote(pdbresser)+", "+quote(resType)+", "+secStructType+", "+secStructId+", "+scopId+", "+sccs+", "+sunid+", "+orderIn+", "+domainType+", "+domainNumReg+", "+allRsa+", "+scRsa+", "+consurfhsspScore+","+consurfhsspColor+","+ecId+","+csaNums+","+csaChemFuncs+","+csaEvids+")";
			//System.out.println(sql);
			stmt = conn.createStatement();
			stmt.executeUpdate(sql);
			stmt.close();
		}			
	}

	/**
	 * Write residue info to given db, using our db graph OWL format, 
	 * i.e. tables: residue_info
	 * Only considers standard amino acids 
	 * @param conn
	 * @param db
	 * @throws SQLException
	 */
	public void writeToDbFast(MySQLConnection conn, String db) throws SQLException, IOException {

		Statement stmt;
		String sql = "";
		
		conn.setSqlMode("NO_UNSIGNED_SUBTRACTION,TRADITIONAL");
		
		PrintStream resOut = new PrintStream(new FileOutputStream(getPdbCode()+chainCode+"_residues.txt"));
		
		for (Residue residue:this) {
			if (!(residue instanceof AaResidue)) continue;
			int resser = residue.getSerial();
			String resType = String.valueOf(((AaResidue)residue).getAaType().getOneLetterCode());
			String pdbresser = residue.getPdbSerial();
			
			String secStructType = "\\N";
			String secStructId = "\\N";
			if (secondaryStructure != null) {
				if (secondaryStructure.getSecStrucElement(resser) != null) {
					secStructType = String.valueOf(secondaryStructure.getSecStrucElement(resser).getType());
					secStructId = secondaryStructure.getSecStrucElement(resser).getId();
				}
			}
			
			String scopId = "\\N";
			String sccs = "\\N";
			String sunid = "\\N";
			String orderIn = "\\N";
			String domainType = "\\N";
			String domainNumReg = "\\N";
			if (scop != null) {
				if (scop.getScopRegion(resser)!=null) {
					ScopRegion sr = scop.getScopRegion(resser);
					scopId = sr.getSId();
					sccs = sr.getSccs();
					sunid = String.valueOf(sr.getSunid());
					orderIn = String.valueOf(sr.getOrder());
					domainType = String.valueOf(sr.getDomainType());
					domainNumReg = String.valueOf(sr.getNumRegions());
				}
			}
			
			Double allRsa = getAllRsaFromResSerial(resser);
			Double scRsa = getScRsaFromResSerial(resser);
			
			// support for consurf parsing now discontinued: web site changed format, need a local executable or similar
			Double consurfhsspScore = null;//((AaResidue)residue).getConsurfScore();
			Integer consurfhsspColor = null;//((AaResidue)residue).getConsurfColor();
			
			String ecId = "\\N";
			if (ec != null) {
				if (ec.getECRegion(resser) != null) {
					ecId = ec.getECNum(resser);
				}
			}
			
			String csaNums = "\\N";
			String csaChemFuncs = "\\N";
			String csaEvids = "\\N";
			if (catalSiteSet != null) {
				if (catalSiteSet.getCatalSite(resser) != null) {
					csaNums = catalSiteSet.getCatalSiteNum(resser);
					csaChemFuncs = catalSiteSet.getCatalSiteChemFunc(resser);
					csaEvids = catalSiteSet.getCatalSiteEvid(resser);
				}
			}
			
			resOut.println(getPdbCode()+"\t"+chainCode+"\t"+(pdbChainCode.equals(PdbAsymUnit.NULL_CHAIN_CODE)?"-":pdbChainCode)+"\t"+resser+"\t"+pdbresser+"\t"+resType+"\t"+secStructType+"\t"+secStructId+"\t"+scopId+"\t"+sccs+"\t"+sunid+"\t"+orderIn+"\t"+domainType+"\t"+domainNumReg+"\t"+allRsa+"\t"+scRsa+"\t"+consurfhsspScore+"\t"+consurfhsspColor+"\t"+ecId+"\t"+csaNums+"\t"+csaChemFuncs+"\t"+csaEvids);
			
		}
		resOut.close();
		sql = "LOAD DATA LOCAL INFILE '"+getPdbCode()+chainCode+"_residues.txt' INTO TABLE "+db+".pdb_residue_info (pdb_code, chain_code, pdb_chain_code, res_ser, pdb_res_ser, res_type, sstype, ssid, scop_id, sccs, sunid, order_in, domain_type, domain_num_reg, all_rsa, sc_rsa, consurfhssp_score, consurfhssp_color, ec, csa_site_nums, csa_chem_funcs, csa_evid);";
		//System.out.println(sql);
		stmt = conn.createStatement();
		stmt.executeUpdate(sql);
		stmt.close();
		File fileToDelete = new File(getPdbCode()+chainCode+"_residues.txt");
		if (fileToDelete.exists()) {
			fileToDelete.delete();
		}
	}

	private static String quote(String s) {
		return ("'"+s+"'");
	}
	
	/**
	 * Restricts thisPdb object to residues that belong to the given sunid
	 * Can only be used if SCOP annotation is loaded. Check it with {@link #hasScop()}
	 * @param sunid
	 */
	public void restrictToScopDomain (int sunid) {
		Vector<ScopRegion> scopRegions = this.scop.getScopRegions(sunid);
		if (scopRegions.size()!=0) {
			this.sid = scopRegions.get(0).getSId();
			if (scopRegions.get(0).getDomainType() != ScopRegion.DomainType.WHOLECHAIN) {
				restrictToScopRegions(scopRegions);
				
				Iterator<ScopRegion> allScopRegions = this.scop.getIterator();
				while (allScopRegions.hasNext()) {
					ScopRegion scopRegion = allScopRegions.next();
					if (!scopRegion.getSId().equals(sid)) {
						allScopRegions.remove();
					}					
				}
			}
		}
	}
	
	/**
	 * Restricts thisPdb object to residues that belong to the given sid
	 * Can only be used if SCOP annotation is loaded. Check it with {@link #hasScop()} 
	 * @param sid
	 */
	public void restrictToScopDomain (String sid) {

		Vector<ScopRegion> scopRegions = this.scop.getScopRegions(sid);
		if (scopRegions.size()!=0) {
			this.sid = sid;
			if (scopRegions.get(0).getDomainType() != ScopRegion.DomainType.WHOLECHAIN) {
				restrictToScopRegions(scopRegions);
				
				Iterator<ScopRegion> allScopRegions = this.scop.getIterator();
				while (allScopRegions.hasNext()) {
					ScopRegion scopRegion = allScopRegions.next();
					if (!scopRegion.getSId().equals(sid)) {
						allScopRegions.remove();
					}					
				}
			}
		}
	}
	
	/**
	 * Restricts thisPdb object to residues that belong to the given sids
	 * Can only be used after calling checkScop() 
	 * @param sid
	 */
	public void restrictToScopDomains (String[] sids) {
		
		TreeSet<String> keepSids = new TreeSet<String>();
		
		Vector<ScopRegion> scopRegions = new Vector<ScopRegion>();
		for (String sid: sids) {
			if (!keepSids.contains(sid)) {
				Vector<ScopRegion> curScopRegions = this.scop.getScopRegions(sid);
				if (curScopRegions.size() != 0) {
					scopRegions.addAll(curScopRegions);
					keepSids.add(sid);
				}
			}
		}
		
		if (scopRegions.size()!=0) {
			Iterator<String> it = keepSids.iterator();
			this.sid = it.next();
			while (it.hasNext()) {
				this.sid += ","+it.next(); 
			}
			
			if (scopRegions.get(0).getDomainType() != ScopRegion.DomainType.WHOLECHAIN) {
				restrictToScopRegions(scopRegions);
				
				Iterator<ScopRegion> allScopRegions = this.scop.getIterator();
				while (allScopRegions.hasNext()) {
					ScopRegion scopRegion = allScopRegions.next();
					if (!keepSids.contains(scopRegion.getSId())) {
						allScopRegions.remove();
					}					
				}
			}
		}
	}
	
	/**
	 * Restricts thisPdb object to residues that belong to the given sunids
	 * Can only be used after calling checkScop() 
	 * @param sid
	 */
	public void restrictToScopDomains (int[] sunids) {
		
		TreeSet<String> keepSids = new TreeSet<String>();
		TreeSet<Integer> keepSunids = new TreeSet<Integer>();
		
		Vector<ScopRegion> scopRegions = new Vector<ScopRegion>();
		for (int sunid: sunids) {
			if (!keepSunids.contains(sunid)) {
				Vector<ScopRegion> curScopRegions = this.scop.getScopRegions(sunid);
				if (curScopRegions.size() != 0) {
					scopRegions.addAll(curScopRegions);
					keepSunids.add(sunid);
					keepSids.add(curScopRegions.get(0).getSId());
				}
			}
		}
		
		if (scopRegions.size()!=0) {
			Iterator<String> it = keepSids.iterator();
			this.sid = it.next();
			while (it.hasNext()) {
				this.sid += ","+it.next(); 
			}
			
			if (scopRegions.get(0).getDomainType() != ScopRegion.DomainType.WHOLECHAIN) {
				restrictToScopRegions(scopRegions);
				
				Iterator<ScopRegion> allScopRegions = this.scop.getIterator();
				while (allScopRegions.hasNext()) {
					ScopRegion scopRegion = allScopRegions.next();
					if (!keepSids.contains(scopRegion.getSId())) {
						allScopRegions.remove();
					}					
				}
			}
		}
	}
	
	/**
	 * Restricts this PdbChain object to residues within the given ScopRegions 
	 * @param scopRegions
	 */
	private void restrictToScopRegions (Vector<ScopRegion> scopRegions) {
		IntervalSet intervSet = new IntervalSet();
		Iterator<ScopRegion> it = scopRegions.iterator();
		while(it.hasNext()) {
			ScopRegion scopRegion = it.next();
			intervSet.add(scopRegion.getInterval());
		}
		restrictToIntervalSet(intervSet);
	}
	
	/**
	 * Restricts this PdbChain object to residues within the given IntervalSet
	 * @param intervSet a set of internal residue serials
	 */
	public void restrictToIntervalSet(IntervalSet intervSet) {
		
		// getting list of the residue serials to keep
		TreeSet<Integer> resSersToKeep = intervSet.getIntegerSet();

		// removing residues
		Iterator<Residue> resIt = iterator();
		while (resIt.hasNext()) {
			int resser = resIt.next().getSerial();
			if (!resSersToKeep.contains(resser)) {
				resIt.remove();
				if (catalSiteSet != null) {
					catalSiteSet.removeCatalSiteRes(resser);
				}
			}
		}
		// reinitialise the maps
		initialiseMaps();
			
		// setting sequence to scop sequence and obsLength and fullLength respectively
		Iterator<Interval> regionsToKeep = intervSet.iterator();
		String newSequence = "";
		while (regionsToKeep.hasNext()) {
			Interval region = regionsToKeep.next();
			newSequence += sequence.getSeq().substring((region.beg-1),region.end);
		}
		sequence = new Sequence(sequence.getName()+"_partial",newSequence);
	}
	
	/**
	 * Returns the number of residues with BSA (buried surface accessibility) above
	 * given BSA value. 
	 * @param value
	 * @return
	 */
	public int getNumResiduesWithBsaAbove(double value) {
		int count=0;
		for (Residue residue:this) {
			if ((residue instanceof AaResidue)) {
				if (residue.getBsa()>value) {
					count++;
				}
			}
		}
		return count;
	}
	
	/**
	 * Returns 2 lists of residues as a {@link InterfaceRimCore} object: core residues are those 
	 * with bsa>0 and with bsa/asa ratio above the given cut-off, rim those with bsa>0 and with 
	 * bsa/asa below the cut-off.
	 * Additionally residues will be considered at interface surface only if their ASA is above 
	 * given minAsaForSurface 
	 * @param bsaToAsaCutoff
	 * @param minAsaForSurface
	 * @return
	 */
	public InterfaceRimCore getRimAndCore(double bsaToAsaCutoff, double minAsaForSurface) {
		List<Residue> core = new ArrayList<Residue>();
		List<Residue> rim = new ArrayList<Residue>();
		for (Residue residue:this) {
			if (residue.getAsa()>minAsaForSurface && residue.getBsa()>0) {
				if (residue.getBsaToAsaRatio()<bsaToAsaCutoff) {
					rim.add(residue);
				} else {
					core.add(residue);
				}
			}
		}
		return new InterfaceRimCore(rim,core,bsaToAsaCutoff,minAsaForSurface);
	}
	
	/**
	 * Returns 2 list of residues as a {@link InterfaceRimCore} object (see {@link #getRimAndCore(double)})
	 * The core is required to have a minimum of minNumResidues. If the minimum is not 
	 * reached with the bsaToAsaSoftCutoff, then the cutoff is relaxed in relaxationStep steps 
	 * until reaching the bsaToAsaHardCutoff.
	 * Residues will be considered at interface surface only if their ASA is above 
	 * given minAsaForSurface  
	 * @param bsaToAsaSoftCutoff
	 * @param bsaToAsaHardCutoff
	 * @param relaxationStep
	 * @param minNumResidues
	 * @param minAsaForSurface
	 * @return
	 */
	public InterfaceRimCore getRimAndCore(double bsaToAsaSoftCutoff, double bsaToAsaHardCutoff, double relaxationStep, int minNumResidues, double minAsaForSurface) {
		InterfaceRimCore rimCore = null;
		// we introduce a margin of relaxationSte*0.10 to be sure we do go all the way down to bsaToAsaHardCutoff (necessary because of rounding)
		for (double cutoff=bsaToAsaSoftCutoff;cutoff>=bsaToAsaHardCutoff-relaxationStep*0.10;cutoff-=relaxationStep) {
			rimCore = getRimAndCore(cutoff, minAsaForSurface);
			//System.out.printf("cutoff %4.2f, core size: %d\n",cutoff,rimCore.getCoreSize());
			if (rimCore.getCoreSize()>=minNumResidues) {
				break;
			}
		}
		return rimCore;
	}
	
	/**
	 * Returns 2 lists of residues as a {@link InterfaceRimCore} object.
	 * Following the Chakrabarti definition (see Chakrabarti, Janin Proteins 2002):
	 * core residues have at least one atom fully buried (bsa/asa=1) and rim residues are all the
	 * rest still with bsa>0 but bsa/asa<1 (all atoms partially accessible)
	 * Residues will be considered at interface surface only if their ASA is above 
	 * given minAsaForSurface   
	 * @param minAsaForSurface
	 * @return
	 */
	public InterfaceRimCore getRimAndCoreChakrabarti(double minAsaForSurface) {
		InterfaceRimCore rimcore = new InterfaceRimCore();
		for (Residue residue:this) {
			if (residue.getAsa()>minAsaForSurface && residue.getBsa()>0) {
				boolean iscore = false;
				for (Atom atom:residue) {
					if (atom.getBsa()/atom.getAsa()>0.999) {
						iscore = true;
						break;
					}
				}
				if (iscore) rimcore.addCoreResidue(residue);
				else rimcore.addRimResidue(residue);
			}
		}
		return rimcore;
	}
	
	/**
	 * Returns 2 lists of residues as a {@link InterfaceRimCore} object. 
	 * Following the Levy definition (see Levy JMB 2010):
	 * core residues have rASA(c) below 25% while rASA(u) was above 25%, 
	 * rim residues have rASA(c) above 25%
	 * The support residues, a third class of residues introduced in Levy JMB 2010, are not included
	 * in either the core or the rim sets here. Support residues are those with rASA(u)<25%  
	 * @param rASAcutoff
	 * @return
	 */
	public InterfaceRimCore getRimAndCoreLevy(double rASAcutoff) {
		InterfaceRimCore rimcore = new InterfaceRimCore();
		for (Residue residue:this) {
			if (residue.getBsa()>0) {
				if (residue instanceof AaResidue) {
					AaResidue aares = (AaResidue) residue;
					// rBSA = rASA(u)-rASA(c) --> rASA(c) = rASA(u) - rBSA
					double rbsa = aares.getBsa()/aares.getAaType().getAsaInExtTripept();
					double rasau = aares.getAsa()/aares.getAaType().getAsaInExtTripept();
					double rasac = rasau -rbsa;
					if ( rasac<rASAcutoff && rasau>rASAcutoff) {  
						rimcore.addCoreResidue(residue);
					} else {
						rimcore.addRimResidue(residue);
					}
				}
			}
		}
		return rimcore;
	}
	
	/**
	 * Returns a list of all surface residues, i.e. all residues whose 
	 * ASA is above the given minAsaForSurface
	 * @param minAsaForSurface
	 * @return
	 */
	public List<Residue> getSurfaceResidues(double minAsaForSurface) {
		List<Residue> surfResidues = new ArrayList<Residue>();
		for (Residue res:this) {
			if (res.getAsa()>minAsaForSurface) surfResidues.add(res);
		}
		return surfResidues;
	}
	
	/**
	 * Mirror this PdbChain structure by inverting through the origin.
	 */
	public void mirror() {
		this.bounds = null; // we must reset bounds whenever the coordinates are changed
		for (Residue residue:this) {
			for (Atom atom:residue) {
				Point3d coords = atom.getCoords();
				coords.x *= -1;
				coords.y *= -1;
				coords.z *= -1;				
			}
		}
	}

	/**
	 * Transforms (rotation+translation) this structure in place as indicated by the given matrix. 
	 * @param m the rotation/translation matrix
	 */
	public void transform(Matrix4d m) {	
		this.bounds = null; // we must reset bounds whenever the coordinates are changed
		for (Residue residue:this) {
			for (Atom atom:residue) {
				Point3d coords = atom.getCoords();
				m.transform(coords);
			}
		}
	}
	
	/**
	 * Transforms (rotation+translation) this structure as indicate by the given matrix placing the
	 * result into the given pdb parameter (its existing contents will be wiped out). 
	 * @param m the rotation/translation matrix
	 * @param pdb
	 */
	public void transform(Matrix4d m, PdbChain pdb) {
		pdb = this.copy(this.parent);
		pdb.bounds = null;
		pdb.transform(m);
	}
	
	/**
	 * Translates this PdbChain to the given unit cell (direction).
	 * e.g. doCrystalTranslation(new Point3i(1,1,1)) will translate this PdbChain to 
	 * crystal cell (1,1,1), considering always this PdbChain's cell to be (0,0,0)
	 * The bounds are translated too so there is no need to recalculate them later
	 * @param direction
	 */
	public void doCrystalTranslation(Point3i direction) {
		Matrix4d m = parent.getCrystalCell().getTransform(direction);
		if (bounds!=null) {
			bounds.translate(new Vector3d(m.m03,m.m13,m.m23));
		}
		for (Residue residue:this) {
			for (Atom atom:residue) {
				Point3d coords = atom.getCoords();
				m.transform(coords);
			}
		}
	}
	
	/**
	 * Tells wheter given PdbChain's centre of mass is closer to this' centre of mass than twice 
	 * the maximum diagonal dimension of the unit cell.
	 * @param pdb
	 * @return
	 */
	public boolean isCrystalNeighbor(PdbChain pdb) {
		Point3d thisCoM  = this.getCentroid();
		Point3d otherCoM = pdb.getCentroid();
		if (thisCoM.distance(otherCoM)<2.0*parent.getCrystalCell().getMaxDimension()) {
			return true;
		}
		return false;
	}
	
	/**
	 * Returns an integer triplet indicating the crystal coordinates of this PdbChain with respect to the given one.
	 * e.g. if given PdbChain is a translation of this to adjacent cell (-1,0,0) then the triplet (-1,0,0) is returned  
	 * @param pdb
	 * @return
	 */
	public Point3i getCrystalSeparation(PdbChain pdb) {
		Point3d thisCoM  = this.getCentroid();
		Point3d otherCoM = pdb.getCentroid();
		parent.getCrystalCell().getCrystalFromOrthCoords(thisCoM);
		parent.getCrystalCell().getCrystalFromOrthCoords(otherCoM);
		double asep = otherCoM.x-thisCoM.x;
		double bsep = otherCoM.y-thisCoM.y;
		double csep = otherCoM.z-thisCoM.z;
		//System.out.printf("a: %5.2f ",asep);
		//System.out.printf("b: %5.2f ",bsep);
		//System.out.printf("c: %5.2f \n",csep);
		return new Point3i((int)asep, (int)bsep, (int)csep);
	}
	
	/**
	 * Transform the coordinates of this structure translating them to the given center
	 * and rotating them so that the given axis aligns with the z-axis
	 * @param center
	 * @param axis
	 */
	public void transformToCenterAndAxis(Point3d center, Vector3d axis) {
		this.bounds = null; // we must reset bounds whenever the coordinates are changed
		// finding the rotation matrix to align z axis to the given inertia axis
		Vector3d r = new Vector3d();
		Vector3d k = new Vector3d(0,0,1);
		r.cross(axis, k); // this is the axis of rotation
		double alpha = axis.angle(k); // this is the angle to rotate
		AxisAngle4d axisAngle = new AxisAngle4d(r, alpha);
		// note that the matrix needs to be initialised to the unit matrix otherwise setRotation() doesn't work properly
		Matrix4d rot = new Matrix4d(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1); 
		rot.setRotation(axisAngle);
		for (Residue residue: this) {
			for (Atom atom:residue) {
				Point3d coords = atom.getCoords();
				// translate to new origin
				coords.sub(center);
				// rotate so that z axis is the given axis
				rot.transform(coords);
			}
		}
	}
	
	/**
	 * Rotates this structure around rotAxis with the given rotAngle 
	 * @param rotAxis the vector around which the rotation will be performed
	 * @param rotAngle the rotation angle in radians
	 */
	public void rotate(Vector3d rotAxis, double rotAngle) {
		this.bounds = null; // we must reset bounds whenever the coordinates are changed
		AxisAngle4d axisAngle = new AxisAngle4d(rotAxis, rotAngle);
		// note that the matrix needs to be initialised to the unit matrix otherwise setRotation() doesn't work properly
		Matrix4d rot = new Matrix4d(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1); 
		rot.setRotation(axisAngle);
		for (Residue residue: this) {
			for (Atom atom:residue) {
				rot.transform(atom.getCoords());
			}
		}
	}

	
	/**
	 * Moves this structure such that the center of mass is at the origin using all atoms
	 */
	public void moveToOrigin() {
		this.bounds = null; // we must reset bounds whenever the coordinates are changed
		Point3d center = getCentroid();
		for (Residue residue: this) {
			for (Atom atom:residue) {
				atom.getCoords().sub(center);
			}
		}
	}
	
	/**
	 * Moves this structure such that the center of mass of the given subset of residues is at the origin using only CA atoms
	 */
	public void moveToOrigin(TreeSet<Integer> residues) {
		this.bounds = null; // we must reset bounds whenever the coordinates are changed
		Vector3d sumVector = new Vector3d();
		int numVectors = 0;
		for(int resser:residues) {
			Point3d coords = getAtomCoord(resser, "CA");
			sumVector.add(coords);
			numVectors++;			
		}

		sumVector.scale(1.0/numVectors);
		//System.out.println(sumVector);
		
		for (Residue residue: this) {
			for (Atom atom:residue) {
				atom.getCoords().sub(sumVector);
			}
		}
	}
	
	/**
	 * Returns the centroid (average coordinate of all atoms) of this chain 
	 * (i.e. center of mass, but disregarding atom masses)
	 * Atoms of all kinds of residues are considered
	 * @return
	 */
	public Point3d getCentroid() {
		Vector3d sumVector = new Vector3d();
		int numAtoms = 0;
		for (Residue residue:this) {
			for (Atom atom:residue) {
				Point3d coords = atom.getCoords();
				sumVector.add(coords);
				numAtoms++;
			}
		}
		sumVector.scale(1.0/numAtoms);
		return new Point3d(sumVector);
	}
	
	/**
	 * Calculates molecular weight of this PDB chain for all atoms of all residues
	 * @return
	 */
	public double getMass() {
		double mass = 0;
		for (Residue res:this) {
			for (Atom atom:res) {
				mass+=atom.getType().getAtomicMass();	
			}
		}
		return mass;
	}
	
	/**
	 * Gets the phi angle in degrees for given residue serial
	 * @param i
	 * @return the phi angle or NaN if there are no coordinates for given i or i-1
	 */
	public double getPhiAngle(int i) { 
		if (!hasCoordinates(i-1, "C") || !hasCoordinates(i, "N") || !hasCoordinates(i, "CA") || !hasCoordinates(i, "C")) {
			return Double.NaN;
		}
		Point3d Ciminus1 = getAtomCoord(i-1, "C");
		Point3d Ni = getAtomCoord(i, "N");
		Point3d CAi = getAtomCoord(i, "CA");
		Point3d Ci = getAtomCoord(i, "C");
		return getTorsionAngle(Ciminus1, Ni, CAi, Ci);	
	}
	
	/**
	 * Gets the psi angle in degrees for the given residue serial
	 * @param i
	 * @return the psi angle or NaN if there are no coordinates for given i or i+1
	 */
	public double getPsiAngle(int i) {
		if (!hasCoordinates(i, "N") || !hasCoordinates(i, "CA") || !hasCoordinates(i, "C") || !hasCoordinates(i+1,"N")) {
			return Double.NaN;
		}
		Point3d Ni = getAtomCoord(i, "N");
		Point3d CAi = getAtomCoord(i, "CA");
		Point3d Ci = getAtomCoord(i, "C");
		Point3d Niplus1 = getAtomCoord(i+1, "N");
		return getTorsionAngle(Ni, CAi, Ci, Niplus1);
	}
	
	/**
	 * Gets the omega angle in degrees for the given residue serial
	 * @param i
	 * @return the omega angle or NaN if there are no coordinates for given i or i+1
	 */
	public double getOmegaAngle(int i) {
		if (!hasCoordinates(i, "CA") || !hasCoordinates(i, "C") || !hasCoordinates(i+1, "N") || !hasCoordinates(i+1,"CA")) {
			return Double.NaN;
		}
		Point3d CAi = getAtomCoord(i, "CA");
		Point3d Ci = getAtomCoord(i, "C");
		Point3d Niplus1 = getAtomCoord(i+1, "N");
		Point3d CAiplus1 = getAtomCoord(i+1, "CA");
		return getTorsionAngle(CAi, Ci, Niplus1, CAiplus1);
	}
	
	/**
	 * Gets the torsion angle in degrees for the 4 given atoms 
	 * See http://en.wikipedia.org/wiki/Dihedral_angle for how to calculate it.
	 * @param atom1
	 * @param atom2
	 * @param atom3
	 * @param atom4
	 * @return
	 */
	private double getTorsionAngle(Point3d atom1, Point3d atom2, Point3d atom3, Point3d atom4) {
		Vector3d b1 = new Vector3d();
		Vector3d b2 = new Vector3d();
		Vector3d b3 = new Vector3d();
		Vector3d b1xb2 = new Vector3d();
		Vector3d b2xb3 = new Vector3d();
		Vector3d b1scalelengthb2 = new Vector3d();
		b1.sub(atom2, atom1);
		b2.sub(atom3, atom2);
		b3.sub(atom4, atom3);
		b1xb2.cross(b1, b2);
		b2xb3.cross(b2, b3);
		b1scalelengthb2.scale(b2.length(), b1);
		return Math.toDegrees(Math.atan2(b1scalelengthb2.dot(b2xb3), b1xb2.dot(b2xb3)));		
	}
	
	/**
	 * Gets all phi/psi angles in degrees for this PdbChain structure 
	 * @return residue serials as keys, values arrays of 2 doubles: first phi angle, second psi angle 
	 */
	public TreeMap<Integer, double[]> getAllPhiPsi() {
		TreeMap<Integer, double[]> phipsi = new TreeMap<Integer, double[]>();
		for (int resser: getAllResSerials()) {
			double[] angles = {getPhiAngle(resser), getPsiAngle(resser)};
			phipsi.put(resser, angles);
		}
		return phipsi;
	}
	
	/**
	 * Prints in 2 columns all phi, psi angles for this protein structure.
	 * Use to generate a Ramachandran plot.
	 */
	public void printAllPhiPsi() {
		TreeMap<Integer, double[]> phipsi = getAllPhiPsi();
		for (int resser:phipsi.keySet()) {
			System.out.printf("%7.2f %7.2f\n",phipsi.get(resser)[0],phipsi.get(resser)[1]);
		}
	}
	
	/**
	 * Returns a TreeMap of residue serials to standard amino acid unobserved residue's one letter codes 
	 * @return
	 */
	public TreeMap<Integer,Character> getUnobservedResidues() {
		TreeMap<Integer,Character> unobserved = new TreeMap<Integer,Character>();

		// detect all unobserved residues
		for(int i = 1; i <= getFullLength(); ++i) {
			if(!residues.containsKey(i)) {
				if (sequence.getSeq().charAt(i-1)!=AminoAcid.XXX.getOneLetterCode()) {
					unobserved.put(i,sequence.getSeq().charAt(i-1));
				}
			}
		}
		return unobserved;
		
	}

	/**
	 * Deep copies this PdbChain object
	 * @return
	 */
	public PdbChain copy(PdbAsymUnit parent) {
		PdbChain newPdb = new PdbChain();
		newPdb.parent = parent;
		newPdb.pdbChainCode = this.pdbChainCode;
		newPdb.chainCode = this.chainCode;
		newPdb.sid = this.sid;
		newPdb.targetNum = this.targetNum;
		newPdb.caspModelNum = this.caspModelNum;
		newPdb.caspAuthorStr = this.caspAuthorStr;
		newPdb.caspMethodStr = this.caspMethodStr;
		if (this.caspParents!=null) {
			newPdb.caspParents = new String[this.caspParents.length];
			for (int i=0;i<newPdb.caspParents.length;i++){
				newPdb.caspParents[i] = this.caspParents[i];
			}
		}
		newPdb.groupNum = this.groupNum;
		newPdb.hasASA = this.hasASA;
		newPdb.hasBfactors = this.hasBfactors;
		newPdb.isNonPolyChain = this.isNonPolyChain;
		newPdb.hasAltCodes = this.hasAltCodes;
				
		newPdb.sequence = this.sequence;
		if (secondaryStructure!=null) {
			newPdb.secondaryStructure = this.secondaryStructure.copy();
		}
		if (this.scop!=null) newPdb.scop = this.scop.copy();
		if (this.ec!=null) newPdb.ec = this.ec.copy();
		if (this.catalSiteSet!=null) newPdb.catalSiteSet = this.catalSiteSet.copy();
		// residues
		newPdb.residues = new TreeMap<Integer, Residue>();
		for (int resser:this.residues.keySet()) {
			newPdb.residues.put(resser,this.residues.get(resser).copy(newPdb));
		}
		// atomser2atom
		newPdb.atomser2atom = new TreeMap<Integer, Atom>();
		for (Residue residue:newPdb.residues.values()) {
			for (Atom atom:residue.getAtoms()) {
				newPdb.atomser2atom.put(atom.getSerial(), atom);
			}
		}
		// resser2pdbresser and pdbresser2resser
		if (resser2pdbresser!=null) {
			newPdb.resser2pdbresser = new TreeMap<Integer, String>();
			newPdb.pdbresser2resser = new TreeMap<String, Integer>();
			for (int resser:this.resser2pdbresser.keySet()) {
				newPdb.resser2pdbresser.put(resser, this.resser2pdbresser.get(resser));
			}
			for (String pdbresser:this.pdbresser2resser.keySet()) {
				newPdb.pdbresser2resser.put(pdbresser, this.pdbresser2resser.get(pdbresser));
			}
		}
		
		if (this.bounds!=null) {
			newPdb.bounds = new BoundingBox(bounds.xmin, bounds.xmax, bounds.ymin, bounds.ymax, bounds.zmin, bounds.zmax);
		}
		
		return newPdb;
	}
	
	

	/*--------------------------------------- static methods -----------------------------------------*/
	
	/**
	 * Loads a pdb structure where arg can be a pdbcode+chaincode or a pdb file name.
	 * If something goes wrong, prints an error message and exits.
	 * @param arg a pdbcode+chaincode (e.g. 1tdrB) or a pdb file name
	 * @return the newly created pdb object
	 */
	public static PdbChain readStructureOrExit(String arg) {
		return readFromFileOrPdbCode(arg, true, true);
	}
	
	/**
	 * Loads a pdb structure where arg can be a pdbcode+chaincode or a pdb file name.
	 * If something goes wrong, prints an error message and returns null;
	 * @param arg a pdbcode+chaincode (e.g. 1tdrB) or a pdb file name
	 * @return the newly created pdb object or null
	 */	
	public static PdbChain readStructureOrNull(String arg) {
		return readFromFileOrPdbCode(arg, false, true);		
	}
	
	/**
	 * Not sure why in this method the 'exit' parameter to readFromFileOrPdbCode was set to true
	 * @param arg
	 * @param chain
	 * @return
	 */
	public static PdbChain readStructureOrNull(String arg, String chain) {
		// Not sure why in the following call the 'exit' parameter was set to true,
		// I'm setting it to false beacuse otherwise it would cause the program to exit.
		return readFromFileOrPdbCode(arg, chain, false, true);		
	}	

	/**
	 * Loads a pdb structure given a pdbcode+chaincode or a pdb file name and chain code.
	 * Common exceptions are caught internally. The behaviour in case of an error is
	 * specified by the parameters <code>exit</code> and <code>silent</code>.
	 * Parameter chain code is only used if first parameter is a file, otherwise ignored.
	 * @param arg a pdbcode+chaincode (e.g. 1tdrB) or a pdb file name
	 * @param chain a chain code or null (first chain code in file)
	 * @param exit if true, system.exit(1) is called on error
	 * @param silent if false, error messages will be printed
	 * @return the structure object
	 */
	public static PdbChain readFromFileOrPdbCode(String arg, String chain, boolean exit, boolean silent) {
		PdbChain pdb = null;
		
		// check if argument is a filename
		File inFile = new File(arg);
		if(inFile.canRead()) {
			if(!silent) System.out.println("Reading file " + arg);
			try {
				PdbAsymUnit fullpdb = new PdbAsymUnit(new File(arg));
				if(chain == null) {
					//String[] chains = parser.getChains();
					if(!silent) System.out.println("Loading chain " + fullpdb.getFirstChain().getPdbChainCode());
					pdb = fullpdb.getFirstChain();
				} else {
					if(!silent) System.out.println("Loading chain " + chain);
					pdb = fullpdb.getChain(chain);
				}
			} catch (PdbLoadException e) {
				if(!silent) System.err.println("Error loading file " + arg + ":" + e.getMessage());
			} catch (FileFormatException e) {
				if(!silent) System.err.println("Error loading file " + arg + ":" + e.getMessage());
			} catch (IOException e) {
				if(!silent) System.err.println("Error loading file " + arg + ":" + e.getMessage());
			}
		} else {
			// check if argument is a pdb code
			if(arg.length() < 4 || arg.length() > 5) {
				if(!silent) System.err.println(arg + "is neither a valid file name nor a valid pdb(+chain) code");
				if(exit) System.exit(1);
			} else {
				String pdbCode = arg.substring(0,4);
				String chainCode = arg.substring(4,5);
				try {
					if(!silent) System.out.println("Loading pdb code " + pdbCode);
					try {
						PdbAsymUnit fullpdb = new PdbAsymUnit(pdbCode, new MySQLConnection(), PdbaseParser.DEFAULT_PDBASE_DB);
						
						if(chainCode.length() == 0) {

							chainCode = fullpdb.getFirstChain().getPdbChainCode();
						}
						
						if(!silent) System.out.println("Loading chain " + chainCode);
						pdb = fullpdb.getChain(chainCode);

					} catch (PdbLoadException e) {
						if(!silent) System.err.println("Error loading pdb structure:" + e.getMessage());
						if(exit) System.exit(1);
					}

				} catch (SQLException e) {
					if(!silent) System.err.println("Database error: " + e.getMessage());
					if(exit) System.exit(1);
				} catch (PdbCodeNotFoundException e) {
					if(!silent) System.err.println("Pdb code " + pdbCode + " not found in database.");
					if(exit) System.exit(1);
				}
			}
		}
		return pdb;	
	}
	
	
	/**
	 * Loads a pdb structure given a pdbcode+chaincode or a pdb file name.
	 * Common exceptions are caught internally. The behvaiour in case of an error is
	 * specified by the parameters <code>exit</code> and <code>silent</code>.
	 * @param arg a pdbcode+chaincode (e.g. 1tdrB) or a pdb file name
	 * @param exit if true, system.exit(1) is called on error
	 * @param silent if false, error messages will be printed
	 * @return the structure object
	 */
	public static PdbChain readFromFileOrPdbCode(String arg, boolean exit, boolean silent) {
		return readFromFileOrPdbCode(arg, null, exit, silent);
	}



}


