package owl.mutanom.core;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;

import owl.core.connections.CSAConnection;
import owl.core.features.CsaFeature;
import owl.core.features.Feature;
import owl.core.features.FeatureType;
import owl.core.features.InvalidFeatureCoordinatesException;
import owl.core.features.OverlappingFeatureException;
import owl.core.runners.NaccessRunner;
import owl.core.sequence.alignment.PairwiseSequenceAlignment;
import owl.core.sequence.alignment.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;
import owl.core.structure.PdbChain;
import owl.core.structure.features.CatalSiteSet;
import owl.core.structure.features.CatalyticSite;
import owl.core.structure.graphs.RIGraph;
import owl.mutanom.output.PyMolScriptMaker;
import owl.mutanom.output.PyMolScriptMaker.Color;



/**
 * A structural region in a protein.
 * @author stehr
 *
 */
public class Substructure {

	/*--------------------------- type definitions --------------------------*/
	public enum SubstructureType {XRAY, NMR, PREDICTION, OTHER};
	
	/*------------------------------ constants ------------------------------*/
	public static final File LOCAL_SIFTS_FILE = new File("/project/StruPPi/henning/projects/mutanom/analysis/data/sifts/pdb_chain_uniprot.lst");
	public static final String PDB_FILE_URL = "file:///project/StruPPi/BiO/DBd/PDB-REMEDIATED/data/structures/unzipped/all/pdb/pdb%s.ent";
	public static final String RENUM_SCRIPT = "/project/StruPPi/CASP8/scripts/renumpdb.py %d %s %s %s";
	public static final String PDB2UNIPROT_SQL = "SELECT sp_id from mutanom.pdb2uniprot where pdb=\"%s\" and chain=\"%s\";";
	public static final String PDB2UNIPROT_SQL_OFFSET = "SELECT sp_beg-pdb_beg as offset from mutanom.pdb2uniprot where pdb=\"%s\" and chain=\"%s\";";
	public static final int ALIGNMENT_UNDEFINED = Integer.MAX_VALUE;
	public static final String ALIGNMENT_UNDEFINED_STR = null;
	
	
	public static final String DEFAULT_RIG_TYPE = "Ca";
	public static final double DEFAULT_RIG_CUTOFF = 8.0;
	
	private static final String NACCESS_EXECUTABLE = "/project/StruPPi/bin/naccess";
	private static final String NACCESS_PARAMETERS = "";
	public static final double EXPOSURE_CUTOFF = 5.0; // everything above this cutoff is considered exposed
	/*--------------------------- member variables --------------------------*/
	int start;	// start position in Uniprot sequence
	int end;	// end position in Uniprot sequence
	int offset;	// offset between residue numbers in uniprot sequence and residue numbers in Pdbase sequence
				// this is a quick workaround for the problem that renumbered PDB files are not being processed
				// correctly by the PdbChain class. So instead, we load known structures from Pdbase and apply this
				// offset to locate the mutation in the structure. In the long run this should be replaced by
				// a more general solution using an alignment. Update: Now we use the alignment, offset should not be used anymore!
	SubstructureType type;
	File pdbFile;
	
	// additional data which can be optionally loaded to enable further analysis
	// each method which requires one of these is responsible for checking that the data has been loaded
	PdbChain pdb; // the actual coordinates (to be loaded with loadPdb...)
	RIGraph graph; // the accompanying graph (to be loaded with loadGraph())
	Map<FeatureType, Collection<Feature>> features; // to be loaded with loadFeatures()
	
	// alignment between uniprot sequence and pdb sequence (to be loaded with initializeAlignment...)
	Map<Integer,Integer> cif2uni;	// stores the mapping of residue numbers from cif to uniprot
	Map<Integer,Integer> uni2cif;	// stores the mapping of residue numbers from uniprot to cif
	Map<String,Integer> pdb2uni;	// stores the mapping of residue numbers from pdbfile to uniprot
	Map<Integer,String> uni2pdb;	// stores the mapping of residue numbers from uniprot to pdbfile	
	boolean alignmentInitialized = false;
	
	// for experimental structures
	String pdbCode;
	String chainCode;	// sometimes multiple identical chains match, we ignore this for the moment
	
	// for predicted structures
	String templatesString;		// workaround to store template information as string before structured storage is implemented
	Collection<Template> templates;
	
	/*----------------------------- constructors ----------------------------*/

//	/**
//	 * Create a new substructure object
//	 */
//	public Substructure(int start, int end, SubstructureType type, File pdbFile) {
//		this.start = start;
//		this.end = end;
//		this.type = type;
//		this.pdbFile = pdbFile;
//		this.pdbCode = null;
//		this.offset = ALIGNMENT_UNDEFINED;
//		pdb = null;
//	}
	
	/**
	 * Create a new substructure object for a structure from the PDB.
	 * @param start
	 * @param end
	 * @param type
	 * @param pdbFile
	 */
	public Substructure(int start, int end, SubstructureType type, String pdbCode, String chainCode) {
		this.start = start;
		this.end = end;
		this.type = type;
		this.pdbCode = pdbCode;
		this.chainCode = chainCode;
		this.pdbFile = null;
		this.offset = ALIGNMENT_UNDEFINED;
		pdb = null;
	}

	/**
	 * Initializes the alignment from uniprot to pdb sequence. Pdb structures have to be loaded previously.
	 * Predicted structures have to be numbered from 1 (so not renumbered according to Uniprot sequence).
	 * The full PDB sequence is used.
	 */
	public void initializeAlignment(String uniprotSeq, String uniprotId) { 
		// get alignment end points
		int uniBeg, uniEnd, cifBeg, cifEnd;
		if(this.type == SubstructureType.PREDICTION) {
			uniBeg = this.start;
			uniEnd = this.end;
			cifBeg = 1;
			//cifEnd = this.end - this.start + 1;	// wrong because there may be gaps in the Uniprot2Pdb alignment
			cifEnd = this.pdb.getFullLength();
		} else {
			// dirty workaround: this structure (1azs) is from a different Uniprot entry with almost identical sequence,
			// so we use it even it is not in SIFTS
			if(uniprotId.equals("P63092")) {
				uniBeg = 1;
				uniEnd = 394;
				cifBeg = 1;
				cifEnd = 394;			
			} else {
				try {
					SiftsMapping sm = new SiftsMapping(this.pdbCode, this.chainCode, LOCAL_SIFTS_FILE);
					uniBeg = sm.uniBeg;
					uniEnd = sm.uniEnd;
					cifBeg = sm.cifBeg;
					cifEnd = sm.cifEnd;
					if(!sm.uniprotId.equals(uniprotId)) {
						System.err.println("Warning: UniprotId " + uniprotId + " does not match the one found in SIFTS record " + sm.uniprotId);
					}
				} catch(IOException e) {
					System.err.println("IOException while retrieving SIFTS mapping: " + e.getMessage());
					return;
				} catch (SiftsMappingNotFoundException e) {
					System.err.println("Could not find " + pdbCode + " " + chainCode + " in SIFTS file");
					return;
				}
			}
		}		
		// align uniprot to (full) pdb sequence
		Map<Integer,Integer> uni2cif = new HashMap<Integer, Integer>();
		Map<Integer,Integer> cif2uni = new HashMap<Integer, Integer>();
		Map<Integer,String> uni2pdb = new HashMap<Integer, String>();
		Map<String,Integer> pdb2uni = new HashMap<String, Integer>();
		String pdbSeq = this.pdb.getSequence().getSeq().substring(cifBeg-1, cifEnd);
		String uniSeq = uniprotSeq.substring(uniBeg-1, uniEnd);
		try {
			PairwiseSequenceAlignment al = new PairwiseSequenceAlignment(uniSeq, pdbSeq,"Uniprot","Pdb");
			if(al.getPercentIdentity() < 95) {
				System.err.println("Warning: Alignment between " + uniprotId + " and " + pdbCode+chainCode + " has sequence identity < 95%");
				// DEBUG: al.printAlignment();
			}
			for (int i = 0; i < al.getMapping1To2().length; i++) {
				int uni = i + uniBeg;
				int cif = al.getMapping1To2()[i] + cifBeg;
				// TODO: check for unobserved
				if(cif > 0) {
					uni2cif.put(uni, cif);
					String pdb = this.pdb.getPdbResSerFromResSer(cif);
					if(pdb != null) uni2pdb.put(uni, pdb);
				}
			}
			for (int i = 0; i < al.getMapping2To1().length; i++) {
				int cif = i + cifBeg;
				int uni = al.getMapping2To1()[i] + uniBeg;
				// TODO: check for unobserved
				if(uni > 0) {
					cif2uni.put(cif, uni);
					String pdb = this.pdb.getPdbResSerFromResSer(cif);
					if(pdb != null) pdb2uni.put(pdb, uni);
				}
			}			
			this.uni2cif = uni2cif;
			this.cif2uni = cif2uni;
			this.uni2pdb = uni2pdb;
			this.pdb2uni = pdb2uni;
			alignmentInitialized = true;
		} catch (PairwiseSequenceAlignmentException e) {
			System.err.println("Could not align sequences: " + e.getMessage());
			return;
		}
	}
	
	/*-------------------------- getters and setters ------------------------*/
	
//	/**
//	 * Returns the offset between residue numbers in the Uniprot sequence and the Pdbase sequence.
//	 * Returns OFFSET_UNDEFINED if offset hasn't been set.
//	 */
//	public int getOffset() {
//		return this.offset;
//	}

	/**
	 * Sets the offset between residue numbers in the Uniprot sequence and the Pdbase sequence
	 */
	public void setOffset(int offset) {
		this.offset = offset;
	}
	
	public void setTemplatesString(String s) {
		this.templatesString = s;
	}
	
	public String getTemplatesStr() {
		return this.templatesString;
	}
	
	/**
	 * Gets the starting position of this substructure in the reference sequence.
	 * @return the starting position (in amino acids)
	 */
	public int getBegPos() {
		return this.start;
	}
	
	/**
	 * Gets the end position of this substructure in the reference sequence.
	 * @return the end position (in amino acids)
	 */	
	public int getEndPos() {
		return this.end;
	}
	
	public String getRange() {
		return this.start + "-" + this.end;
	}
	
	public SubstructureType getType() {
		return this.type;
	}
	
	public char getTypeChar() {
		char c;
		if(this.type != null) {
			switch(this.type) {
			case XRAY: c = 'X'; break;
			case NMR: c = 'N'; break;
			case PREDICTION: c = 'P'; break;
			default: c = 'O';
			}
		} else c = 'U';	// unknown
		return c;
	}
	
	/**
	 * Returns the four letter pdb code of this substructure or null if none is specified.
	 * @return the four letter pdb code or null if none is specified
	 */
	public String getPdbCode() {
		if(pdbCode != null && pdbCode.length() >= 4) {
			return pdbCode.substring(0,4);
		} else {
			return null;
		}
	}
	
	/**
	 * Returns the one letter chain code of this substructure or null if none is specified.
	 * @return the one letter pdb chain code or null if none is specified
	 */
	public String getChainCode() {
		if(this.chainCode != null && this.chainCode.length() > 0) {
			return this.chainCode.substring(0,1);
		} else {
			return null;
		}
	}
	
	/**
	 * Returns the pdb object associated with this substructure.
	 * The Pdb object can be loaded with loadPdb...
	 * @return the pdb object or null if it has not yet been loaded
	 */
	public PdbChain getPdb() {
		return this.pdb;
	}
	
	/**
	 * Returns the RIGraph object associated with this substructure.
	 * The PdbChain object can be loaded with loadGraph().
	 * @return the graph object or null if it has not yet been loaded
	 */	
	public RIGraph getGraph() {
		return this.graph;
	}
	
	/*-------------------------- implemented methods ------------------------*/
	
	public boolean addFeature(Feature feature) throws InvalidFeatureCoordinatesException, OverlappingFeatureException {
		boolean result = false;
		FeatureType ft = feature.getType();	// will throw a null pointer exception if feature == null
		if(this.features == null) {
			this.features = new HashMap<FeatureType, Collection<Feature>>();
		}
		Collection<Feature> fc = this.features.get(ft);
		if(fc==null) {
			fc = new LinkedList<Feature>();
			features.put(ft, fc);
		}
		// TODO: check for invalid coordinates and overlapping features
		result = fc.add(feature);
		return result;
	}
	
	public Collection<Feature> getFeatures() {
		Collection<Feature> result = new LinkedList<Feature>();
		for(FeatureType type:this.features.keySet()) {
			for(Feature f:this.features.get(type)) {
				result.add(f);
			}
		}
		return result;
	}
	
	public Collection<Feature> getFeaturesOfTypeForPosition(FeatureType featureType, int position) {
		Collection<Feature> result = new LinkedList<Feature>(); 
			for(Feature f:this.features.get(featureType)) {
				if(f.getIntervalSet().getIntegerSet().contains(position)) result.add(f);
			}
		return result;
	}
	
	/*---------------------------- public methods ---------------------------*/
	public String toString() {
		return String.format("%4d-%4d (%s) [%s]", start, end, getTypeChar(), this.type == SubstructureType.PREDICTION?templatesString:getPdbCode()+getChainCode());
	}
	
	/**
	 * @return true iff the coordinates for this substructure have been loaded.
	 */
	public boolean isPdbLoaded() {
		return pdb != null;
	}
	
	/**
	 * @return true iff the alignment has been initialized.
	 */
	public boolean isAlignmentInitialized() {
		return alignmentInitialized;
	}
	
	/**
	 * @return true iff a RIG for this substructure has been loaded (requires PdbChain to be loaded first)
	 */
	public boolean isGraphLoaded() {
		return this.graph != null;
	}
	
	/**
	 * Loads the coordinates from the given pdb file into memory and stores the
	 * filename in the pdbFile member.
	 */
	public void loadPdbFromFile(File pdbFile) {
		// DEBUG: System.out.println("Loading file " + pdbFile);
		this.pdb = PdbChain.readStructureOrNull(pdbFile.toString(), this.chainCode);
		if(this.pdb != null) this.pdbFile = pdbFile; else System.err.println("Error loading PDB file " + pdbFile);
	}
	
	/**
	 * Loads the coordinates from the pdb file into memory using the pdbCode member.
	 */
	public void loadPdbFromDb() {
		if(this.pdbCode != null && this.chainCode != null) {
			this.pdb = PdbChain.readStructureOrNull(this.pdbCode+this.chainCode);
			// DEBUG:
//			if(this.pdbCode.equals("2ovr")) {
//				System.out.println("Loading " + this.pdbCode + " " + this.chainCode);
//				System.out.println(this.pdb.getSequence());
//			}
			// DEBUG:
//			if(this.pdbCode.equals("2o8b")) {
//				System.out.println("Loading " + this.pdbCode + " " + this.chainCode);
//				System.out.println("Cif chain code: " + this.pdb.getChainCode());
//				System.out.println("Pdb chain code: " + this.pdb.getPdbChainCode());
//			}

			
		}
	}
	
	/**
	 * Loads the residue interaction graph for this substructure.
	 * Requires PdbChain to be loaded first.
	 */
	public void loadGraph(String contactType, double cutoff) {
		if(!isPdbLoaded()) {
			System.err.println("Could not load RIG because structure is not available");
			return;
		}
		this.graph = this.pdb.getRIGraph(contactType, cutoff);
	}
	
	/**
	 * Loads structural/functional features for this substructure from external databases.
	 * Features include (known or predicted): 
	 * - ATP binding sites
	 * - phosphorilation sites
	 * - catalytic sites
	 * - metal binding sites
	 * - DNA binding site?
	 * @return the number of features loaded
	 */
	public int loadFeatures() {
		if(this.pdb == null) {
			System.err.println("Could not load features: structure not loaded");
			return 0;
		}
		int featuresLoaded = 0;
		if(this.features == null) this.features = new HashMap<FeatureType, Collection<Feature>>();
		// load catalytic sites from CSA
		try {
			CSAConnection.parseCSA(this.pdb,CatalSiteSet.LATEST_VERSION, false);
			if(this.pdb.getCSA() != null) {
				Iterator<CatalyticSite> it = this.pdb.getCSA().getIterator();
				while(it.hasNext()) {
					CatalyticSite cs = it.next();
					try {
						if(this.addFeature(new CsaFeature(cs))) featuresLoaded++;
					} catch (InvalidFeatureCoordinatesException e) {
					} catch (OverlappingFeatureException e) {
					}
				}
			}
		} catch (IOException e) {
			System.err.println("Error reading from Catalytic Site Atlas: " + e.getMessage());
		}
		// load features from Uniprot
		
		
		// load features from ProSite

		//System.out.println("Features found for this substructure:");
		for(Feature f:this.getFeatures()) {
			System.out.println(f.toString());
		}
		return featuresLoaded;
	}
	
	/**
	 * If substructure type is experimental (Xray or NMR), download the PDB file to the given directory.
	 * If type is different or pdbCode is not specified, does nothing.
	 * @return a handle for the downloaded file or null if nothing was downloaded
	 */
	public File downloadStructureFromPdb(File newFile) {
		String pdbCode = this.pdbCode;
		if((this.type == SubstructureType.XRAY || this.type == SubstructureType.NMR) && pdbCode != null) {
			try {
				System.out.println("Writing file " + newFile);
				PrintWriter out = new PrintWriter(newFile);
				URL url = new URL(String.format(PDB_FILE_URL, pdbCode));
				BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
				String line;
				while((line = in.readLine()) != null) {
					out.println(line);
				}
				out.close();
				in.close();
				return newFile;
			} catch (MalformedURLException e) {
				System.err.println("Error downloading PDB file for " + pdbCode + ": " + e.getMessage());
			} catch (IOException e) {
				System.err.println("Error downloading PDB file for " + pdbCode + ": " + e.getMessage());
			}
		}
		return null;
	}
	
	/**
	 * Returns the relative surface accessibility for the given position or NaN if something went wrong.
	 * Structure needs to be loaded previously (with correct residue numbering). If rsa can not be
	 * retrieved from PDB object initially, tries to run NACCESS.
	 * @param pos the position for which to evaluate the surface accessibility.
	 * @return the relative surface accessibility of the given position or Double.NaN
	 */
	public double getSurfaceAccessibility(int pos) {
		Double rsa = Double.NaN;
		if(isPdbLoaded()) {
			// assuming that numbering of residues is correct
			if(pdb.getSurfaceAccessibilities() == null) {
				// If results are not cached, run NACCESS
				try {
					System.out.println("Running NACCESS on " + this.getPdbCode() + " " + this.getChainCode());
					NaccessRunner nar = new NaccessRunner(new File(NACCESS_EXECUTABLE), NACCESS_PARAMETERS);
					nar.runNaccess(pdb);
					//pdb.runNaccess();
				} catch (IOException e) {
					System.err.println("Error running NACCESS: " + e.getMessage());
					//System.exit(1);
				}
			}
			rsa = pdb.getAllRsaFromResSerial(pos);
			if(rsa==null) System.out.println("No surface accessibility value for position " + pos);
		} else {
			System.err.println("Structure for target " + this.getPdbCode() + " is not loaded");
		}
		return rsa==null?Double.NaN:rsa.doubleValue();
	}
	
	/**
	 * Returns true iff the residue at the given position has a relative surface accessibility above
	 * the threshold EXPOSURE_CUTOFF
	 * @param pos
	 * @return true if the residue is exposed, false otherwise
	 */
	public boolean isExposed(int pos) {
		double rsa = getSurfaceAccessibility(pos);
		return rsa > EXPOSURE_CUTOFF;
	}
	
	/**
	 * Returns the default name for a pdb file associated with this substructure.
	 * @return the default pdb file name
	 */
	public String getDefaultPdbFileName(String geneName) {
		if(this.type == SubstructureType.PREDICTION) {
			//return geneName + "_" + this.getRange() + "_pred.renum.pdb";
			return geneName + "_" + this.getRange() + ".pdb";
		} else return null;
// obsolete:
//		else {
//			String s = geneName + "_" + this.getRange() + "_" + this.getPdbCode() + ".renum.pdb";
//			if(this.getPdbCode() == null || this.getPdbCode().length() == 0) System.err.println("Warning: Null PDB code in file name " + s);
//			return s;
//		}
	}
	
	/*---------------------------- static methods ---------------------------*/
	public static SubstructureType parseType(String s) {
		SubstructureType sst = null;
		char c = s.charAt(0);
		switch(c) {
			case 'X': sst = SubstructureType.XRAY; break;
			case 'N': sst = SubstructureType.NMR; break;
			case 'P': sst = SubstructureType.PREDICTION; break;
			default: sst = SubstructureType.OTHER; break;
		}
		return sst;
	}
	
	/**
	 * Renumber a pdb text file by adding the given offset to each residue number.
	 * @param pdbFile
	 * @param offset
	 */
	public static void renumberPdb(File inFile, File outFile, int offset, char chain) {
		String cmd = String.format(RENUM_SCRIPT, offset+1, inFile, outFile, chain);
		try {
			Runtime.getRuntime().exec(cmd);
		} catch (IOException e) {
			System.err.println("Error executing command line: " + cmd);
		}
	}

// DELETE:
//	/**
//	 * Queries the official PDB to UniProt mapping and returns the offset by which
//	 * the pdb has to be shiftet in order to match the Uniprot residue numbering.
//	 * Opens a database connection using default parameters ~/.my.cnf
//	 * @param pdbFile
//	 * @param chain
//	 * @return
//	 * See also: Gene.getOffsetFromSifts()
//	 */
//	public static int getPdb2UniprotOffset(String pdbCode, String chainCode) {
//		int offset = 0;
//		String query = String.format(PDB2UNIPROT_SQL_OFFSET, pdbCode.toLowerCase(), chainCode.toUpperCase());
//		try {
//			MySQLConnection conn = new MySQLConnection();
//			offset = conn.getIntFromDb(query);
//			conn.close();
//		} catch (SQLException e) {
//			System.err.println("Error executing database query " + query + ": " + e.getMessage());			
//		}
//		return offset;
//	}
//	
//	/**
//	 * Queries the official PDB to UniProt mapping and returns the UniprotID
//	 * maching the given pdb and chain code.
//	 * Opens a database connection using default parameters in ~/.my.cnf
//	 * @param pdbFile
//	 * @param chain
//	 * @return the uniprot id for the given pdb and chain code or null if none was found.
//	 */
//	public static String getPdb2Uniprot(String pdbCode, String chainCode) {
//		String uniprotID = null;
//		String query = String.format(PDB2UNIPROT_SQL, pdbCode.toLowerCase(), chainCode.toUpperCase());
//		try {
//			MySQLConnection conn = new MySQLConnection();
//			uniprotID = conn.getStringFromDb(query);
//			conn.close();
//		} catch (SQLException e) {
//			System.err.println("Error executing database query " + query + ": " + e.getMessage());
//		}
//		return uniprotID;
//	}
	
	/**
	 * Returns the cif residue number in the loaded pdb structure corresponding to the given position in the Uniprot sequence.
	 * @return the residue number in the pdb structure or ALIGNMENT_UNDEFINED if no such residue exists
	 */
	public int mapUniprotResser2Cif(int uniprotResser) {
		if(!this.alignmentInitialized || this.uni2cif == null) {
			System.err.println("Error: Alignment not initialized for " + this.toString());
			return ALIGNMENT_UNDEFINED;
		}
		Integer cifResser = uni2cif.get(uniprotResser);
		if(cifResser == null) {
			return ALIGNMENT_UNDEFINED;
		} else {
			return cifResser.intValue();
		}
	}
	
	/**
	 * Returns the residue number in the uniprot sequence corresponding to the given cif residue number in the loaded pdb structure.
	 * @return the residue number in the uniprot sequence or ALIGNMENT_UNDEFINED if no such residue exists
	 */
	public int mapCifResser2Uniprot(int cifResser) {
		if(!this.alignmentInitialized || this.cif2uni == null) {
			System.err.println("Error: Alignment not initialized for " + this.toString());
			return ALIGNMENT_UNDEFINED;
		}
		Integer uniprotResser = cif2uni.get(cifResser);
		if(uniprotResser == null) {
			return ALIGNMENT_UNDEFINED;
		} else {
			return uniprotResser.intValue();
		}
	}
	
	/**
	 * Returns the pdbfile residue number in the loaded pdb structure corresponding to the given position in the Uniprot sequence.
	 * @return the residue number (possibly with insertion code) in the pdb structure or ALIGNMENT_UNDEFINED_STR if no such residue exists
	 */
	public String mapUniprotResser2Pdb(int uniprotResser) {
		if(!this.alignmentInitialized || this.uni2pdb == null) {
			System.err.println("Error: Alignment not initialized for " + this.toString());
			return ALIGNMENT_UNDEFINED_STR;
		}
		String pdbResser = uni2pdb.get(uniprotResser);
		if(pdbResser == null) {
			return ALIGNMENT_UNDEFINED_STR;
		} else {
			return pdbResser;
		}
	}
	
	/**
	 * Returns the residue number in the uniprot sequence corresponding to the given pdbfile residue number in the loaded pdb sequence.
	 * @return the residue number in the uniprot sequence or ALIGNMENT_UNDEFINED if no such residue exists
	 */
	public int mapPdbResser2Uniprot(String pdbResser) {
		if(!this.alignmentInitialized || this.pdb2uni == null) {
			System.err.println("Error: Alignment not initialized for " + this.toString());
			return ALIGNMENT_UNDEFINED;
		}
		Integer uniprotResser = pdb2uni.get(pdbResser);
		if(uniprotResser == null) {
			return ALIGNMENT_UNDEFINED;
		} else {
			return uniprotResser.intValue();
		}
	}
	
	public void writePseWithAllMutations(File outDir, String baseName, File pdbFileDir, Collection<Mutation> mutations) throws IOException {
		File sessionFile = new File(outDir, baseName + "_" + this.getRange() + ".pse");
		String geneName = baseName;
		String pdbFileName = this.getDefaultPdbFileName(geneName);		
		String objName = geneName + "_" + this.getRange();
		PyMolScriptMaker psm = new PyMolScriptMaker(true);
		psm.loadDefaultSettings();
		psm.load(new File(pdbFileDir, pdbFileName), objName);
		for(Mutation m:mutations) {
			String selName = String.format("%s%d", m.before.getOneLetterCode(), m.position);			
			psm.highlightResidue(objName, m.position, this.getChainCode().charAt(0), Color.RED, selName, null);
		}
		
		psm.saveSession(sessionFile);
		psm.executeAndClose();
	}
	
}
