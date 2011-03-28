package owl.mutanom.core;

import java.io.*;
import java.sql.*;
import java.net.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.connections.CSAConnection;
import owl.core.connections.NoMatchFoundException;
import owl.core.connections.PDomainsConnection;
import owl.core.connections.PhosphoSiteConnection;
import owl.core.connections.PrositeHit;
import owl.core.connections.PrositeScanner;
import owl.core.connections.UniProtConnection;
import owl.core.features.CsaFeature;
import owl.core.features.Feature;
import owl.core.features.FeatureType;
import owl.core.features.GeneralFeature;
import owl.core.features.HasFeatures;
import owl.core.features.InvalidFeatureCoordinatesException;
import owl.core.features.OverlappingFeatureException;
import owl.core.features.PhosphoSitePlusFeature;
import owl.core.features.PrositeFeature;
import owl.core.features.ProteinModificationType;
import owl.core.features.StructuralDomainFeature;
import owl.core.features.StructuralDomainType;
import owl.core.features.UniprotFeature;
import owl.core.sequence.Sequence;
import owl.core.sequence.alignment.PairwiseSequenceAlignment;
import owl.core.sequence.alignment.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;
import owl.core.structure.AminoAcid;
import owl.core.structure.features.CatalSiteSet;
import owl.core.structure.features.CatalyticSite;
import owl.core.util.FileFormatException;
import owl.core.util.Interval;
import owl.core.util.IntervalSet;
import owl.core.util.MySQLConnection;
import owl.core.util.Pair;
import owl.mutanom.core.Substructure.SubstructureType;


import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseType;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.dbx.pdb.Pdb;
import uk.ac.ebi.kraken.interfaces.uniprot.features.HasFeatureDescription;

/**
 * Gathers information about a given gene, such as:
 * - UniprotID
 * - Genetic/Protein Sequence
 * - Known and predicted structural regions
 * - Structural Features
 * - Mutations
 * - Information about the given mutations, such as:
 *   - clustering index
 *   - surface accessibility
 *   - proximity to catalytic and ATP binding sites
 * and provides methods to output and visualize this information as:
 * - text
 * - pymol script/3D image
 * - sequence scheme (perl sequence rendering)
 * @author stehr
 */
public class Gene implements HasFeatures {

	/*--------------------------- type definitions --------------------------*/
	public enum SeqType {DNA,PROTEIN};
	public enum SubcellularLocalization{NUCLEUS, CYTOPLASM, MEMBRANE, TRANSMEMBRANE};
	
	/*------------------------------ constants ------------------------------*/
	// old
	//public static final File TARGET_FILE = new File("/project/PyMol/targets/targets_top20_20090320.mod.txt");
	//public static final File TARGET_FILE = new File("/project/PyMol/targets/targets_top29_20091106.mod.txt");
//	public static final File TARGET_FILE = new File("/project/PyMol/targets/targets_top29_20100302.mod.txt");
	//public static final File SUBSTRUCTURE_FILE = new File("/project/PyMol/targets/substructures_top20_20090320.mod.txt");
	//public static final File SUBSTRUCTURE_FILE = new File("/project/PyMol/targets/substructures_top29_20091106.mod.txt");
//	public static final File TARGET_FILE_11 = new File("/project/PyMol/targets/targets_top11_20100309.mod.txt");
//	public static final File SUBSTRUCTURE_FILE = new File("/project/PyMol/targets/substructures_top29_2010-02-25.mod.txt");	
//	public static final File MUTATION_FILE = new File("/project/PyMol/targets/mutations_top20_20090504.txt");
//	public static final File ANNOTATION_FILE = new File("/project/PyMol/targets/annotations_top20_20091007.txt");

//	public static final String MUTATION_TABLE = "mutanom.missense_mutation";
//	public static final String MUTATION_QUERY = "SELECT DISTINCT mut_aa FROM " + MUTATION_TABLE + " WHERE gene=\"%s\"";
//	public static final String MUTATION_TABLE2 = "mutanom.cosmic43_mis_4tis";
//	public static final String MUTATION_TABLE3 = "mutanom.cosmic43_mis_5tis";	
//	public static final String MUTATION_QUERY2 = "SELECT DISTINCT mut_aa FROM " + MUTATION_TABLE2 + " WHERE gene_name=\"%s\"";	
//	public static final String MUTATION_QUERY3 = "SELECT DISTINCT mut_aa FROM " + MUTATION_TABLE3 + " WHERE gene_name=\"%s\"";	

//	public static final String ENSP_QUERY = "SELECT ensp FROM mutanom.sp_with_snps WHERE sp_id = \"%s\" LIMIT 1";
//	public static final String SNP_QUERY = "SELECT pos, aa_ref, aa_mut FROM mutanom.snp WHERE ensp = \"%s\" AND aa_ref != aa_mut";	
//	
//	public static final String ENSG_QUERY = "SELECT ensg FROM uniprot_15_8.sp2ens_fromUniprot where sp_id = \"%s\" LIMIT 1";
//	public static final String ENST_QUERY = "SELECT enst FROM mutanom.ens47_sp2ens where sp_id = \"%s\"";	// adding LIMIT 1 will results in the same set of SNPs	
//	public static final String SNP_QUERY2 = "SELECT beg_aa, mut FROM mutanom.ensvar56_top21 WHERE ensg=\"%s\" AND mut RLIKE \"^[A-Z]/[A-Z]$\" AND beg_aa = end_aa";
//	public static final String SNP_QUERY3 = "SELECT beg_aa, mut FROM mutanom.ensvar56_top21 WHERE enst in (%s) AND mut RLIKE \"^[A-Z]/[A-Z]$\" AND beg_aa = end_aa group by mut,beg_aa";
	
	// This should not be used anymore, but we keep it here in case we want to compare
	//public static final File ANNOTATION_FILE = new File("/project/PyMol/targets/annotations_top20_20091007.txt");
	public static final File ANNOTATION_FILE = new File("/project/StruPPi/henning/projects/mutanom/analysis/data/manual/pik3ca.annotation.txt");
	
	// new

	// configurables, TODO: load from config file or environment variable (through static method)
	public static final String BASE_DIR =				"/project/StruPPi/henning/projects/mutanom/analysis/";
	public static final String DATABASE = 				"mutanom3"; 
	
	// database
	public static final String SUBSTRUCTURE_TABLE = DATABASE + ".substructures";
	public static final String SUBSTRUCTURE_QUERY = "SELECT pos_beg,pos_end,pdb_code,chain_code,type FROM " + SUBSTRUCTURE_TABLE + " WHERE gene_name=\"%s\"";
	
	public static final String MUTATION_TABLE = DATABASE + ".mutations";
	public static final String MUTATION_QUERY = "SELECT DISTINCT CONCAT(mut_dna, ';', mut_aa) FROM " + MUTATION_TABLE + " WHERE gene_name=\"%s\"";
	public static final String MUT_ID_QUERY = "SELECT mut_id FROM " + MUTATION_TABLE + " WHERE gene_name=\"%s\" AND mut_dna=\"%s\" AND mut_aa=\"%s\"";
	
	public static final String SNP_TABLE = 	DATABASE + ".snps";
	public static final String SNP_QUERY = 	"SELECT dbsnp_id, SUBSTR(mut_aa,1,1) AS ref_aa, SUBSTR(mut_aa,3,1) AS mut_aa, pos_aa " +
											"FROM " + SNP_TABLE + " WHERE mut_aa RLIKE \"^[A-Z]/[A-Z]$\"" +
																	" AND SUBSTR(mut_aa,1,1) != SUBSTR(mut_aa,3,1)" +
																	" AND clin=0" + 
																	" AND gene_name=\"%s\";";	

	public static final String ONC_SUP_TABLE = DATABASE + ".genes_onc_sup";
	public static final String ONC_SUP_QUERY = "SELECT class FROM " + ONC_SUP_TABLE + " WHERE gene_name=\"%s\"";
	
	
	// paths (normally, these should be File objects but in this case they are templates which still beed to be formatted)
	public static final String COSMIC_SEQ_PATH   = BASE_DIR + "data/cosmic/seq/%s.aa.fasta";
	public static final String COSMIC_CDNA_PATH  = BASE_DIR + "data/cosmic/seq/%s.cdna.fasta";
	public static final String UNIPROT_SEQ_PATH  = BASE_DIR + "data/uniprot/seq/%s.fasta";	
	
	// URLs
	public static final String COSMIC_SEQ_URL    = "ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/fasta_files/%s/%s_protein.txt";
	public static final String COSMIC_CDNA_URL   = "ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/fasta_files/%s/%s_cdna.txt";
	public static final String UNIPROT_SEQ_URL   = "http://www.uniprot.org/uniprot/%s.fasta";
	public static final String COSMIC_GENE_URL   = "http://www.sanger.ac.uk/perl/genetics/CGP/cosmic?action=gene&ln=%s";
	
	// misc (TODO: move the scripts to BASE_DIR/scripts)
	public static final String RENDER_CMD = "/bin/sh /project/StruPPi/henning/projects/mutanom/scripts/render_targets.sh %s %d %s";
	public static final String RULER_CMD = "/bin/sh /project/StruPPi/henning/projects/mutanom/scripts/render_ruler.sh %d %d %s";
	public static final String COSMIC_UNIPROT_REGEXP = "http://www.expasy.org/uniprot/([A-Z0-9]+)";
		
	// parameters
	public static final boolean LOAD_UNIQUE = false; // if true, discard mutations with identical positions (but different aa exchange)
	public static final boolean USE_ONLINE_CSA = false; // whether to load catalytic sites from online CSA (otherwise local file)
	
	/*--------------------------- member variables --------------------------*/
	protected String geneName;
	protected String uniProtID;
	protected String cosmicSeq;							// protein sequence from COSMIC
	protected String cosmicCdnaSeq; 					// cDNA sequence from COSMIC
	protected String uniprotSeq;						// protein sequence from UniProt
	protected PairwiseSequenceAlignment cosmic2Uniprot;
	protected int length;								// reference length (from Cosmic)
	protected Collection<Substructure> subStructures;	// loaded substructures (see loadSubstructures())
	protected Collection<Mutation> mutations;			// loaded mutations (see loadMutations())
	protected Collection<SNP> snps;						// loaded SNPs (see loadSNPs())
	
	Map<FeatureType, Collection<Feature>> features; 	// functional residues (see loadFeatures())
	
	/*----------------------------- constructors ----------------------------*/
	public Gene(String geneName) {
		this.geneName = geneName;
		this.uniProtID = null;
	}
	
	public Gene(String geneName, String uniProtID) {
		this.geneName = geneName;
		this.uniProtID = uniProtID;
	}
	
	/*-------------------------- getters and setters ------------------------*/
	public String getGeneName() {
		return this.geneName;
	}
	
	public void setUniprotId(String uniprotId) {
		this.uniProtID = uniprotId;
	}
	
	public String getUniprotId() {
		return this.uniProtID;
	}
	
	public void setLength(int length) {
		this.length = length;
	}
	
	public int getLength() {
		return this.length;
	}
	
	public Collection<Substructure> getSubstructures() {
		return this.subStructures;
	}
	
	public Collection<Mutation> getMutations() {
		return this.mutations;
	}
	
	public Collection<SNP> getSNPs() {
		return this.snps;
	}
	
	/**
	 * Returns the uniprot sequence associated with this gene. Uniprot sequence can
	 * be initialized using loadUniprotSequence()
	 * @return the uniprot sequence or null if sequence has not been loaded
	 */
	public String getUniprotSeq() {
		return this.uniprotSeq;
	}
	
	public String getCosmicSeq() {
		return this.cosmicSeq;
	}
	
	public String getCosmicCdnaSeq() {
		return this.cosmicCdnaSeq;
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
		if(this.features != null) {
			for(FeatureType type:this.features.keySet()) {
				for(Feature f:this.features.get(type)) {
					result.add(f);
				}
			}
		}
		return result;
	}
	
	public Collection<Feature> getFeaturesOfTypeForPosition(FeatureType featureType, int position) {
		Collection<Feature> result = new LinkedList<Feature>();
		if(this.features != null) {
			if(this.features.get(featureType) != null) {
				for(Feature f:this.features.get(featureType)) {
					if(f.getIntervalSet().getIntegerSet().contains(position)) result.add(f);
				}
			}
		}
		return result;
	}	
	
	public Collection<FeatureType> getFeatureTypes() {
		return this.features.keySet();
	}

	public Collection<Feature> getFeaturesForPositon(int position) {
		Collection<Feature> result = new LinkedList<Feature>();
		if(this.features != null) {
			for(FeatureType type:this.features.keySet()) {
				for(Feature f:this.features.get(type)) {
					if(f.getIntervalSet().getIntegerSet().contains(position)) result.add(f);
				}
			}
		}
		return result;
	}

	public Collection<Feature> getFeaturesOfType(FeatureType featureType) {
		Collection<Feature> result = this.features.get(featureType);
		return result!=null?result:new LinkedList<Feature>();
	}

	/**
	 * Searches for features at the given position and if a feature was found, returns the functional class
	 * such as Phosphorilation, ATP binding site, etc. independent of the source (which corresponds to FeatureType).
	 * The set alreadyObserved, which has to be previously initialized, stores features which have already
	 * been reported, such that no feature will be reported twice.
	 * @param gene the gene in which to search for features
	 * @param pos the position at which to search for features
	 * @param alreadyObserved a set which stores all features which have already been found
	 * @return the identified functional class of the feature
	 */
	public Pair<Feature, String> getFirstFeatureAndFunctionalClass(int pos, HashSet<Feature> alreadyObserved) {
			// check whether this position matches a feature
			for(Feature f:this.getFeaturesOfTypeForPosition(FeatureType.UNIPROT, pos)) {
				if(!alreadyObserved.contains(f)) {	// avoid double counting of the same feature
					alreadyObserved.add(f);
					UniprotFeature uf = (UniprotFeature) f;
					if(uf.getUniprotTypeName().equals("ACT_SITE")) return new Pair<Feature, String>(f,"ACT_SITE");
					if(uf.getUniprotTypeName().equals("NP_BIND")) {
						if(uf.getDescription().startsWith("ATP")) return new Pair<Feature, String>(f,"ATP_BIND");
						else if(uf.getDescription().startsWith("GTP")) return new Pair<Feature, String>(f,"GTP_BIND");
						else System.err.println("Unknown NP_BIND feature: " + uf.getDescription()); break;
					}
					if(uf.getUniprotTypeName().equals("MOD_RES")) {
						if(uf.getDescription().startsWith("Phospho")) return new Pair<Feature, String>(f,"PHO_RES");
						else return new Pair<Feature, String>(f,"MOD_RES");
					}
					if(uf.getUniprotTypeName().equals("CARBOHYD")) return new Pair<Feature, String>(f,"GLYC_RES");
					if(uf.getUniprotTypeName().equals("DNA_BIND")) return new Pair<Feature, String>(f,"DNA_BIND");
				}
			}
		
			for(Feature f:this.getFeaturesOfTypeForPosition(FeatureType.CSA, pos)) {
				if(!alreadyObserved.contains(f)) {	// avoid double counting of the same feature
					alreadyObserved.add(f);
					return new Pair<Feature, String>(f,"ACT_SITE"); 
				}
			}
		
			for(Feature f:this.getFeaturesOfTypeForPosition(FeatureType.GENERAL, pos)) {
				if(!alreadyObserved.contains(f)) {	// avoid double counting of the same feature
					alreadyObserved.add(f);
					GeneralFeature gf = (GeneralFeature) f;
					if(gf.getDescription().equals("ACT_SITE")) return new Pair<Feature, String>(f,"ACT_SITE");
					if(gf.getDescription().equals("ATP_BIND")) return new Pair<Feature, String>(f,"ATP_BIND");
					if(gf.getDescription().equals("GTP_BIND")) return new Pair<Feature, String>(f,"GTP_BIND"); 
					if(gf.getDescription().equals("ATY_RES")) return new Pair<Feature, String>(f,"MOD_RES"); 
					if(gf.getDescription().equals("MOD_RES")) return new Pair<Feature, String>(f,"MOD_RES"); 
					if(gf.getDescription().equals("UBQ_RES")) return new Pair<Feature, String>(f,"UBQ_RES"); 
					if(gf.getDescription().equals("PHO_RES")) return new Pair<Feature, String>(f,"PHO_RES"); 
					if(gf.getDescription().equals("CARBOHYD")) return new Pair<Feature, String>(f,"GLYC_RES"); 
					if(gf.getDescription().equals("DNA_BIND")) return new Pair<Feature, String>(f,"DNA_BIND"); 
				}
			}
		
			for(Feature f:this.getFeaturesOfTypeForPosition(FeatureType.PHOSPHOSITE, pos)) {
				if(!alreadyObserved.contains(f)) {	// avoid double counting of the same feature
					alreadyObserved.add(f);
					if(((PhosphoSitePlusFeature)f).getModType() == ProteinModificationType.PHOSPHORYLATION) return new Pair<Feature, String>(f,"PHO_RES");
					else if(((PhosphoSitePlusFeature)f).getModType() == ProteinModificationType.UBIQUITINATION) return new Pair<Feature, String>(f,"UBQ_RES");
					else return new Pair<Feature, String>(f,"MOD_RES");
				}
			}
			return new Pair<Feature, String>(null,"NONE");						
	}
	
	/*---------------------------- public methods ---------------------------*/
	public String toString() {
		String s;
		s = geneName + "\t" + this.uniProtID + "\t" + length;
		if(this.subStructures != null) {
			s += "\t[";
			for(Substructure ss: this.subStructures) {
				s += ss.toString() + " ";
			}
			s += "]";
		}
		if(this.mutations != null) {
			s += "\t[";
			for(Mutation m:this.mutations) {
				s += m.toString() + " ";
			}
			s += "]";
		}
		if(this.snps != null) {
			s += "\t{";
			for(Mutation m:this.snps) {
				s += m.toString() + " ";
			}
			s += "}";
		}
		return s;
	}
	
	/**
	 * Returns whether this gene is an Oncogene.
	 * Relies on proper annotation in the database.
	 * @param conn
	 * @return true, if the gene is annotated as an oncogene, false otherwise
	 */
	public boolean isOncogene(MySQLConnection conn) {
		String type = conn.getStringFromDb(String.format(ONC_SUP_QUERY,this.getGeneName()));
		if(type != null && type.equals("O")) {
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * Returns whether this gene is a Tumor Suppressor gene.
	 * Relies on proper annotation in the database.
	 * @param conn
	 * @return true, if the gene is annotated as a tumor suppressor, false otherwise
	 */	
	public boolean isTumorSuppressor(MySQLConnection conn) {
		String type = conn.getStringFromDb(String.format(ONC_SUP_QUERY,this.getGeneName()));
		if(type != null && type.equals("S")) {
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * Returns the annotation symbol for this gene in the Oncogene/Tumor Suppressor table.
	 * In most cases, the methods isOncogene() and isTumorSuppressor() should be used instead.
	 * @param conn
	 * @return the annotation symbol for this gene in the database or null if no annotation was found
	 */		
	public String getOncSupAnnotation(MySQLConnection conn) {
		return conn.getStringFromDb(String.format(ONC_SUP_QUERY,this.getGeneName()));
	}
	
	/**
	 * Loads the substructure information for this gene from the default database table.
	 * All previously loaded substructures will be overwritten.
	 * @param onlyPredicted if true, only predicted substructures will be loaded
	 * @return the number of substructures found
	 * @throws SQLException 
	 */
	public int loadSubstructures(MySQLConnection conn, boolean onlyPredicted) throws SQLException {
		Collection<Substructure> ssl = new LinkedList<Substructure>();
		
		Statement s = conn.createStatement();
		ResultSet rs  = s.executeQuery(String.format(SUBSTRUCTURE_QUERY, this.getGeneName()));
			while(rs.next()) {
				// Result: pos_beg, pos_end, pdb_code, chain_code, type
				int posBeg = rs.getInt(1);
				int posEnd = rs.getInt(2);
				String pdbCode = rs.getString(3);
				String chainCode = rs.getString(4);
				String typeStr = rs.getString(5);
				SubstructureType type;
				if(typeStr.toUpperCase().startsWith("X")) type = SubstructureType.XRAY; 
					else if(typeStr.toUpperCase().startsWith("N")) type = SubstructureType.NMR;
					else if(typeStr.toUpperCase().startsWith("P")) type = SubstructureType.PREDICTION;
					else type = SubstructureType.OTHER;
				Substructure ss = new Substructure(posBeg, posEnd, type, pdbCode, chainCode);
				if(!onlyPredicted || ss.getType() == SubstructureType.PREDICTION) {
					ssl.add(ss);
				}
			}		
		this.subStructures = ssl;
		return ssl.size();
	}

// OBSOLETE:
//	/**
//	 * Loads the substructure information for this gene from the default file (Google Docs table).
//	 * All previously loaded substructures will be overwritten.
//	 * @param onlyExperimental if true, predicted substructures will not be loaded
//	 * @return the number of substructures found
//	 */
//	public int loadSubstructures_old(boolean onlyExperimental) {
//		Collection<Substructure> ssl = new LinkedList<Substructure>();
//		try {
//			BufferedReader in = new BufferedReader(new FileReader(SUBSTRUCTURE_FILE));
//			String line = in.readLine(); // skip header line
//			while((line=in.readLine()) != null) {
//				String[] fields = line.split("\t");
//				String gene = fields[0];
//				String range = fields[1];
//				String typeStr = fields[3];
//				String pdbCodeRaw = fields[4];
//				String templatesStr = null;
//				if(fields.length > 5) templatesStr = fields[5];
//				Substructure.SubstructureType type = Substructure.parseType(typeStr);
//				String[] startEnd = range.split("-");
//				int start = Integer.parseInt(startEnd[0]);
//				int end = Integer.parseInt(startEnd[1]);
//				String pdbCode = null;
//				String chainCode = null;
//				if(pdbCodeRaw.length() >= 4) pdbCode = pdbCodeRaw.substring(0, 4);
//				if(pdbCodeRaw.length() > 4 && Character.isLetter(pdbCodeRaw.charAt(4))) chainCode = pdbCodeRaw.substring(4, 5);
//				if((type == SubstructureType.XRAY || type == SubstructureType.NMR) && pdbCodeRaw.length() < 5) {
//					System.err.println("Warning. Expected five letter PDB code but found: " + pdbCodeRaw);
//				}
//				if(type == SubstructureType.PREDICTION) chainCode = "A";
//				if(gene.equals(this.geneName)) {
//					Substructure ss = new Substructure(start, end, type, pdbCode, chainCode);
//					if(ss.getType() == SubstructureType.PREDICTION) {
//						ss.setTemplatesString(templatesStr);
//						//ss.setOffset(ss.start-1);	// TODO: Assuming file has been numbered according to Uniprot sequence, (as is 
//													//       the case for our predictions), otherwise results will be undefined.
//					} else {
//						// determine offset between Uniprot sequence and Pdbase sequence (now obsoleted by internal alignment)
////						int offset = getOffsetFromSifts(pdbCode, chainCode);
////						ss.setOffset(offset);			
//					}
//					if(type != SubstructureType.PREDICTION || onlyExperimental == false) {
//						ssl.add(ss);
//						//System.out.print(ss + " ");
//					}
//				}				
//			}
//			in.close();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		this.subStructures = ssl;
//		return ssl.size();
//	}
	
	/**
	 * Returns true if the method to load substructures for this gene has previously been called
	 * (not neccessarily that substructures were actually loaded), and false otherwise.
	 * @return true if substructure have been loaded, false otherwise
	 */
	public boolean areSubstructuresLoaded() {
		return this.subStructures != null;
	}
	
	/**
	 * Loads structural domain definitions from pDomains service (http://pdomains.sdsc.edu) and stores them as StructuralDomainFeatures
	 * @param domType the method to use to assign structural domain definitions
	 * @return the number of domains loaded
	 */
	public int loadStructuralDomains(StructuralDomainType domainType) {
		int domainsLoaded = 0;
		if(this.features == null) this.features = new HashMap<FeatureType, Collection<Feature>>();
		Collection<Feature> domainCollection = new LinkedList<Feature>();
		System.out.print(" pDomains\t");
		for(Substructure ss:this.getSubstructures()) {
			if(ss.type == SubstructureType.PREDICTION) {
				// assume that structures are single domain
				System.out.println("Assuming structure " + ss.getPdbCode() + " is single domain.");
				IntervalSet int1 = IntervalSet.createFromInterval(new Interval(ss.getBegPos(),ss.getEndPos()));
				String domName = ss.getPdbCode() + "A1";
				domainCollection.add(new StructuralDomainFeature(domainType,domName,int1));
				domainsLoaded++;
			} else  {	// this does not work for structures not in the PDB
//				if(ss.getPdbCode().equals("2rd0")) {
//					System.out.println("DEBUG " + ss.getPdbCode() + ":" + ss.getBegPos() + "-" + ss.getEndPos());
//				}
				// workarounds for structures which are missing or wrongly annotated in pDomains
				if(ss.getPdbCode().equals("3ny5") && ss.getChainCode().equals("A")) {
					IntervalSet newIntervals = IntervalSet.createFromInterval(new Interval(ss.getBegPos(),ss.getEndPos()));
					domainCollection.add(new StructuralDomainFeature(domainType,"3NY5A1",newIntervals));
					domainsLoaded++;
				} else 
				if(ss.getPdbCode().equals("2x2u") && ss.getChainCode().equals("A")) {
					IntervalSet int1 = IntervalSet.createFromInterval(new Interval(29,154));
					domainCollection.add(new StructuralDomainFeature(domainType,"2X2UA1",int1));
					IntervalSet int2 = IntervalSet.createFromInterval(new Interval(155,274));
					domainCollection.add(new StructuralDomainFeature(domainType,"2X2UA2",int2));					
					domainsLoaded+=2;
				} else
				if(ss.getPdbCode().equals("3na3") && ss.getChainCode().equals("A")) {
					IntervalSet int1 = IntervalSet.createFromInterval(new Interval(3,208));
					domainCollection.add(new StructuralDomainFeature(domainType,"3NA3A1",int1));
					IntervalSet int2 = IntervalSet.createFromInterval(new Interval(209,299));
					domainCollection.add(new StructuralDomainFeature(domainType,"3NA3A2",int2));					
					domainsLoaded+=2;
				} else
				// the following (and more) were wrong in pDomains before we switched to the NCBI method
//				if(ss.getPdbCode().equals("2x1w") && ss.getChainCode().equals("L")) {
//					domainCollection.add(new StructuralDomainFeature(domainType,"2X1WL1",IntervalSet.createFromInterval(new Interval(123,219))));
//					domainCollection.add(new StructuralDomainFeature(domainType,"2X1WL2",IntervalSet.createFromInterval(new Interval(220,326))));
//					domainsLoaded+=2;
//				} else 
//				if(ss.getPdbCode().equals("3kvq") && ss.getChainCode().equals("A")) {
//					domainCollection.add(new StructuralDomainFeature(domainType,"3KVQA1",IntervalSet.createFromInterval(new Interval(667,756))));
//					domainsLoaded++;
//				} else				
				{
					PDomainsConnection pDomains = new PDomainsConnection();
					Map<String, IntervalSet> domainMap = pDomains.getDomains(ss.getPdbCode(), ss.getChainCode(), domainType);
					if(domainMap.size() > 0) {
						for(String domName:domainMap.keySet()) {
							IntervalSet intervals = domainMap.get(domName);	// intervals in pdb (cif) coordinates
							IntervalSet newIntervals = new IntervalSet();	// new intervals in Uniprot coordinates
							for(Interval domInt:intervals) {
								if(domInt.end > ss.getPdb().getFullLength()) {
									System.err.println("Error: Domain interval " + domInt + " is outside of protein length (" + ss.getPdb().getFullLength() + ") for " + ss.getPdbCode()+ss.getChainCode());									
								}
								int begIdx = domInt.beg;
								int endIdx = domInt.end;
								int domStart = Substructure.ALIGNMENT_UNDEFINED;
								while(domStart == Substructure.ALIGNMENT_UNDEFINED && begIdx < domInt.end) {		// if starting residue is unobserved, try next
									domStart = ss.mapCifResser2Uniprot(begIdx++);
								}
								int domEnd = Substructure.ALIGNMENT_UNDEFINED;
								while(domEnd == Substructure.ALIGNMENT_UNDEFINED && endIdx > domInt.beg) {		// // if final residue is unobserved, try previous
									domEnd = ss.mapCifResser2Uniprot(endIdx--);
								}
								if(domStart == Substructure.ALIGNMENT_UNDEFINED || domEnd == Substructure.ALIGNMENT_UNDEFINED) {
									System.err.println("Warning: Domain " + domInt.beg + "-" + domInt.end + " in " + ss.getPdbCode()+ss.getChainCode() + " could not be mapped to Uniprot sequence.");
								} else {
									// add to intervals
									newIntervals.add(new Interval(domStart,domEnd));
								}
							}
							// create feature and add to domainCollection
							domainCollection.add(new StructuralDomainFeature(domainType,domName,newIntervals));
							// DEBUG: System.out.println(new StructuralDomainFeature(domainType,domName,newIntervals));
							domainsLoaded++;
						}
					} else {
						System.err.println("No domains found for " + ss.getPdbCode() + ss.getChainCode());
					}
				}
			}
		}
		// add domain collection to features
		features.put(FeatureType.SDOMAIN, domainCollection);
		System.out.print("(" + domainsLoaded + ")\t");
		for(Feature f: domainCollection) {
			System.out.print("[" + f + "] ");
		}
		System.out.println();
		return domainsLoaded;
	}
	
	/**
	 * Load selected annotations for this gene from online Uniprot database and store as features.
	 * @return the number of features loaded
	 */
	public int loadUniprotFeatures() {
		int featuresLoaded = 0;
		if(this.features == null) this.features = new HashMap<FeatureType, Collection<Feature>>();
		
		// load features from Uniprot
		UniProtConnection uc = new UniProtConnection();
		System.out.print(" Uniprot\t");
		UniProtEntry upe;
		try {
			upe = uc.getEntry(this.uniProtID);
		} catch (NoMatchFoundException e1) {
			System.err.println("Uniprot entry not found: " + this.uniProtID);
			return 0;
		}
		for(uk.ac.ebi.kraken.interfaces.uniprot.features.Feature hit:upe.getFeatures()) {
			String key = hit.getType().getName();
			if(key.equals("NP_BIND") || key.equals("ACT_SITE") ||	// ATP/GTP binding, enzyme active site
			   key.equals("MOD_RES") || key.equals("CARBOHYD")) {	// post-transl. mod., glycosilation site
				String description = ((HasFeatureDescription) hit).getFeatureDescription().getValue();
				int begin = hit.getFeatureLocation().getStart();
				int end = hit.getFeatureLocation().getEnd();
				try {
					//System.out.print(".");
					if(this.addFeature(new UniprotFeature(begin, end, key, description))) featuresLoaded++;
				} catch (InvalidFeatureCoordinatesException e) {
				} catch (OverlappingFeatureException e) {
				}
			}
		}
		System.out.println("(" + featuresLoaded + ")");
		return featuresLoaded;
	}
	
	/**
	 * Load selected annotations for this gene using an online Prosite scan and store hits as features.
	 * @return the number of features loaded
	 */
	public int loadPrositeFeatures() {
		int featuresLoaded = 0;
		if(this.features == null) this.features = new HashMap<FeatureType, Collection<Feature>>();
				
		// load features from ProSite
		PrositeScanner ps = new PrositeScanner();
		System.out.print(" Prosite\t");
		ps.scanSeq(this.getUniprotId());
		for(PrositeHit hit:ps.getLatestHits()) {
			if(hit.signatureAc.equals("PS00107") || hit.signatureAc.equals("PS00383")    // Kinase ATP binding region, tyrosine phosphatase active site
			|| hit.signatureAc.equals("PS00109") || hit.signatureAc.equals("PS00108")) { // Protein kinase active-site
				try {
					//System.out.print(".");
					if(this.addFeature(new PrositeFeature(hit))) featuresLoaded++;
				} catch (InvalidFeatureCoordinatesException e) {
				} catch (OverlappingFeatureException e) {
				}
			}
		}
		System.out.println("(" + featuresLoaded + ")");
		
		// print features (for debugging)
//		for(Feature f:this.getFeatures()) {
//			System.out.println(f.toString());
//		}
		return featuresLoaded;
	}	
	
	/**
	 * Load selected annotations for this gene by parsing a PhosphoSitePlus html file.
	 * @param htmlFile an HTML file downloaded from PhosphoSitePlus.
	 * @return the number of features loaded
	 */
	public int loadPhosphoSitePlusFeatures(File htmlFile) {
		int featuresLoaded = 0;
		if(this.features == null) this.features = new HashMap<FeatureType, Collection<Feature>>();
		System.out.print(" PhosphoSite\t");
		HashSet<Integer> alreadyAdded = new HashSet<Integer>();
		for(PhosphoSitePlusFeature m:PhosphoSiteConnection.getModifications(htmlFile)) {
			// check amino acid
			if(m.getPosition() > this.uniprotSeq.length()) {
				System.err.println("Invalid feature index " + m + " (seqlen=" + this.uniprotSeq.length() + "). Skipping.");
				continue;
			}
			if(m.getAminoAcid() != AminoAcid.getByOneLetterCode(this.uniprotSeq.charAt(m.getPosition()-1))) {
				System.err.println("Invalid amino acid " + m + " (seq[i]=" + this.uniprotSeq.charAt(m.getPosition()-1) + "). Skipping.");
				continue;				
			}
			try {
				if(!alreadyAdded.contains(m.getPosition()))	{ // filter out identical positions
					if(this.addFeature(m)) {
						//System.out.print(".");
						//System.out.print(m + " ");
						alreadyAdded.add(m.getPosition());
						featuresLoaded++;
					}
				}
			} catch (InvalidFeatureCoordinatesException e) {
				e.printStackTrace();
			} catch (OverlappingFeatureException e) {
				e.printStackTrace();
			}
		}
		System.out.println("(" + featuresLoaded + ")");
		return featuresLoaded;		
	}
	
	public int loadCSAFeatures() {
		int featuresLoaded = 0;
		if(this.features == null) this.features = new HashMap<FeatureType, Collection<Feature>>();
		System.out.print(" CatSites(CSA)\t");
		if(this.subStructures == null) {
			System.err.println("Error loading features from CSA. Substructure not loaded for " + this.geneName);
			return 0;
		}
		
		for(Substructure ss:this.getSubstructures()) {
			if(ss.getType() != SubstructureType.PREDICTION) {	// otherwise we don't have a PDB code
				try {
					if(CSAConnection.parseCSA(ss.pdb, CatalSiteSet.LATEST_VERSION, USE_ONLINE_CSA) > 0) {
						System.err.println("Warning: Residue mismatches detected when loading CSA annotation for " + this.geneName + " " + ss.pdbCode);
					}
					if(ss.pdb.getCSA() == null) {
						//System.err.println("Could not find CSA annotations for " + this.geneName + " " + ss.pdbCode);
					} else {
						Iterator<CatalyticSite> it = ss.pdb.getCSA().getIterator();
						while(it.hasNext()) {
							CatalyticSite cs = it.next();
							// modify set to remap to Uniprot sequence
							Set<Integer> bakSet = new HashSet<Integer>();
							Map<Integer, String> newSet = new HashMap<Integer,String>();
							Set<Integer> resSet = cs.getRes();
							bakSet.addAll(resSet);	// backup copy of original set
							for(int pos: resSet) {
								newSet.put(ss.mapCifResser2Uniprot(pos),cs.getChemFuncFromResSerial(pos));
							}
							// remove all old residues
							for(int pos:bakSet) {
								cs.remRes(pos);
							}
							// put back residues with new number
							for(int pos:newSet.keySet()) {
								cs.addRes(pos, newSet.get(pos));
							}
							// add feature
							try {
								if(this.addFeature(new CsaFeature(cs))) {
									//System.out.print(cs + " ");
									//System.out.print(".");
									featuresLoaded++;
								}
							} catch (InvalidFeatureCoordinatesException e) {
								e.printStackTrace();
							} catch (OverlappingFeatureException e) {
								e.printStackTrace();
							}
						}
					}
				} catch (IOException e) {
					System.err.println("Error accessing CSA: " + e.getMessage());
				}
			}
		}
		System.out.print("(" + featuresLoaded + ")\t");
		for(Feature f: this.getFeaturesOfType(FeatureType.CSA)) {
			System.out.print("[" + f + "] ");
		}
		System.out.println();
		return featuresLoaded;	
	}
	
	public int loadManualFeatures() {
		int featuresLoaded = 0;
		if(this.features == null) this.features = new HashMap<FeatureType, Collection<Feature>>();
		System.out.print(" ManualAnno\t");
		try {
			BufferedReader in = new BufferedReader(new FileReader(ANNOTATION_FILE));
			String line;
			while((line = in.readLine()) != null) {
				if(line.startsWith(this.getGeneName())) {
					String[] fields = line.split("\t");
					String posStr = fields[1];
					String cat = fields[2];
					//String comment = fields[3];
					if(cat.equals("ACT_SITE") || cat.equals("ATP_BIND") || cat.equals("GTP_BIND")  
											  || cat.equals("PHO_RES")  || cat.equals("MOD_RES") || cat.equals("UBQ_RES")
											  || cat.equals("CARBOHYD") || cat.equals("DNA_BIND") || cat.equals("ATY_RES")) {
						IntervalSet pos = Interval.getIntervals(Interval.parseSelectionString(posStr));
						Feature f = new GeneralFeature(pos, cat);
						try {
							if(this.addFeature(f)) {
								featuresLoaded++;
							}
						} catch (InvalidFeatureCoordinatesException e) {
							e.printStackTrace();
						} catch (OverlappingFeatureException e) {
							e.printStackTrace();
						}
					}
				}
			}
		} catch (IOException e) {
			System.err.println("Error reading from file " + ANNOTATION_FILE + ":" + e.getMessage());
		}
		
		System.out.println("(" + featuresLoaded + ")");
		return featuresLoaded;
	}
	
	/**
	 * Returns a collection of the loaded substructures which have mutations.
	 * @return
	 */
	public Collection<Substructure> getSubstructuresWithMutations() {
		HashSet<Substructure> newSet = new HashSet<Substructure>();
		for(Mutation m:this.mutations) {
			Substructure ss = this.getSubstructure(m.position);
			if(ss != null) newSet.add(ss);
		}
		return newSet;
	}
	
	/**
	 * Of the currently loaded mutations returns the ones which occur in known structures and at observed
	 * positions. Returns null if mutations are not loaded or substructures are not loaded.
	 * @param missenseOnly if true, consider only missense mutations
	 * @param experimentalOnly if true, consider only XRAY and NMR stuctures.
	 * @return
	 */
	public Collection<Mutation> getMutationsInStructures(boolean missenseOnly, boolean experimentalOnly) {
		Collection<Mutation> strMut = new LinkedList<Mutation>();
		if(this.mutations == null) {
			return null;
		}
		if(!this.areSubstructuresLoaded()) {
			return null;
		}
		if(this.mutations.size() > 0) {
			for(Mutation m:this.getMutations()) {
				if(!missenseOnly || m.getType() == Mutation.MutType.MISSENSE) {
					if(this.getSubstructure(m.position) != null) {
						Substructure ss = this.getSubstructure(m.position);
						if(!experimentalOnly || ss.getType() == Substructure.SubstructureType.NMR || ss.getType() == Substructure.SubstructureType.XRAY) {
							if(ss.pdb.hasCoordinates(ss.mapUniprotResser2Cif(m.position))) strMut.add(m);							
						}
					}
				}
			}
		}
		return strMut;
	}

// OBSOLETE:
//	/**
//	 * Returns the offset between Uniprot sequence and Pdbase sequence read from the official online SIFTS mapping
//	 * @param pdbCode
//	 * @param chainCode
//	 * @return
//	 * @throws IOException 
//	 */
//	public static int getOffsetFromSifts(String pdbCode, String chainCode) throws IOException {
//		int offset = Substructure.OFFSET_UNDEFINED;
//		try {
//			URL url = new URL(SIFTS_URL);
//			BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
//			String line;
//			while((line = in.readLine()) != null) {
//				if(line.startsWith(pdbCode)) {
//					//DEBUG: System.out.println(line);
//					String[] fields = line.split("\t");
//					//DEBUG: System.out.println(fields);
//					String pdb = fields[0];
//					String chain = fields[1];
//					//String uniprotId = fields[2];
//					int cifBeg = Integer.parseInt(fields[4]);
//					int cifEnd = Integer.parseInt(fields[5]);				
//					//int pdbBeg = Integer.parseInt(fields[6]);
//					//int pdbEnd = Integer.parseInt(fields[7]);				
//					int uniBeg = Integer.parseInt(fields[8]);
//					int uniEnd = Integer.parseInt(fields[9]);
//					if(pdb.equals(pdbCode) && chain.equals(chainCode)) {
//						in.close();
//						offset = uniBeg - cifBeg;
//						if(uniEnd-uniBeg != cifEnd-cifBeg) {
//							System.err.println("Warning: SIFTS mapping reports different lenght of sequence between Uniprot and PDBase for " + pdbCode+" "+chainCode);
//						}
//						return offset;
//					}
//				}
//			}
//		} catch (MalformedURLException e) {
//			e.printStackTrace();
//		}
//		System.err.println("Warning: " + pdbCode+" "+chainCode + " not found in SIFTS mapping at " + SIFTS_URL);
//		return offset;
//	}
	
	/**
	 * Loads mutation information for this gene from the database.
	 * @param onlyMissense if true, only missense mutations are loaded
	 * @return the number of mutations loaded.
	 */
	public int loadMutations(MySQLConnection conn, boolean onlyMissense) {
		Collection<Mutation> muts = new LinkedList<Mutation>();
		TreeSet<Integer> positions = new TreeSet<Integer>();
		boolean unique = LOAD_UNIQUE;	// if true, discard mutations with identical positions
		try {
			Statement s = conn.createStatement();
			ResultSet rs = s.executeQuery(String.format(MUTATION_QUERY,this.geneName));
			while(rs.next()) {
				String[] fields = rs.getString(1).split(";");
				String dnaMutStr = fields[0];
				String mutStr = fields[1];
				String query = String.format(MUT_ID_QUERY, this.geneName, dnaMutStr, mutStr);
				String mutIdStr = conn.getStringFromDb(query);
				//System.out.println(query);
				//System.exit(1);
				if(mutStr.startsWith("p.")) {
					String mutSubStr = mutStr.substring(2);
					Mutation m = Mutation.parseMutation(mutSubStr);
					m.setDnaMutStr(dnaMutStr);
					m.setMutId(Integer.parseInt(mutIdStr));
					m.position = this.mapCosmicPos2Uniprot(m.position);
					if(m != null && m.position > 0 && (!onlyMissense || m.isMissense())) {
						if(!unique || !positions.contains(m.position)) {
							positions.add(m.position);
							muts.add(m);
						}
					} else {
						System.err.println("Error parsing mutation:" + mutSubStr);
					}
				} else {
					System.err.println("Error parsing mutation. Expected 'p.', found: " + mutStr);
				}
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		this.mutations = muts;
		return muts.size();
	}
	
	// alternative method name to keep old code working:
//	/**
//	 * See: {@link loadMutationsFromDb()}
//	 */
//	public int loadMutationsFromDb(MySQLConnection conn) {
//		return loadMutations(conn);
//	}
	
	// Keep only for the unlikely case that we want to reproduce some old data:
//	/**
//	 * Loads the mutations associated with this genes from the default file (Google Docs table).
//	 * All previously loaded mutations will be overwritten.
//	 * @return the number of mutations found
//	 */
//	public int loadMutations_old() {
//		Collection<Mutation> muts = new LinkedList<Mutation>();
//		TreeSet<Integer> positions = new TreeSet<Integer>();
//		boolean unique = LOAD_UNIQUE;	// if true, discard mutations with identical positions
//		try {
//			BufferedReader in = new BufferedReader(new FileReader(MUTATION_FILE));
//			String line = in.readLine(); // skip header line
//			while((line=in.readLine()) != null) {
//				String[] fields = line.split("\t");
//				String gene = fields[0];
//				String mutationStr = fields[1];
//				String mutationConseq = fields[7];
//				String mutationConseqDetails = fields[8];
//				if(gene.equals(this.geneName)) {
//					Mutation m = Mutation.parseMutation(mutationStr);
//					if(m != null) {
//						m.setParent(this);
//						m.setComment(mutationConseq);
//						m.setCommentDetails(mutationConseqDetails);
//						m.setMutStr(mutationStr);
//						if(!unique || !positions.contains(m.position)) {
//							positions.add(m.position);
//							muts.add(m);
//						}
//					} else {
//						System.err.println("Failed to parse mutation string: " + mutationStr + " in file " + MUTATION_FILE);
//					}
//				}				
//			}
//			in.close();
//			this.mutations = muts;
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		return muts.size();
//	}
	
// OBSOLETE:
//	/**
//	 * Loads mutation information for this gene from the default database (mutanom.missense_mutation).
//	 * @return the number of mutations loaded.
//	 */
//	public int loadMutationsFromDb(MySQLConnection conn) {
//		Collection<Mutation> muts = new LinkedList<Mutation>();
//		TreeSet<Integer> positions = new TreeSet<Integer>();
//		boolean unique = LOAD_UNIQUE;	// if true, discard mutations with identical positions
//		try {
//			Statement s = conn.createStatement();
//			ResultSet rs = s.executeQuery(String.format(MUTATION_QUERY,this.geneName));
//			while(rs.next()) {
//				String mutStr = rs.getString(1);
//				if(mutStr.startsWith("p.")) {
//					String mutSubStr = mutStr.substring(2);
//					Mutation m = Mutation.parseMutation(mutSubStr);
//					m.position = this.mapCosmicPos2Uniprot(m.position);
//					if(m != null && m.position > 0) {
//						if(!unique || !positions.contains(m.position)) {
//							positions.add(m.position);
//							muts.add(m);
//						}
//					} else {
//						System.err.println("Error parsing mutation:" + mutSubStr);
//					}
//				} else {
//					System.err.println("Error parsing mutation. Expected 'p.', found: " + mutStr);
//				}
//			}
//		} catch (SQLException e) {
//			e.printStackTrace();
//		}
//		this.mutations = muts;
//		return muts.size();
//	}
//	
//	/**
//	 * Loads mutation information for this gene from the database mutanom.cosmic43_mis_4tis;
//	 * @return the number of mutations loaded.
//	 */
//	public int loadMutationsFromDb2(MySQLConnection conn) {
//		Collection<Mutation> muts = new LinkedList<Mutation>();
//		TreeSet<Integer> positions = new TreeSet<Integer>();
//		boolean unique = LOAD_UNIQUE;	// if true, discard mutations with identical positions
//		try {
//			Statement s = conn.createStatement();
//			ResultSet rs = s.executeQuery(String.format(MUTATION_QUERY2,this.geneName));
//			while(rs.next()) {
//				String mutStr = rs.getString(1);
//				if(mutStr.startsWith("p.")) {
//					String mutSubStr = mutStr.substring(2);
//					Mutation m = Mutation.parseMutation(mutSubStr);
//					m.position = this.mapCosmicPos2Uniprot(m.position);
//					if(m != null && m.position > 0) {
//						if(!unique || !positions.contains(m.position)) {
//							positions.add(m.position);
//							muts.add(m);
//						}
//					} else {
//						System.err.println("Error parsing mutation:" + mutSubStr);
//					}
//				} else {
//					System.err.println("Error parsing mutation. Expected 'p.', found: " + mutStr);
//				}
//			}
//		} catch (SQLException e) {
//			e.printStackTrace();
//		}
//		this.mutations = muts;
//		return muts.size();
//	}
//
//	/**
//	 * Loads mutation information for this gene from the database mutanom.cosmic43_mis_5tis;
//	 * @return the number of mutations loaded.
//	 */
//	public int loadMutationsFromDb3(MySQLConnection conn) {
//		Collection<Mutation> muts = new LinkedList<Mutation>();
//		TreeSet<Integer> positions = new TreeSet<Integer>();
//		boolean unique = LOAD_UNIQUE;	// if true, discard mutations with identical positions
//		try {
//			Statement s = conn.createStatement();
//			ResultSet rs = s.executeQuery(String.format(MUTATION_QUERY3,this.geneName));
//			while(rs.next()) {
//				String mutStr = rs.getString(1);
//				if(mutStr.startsWith("p.")) {
//					String mutSubStr = mutStr.substring(2);
//					Mutation m = Mutation.parseMutation(mutSubStr);
//					m.position = this.mapCosmicPos2Uniprot(m.position);
//					if(m != null && m.position > 0) {
//						if(!unique || !positions.contains(m.position)) {
//							positions.add(m.position);
//							muts.add(m);
//						}
//					} else {
//						System.err.println("Error parsing mutation:" + mutSubStr);
//					}
//				} else {
//					System.err.println("Error parsing mutation. Expected 'p.', found: " + mutStr);
//				}
//			}
//		} catch (SQLException e) {
//			e.printStackTrace();
//		}
//		this.mutations = muts;
//		return muts.size();
//	}
	
	/**
	 * Verifies that all missense mutations for which structures are known have the correct residue at the respective position
	 * in the structure. Return true if this is the case, false otherwise. Can also be used for SNPs with loadSnpsAsMutations().
	 * @return true iff all mutations which can be found in a loaded structure have the correct residue at the respective position 
	 */
	public boolean checkMutations() {
		boolean result = true;
		if(this.mutations == null) {
			System.err.println("Error: No mutations loaded for " + this.geneName);
			return false;
		}
		if(this.subStructures == null) {
			System.err.println("Error: No substructures loaded for " + this.geneName);
			return false;			
		}
		for(Mutation m:this.mutations) {
				//System.out.println(m);
				Substructure ss = this.getSubstructure(m.position);
				if(ss != null) {
					if(!ss.isPdbLoaded()) {
						System.err.println("Error: No pdb loaded for " + this.geneName + " and " + m.toString());
						return false;
					}
					if(!ss.isAlignmentInitialized()) {
						System.err.println("Error: Alignment not initialized for " + this.geneName + " and " + ss.toString());
						return false;
					}
					int uniPos = m.position;
					int pdbPos = ss.mapUniprotResser2Cif(uniPos);
					String pdbSeq = ss.pdb.getSequence();
					AminoAcid mutAA = m.before;
					if(pdbSeq.length() < pdbPos) {
						System.err.printf("Error: Mutation %s%d at mapped position %d can not be found in sequence of length %d in %s %s (offset=%d)\n", mutAA.getOneLetterCode(), uniPos, pdbPos, pdbSeq.length(), ss.getPdbCode(), ss.getChainCode(), ss.offset);
					} else {
						AminoAcid pdbAA = AminoAcid.getByOneLetterCode(pdbSeq.charAt(pdbPos-1));
						if(mutAA != pdbAA) {
							System.out.printf("Warning: Mutation %s%d does not match %s%d in %s%s\n", mutAA.getOneLetterCode(), uniPos, pdbAA.getOneLetterCode(), pdbPos, ss.getPdbCode(), ss.getChainCode());
							result = false;
						}
						if(!ss.pdb.hasCoordinates(pdbPos)) {
							System.out.println("Warning: Residue at mutated position " + pdbPos + " has no coordinates");
						}
					}
				}
		}
		return result;
	}
	
	/**
	 * Loads SNP data for this gene from the database.
	 * Previously loaded SNPs will be overwritten.
	 * @return the number of snps found.
	 */
	public int loadSNPs(MySQLConnection conn) {
		Collection<SNP> snps = new LinkedList<SNP>();
		TreeSet<Integer> positions = new TreeSet<Integer>();
		boolean unique = LOAD_UNIQUE;	// if true, discard mutations with identical positions
		try {
			Statement s = conn.createStatement();
			ResultSet rs = s.executeQuery(String.format(SNP_QUERY,this.geneName));
			while(rs.next()) {
				String snpIdStr = rs.getString(1);
				int snpId = Integer.parseInt(snpIdStr.substring(2));
				String refAA = rs.getString(2);
				String mutAA = rs.getString(3);
				int pos = rs.getInt(4);
				//double freq = rs.getDouble(5);
				double freq = Double.NaN;
				SNP snp = new SNP(AminoAcid.getByOneLetterCode(refAA.charAt(0)), AminoAcid.getByOneLetterCode(mutAA.charAt(0)), pos, snpId, freq);
				if(!unique || !positions.contains(snp.position)) {
					positions.add(snp.position);
					snps.add(snp);
				}
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		this.snps = snps;
		return snps.size();		
	}
	
	/**
	 * Loads SNP data for this gene from the database but store in this.mutations.
	 * This allows to apply methods only implemented for mutations to be used for
	 * SNPs (e.g. for testing purposes). Note: SNP is a subclass of Mutation.
	 * Previously loaded mutations will be overwritten.
	 * @return the number of snps found.
	 */
	public int loadSNPsAsMutations(MySQLConnection conn) {
		Collection<Mutation> snps = new LinkedList<Mutation>();
		TreeSet<Integer> positions = new TreeSet<Integer>();
		boolean unique = LOAD_UNIQUE;	// if true, discard mutations with identical positions
		try {
			Statement s = conn.createStatement();
			ResultSet rs = s.executeQuery(String.format(SNP_QUERY,this.geneName));
			while(rs.next()) {
				String snpIdStr = rs.getString(1);
				int snpId = Integer.parseInt(snpIdStr.substring(2));
				String refAA = rs.getString(2);
				String mutAA = rs.getString(3);
				int pos = rs.getInt(4);
				//double freq = rs.getDouble(5);
				double freq = Double.NaN;
				SNP snp = new SNP(AminoAcid.getByOneLetterCode(refAA.charAt(0)), AminoAcid.getByOneLetterCode(mutAA.charAt(0)), pos, snpId, freq);
				if(!unique || !positions.contains(snp.position)) {
					positions.add(snp.position);
					snps.add(snp);
				}
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		this.mutations = snps;
		return snps.size();		
	}
		
	// alternative method name to keep old code working:
//	/**
//	 * See: {@link loadSNPS()}
//	 */
//	public int loadSnpsFromDb(MySQLConnection conn) {
//		return loadMutations(conn);
//	}	
	
	
// OBSOLETE:
//	/**
//	 * Loads SNP data for this gene from the default database (mutanom.snp).
//	 * Previously loaded SNPs will be overwritten.
//	 * @return the number of snps found.
//	 */
//	public int loadSnpsFromDb(MySQLConnection conn) {
//		Collection<Mutation> snps = new LinkedList<Mutation>();
//		TreeSet<Integer> positions = new TreeSet<Integer>();
//		boolean unique = LOAD_UNIQUE;	// if true, discard mutations with identical positions
//		try {
//			Statement s = conn.createStatement();
//			String ensp = conn.getStringFromDb(String.format(ENSP_QUERY,this.uniProtID));
//			if(ensp != null) {
//				ResultSet rs = s.executeQuery(String.format(SNP_QUERY,ensp));
//				while(rs.next()) {
//					int pos = rs.getInt(1);
//					String refAA = rs.getString(2);
//					String mutAA = rs.getString(3);
//					Mutation m = new Mutation(AminoAcid.getByThreeLetterCode(refAA), AminoAcid.getByThreeLetterCode(mutAA), pos);
//					if(!unique || !positions.contains(m.position)) {
//						positions.add(m.position);
//						snps.add(m);
//					}
//				}
//			}
//		} catch (SQLException e) {
//			e.printStackTrace();
//		}
//		this.snps = snps;
//		return snps.size();		
//	}
//	
//	/**
//	 * Loads SNP data for this gene from the an alternative database (mutanom.ensvar56_top21).
//	 * Previously loaded SNPs will be overwritten.
//	 * @return the number of snps found.
//	 */
//	public int loadSnpsFromDb2(MySQLConnection conn) {
//		Collection<Mutation> snps = new LinkedList<Mutation>();
//		TreeSet<Integer> positions = new TreeSet<Integer>();
//		boolean unique = LOAD_UNIQUE;	// if true, discard mutations with identical positions
//		try {
//			String ret = conn.getStringFromDb(String.format(ENSG_QUERY,this.uniProtID));
//			if(ret != null) {
//				String[] ensgs = ret.split("; ");
//				String ensg = ensgs[0];
//				Statement s = conn.createStatement();
//				ResultSet rs = s.executeQuery(String.format(SNP_QUERY2,ensg));
//				int mismatches = 0;
//				while(rs.next()) {
//					int pos = rs.getInt(1);
//					String mutStr = rs.getString(2);
//					if(this.getUniprotSeq().length() < pos) {
//						//System.out.print(pos + ">" + g.getUniprotSeq().length() + " ");
//						mismatches++;
//					} else if (this.getUniprotSeq().charAt(pos-1) != mutStr.charAt(0)) {
//						//System.out.print(g.getUniprotSeq().charAt(pos-1) +"!="+ mutStr.charAt(0)+ " ");
//						mismatches++;
//					} else {
//						AminoAcid refAA = AminoAcid.getByOneLetterCode(mutStr.charAt(0));
//						AminoAcid mutAA = AminoAcid.getByOneLetterCode(mutStr.charAt(2));
//						Mutation m = new Mutation(refAA, mutAA, pos);
//						if(!unique || !positions.contains(m.position)) {
//							positions.add(m.position);
//							snps.add(m);
//						}
//					}
//				}
//			}
//		} catch (SQLException e) {
//			e.printStackTrace();
//		}	
//		this.snps = snps;
//		return snps.size();		
//	}
//	
//	/**
//	 * Loads SNP data for this gene from the an alternative database (mutanom.ensvar56_top21)
//	 * using an alternative mapping of Uniprot to Ensemble (mutanom.ens47_sp2ens). This may
//	 * map multiple ENSTs to one Uniprot entry, so we scan them all for mutations.
//	 * Previously loaded SNPs will be overwritten.
//	 * @return the number of snps found.
//	 */
//	public int loadSnpsFromDb3(MySQLConnection conn) {
//		Collection<Mutation> snps = new LinkedList<Mutation>();
//		TreeSet<Integer> positions = new TreeSet<Integer>();
//		boolean unique = LOAD_UNIQUE;	// if true, discard mutations with identical positions
//		try {
//			String[] ensts = conn.getStringsFromDb(String.format(ENST_QUERY,this.uniProtID));
//			if(ensts != null) {
//				String inClause = "\"" + ensts[0] + "\"";
//				for (int i = 1; i < ensts.length; i++) {
//					inClause += ",\"" + ensts[i] + "\"";
//				}
//				Statement s = conn.createStatement();
//				ResultSet rs = s.executeQuery(String.format(SNP_QUERY3,inClause));
//				int mismatches = 0;
//				while(rs.next()) {
//					int pos = rs.getInt(1);
//					String mutStr = rs.getString(2);
//					if(this.getUniprotSeq().length() < pos) {
//						//System.out.print(pos + ">" + g.getUniprotSeq().length() + " ");
//						mismatches++;
//					} else if (this.getUniprotSeq().charAt(pos-1) != mutStr.charAt(0)) {
//						//System.out.print(g.getUniprotSeq().charAt(pos-1) +"!="+ mutStr.charAt(0)+ " ");
//						mismatches++;
//					} else {
//						AminoAcid refAA = AminoAcid.getByOneLetterCode(mutStr.charAt(0));
//						AminoAcid mutAA = AminoAcid.getByOneLetterCode(mutStr.charAt(2));
//						Mutation m = new Mutation(refAA, mutAA, pos);
//						if(!unique || !positions.contains(m.position)) {
//							positions.add(m.position);
//							snps.add(m);
//						}
//					}
//				}
//			}
//		} catch (SQLException e) {
//			e.printStackTrace();
//		}	
//		this.snps = snps;
//		return snps.size();		
//	}

	
// OBSOLETE:	
//	public int loadSNPsAsMutations(MySQLConnection conn) {
//		Collection<Mutation> snps = new LinkedList<Mutation>();
//		TreeSet<Integer> positions = new TreeSet<Integer>();
//		boolean unique = LOAD_UNIQUE;	// if true, discard mutations with identical positions
//		try {
//			Statement s = conn.createStatement();
//			String ensp = conn.getStringFromDb(String.format(ENSP_QUERY,this.uniProtID));
//			if(ensp != null) {
//				ResultSet rs = s.executeQuery(String.format(SNP_QUERY,ensp));
//				while(rs.next()) {
//					int pos = rs.getInt(1);
//					String refAA = rs.getString(2);
//					String mutAA = rs.getString(3);
//					Mutation m = new Mutation(AminoAcid.getByThreeLetterCode(refAA), AminoAcid.getByThreeLetterCode(mutAA), pos);
//					if(!unique || !positions.contains(m.position)) {
//						positions.add(m.position);
//						snps.add(m);
//					}
//				}
//			}
//		} catch (SQLException e) {
//			e.printStackTrace();
//		}
//		this.mutations = snps;
//		return mutations.size();			
//	}
//	
//	/**
//	 * Load SNPs from alternative dataset taken from EnsemblVariation56 (dbSNP130), mapped through ENSG ID. 
//	 * @param conn
//	 * @return
//	 */
//	public int loadSNPsAsMutations2(MySQLConnection conn) {
//		Collection<Mutation> snps = new LinkedList<Mutation>();
//		TreeSet<Integer> positions = new TreeSet<Integer>();
//		boolean unique = LOAD_UNIQUE;	// if true, discard mutations with identical positions
//		try {
//			String ret = conn.getStringFromDb(String.format(ENSG_QUERY,this.uniProtID));
//			if(ret != null) {
//				String[] ensgs = ret.split("; ");
//				String ensg = ensgs[0];
//				Statement s = conn.createStatement();
//				ResultSet rs = s.executeQuery(String.format(SNP_QUERY2,ensg));
//				int mismatches = 0;
//				while(rs.next()) {
//					int pos = rs.getInt(1);
//					String mutStr = rs.getString(2);
//					if(this.getUniprotSeq().length() < pos) {
//						//System.out.print(pos + ">" + g.getUniprotSeq().length() + " ");
//						mismatches++;
//					} else if (this.getUniprotSeq().charAt(pos-1) != mutStr.charAt(0)) {
//						//System.out.print(g.getUniprotSeq().charAt(pos-1) +"!="+ mutStr.charAt(0)+ " ");
//						mismatches++;
//					} else {
//						AminoAcid refAA = AminoAcid.getByOneLetterCode(mutStr.charAt(0));
//						AminoAcid mutAA = AminoAcid.getByOneLetterCode(mutStr.charAt(2));
//						Mutation m = new Mutation(refAA, mutAA, pos);
//						if(!unique || !positions.contains(m.position)) {
//							positions.add(m.position);
//							snps.add(m);
//						}
//					}
//				}
//
//			}
//		} catch (SQLException e) {
//			e.printStackTrace();
//		}	
//		this.mutations = snps;
//		return mutations.size();			
//	}	
//	
//	
//public int loadSNPsAsMutations3(MySQLConnection conn) {
//	Collection<Mutation> snps = new LinkedList<Mutation>();
//	TreeSet<Integer> positions = new TreeSet<Integer>();
//	boolean unique = LOAD_UNIQUE;	// if true, discard mutations with identical positions
//	try {
//		String[] ensts = conn.getStringsFromDb(String.format(ENST_QUERY,this.uniProtID));
//		if(ensts != null) {
//			String inClause = "\"" + ensts[0] + "\"";
//			for (int i = 1; i < ensts.length; i++) {
//				inClause += ",\"" + ensts[i] + "\"";
//			}
//			Statement s = conn.createStatement();
//			ResultSet rs = s.executeQuery(String.format(SNP_QUERY3,inClause));
//			int mismatches = 0;
//			while(rs.next()) {
//				int pos = rs.getInt(1);
//				String mutStr = rs.getString(2);
//				if(this.getUniprotSeq().length() < pos) {
//					//System.out.print(pos + ">" + g.getUniprotSeq().length() + " ");
//					mismatches++;
//				} else if (this.getUniprotSeq().charAt(pos-1) != mutStr.charAt(0)) {
//					//System.out.print(g.getUniprotSeq().charAt(pos-1) +"!="+ mutStr.charAt(0)+ " ");
//					mismatches++;
//				} else {
//					AminoAcid refAA = AminoAcid.getByOneLetterCode(mutStr.charAt(0));
//					AminoAcid mutAA = AminoAcid.getByOneLetterCode(mutStr.charAt(2));
//					Mutation m = new Mutation(refAA, mutAA, pos);
//					if(!unique || !positions.contains(m.position)) {
//						positions.add(m.position);
//						snps.add(m);
//					}
//				}
//			}
//		}
//	} catch (SQLException e) {
//		e.printStackTrace();
//	}	
//	this.mutations = snps;
//	return snps.size();	
//}
	
	
//	public boolean hasKnownStructures() {
//		return false;
//	}
//	
//	public boolean hasAtpBindingSite() {
//		return false;
//	}
//	
//	public boolean hasDnaBindingSite() {
//		return false;
//	}
//	
//	public boolean hasCatalyticSite() {
//		return false;
//	}
//	
//	public boolean hasTransmembraneRegion() {
//		return false;
//	}
	
	/**
	 * Returns the first substructure found to be covering the given position in this gene's sequence or
	 * null if no such structure exists.
	 * @param pos the residue number for which the substructure is to be returned.
	 * @return the first substructure covering the the given positition
	 */
	public Substructure getSubstructure(int pos) {
		Substructure ss = null;
		for(Substructure s:this.subStructures) {
			if(s != null && pos >= s.start && pos <= s.end) ss = s;
		}
		return ss;
	}
	
	/**
	 * Writes a data file with substructure information.
	 * @param outFile the file to be written
	 * TODO: @param restrictToType if not null, only substructures of the given type will be written (e.g. SubstructureType.XRAY)
	 */
	public void writeSubstructureData(File outFile, SubstructureType restrictToType) {
		try {
			PrintWriter out = new PrintWriter(outFile);
			if (this.subStructures != null) {
				
				for (Substructure ss : this.subStructures) {
					if(restrictToType==null || ss.getType() == restrictToType) {
						out.printf("%s\t%s\t%s\n", ss.getRange(), ss.getTypeChar(),
								ss.getType()==SubstructureType.PREDICTION?ss.getTemplatesStr():ss.getPdbCode()+ss.getChainCode());
					}
				}
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Writes a data file with mutation information.
	 * @param outFile the file to be written
	 * TODO: @param restrictToType if not null, only mutations in substructures of the given type will be written (e.g. SubstructureType.XRAY)
	 */
	public void writeMutationData(File outFile, SubstructureType restrictToType) {
		try {
			PrintWriter out = new PrintWriter(outFile);
			if (this.mutations != null) {
				
				for (Mutation m : this.mutations) {
					Substructure ss = this.getSubstructure(m.position);
					if(restrictToType==null || (ss != null && ss.getType() == restrictToType) && ss.pdb.hasCoordinates(ss.mapUniprotResser2Cif(m.position))) {
						String mutStr = m.getDnaMutStr() == null?m.toString():m.getDnaMutStr();
						out.printf("%d\t%s\t%s\t%s\n", m.position, mutStr, m.getComment(), m.getCommentDetails());
					}
				}
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}

	/**
	 * Writes a data file with mutation information.
	 * @param outFile the file to be written
	 * TODO: @param restrictToType if not null, only mutations in substructures of the given type will be written (e.g. SubstructureType.XRAY)
	 */
	public void writeSnpData(File outFile, SubstructureType restrictToType) {
		if(this.snps != null) {
			try {
				PrintWriter out = new PrintWriter(outFile);
				for (Mutation m : this.snps) {
					Substructure ss = this.getSubstructure(m.position);
					if(restrictToType==null || (ss != null && ss.getType() == restrictToType)) {
						out.printf("%d\t%s\n", m.position, m.toString());
					}
				}
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	public void renderRuler(String outFileName, int width) {
		try {
			String cmd = String.format(RULER_CMD, this.length,width,outFileName);
			System.out.println(cmd);
			Process p = Runtime.getRuntime().exec(cmd);
			
			BufferedReader pout = new BufferedReader(new InputStreamReader(p.getInputStream()));
			String line;
			while((line=pout.readLine())!=null) {
				System.err.println(line);
			}			
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}
	
	public void renderStructuralRegions(String outFileName, int width) {
		try {
			File tmpFile = File.createTempFile(this.geneName, ".txt");
			//tmpFile.deleteOnExit();
			PrintWriter out = new PrintWriter(tmpFile);
			out.println(this.geneName);
			out.println(this.length);
			if(this.subStructures != null) {
				for(Substructure ss:this.subStructures) {
					out.printf("%s %s %s\n", ss.getRange(), ss.getTypeChar(), "?");
				}
			}
			out.close();
			String cmd = String.format(RENDER_CMD, tmpFile.toString(),width,outFileName);
			System.out.println(cmd);
			Process p = Runtime.getRuntime().exec(cmd);
			
//			PrintWriter png = new PrintWriter(new File(outFileName));
			BufferedReader pout = new BufferedReader(new InputStreamReader(p.getInputStream()));
			String line;
			while((line=pout.readLine())!=null) {
//				png.println(line);
				System.err.println(line);
			}
//			png.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
		
	/**
	 * Loads the protein sequence for this gene name from COSMIC.
	 * offline: load from downloaded COSMIC data (for reproducibility)
	 * online:  load from current online COSMIC (for latest data)
	 * Sequence is stored in a member variable.
	 * @param online if true, work in online mode, otherwise in offline mode
	 * @return the downloaded sequence or null on error
	 */
	public String loadCosmicSequence(boolean online) {
		String newSeq = null;
		if(this.geneName == null) {
			System.err.println("Error loading COSMIC sequence. Gene name is null.");
		} else {
			File seqFile = null;	// either a local sequence file or a downloaded temp file
			if(online) {
				URL url = null;
				try {
					seqFile = File.createTempFile("cosmic", ".fa");
					seqFile.deleteOnExit();
					PrintWriter out = new PrintWriter(seqFile);
					url = new URL(String.format(COSMIC_SEQ_URL,this.geneName.substring(0, 1), this.geneName));
					BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
					String line;
					while ((line = in.readLine()) != null) {
						out.println(line);
					}
					in.close();
					out.close();
				} catch (MalformedURLException e) {
					e.printStackTrace();
				} catch (IOException e) {
					System.err.println("Error connecting to url " + url + ": " + e.getMessage());
				}
			} else {
				seqFile = new File(String.format(COSMIC_SEQ_PATH, this.geneName));
			}
			// now seqFile should exist
			if(seqFile != null) {
				Sequence seq = new Sequence();
				try {
					seq.readFromFastaFile(seqFile);
					// tests
					if(!seq.getName().equals(this.geneName)) {
						System.err.println("Warning: Expected gene name " + this.geneName + " in fasta header but found " + seq.getName());
					}				
					if(this.length > 0 && this.length != seq.getLength() - 1) {	// cosmic seqs contains stop codon
						System.err.println("Warning: Length of sequence downloaded from COSMIC (" + (seq.getLength() - 1) + ") does not match annotated length (" + this.length + ") for " + this.geneName);				
					}
					this.cosmicSeq = seq.getSeq().substring(0, seq.getSeq().length()-1);
					newSeq = seq.getSeq().substring(0, seq.getSeq().length()-1);
				} catch (FileFormatException e) {
					System.err.println("Error reading from sequence file " + seqFile + ": " + e.getMessage());
				} catch (IOException e) {
					System.err.println("Error reading from sequence file " + seqFile + ": " + e.getMessage());
				}
			}
		}
		return newSeq;		
	}
	
	/**
	 * Loads the cDNA sequence for this gene name from COSMIC.
	 * offline: load from downloaded COSMIC data (for reproducibility)
	 * online:  load from current online COSMIC (for latest data)
	 * Sequence is stored in a member variable.
	 * @param online if true, work in online mode, otherwise in offline mode
	 * @return the obtained sequence or null on error
	 */
	public String loadCosmicCdnaSequence(boolean online) {
		String newSeq = null;
		if(this.geneName == null) {
			System.err.println("Error loading COSMIC sequence. Gene name is null.");
		} else {
			File seqFile = null;	// either a local sequence file or a downloaded temp file
			if(online) {
				URL url = null;
				try {
					seqFile = File.createTempFile("cosmic", ".cdna.fa");
					seqFile.deleteOnExit();
					PrintWriter out = new PrintWriter(seqFile);
					url = new URL(String.format(COSMIC_CDNA_URL,this.geneName.substring(0, 1), this.geneName));
					BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
					String line;
					while ((line = in.readLine()) != null) {
						out.println(line);
					}
					in.close();
					out.close();
				} catch (MalformedURLException e) {
					e.printStackTrace();
				} catch (IOException e) {
					System.err.println("Error connecting to url " + url + ": " + e.getMessage());
				}
			} else {
				seqFile = new File(String.format(COSMIC_CDNA_PATH, this.geneName));
			}
			// now seqFile should exist
			if(seqFile != null) {
				Sequence seq = new Sequence();
				try {
					seq.readFromFastaFile(seqFile);
					// tests
					if(!seq.getName().equals(this.geneName)) {
						System.err.println("Warning: Expected gene name " + this.geneName + " in fasta header but found " + seq.getName());
					}				
					if(this.length > 0 && this.length != seq.getLength() - 1) {	// cosmic seqs contains stop codon
						System.err.println("Warning: Length of sequence downloaded from COSMIC (" + (seq.getLength() - 1) + ") does not match annotated length (" + this.length + ") for " + this.geneName);				
					}
					this.cosmicCdnaSeq = seq.getSeq().substring(0, seq.getSeq().length()-1);
					newSeq = seq.getSeq().substring(0, seq.getSeq().length()-1);
				} catch (FileFormatException e) {
					System.err.println("Error reading from sequence file " + seqFile + ": " + e.getMessage());
				} catch (IOException e) {
					System.err.println("Error reading from sequence file " + seqFile + ": " + e.getMessage());
				}
			}
		}
		return newSeq;
	}
	
	/**
	 * Loads the UniprotId for this gene from the COSMIC website and sets the member variable accordingly.
	 * @return the Uniprot id or null if either the gene was not found in COSMIC or no Uniprot ID was found.
	 */
	public String loadUniprotIdFromCosmic() {
		String uniProtId = null;		
		try {
			URL url = new URL(String.format(COSMIC_GENE_URL, this.getGeneName()));
			BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
			Pattern p = Pattern.compile(COSMIC_UNIPROT_REGEXP);
			String line;
			while((line = in.readLine()) != null) {
				Matcher m = p.matcher(line);
				if(m.find()) {
					uniProtId = m.group(1);	
					this.uniProtID = uniProtId;
					break;
				}
			}
			in.close();
		} catch(IOException e) {
			System.err.println("Gene " + this.getGeneName() + " not found in COSMIC");
		}
		return uniProtId;
	}
	
	/**
	 * Returns a set of residues for which (according to Uniprot) known structures are available.
	 * @param onlyXRay if true, consider only xray structures as known (otherwise also NMR)
	 * @param maxResolution maximum resolution for xray structures to be considered
	 * @param minLength minmum length for xray structures to be considered (including unobserved resdiues)
	 * @return a set of residue numbers for which structures are known (according to Uniprot)
	 */
	public TreeSet<Integer> getKnownRegionsFromUniprot(boolean onlyXRay, double maxResolution, int minLength) {
				
		IntervalSet knownRegions = new IntervalSet();
		if(this.uniProtID == null) {
			System.err.println("Error loading substructures from Uniprot. UniprotID is null.");
			return null;
		}
		// get Uniprot entry
		UniProtConnection upc = new UniProtConnection();
		UniProtEntry entry;
		try {
			entry = upc.getEntry(this.uniProtID);
			
			// get PDB cross references
			Collection<Pdb> refs = entry.getDatabaseCrossReferences(DatabaseType.PDB);
			for(Pdb ref:refs) {
				//System.out.println(ref);
				//String pdbCode = ref.getPdbAccessionNumber().getValue();
				String chainsStr = ref.getPdbChains().getValue();
				String[] chains = chainsStr.split(",");
				String[] fields = chains[0].trim().split("=");
				//String chain = fields[0].split("/")[0];
				if(fields.length < 2) {
					System.err.println("Error parsing xref from Uniprot for " + this.uniProtID + ". Sequence range not found.");
					continue;	// does this still work in try/catch block?
				}
				String[] begEnd = fields[1].split("-");
				int beg = Integer.parseInt(begEnd[0]);
				int end = Integer.parseInt(begEnd[1]);	    		
				String method = ref.getPdbMethod().getValue();
				String resStr = ref.getPdbResolution().getValue();
				double res = resStr.equals("-")?Double.NaN:Double.parseDouble(resStr.split(" ")[0]);
				//System.out.printf("%s:%s:%.2f:%s:%d-%d\n", pdbCode, method, res, chain, beg, end);
				if(method.equals("X-ray") || (!onlyXRay && method.equals("NMR"))) {
					if(!method.equals("X-ray") || res <= maxResolution) {
						if(end-beg+1 >= minLength) {
							knownRegions.add(new Interval(beg, end));
						}
					}
				}
			}
		} catch(NoMatchFoundException e) {
			System.err.println("Error loading substructures from Uniprot. UniprotID " + this.getUniprotId() + " not found.");
		}
		TreeSet<Integer> knownResidues = knownRegions.getIntegerSet();
		return knownResidues;
	}
	
	/**
	 * Loads the protein sequence for this gene from UniProt based on the Uniprot identifyer.
	 * offline: load from downloaded Uniprot data (for reproducibility)
	 * online:  load from current online UniProt (for latest data) 
	 * The sequence is stored as a member variable.
	 * @param online if true, work in online mode, otherwise in offline mode
	 * @return the obtained sequence or null on error
	 */
	public String loadUniprotSequence(boolean online) {
		String newSeq = null;
		if(this.uniProtID == null) {
			System.err.println("Could not load UniProt sequence. UniprotID for " + this.geneName + " is null");
		} else {
			File seqFile = null;	// either a local sequence file or a downloaded temp file
			if(online) {			
				URL url = null;
				try {
					seqFile = File.createTempFile("uniprot", ".fa");
					seqFile.deleteOnExit();
					PrintWriter out = new PrintWriter(seqFile);
					url = new URL(String.format(UNIPROT_SEQ_URL,this.uniProtID));
					BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
					String line;
					while ((line = in.readLine()) != null) {
						out.println(line);
					}
					in.close();
					out.close();
				} catch (MalformedURLException e) {
					e.printStackTrace();				
				} catch (IOException e) {
					System.err.println("Error connecting to url " + url + ": " + e.getMessage());
				}
			} else {
				seqFile = new File(String.format(UNIPROT_SEQ_PATH,this.uniProtID));
			}
			// now seqFile should exist
			if(seqFile != null && seqFile.canRead()) {
				Sequence seq = new Sequence();
				try {
					seq.readFromFastaFile(seqFile);
					// tests
					if(seq.getName().indexOf(this.uniProtID) < 0) {
						System.err.println("Warning: Expected UniProtID " + this.uniProtID + " in fasta header but found " + seq.getName());
					}
					if(seq.getName().indexOf(this.uniProtID) < 0) {
						System.err.println("Warning: Expected gene name " + this.geneName + " in fasta header but found " + seq.getName());
					}				
					if(this.length > 0 && this.length != seq.getLength()) {
						System.err.println("Warning: Length of sequence downloaded from UniProt (" + seq.getLength() + ") does not match annotated length (" + this.length + ")");				
					}
					this.uniprotSeq = seq.getSeq();
					newSeq = seq.getSeq();
				} catch (FileFormatException e) {
					System.err.println("Error reading from sequence file " + seqFile + ": " + e.getMessage());
				} catch (IOException e) {
					System.err.println("Error reading from sequence file " + seqFile + ": " + e.getMessage());
				}
			}
		}
		return newSeq;
	}
	
	/**
	 * Loads an alignment from Cosmic to Uniprot sequence.
	 */
	public void loadCosmic2UniprotAlignment() {
		if(this.uniprotSeq == null) {
			System.err.println("Error creating cosmic2Uniprot alignment. Uniprot sequence not loaded.");
			return;
		}
		if(this.cosmicSeq == null) {
			System.err.println("Error creating cosmic2Uniprot alignment. Cosmic sequence not loaded.");
			return;
		}
		try {
			this.cosmic2Uniprot = new PairwiseSequenceAlignment(this.cosmicSeq, this.uniprotSeq, "Cosmic", "Uniprot");
		} catch (PairwiseSequenceAlignmentException e) {
			System.err.println("Could not load Cosmic to Uniprot alignment for " + this.geneName);
		}
	}
	
	/**
	 * Based on the preloaded alignment between cosmic sequence and uniprot sequence, return the uniprot position
	 * given a cosmic position. Positions are counted from 1 to length.
	 * @param cosmicPos the position in the cosmic sequence
	 * @return the corresponding position in the uniprot sequence or -1 if the given position maps to a gap
	 */
	public int mapCosmicPos2Uniprot(int cosmicPos) {
		if(this.cosmic2Uniprot == null) this.loadCosmic2UniprotAlignment();
		return this.cosmic2Uniprot.getMapping1To2()[cosmicPos-1]+1;
	}

	/**
	 * Based on the preloaded alignment between cosmic sequence and uniprot sequence, return the cosmic position
	 * given a uniprot position. Positions are counted from 1 to length.
	 * @param uniprotPos the position in the uniprot sequence
	 * @return the corresponding position in the cosmic sequence or -1 if the given position maps to a gap
	 */
	public int mapUniprotPos2Cosmic(int uniprotPos) {
		if(this.cosmic2Uniprot == null) this.loadCosmic2UniprotAlignment();
		return this.cosmic2Uniprot.getMapping2To1()[uniprotPos-1]+1;
	}
	
	/**
	 * Aligns Uniprot vs. Cosmic sequence for this gene.
	 * @return the %identity
	 */
	public double compareUniprotVsCosmicSequence() {	
		return 0;
	}
	
	public void visualizeMutations(File outFile) {
	}
	
	public void writePseWithAllMutations(File outDir, String baseName, File pdbFileDir) throws IOException {
		for(Substructure ss:this.getSubstructures()) {
			ss.writePseWithAllMutations(outDir, baseName, pdbFileDir, this.getMutations());
		}
	}
		
	public void visualizeMutationsAndFeatures(File outFile) {
		
	}
	
	/*--------------------------------- main --------------------------------*/
	
	public static void highlightResidue(PrintWriter out, int pos, String color) {
		out.println("select resi " + pos);
		out.println("show sticks, sele");
		out.println("color " + color + ", sele");
		out.println("show spheres, sele and name ca");	
	}
	
	public static void showMLH1(MySQLConnection conn) {
		Gene g = new Gene("MLH1","P40692");
		g.loadMutations(conn, true);	// only missense
		g.loadSNPs(conn);
				
		String pdbFileName = "/project/PyMol/pdbs_final/MLH1_1-340_pred.pdb";
		File scriptFile = new File("/project/StruPPi/henning/projects/mutanom/scripts/test/MLH1.pml");
		int start = 1;
		int end = 340;
		int offset = 0;
		try {
			PrintWriter out = new PrintWriter(scriptFile);
			out.println("load " + pdbFileName);
			out.println("as cartoon");
			out.println("orient");
			out.println("set sphere_scale, 1.2");
			out.println("set side_chain_helper, 1");
			out.println("set sphere_transparency, 0.1");
			for(Mutation m:g.snps) {
				if(m.position >= start && m.position <= end) {
					highlightResidue(out, m.position+offset, "yellow");
				}
			}	
			for(Mutation m:g.mutations) {
				System.out.println(m);
				if(m.position >= start && m.position <= end) {
					System.out.print(".");
					highlightResidue(out, m.position+offset, "red");
				}
			}
			out.close();
			System.out.println();
		} catch (IOException e) {
			e.printStackTrace();
		} 
	}
	
	public static void showERBB2(MySQLConnection conn) {
		Gene g = new Gene("ERBB2","P04626");
		g.loadMutations(conn, true);
		g.loadSNPs(conn);
		
		String pdbFileName = "/project/PyMol/pdbs_final/ERBB2_720-987_pred.pdb";
		File scriptFile = new File("/project/StruPPi/henning/projects/mutanom/scripts/test/ERBB2.pml");
		int start = 720;
		int end = 987;
		int offset = 0;
		try {
			PrintWriter out = new PrintWriter(scriptFile);
			out.println("load " + pdbFileName);
			out.println("as cartoon");
			out.println("orient");
			out.println("set sphere_scale, 1.2");
			out.println("set sphere_transparency, 0.1");
			for(Mutation m:g.mutations) {
				if(m.position >= start && m.position <= end) {
					highlightResidue(out, m.position+offset, "red");
					System.out.print(".");
				}
			}
			for(Mutation m:g.snps) {
				if(m.position >= start && m.position <= end) {
					highlightResidue(out, m.position+offset, "yellow");
				}
			}			
			System.out.println();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		} 

	}
	
	public static void main(String[] args) {
		MySQLConnection conn;
		try {
			conn = new MySQLConnection();
			showMLH1(conn);
			showERBB2(conn);
			System.out.println("Done.");
		} catch (SQLException e) {
			e.printStackTrace();
		}

	}

	
}
