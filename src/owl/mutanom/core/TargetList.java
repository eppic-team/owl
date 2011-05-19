
package owl.mutanom.core;

//import java.awt.BasicStroke;
//import java.awt.Graphics2D;
//import java.awt.Rectangle;
import java.io.*;
import java.sql.SQLException;
import java.util.*;

import owl.core.connections.NoMatchFoundException;
import owl.core.connections.PrositeHit;
import owl.core.connections.PrositeScanner;
import owl.core.connections.UniProtConnection;
import owl.core.features.Feature;
import owl.core.features.FeatureType;
import owl.core.features.GeneralFeature;
import owl.core.features.PhosphoSitePlusFeature;
import owl.core.features.ProteinModificationType;
import owl.core.features.StructuralDomainType;
import owl.core.features.UniprotFeature;
import owl.core.sequence.Sequence;
import owl.core.sequence.alignment.PairwiseSequenceAlignment;
import owl.core.sequence.alignment.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;
import owl.core.structure.AminoAcid;
import owl.core.structure.AaResidue;
import owl.core.structure.graphs.RIGNode;
import owl.core.util.Goodies;
import owl.core.util.Interval;
import owl.core.util.MySQLConnection;
import owl.core.util.Pair;
import owl.mutanom.FoldXRunner;
import owl.mutanom.core.Mutation.MutType;
import owl.mutanom.core.Substructure.SubstructureType;
import owl.mutanom.output.PyMolScriptMaker;
import owl.mutanom.output.SvgGenerator;
import owl.mutanom.output.PyMolScriptMaker.Color;

//import org.apache.batik.dom.GenericDOMImplementation;
//import org.apache.batik.svggen.SVGGraphics2D;
//import org.w3c.dom.DOMImplementation;
//import org.w3c.dom.Document;




import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
//import uk.ac.ebi.kraken.interfaces.uniprot.features.ActSiteFeature;
//import uk.ac.ebi.kraken.interfaces.uniprot.features.BindingFeature;
//import uk.ac.ebi.kraken.interfaces.uniprot.features.Feature;
//import uk.ac.ebi.kraken.interfaces.uniprot.features.ModResFeature;
//import uk.ac.ebi.kraken.interfaces.uniprot.features.NpBindFeature;
//import uk.ac.ebi.kraken.interfaces.uniprot.features.RegionFeature;

/**
 * A list of targets genes to be analyzed. Provides methods to load and analyze a target list.
 * @author stehr
 *
 */
public class TargetList {
	/*------------------------------ constants ------------------------------*/
	public static final String GENE_FILE_NAME = 		"genes.txt"; // for supplementary material
	
	public static final String DATABASE = 				"mutanom3"; 
	
	public static final String TARGETS_TABLE = 			DATABASE + ".targets";	
	public static final String GENE2SP_TABLE = 			DATABASE + ".gene2sp";
	public static final String GENE2NM_TABLE = 			DATABASE + ".gene2nm";
	public static final String GENE2NP_TABLE = 			DATABASE + ".gene2np";
	public static final String PRED_TARGETS_TABLE =		DATABASE + ".targets_predicted";

	
	public static final String QUERY_LOAD_TARGETS = 	"SELECT * FROM " + TARGETS_TABLE;
	public static final String QUERY_GENE2SP = 			"SELECT sp_id FROM " + GENE2SP_TABLE + " WHERE gene_name='%s'";
	public static final String QUERY_GENE2NM = 			"SELECT nm FROM " + GENE2NM_TABLE + " WHERE gene_name='%s'";
	public static final String QUERY_GENE2NP = 			"SELECT np FROM " + GENE2NP_TABLE + " WHERE gene_name='%s'";
	public static final String QUERY_LOAD_PREDICTED = 	"SELECT * FROM " + PRED_TARGETS_TABLE;
	
	public static final boolean LOAD_ONLY_MISSENSE = true;	// if true, only missense mutations are loaded (excluding silent mutations)
	
	/*--------------------------- member variables --------------------------*/
	Collection<Gene> targets;	// collection of target genes
	
	/*----------------------------- constructors ----------------------------*/
	public TargetList() {
		targets = new LinkedList<Gene>(); 
	}
	
	/*---------------------------- public methods ---------------------------*/
	public void addTarget(Gene gene) {
		this.targets.add(gene);
	}
	
	/**
	 * Removes the gene with the given name from this target list.
	 * @param geneName
	 * @return the gene which has been removed or null if removal was not successfull
	 */
	public Gene removeTarget(String geneName) {
		Gene toBeRemoved = null;
		for(Gene g:targets) {
			if(g.getGeneName().equals(geneName)) {
				toBeRemoved = g;
			}
		}
		if(toBeRemoved != null) targets.remove(toBeRemoved);
		return toBeRemoved;
	}
	
	/**
	 * Returns the subset of targets which are annotated as Oncogenes.
	 * Relies on proper annotation in the database.
	 * @return A targetlist containing only the oncogenes
	 */
	public TargetList getOncogenes(MySQLConnection conn) {
		TargetList newList = new TargetList();
		for(Gene g:this.targets) {
			if(g.isOncogene(conn)) {
				newList.addTarget(g);
			}
		}
		return newList;
	}
	
	/**
	 * Returns the subset of targets which are annotated as Tumor Suppressors.
	 * Relies on proper annotation in the database.
	 * @return A TargetList containing only the tumor suppressors
	 */	
	public TargetList getTumorSuppressors(MySQLConnection conn) {
		TargetList newList = new TargetList();
		for(Gene g:this.targets) {
			if(g.isTumorSuppressor(conn)) {
				newList.addTarget(g);
			}
		}
		return newList;
	}
	
	/**
	 * Returns the subset of targets which are neither annotated as Oncogene nor as Tumor Suppressor.
	 * Relies on proper annotation in the database.
	 * @return A TargetList containing only the gene neither annotated as oncogenes nor as tumor suppressors
	 */	
	public TargetList getUnAnnotated(MySQLConnection conn) {
		TargetList newList = new TargetList();
		for(Gene g:this.targets) {
			if(!g.isTumorSuppressor(conn) && !g.isOncogene(conn)) {
				newList.addTarget(g);
			}
		}
		return newList;
	}	
	
	/**
	 * Returns the subset of targets which are annotated with the given subcellular localization.
	 * Relies on proper annotation in database table 'mutanom.gene'.
	 * @return A targetlist containing only the specified genes
	 */
	public TargetList getLocalizedList(MySQLConnection conn, Gene.SubcellularLocalization loc) {
		TargetList newList = new TargetList();
		String table = "mutanom.gene";
		String locStr = "";
		switch(loc) {
			case NUCLEUS: locStr = "N"; break;
			case CYTOPLASM: locStr = "C"; break;
			case MEMBRANE: locStr = "M"; break;
			case TRANSMEMBRANE: locStr = "T"; break;
		}
		String sql = "SELECT COUNT(*) FROM " + table + " WHERE name='%s' AND loc='%s'";
		for(Gene g:this.targets) {
			int rightLocation = conn.getIntFromDb(String.format(sql,g.getGeneName(),locStr));
			if(rightLocation != 0) newList.addTarget(g);
		}
		return newList;
	}	
	
	/**
	 * Returns a new target list of the same size where the members are sampled at random with replacements
	 * from the original list. A population of such samples is used for bootstrapping analyses.
	 * TODO: Is it a problem that there may be multiple copies of the same target in the sample?
	 */
	public TargetList getRandomSample() {
		int n = this.targets.size();
		Gene[] targetArray = new Gene[n];
		targetArray = this.targets.toArray(targetArray);
		TargetList newList = new TargetList();
		Random r = new Random();
		for (int i = 0; i < this.targets.size(); i++) {
			int pick = r.nextInt(this.targets.size());
			newList.addTarget(targetArray[pick]);
		}
		return newList;
	}
	
	/**
	 * Returns a new target list which excludes the given target from the list.
	 */
	public TargetList getOneLess(Gene toExclude) {
		TargetList newList = new TargetList();
		for(Gene g:this.targets) {
			if(g!=toExclude) {
				newList.addTarget(g);
			}
		}
		return newList;
	}
	
	/**
	 * Removes the first substructure of any gene within this list matching the given pdb code and chain code
	 * TODO: not tested!
	 * @param pdbAndChainCode
	 * @return
	 */
	public Substructure removeSubstructure(String pdbCode, String chainCode) {
		Substructure toRemove = null;
		for(Gene g:this.targets) {
			for(Substructure ss: g.subStructures) {
				if((ss.getPdbCode().equals(pdbCode) && ss.getChainCode().equals(chainCode))) {
					toRemove = ss;
					break;
				}
			}
			if(toRemove != null) {
				if(g.subStructures.remove(toRemove)) {
					return toRemove;
				} else {
					return null;
				}
			}
		}
		return null;
	}
	
	public Collection<Gene> getTargets() {
		return targets;
	}
	
	public void report() {
		int c = 1;
		for(Gene g:this.targets) {
			System.out.printf("%02d:\t", c++);
			System.out.println(g);
		}
	}
	
	/**
	 * Returns the number of loaded mutations summed over all genes.
	 * @return number of mutations found
	 */
	public int getNumberOfMutations() {
		int num = 0;
		for(Gene g:this.targets) {
			num += g.getMutations().size();
		}
		return num;
	}
	
	/**
	 * Returns the number of mutated positions summed over all genes.
	 * @returnnumber of mutated positions found
	 */
	public int getNumberOfMutatedPositions() {
		int num = 0;
		for(Gene g:this.targets) {
			TreeSet<Integer> positions = new TreeSet<Integer>();
			for(Mutation m:g.getMutations()) {
				positions.add(m.position);
			}
			num += positions.size();
		}
		return num;		
	}
	
	/**
	 * Prints statistics about this target list (# genes, #substructures, #mutations)
	 */
	public void printStats(MySQLConnection conn) {
		int numGenes = this.getTargets().size();
		int numSubStructures = 0;
		int numMutations = 0;
		int numMutationsWithStructure = 0;
		int numSNPs = 0;
		int numSNPsWithStructure = 0;
		for(Gene g:this.getTargets()) {
			int str = 0;
			int mut = 0;
			int snp = 0;
			int strSnp = 0;
			int funcSites = 0;	// all features except domains 
			str = g.getSubstructures()==null?0:g.getSubstructures().size();
			mut = g.getMutations()==null?0:g.getMutations().size();
			snp = g.getSNPs()==null?0:g.getSNPs().size();
			numSubStructures += str;
			numMutations += mut;
			numSNPs += snp;
			
			// structural mutations
			int strMut = 0;			// structural mutations
			//int strMutPos = 0; // distinct mutated positions (in structural regions)
			int strMutMatch = 0; 	// structural mutations matching the amino acid of the PDB sequence
			//int strMutMatchPos = 0; // distinct mutated positions (in structural regions) where amino acid matches
			if(mut > 0) {
				TreeSet<Integer> mutatedPositions = new TreeSet<Integer>();
				TreeSet<Integer> mutatedMatchPositions = new TreeSet<Integer>();
				for(Mutation m:g.getMutations()) {
					if(m.type == MutType.MISSENSE && g.getSubstructure(m.position) != null) {
						Substructure ss = g.getSubstructure(m.position);
						if(ss.pdb.containsStdAaResidue(ss.mapUniprotResser2Cif(m.position))) {
							strMut++;
							//int cosmicPos = g.mapUniprotPos2Cosmic(m.position);
							//if(m.getWtAA() == AminoAcid.getByOneLetterCode(g.getCosmicSeq().charAt(cosmicPos-1))) {
							//if(m.getWtAA() == AminoAcid.getByOneLetterCode(g.getUniprotSeq().charAt(m.position-1))) {
							if(m.getWtAA() == ((AaResidue)ss.pdb.getResidue(ss.mapUniprotResser2Cif(m.position))).getAaType()) {
								strMutMatch++;
								mutatedMatchPositions.add(m.position);
							}
							mutatedPositions.add(m.position);
						}
					}
				}
				//strMutPos = mutatedPositions.size();
				//strMutMatchPos = mutatedMatchPositions.size();
			}
			numMutationsWithStructure += strMut;
			
			// structural SNPs
			if(snp > 0) {
				for(Mutation m:g.getSNPs()) {
					if(g.getSubstructure(m.position) != null) {
						Substructure ss = g.getSubstructure(m.position);
						if(ss.pdb.containsStdAaResidue(ss.mapUniprotResser2Cif(m.position))) {
							strSnp++;
						}
					}
				}
			}
			numSNPsWithStructure += strSnp;
			
			// features
			for(Feature f: g.getFeatures()) {
				if(f.getType() != FeatureType.SDOMAIN) funcSites ++;
			}
			
			// sequence coverage
			int len = g.getLength();	// total (cosmic) length
			int res = 0;				// total residues covered by structural regions
			int strRes = 0;				// observed residues in structural regions
			//double cov = 0.0;			// total sequence coverage in %
			double strCov = 0.0;		// sequence coverage by observed residues in %			
			if(str > 0) {
				for(Substructure ss:g.getSubstructures()) {
					res += ss.getEndPos()-ss.getBegPos()+1;
					strRes += ss.getPdb().getStdAaObsLength();
				}
				//cov = 100.0 * res / len;
				strCov = 100.0 * strRes / len;
			}
			
			// oncogene/tumor suppressor
			String type = g.getOncSupAnnotation(conn);
			if(type == null) type="-";
			
			// pdbs
			String pdbStr = "";
			for(Substructure ss:g.getSubstructures()) {
				pdbStr += ", " + ss.getPdbCode().toUpperCase()+ss.getChainCode().toUpperCase()+" (" + ss.getRange() + ")";
			}
			pdbStr = pdbStr.substring(2);
			
			// print stats per gene
			//System.out.printf("%s\tSubStr: %d\tMut: %3d\tMutPos: %3d\tMutMatch: %3d\tMutMatchPos: %3d\tSnp: %3d\tFunc:%3d\tLen: %4d\tRes: %4d\tStrRes:%3d\tCov:%4.0f\tStrCov:%4.0f\tOncSup: %s\n", g.getGeneName(), str, strMut, strMutPos, strMutMatch, strMutMatchPos, strSnp, funcSites, len, res, strRes, cov, strCov, type);
			//System.out.printf("%s,SubStr:,%d,Mut: ,%3d,MutPos: ,%3d,MutMatch: ,%3d,MutMatchPos: ,%3d,Snp: ,%3d,Func:,%3d,Len: ,%4d,Res: ,%4d,StrRes:,%3d,Cov:,%4.0f,StrCov:,%4.0f,OncSup: ,%s\n", g.getGeneName(), str, strMut, strMutPos, strMutMatch, strMutMatchPos, strSnp, funcSites, len, res, strRes, cov, strCov, type);
			System.out.printf("%s;%d;%d;%d;%d;%6.2f;%s\n", g.getGeneName(), len, strMut, strSnp, funcSites, strCov, pdbStr);

		}
		System.out.println("Number of target genes: " + numGenes);
		System.out.println("Number of substructures: " + numSubStructures);
		System.out.println("Number of mutations: " + numMutations);
		System.out.println("Number of missense mutations with structure: " + numMutationsWithStructure);
		System.out.println("Number of SNPs: " + numSNPs);
		System.out.println("Number of SNPs with structure: " + numSNPsWithStructure);		
	
		// specialized stats
	
	}
	
	public void printUniProtIds() {
		for(Gene g: this.targets) {
			System.out.println(g.geneName + "\t" + g.uniProtID);
		}
	}
	
	/**
	 * Loads the substructure information for all target genes in the list
	 * @throws SQLException 
	 */
	public void loadSubstructures(MySQLConnection conn, boolean onlyPredicted) throws SQLException {
		for(Gene g:this.targets) {
			//System.out.print(g.getGeneName() + "\t");
			g.loadSubstructures(conn, onlyPredicted);
			//System.out.println();
		}
	}
	
	/**
	 * For previously loaded substructures, load the corresponding PDB objects. If substructure
	 * is a known structure, load from Pdbase, otherwise load from a PDB files with default name
	 * from the given directory. Also initializes the alignments from Uniprot to PDB sequence
	 */
	public void loadPdbsAndAlignments(File pdbDir) {
		for(Gene g:this.targets) {
			for(Substructure ss:g.getSubstructures()) {
				if(ss.getType() == Substructure.SubstructureType.PREDICTION) {
					File pdbFile = new File(pdbDir, ss.getDefaultPdbFileName(g.geneName));
					ss.loadPdbFromFile(pdbFile);					
				} else {
					ss.loadPdbFromDb();
				}
				if(g.getUniprotSeq() == null) {
					System.err.println("Cannot initialize alignment, Uniprot sequence for gene " + g.geneName + " not loaded.");
				} else {
					ss.initializeAlignment(g.getUniprotSeq(), g.getUniprotId());
				}
			}
		}		
	}
	
	/**
	 * Load RIGs and sequence/structure features such as catalytic sites for all substructures.
	 * @param useManFuncAnno TODO
	 * @return the number of features loaded
	 */
	public int loadFeaturesAndGraphsAndDomains(File phosphoSiteHtmlDir, String contactType, double cutoff, StructuralDomainType domType, boolean useManFuncAnno) {
		int featuresLoaded = 0;
		int graphsLoaded = 0;
		int domainsLoaded = 0;
		System.out.println("Loading features...");
		for(Gene g:this.targets) {
			System.out.println(g.geneName);
			// Uniprot
			featuresLoaded += g.loadUniprotFeatures();
			// PhoshoSitePlus
			File htmlFile = new File(phosphoSiteHtmlDir, g.getGeneName() + ".html");
			featuresLoaded += g.loadPhosphoSitePlusFeatures(htmlFile);
			// ProSite
			//featuresLoaded += g.loadPrositeFeatures();
			// CSA
			featuresLoaded += g.loadCSAFeatures();
			// Manual annotations
			featuresLoaded += g.loadManualFeatures();
			// Graphs (for proximity calculation)
			for(Substructure ss:g.getSubstructures()) {
				//System.out.println(ss.getPdbCode() + ss.getChainCode());
				ss.loadGraph(contactType, cutoff);
				if(ss.isGraphLoaded()) graphsLoaded++;
				//featuresLoaded += ss.loadFeatures();
			}
			// Domains (for cluster analysis)
			domainsLoaded += g.loadStructuralDomains(domType);
		}
		System.out.println();
		System.out.println(graphsLoaded + " graphs loaded.");
		System.out.println(featuresLoaded + " features loaded.");
		System.out.println(domainsLoaded + " domains loaded.");
		return featuresLoaded;
	}
	
	/**
	 * Prints statistics about the types of functional site annotations
	 */
	public void printFunctionalSiteStats() {
		int act = 0;
		int atp = 0;
		int gtp = 0;
		int pho = 0;
		int ubq = 0;
		int ned = 0;
		int sum = 0;
		int met = 0;
		int aty = 0;
		int gly = 0;
		int nit = 0;
		int dom = 0;
		int mod = 0;
		int unk = 0;
		int total = 0;
		for(Gene g: this.getTargets()) {
			for(Feature f: g.getFeatures()) {
				total++;
				switch(f.getType()) {
				case GENERAL:
					if(f.getDescription().startsWith("UBQ_RES")) ubq++; else
					if(f.getDescription().startsWith("PHO_RES")) pho++; else
					if(f.getDescription().startsWith("ACT_SITE")) act++; else
					if(f.getDescription().startsWith("ATP_BIND")) atp++; else
					if(f.getDescription().startsWith("GTP_BIND")) gtp++; else
					if(f.getDescription().startsWith("ATY_RES")) aty++; else
					if(f.getDescription().startsWith("MOD_RES")) mod++; else	
					{System.out.println(g+" "+f); unk++;}
					break;
				case SDOMAIN: dom++; break;
				case UNIPROT:
					UniprotFeature f2 = (UniprotFeature) f;
					String upType = f2.getUniprotTypeName();
					if(upType.startsWith("ACT_SITE")) act++; else
					if(upType.startsWith("CARBOHYD")) gly++; else	
					if(upType.startsWith("NP_BIND")) {
						if(f2.getDescription().startsWith("ATP")) atp++; else
						if(f2.getDescription().startsWith("GTP")) gtp++; else
						{System.out.println(f); unk++;}
					} else
					if(upType.startsWith("MOD_RES")) {
						if(f2.getDescription().indexOf("Phospho") >= 0) pho++; else
						if(f2.getDescription().indexOf("methyl") >= 0) met++; else
						if(f2.getDescription().indexOf("acetyl") >= 0) aty++; else
						if(f2.getDescription().indexOf("nitroso") >= 0) {nit++;} else
						{System.out.println(f); unk++;}
					} else
					{System.out.println(f); unk++;}
					break;
				case CSA: act++; break;
				case PHOSPHOSITE: 
					PhosphoSitePlusFeature f3 = (PhosphoSitePlusFeature) f;
					switch(f3.getModType()) {
					case PHOSPHORYLATION: pho++; break;
					case UBIQUITINATION: ubq++; break;
					case ACETYLATION: aty++; break;
					case METHYLATION: met++; break;
					case NEDDYLATYION: {ned++;}break;
					case SUMOYLATION: {sum++;}break;
					default: System.out.println(f); unk++;
					}
					break;
				default: System.out.println(f); unk++;
				}
			}
		}
		System.out.println();
		System.out.println("Feature types:");
		System.out.printf("act=%d atp=%d gtp=%d pho=%d met=%d aty=%d gly=%d nit=%d ubq=%d ned=%d sum=%d dom=%d unk=%d\n", act, atp, gtp, pho, met, aty, gly, nit, ubq, ned, sum, dom, unk);
		System.out.printf("sum=%d total=%d\n", act+ atp+ gtp+ pho+ met+ aty+ gly+ nit+ ubq+ ned+ +sum +unk, total-dom);
	}
	
	/**
	 * Check for each observed residue whether it coincides with or is proximal to a functional site
	 * as annotated by Uniprot and other databases.
	 * @param silent if silent=false, print stats of hits and proximal hits
	 * @param conn
	 * @throws SQLException 
	 */
	public double countFeaturesAllResidues(boolean silent, PrintStream out) {

		int actHits = 0, actProx = 0;	// Uniprot+CSA+Manual
		int modHits = 0, modProx = 0;	// Uniprot+PhosphoSite+Manual (all except phosphorilation and glycosilation)
		int phoHits = 0, phoProx = 0;	// Uniprot+PhosphoSite+Manual (phosphorilation sites)
		int ubqHits = 0, ubqProx = 0;	// PhosphoSite+Manual (ubiquitination sites)		
		int atpHits = 0, atpProx = 0;	// Uniprot+Prosite+Manual (atp binding sites)
		int gtpHits = 0, gtpProx = 0;	// Uniprot+Prosite+Manual (gtp binding sites)
		int glycHits = 0, glycProx = 0;	// Uniprot+Manual
		int dnaHits = 0, dnaProx = 0;   // Uniprot+Manual
		int numMut = 0;
		
		for(Gene g:this.getTargets()) {
			if(!g.areSubstructuresLoaded()) {
				System.err.println("Skipping gene " + g.getGeneName() + ": substructures not loaded.");
			}
			for(Substructure ss:g.getSubstructures()) {
				if(ss.getType() == SubstructureType.PREDICTION) {
					System.out.println("Skipping predicted structure " + g.getGeneName() + "_" + ss.getRange());
					continue;
				}
				if(!ss.isPdbLoaded()) {
					System.err.println("Skipping " + g.getGeneName() + " " + ss.getRange() + " " + ss.getPdbCode()+ss.getChainCode() + ": Structure not loaded.");
				} else {
					if(!silent) System.out.print(".");				
					//String pdbCode = ss.getPdbCode();
					String cifChain = ss.pdb.getChainCode();
					String pdbChain = ss.pdb.getPdbChainCode(); // = ss.getChainCode();
					if(!pdbChain.equals(cifChain)) {
						System.err.println("Warning: pdb chain code " + pdbChain + " != cif chain code " + cifChain);
					}
					for(int p: ss.pdb.getAllStdAaResSerials()) {	// only observed residues
						//int pdbPos = Integer.parseInt(ss.pdb.getResidue(p).getPdbSerial());	// cif -> pdb serial
						int cifPos = p;
						int uniprotPos = ss.mapCifResser2Uniprot(p);
						//AminoAcid wtAa = ss.pdb.getResidue(p).getAaType();		
						// check for functional sites
						String fsHit = null;
						String fsNear = null;
						// Remember observed features to prevent overcounting of the same feature
						HashSet<Feature> alreadyObserved = new HashSet<Feature>();
						// check whether this position matches a feature
						for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.UNIPROT, uniprotPos)) {
							alreadyObserved.add(f);
							UniprotFeature uf = (UniprotFeature) f;
							if(uf.getUniprotTypeName().equals("ACT_SITE")) {actHits++; fsHit="ACT_SITE"; break;}
							if(uf.getUniprotTypeName().equals("NP_BIND")) {
								if(uf.getDescription().startsWith("ATP")) {atpHits++; fsHit="NP_BIND"; break;}
								else if(uf.getDescription().startsWith("GTP")) {gtpHits++; fsHit="NP_BIND"; break;}
								else System.err.println("Unknown NP_BIND feature: " + uf.getDescription()); break;
							}
							if(uf.getUniprotTypeName().equals("MOD_RES")) {
								if(uf.getDescription().startsWith("Phospho")) phoHits++;
								else modHits++; 
								fsHit="MOD_RES"; break;
							}
							if(uf.getUniprotTypeName().equals("CARBOHYD")) {glycHits++;	fsHit="CARBOHYD"; break;}
							if(uf.getUniprotTypeName().equals("DNA_BIND")) {dnaHits++;	fsHit="DNA_BIND"; break;}
						}
						if(fsHit == null) {
							for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.CSA, uniprotPos)) {
								alreadyObserved.add(f);
								actHits++;
								fsHit="ACT_SITE";
								break;
							}
						}
						if(fsHit == null) {
							for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.GENERAL, uniprotPos)) {
								alreadyObserved.add(f);
								GeneralFeature gf = (GeneralFeature) f;
								if(gf.getDescription().equals("ACT_SITE")) {actHits++; fsHit="ACT_SITE"; break;}
								if(gf.getDescription().equals("ATP_BIND")) {atpHits++; fsHit="ATP_BIND"; break;}
								if(gf.getDescription().equals("GTP_BIND")) {gtpHits++; fsHit="GTP_BIND"; break;}
								if(gf.getDescription().equals("ATY_RES")) {modHits++; fsHit="MOD_RES"; break;}								
								if(gf.getDescription().equals("MOD_RES")) {modHits++; fsHit="MOD_RES"; break;}
								if(gf.getDescription().equals("UBQ_RES")) {ubqHits++; fsHit="UBQ_RES"; break;}
								if(gf.getDescription().equals("PHO_RES")) {phoHits++; fsHit="PHO_RES"; break;}
								if(gf.getDescription().equals("CARBOHYD")) {glycHits++;	fsHit="CARBOHYD"; break;}
								if(gf.getDescription().equals("DNA_BIND")) {dnaHits++;	fsHit="DNA_BIND"; break;}
							}
						}
						if(fsHit == null) {
							for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.PHOSPHOSITE, uniprotPos)) {
								alreadyObserved.add(f);
								fsHit="MOD_RES";
								if(((PhosphoSitePlusFeature)f).getModType() == ProteinModificationType.PHOSPHORYLATION) phoHits++;
								else if(((PhosphoSitePlusFeature)f).getModType() == ProteinModificationType.UBIQUITINATION) ubqHits++;
								else modHits++;
								break;
							}
						}
						if(ss.graph == null) System.err.println("Error: graph not loaded for " + g.geneName + " " + ss.getRange());
						RIGNode mutNode = ss.graph.getNodeFromSerial(cifPos);
						if(mutNode != null) {	// residue observed
							numMut++;
							if(fsHit == null) {	// only check proximal hit if no direct hit found
								for(RIGNode n:ss.graph.getNeighbors(mutNode)) {
									int pdbNbPos = n.getResidueSerial();
									int uniNbPos = ss.mapCifResser2Uniprot(pdbNbPos);
									if(fsNear == null) { // only check if no prox hit found yet
										// check whether this position matches a feature
										for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.UNIPROT, uniNbPos)) {
											if(!alreadyObserved.contains(f)) {	// avoid double counting of the same feature
												alreadyObserved.add(f);
												UniprotFeature uf = (UniprotFeature) f;
												if(uf.getUniprotTypeName().equals("ACT_SITE")) {actProx++; fsNear="ACT_SITE";break;}
												if(uf.getUniprotTypeName().equals("NP_BIND")) {
													if(uf.getDescription().startsWith("ATP")) {atpProx++; fsNear="NP_BIND"; break;}
													else if(uf.getDescription().startsWith("GTP")) {gtpProx++; fsNear="NP_BIND"; break;}
													else System.err.println("Unknown NP_BIND feature: " + uf.getDescription()); break;												
												}
												if(uf.getUniprotTypeName().equals("MOD_RES")) {
													if(uf.getDescription().startsWith("Phospho")) phoProx++;
													else modProx++; 
													fsNear="MOD_RES"; break;
												}
												if(uf.getUniprotTypeName().equals("CARBOHYD")) {glycProx++; fsNear="CARBOHYD"; break;}
												if(uf.getUniprotTypeName().equals("DNA_BIND")) {dnaProx++; fsNear="DNA_BIND"; break;}
											}
										}
									}
									if(fsNear == null) {
										for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.CSA, uniNbPos)) {
											if(!alreadyObserved.contains(f)) {	// avoid double counting of the same feature
												actProx++;
												fsNear="ACT_SITE";
												break;
											}
										}
									}
									if(fsNear == null) {
										for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.GENERAL, uniNbPos)) {
											if(!alreadyObserved.contains(f)) {	// avoid double counting of the same feature
												alreadyObserved.add(f);
												GeneralFeature gf = (GeneralFeature) f;
												if(gf.getDescription().equals("ACT_SITE")) {actProx++; fsNear="ACT_SITE";break;}
												if(gf.getDescription().equals("ATP_BIND")) {atpProx++; fsNear="ATP_BIND"; break;}
												if(gf.getDescription().equals("GTP_BIND")) {gtpProx++; fsNear="GTP_BIND"; break;}
												if(gf.getDescription().equals("MOD_RES")) {modProx++; fsNear="MOD_RES"; break; }
												if(gf.getDescription().equals("ATY_RES")) {modProx++; fsNear="MOD_RES"; break; }
												if(gf.getDescription().equals("UBQ_RES")) {ubqProx++; fsNear="UBQ_RES"; break; }
												if(gf.getDescription().equals("PHO_RES")) {phoProx++; fsNear="PHO_RES"; break; }
												if(gf.getDescription().equals("CARBOHYD")) {glycProx++; fsNear="CARBOHYD"; break;}
												if(gf.getDescription().equals("DNA_BIND")) {dnaProx++; fsNear="DNA_BIND"; break;}
											}
										}
									}
									if(fsNear == null) {
										for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.PHOSPHOSITE, uniNbPos)) {
											if(!alreadyObserved.contains(f)) {	// avoid double counting of the same feature
												fsNear="MOD_RES";
												if(((PhosphoSitePlusFeature)f).getModType() == ProteinModificationType.PHOSPHORYLATION) phoProx++;
												else if(((PhosphoSitePlusFeature)f).getModType() == ProteinModificationType.UBIQUITINATION) ubqProx++;
												else modProx++;
												break;
											}
										}
									}
								}
							}
						}						
					}
				}
			}
		}
		int totalHits = actHits+phoHits+ubqHits+modHits+atpHits+gtpHits+glycHits+dnaHits;
		int totalProx = actProx+phoProx+ubqProx+modProx+atpProx+gtpProx+glycProx+dnaProx;
		if(!silent) System.out.println();
		if(!silent) System.out.printf("Active sites:          %3d hits %3d proximal %3d total\n", actHits, actProx, actHits+actProx);
		if(!silent) System.out.printf("Phosphorilation sites: %3d hits %3d proximal %3d total\n", phoHits, phoProx, phoHits+phoProx);
		if(!silent) System.out.printf("Ubiquitination sites:  %3d hits %3d proximal %3d total\n", ubqHits, ubqProx, ubqHits+ubqProx);
		if(!silent) System.out.printf("Glycosilation sites:   %3d hits %3d proximal %3d total\n", glycHits, glycProx, glycHits+glycProx);		
		if(!silent) System.out.printf("Other modifications:   %3d hits %3d proximal %3d total\n", modHits, modProx, modHits+modProx);
		if(!silent) System.out.printf("ATP binding sites:     %3d hits %3d proximal %3d total\n", atpHits, atpProx, atpHits+atpProx);
		if(!silent) System.out.printf("GTP binding sites:     %3d hits %3d proximal %3d total\n", gtpHits, gtpProx, gtpHits+gtpProx);		
		if(!silent) System.out.printf("DNA binding sites:     %3d hits %3d proximal %3d total\n", dnaHits, dnaProx, dnaHits+dnaProx);		
		if(!silent) System.out.printf("Total:                 %3d hits %3d proximal %3d total\n", totalHits, totalProx, totalHits+totalProx);
		if(!silent) System.out.printf("Number of mutations:   %3d\n", numMut);
		if(!silent) System.out.printf("all=[%3d %3d %3d %3d %3d %3d %3d];\n", actHits+actProx, phoHits+phoProx, ubqHits+ubqProx, glycHits+glycProx+modHits+modProx, atpHits+atpProx, gtpHits+gtpProx, numMut);
		double fracHits = 1.0 * totalHits / numMut;
		double fracProx = 1.0 * totalProx / numMut;
		double fracTotal = 1.0 * (totalHits+totalProx) / numMut;
		//double fracTotal = fracHits;
		if(!silent) System.out.printf("Fraction of hits:       %6.4f Prox: %6.4f Total: %6.4f\n", fracHits, fracProx, fracTotal);
		if(!silent && out != null) out.printf("%f\t%f\t%f\t%e\n", fracTotal, fracTotal, fracTotal, 0.0); // fraction, jackknife, p-value
		if(!silent) System.out.println();
		return fracTotal;
	}
	
//	/**
//	 * Loads the mutation information for all target genes in the list
//	 */
//	public void loadMutations() {
//		for(Gene g:this.targets) {
//			g.loadMutations();
//		}
//	}	
	
	/**
	 * Verifies that all mutations for which structures are known have the correct residue at the respective position
	 * in the structure. Prints error messages otherwise.
	 */
	public void checkMutations() {
		for(Gene g:this.targets) {
			System.out.println("-- " + g.geneName + " --");
			g.checkMutations();
		}		
	}
	
	/**
	 * Loads the mutation information for all target genes from the default database table.
	 */
	public void loadMutationsFromDb(MySQLConnection conn) {
		for(Gene g:this.targets) {
			g.loadMutations(conn, LOAD_ONLY_MISSENSE);
		}
	}
	
//	/**
//	 * Loads the mutation information for all target genes from the alternate database table.
//	 */
//	public void loadMutationsFromDb2(MySQLConnection conn) {
//		for(Gene g:this.targets) {
//			g.loadMutationsFromDb(conn);
//		}
//	}
	
	/**
	 * Load SNP data for all genes from the database.
	 */
	public void loadSnpsFromDb(MySQLConnection conn) {
		for(Gene g:this.targets) {
			g.loadSNPs(conn);
		}		
	}
	
//	public void loadSnpsFromDb2(MySQLConnection conn) {
//		for(Gene g:this.targets) {
//			g.loadSnpsFromDb2(conn);
//		}		
//	}
	
	/**
	 * Load SNP data as mutations for all genes.
	 */
	public void loadSNPsAsMutations(MySQLConnection conn) {
		for(Gene g:this.targets) {
			g.loadSNPsAsMutations(conn);
		}		
	}
	
//	public void loadSNPsAsMutations2(MySQLConnection conn) {
//		for(Gene g:this.targets) {
//			g.loadSNPsAsMutations(conn);
//		}		
//	}
	
	/**
	 * Remove all mutations which are not missense mutations
	 * @return the number of removed mutations or -1 if something went wrong
	 */
	public int filterMutationsMissense() {
		int numRemoved = 0;
		for(Gene g:this.targets) {
			if(g.mutations == null) {
				System.err.println("Error: No mutations loaded for " + g.geneName);
				continue;
			}
			HashSet<Mutation> toRemove = new HashSet<Mutation>();
			for(Mutation m:g.mutations) {
				if(!m.isMissense()) toRemove.add(m);
			}
			for(Mutation m:toRemove) {
				if(g.mutations.remove(m)) numRemoved++;
			}
		}
		return numRemoved; 
	}	
	
	/**
	 * Remove all SNPs which are not missense mutations
	 * @return the number of removed mutations or -1 if something went wrong
	 */
	public int filterSNPsMissense() {
		int numRemoved = 0;
		for(Gene g:this.targets) {
			if(g.snps == null) {
				System.err.println("Error: No snps loaded for " + g.geneName);
				continue;
			}
			HashSet<Mutation> toRemove = new HashSet<Mutation>();
			for(Mutation m:g.snps) {
				if(!m.isMissense()) toRemove.add(m);
			}
			for(Mutation m:toRemove) {
				if(g.snps.remove(m)) numRemoved++;
			}
		}
		return numRemoved; 
	}
	
	/**
	 * Removes all substructures which are not X-ray structures (for energy calculations)
	 * @return the number of structures removed or -1 if something went wrong
	 */
	public int filterSubstructuresXRay() {
		int numRemoved = 0;
		for(Gene g:this.targets) {
			if(!g.areSubstructuresLoaded()) {
				System.err.println("Skipping gene " + g.getGeneName() + ". Substructures not loaded.");
				continue;
			}
			LinkedList<Substructure> toRemove = new LinkedList<Substructure>();
			for(Substructure ss:g.getSubstructures()) {
				if(ss.getType() != SubstructureType.XRAY) toRemove.add(ss);
			}
			for(Substructure ss:toRemove) {
				if(g.getSubstructures().remove(ss)) numRemoved++;
			}
		}
		return 0;
	}
	
	/**
	 * Remove all mutations without known substructure
	 * @return the number of removed mutations or -1 if something went wrong
	 */
	public int filterMutationsKnownStructure() {
		int numRemoved = 0;
		for(Gene g:this.targets) {
			if(g.mutations == null) {
				System.err.println("Error: No mutations loaded for " + g.geneName);
				continue;
			}
			if(g.subStructures == null) {
				System.err.println("Error: No substructures loaded for " + g.geneName);
				continue;
			}
			HashSet<Mutation> toRemove = new HashSet<Mutation>();
			for(Mutation m:g.mutations) {
				if(g.getSubstructure(m.position) == null) toRemove.add(m);
			}
			for(Mutation m:toRemove) {
				if(g.mutations.remove(m)) numRemoved++;
			}
		}
		return numRemoved; 
	}
	
	/**
	 * Remove all SNPs without known substructure
	 * @return the number of removed mutations or -1 if something went wrong
	 */
	public int filterSNPsKnownStructure() {
		int numRemoved = 0;
		for(Gene g:this.targets) {
			if(g.snps == null) {
				System.err.println("Error: No SNPs loaded for " + g.geneName);
				continue;
			}
			if(g.subStructures == null) {
				System.err.println("Error: No substructures loaded for " + g.geneName);
				continue;
			}
			HashSet<Mutation> toRemove = new HashSet<Mutation>();
			for(Mutation m:g.snps) {
				if(g.getSubstructure(m.position) == null) toRemove.add(m);
			}
			for(Mutation m:toRemove) {
				if(g.snps.remove(m)) numRemoved++;
			}
		}
		return numRemoved; 
	}
	
	/**
	 * Remove all mutations without substructure or where the position in the substructure is not observed
	 * @return the number of removed mutations or -1 if something went wrong
	 */
	public int filterMutationsObserved() {
		int numRemoved = 0;
		for(Gene g:this.targets) {
			if(g.mutations == null) {
				System.err.println("Error: No mutations loaded for " + g.geneName);
				continue;
			}
			if(g.subStructures == null) {
				System.err.println("Error: No substructures loaded for " + g.geneName);
				continue;			
			}
			// mark mutations for deletion
			HashSet<Mutation> toRemove = new HashSet<Mutation>();
			for(Mutation m:g.mutations) {
				Substructure ss = null;
				if((ss = g.getSubstructure(m.position)) == null) {
					toRemove.add(m);
				} else {
					// mapping uniprot codes to pdb codes (use offset)
					int uniPos = m.position;
					int pdbPos = ss.mapUniprotResser2Cif(uniPos);
					if(pdbPos == Substructure.ALIGNMENT_UNDEFINED) {
						System.err.println("Warning: Unable to map residue number " + uniPos + " from Uniprot to Pdb for " + g.geneName + " " + ss);
					}
					if(ss.pdb == null) {
						System.err.println("Warning: No PdbChain object loaded for " + g.geneName + " " + ss.toString());
						toRemove.add(m);
					} else {
						if(!ss.pdb.containsStdAaResidue(pdbPos)) toRemove.add(m);
					}
				}
			}
			// actually delete mutations
			for(Mutation m:toRemove) {
				if(g.mutations.remove(m)) numRemoved++;
			}
		}
		return numRemoved; 
	}
	
	/**
	 * Remove all mutations which are in regions with predicted structure
	 * @return the number of removed mutations or -1 if something went wrong
	 */
	public int filterMutationsPredicted() {
		int numRemoved = 0;
		for(Gene g:this.targets) {
			if(g.mutations == null) {
				System.err.println("Error: No mutations loaded for " + g.geneName);
				continue;
			}
			if(g.subStructures == null) {
				System.err.println("Error: No substructures loaded for " + g.geneName);
				continue;
			}
			HashSet<Mutation> toRemove = new HashSet<Mutation>();
			for(Mutation m:g.mutations) {
				Substructure ss = g.getSubstructure(m.position);
				if(ss != null && ss.getType() == SubstructureType.PREDICTION) toRemove.add(m);
			}
			for(Mutation m:toRemove) {
				if(g.mutations.remove(m)) numRemoved++;
			}
		}
		return numRemoved; 		
	}
	
	/*--------------------------- abstract methods --------------------------*/
	/**
	 * Randomizes the current set of mutations.
	 * Modifies the currently loaded set of mutations (for all genes) such
	 * that within each structural domain the number of mutations stays the
	 * same but the positions are randomly chosen from the set of observed
	 * residues of the respective domain. Structures and Domains have to be
	 * loaded previously.
	 * TODO: Should this be a method of class Gene?
	 * @param withReplacement if false, makes sure that each position is only
	 * chosen once.
	 */
	public void randomizeMutations(boolean withReplacement) {
		for(Gene g:this.getTargets()) {
			//System.out.println(g.getGeneName());
			//System.out.println(g.mutations);
			//LinkedList<Mutation> newMutations = new LinkedList<Mutation>();
			HashSet<Mutation> newMutations = new HashSet<Mutation>();
			if(g.getFeatures().size() == 0) System.err.println("Warning. No features found when randomizing mutations for " + g.getGeneName());
			for(Feature f:g.getFeatures()) {
				if(f.getType() == FeatureType.SDOMAIN) {
					TreeSet<Integer> intSet = f.getIntervalSet().getIntegerSet();	// set of positions in current domain
					Substructure ss = g.getSubstructure(intSet.first());	// this should never be null because of the way we load domains
					if(ss == null) System.err.println("Error: Substructure for " + f.toString() + " should never be null");
					// count number of observed mutations in this domain
					int numMutations = 0;
					for(Mutation m:g.getMutations()) {
						if(intSet.contains(m.position) && ss.pdb.containsStdAaResidue(ss.mapUniprotResser2Cif(m.position))) {
							if(ss != g.getSubstructure(m.position)) System.err.printf("Assertion failed: ss(%d)=%s%s != ss(%s)\n", m.position, ss.getPdbCode(),ss.getChainCode(), f.toString());
							numMutations++;
						}
					}
					// get random positions in domain
					int c = 0;
					Set<Integer> resSet = new TreeSet<Integer>();	// this will be the set of allowed positions from which to select at random
					for(int i: intSet) {
						int pdbPos = ss.mapUniprotResser2Cif(i);
						if(ss.pdb.containsStdAaResidue(pdbPos)) resSet.add(pdbPos);
					}
					Integer[] resArr = new Integer[resSet.size()];	// the same in array form so that random indices can be drawn
					resArr = resSet.toArray(resArr);
					Random rand = new Random();
					Set<Integer> rndIdx = new HashSet<Integer>();	// remember indices for case w/o replacement
					if(resArr.length <= 0) System.err.printf("Error: resArr.length=%d in %s #mut=%d\n", resArr.length, f.toString(), numMutations);
					// we had problems that this loop would not terminate
					while(c < numMutations) {
						int r = rand.nextInt(resArr.length);		
						if(withReplacement || !rndIdx.contains(r)) {
							rndIdx.add(r);
							int rndPdbPosition = resArr[r]; // convert random index to pdb position 
							int uniPos = ss.mapCifResser2Uniprot(rndPdbPosition); // convert to uniprot
							Mutation newMut = new Mutation(AminoAcid.XXX, AminoAcid.XXX, uniPos);
							newMutations.add(newMut);
							c++;
						}
					}					
				}
			}
			
// This is the old way of randomizing within whole substructures, after testing the new way, this can be deleted			
//			for(Substructure ss:g.getSubstructuresWithMutations()) {
//				int numMutations = 0;				
//				for(Mutation m:g.getMutations()) {
//					// count number of mutations in this substructure
//					if(g.getSubstructure(m.position) == ss) {
//						if(ss.pdb.hasCoordinates(ss.mapUniprotResser2Pdb(m.position))) {
//							numMutations++;
//						}
//					}
//				}
//				// get random positions in observed substructure
//				int c = 0;
//				Set<Integer> resSet = ss.pdb.getAllSortedResSerials();
//				Integer[] resArr = new Integer[resSet.size()];
//				resArr = resSet.toArray(resArr);
//				Random rand = new Random();
//				Set<Integer> rndIdx = new HashSet<Integer>();	// remember indices for case w/o replacement
//				while(c < numMutations) {
//					int r = rand.nextInt(resArr.length);		
//					if(withReplacement || !rndIdx.contains(r)) {
//						rndIdx.add(r);
//						int rndPdbPosition = resArr[r]; // convert random index to pdb position 
//						int uniPos = ss.mapPdbResser2Uniprot(rndPdbPosition); // convert to uniprot
//						Mutation newMut = new Mutation(AminoAcid.XXX, AminoAcid.XXX, uniPos);
//						newMutations.add(newMut);
//						c++;
//					}
//				}
//			}
				
			g.mutations = newMutations;
			//System.out.println(newMutations);
		}
	}
	
	/**
	 * Writes a text file with information about each gene, and data files for each gene
	 * with information about substructures and mutations.
	 * @param outDir
	 */
	public void writeGeneDataFiles(File outDir, SubstructureType restrictToType) {
		try {
			File outFile = new File(outDir,GENE_FILE_NAME);
			System.out.println("Writing " + outFile);
			PrintWriter out = new PrintWriter(outFile);
			File substructureFile;
			File mutationFile;
			File snpFile;
			for(Gene g:this.targets) {
				out.printf("%s\t%d\n", g.geneName, g.length);
				substructureFile = new File(outDir, g.geneName + ".substructures.txt");
				mutationFile = new File(outDir, g.geneName + ".mutations.txt");
				snpFile = new File(outDir, g.geneName + ".snps.txt");
				System.out.println("Writing " + substructureFile);
				g.writeSubstructureData(substructureFile, restrictToType);
				System.out.println("Writing " + mutationFile);
				g.writeMutationData(mutationFile, restrictToType);
				System.out.println("Writing " + snpFile);
				g.writeSnpData(snpFile, restrictToType);
			}				
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Writes the data used for this analysis to files (e.g. for FTP download)
	 * @param outDir
	 * @param append if true, skip header and append to outfile files
	 */
	public void writeOriginalDataFiles(File outDir, boolean append) {
		
		System.out.println("Writing original data to " + outDir);
		File geneFile = new File(outDir, "genes.txt");
		File mutFile = new File(outDir, "mutations.txt");
		File snpFile = new File(outDir, "snps.txt");
		File annotFile = new File(outDir, "annotations.txt");
		ArrayList<Mutation> muts;
		
		// genes
		try{
			System.out.println("--- Genes ---");
			PrintWriter out = new PrintWriter(new FileWriter(geneFile, append));
			if(!append) out.printf("%s\t%s\t%s\t%s\n", "#name", "unip_id", "len_aa", "class");			
			for(Gene g:this.getTargets()) {
				out.printf("%s\t%s\t%d\t%s\n", g.getGeneName(), g.getUniprotId(), g.getLength(), "Onc");
			}
			out.close();
		} catch(IOException e) {
			System.err.println("Error writing to file " + geneFile + ": " + e.getMessage());
		}

		// mutations
		try{
			System.out.println("--- Mutations ---");
			PrintWriter out = new PrintWriter(new FileWriter(mutFile, append));
			if(!append) out.printf("%s\t%s\t%s\t%s\t%s\t%s\n", "#gene", "pos_aa", "aa_wt","aa_mut", "cosmic_id", "mut_dna");
			for(Gene g:this.getTargets()) {
				muts = new ArrayList<Mutation>(g.getMutations());
				Collections.sort(muts, new Comparator<Mutation>() {
					public int compare(Mutation m1, Mutation m2) {
						if(m1.position > m2.position) return 1;
						if(m1.position < m2.position) return -1;
						return m1.after.compareTo(m2.after);
					}
				});
				for(Mutation m:muts) {
					if(m.type == MutType.MISSENSE && g.getSubstructure(m.position) != null) {
						Substructure ss = g.getSubstructure(m.position);
						if(ss.pdb.containsStdAaResidue(ss.mapUniprotResser2Cif(m.position))) {
							out.printf("%s\t%d\t%s\t%s\t%d\t%s\n", g.getGeneName(), m.position, m.before.getThreeLetterCode(), m.after.getThreeLetterCode(), m.mutId, m.dnaMutStr);
						}
					}
				}
			}
			out.close();
		} catch(IOException e) {
			System.err.println("Error writing to file " + mutFile + ": " + e.getMessage());
		}

		// snps
		try{
			System.out.println("--- SNPs ---");
			PrintWriter out = new PrintWriter(new FileWriter(snpFile, append));
			if(!append) out.printf("%s\t%s\t%s\t%s\t%s\n", "#gene", "pos_aa", "aa_wt","aa_mut", "rsId");
			for(Gene g:this.getTargets()) {
				muts = new ArrayList<Mutation>(g.getSNPs());
				Collections.sort(muts, new Comparator<Mutation>() {
					public int compare(Mutation m1, Mutation m2) {
						if(m1.position > m2.position) return 1;
						if(m1.position < m2.position) return -1;
						return m1.after.compareTo(m2.after);
					}
				});
				for(Mutation m:muts) {
					SNP snp = (SNP) m;
					if(g.getSubstructure(m.position) != null) {
						out.printf("%s\t%d\t%s\t%s\t%d\n", g.getGeneName(), m.position, m.before.getThreeLetterCode(), m.after.getThreeLetterCode(), snp.rsId);
					}
				}
			}
			out.close();
		} catch(IOException e) {
			System.err.println("Error writing to file " + snpFile + ": " + e.getMessage());
		}

//		// annotations
//		try{
//			System.out.println("--- Annotations ---");
//			PrintWriter out = new PrintWriter(new FileWriter(annotFile));
//			out.printf("%s\t%s\t%s\t%s\n","#gene", "pos", "source", "descr");
//			for(Gene g:this.getTargets()) {			
//				for(Feature f: g.getFeatures()) {
//					if(f.getType() != FeatureType.SDOMAIN) {
//						String type = f.getType().toString();
//						if(type.equals("GENERAL")) type = "MANUAL";
//						out.printf("%s\t%s\t%s\t%s\n",g.getGeneName(), f.getIntervalSet(), type, f.getDescription());
//					}
//				}
//			}
//			out.close();
//		} catch(IOException e) {
//			System.err.println("Error writing to file " + annotFile + ": " + e.getMessage());
//		}
		
		// annotations
		try{
			System.out.println("--- Annotations ---");
			PrintWriter out = new PrintWriter(new FileWriter(annotFile, append));
			if(!append) out.printf("%s\t%s\t%s\t%s\n", "#gene", "location", "type", "source");
			for(Gene g:this.getTargets()) {
				if(!g.areSubstructuresLoaded()) {
					System.err.println("Skipping gene " + g.getGeneName() + ": substructures not loaded.");
				}
				HashSet<Feature> alreadyObserved = new HashSet<Feature>();
				for(Substructure ss:g.getSubstructures()) {
					//				if(ss.getType() == SubstructureType.PREDICTION) {
					//					System.out.println("Skipping predicted structure " + g.getGeneName() + "_" + ss.getRange());
					//					continue;
					//				}
					if(!ss.isPdbLoaded()) {
						System.err.println("Skipping " + g.getGeneName() + " " + ss.getRange() + " " + ss.getPdbCode()+ss.getChainCode() + ": Structure not loaded.");
					} else {
						String cifChain = ss.pdb.getChainCode();
						String pdbChain = ss.pdb.getPdbChainCode(); // = ss.getChainCode();
						if(!pdbChain.equals(cifChain)) {
							System.err.println("Warning: pdb chain code " + pdbChain + " != cif chain code " + cifChain);
						}
						for(int cifPos: ss.pdb.getAllStdAaResSerials()) {	// only observed residues
							//int pdbPos = Integer.parseInt(ss.pdb.getResidue(p).getPdbSerial());	// cif -> pdb serial
							int uniprotPos = ss.mapCifResser2Uniprot(cifPos);
							//AminoAcid wtAa = ss.pdb.getResidue(p).getAaType();
							// check for functional sites
							Pair<Feature, String> featureAndFunction = g.getFirstFeatureAndFunctionalClass(uniprotPos, alreadyObserved);
							Feature f = featureAndFunction.getFirst();
							String fsHit = featureAndFunction.getSecond();
							if(!fsHit.equals("NONE")) {
								String type = f.getType().toString();
								if(type.equals("GENERAL")) type = "MANUAL";
								out.printf("%s\t%s\t%s\t%s\n", g.getGeneName(), f.getIntervalSet(), fsHit, type);
							}
						}
					}
				}
			}
			out.close();
		} catch(IOException e) {
			System.err.println("Error writing to file " + annotFile + ": " + e.getMessage());
		}

	}
	
	/**
	 * Renders png images with sequence rulers for the genes in this list.
	 * @param outDir where the output files will be written to
	 */
	public void renderRulers(File outDir) {
		// rendering images
		File pngFile;
		for(Gene g:this.targets) {
			pngFile = new File(outDir, g.geneName + ".png");
			g.renderRuler(pngFile.toString(), 800);
		}		
	}
	
	/**
	 * Renders a graphical overview for the targets in the list.
	 * @param outDir where the output files will be written to
	 */
	public void renderOverview(File outDir) {
		try {
			File htmlFile = new File(outDir, "index.html");
			PrintWriter out = new PrintWriter(htmlFile);
			out.println("<html>");
			
			// rendering images
			File pngFile;
			for(Gene g:this.targets) {
				pngFile = new File(outDir, g.geneName + ".png");
				g.renderStructuralRegions(pngFile.toString(), 800);
				out.println("<h3>" + g.geneName + "</h3>");
				out.println("<img src=\"" + pngFile.toString() + "\" /><br />");
			}

			out.println("</html>");			
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/* ---------------------- TODO ------------------------ */
		
	/**
	 * Loads sequences for all genes in this list from COSMIC.
	 * Accession is by gene name, so names have to match COSMIC notation.
	 * @param outFile if not null write all sequences in multi-fasta format to this file
	 * @param online if true, load from COSMIC ftp server, otherwise use local files defined in Gene.COSMIC_SEQ_PATH
	 */
	public void loadCosmicSequences(File outFile, boolean online) {
		PrintStream out = null;
		try {
			if (outFile != null) {
				System.out.println("Writing Cosmic sequences to " + outFile);
				out = new PrintStream(outFile);
			}
			for (Gene g : this.targets) {
				String s = g.loadCosmicSequence(online);
				//System.out.println(g.geneName + ":\t" + s);
				if(out != null) {
					Sequence seq = new Sequence(g.geneName, s);
					seq.writeToPrintStream(out);
				}
			}
			if(out != null) out.close();
		} catch (IOException e) {
			System.err.println("Error writing to file " + outFile + ": " + e.getMessage());
		}		
	}
	
	/**
	 * Loads sequences for all genes in this list from the Uniprot website.
	 * Accession is by UniprotID, so annotated Ids are assumed to be correct.
	 * @param outFile if not null write all sequences in multi-fasta format to this file
	 * @param online if true, load from online UniProt, otherwise use local files defined in Gene.UNIPROT_SEQ_PATH
	 */
	public void loadUniprotSequences(File outFile, boolean online) {
		PrintStream out = null;
		try {
			if (outFile != null) {
				System.out.println("Writing Uniprot sequences to " + outFile);
				out = new PrintStream(outFile);
			}
			for (Gene g : this.targets) {
				String s = g.loadUniprotSequence(online);
				//System.out.println(g.geneName + ":\t" + s);
				if(out != null) {
					Sequence seq = new Sequence(g.geneName, s);
					seq.writeToPrintStream(out);
				}
			}
			if(out != null) out.close();
		} catch (IOException e) {
			System.err.println("Error writing to file " + outFile + ": " + e.getMessage());
		}
	}
	
	/**
	 * Aligns Uniprot vs. Cosmic sequence for each gene in the list.
	 * Reports only mismatches.
	 * @param online if true, load sequeneces from online UniProt and COSMIC, otherwise use local files defined in Gene.*_SEQ_PATH
	 * @return true if all sequences match, false otherwise
	 */
	public boolean compareUniprotVsCosmicSequence(boolean online) {		
		boolean allEqual = true;
		for(Gene g:this.targets) {
			if(g.cosmicSeq == null) g.loadCosmicSequence(online);
			if(g.uniprotSeq == null) g.loadUniprotSequence(online);
			if(!g.cosmicSeq.equals(g.uniprotSeq)) {
				allEqual = false;
				System.out.println("Warning for gene " + g.geneName);
				System.out.println("Cosmic sequence");
				System.out.println(g.cosmicSeq);
				System.out.println("does not match Uniprot sequence");
				System.out.println(g.uniprotSeq);
				
				// align cosmic vs. uniprot				
				try {
					PairwiseSequenceAlignment pa = new PairwiseSequenceAlignment(g.cosmicSeq, g.uniprotSeq, "cosmic", "uniprot");
					pa.printReport();
				} catch (PairwiseSequenceAlignmentException e) {
					System.err.println("Error. Could not create sequence alignment: " + e.getMessage());
				}
			}
		}		
		return allEqual;
	}

// OBSOLETE:
//	/**
//	 * Downloads PDB files for the experimental substructures of genes in this list,
//	 * trying to do renumbering of residues according to Uniprot sequence.
//	 * @param targetDir
//	 */
//	public void downloadPdbsForSubstructures(File targetDir) {
//		File pdbFile = null;
//		for(Gene g:this.targets) {
//			for(Substructure s:g.getSubstructures()) {
//				pdbFile = new File(targetDir, s.getDefaultPdbFileName(g.getGeneName()));
//				File newFile = s.downloadStructureFromPdb(pdbFile);
//				if(newFile != null) {
//					// renumber
//					// TODO: this has to be improved because for some structures, a simple offset is not enough
//					int offset = Substructure.getPdb2UniprotOffset(s.getPdbCode(), s.getChainCode());
//					File pdbRenum = new File(targetDir, g.geneName + "_" + s.getRange() + "_" + s.getPdbCode() + ".renum.pdb");
//					Substructure.renumberPdb(pdbFile, pdbRenum, offset, s.getChainCode().charAt(0));
//				}
//			}
//		}
//	}
	
	/**
	 * Checks for each substructure whether there exists a pdb file.
	 */
	public void checkPdbFiles() {
		//File pdbDir = new File("/project/PyMol/pdbs_final");
		File pdbDir = new File("/home/web/lappe/stehr/mutanom/pdb/");
		String fileName;
		for(Gene g:this.targets) {
			for(Substructure s:g.getSubstructures()) {
				if(s.getType() == SubstructureType.PREDICTION) {
					fileName = g.geneName + "_" + s.getRange() + "_" + "pred.renum.pdb";
				} else {
					fileName = g.geneName + "_" + s.getRange() + "_" + s.getPdbCode() + ".renum.pdb";
				}
				if(!new File(pdbDir, fileName).canRead()) {
					System.out.println("Missing pdb file: " + fileName);
				}
			}
		}
	}
	
	/**
	 * Generated a pymol session and png preview for each mutation
	 * @param outDir
	 */
	public void generatePngsAndPyMolSessions(File outDir, File pdbDir) {		
			for(Gene g:this.targets) {
				for(Mutation m:g.getMutations()) {
					String baseName = m.getFileBasename();
					try {
						System.out.println("Writing files for " + baseName);
						m.writePngAndPyMolSession(outDir, baseName, pdbDir);
					} catch (IOException e) {
						System.err.println("Error writing PyMol session files and/or PNG previews for " + baseName + ": " + e.getMessage());
					}
				}
			}
			System.out.println("Warning: PyMol jobs may be still running");
	}
	
	public void generatePseWithAllMutations(File outDir, File pdbDir) {
		for(Gene g:this.targets) {
			String baseName = g.getGeneName();
			try {
				System.out.println("Writing session file for " + baseName);
				g.writePseWithAllMutations(outDir, baseName, pdbDir);
			} catch (IOException e) {
				System.err.println("Error writing PyMol session files and/or PNG previews for " + baseName + ": " + e.getMessage());
			}
		}
	}
	
	/**
	 * Generates for each gene a PyMol session (.pse) file with all substructures as separate objects and mutations and SNPs mapped onto the structures
	 * @param sessionFile the output file which is being written
	 * @param pdbDir directory where the pdb files for all substructures are found
	 * @throws IOException
	 */
	public void generatePseWithAllMutationsAndSNPs(File outDir, File pdbDir) throws IOException {
		PyMolScriptMaker psm;
		File sessionFile;
		for(Gene g:this.targets) {
			System.out.println(g.geneName);
			sessionFile = new File(outDir,g.geneName + ".pse");
			psm = new PyMolScriptMaker(true);
			psm.loadDefaultSettings();
			for(Substructure ss:g.subStructures) {
				String pdbFileName = ss.getDefaultPdbFileName(g.geneName);		
				String objName = g.geneName + "_" + ss.getRange();
				psm.load(new File(pdbDir, pdbFileName), objName);
//				System.out.println();
//				for(Mutation m:g.getSNPs()) {
//					System.out.print(m + " ");
//					int pdbPos = ss.mapUniprotResser2Pdb(m.position);
//					String selName = String.format("s_%s%d", m.before.getOneLetterCode(), pdbPos);			
//					//selName = g.geneName + "_snps";
//					psm.highlightResidue(objName, m.position, ss.getChainCode().charAt(0), Color.YELLOW, selName);
//				}
//				for(Feature f:g.getFeatures()) {
//					//System.out.print(f + " ");
//					for(int pos:f.getIntervalSet().getIntegerSet()) {
//						//int pos = f.getIntervalSet().getIntegerSet().first();
//						int pdbPos = ss.mapUniprotResser2Pdb(pos);
//						String selName = String.format("f%d_%s", pdbPos, f.getDescription());
//						System.out.println(f);
//						//selName = g.geneName + "_features";
//						psm.highlightResidue(objName, pdbPos, ss.getChainCode().charAt(0), Color.YELLOW, selName);
//					}
//				}
				for(Mutation m:g.getMutations()) {
					int pdbPos = ss.mapUniprotResser2Cif(m.position);
					if(pdbPos != Substructure.ALIGNMENT_UNDEFINED) {
					//System.out.print(m + " ");
					String selName = String.format("m_%s%d", m.before.getOneLetterCode(), pdbPos);
					//selName = "mutations";
					psm.highlightResidue(objName, pdbPos, ss.getChainCode().charAt(0), Color.RED, selName, null);
					}
				}
			}
			psm.saveSession(sessionFile);
			psm.executeAndClose();
		}
	}
	
	/**
	 * Generates for each gene a PyMol session (.pse) file with all substructures as separate objects\
	 * and mutations, SNPs, functional sites and domains mapped onto the structures.
	 * @param outDir the output directory where the pse files will be written to
	 * @param pdbDir directory with the pdb files for the substructures, name format: 1xyzA.pdb, numbering=cif
	 * @throws IOException
	 */
	public void generatePseWithMutSnpFuncAndDomains(File outDir, File pdbDir) throws IOException {
		PyMolScriptMaker.Color[] domainColors = {Color.WHEAT, Color.PALEGREEN, Color.LIGHTBLUE, Color.LIGHTPINK, Color.PALEYELLOW, Color.PALECYAN, Color.LIGHTORANGE, Color.BLUEWHITE};
		PyMolScriptMaker psm;
		File sessionFile;
		for(Gene g:this.targets) {
			System.out.println(g.geneName);
			sessionFile = new File(outDir,g.geneName + ".pse");
			psm = new PyMolScriptMaker(true);
			psm.loadDefaultSettings();
			for(Substructure ss:g.subStructures) {
				String chain = ss.getPdb().getChainCode();	// cif chain code
				String pdbChain = ss.getPdb().getPdbChainCode(); // pdb chain code
				String pdbFileName = ss.getPdbCode() + pdbChain + ".pdb";	
				//String objName = g.geneName + "_" + ss.getRange() + "_" + ss.getPdbCode() + chain;
				String objName = g.geneName + "_" + ss.getPdbCode() + pdbChain;
				psm.load(new File(pdbDir, pdbFileName), objName);
				// DOM			
				int i = 0;
				for(Feature f:g.getFeaturesOfType(FeatureType.SDOMAIN)) {
					for(Interval intv:f.getIntervalSet()) {
						String selName = "d_" + f.getDescription();
						psm.highlightInterval(objName, intv, chain.charAt(0), domainColors[i], selName);
					}
					i++;
				}
				// SNP
				String snpGroup = ss.getPdbCode() + "_snp";
				for(Mutation m:g.getSNPs()) {
					int pdbPos = ss.mapUniprotResser2Cif(m.position);
					if(pdbPos != Substructure.ALIGNMENT_UNDEFINED) {
						String selName = String.format("s_%s%d", m.before.getOneLetterCode(), pdbPos);			
						psm.highlightResidue(objName, m.position, chain.charAt(0), Color.YELLOW, selName, snpGroup);
					}
				}
				// FUNC
				String funcGroup = ss.getPdbCode() + "_func";
				for(Feature f:g.getFeatures()) {
					if(f.getType() != FeatureType.SDOMAIN) {
						for(int pos:f.getIntervalSet().getIntegerSet()) {
							int pdbPos = ss.mapUniprotResser2Cif(pos);
							if(pdbPos != Substructure.ALIGNMENT_UNDEFINED) {
								String selName = String.format("f%d_%s", pdbPos, f.getDescription());
								psm.highlightResidue(objName, pdbPos, chain.charAt(0), Color.BLUE, selName, funcGroup);
							}
						}
					}
				}
				// MUT
				String mutGroup = ss.getPdbCode() + "_mut";				
				for(Mutation m:g.getMutations()) {
					int pdbPos = ss.mapUniprotResser2Cif(m.position);
					if(pdbPos != Substructure.ALIGNMENT_UNDEFINED) {
						String selName = String.format("m_%s%d", m.before.getOneLetterCode(), pdbPos);
						psm.highlightResidue(objName, pdbPos, chain.charAt(0), Color.RED, selName, mutGroup);
					}
				}
			}
			psm.saveSession(sessionFile);
			psm.executeAndClose();
			// wait a bit because too many open pymol processes may cause memory shortage
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * Prints a list of prosite motifs most frequently found in this set of genes.
	 * Motifs are counted only once per gene to avoid bias.
	 */
	public void findFrequentPrositeMotifs() {
		HashMap<String, Integer> motifCount = new HashMap<String, Integer>();
		PrositeScanner ps = new PrositeScanner();
		// count motifs
		for(Gene g:this.targets) {
			System.out.println("Scan prosite for " + g.geneName);
			ps.scanSeq(g.getUniprotId());
			HashSet<String> localHits = new HashSet<String>();
			for(PrositeHit hit:ps.getLatestHits()) {
				// makes sure that every motif is only counted once per gene
				localHits.add(hit.signatureAc);
			}
			for(String motif:localHits) {
				if(motifCount.containsKey(motif)) {
					int oldCount = motifCount.get(motif);
					int newCount = oldCount + 1;
					motifCount.put(motif, newCount);
				} else {
					motifCount.put(motif,1);
				}
			}
		}
		
		// print motifs
		LinkedHashMap<String,Integer> sm = Goodies.sortMapByValue(motifCount, false);
		System.out.println();
		System.out.println("Most frequent prosite motifs in this set:");
		for(String motifAc:sm.keySet()) {
			int count = sm.get(motifAc);
			String descr = PrositeScanner.getPatternDescription(motifAc);
			System.out.printf("%3d %s %s\n", count, motifAc, descr);
		}
	}
	
	/**
	 * Prints a list of uniprot sequence features frequently found in this set of genes.
	 * Motifs are counted only once per gene to avoid bias.
	 */
	public void findFrequentUniprotMotifs() {
		HashMap<String, Integer> motifCount = new HashMap<String, Integer>();
		UniProtConnection uc = new UniProtConnection();
		// count motifs
		for(Gene g:this.targets) {
			System.out.print("Scan Uniprot for " + g.geneName);
			UniProtEntry upe;
			try {
				upe = uc.getEntry(g.uniProtID);
				HashSet<String> localHits = new HashSet<String>();
				for(uk.ac.ebi.kraken.interfaces.uniprot.features.Feature hit:upe.getFeatures()) {
					// makes sure that every motif is only counted once per gene
					String key = hit.getType().getName() + " " + hit.getType().getDisplayName();
					System.out.print(".");
	//				if(hit.getType() == FeatureType.NP_BIND) {
	//					System.out.print(".");
	//					key = ((HasFeatureDescription) hit).getFeatureDescription().getValue();
	//					localHits.add(key);
	//				}
					localHits.add(key);
				}
				System.out.println();
				for(String motif:localHits) {
					if(motifCount.containsKey(motif)) {
						int oldCount = motifCount.get(motif);
						int newCount = oldCount + 1;
						motifCount.put(motif, newCount);
					} else {
						motifCount.put(motif,1);
					}
				}
			} catch (NoMatchFoundException e) {
				System.err.println("No UniProt entry found for" +  g.uniProtID);
			}
		}
		
		// print features
		LinkedHashMap<String,Integer> sm = Goodies.sortMapByValue(motifCount, false);
		System.out.println();
		System.out.println("Most frequent UniProt feature types in "+ this.getTargets().size() + " proteins:");
		for(String motifAc:sm.keySet()) {
			int count = sm.get(motifAc);
			System.out.printf("%3d %s\n", count, motifAc);
		}
	}
	
	/**
	 * Print the pdb codes of all substructures in this target list
	 * @param printChainCodes if true, print chain codes as well
	 */
	public void printPdbCodes(boolean printChainCodes) {
		for(Gene g:this.targets) {
			for(Substructure ss:g.getSubstructures()) {
				String code = ss.getPdbCode();
				if(printChainCodes) code += ss.getChainCode();
				System.out.println(code);
			}
		}
	}
	
	/**
	 * Writes the file individuals_list.txt for each substructure in this list. 
	 * The file conforms to the FoldX format for evaluating mutations using <BuildModel>.
	 */
	public void writeIndividualLists(File outDir) {
		try {
		for(Gene g:this.getTargets()) {
			Collection<Substructure> substructures = g.getSubstructures();
			if(substructures == null) {
				System.out.println("No substructures with mutations found for " + g.getGeneName());
			} else {
				for(Substructure ss:substructures) {
					if(ss == null) System.err.println("Unexpected error: ss == null");
					String pdbCode = ss.getPdbCode();
					String chainCode = ss.getChainCode();
					String pdbChainCode = ss.pdb.getPdbChainCode();
					if(ss.getType() != SubstructureType.PREDICTION) {
						if(ss.pdb != null) {
							//File outFile = new File(outDir, "/" + pdbCode + "/" + FoldXRunner.INDIVIDUAL_LIST_FILE);
							File outFile = new File(outDir, pdbCode + ".mut");	// original pdb file with whole complex
							File outFile2 = new File(outDir, pdbCode + chainCode + ".mut"); // file from pdbase with single chain
							PrintWriter out = new PrintWriter(outFile);
							PrintWriter out2 = new PrintWriter(outFile2);
							for(Mutation m:g.mutations) {
								if(g.getSubstructure(m.position) == ss) {
									int posUni = m.position;				// position in uniprot sequence
									int posCif = ss.uni2cif.get(posUni);	// position in cif structure (pdbase)
									int posPdb = Integer.parseInt(ss.pdb.getPdbResSerFromResSer(posCif));	// position in original pdb file
									out.println(m.before.getOneLetterCode() + pdbChainCode + posPdb + m.after.getOneLetterCode()+";");
									out2.println(m.before.getOneLetterCode() + chainCode + posCif + m.after.getOneLetterCode()+";");
									if(!chainCode.equals(pdbChainCode)) System.out.println(chainCode + " != " + pdbChainCode);
								}
							}
							out.close();
							out2.close();
							System.out.println("Wrote file " + outFile);
							System.out.println("Wrote file " + outFile2);
						} else {
							System.err.println("No structure loaded for " + g.getGeneName() + " " + ss.getRange() + " " + ss.getPdbCode()+ss.getChainCode());
						}
					} else {
						System.out.println("Skipping predicted substructure " + g.getGeneName() + " " + ss.getRange());
					}
				}
			}
		}
		} catch(FileNotFoundException e) {
			System.err.println("Error writing to " + outDir + ": " + e.getMessage());
		}
	}
	
	/**
	 * Writes a shell script which can be called to run a comprehensive mutational
	 * screening (all positions and all possible target amino acids) on the cluster.
	 * The results will be written to the defined outDir.
	 * @throws FileNotFoundException 
	 */
	public void writeFoldXMutScanMasterShellScript(File scriptFile, File inDir, File outDir, boolean useSingleChainPdbFiles) {
		try {
			FoldXRunner.writeMutScanMasterShellScript(scriptFile, this, inDir, outDir, useSingleChainPdbFiles);
		} catch (FileNotFoundException e) {
			System.err.println("Error writing FoldXMutScanMasterScriptFile:" + e.getMessage());
		}
	}
	
	/**
	 * Parses the NACCESS output files for the mutanom target's pdb structures and writes the results to db.
	 * Creates the db table if not exists. Expects the residue table to be created previously using
	 * createAllPositionsDbTable().
	 * @param conn
	 * @param resTable the residue table which has be created previously
	 * @param outTable the table where the results will be written to (will be created if necessary)
	 * @throws SQLException 
	 */
	public void loadExposureToDb(MySQLConnection conn, String resTable, String outTable, File rsaDir) throws IOException, SQLException {
		
		// make sure that residue table exists
		int residues = conn.getIntFromDb("SELECT COUNT(*) FROM " + resTable);
		if(residues <= 0) {
			System.err.println("Error. Table " + resTable + " not found or is empty.");
			return;
		}
		
		// create output table
		System.out.println("Creating table (if not exists) " + outTable + "...");
		String sqlCreate = "CREATE TABLE IF NOT EXISTS " + outTable + "(res_idx int, rsa_sc float, rsa_cx float)";
		conn.executeSql(sqlCreate);
		
		System.out.print("Loading NACCESS results");
		
		int naccessResults = 0;
		for(Gene g:this.getTargets()) {
			if(!g.areSubstructuresLoaded()) {
				System.err.println("Skipping gene " + g.getGeneName() + ": substructures not loaded.");
			}
			for(Substructure ss:g.getSubstructures()) {
				//if(ss.getType() == SubstructureType.PREDICTION) {
				//	System.out.println("Skipping predicted structure " + g.getGeneName() + "_" + ss.getRange());
				//	continue;
				//}
				if(!ss.isPdbLoaded()) {
					System.err.println("Skipping " + g.getGeneName() + " " + ss.getRange() + " " + ss.getPdbCode()+ss.getChainCode() + ": Structure not loaded.");
				} else {
					System.out.print(".");				
					String pdbCode = ss.getPdbCode();
					String cifChain = ss.pdb.getChainCode();
					String pdbChain = ss.pdb.getPdbChainCode(); // = ss.getChainCode();
					if(!pdbChain.equals(cifChain)) {
						//System.err.println("Warning: pdb chain code " + pdbChain + " != cif chain code " + cifChain);
					}
					//System.out.print(pdbCode);
					//File rsaScFile = new File(rsaDir, pdbCode + pdbChain + ".rsa");
					File rsaCxFile = new File(rsaDir, pdbCode + ".rsa");
					String line;
					//TreeMap<String,Double> rsasSc = new TreeMap<String, Double>();
					TreeMap<String,Double> rsasCx = new TreeMap<String, Double>();
					BufferedReader in;
					// read single chain files
//					in = new BufferedReader(new FileReader(rsaScFile));
//					while((line = in.readLine()) != null) {
//						if(line.startsWith("RES")) {
//							int resNum = Integer.parseInt(line.substring(9, 13).trim());
//							String res = line.substring(4,7);
//							char chain = line.charAt(8);
//							double rsa = Double.parseDouble(line.substring(23,28));
//							String key = String.format("%s%04d%s",chain,resNum,res);
//							rsasSc.put(key, rsa);
//						}
//					}
//					in.close();
					// read complex file
					in = new BufferedReader(new FileReader(rsaCxFile));
					while((line = in.readLine()) != null) {
						if(line.startsWith("RES")) {
							int resNum = Integer.parseInt(line.substring(9, 13).trim());
							String res = line.substring(4,7);
							char chain = line.charAt(8);
							double rsa = Double.parseDouble(line.substring(23,28));
							String key = String.format("%s%04d%s",chain,resNum,res);
							rsasCx.put(key, rsa);
						}
					}
					in.close();
					// write rsa values to database
					for(int p: ss.pdb.getAllStdAaResSerials()) {	// only observed residues
						int pdbPos = Integer.parseInt(ss.pdb.getResidue(p).getPdbSerial());	// cif -> pdb serial
						int cifPos = p;
						AminoAcid wtAa = ((AaResidue)ss.pdb.getResidue(p)).getAaType();		
//						String keySc = String.format("%s%04d%s",cifChain,cifPos,wtAa.getThreeLetterCode());
//						if(!rsasSc.containsKey(keySc)) {
//							System.err.println("Could not find RSA entry for " + pdbCode + " " + keySc);
//						} else {
//							double rsaSc = rsasSc.get(keySc);
//							String sql = "UPDATE " + table + " SET rsa_sc=%f WHERE pdb_code='%s' AND pdb_chain='%s' AND cif_pos=%d";
//							conn.executeSql(String.format(sql, rsaSc, pdbCode, pdbChain, cifPos));
//						}
						String keyCx = String.format("%s%04d%s",pdbChain,pdbPos,wtAa.getThreeLetterCode());
						if(!rsasCx.containsKey(keyCx)) {
							System.err.println("Could not find RSA entry for " + pdbCode + " " + keyCx);
						} else {
							double rsaCx = rsasCx.get(keyCx);
							String query = "SELECT idx FROM " + resTable + " WHERE pdb_code='%s' AND pdb_chain='%s' AND cif_pos=%d";
							int resIdx = conn.getIntFromDb(String.format(query, pdbCode, pdbChain, cifPos));
							if(resIdx < 0) {
								String err = "Could not find entry for pdb_code='%s', pdb_chain='%s', cif_pos=%d in table " + resTable + ". Skipping.";
								System.err.println(String.format(err, pdbCode, pdbChain, cifPos));
								System.err.println(String.format(query, pdbCode, pdbChain, cifPos));
							} else {
								String sql = "INSERT INTO " + outTable + "(res_idx, rsa_cx) VALUES (%d, %f)";
								conn.executeSql(String.format(sql, resIdx, rsaCx));
								naccessResults++;
								// old:
								//String sql = "UPDATE " + table + " SET rsa_cx=%f WHERE pdb_code='%s' AND pdb_chain='%s' AND cif_pos=%d";
								//conn.executeSql(String.format(sql, rsaCx, pdbCode, pdbChain, cifPos));
							}
						}	
					}
				}
			}
		}
		System.out.println();
		if(residues != naccessResults) {
			System.err.println("Warning: Only " + naccessResults + " NACCESS results founds for " + residues + " in residue table.");
		}
	}
	
	/**
	 * Check for each observed residue whether it coincides with or is proximal to a functional site
	 * as annotated by Uniprot and other databases.
	 * @param conn
	 * @throws SQLException 
	 */
	public void writeProximitiesToDb(MySQLConnection conn) throws SQLException {
		System.out.print("Writing proximities to functional sites to database");
		
		int actHits = 0, actProx = 0;
		int modHits = 0, modProx = 0;
		int atpHits = 0, atpProx = 0;
		int glycHits = 0, glycProx = 0;
		int numMut = 0;
		
		for(Gene g:this.getTargets()) {
			if(!g.areSubstructuresLoaded()) {
				System.err.println("Skipping gene " + g.getGeneName() + ": substructures not loaded.");
			}
			for(Substructure ss:g.getSubstructures()) {
				if(ss.getType() == SubstructureType.PREDICTION) {
					System.out.println("Skipping predicted structure " + g.getGeneName() + "_" + ss.getRange());
					continue;
				}
				if(!ss.isPdbLoaded()) {
					System.err.println("Skipping " + g.getGeneName() + " " + ss.getRange() + " " + ss.getPdbCode()+ss.getChainCode() + ": Structure not loaded.");
				} else {
					System.out.print(".");				
					String pdbCode = ss.getPdbCode();
					String cifChain = ss.pdb.getChainCode();
					String pdbChain = ss.pdb.getPdbChainCode(); // = ss.getChainCode();
					if(!pdbChain.equals(cifChain)) {
						System.err.println("Warning: pdb chain code " + pdbChain + " != cif chain code " + cifChain);
					}
					for(int p: ss.pdb.getAllStdAaResSerials()) {	// only observed residues
						//int pdbPos = Integer.parseInt(ss.pdb.getResidue(p).getPdbSerial());	// cif -> pdb serial
						int cifPos = p;
						int uniprotPos = ss.mapCifResser2Uniprot(p);
						//AminoAcid wtAa = ss.pdb.getResidue(p).getAaType();		
						// check for functional sites
						String fsHit = null;
						String fsNear = null;
						// Remember observed features to prevent overcounting of the same feature
						HashSet<Feature> alreadyObserved = new HashSet<Feature>();
						// check whether this position matches a feature
						for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.UNIPROT, uniprotPos)) {
							alreadyObserved.add(f);
							UniprotFeature uf = (UniprotFeature) f;
							if(uf.getUniprotTypeName().equals("ACT_SITE")) {actHits++; fsHit="ACT_SITE";}
							if(uf.getUniprotTypeName().equals("NP_BIND")) {atpHits++; fsHit="NP_BIND";}
							if(uf.getUniprotTypeName().equals("MOD_RES")) {modHits++; fsHit="MOD_RES";}
							if(uf.getUniprotTypeName().equals("CARBOHYD")) {glycHits++;	fsHit="CARBOHYD";}					
						}	
						if(ss.graph == null) System.err.println("Error: graph not loaded for " + g.geneName + " " + ss.getRange());
						RIGNode mutNode = ss.graph.getNodeFromSerial(cifPos);
						if(mutNode != null) {	// residue observed
							numMut++;
							if(fsHit == null) {	// only check proximal hit if no direct hit found
								for(RIGNode n:ss.graph.getNeighbors(mutNode)) {
									int pdbNbPos = n.getResidueSerial();
									int uniNbPos = ss.mapCifResser2Uniprot(pdbNbPos);
									// check whether this position matches a feature
									for(owl.core.features.Feature f:g.getFeaturesOfTypeForPosition(FeatureType.UNIPROT, uniNbPos)) {
										if(!alreadyObserved.contains(f)) {	// avoid double counting of the same feature
											alreadyObserved.add(f);
											UniprotFeature uf = (UniprotFeature) f;
											if(uf.getUniprotTypeName().equals("ACT_SITE")) {actProx++; fsNear="ACT_SITE";}
											if(uf.getUniprotTypeName().equals("NP_BIND")) {atpProx++; fsNear="NP_BIND";}
											if(uf.getUniprotTypeName().equals("MOD_RES")) {modProx++; fsNear="MOD_RES";}
											if(uf.getUniprotTypeName().equals("CARBOHYD")) {glycProx++; fsNear="CARBOHYD";}
										}
									}						
								}
							}
						}						
						// write hits/proximities to database
						String table = "mutanom.all_observed";
						if(fsHit != null) {
							String sqlHit = "UPDATE " + table + " SET fs_hit='%s' WHERE pdb_code='%s' AND pdb_chain='%s' AND cif_pos=%d";
							conn.executeSql(String.format(sqlHit, fsHit, pdbCode, pdbChain, cifPos));
						}
						if(fsNear != null) {
							String sqlNear = "UPDATE " + table + " SET fs_near='%s' WHERE pdb_code='%s' AND pdb_chain='%s' AND cif_pos=%d";
							conn.executeSql(String.format(sqlNear, fsNear, pdbCode, pdbChain, cifPos));
						}
					}
				}
			}
		}
		System.out.println();
		System.out.printf("actHits = %4d, actProx = %4d, sum = %4d\n",actHits,actProx,actHits+actProx);
		System.out.printf("modHits = %4d, modProx = %4d, sum = %4d\n",modHits,modProx,modHits+modProx);
		System.out.printf("atpHits = %4d, atpProx = %4d, sum = %4d\n",atpHits,atpProx,atpHits+atpProx);
		System.out.printf("glycHits = %3d, glycProx = %3d, sum = %4d\n",glycHits,glycProx,glycHits+glycProx);
		int all = actHits+actProx+modHits+modProx+atpHits+atpProx+glycHits+glycProx;
		double frac = 1.0 * all / numMut;
		System.out.printf("total:    %4d            %4d        %4d\n",actHits+modHits+atpHits+glycHits, actProx+modProx+atpProx+glycProx, all);		
		System.out.printf("numMut = %d\n",numMut);
		System.out.printf("fract  = %6.4f\n",frac);	
	}
	
	/**
	 * Collect the results of a FoldX cluster job submitted with the script created by
	 * writeFoldXMutScanMasterShellScript(). Creates the db table if not exists.
	 * Expects the residue table to be created previously using createAllPositionsDbTable().
	 * @param resTable the residue table which has be created previously
	 * @param outTable the table where the results will be written to (will be created if necessary)
	 * @param resultDir the directory where the results of the cluster job were written to
	 * @throws SQLException 
	 */
	public void collectFoldXResults(MySQLConnection conn, String resTable, String outTable, File resultDir, boolean singleChain) throws SQLException {
		
		boolean useSingleChainPdbFiles = singleChain;	// if false, results for original pdb files will be loaded, otherwise for single chain files
		
		int numFiles = 0;
		int filesNotFound = 0;
		int maxNumMut = 0;
		int minNumMut = Integer.MAX_VALUE;
		int sumNumMut = 0;
		
		// make sure that residue table exists
		int numRes = conn.getIntFromDb("SELECT COUNT(*) FROM " + resTable);
		if(numRes <= 0) {
			System.err.println("Error. Table " + resTable + " not found or is empty.");
			return;
		}
		
		// create output table
		System.out.println("Creating table (if not exists) " + outTable + "...");
		String sqlCreate = "CREATE TABLE IF NOT EXISTS " + outTable + "(res_idx int, aa_mut char(1), ddg_sc float, ddg_cx float)";
		conn.executeSql(sqlCreate);
				
		System.out.print("Writing results to database");
		
		for(Gene g:this.getTargets()) {
			if(!g.areSubstructuresLoaded()) {
				System.err.println("Skipping gene " + g.getGeneName() + ": substructures not loaded.");
			}
			for(Substructure ss:g.getSubstructures()) {
				if(!ss.isPdbLoaded()) {
					System.err.println("Skipping " + g.getGeneName() + " " + ss.getRange() + " " + ss.getPdbCode()+ss.getChainCode() + ": Structure not loaded.");
				} else {
					System.out.print(".");
					String pdbCode = ss.getPdbCode();
					String cifChain = ss.pdb.getChainCode();
					String pdbChain = ss.pdb.getPdbChainCode(); // = ss.getChainCode();
					if(!pdbChain.equals(cifChain)) {
						System.err.println("Warning: pdb chain code " + pdbChain + " != cif chain code " + cifChain);
					}
					for(int p: ss.pdb.getAllStdAaResSerials()) {	// only observed residues
						String pdbPos = "?";
						String cifPos = "?";
						pdbPos = ss.pdb.getResidue(p).getPdbSerial();	// cif -> pdb serial
						cifPos = Integer.toString(p);
						AminoAcid wtAa = ((AaResidue)ss.pdb.getResidue(p)).getAaType();
						String fileName = null;
						if(useSingleChainPdbFiles) {
							fileName = FoldXRunner.getJobName(pdbCode, pdbChain, cifPos) + ".txt";
						} else {
							fileName = FoldXRunner.getJobName(pdbCode, pdbChain, pdbPos) + ".txt";
						}
						try {
							BufferedReader in = new BufferedReader(new FileReader(new File(resultDir, fileName)));
							String line;
							int numMut = 0;
							while((line = in.readLine()) != null) {
								String[] fields = line.split("\t");
								if(fields.length < 23) {
									System.err.println("Error parsing " + fileName + " line " + numMut + ":" + line);
									System.exit(1);
								}
								String mutStr = fields[0].trim();
								double ddg = Double.parseDouble(fields[1]);
								String basename = (mutStr.split("\\."))[0];
								int mutNum = Integer.parseInt((basename.split("_"))[1]);	// e.g. mut_12.pdb -> 12 
								if(wtAa.getNumber() == mutNum) {
									if(Math.abs(ddg) > 2.0) {
										System.err.printf("Warning: Energy difference for WT (%s) should be near zero but is %f in %s:\n", wtAa, ddg, fileName);
										System.err.println(line);
										//System.exit(1);
									}
								} else {
									numMut++;
									AminoAcid mutAa = AminoAcid.getByNumber(mutNum);
									try {
										writeMutScanDbEntry(conn, resTable, outTable, pdbCode, cifChain, Integer.parseInt(cifPos), pdbChain, Integer.parseInt(pdbPos), wtAa, mutAa, ddg, useSingleChainPdbFiles);
									} catch (SQLException e) {
										System.err.println("Error writing to database: " + e.getMessage());
									}
								}
							}
							numFiles++;
							sumNumMut += numMut;
							minNumMut = Math.min(minNumMut, numMut);
							maxNumMut = Math.max(maxNumMut, numMut);
							in.close();
						} catch (IOException e) {
							System.err.println("Could not read from file " + fileName + ": " + e.getMessage());
							filesNotFound++;
						}
					}			
				}
			}
		}
		System.out.println();
		System.out.println("Total residues:   " + numRes);
		System.out.println("Results read:     " + numFiles);
		System.out.println("Files not found:  " + filesNotFound);
		System.out.println("Total Mutations:  " + sumNumMut + " (should be " + numFiles + " * 19 = " + (numFiles * 19) + ")");
		System.out.println("Min Num Mutations: " + minNumMut);
		System.out.println("Max Num Mutations: " + maxNumMut);
		System.out.println("Avg Num Mutations: " + 1.0 * sumNumMut / numFiles);
	}
	
	public static void writeMutScanDbEntry(MySQLConnection conn, String resTable, String outTable, String pdbCode, String cifChain, int cifPos, String pdbChain, int pdbPos, AminoAcid aaWt, AminoAcid aaMut, double deltaDeltaG, boolean singleChain) throws SQLException {
		String sql;
		
		String query = "SELECT idx FROM " + resTable + " WHERE pdb_code='%s' AND pdb_chain='%s' AND cif_pos=%d";
		int resIdx = conn.getIntFromDb(String.format(query, pdbCode, pdbChain, cifPos));
		if(resIdx < 0) {
			String err = "Could not find entry for pdb_code='%s', pdb_chain='%s', cif_pos=%d in table " + resTable + ". Skipping.";
			System.err.println(String.format(err, pdbCode, pdbChain, cifPos));
			System.err.println(String.format(query, pdbCode, pdbChain, cifPos));
		} else {
			if(singleChain) {
				sql = "INSERT INTO " + outTable + "(res_idx, aa_mut, ddg_sc) VALUES (%d, '%s', %f)";
			} else {
				sql = "INSERT INTO " + outTable + "(res_idx, aa_mut, ddg_cx) VALUES (%d, '%s', %f)";			
			}
			conn.executeSql(String.format(sql, resIdx, Character.toString(aaMut.getOneLetterCode()), deltaDeltaG));
		}
	}
	
	public static void writeMutScanDbEntry_old(MySQLConnection conn, String dbTable, String pdbCode, String cifChain, int cifPos, String pdbChain, int pdbPos, AminoAcid aaWt, AminoAcid aaMut, double deltaDeltaG, boolean singleChain) throws SQLException {
		String sql;
		if(singleChain) {
			//sql = String.format("INSERT INTO %s(pdb_code, cif_chain, cif_pos, pdb_chain, pdb_pos, aa_wt, aa_mut, ddg_sc) VALUES('%s','%s',%d,'%s',%d,'%s','%s',%f)",
			//		dbTable, pdbCode, cifChain, cifPos, pdbChain, pdbPos, aaWt.getOneLetterCode(), aaMut.getOneLetterCode(), deltaDeltaG);
			sql = String.format("UPDATE %s SET ddg_sc=%f WHERE pdb_code='%s' AND cif_chain='%s' AND cif_pos=%d AND aa_wt='%s' AND aa_mut='%s'",
					dbTable, deltaDeltaG, pdbCode, cifChain, cifPos, aaWt.getOneLetterCode(), aaMut.getOneLetterCode());		
		} else {			
			//sql = String.format("UPDATE %s SET ddg_cx=%f WHERE pdb_code='%s' AND cif_chain='%s' AND cif_pos=%d AND aa_wt='%s' AND aa_mut='%s'",
			//		dbTable, deltaDeltaG, pdbCode, cifChain, cifPos, aaWt.getOneLetterCode(), aaMut.getOneLetterCode());
						
			sql = String.format("INSERT INTO %s(pdb_code, cif_chain, cif_pos, pdb_chain, pdb_pos, aa_wt, aa_mut, ddg_cx) VALUES('%s','%s',%d,'%s',%d,'%s','%s',%f)",
					dbTable, pdbCode, cifChain, cifPos, pdbChain, pdbPos, aaWt.getOneLetterCode(), aaMut.getOneLetterCode(), deltaDeltaG);
		}
		conn.executeSql(sql);
	}
	
	/**
	 * Collects SNP data for this target list from dbSNP xml files. Creates the destination
	 * table in the database (if not exists). Only parses SNPs whose NM (base) number matches the gene.
	 * @param conn the database connection for writing the results and querying nm numbers
	 * @param snpTable the database table where the results will be written to
	 * @param xmlDir xml files with SNPs are expected in this directory
	 * @throws SQLException 
	 */
	public void collectNcbiSnpsFromXmlFiles(MySQLConnection conn, String snpTable, File xmlDir) throws SQLException {

		String createSql= "CREATE TABLE IF NOT EXISTS " + snpTable + "(dbsnp_id varchar(15), gene_name varchar(30), mut_dna varchar(50), mut_aa varchar(3), pos_aa int, freq float)";
		conn.executeSql(createSql);
		
		DbSnpConnection dbSnp = new DbSnpConnection();
		
		for(Gene g:this.targets) {
			String nm = conn.getStringFromDb(String.format(QUERY_GENE2NM, g.getGeneName()));
			String np = conn.getStringFromDb(String.format(QUERY_GENE2NP, g.getGeneName()));			
			System.out.print(g.getGeneName() + "\t" + nm + "\t" + np + "\t");
			File xmlFile = new File(xmlDir, "s_" + g.getGeneName() + ".xml");
			if(!xmlFile.canRead()) {
				System.out.println("File not found: " + xmlFile.getPath());
			} else {
				Collection<SNP> snps = dbSnp.parseSnpsFromXmlFile(xmlFile, np, nm);
				System.out.println(snps.size() + " snps found.");
				for(SNP snp:snps) {
					String mutStr = snp.getWtAA().getOneLetterCode() + "/" + snp.getMutAA().getOneLetterCode();
					conn.executeSql(String.format("INSERT INTO " + snpTable + "(dbsnp_id, gene_name, mut_dna, mut_aa, pos_aa) VALUES('%s','%s','%s','%s',%d)", "rs" + snp.rsId, g.getGeneName(), snp.getDnaMutStr(), mutStr, snp.getPos()));
				}		
			}
		}
	}
		
	/**
	 * Creates a table with all observed residues in all loaded substructrues.
	 * @param conn an active database connection
	 * @param table the database and table name
	 * @throws SQLException 
	 */
	public void createAllPositionsDbTable(MySQLConnection conn, String table) throws SQLException {
		String createSql= "CREATE TABLE IF NOT EXISTS " + table + "(idx int primary key auto_increment, pdb_code char(4), type char(1), cif_chain char(1), cif_pos int, pdb_chain char(1), pdb_pos int, aa_wt char(1))";
		conn.executeSql(createSql);
		String sql = "INSERT INTO " + table + "(pdb_code, type, cif_chain, cif_pos, pdb_chain, pdb_pos, aa_wt) VALUES('%s','%s','%s',%d,'%s',%d,'%s')";
		for(Gene g:this.targets) {
			for(Substructure ss:g.getSubstructures()) {
				String pdbCode = ss.pdb.getPdbCode();
				String cifChain = ss.pdb.getChainCode();
				String pdbChain = ss.pdb.getPdbChainCode();
				for(int cifPos:ss.pdb.getAllStdAaResSerials()) {
					AminoAcid wtAa = ((AaResidue)ss.pdb.getResidue(cifPos)).getAaType();
					int pdbPos = Integer.parseInt(ss.pdb.getResidue(cifPos).getPdbSerial());
					String type = Character.toString(ss.getTypeChar());
					conn.executeSql(String.format(sql, pdbCode, type, cifChain,cifPos,pdbChain,pdbPos,wtAa.getOneLetterCode()));
				}
			}
		}
	}
	
	/**
	 * Simply print out a list with the lengths of all substructures in this list.
	 */
	public void printTargetsLengths() {
		for(Gene g:this.targets) {
			for(Substructure ss:g.getSubstructures()) {
				System.out.println(ss.end-ss.start+1);
			}
		}
	}
	
	/*---------------------------- static methods ---------------------------*/
	
	/**
	 * @param online if true, load sequences from online Uniprot/COSMIC, otherwise from local files defined in Gene.*_SEQ_PATH
	 */
	public static TargetList get2PredictedStructures(boolean online) {
		TargetList newTargetList = new TargetList();
		
//		Gene mlh1 = new Gene("MLH1","P40692");
//		mlh1.setLength(756);
//		mlh1.loadUniprotSequence(online);
//		if(mlh1.getLength() != mlh1.getUniprotSeq().length()) {
//			System.err.println("Warning: Gene length does not match length of uniprot sequence");
//		}
//		newTargetList.addTarget(mlh1);
		
		Gene erbb2 = new Gene("ERBB2","P04626");
		erbb2.setLength(1255);
		erbb2.loadUniprotSequence(online);
		if(erbb2.getLength() != erbb2.getUniprotSeq().length()) {
			System.err.println("Warning: Gene length does not match length of uniprot sequence");
		}
		newTargetList.addTarget(erbb2);
		
		return newTargetList;
	}
	
	
	/**
	 * Returns a target list with the modelled targets from db table 'targets_predicted'
	 * @param online whether to load sequence from online Uniprot
	 * @return
	 */
	public static TargetList loadPredictedTargets(MySQLConnection conn, boolean online) {
		TargetList newTargetList = new TargetList();
		
		// get genes from database
		String[] genes = conn.getStringsFromDb(QUERY_LOAD_PREDICTED);
		if(genes == null || genes.length == 0) {
			System.err.println("Error. Could not load target list.");
			System.err.println(QUERY_LOAD_PREDICTED);
		} else {
			for(String geneName:genes) {
				String uniprotId = conn.getStringFromDb(String.format(QUERY_GENE2SP, geneName));
				Gene g = new Gene(geneName, uniprotId);
				g.loadUniprotSequence(online);	
				g.loadCosmicSequence(online);
				//int u = g.getUniprotSeq().length();
				//int c = g.getCosmicSeq().length();
				// TODO: Move the following line to a seperate data check function
				g.setLength(g.getCosmicSeq().length());
				newTargetList.addTarget(g);
				//System.out.print(".");
			}
		}		
		return newTargetList;
	}
	
// OBSOLETE:	
//	/**
//	 * Loads the 11 targets which are included in the results for the paper
//	 * @param online if true, load sequences from online Uniprot/COSMIC, otherwise from local files defined in Gene.*_SEQ_PATH
//	 * @return the target list with the 11 targets
//	 */
//	public static TargetList getTargetListForPublication(boolean online) {
//		TargetList newTargetList = new TargetList();
//		try {
//			BufferedReader in = new BufferedReader(new FileReader(Gene.TARGET_FILE_11));
//			String line = in.readLine(); // skip header line
//			while((line=in.readLine()) != null) {
//				String[] fields = line.split("\t");
//				String name = fields[0];
//				String uniProtId = fields[2];
//				int length = Integer.parseInt(fields[3]);
//				Gene gene = new Gene(name, uniProtId);
//				gene.setLength(length);
//				gene.loadUniprotSequence(online);
//				if(gene.getLength() != gene.getUniprotSeq().length()) {
//					System.err.println("Warning: Gene length does not match length of uniprot sequence");
//				}
//				newTargetList.addTarget(gene);
//			}
//			in.close();
//		} catch (FileNotFoundException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		return newTargetList;		
//	}

// OBSOLETE
//	/**
//	 * Old targets list with which results for the Oncogene version of the paper were produced
//	 * @param online if true, load sequences from online Uniprot/COSMIC, otherwise from local files defined in Gene.*_SEQ_PATH
//	 * @return
//	 */
//	public static TargetList getTop20TargetList(boolean online) {
//		TargetList newTargetList = new TargetList();
//		try {
//			BufferedReader in = new BufferedReader(new FileReader(Gene.TARGET_FILE));
//			String line = in.readLine(); // skip header line
//			while((line=in.readLine()) != null) {
//				String[] fields = line.split("\t");
//				String name = fields[0];
//				String uniProtId = fields[2];
//				int length = Integer.parseInt(fields[3]);
//				Gene gene = new Gene(name, uniProtId);
//				gene.setLength(length);
//				gene.loadUniprotSequence(online);
//				if(gene.getLength() != gene.getUniprotSeq().length()) {
//					System.err.println("Warning: Gene length does not match length of uniprot sequence");
//				}
//				newTargetList.addTarget(gene);
//			}
//			in.close();
//		} catch (FileNotFoundException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		return newTargetList;
//	}
	
	/**
	 * Top29 target list for selecting genes for the Cancer Research version of the paper.
	 * @param online if true, load sequences from online Uniprot/COSMIC, otherwise from local files defined in Gene.*_SEQ_PATH
	 * @return
	 */
	public static TargetList loadTop29TargetList(MySQLConnection conn, boolean online) {
		TargetList newTargetList = new TargetList();
		
		// get genes from database
		String sql = "SELECT * FROM mutanom2.targets_top29";
		String[] genes = conn.getStringsFromDb(sql);
		if(genes == null || genes.length == 0) {
			System.err.println("Error. Could not load top29 target list.");
			System.err.println(sql);
		} else {
			for(String geneName:genes) {
				String uniprotId = conn.getStringFromDb("SELECT sp_id FROM mutanom2.gene2sp WHERE gene_name='" + geneName + "';");
				Gene g = new Gene(geneName, uniprotId);
				g.loadUniprotSequence(online);	
				g.loadCosmicSequence(online);
				int u = g.getUniprotSeq().length();
				int c = g.getCosmicSeq().length();
				// TODO: Move the following line to a seperate data check function
				if(u != c) System.err.println(String.format("Warning: %s cosmic len (%d) != uniprot len (%d)", geneName, c, u));
				g.setLength(g.getCosmicSeq().length());
				newTargetList.addTarget(g);
			}
		}		
		return newTargetList;
	}
	
	/**
	 * Top18 target list for selecting genes for the Cancer Research version of the paper.
	 * @param online if true, load sequences from online Uniprot/COSMIC, otherwise from local files defined in Gene.*_SEQ_PATH
	 * @return
	 */
	public static TargetList loadTop18TargetList(MySQLConnection conn, boolean online) {
		TargetList newTargetList = new TargetList();
		
		// get genes from database
		String sql = "SELECT * FROM mutanom2.targets_top18";
		String[] genes = conn.getStringsFromDb(sql);
		if(genes == null || genes.length == 0) {
			System.err.println("Error. Could not load top18 target list.");
			System.err.println(sql);
		} else {
			for(String geneName:genes) {
				String uniprotId = conn.getStringFromDb("SELECT sp_id FROM mutanom2.gene2sp WHERE gene_name='" + geneName + "';");
				Gene g = new Gene(geneName, uniprotId);
				g.loadUniprotSequence(online);	
				g.loadCosmicSequence(online);
				int u = g.getUniprotSeq().length();
				int c = g.getCosmicSeq().length();
				// TODO: Move the following line to a seperate data check function
				if(u != c) System.err.println(String.format("Warning: %s cosmic len (%d) != uniprot len (%d)", geneName, c, u));
				g.setLength(g.getCosmicSeq().length());
				newTargetList.addTarget(g);
			}
		}		
		return newTargetList;
	}
	
	/**
	 * Load targets from table 'targets' in the database.
	 * @param online if true, load sequences from online Uniprot/COSMIC, otherwise from local files defined in Gene.*_SEQ_PATH
	 * @return
	 */
	public static TargetList loadTargets(MySQLConnection conn, boolean online) {
		TargetList newTargetList = new TargetList();
		
		// get genes from database
		String[] genes = conn.getStringsFromDb(QUERY_LOAD_TARGETS);
		if(genes == null || genes.length == 0) {
			System.err.println("Error. Could not load target list.");
			System.err.println(QUERY_LOAD_TARGETS);
		} else {
			for(String geneName:genes) {
				String uniprotId = conn.getStringFromDb(String.format(QUERY_GENE2SP, geneName));
				Gene g = new Gene(geneName, uniprotId);
				g.loadUniprotSequence(online);	
				g.loadCosmicSequence(online);
				//int u = g.getUniprotSeq().length();
				//int c = g.getCosmicSeq().length();
				// TODO: Move the following line to a seperate data check function
				g.setLength(g.getCosmicSeq().length());
				newTargetList.addTarget(g);
				//System.out.print(".");
			}
		}		
		return newTargetList;
	}
	
	/**
	 * Get a target list consisting of a single gene by its gene name.
	 * @param conn
	 * @param online
	 * @param geneName
	 * @return
	 */
	public static TargetList loadSingletonList(MySQLConnection conn, boolean online, String geneName) {
		TargetList newTargetList = new TargetList();
		String uniprotId = conn.getStringFromDb("SELECT sp_id FROM mutanom2.gene2sp WHERE gene_name='" + geneName + "';");
		Gene g = new Gene(geneName, uniprotId);
		g.loadUniprotSequence(online);	
		g.loadCosmicSequence(online);
		int u = g.getUniprotSeq().length();
		int c = g.getCosmicSeq().length();
		// TODO: Move the following line to a seperate data check function
		if(u != c) System.err.println(String.format("Warning: %s cosmic len (%d) != uniprot len (%d)", geneName, c, u));
		g.setLength(g.getCosmicSeq().length());
		newTargetList.addTarget(g);
		return newTargetList;
	}
	
	/**
	 * Connects to database and returns connection object or null if connection failed.
	 * Settings for user name, password, and server are taken from ~/.my.cnf
	 * @return the connection object or null
	 */
	public static MySQLConnection connectToDb()	{
		MySQLConnection conn = null;
		try {
			conn = new MySQLConnection();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return conn;
	}
	
	/*--------------------------------- main --------------------------------*/
	public static void main(String[] args) throws SQLException {
		
		// supress output of JAligner
//		File trashLogFile = null;
//		try {
//			trashLogFile = File.createTempFile("logger", ".trash");
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		trashLogFile.deleteOnExit();
//		System.setProperty("java.util.logging.config.file",trashLogFile.getAbsolutePath());
		
		MySQLConnection conn = TargetList.connectToDb();
		
		// load data
		boolean online = false;
		TargetList tl = TargetList.loadTop29TargetList(conn, online);
		tl.loadSubstructures(conn, false);	// not only predicted
		//tl.loadMutations();
		tl.loadMutationsFromDb(conn);		// all cosmic mutations
		//tl.loadSnpsFromDb(conn);
		
		// load sequences
		//File cosmicFile = new File("/project/PyMol/targets/seq/top19_cosmic.fa");
		//File uniprotFile = new File("/project/PyMol/targets/seq/top19_uniprot.fa");
		//tl.loadUniprotSequences(uniprotFile);
		//tl.loadCosmicSequences(cosmicFile);
		//tl.compareUniprotVsCosmicSequence();

		// download structures
		//tl.report();
		//File outDir = new File("/project/PyMol/test");
		//tl.downloadPdbsForSubstructures(outDir);
		
		// output data
		//File outDir = new File("/home/web/lappe/stehr/mutanom/v4/gene");
		//File outDir = new File("/project/StruPPi/henning/projects/mutanom/meeting2");
		//tl.renderOverview(outDir);
		//tl.renderRulers(outDir);
		//tl.writeGeneDataFiles(outDir);
		//File pseDir = new File("/home/web/lappe/stehr/mutanom/pse/");
		//File pdbDir = new File("/home/web/lappe/stehr/mutanom/pdb/");
		//File pseDir = new File("/project/StruPPi/henning/projects/mutanom/meeting2");
		//File pdbDir = new File("/project/StruPPi/henning/projects/mutanom/meeting2");
		//tl.generatePngsAndPyMolSessions(pseDir, pdbDir);
		//tl.generatePseWithAllMutations(pseDir, pdbDir);
		//tl.report();
		//tl.checkPdbFiles();
		//tl.printPdbCodes(false);
		//tl.printPdbCodes(true);
		//tl.loadPdbsAndAlignments(new File("/project/PyMol/pdbs_download/pred_not_renum/"));	// pdbDir only for predicted pdbs
		tl.printTargetsLengths();
		System.out.println("Done.");
	}

	public void writeSeqOverviewForPrint(File outFile) {
        
		SvgGenerator.writeSeqOverviewForPrint(outFile, this);

		
	}	

}
