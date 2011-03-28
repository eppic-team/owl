package owl.mutanom;

import java.io.IOException;
import java.sql.SQLException;
import java.util.Collection;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


import owl.core.connections.NoMatchFoundException;
import owl.core.connections.SiftsConnection;
import owl.core.connections.UniProtConnection;
import owl.core.connections.UniProtPdbRef;
import owl.core.features.SiftsFeature;
import owl.core.util.Interval;
import owl.core.util.MySQLConnection;
import owl.mutanom.core.Gene;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;

/**
 * Creates the database tables 'gene2sp', 'genes' and 'pdbs' from the COSMIC table
 * and some online resources (Uniprot, Cosmic, SIFTS, local Pdbase).
 * 
 * Input: table with mutations downloaded from COSMIC
 * Output: Database tables 'gene2sp', 'genes', 'pdbs'
 * Dependencies:
 * - Cosmic website  (for obtaining sequences)
 * - online Uniprot  (for sequences, pdb x-refs, updated ids
 * - online SIFTS    (chain2uniprot from ftp, for pdb x-refs) 
 * - local pdbase db (for querying experimental method for SIFTS entries) TODO: query online
 * 
 * Detailed explanation:
 * This class contains tools for preselecting target genes from the COSMIC database.
 * Method are provided for creating the basic tables in the database. The main result
 * is the table 'genes' which contains upper bounds for the number of mutations in
 * structural regions. Based in this we can make a pre-selection of genes for which
 * we will then manually select substructures and put them into table 'substructures'.
 * Automatic selection would be too troublesome to implement in a robust way for the
 * moment. Based on the selected substructures, we can then calculate the actual
 * numbers of mutations and SNPs in structural (and observed) regions as well as the
 * coverage of the sequence with structural regions which we will beed to select the
 * genes finally used in the analysis. See: GeneInfo.java
 * 
 * Note: This class should use online resources as much as possible to get the latest
 * data and to be self-contained. Subsequent analysis steps should then use only
 * downloaded data to be reproducable.
 * 
 * @author stehr
 */
public class GenePreSelectionTools {
	
	/*------------------------------ constants ------------------------------*/
	// old
//	public static final String QUERY_GENES_43 = "select gene_name, count(distinct mut_aa) as cnt from mutanom.cosmic43_mis_4tis group by gene_name having cnt > 2 order by cnt desc;";
//	public static final String QUERY_MUTATIONS_43_4 = "SELECT DISTINCT mut_aa FROM mutanom.cosmic43_mis_4tis WHERE gene_name=\"%s\"";
//	public static final String QUERY_MUTATIONS_43_5 = "SELECT DISTINCT mut_aa FROM mutanom.cosmic43_mis_5tis WHERE gene_name=\"%s\"";
//
//	public static final String QUERY_GENES_44 = "select gene_name, count(distinct mut_aa) as cnt from mutanom.cosmic44_mis_5tis group by gene_name having cnt > 2 order by cnt desc;";
//	public static final String QUERY_GENES_44_ALL = "select gene_name, count(distinct mut_aa) as cnt from mutanom.cosmic44_mis_5tis group by gene_name order by cnt desc;";
//	public static final String QUERY_MUTATIONS_44_4 = "SELECT DISTINCT mut_aa FROM mutanom.cosmic44_mis_4tis WHERE gene_name=\"%s\"";
//	public static final String QUERY_MUTATIONS_44_5 = "SELECT DISTINCT mut_aa FROM mutanom.cosmic44_mis_5tis WHERE gene_name=\"%s\"";
//
//	public static final String INSERT_GENE2SP_NULL = "INSERT INTO mutanom.gene2sp VALUES(\"%s\",null);";
//	public static final String QUERY_GENE2SP_44 = "SELECT sp_id FROM mutanom.cosmic44_gene2sp WHERE gene=\"%s\";";	
//	public static final String INSERT_GENE2SP_44 = "INSERT INTO mutanom.cosmic44_gene2sp VALUES(\"%s\",\"%s\");";
//	public static final String INSERT_GENE2SP_NULL_44 = "INSERT INTO mutanom.cosmic44_gene2sp VALUES(\"%s\",null);";
	
	// database
	public static final String DATABASE = 		"mutanom3";
	public static final String COSMIC_TABLE = 	DATABASE + ".mut_cosmic49_mis_8tis";	// cosmic data
	public static final String GENE2SP_TABLE = 	DATABASE + ".gene2sp";					// gene -> uniprot mapping
	public static final String GENES_TABLE = 	DATABASE + ".genes";					// mutation statistics for gene pre-selection
	public static final String PDBS_TABLE = 	DATABASE + ".pdbs"; 					// uniprot -> pdb mapping
	
	public static final String QUERY_ALL_GENES = "SELECT gene_name, count(distinct mut_dna) as cnt_mut_dna, count(distinct mut_aa) as cnt_mut_aa FROM " + COSMIC_TABLE + " group by gene_name having cnt_mut_aa >= %d order by cnt_mut_aa desc;";	
	public static final String QUERY_MUTATIONS = "SELECT DISTINCT CONCAT(mut_dna, \";\", mut_aa) FROM " + COSMIC_TABLE + " WHERE gene_name=\"%s\"";

	public static final String CREATE_GENE2SP_TABLE = "CREATE TABLE IF NOT EXISTS " + GENE2SP_TABLE + "(gene_name varchar(30), sp_id char(6), source varchar(20));";
	public static final String INSERT_GENE2SP = "INSERT INTO " + GENE2SP_TABLE + " VALUES(\"%s\",\"%s\", \"%s\");";	
	public static final String QUERY_GENE2SP = "SELECT sp_id FROM " + GENE2SP_TABLE + " WHERE gene_name=\"%s\";";	

	public static final String CREATE_GENES_TABLE = "CREATE TABLE IF NOT EXISTS " + GENES_TABLE + "(gene_name varchar(30), sp_id char(6), sp_len int, cos_len int, distinct_dna_mut int, distinct_mut int, struc_mut int, struc_mut_pos int, str_res int, struc_regions varchar(255));";
	public static final String INSERT_GENES = "INSERT INTO " + GENES_TABLE + "(gene_name, sp_id, sp_len, cos_len, distinct_dna_mut, distinct_mut , struc_mut, struc_mut_pos, str_res, struc_regions) VALUES(\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,\"%s\")";
	public static final String QUERY_CANDIDATE_GENES = "SELECT gene_name FROM " + GENES_TABLE + " WHERE struc_mut >= %d";
	
	public static final String PDBASE_GET_KEY = "SELECT entry_key FROM pdbase.mms_entry WHERE id='%s'";
	public static final String PDBASE_GET_METHOD = "SELECT substring_index(method,',',1) FROM pdbase.exptl WHERE entry_key=%d"; 
	public static final String PDBASE_GET_RESOL = "SELECT ls_d_res_high FROM pdbase.refine WHERE entry_key=%d";
	
	// parameters (for estimation of the number of structural mutations)
	public static final boolean XRAY_ONLY = 		true;		// if true, only consider xray structures (otherwise also NMR)
	public static final double 	MAX_RESOLUTION = 	5.0;		// maximum resolution for xray structures to be considered
	public static final int 	MIN_PDB_LENGTH = 	30;			// minimum length for structures to be considered
	
	// parameters for obtaining PDBs for given Uniprot entries
	public static final boolean USE_SIFTS = 		true;		// if true, get PDB info from SIFTS (recommended), otherwise from UniProt
	
	// parameters for both
	public static final int 	MIN_NUM_MUT = 		5; 			// genes with less distinct aa mutations will be ignored
	
	// aliases
	public static final boolean ONLINE = true;		// whether to use online databases for loading Uniprot and COSMIC sequence
	public static final boolean OFFLINE = false;	// whether to use online databases for loading Uniprot and COSMIC sequence
	
	/**
	 * Query table cosmic_mis_5tis and for all genes with more than 3 distinct aa mutations get the UniprotID wither from COSMIC or
	 * by querying Uniprot with the Cosmic sequence. Write the results to database table gene2sp.
	 * @param conn an active database connection
	 * @throws SQLException on database error
	 */
	public static void writeGene2SpTable(MySQLConnection conn) throws SQLException {
		
		createGene2SpTable(conn);	// creates the table for mapping gene names to Uniprot IDs (if not exists)
		
		String[] geneNames = conn.getStringsFromDb(String.format(QUERY_ALL_GENES, MIN_NUM_MUT));
		UniProtConnection upc = new UniProtConnection();
		for(String geneName:geneNames) {
			System.out.println(geneName);
			String sql = null;
			String source = null;
			Gene g = new Gene(geneName);
			g.loadCosmicSequence(ONLINE);	// use online Uniprot
			g.loadUniprotIdFromCosmic();	// uses online COSMIC
			if(g.getUniprotId() == null) {
				String uniprotId = upc.getHumanEntryBySequence(g.getCosmicSeq()).getPrimaryUniProtAccession().getValue();
				g.setUniprotId(uniprotId);
				System.out.println(uniprotId);
				source = "blast";
			} else {
				source = "cosmic";
			}				
			sql = String.format(INSERT_GENE2SP, geneName, g.getUniprotId(), source);
			conn.executeSql(sql);
		}
	}
	
	/**
	 * Report for each gene the number of mutations in known regions. This is the cutoff we use for selecting genes for analysis.
	 * Note that this does not consider unobserved residues, so it is rather an upper bound on the number of mutations.
	 * @param conn
	 * @throws SQLException
	 */
	public static void writeGenesTable(MySQLConnection conn) throws SQLException {
		
		createGenesTable(conn);	// creates the table 'genes' which contains the number of cosmic mutations per gene (if not exists)
		
		//int minNumStrMut = 7;			// minimum number of mutations in structural regions for a gene to be selected 

		Pattern p = Pattern.compile("\\d+");	// position information in mutation string
		
		// load candidate genes
		String[] geneNames = conn.getStringsFromDb(String.format(QUERY_ALL_GENES, MIN_NUM_MUT));
		System.out.printf("Gene\tUniprot\tlenSp\tlenCos\tnDnaMut\tnumMut\tstrMut\tmutPos\tstrRes\tKnown regions\n");
		for(String geneName:geneNames) {
			int numStrMut = 0;
			Gene g = new Gene(geneName);			
			g.setUniprotId(conn.getStringFromDb(String.format(QUERY_GENE2SP, geneName)));
			if(g.getUniprotId() != null) {
				// load sequences (to get length information)
				g.loadCosmicSequence(ONLINE);
				g.loadUniprotSequence(ONLINE);
				// query uniprot for (xray structure) cross references
				TreeSet<Integer> knownResidues = g.getKnownRegionsFromUniprot(XRAY_ONLY, MAX_RESOLUTION, MIN_PDB_LENGTH);
				// count mutations within known structures
				TreeSet<Integer> strMutPos = new TreeSet<Integer>();
				TreeSet<String> distinctAaMut = new TreeSet<String>();
				String[] mutStrings = conn.getStringsFromDb(String.format(QUERY_MUTATIONS, g.getGeneName()));
				if(mutStrings != null) {
					for(String mutStr:mutStrings) {
						String[] fields = mutStr.split(";");
						//String dnaMut = fields[0];
						String aaMut = fields[1];
						distinctAaMut.add(aaMut);
						Matcher m = p.matcher(aaMut);
						if(m.find()) {
							int pos = Integer.parseInt(m.group());
							if(knownResidues.contains(pos)) {
								strMutPos.add(pos);
								numStrMut++;
							}
						} else {
							System.err.println("Could not parse mutation string " + aaMut);
						}
					}
				}
				// report results
				//if(numStrMut >= minNumStrMut) {
				conn.executeSql(String.format(INSERT_GENES,                   g.getGeneName(), g.getUniprotId(), g.getUniprotSeq().length(), g.getCosmicSeq().length(), mutStrings.length, distinctAaMut.size(), numStrMut, strMutPos.size(), knownResidues.size(), Interval.createSelectionString(knownResidues)));
				System.out.printf("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", g.getGeneName(), g.getUniprotId(), g.getUniprotSeq().length(), g.getCosmicSeq().length(), mutStrings.length, distinctAaMut.size(), numStrMut, strMutPos.size(), knownResidues.size(), Interval.createSelectionString(knownResidues));
				//}
			} else {
				System.out.printf("%s\tno UniprotID found\n", g.getGeneName());
			}		
		}
	}
	
	
	
//	/**
//	 * Report for each gene the number of distinct mutations in known regions. This is the cutoff we use for selecting genes for analysis.
//	 * @param conn
//	 * @throws SQLException
//	 */
//	public static void countMutations_old(MySQLConnection conn) throws SQLException {
//		
//		boolean xrayOnly = true;		// if true, only consider xray structures (otherwise also NMR)
//		double maxResolution = 5.0;		// maximum resolution for xray structures to be considered
//		int minLength = 30;				// minimum length for structures to be considered
//		
//		Pattern p = Pattern.compile("\\d+");	// position information in mutation string
//		
//		System.out.println("Number of mutations in known structures (4 tissues / 5tissues):");
//		System.out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Gene", "UnipID", "Cos43_4", "Cos43_5", "Cos44", "Cos44_5", "Known regions");
//		// load potential genes from database (those with >2 missense mutations)
//		String[] geneNames = conn.getStringsFromDb(QUERY_ALL_GENES);		
//		//String geneName = "TP53";
//		for(String geneName:geneNames) {
//			int mut4tissues = 0;
//			int mut5tissues = 0;
//			int mut44_4tissues = 0;
//			int mut44_5tissues = 0;
//			// get UniprotID
//			Gene g = new Gene(geneName);			
//			g.setUniprotId(conn.getStringFromDb(String.format(QUERY_GENE2SP_44, geneName)));
//			//g.loadUniprotIdFromCosmic();
//			if(g.getUniprotId() != null) {			
//				// query uniprot for (xray structure) cross references
//				TreeSet<Integer> knownResidues = g.getKnownRegionsFromUniprot(xrayOnly, maxResolution, minLength);
//				// count mutations within known structures
//				String[] mutStrings = conn.getStringsFromDb(String.format(QUERY_MUTATIONS_43_4, g.getGeneName()));
//				if(mutStrings != null) {
//					for(String mutStr:mutStrings) {
//						Matcher m = p.matcher(mutStr);
//						if(m.find()) {
//							int pos = Integer.parseInt(m.group());
//							if(knownResidues.contains(pos)) {
//								mut4tissues++;
//							}
//						}
//					}
//				}
//				String[] mutStrings2 = conn.getStringsFromDb(String.format(QUERY_MUTATIONS_43_5, g.getGeneName()));
//				if(mutStrings != null) {
//				for(String mutStr:mutStrings2) {
//					Matcher m = p.matcher(mutStr);
//					if(m.find()) {
//						int pos = Integer.parseInt(m.group());
//						if(knownResidues.contains(pos)) {
//							mut5tissues++;
//						}
//					}
//				}
//				}
//				String[] mutStrings3 = conn.getStringsFromDb(String.format(QUERY_MUTATIONS_44_4, g.getGeneName()));
//				if(mutStrings != null) {
//				for(String mutStr:mutStrings3) {
//					Matcher m = p.matcher(mutStr);
//					if(m.find()) {
//						int pos = Integer.parseInt(m.group());
//						if(knownResidues.contains(pos)) {
//							mut44_4tissues++;
//						}
//					}
//				}
//				}
//				String[] mutStrings4 = conn.getStringsFromDb(String.format(QUERY_MUTATIONS_44_5, g.getGeneName()));
//				if(mutStrings != null) {
//				for(String mutStr:mutStrings4) {
//					Matcher m = p.matcher(mutStr);
//					if(m.find()) {
//						int pos = Integer.parseInt(m.group());
//						if(knownResidues.contains(pos)) {
//							mut44_5tissues++;
//						}
//					}
//				}
//				}
//				// print results
//				if(mut4tissues > 0 || mut5tissues > 0 || mut44_4tissues > 0 || mut44_5tissues > 0) 
//				System.out.printf("%s\t%s\t%d\t%d\t%d\t%d\t%s\n", g.getGeneName(), g.getUniprotId(), mut4tissues, mut5tissues, mut44_4tissues, mut44_5tissues, Interval.createSelectionString(knownResidues));
//			} else {
//				System.out.printf("%s\tno UniprotID found\n", g.getGeneName());
//			}
//			//break;
//		}
//	}
	
	/**
	 * Loads information about associated PDB structures for each gene and writes them to database table 'pdbs'.
	 * This table lists all structures referenced in Uniprot, the actual substructures for the analysis have to be
	 * selected and inserted into table 'substructures'.
	 * @param conn an active database connection
	 * @param table the database and table name, if null output to screen
	 * @throws SQLException if an error occurs when loading from database
	 * @throws IOException if an error occurs when loading from SIFTS
	 */
	public static void writePdbsTable(MySQLConnection conn) throws SQLException, IOException {
		
		UniProtPdbRef.createDbTable(conn, PDBS_TABLE);
		
		UniProtConnection uniprot = new UniProtConnection();	// online
		SiftsConnection sifts = new SiftsConnection();			// online
		
		String[] genes = conn.getStringsFromDb(String.format(QUERY_CANDIDATE_GENES,MIN_NUM_MUT));
		System.out.printf("%s\t%s\t%s\n", "Gene", "UniProt", "SIFTS");
		for(String geneName:genes) {
			Gene g = new Gene(geneName);
			String uniprotId = conn.getStringFromDb(String.format(QUERY_GENE2SP, geneName));
			g.setUniprotId(uniprotId);
			int numUni = 0;
			int numSifts = 0;
			if(g.getUniprotId() != null) {
				try {
					UniProtEntry entry = uniprot.getEntry(uniprotId);
					Collection<UniProtPdbRef> pdbRefs = UniProtConnection.getPdbRefs(entry);
					numUni = pdbRefs.size();
					if(!USE_SIFTS) {
						for(UniProtPdbRef ref:pdbRefs) {
							if(PDBS_TABLE != null) {
								ref.writeToDb(conn, PDBS_TABLE);
							} else {
								System.out.println(ref.toString());
							}
						}
					}
				} catch (NoMatchFoundException e) {
					System.err.println("Error. No UniProt entry found for " + geneName + " " + uniprotId);
					continue;
				}
				try {
					Collection<SiftsFeature> sfs = sifts.getUniprot2PdbMappings(g.getUniprotId());
					numSifts = sfs.size();
					if(USE_SIFTS) {
						for(SiftsFeature sf:sfs) {
							// construct a UniProtPdbRef object
							String pdbCode = sf.getPdbCode();
							String[] chains = new String[1];
							chains[0] = sf.getPdbChainCode();
							int beg = sf.getIntervalSet().first().beg;
							int end = sf.getIntervalSet().first().end;
							// getting experimental method and resolution
							int entryKey = conn.getIntFromDb(String.format(PDBASE_GET_KEY, pdbCode));
							String method = conn.getStringFromDb(String.format(PDBASE_GET_METHOD, entryKey));
							if(method == null) {
								System.err.println("Error: No experimental method found for " + pdbCode);
							} else
							switch(Character.toUpperCase(method.charAt(0))) {
								case 'S': method = "N"; break;
								case 'X': method = "X"; break;
								default: System.err.println("Warning: Unknown experimental method " + method + " found in Pdbase for " + pdbCode);
							}
							String resolStr = conn.getStringFromDb(String.format(PDBASE_GET_RESOL, entryKey));
							double resol = resolStr==null?-0:Double.parseDouble(resolStr);
							UniProtPdbRef ref = new UniProtPdbRef(geneName,uniprotId,pdbCode,chains,method,resol,beg,end);
							if(PDBS_TABLE != null) {
								ref.writeToDb(conn, PDBS_TABLE);
							} else {
								System.out.println(ref.toString());
							}												
						}
					}
				} catch (NoMatchFoundException e) {
					System.err.println("Warning: No SIFTS record found for " + geneName + " " + uniprotId);
				}
			}
			System.out.printf("%s\t%d\t%d\n", geneName, numUni, numSifts);
		}
	}
		
	/**
	 * Print information about each gene which helps to find the right Ensemble transcript. We are pretty sure
	 * we know the right Ensembl gene but it is not clear which is the right transcript.
	 */
	public static void findTranscript(MySQLConnection conn) {
		//String[] geneNames = conn.getStringsFromDb(QUERY_ALL_GENES);
		String[] geneNames = conn.getStringsFromDb(String.format(QUERY_CANDIDATE_GENES,7));
		
		for(String geneName:geneNames) {
			System.out.println(geneName);
			// get Uniprot ID
			Gene g = new Gene(geneName);			
			g.setUniprotId(conn.getStringFromDb(String.format(QUERY_GENE2SP, geneName)));
			g.loadCosmicSequence(ONLINE);
			
			// get ENSG
			String[] ensgs = conn.getStringsFromDb("SELECT DISTINCT ensg FROM mutanom2.sp2ens_top148 WHERE sp_id=\"" + g.getUniprotId() + "\"");
			if(ensgs != null) {
			System.out.println("  ENSGs found: " + ensgs.length);
			for(String ensg:ensgs) {
				System.out.println("  " + ensg);
				String[] ensts = conn.getStringsFromDb("SELECT DISTINCT enst FROM mutanom2.sp2ens_top148 WHERE ensg=\"" + ensg + "\"");
				if(ensts != null) {
					System.out.println("  ENSTs found: " + ensts.length);
					for(String enst:ensts) {
						System.out.println("    " + enst);
						String[] snps = conn.getStringsFromDb("SELECT DISTINCT CONCAT(mut_aa,beg_aa) FROM mutanom2.snps_ensvar59_top29 WHERE mut_aa rlike \"^[A-Z]/[A-Z]$\" AND beg_aa = end_aa AND enst=\"" + enst + "\"");
						int matches = 0;
						if(snps != null) {
						for(String snp:snps) {
							char from = snp.charAt(0);
							//char to = snp.charAt(2);
							int pos = Integer.parseInt(snp.substring(3));
							if(g.getCosmicSeq().length() >= pos && g.getCosmicSeq().charAt(pos-1)==from) matches++;
						}
						}
						System.out.println("    SNPs found: " + (snps==null?0:snps.length) + " Matching: " + matches);
						//String[] snps2 = conn.getStringsFromDb("SELECT DISTINCT CONCAT(mut_aa,beg_aa) FROM mutanom2.snps_ensvar59 WHERE mut_aa rlike \"[A-Z]/[A-Z]\" AND beg_aa = end_aa AND enst=\"" + enst + "\" AND genot_freq >= 0.01 and genot_freq < 0.5");
						//System.out.println("    SNPs found: " + (snps2==null?0:snps2.length) + " (0.01 <= f < 0.5)");			
					}
				} else {
					System.out.println("  ENSTs found: 0");
				}
			}
			} else {
				System.out.println("  ENSGs found: 0");
			}
		}	
	}
	
	/**
	 * Print statistics about the SNPs per gene.
	 */
	public static void snpStats(MySQLConnection conn) {
		String[] geneNames = conn.getStringsFromDb(String.format(QUERY_CANDIDATE_GENES,5));
		
		System.out.println("Gene\tsnpper\tmtch_co\tmtch_sp\tenst\tensg\tcommon");
		for(String geneName:geneNames) {
			System.out.print(geneName + "\t");
			// get Uniprot ID
			Gene g = new Gene(geneName);			
			g.setUniprotId(conn.getStringFromDb(String.format(QUERY_GENE2SP, geneName)));
			g.loadCosmicSequence(ONLINE);
			g.loadUniprotSequence(ONLINE);
			
			// get SNPper SNPs
			String[] snps = conn.getStringsFromDb("SELECT DISTINCT CONCAT(mut_aa,pos_aa) FROM mutanom2.snps_snpper_top29 WHERE mut_aa rlike \"^[A-Z]/[A-Z]$\" AND substr(mut_aa,1,1)!=substr(mut_aa,3,1) AND gene_name=\"" + geneName + "\"");
			System.out.print((snps==null?0:snps.length) + "\t");
			// cosmic matches
			int matches = 0;
			if(snps!=null) {
				for(String snp:snps) {
					char from = snp.charAt(0);
					//char to = snp.charAt(2);
					int pos = Integer.parseInt(snp.substring(3));
					if(g.getCosmicSeq().length() >= pos && g.getCosmicSeq().charAt(pos-1)==from) matches++;
				}
			}
			System.out.print(matches + "\t");
			// uniprot matches
			matches = 0;
			if(snps!=null) {
				for(String snp:snps) {
					char from = snp.charAt(0);
					//char to = snp.charAt(2);
					int pos = Integer.parseInt(snp.substring(3));
					if(g.getUniprotSeq().length() >= pos && g.getUniprotSeq().charAt(pos-1)==from) matches++;
				}
			}
			System.out.print(matches + "\t");
			
			//String[] snps2 = conn.getStringsFromDb("SELECT DISTINCT dbsnp_id FROM mutanom2.snps_snpper_top29 WHERE mut_aa rlike \"^[A-Z]/[A-Z]$\" AND substr(mut_aa,1,1)!=substr(mut_aa,3,1) AND gene_name=\"" + geneName + "\"");
			//System.out.print((snps2==null?0:snps2.length) + "\t");			
			
			// get Ensembl SNPs
			String[] ensgs = conn.getStringsFromDb("SELECT DISTINCT ensg FROM mutanom2.sp2ens_top148 WHERE sp_id=\"" + g.getUniprotId() + "\"");
			if(ensgs == null) {
				System.out.print("No ENSG found!");
			} else {
				String[] ensts = conn.getStringsFromDb("SELECT DISTINCT enst FROM mutanom2.sp2ens_top148 WHERE ensg=\"" + ensgs[0] + "\"");
				if(ensts == null) {
					System.out.print("No ENST found!");
				} else {
						String[] snps3 = conn.getStringsFromDb("SELECT DISTINCT snp_id FROM mutanom2.snps_ensvar59_top29 WHERE mut_aa rlike \"[A-Z]/[A-Z]\" AND beg_aa = end_aa AND enst=\"" + ensts[0] + "\"");
						System.out.print((snps3==null?0:snps3.length) + "\t");
						//String[] snps4 = conn.getStringsFromDb("SELECT DISTINCT mut_aa FROM mutanom2.snps_ensvar59 WHERE mut_aa rlike \"[A-Z]/[A-Z]\" AND beg_aa = end_aa AND enst=\"" + ensts[0] + "\" AND genot_freq >= 0.01 and genot_freq < 0.5");
						//System.out.print(snps4==null?0:snps4.length);			
						String[] snps5 = conn.getStringsFromDb("SELECT DISTINCT snp_id FROM mutanom2.snps_ensvar59_top29 WHERE mut_aa rlike \"[A-Z]/[A-Z]\" AND beg_aa = end_aa AND ensg=\"" + ensgs[0] + "\"");
						System.out.print((snps5==null?0:snps5.length) + "\t");
						// SNP ids which are in common between snpper and ensvar (by ensg)
						String[] snps6= conn.getStringsFromDb("select distinct s.dbsnp_id from mutanom2.snps_ensvar59_top29 as e inner join mutanom2.snps_snpper_top29 as s on e.snp_id=s.dbsnp_id where s.pos_aa = e.beg_aa and s.mut_aa = e.mut_aa and ensg=\""+ ensgs[0] + "\" and gene_name=\"" + geneName + "\";");
						System.out.print((snps6==null?0:snps6.length) + "\t");
						System.out.print(ensgs[0] + "\t");
						System.out.print(ensts[0] + "\t");
				}
			}
			System.out.println(".");
		}	
	}	
	
	/**
	 * Creates the gene2sp table (if not exists).
	 * @param conn
	 * @throws SQLException
	 */
	public static void createGene2SpTable(MySQLConnection conn) throws SQLException {
		conn.executeSql(CREATE_GENE2SP_TABLE);
	}
	
	/**
	 * Creates the genes table (if not exists).
	 * @param conn
	 * @throws SQLException
	 */
	public static void createGenesTable(MySQLConnection conn) throws SQLException {
		conn.executeSql(CREATE_GENES_TABLE);
	}
	
	/**
	 * Some Uniprot ids were found to be wrong or outdated. We are fixing this here (semi-) manually.
	 * TODO: Automatically query Uniprot for updated IDs.
	 * @param conn
	 * @throws SQLException
	 */
	public static void updateUniprotIds(MySQLConnection conn) throws SQLException {
		conn.executeSql("UPDATE " + GENE2SP_TABLE + " SET sp_id=\"Q5VST9\" WHERE sp_id=\"Q96AA2\"");	// outdated id
		conn.executeSql("UPDATE " + GENE2SP_TABLE + " SET sp_id=\"P08F94\" WHERE sp_id=\"Q8TCZ9\""); 	// outdated id
		conn.executeSql("UPDATE " + GENE2SP_TABLE + " SET sp_id=\"Q70Z35\" WHERE sp_id=\"Q9H805\"");	// outdated id
		conn.executeSql("UPDATE " + GENE2SP_TABLE + " SET sp_id=\"P16473\" WHERE gene_name=\"TSHR\"");	// was Q59GA2 (wrong blast hit)
	}
	
	public static void main(String[] args) throws SQLException, IOException {
		
		// TODO:
		// Usage: GeneSelection <command> <database>
		// Commands:
		// -u find uniprot ids and write to table 'gene2sp' (creates table first)
		// -g count mutations and write to table 'genes' (creates table first)
		// -s write number of SNPs per gene to screen
		// -p load PDBs for per gene to table 'pdbs'
		// -t write Ensembl transcript ids for each UniprotID to screen
		
		MySQLConnection conn = new MySQLConnection();
		
		// --- Creating database tables ---
		
		//writeGene2SpTable(conn);	// writes uniprot Ids to database table 'gene2sp' (created if not exists)
		
		// updateUniprotIds(conn);	// manually update some outdated uniprot ids with the latest ones
		
		//writeGenesTable(conn);	// reads uniprot ids from database and counts # of mutations in known regions + writes to 'genes' table
				
		//writePdbsTable(conn);		// loads all PDB cross references for each gene to database table 'pdbs'
		
		// --- Information methods which use the tables from above ---
		
		// snpStats(conn);				// prints a table with the number of total missense SNPs found for each of the top29 genes
		
		//findTranscript(conn);		// shows which ensembl ENSG and ENST records are founds for each UniprotID
		
		conn.close();
		System.out.println("done.");
	}
	
}
