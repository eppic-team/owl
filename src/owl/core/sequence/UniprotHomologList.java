package owl.core.sequence;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.xml.sax.SAXException;

import owl.core.connections.EmblWSDBfetchConnection;
import owl.core.connections.NoMatchFoundException;
import owl.core.connections.UniProtConnection;
import owl.core.runners.TcoffeeError;
import owl.core.runners.TcoffeeRunner;
import owl.core.runners.blast.BlastError;
import owl.core.runners.blast.BlastHit;
import owl.core.runners.blast.BlastHitList;
import owl.core.runners.blast.BlastRunner;
import owl.core.runners.blast.BlastXMLParser;
import owl.core.sequence.alignment.AlignmentConstructionError;
import owl.core.sequence.alignment.MultipleSequenceAlignment;
import owl.core.util.FileFormatError;
import owl.core.util.Goodies;
import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseType;
import uk.ac.ebi.kraken.interfaces.uniprot.NcbiTaxonomyId;
import uk.ac.ebi.kraken.interfaces.uniprot.Organelle;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.dbx.embl.Embl;
import uk.ac.ebi.kraken.uuw.services.remoting.EntryIterator;

/**
 * Class to store a set of uniprot homologs of a given sequence.
 * It contains methods to blast against a uniprot database to get the homolog
 * list and to retrieve data from uniprot (taxonomy, dbrefs) and embl cds sequences.
 * 
 * @see UniprotHomolog
 * 
 * @author duarte_j
 *
 */
public class UniprotHomologList implements Iterable<UniprotHomolog>{

	/*------------------------ constants --------------------------*/
	
	private static final String BLASTOUT_SUFFIX = "blast.out.xml";
	private static final String FASTA_SUFFIX = ".fa";
	private static final String BLAST_BASENAME = "homSearch";
	
	private static final int 	BLAST_OUTPUT_TYPE = 7;  // xml output
	private static final boolean BLAST_NO_FILTERING = true;
	private static final String UNIPROT_VER_FILE = "reldate.txt";
	
	private static final String  TCOFFEE_ALN_OUTFORMAT = "fasta";
	
	private static final boolean DEBUG = false;
	
	/*-------------------------- members --------------------------*/
	
	private Sequence seq; 							 // the sequence to which the homologs refer
	private List<UniprotHomolog> list; 				 // the list of homologs
	private Map<String,List<UniprotHomolog>> lookup; // to speed up searches (uniprot ids to Homologs lists) 
													 // it's a list because blast can hit a single uniprot in multiple regions (for us
													 // that's multiple BlastHits)
	private double idCutoff; 						 // the identity cutoff (see restrictToMinId() )
	private String uniprotVer;						 // the version of uniprot used in blasting, read from the reldate.txt uniprot file
	
	
	
	public UniprotHomologList(String tag, String sequence) {
		this.seq = new Sequence(tag, sequence);
		this.idCutoff = 0.0; // i.e. no filter
	}
	
	public UniprotHomologList(Sequence seq) {
		this.seq = seq;
		this.idCutoff = 0.0; // i.e. no filter
	}
	
	/**
	 * Performs a blast search to populate this list of homologs.
	 * All blast output files will be removed on exit.
	 * @param blastBinDir
	 * @param blastDbDir
	 * @param blastDb
	 * @param blastNumThreads
	 * @param cacheFile a file with the cached xml blast output file, if null blast will be always run
	 * @throws IOException
	 * @throws BlastError
	 */
	public void searchWithBlast(String blastBinDir, String blastDbDir, String blastDb, int blastNumThreads, File cacheFile) throws IOException, BlastError {
		File outBlast = null;
		boolean fromCache = false;
		if (cacheFile!=null && cacheFile.exists()) {
			outBlast = cacheFile;
			fromCache = true;
			System.out.println("Reading blast results from cache file "+cacheFile);
		} else {
			outBlast = File.createTempFile(BLAST_BASENAME,BLASTOUT_SUFFIX);
			File inputSeqFile = File.createTempFile(BLAST_BASENAME,FASTA_SUFFIX);
			if (!DEBUG) {
				outBlast.deleteOnExit();
				inputSeqFile.deleteOnExit();
			}
			seq.writeToFastaFile(inputSeqFile);
			BlastRunner blastRunner = new BlastRunner(blastBinDir, blastDbDir);
			blastRunner.runBlastp(inputSeqFile, blastDb, outBlast, BLAST_OUTPUT_TYPE, BLAST_NO_FILTERING, blastNumThreads);
			this.uniprotVer = readUniprotVer(blastDbDir);
			try {
				Goodies.copyFile(outBlast, cacheFile);
			} catch (IOException e) {
				System.err.println("Couldn't write the blast cache file "+cacheFile);
				System.err.println(e.getMessage());
			}
		}
		
		BlastHitList blastList = null;
		try {
			BlastXMLParser blastParser = new BlastXMLParser(outBlast);
			blastList = blastParser.getHits();
		} catch (SAXException e) {
			if (fromCache) {
				throw new IOException("Cache file "+cacheFile+" does not comply with blast XML format.");
			} else {
				// if this happens it means that blast doesn't format correctly its XML, i.e. has a bug
				System.err.println("Unexpected error: "+e.getMessage());
				System.exit(1);
			}
		}

		if (fromCache) {
			if (!blastList.getQueryId().equals(seq.getName())) {
				throw new IOException("Query id from cache file "+cacheFile+" does not match the id from the sequence: "+seq.getName());
			}
			this.uniprotVer = readUniprotVer(cacheFile.getParent());
			String uniprotVerFromBlastDbDir = readUniprotVer(blastDbDir);
			if (!uniprotVerFromBlastDbDir.equals(uniprotVer)) {
				throw new IOException("Uniprot version from blast db dir "+blastDbDir+" does not match version in cache dir "+cacheFile.getParent());
			}
		}

		
		this.list = new ArrayList<UniprotHomolog>();
		for (BlastHit hit:blastList) {
			String sid = hit.getSubjectId();
			Matcher m = Sequence.DEFLINE_PRIM_ACCESSION_REGEX.matcher(sid);
			if (m.matches()) {
				String uniId = m.group(1);
				list.add(new UniprotHomolog(hit,uniId));
			} else {
				System.err.println("Could not find uniprot id in subject id "+sid);
			}
		}
		initialiseMap();
	}
	
	//TODO write a searchithPSIBlast method

	private void initialiseMap() {
		this.lookup = new HashMap<String, List<UniprotHomolog>>();
		for (UniprotHomolog hom:this) {
			if (lookup.containsKey(hom.getUniId())) {
				lookup.get(hom.getUniId()).add(hom);	
			} else {
				List<UniprotHomolog> list = new ArrayList<UniprotHomolog>();
				list.add(hom);
				lookup.put(hom.getUniId(), list);
			}
			
		}
	}
	
	private String readUniprotVer(String blastDbDir) {
		String ver = "unknown";
		File uniprotVerFile = new File(blastDbDir,UNIPROT_VER_FILE);
		try {
			
			BufferedReader br = new BufferedReader(new FileReader(uniprotVerFile));
			String line;
			Pattern p = Pattern.compile("^UniProt\\sKnowledgebase\\sRelease\\s([\\d._]+)\\s.*");
			while ((line=br.readLine())!=null){
				Matcher m = p.matcher(line);
				if (m.matches()) {
					ver = m.group(1);
					break;
				}
			}
			br.close();
		} catch(IOException e) {
			System.err.println("Warning: couldn't read uniprot version from file "+uniprotVerFile);
		}
		return ver;
	}
	
	/**
	 * Retrieves from UniprotKB the sequence, taxonomy and EMBL CDS ids data,
	 * by using the remote Uniprot API
	 */
	public void retrieveUniprotKBData() {
		UniProtConnection uniprotConn = new UniProtConnection();
		if (!uniprotConn.getVersion().equals(this.uniprotVer)){
			System.err.println("Warning! Uniprot version used for blast ("+uniprotVer+") and uniprot version being queried with api ("+uniprotConn.getVersion()+")don't match!");
		}
		List<String> uniprotIds = new ArrayList<String>();
		for (UniprotHomolog hom:this) {
			uniprotIds.add(hom.getUniId());
		}
		EntryIterator<UniProtEntry> entries = uniprotConn.getMultipleEntries(uniprotIds);

		for (UniProtEntry entry:entries) {
			List<UniprotHomolog> homs = this.getHomolog(entry.getPrimaryUniProtAccession().getValue());
			for (UniprotHomolog hom:homs) {
				hom.setUniprotSeq(new Sequence(hom.getUniId(),entry.getSequence().getValue()));
				List<String> taxIds = new ArrayList<String>();
				for(NcbiTaxonomyId ncbiTaxId:entry.getNcbiTaxonomyIds()) {
					taxIds.add(ncbiTaxId.getValue());
				}
				hom.setTaxIds(taxIds);
				Collection<Embl> emblrefs = entry.getDatabaseCrossReferences(DatabaseType.EMBL);
				List<String> emblCdsIds = new ArrayList<String>();
				Set<String> tmpEmblCdsIdsSet = new TreeSet<String>();
				for(Embl ref:emblrefs) {
					String emblCdsIdWithVer = ref.getEmblProteinId().getValue();
					if (!emblCdsIdWithVer.equals("-")) { // for non annotated genomic dna cds sequences the identifier is '-', we ignore them
						String emblCdsId = emblCdsIdWithVer.substring(0, emblCdsIdWithVer.lastIndexOf("."));
						//emblCdsIds.add(emblCdsId);
						tmpEmblCdsIdsSet.add(emblCdsId);
					}
				}
				emblCdsIds.addAll(tmpEmblCdsIdsSet); // we use the set to be sure there are no duplicates (it does happen sometimes)
				hom.setEmblCdsIds(emblCdsIds);
				
				List<Organelle> orglls = entry.getOrganelles();
				if (orglls.size()>0) {
					hom.setGeneEncodingOrganelle(orglls.get(0).getType());
					if (orglls.size()>1) {
						for (Organelle orgll:orglls){ 
							if (!orgll.getType().equals(hom.getGeneEncodingOrganelle())) {
								System.err.println("Warning! Different gene encoding organelles for Uniprot "+hom.getUniId());
							}
						}
					}
				}
			}
		}
	}
	
	/**
	 * Retrieves from EMBL DB fetch web service the EMBL CDS sequences
	 * @param cacheFile a FASTA file containing the sequences to retrieve. If present and if
	 * it contains ALL required sequences then they are read from cacheFile. If null or file
	 * does not exist or file older than {@value #EmblWSDBfetchConnect.MAX_CACHE_AGE} then the sequences are 
	 * retrieved from EMBL DB fetch
	 * @throws IOException
	 */
	public void retrieveEmblCdsSeqs(File cacheFile) throws IOException {
		List<String> allIds = new ArrayList<String>();
		
		for (UniprotHomolog hom:this) {
			allIds.addAll(hom.getEmblCdsIds());
		}
		List<Sequence> allSeqs = null;
		try {
			allSeqs = EmblWSDBfetchConnection.fetchEMBLCDS(allIds, cacheFile);
		} catch (NoMatchFoundException e) {
			// this is unlikely to happen here, that's why we don't write a better error message
			System.err.println("Couldn't retrieve EMBL CDS sequences for some EMBL cds ids"); 
		}
		
		// we put the list (containing all the sequences from all the homologs) in a lookup table
		// so that we can then retrieve the ones corresponding to each homolog below
		Map<String,Sequence> lookup = new HashMap<String,Sequence>();
		for (Sequence seq:allSeqs) {
			lookup.put(seq.getSecondaryAccession(), seq);
		}
		 
		for (UniprotHomolog hom:this) {
			List<Sequence> seqs = new ArrayList<Sequence>();
			for (String emblCdsId:hom.getEmblCdsIds()) {
				seqs.add(lookup.get(emblCdsId));
			}
			hom.setEmblCdsSeqs(seqs);
		}		
	}
	
	/**
	 * Gets the Homolog given a uniprot ID
	 * @param uniprotId
	 * @return
	 */
	public List<UniprotHomolog> getHomolog(String uniprotId) {
		return this.lookup.get(uniprotId);
	}
	
	public Iterator<UniprotHomolog> iterator() {
		return this.list.iterator();
	}

	/**
	 * Write to the given file the query sequence and all the homolog sequences (full 
	 * uniprot sequences) in fasta format.
	 * @param outFile
	 * @throws FileNotFoundException
	 */
	public void writeToFasta(File outFile) throws FileNotFoundException {
		PrintWriter pw = new PrintWriter(outFile);
		
		int len = 80;

		pw.println(MultipleSequenceAlignment.FASTAHEADER_CHAR + seq.getName());
		for(int i=0; i<seq.getSeq().length(); i+=len) {
			pw.println(seq.getSeq().substring(i, Math.min(i+len,seq.getSeq().length())));
		}
		
		for(UniprotHomolog hom:this) {
			
			String sequence = hom.getUniprotSeq().getSeq();
			pw.println(MultipleSequenceAlignment.FASTAHEADER_CHAR + hom.getBlastHit().getSubjectId());
			for(int i=0; i<sequence.length(); i+=len) {
				pw.println(sequence.substring(i, Math.min(i+len,sequence.length())));
			}
		}
		pw.println();
		pw.close();
	}
	
	/**
	 * Runs t_coffee to align all sequences of homologs and the query sequence
	 * returning a MultipleSequenceAlignment object
	 * @param tcoffeeBin
	 * @param veryFast whether to use t_coffee's very fast alignment (and less accurate) mode
	 * @return
	 * @throws IOException
	 * @throws TcoffeeError 
	 */
	public MultipleSequenceAlignment getTcoffeeAlignment(File tcoffeeBin, boolean veryFast) throws IOException, TcoffeeError {
		File homologSeqsFile = File.createTempFile("homologs.", ".fa");
		File outTreeFile = File.createTempFile("homologs.", ".dnd");
		File alnFile = File.createTempFile("homologs.",".aln");
		File tcoffeeLogFile = File.createTempFile("homologs.",".tcoffee.log");
		if (!DEBUG) {
			homologSeqsFile.deleteOnExit();
			alnFile.deleteOnExit();
			tcoffeeLogFile.deleteOnExit();
			outTreeFile.deleteOnExit(); 
		}
		this.writeToFasta(homologSeqsFile);
		TcoffeeRunner tcr = new TcoffeeRunner(tcoffeeBin);
		tcr.runTcoffee(homologSeqsFile, alnFile, TCOFFEE_ALN_OUTFORMAT, outTreeFile, null, tcoffeeLogFile, veryFast);

			MultipleSequenceAlignment aln = null;
		
		try {
			aln = new MultipleSequenceAlignment(alnFile.getAbsolutePath(), MultipleSequenceAlignment.FASTAFORMAT);
		} catch (FileFormatError e) {
			System.err.println("Unexpected error, output file of tcoffee "+alnFile+" does not seem to be in the right format.");
			System.err.println("Error: "+e.getMessage());
			System.exit(1);
		} catch (AlignmentConstructionError e) {
			System.err.println("Unexpected error, output file of tcoffee "+alnFile+" seems to contain .");
			System.err.println("Error: "+e.getMessage());
			System.exit(1);
		}
		
		return aln;
	}
 
	/**
	 * Removes homologs below the given sequence identity value.
	 * @param idCutoff
	 */
	public void restrictToMinId(double idCutoff) {
		this.idCutoff = idCutoff;
		List<String> toRemove = new ArrayList<String>();
		Iterator<UniprotHomolog> it = list.iterator();
		while (it.hasNext()) {
			UniprotHomolog hom = it.next();
			if ((hom.getBlastHit().getPercentIdentity()/100.0)<=idCutoff) {
				it.remove();
				toRemove.add(hom.getUniId());
			}
		}
		// updating also lookup table
		// 1st we go through toRemove uniIds and search each list for ids below cutoff
		// thus it can happen that a uniId deletes all members of a list in its first appearance
		for (String uniId:toRemove) {
			List<UniprotHomolog> homs = lookup.get(uniId);
			Iterator<UniprotHomolog> it2 = homs.iterator();
			while (it2.hasNext()) {
				UniprotHomolog hom = it2.next();
				if ((hom.getBlastHit().getPercentIdentity()/100.0)<=idCutoff) {
					it2.remove();
				}
			}
		}
		// after the purge, now we get rid of uniIds mapping to empty list
		Iterator<List<UniprotHomolog>> it3 = lookup.values().iterator(); 
		while (it3.hasNext()) {
			List<UniprotHomolog> homs = it3.next();
			if (homs.isEmpty()) {
				it3.remove();
			}
		}

	}
	
	/**
	 * Returns the number of homologs in this list
	 * @return
	 */
	public int size() {
		return list.size();
	}
	
	/**
	 * Returns the sequence identity cutoff, see {@link #restrictToMinId(double)}
	 * @return
	 */
	public double getIdCutoff() {
		return idCutoff;
	}
	
	/**
	 * Gets the uniprot version used for blasting
	 * @return
	 */
	public String getUniprotVer() {
		return uniprotVer;
	}
	
	public void checkEmblCDSMatching() {
		for (UniprotHomolog hom:this) {
			hom.checkEmblCDSMatching();
		}		
	}
	
	public int getNumHomologsWithCDS() {
		int count = 0;
		for (UniprotHomolog hom:this) {
			if (hom.hasCDS()) {
				count++;
			}
		}				
		return count;
	}
	
	public int getNumHomologsWithValidCDS() {
		int count = 0;
		for (UniprotHomolog hom:this) {
			if (hom.getRepresentativeCDS()!=null) {
				count++;
			}
		}
		return count;
	}

}
