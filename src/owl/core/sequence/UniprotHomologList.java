package owl.core.sequence;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;

import org.xml.sax.SAXException;

import owl.core.connections.EmblWSDBfetchConnection;
import owl.core.connections.NoMatchFoundException;
import owl.core.connections.UniProtConnection;
import owl.core.runners.blast.BlastError;
import owl.core.runners.blast.BlastHit;
import owl.core.runners.blast.BlastHitList;
import owl.core.runners.blast.BlastRunner;
import owl.core.runners.blast.BlastXMLParser;
import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseType;
import uk.ac.ebi.kraken.interfaces.uniprot.NcbiTaxonomyId;
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

	private static final boolean DEBUG = false;
	
	/*-------------------------- members --------------------------*/
	
	private Sequence seq; // the sequence to which the homologs refer
	private List<UniprotHomolog> list; // the list of homologs
	private Map<String,UniprotHomolog> lookup; // to speed up searches (uniprot ids to Homologs)
	
	
	public UniprotHomologList(String tag, String sequence) {
		seq = new Sequence(tag, sequence);
	}
	
	public UniprotHomologList(Sequence seq) {
		this.seq = seq;
	}
	
	/**
	 * Performs a blast search to populate this list of homologs.
	 * All blast output files will be removed on exit.
	 * @param blastBinDir
	 * @param blastDbDir
	 * @param blastDb
	 * @param blastNumThreads
	 * @throws IOException
	 * @throws BlastError
	 */
	public void searchWithBlast(String blastBinDir, String blastDbDir, String blastDb, int blastNumThreads) throws IOException, BlastError {
		File outBlast = File.createTempFile(BLAST_BASENAME,BLASTOUT_SUFFIX);
		File inputSeqFile = File.createTempFile(BLAST_BASENAME,FASTA_SUFFIX);
		if (!DEBUG) {
			outBlast.deleteOnExit();
			inputSeqFile.deleteOnExit();
		}
		seq.writeToFastaFile(inputSeqFile);
		BlastRunner blastRunner = new BlastRunner(blastBinDir, blastDbDir);
		blastRunner.runBlastp(inputSeqFile, blastDb, outBlast, BLAST_OUTPUT_TYPE, BLAST_NO_FILTERING, blastNumThreads);
		
		BlastHitList blastList = null;
		try {
			BlastXMLParser blastParser = new BlastXMLParser(outBlast);
			blastList = blastParser.getHits();
		} catch (SAXException e) {
			// if this happens it means that blast doesn't format correctly its XML, i.e. has a bug
			System.err.println("Unexpected error: "+e.getMessage());
			System.exit(1);
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
		this.lookup = new HashMap<String, UniprotHomolog>();
		for (UniprotHomolog hom:this) {
			lookup.put(hom.getUniId(), hom);
		}
	}
	
	/**
	 * Retrieves from UniprotKB the taxonomy and EMBL CDS ids data,
	 * by using the remote Uniprot API
	 */
	public void retrieveUniprotKBData() {
		UniProtConnection uniprotConn = new UniProtConnection();
		List<String> uniprotIds = new ArrayList<String>();
		for (UniprotHomolog hom:this) {
			uniprotIds.add(hom.getUniId());
		}
		EntryIterator<UniProtEntry> entries = uniprotConn.getMultipleEntries(uniprotIds);

		for (UniProtEntry entry:entries) {
			UniprotHomolog hom = this.getHomolog(entry.getPrimaryUniProtAccession().getValue()); 
			List<String> taxIds = new ArrayList<String>();
			for(NcbiTaxonomyId ncbiTaxId:entry.getNcbiTaxonomyIds()) {
				taxIds.add(ncbiTaxId.getValue());
			}
			hom.setTaxIds(taxIds);
			Collection<Embl> emblrefs = entry.getDatabaseCrossReferences(DatabaseType.EMBL);
			List<String> emblCdsIds = new ArrayList<String>();
			for(Embl ref:emblrefs) {
				String emblCdsIdWithVer = ref.getEmblProteinId().getValue();
				String emblCdsId = emblCdsIdWithVer.substring(0, emblCdsIdWithVer.lastIndexOf("."));
    			emblCdsIds.add(emblCdsId);
    		}
			hom.setEmblCdsIds(emblCdsIds);
		}
	}
	
	/**
	 * Retrieves from EMBL DB fetch web service the EMBL CDS sequences
	 * @throws IOException
	 */
	public void retrieveEmblCdsSeqs() throws IOException {
		List<String> allIds = new ArrayList<String>();
		
		for (UniprotHomolog hom:this) {
			allIds.addAll(hom.getEmblCdsIds());
		}
		List<Sequence> allSeqs = null;
		try {
			allSeqs = EmblWSDBfetchConnection.fetchEMBLCDS(allIds);
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
	public UniprotHomolog getHomolog(String uniprotId) {
		return this.lookup.get(uniprotId);
	}
	
	public Iterator<UniprotHomolog> iterator() {
		return this.list.iterator();
	}

 
}
