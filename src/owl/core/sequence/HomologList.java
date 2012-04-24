package owl.core.sequence;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.Serializable;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.xml.sax.SAXException;

import owl.core.connections.EmblWSDBfetchConnection;
import owl.core.connections.NoMatchFoundException;
import owl.core.connections.UniProtConnection;
import owl.core.connections.UniprotLocalConnection;
import owl.core.runners.TcoffeeException;
import owl.core.runners.TcoffeeRunner;
import owl.core.runners.blast.BlastException;
import owl.core.runners.blast.BlastHit;
import owl.core.runners.blast.BlastHitList;
import owl.core.runners.blast.BlastRunner;
import owl.core.runners.blast.BlastXMLParser;
import owl.core.sequence.alignment.AlignmentConstructionException;
import owl.core.sequence.alignment.MultipleSequenceAlignment;
import owl.core.sequence.alignment.PairwiseSequenceAlignment;
import owl.core.sequence.alignment.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;
import owl.core.util.FileFormatException;
import owl.core.util.Goodies;
import owl.core.util.Interval;

/**
 * Class to store a set of homologs of a given sequence.
 * It contains methods to blast against a uniprot/uniref database to get the homolog
 * list and to retrieve data from Uniprot (taxonomy, dbrefs), embl cds sequences, 
 * Uniparc sequences...
 * 
 * @see Homolog
 * 
 * @author duarte_j
 *
 */
public class HomologList implements  Serializable {//Iterable<UniprotHomolog>,

	private static final long serialVersionUID = 1L;

	/*------------------------ constants --------------------------*/
	
	private static final String BLASTOUT_SUFFIX = "blast.out.xml";
	private static final String FASTA_SUFFIX = ".fa";
	private static final String BLAST_BASENAME = "homSearch";
	
	private static final int 	BLAST_OUTPUT_TYPE = 7;  // xml output
	private static final boolean BLAST_NO_FILTERING = true;
	private static final String UNIPROT_VER_FILE = "reldate.txt";
	
	private static final String  TCOFFEE_ALN_OUTFORMAT = "fasta";
	
	private static final boolean DEBUG = false;
	
	private static final Log LOGGER = LogFactory.getLog(HomologList.class);

	
	/*-------------------------- members --------------------------*/
	
	private UnirefEntry ref;						// the uniref entry (uniprot/uniparc) to which the homologs refer
	private Interval refInterval;
	private boolean isSubInterval;
	private List<Homolog> list; 					// the list of homologs
	private List<Homolog> subList;			 		// the filtered list of homologs after calling filterToMinIdAndCoverage
	private Map<String,Homolog> lookup;		 		// to speed up searches (uniprot/uniparc ids to Homologs)
													// (used to be lists of homologs as we considered multi-matches of the 
													// same uniprot as different BlastHits, but not anymore since we introduced 
													// BlastHsps)
	private double idCutoff; 						// the identity cutoff (see filterToMinIdAndCoverage() )
	private double qCoverageCutoff;					// the query coverage cutoff (see filterToMinIdAndCoverage() )
	private String uniprotVer;						// the version of uniprot used in blasting, read from the reldate.txt uniprot file
	
	private MultipleSequenceAlignment aln;	  		// the protein sequences alignment

	private int reducedAlphabet;					// the reduced alphabet used to calculate entropies
	private List<Double> entropies;					// entropies for each uniprot reference sequence position
	
	private boolean useUniparc;
	
	
	public HomologList(UnirefEntry ref) {
		this(ref,null);
	}
	
	/**
	 * Create a new UniprotHomologList
	 * @param ref the uniprot entry whose sequence is the reference for this homolog list
	 * @param interv the interval in the uniprot sequence that we actually use (with 
	 * numbering of uniprot seq from 1 to length-1), if null the whole sequence is use 
	 */
	public HomologList(UnirefEntry ref, Interval interv) {
		this.ref = ref;
		if (interv!=null) {
			this.refInterval = interv;
		} else {
			this.refInterval = new Interval(1,ref.getLength());
		}
		if (refInterval.beg==1 && refInterval.end==ref.getLength()) {
			isSubInterval = false;
		} else {
			isSubInterval = true;
		}
		
		this.idCutoff = 0.0; // i.e. no filter

	}
	
	public List<Homolog> getFilteredSubset() {
		return subList;
	}
	
	/**
	 * Performs a blast search based on the reference UniprotEntry to populate this list of homologs.
	 * All blast output files will be removed on exit.
	 * @param blastBinDir
	 * @param blastDbDir
	 * @param blastDb
	 * @param blastNumThreads
	 * @param cacheFile a file with the cached xml blast output file, if null blast will be always run
	 * @throws IOException
	 * @throws BlastException
	 * @throws UniprotVerMisMatchException if uniprot versions of cacheFile given and blastDbDir do not coincide
	 * @throws InterruptedException
	 */
	public void searchWithBlast(String blastBinDir, String blastDbDir, String blastDb, int blastNumThreads, int maxNumSeqs, File cacheFile) throws IOException, BlastException, UniprotVerMisMatchException, InterruptedException {
		File outBlast = null;
		boolean fromCache = false;
		BlastHitList blastList = null;
		
		// if a sub-interval is used, we need to alter the cache file name to contain the subinterval string
		if (cacheFile!=null && isSubInterval) {
			String prefix = cacheFile.getName();
			prefix = prefix.substring(0, prefix.lastIndexOf(".blast.xml"));
			cacheFile = new File(cacheFile.getParent(),prefix+"."+refInterval.beg+"-"+refInterval.end+".blast.xml");
		}
		this.uniprotVer = readUniprotVer(blastDbDir);
		
		if (cacheFile!=null && cacheFile.exists()) {

			outBlast = cacheFile;
			fromCache = true;
			LOGGER.warn("Reading blast results from cache file "+cacheFile);

			try {
				BlastXMLParser blastParser = new BlastXMLParser(outBlast);
				blastList = blastParser.getHits();
				
				// 500 is blast's default, we don't want to check this if we are under default
				if (maxNumSeqs>BlastRunner.BLAST_DEFAULT_MAX_HITS && blastList.size()<maxNumSeqs) { 
					// we are asking for more max hits than present in the file, we have to blast again
					LOGGER.info("Blast cache file exits ("+cacheFile+") but it contains only "+blastList.size()+" hits. Need to re-blast as a max of "+maxNumSeqs+" hits have been requested");
					fromCache = false;
					blastList = null;
				} else {
					// if we do take the cache file we have to do some sanity checks
					if (!blastList.getQueryId().equals(this.ref.getUniId())) {
						throw new IOException("Query id "+blastList.getQueryId()+" from cache file "+cacheFile+
								" does not match the id from the sequence: "+this.ref.getUniId());
					}
					String uniprotVerFromCache = readUniprotVer(cacheFile.getParent());
					if (!uniprotVerFromCache.equals(uniprotVer)) {
						throw new UniprotVerMisMatchException("Uniprot version from blast db dir "+blastDbDir+
								" ("+uniprotVer+") does not match version in cache dir "+cacheFile.getParent()+" ("+uniprotVerFromCache+")");
					}
					if (!blastList.getDb().substring(blastList.getDb().lastIndexOf("/")+1).equals(blastDb)) {
						LOGGER.error("Blast db used in cache file ("+cacheFile+") different from one requested "+blastDb);
						LOGGER.error("Please check the blast cache directory.");
						System.exit(1);
					}
				}
			} catch (SAXException e) {
				throw new IOException("Cache file "+cacheFile+" does not comply with blast XML format. "+e.getMessage());
			}
		} 
		
		if (!fromCache) {
			outBlast = File.createTempFile(BLAST_BASENAME,BLASTOUT_SUFFIX);
			File inputSeqFile = File.createTempFile(BLAST_BASENAME,FASTA_SUFFIX);
			if (!DEBUG) {
				outBlast.deleteOnExit();
				inputSeqFile.deleteOnExit();
			}
			// NOTE: we blast the reference uniprot sequence using only the interval specified
			this.ref.getSeq().getInterval(this.refInterval).writeToFastaFile(inputSeqFile);
			
			BlastRunner blastRunner = new BlastRunner(blastBinDir, blastDbDir);
			blastRunner.runBlastp(inputSeqFile, blastDb, outBlast, BLAST_OUTPUT_TYPE, BLAST_NO_FILTERING, blastNumThreads, maxNumSeqs);

			LOGGER.info("Blasted against "+blastDbDir+"/"+blastDb);
			if (cacheFile!=null) {
				try {
					LOGGER.info("Writing blast cache file "+cacheFile);
					Goodies.copyFile(outBlast, cacheFile);
					cacheFile.setWritable(true, false);
				} catch (IOException e) {
					LOGGER.error("Couldn't write the blast cache file "+cacheFile);
					LOGGER.error(e.getMessage());
				}
			} 
			try {
				BlastXMLParser blastParser = new BlastXMLParser(outBlast);
				blastList = blastParser.getHits();
			} catch (SAXException e) {
				// if this happens it means that blast doesn't format correctly its XML, i.e. has a bug
				LOGGER.fatal("Unexpected error: "+e.getMessage());
				System.exit(1);
			}
		}
		
		this.list = new ArrayList<Homolog>();
		for (BlastHit hit:blastList) {
			String sid = hit.getSubjectId();
			Matcher m = Sequence.DEFLINE_PRIM_ACCESSION_REGEX.matcher(sid);
			if (m.matches()) {
				String uniId = m.group(1);
				UnirefEntry uniref = new UnirefEntry();
				uniref.setUniprotId(uniId);
				list.add(new Homolog(hit,uniref));
			} else {
				Matcher m2 = Sequence.DEFLINE_PRIM_ACCESSION_UNIREF_REGEX.matcher(sid);
				if (m2.matches()) {					
					String uniId = m2.group(1);
					if (uniId.startsWith("UPI")){
						if (useUniparc) {
							UnirefEntry uniref = new UnirefEntry();
							uniref.setUniparcId(uniId);
							list.add(new Homolog(hit,uniref));
						} else {
							LOGGER.warn("Ignoring blast hit "+uniId+" because it is a Uniparc id.");
						}
					}
					else if (uniId.contains("-")) {
						LOGGER.warn("Ignoring blast hit "+uniId+" because it is a Uniprot isoform id.");
					}
					else {	
						UnirefEntry uniref = new UnirefEntry();
						uniref.setUniprotId(uniId);
						list.add(new Homolog(hit,uniref));
					}
				} else {
					LOGGER.error("Could not find uniprot id in subject id "+sid);
				}
			}
		}
		this.subList = list; // initially the subList is the same as the list until filterToMinIdAndCoverage is called
		initialiseMap();
	}
	
	//TODO write a searchithPSIBlast method

	/**
	 * Initialises the lookup map (for speeding up lookups of homologs by homolog identifiers)
	 * Applies only to filtered subset of homologs
	 */
	private void initialiseMap() {
		this.lookup = new HashMap<String, Homolog>();
		for (Homolog hom:subList) {
			lookup.put(hom.getIdentifier(), hom);
		}
	}
	
	public static String readUniprotVer(String blastDbDir) {
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
			LOGGER.warn("Couldn't read uniprot version from file "+uniprotVerFile);
		}
		return ver;
	}
	
	/**
	 * Retrieves from UniprotKB the sequence, taxonomy and EMBL CDS ids data,
	 * by using the remote Uniprot API
	 * @param uniprotConn
	 * @throws UniprotVerMisMatchException 
	 * @throws IOException
	 */
	public void retrieveUniprotKBData(UniProtConnection uniprotConn) throws UniprotVerMisMatchException, IOException {
		
		if (!uniprotConn.getVersion().equals(this.uniprotVer)){
			throw new UniprotVerMisMatchException("Uniprot version used for blast ("+uniprotVer+") and uniprot version being queried with api ("+uniprotConn.getVersion()+") don't match!");
		}
		List<String> uniprotIds = new ArrayList<String>();
		for (Homolog hom:subList) {
			if (hom.isUniprot()) uniprotIds.add(hom.getIdentifier());
		}
		
		List<UnirefEntry> unirefs = uniprotConn.getMultipleUnirefEntries(uniprotIds);

		for (UnirefEntry uniref:unirefs) {
			Homolog hom = this.getHomolog(uniref.getUniprotId());
			hom.getUnirefEntry().setUniprotId(uniref.getUniprotId());
			hom.getUnirefEntry().setNcbiTaxId(uniref.getNcbiTaxId());
			hom.getUnirefEntry().setSequence(uniref.getSequence());
			hom.getUnirefEntry().setTaxons(uniref.getTaxons());			
		}
		
		// now we check if the query did really return all requested uniprot ids 
		HashSet<String> nonreturned = uniprotConn.getNonReturnedIdsLastMultipleRequest();
		
		if (!nonreturned.isEmpty()) {
			Iterator<Homolog> it = subList.iterator(); 
			while (it.hasNext()) {
				Homolog hom = it.next();
				if (!hom.isUniprot()) continue;
				if (nonreturned.contains(hom.getIdentifier())) {
					LOGGER.info("Removing uniprot id "+hom.getIdentifier()+" from homologs because it wasn't returned by the Uniprot connection.");
					it.remove();
				}
			}
			// and update the lookup table 
			initialiseMap();
		}
	}
	
	public void retrieveUniparcData(File cacheFile) throws IOException {
		
		List<String> allIds = new ArrayList<String>();
		for (Homolog hom:subList) {
			if (!hom.isUniprot()) {
				allIds.add(hom.getIdentifier());
			}
		}
		try {
			List<Sequence> allSeqs = EmblWSDBfetchConnection.fetchUniparc(allIds, cacheFile);
			// we put the list (containing all the sequences from all the homologs) in a lookup table
			// so that we can then retrieve the ones corresponding to each homolog below
			Map<String,Sequence> lookup = new HashMap<String,Sequence>();
			for (Sequence seq:allSeqs) {
				lookup.put(seq.getName().substring(0, seq.getName().lastIndexOf(" ")), seq);
			}
			for (Homolog hom:subList) {
				if (hom.isUniprot()) continue;
				if (lookup.containsKey(hom.getIdentifier())) {
					hom.getUnirefEntry().setSequence(lookup.get(hom.getIdentifier()).getSeq());
				}
			}		
		} catch (NoMatchFoundException e) {
			LOGGER.warn("Couldn't retrieve Uniparc sequences");
		}

	}
	
	/**
	 * Retrieves both uniprot and uniparc data from local db
	 * @param uniprotConn
	 * @throws UniprotVerMisMatchException
	 * @throws SQLException
	 */
	public void retrieveUniprotKBData(UniprotLocalConnection uniprotConn) throws UniprotVerMisMatchException, SQLException {
		if (!uniprotConn.getVersion().equals(this.uniprotVer)){
			throw new UniprotVerMisMatchException("Uniprot version used for blast ("+uniprotVer+") and uniprot version being queried from local database ("+uniprotConn.getVersion()+") don't match!");
		}
		List<String> uniIds = new ArrayList<String>();
		for (Homolog hom:subList) {
			uniIds.add(hom.getIdentifier());
		}
		
		List<UnirefEntry> unirefs = uniprotConn.getMultipleUnirefEntries(uniIds);

		for (UnirefEntry uniref:unirefs) {
			Homolog hom = this.getHomolog(uniref.getUniId());
			hom.getUnirefEntry().setId(uniref.getId());
			hom.getUnirefEntry().setUniprotId(uniref.getUniprotId());
			hom.getUnirefEntry().setUniparcId(uniref.getUniparcId());
			hom.getUnirefEntry().setNcbiTaxId(uniref.getNcbiTaxId());
			hom.getUnirefEntry().setSequence(uniref.getSequence());
			hom.getUnirefEntry().setTaxons(uniref.getTaxons());
		}
		
		// now we check if the query did really return all requested uniprot ids 
		HashSet<String> nonreturned = uniprotConn.getNonReturnedIdsLastMultipleRequest();
		
		if (!nonreturned.isEmpty()) {
			Iterator<Homolog> it = subList.iterator(); 
			while (it.hasNext()) {
				Homolog hom = it.next();
				if (nonreturned.contains(hom.getIdentifier())) {
					LOGGER.info("Removing uniprot/uniparc id "+hom.getIdentifier()+" from homologs because it wasn't returned by the Uniprot connection.");
					it.remove();
				}
			}
			// and update the lookup table 
			initialiseMap();
		}
	}
	
	/**
	 * Gets the Homolog given a uniprot/uniparc ID
	 * @param uniprotId
	 * @return
	 */
	public Homolog getHomolog(String identifier) {
		return this.lookup.get(identifier);
	}
	
	/**
	 * Write to the given file the query protein sequence and all the homolog protein 
	 * sequences in FASTA format. Only the subset of entries result of filtering with {@link #filterToMinIdAndCoverage(double, double)}
	 * will be used. If no filtered applied yet then all entries are used.
	 * @param outFile
	 * @param writeQuery if true the query sequence is written as well, if false only homologs 
	 * @throws FileNotFoundException
	 */
	public void writeToFasta(File outFile, boolean writeQuery) throws FileNotFoundException {
		PrintWriter pw = new PrintWriter(outFile);
		
		int len = 80;

		if (writeQuery) {
			pw.println(MultipleSequenceAlignment.FASTAHEADER_CHAR + this.ref.getUniId());
			Sequence refSequence = ref.getSeq().getInterval(refInterval);
			for(int i=0; i<refSequence.getLength(); i+=len) {
				pw.println(refSequence.getSeq().substring(i, Math.min(i+len,refSequence.getLength())));
			}
		}
		
		for(Homolog hom:subList) {
			
			String sequence = hom.getSequence();
			pw.println(MultipleSequenceAlignment.FASTAHEADER_CHAR + hom.getLongSequenceTag());
			for(int i=0; i<sequence.length(); i+=len) {
				pw.println(sequence.substring(i, Math.min(i+len,sequence.length())));
			}
		}
		pw.println();
		pw.close();
	}
	
	/**
	 * Runs t_coffee to align all protein sequences of homologs and the query sequence
	 * returning a MultipleSequenceAlignment object
	 * @param tcoffeeBin
	 * @param veryFast whether to use t_coffee's very fast alignment (and less accurate) mode
	 * @params nThreads number of CPU cores t_coffee should use
	 * @params alnCacheFile
	 * @throws IOException
	 * @throws TcoffeeException 
	 * @throws UniprotVerMisMatchException 
	 */
	public void computeTcoffeeAlignment(File tcoffeeBin, boolean veryFast, int nThreads, File alnCacheFile) throws IOException, TcoffeeException, InterruptedException, UniprotVerMisMatchException {
		File alnFile = null;
		
		if (alnCacheFile!=null && alnCacheFile.exists()) {
			String uniprotVerFromCacheDir = readUniprotVer(alnCacheFile.getParent());
			String uniprotVerFromBlast = this.uniprotVer; // this can be either actually from blast db dir (if blast was run) or read from blast cache dir
			if (!uniprotVerFromBlast.equals(uniprotVerFromCacheDir)) {
				throw new UniprotVerMisMatchException("Uniprot version used for blast "+
						" ("+uniprotVerFromBlast+") does not match version in alignment cache dir "+
						alnCacheFile.getParent()+" ("+uniprotVerFromCacheDir+")");
			}
			
			alnFile = alnCacheFile;
			LOGGER.warn("Reading alignment from cache file " + alnCacheFile+". Won't recompute with t_coffee");
			
		} else {
			// no cache: we compute with t-coffee

			alnFile = File.createTempFile("homologs.",".aln");
			File homologSeqsFile = File.createTempFile("homologs.", ".fa");
			File outTreeFile = File.createTempFile("homologs.", ".dnd");
			File tcoffeeLogFile = File.createTempFile("homologs.",".tcoffee.log");

			this.writeToFasta(homologSeqsFile, true);
			TcoffeeRunner tcr = new TcoffeeRunner(tcoffeeBin);
			tcr.buildCmdLine(homologSeqsFile, alnFile, TCOFFEE_ALN_OUTFORMAT, outTreeFile, null, tcoffeeLogFile, veryFast, nThreads);
			LOGGER.info("Running t_coffee command: " + tcr.getCmdLine());
			tcr.runTcoffee();
			if (!DEBUG) { 
				// note that if the run of tcoffee throws an exception, files are not marked for deletion
				homologSeqsFile.deleteOnExit();
				alnFile.deleteOnExit();
				tcoffeeLogFile.deleteOnExit();
				outTreeFile.deleteOnExit(); 
			}

		}
		
		try {
			aln = new MultipleSequenceAlignment(alnFile.getAbsolutePath(), MultipleSequenceAlignment.FASTAFORMAT);
		} catch (FileFormatException e) {
			throw new IOException(e);
		} catch (AlignmentConstructionException e) {
			throw new IOException(e);
		}

		// if we did pass a cache file but it doesn't exist yet, we have computed the alignment. Now we need to write it to the given cache file
		if (alnCacheFile!=null && !alnCacheFile.exists()) { 
			try {
				writeAlignmentToFile(alnCacheFile);
				LOGGER.info("Writing alignment cache file "+alnCacheFile);
			} catch(FileNotFoundException e) {
				LOGGER.error("Couldn't write alignment cache file "+alnCacheFile);
			}
		}

	}
	
	/**
	 * Returns the protein sequence alignment of query sequence and all homologs
	 * @return
	 */
	public MultipleSequenceAlignment getAlignment() {
		return aln;
	}
	
	public void writeAlignmentToFile(File alnFile) throws FileNotFoundException {
		aln.writeFasta(new PrintStream(alnFile), 80, true);
	}
 
	/**
	 * Creates a subset list of Homologs that have at least the required identities
	 * and query coverage.
	 * Subsequently all methods in this class will refer to the sublist rather than to original unfiltered list.
	 * Upon a subsequent call to this method the a different id/coverage can be chosen and other methods used with those.
	 * @param idCutoff
	 * @param queryCovCutoff
	 */
	public void filterToMinIdAndCoverage(double idCutoff, double queryCovCutoff) {
		this.idCutoff = idCutoff;
		this.qCoverageCutoff = queryCovCutoff;

		this.subList = new ArrayList<Homolog>();

		for (Homolog hom:list) {
			if ((hom.getBlastHit().getTotalPercentIdentity()/100.0)>idCutoff && hom.getBlastHit().getQueryCoverage()>queryCovCutoff) {
				subList.add(hom);
			}
		}
		
		// finally we update the lookup table
		initialiseMap();
	}
	
	/**
	 * Filters the existing subset to the same domain of life (Bacteria, Archaea, Eukaryota) as the reference sequence
	 */
	public void filterToSameDomainOfLife() {
		Iterator<Homolog> it = subList.iterator();
		while (it.hasNext()) {
			Homolog hom = it.next();
			if (!hom.isUniprot()) {
				LOGGER.info("Removing Uniparc homolog "+hom.getIdentifier()+" as no taxonomy info available for Uniparc");
				it.remove();
				continue;
			}
			if (!hom.getUnirefEntry().isInSameDomainOfLife(this.ref)) {
				it.remove();
			}
		}
	}
	
	/**
	 * Removes the redundant (duplicate) sequences in the filtered subset list of Homologs (those remaining after 
	 * calling {@link #filterToMinIdAndCoverage(double, double)}. 
	 * Redundant sequences are any 2 sequences in the homologs list that are 100% identical.
	 * The algorithm for duplicates elimination is based in the fact that any 2 identical sequences 
	 * in the list must have the same id to the query. Thus we can first group them by id to query 
	 * and then calculate a pairwise id matrix within the groups.
	 * It proceeds as follows:
	 * 1) It groups the sequences by sequence id (in percents with 3 decimal figures) to query
	 * 2) If any of the groups have more than one member then the pairwise identities within the group are calculated 
	 *    (all vs all Needleman-Wunsch). From the pairwise identity matrix sequences that are not 100% identity to all the others
	 *    are removed from the group, leaving groups that contain only identical sequences from same species. 
	 * 3) If after this second pruning any group has more than 1 member then a single member is chosen (first one) 
	 */
	public void removeRedundancy() {
		
		// 1) grouping by sequence identity
		Map<String,List<Homolog>> groups = new HashMap<String,List<Homolog>>();
		for (Homolog hom:subList){
			String key = String.format("%6.3f",hom.getPercentIdentity());
			if (groups.containsKey(key)) {
				groups.get(key).add(hom);
			} else {
				List<Homolog> list = new ArrayList<Homolog>();
				list.add(hom);
				groups.put(key, list);
			}
		}
		LOGGER.debug("Number of protein sequence groups for redundancy elimination (based on same identity value): "+groups.size());
		// 2) finding if group members are really identical (all vs all pairwise alignments)
		for (String key:groups.keySet()) {
			// all vs all pairwise alignments
			List<Homolog> list = groups.get(key);
			LOGGER.debug("Size of group "+key+": "+list.size());
			if (list.size()>1) {
				double[][] pairwiseIdMatrix = new double[list.size()][list.size()];
				for (int i=0;i<list.size();i++){
					for (int j=i+1;j<list.size();j++){
						try {
							PairwiseSequenceAlignment aln = new PairwiseSequenceAlignment(list.get(i).getUnirefEntry().getSeq(), list.get(j).getUnirefEntry().getSeq());
							pairwiseIdMatrix[i][j]=aln.getPercentIdentity();
						} catch (PairwiseSequenceAlignmentException e) {
							LOGGER.error("Unexpected error. Couldn't align sequences "+list.get(i).getIdentifier()+" and "+list.get(j).getIdentifier()+" for redundancy removal procedure");
							LOGGER.error(e.getMessage());
						} catch (OutOfMemoryError e) {
							// this happens when very long sequences are used (e.g. 1tki which is a subdomain of the muscular titin protein ~30000 res!)
							// we can continue here, the only effect is that the entries won't be removed
							LOGGER.error("Out of memory while trying to align sequences "+list.get(i).getIdentifier()+" and "+list.get(j).getIdentifier()+" for redundancy removal procedure. Sequences probably too long");
						}
					}
				}
				// mirroring the other side of the matrix
				for (int i=0;i<list.size();i++){
					for (int j=i+1;j<list.size();j++){
						pairwiseIdMatrix[j][i] = pairwiseIdMatrix[i][j]; 
					}
				}

				String matStr = "";
				// if a member is not identical to all the others then we throw it away
				boolean[] hasNoIdenticalPair = new boolean[list.size()];
				for (int i=0;i<list.size()-1;i++){
					boolean hasIdenticalPair = false;
					for (int j=0;j<list.size();j++){
						matStr+=String.format("%5.1f ", pairwiseIdMatrix[i][j]);
						if (pairwiseIdMatrix[i][j]>99.99f) {
							hasIdenticalPair = true;
						}
					}
					if (!hasIdenticalPair) hasNoIdenticalPair[i] = true;
					matStr+="\n";
				}

				LOGGER.debug("Pairwise similarities: \n"+matStr);
				Iterator<Homolog> it = list.iterator();
				int i = 0;
				while (it.hasNext()){
					it.next();
					if (hasNoIdenticalPair[i]) {
						it.remove();
					}
					i++;
				}
				LOGGER.debug("Size of group after elimination of sequences not identical to any of the others: "+list.size());
			}
		}
		// 3) if there are still groups with size>1 then they are really redundant, all sequences except one have to be eliminated
		List<Homolog> toRemove = new ArrayList<Homolog>();
		for (String key:groups.keySet()) {
			List<Homolog> list = groups.get(key);
			if (list.size()>1) {				
				// we remove all but first
				for (int i=1;i<list.size();i++) {
					toRemove.add(list.get(i));
				}
			}
		}
		for (Homolog hom:toRemove){
			this.subList.remove(hom);
			LOGGER.info("Homolog "+hom.getIdentifier()+" removed because it is redundant (another 100% identical sequence exists in the homologs list).");
		}
		if (!toRemove.isEmpty()) {
			LOGGER.info("Number of homologs after redundancy (duplicates) elimination: "+this.getSizeFilteredSubset());
			// finally we update the lookup table
			initialiseMap();
		}
		
	}
	
	/**
	 * Removes from the filtered subset homologs list (those remaining after calling {@link #filterToMinIdAndCoverage(double, double)})
	 * those homologs that are 100% identical to query while covering at least the given minQueryCov. 
	 * The query itself will usually be in this group and removed with this procedure.
	 */
	public void removeIdenticalToQuery(double minQueryCov) {
		boolean identicalsFound = false;
		Iterator<Homolog> it = subList.iterator();
		while (it.hasNext()) {
			Homolog hom = it.next();
			if (hom.getPercentIdentity()>99.99999 && hom.getBlastHit().getQueryCoverage()>minQueryCov) {
				it.remove();
				identicalsFound = true;
				LOGGER.info("Removing "+hom.getIdentifier()+" because it is 100% identical and covers "+
						String.format("%5.1f%%",hom.getBlastHit().getQueryCoverage()*100.0)+" of the query.");				
			}
		}
		// finally we update the lookup table if something changed
		if (identicalsFound) initialiseMap();
	}
	
	/**
	 * Reduces the size of the subset of homologs by skimming it for homologs with same identities 
	 * to query until maxDesiredHomologs is reached.
	 * The procedure is as follows: sequences are grouped by identity to query (rounded to integer percent values) 
	 * and one of each group eliminated at each iteration until the desired number of homologs is reached.
	 * @param maxDesiredHomologs
	 */
	public void skimList(int maxDesiredHomologs) {
		if (subList.size()<=maxDesiredHomologs) {
			return;
		}
		LOGGER.info("List of homologs too long: "+subList.size()+", skimming it.");
		// 1) grouping by sequence identity
		// we make it a treemap simply to have them ordered by ids in the log output
		Map<Integer,List<Homolog>> groups = new TreeMap<Integer,List<Homolog>>(); 
		for (Homolog hom:subList){
			double percentId = hom.getPercentIdentity();
			int key = (int)Math.round(percentId);
			
			if (groups.containsKey(key)) {
				groups.get(key).add(hom);
			} else {
				List<Homolog> list = new ArrayList<Homolog>();
				list.add(hom);
				groups.put(key, list);
			}
		}
		// 2) skimming iteratively
		int countIterations = 0;
		outer:
		while (true) {
			countIterations++;
			LOGGER.info("Skimming round "+countIterations);
			for (int key:groups.keySet()) {
				List<Homolog> group = groups.get(key);
				if (group.size()>1) {
					// remove the last element of the group
					Homolog toRemove = group.get(group.size()-1);
					
					subList.remove(toRemove);					
					group.remove(toRemove);
					
					LOGGER.info("Removed "+toRemove.getIdentifier()+" ("+String.format("%4.1f",toRemove.getPercentIdentity())+"% id to query)");
					if (subList.size()<=maxDesiredHomologs) break outer;
				}
			}
		}
		LOGGER.info("Size of homolog list after skimming: "+subList.size()+" ("+countIterations+" iterations)");

		// and we reinitialise the lookup maps
		initialiseMap();
	}
	
	/**
	 * Returns the number of homologs in the unfiltered list (not filtered)
	 * @return
	 */
	public int getSizeFullList() {
		return list.size();
	}
	
	/**
	 * Returns the number of homologs in the filtered subset i.e. the one 
	 * after calling {@link #filterToMinIdAndCoverage(double, double)}
	 * @return
	 */
	public int getSizeFilteredSubset() {
		return subList.size();
	}
	
	/**
	 * Returns the sequence identity cutoff, see {@link #filterToMinIdAndCoverage(double, double)}
	 * @return
	 */
	public double getIdCutoff() {
		return idCutoff;
	}
	
	/**
	 * Returns the query coverage cutoff, see {@link #filterToMinIdAndCoverage(double, double)}
	 * @return
	 */
	public double getQCovCutoff() {
		return qCoverageCutoff;
	}
	
	/**
	 * Gets the uniprot version used for blasting
	 * @return
	 */
	public String getUniprotVer() {
		return uniprotVer;
	}
	
	/**
	 * Compute the sequence entropies for all reference sequence (uniprot) positions
	 * @param reducedAlphabet
	 */
	public void computeEntropies(int reducedAlphabet) {
		this.reducedAlphabet = reducedAlphabet;
		this.entropies = new ArrayList<Double>(); 
		for (int i=0;i<refInterval.getLength();i++){
			entropies.add(this.aln.getColumnEntropy(this.aln.seq2al(ref.getUniId(),i+1), reducedAlphabet));
		}
	}
	
	public List<Double> getEntropies() {
		return entropies;
	}
	
	public int getReducedAlphabet() {
		return reducedAlphabet;
	}
	
	public void setUseUniparc(boolean useUniparc) {
		this.useUniparc = useUniparc;
	}
	
	/**
	 * Returns a string containing information about the distribution of the sequence entropies for
	 * the alignment of this homolog list, including an entropy value for the distribution.
	 * It divides the sequence entropy values into 6 bins from 0 to max(s), using that distribution 
	 * to calculate an entropy value of it. 
	 * @return
	 */
	public String getAlnVariability() {
		int numBins = 6;
		// the max value for entropy given a reducedAlphabet
		double maxs = Math.log(reducedAlphabet)/Math.log(2);
		double max = Math.max(1,Collections.max(entropies));
		// the bin step given the max and numBins
		double binStep = max/(double)numBins;
		int[] binCounts = new int[numBins];
		int totalLength = entropies.size();
		for (double s:entropies) {
			int i=0;
			for (double binBoundary=binStep;binBoundary<=max;binBoundary+=binStep) {
				if (s>=binBoundary-binStep && s<binBoundary) binCounts[i]++;
				i++;
			}
		}
		
		StringBuffer bf = new StringBuffer();
		
		bf.append("Distribution of entropies: \n");
		for (double binBoundary=binStep;binBoundary<=max;binBoundary+=binStep) {
			bf.append(String.format("%4.2f ", binBoundary));
		}
		bf.append("\n");
		double sumplogp=0.0;
		for (int i=0;i<numBins;i++){
			bf.append(String.format("%4s ",binCounts[i]));
			double prob = (double)binCounts[i]/(double)totalLength; 
			if (prob!=0){ // plogp is defined to be 0 when p=0 (because of limit). If we let java calculate it, it gives NaN (-infinite) because it tries to compute log(0) 
				sumplogp += prob*(Math.log(prob)/Math.log(2));
			}
		}
		double alnVariability = (-1.0)*sumplogp;
		
		bf.append("\n");
		bf.append(String.format("min: %4.2f max: %4.2f max s possible: %4.2f\n",Collections.min(entropies),Collections.max(entropies),maxs));
		bf.append(String.format("Alignment information content: %4.2f\n",alnVariability));
		return bf.toString();
	}
}
