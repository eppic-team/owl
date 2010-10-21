package owl.core.sequence;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
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

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.xml.sax.SAXException;

import owl.core.connections.EmblWSDBfetchConnection;
import owl.core.connections.NoMatchFoundException;
import owl.core.connections.UniProtConnection;
import owl.core.runners.SelectonRunner;
import owl.core.runners.TcoffeeError;
import owl.core.runners.TcoffeeRunner;
import owl.core.runners.blast.BlastError;
import owl.core.runners.blast.BlastHit;
import owl.core.runners.blast.BlastHitList;
import owl.core.runners.blast.BlastRunner;
import owl.core.runners.blast.BlastXMLParser;
import owl.core.sequence.alignment.AlignmentConstructionError;
import owl.core.sequence.alignment.MultipleSequenceAlignment;
import owl.core.sequence.alignment.PairwiseSequenceAlignment;
import owl.core.sequence.alignment.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;
import owl.core.util.FileFormatError;
import owl.core.util.Goodies;
import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseType;
import uk.ac.ebi.kraken.interfaces.uniprot.NcbiTaxon;
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
	
	private static final Log LOGGER = LogFactory.getLog(UniprotHomologList.class);

	
	/*-------------------------- members --------------------------*/
	
	private UniprotEntry ref;						 // the uniprot entry to which the homologs refer
	private List<UniprotHomolog> list; 				 // the list of homologs
	private Map<String,List<UniprotHomolog>> lookup; // to speed up searches (uniprot ids to Homologs lists) 
													 // it's a list because blast can hit a single uniprot in multiple regions (for us
													 // that's multiple BlastHits)
	private double idCutoff; 						 // the identity cutoff (see restrictToMinId() )
	private String uniprotVer;						 // the version of uniprot used in blasting, read from the reldate.txt uniprot file
	
	private MultipleSequenceAlignment aln;	  		// the protein sequences alignment
	private MultipleSequenceAlignment nucAln; 		// the nucleotides alignment

	private int reducedAlphabet;					// the reduced alphabet used to calculate entropies
	private List<Double> entropies;					// entropies for each uniprot reference sequence position
	private List<Double> kaksRatios;				// ka/ks for each uniprot reference sequence position

	private boolean haveCDSData;
	
	public UniprotHomologList(UniprotEntry ref) {
		this.ref = ref;
		this.idCutoff = 0.0; // i.e. no filter
		haveCDSData = false;
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
	 * @throws BlastError
	 * @throws UniprotVerMisMatchException 
	 */
	public void searchWithBlast(String blastBinDir, String blastDbDir, String blastDb, int blastNumThreads, File cacheFile) throws IOException, BlastError, UniprotVerMisMatchException {
		File outBlast = null;
		boolean fromCache = false;
		if (cacheFile!=null && cacheFile.exists()) {
			outBlast = cacheFile;
			fromCache = true;
			LOGGER.warn("Reading blast results from cache file "+cacheFile);
		} else {
			outBlast = File.createTempFile(BLAST_BASENAME,BLASTOUT_SUFFIX);
			File inputSeqFile = File.createTempFile(BLAST_BASENAME,FASTA_SUFFIX);
			if (!DEBUG) {
				outBlast.deleteOnExit();
				inputSeqFile.deleteOnExit();
			}
			// NOTE: we blast the reference uniprot sequence
			this.ref.getUniprotSeq().writeToFastaFile(inputSeqFile);
			BlastRunner blastRunner = new BlastRunner(blastBinDir, blastDbDir);
			blastRunner.runBlastp(inputSeqFile, blastDb, outBlast, BLAST_OUTPUT_TYPE, BLAST_NO_FILTERING, blastNumThreads);
			this.uniprotVer = readUniprotVer(blastDbDir);
			if (cacheFile!=null) {
				try {
					Goodies.copyFile(outBlast, cacheFile);
				} catch (IOException e) {
					System.err.println("Couldn't write the blast cache file "+cacheFile);
					System.err.println(e.getMessage());
				}
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
			if (!blastList.getQueryId().equals(this.ref.getUniprotSeq().getName())) {
				throw new IOException("Query id "+blastList.getQueryId()+" from cache file "+cacheFile+" does not match the id from the sequence: "+this.ref.getUniprotSeq().getName());
			}
			this.uniprotVer = readUniprotVer(cacheFile.getParent());
			String uniprotVerFromBlastDbDir = readUniprotVer(blastDbDir);
			if (!uniprotVerFromBlastDbDir.equals(uniprotVer)) {
				throw new UniprotVerMisMatchException("Uniprot version from blast db dir "+blastDbDir+
						" ("+uniprotVerFromBlastDbDir+") does not match version in cache dir "+cacheFile.getParent()+" ("+uniprotVer+")");
			}
		}

		
		this.list = new ArrayList<UniprotHomolog>();
		for (BlastHit hit:blastList) {
			String sid = hit.getSubjectId();
			Matcher m = Sequence.DEFLINE_PRIM_ACCESSION_REGEX.matcher(sid);
			if (m.matches()) {
				String uniId = m.group(1);
				list.add(new UniprotHomolog(hit,new UniprotEntry(uniId)));
			} else {
				System.err.println("Could not find uniprot id in subject id "+sid);
			}
		}
		initialiseMap();
	}
	
	//TODO write a searchithPSIBlast method

	/**
	 * Initialises the lookup map (for speeding up lookups of homologs by uniprot ids)
	 */
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
	 * @throws UniprotVerMisMatchException 
	 */
	public void retrieveUniprotKBData() throws UniprotVerMisMatchException {
		UniProtConnection uniprotConn = new UniProtConnection();
		if (!uniprotConn.getVersion().equals(this.uniprotVer)){
			throw new UniprotVerMisMatchException("Uniprot version used for blast ("+uniprotVer+") and uniprot version being queried with api ("+uniprotConn.getVersion()+") don't match!");
		}
		List<String> uniprotIds = new ArrayList<String>();
		for (UniprotHomolog hom:this) {
			uniprotIds.add(hom.getUniId());
		}
		EntryIterator<UniProtEntry> entries = uniprotConn.getMultipleEntries(uniprotIds);

		for (UniProtEntry entry:entries) {
			List<UniprotHomolog> homs = this.getHomolog(entry.getPrimaryUniProtAccession().getValue());
			for (UniprotHomolog hom:homs) {
				hom.getUniprotEntry().setUniprotSeq(new Sequence(hom.getUniId(),entry.getSequence().getValue()));
				
				List<NcbiTaxonomyId> ncbiTaxIds = entry.getNcbiTaxonomyIds();
				if (ncbiTaxIds.size()>1) {
					LOGGER.warn("More than one taxonomy id for uniprot entry "+hom.getUniId());
				}
				hom.getUniprotEntry().setTaxId(ncbiTaxIds.get(0).getValue());
				List<String> taxons = new ArrayList<String>();
				for(NcbiTaxon ncbiTaxon:entry.getTaxonomy()) {
					taxons.add(ncbiTaxon.getValue());
				}
				hom.getUniprotEntry().setTaxons(taxons);
				
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
				hom.getUniprotEntry().setEmblCdsIds(emblCdsIds);
				
				List<Organelle> orglls = entry.getOrganelles();
				if (orglls.size()>0) {
					hom.getUniprotEntry().setGeneEncodingOrganelle(orglls.get(0).getType().getValue());
					if (orglls.size()>1) {
						for (Organelle orgll:orglls){ 
							if (!orgll.getType().equals(hom.getUniprotEntry().getGeneEncodingOrganelle())) {
								LOGGER.warn("Different gene encoding organelles for Uniprot "+hom.getUniId());
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
			allIds.addAll(hom.getUniprotEntry().getEmblCdsIds());
		}
		List<Sequence> allSeqs = null;
		try {
			allSeqs = EmblWSDBfetchConnection.fetchEMBLCDS(allIds, cacheFile);
		} catch (NoMatchFoundException e) {
			// this is unlikely to happen here, that's why we don't write a better error message
			LOGGER.warn("Couldn't retrieve EMBL CDS sequences for some EMBL cds ids");
			LOGGER.warn(e.getMessage());
		}
		
		// we put the list (containing all the sequences from all the homologs) in a lookup table
		// so that we can then retrieve the ones corresponding to each homolog below
		Map<String,Sequence> lookup = new HashMap<String,Sequence>();
		for (Sequence seq:allSeqs) {
			lookup.put(seq.getSecondaryAccession(), seq);
		}
		 
		for (UniprotHomolog hom:this) {
			List<Sequence> seqs = new ArrayList<Sequence>();
			for (String emblCdsId:hom.getUniprotEntry().getEmblCdsIds()) {
				if (lookup.containsKey(emblCdsId)) {
					seqs.add(lookup.get(emblCdsId));
				} else { 
					// this will happen when the CDS sequence was not returned by embl dbfetch or the cache file does not have it (it's in the list of missing entries)
					// in either case we don't want the list of emblcs sequences to contain a null
					LOGGER.warn("Sequence for EMBL CDS "+emblCdsId+" of uniprot entry "+hom.getUniId()+" could not be found. Not using it.");
					//TODO should we also remove from the list of embl cds ids this emblCdsId? Not sure if it can cause problems
					// to have an embl cds id without its corresponding sequence
				}
			}
			
			hom.getUniprotEntry().setEmblCdsSeqs(seqs);
		}		
		haveCDSData = true;
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
	 * Write to the given file the query protein sequence and all the homolog protein 
	 * sequences (full uniprot sequences) in fasta format.
	 * @param outFile
	 * @param writeQuery if true the query sequence is written as well, if false only homologs 
	 * @throws FileNotFoundException
	 */
	public void writeToFasta(File outFile, boolean writeQuery) throws FileNotFoundException {
		PrintWriter pw = new PrintWriter(outFile);
		
		int len = 80;

		if (writeQuery) {
			pw.println(MultipleSequenceAlignment.FASTAHEADER_CHAR + this.ref.getUniprotSeq().getName());
			for(int i=0; i<this.ref.getLength(); i+=len) {
				pw.println(ref.getUniprotSeq().getSeq().substring(i, Math.min(i+len,ref.getLength())));
			}
		}
		
		for(UniprotHomolog hom:this) {
			
			String sequence = hom.getUniprotSeq().getSeq();
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
	 * @throws IOException
	 * @throws TcoffeeError 
	 */
	public void computeTcoffeeAlignment(File tcoffeeBin, boolean veryFast) throws IOException, TcoffeeError {
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
		this.writeToFasta(homologSeqsFile, true);
		TcoffeeRunner tcr = new TcoffeeRunner(tcoffeeBin);
		tcr.runTcoffee(homologSeqsFile, alnFile, TCOFFEE_ALN_OUTFORMAT, outTreeFile, null, tcoffeeLogFile, veryFast);
		
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
		
	}
	
	/**
	 * Returns the protein sequence alignment of query sequence and all homologs
	 * @return
	 */
	public MultipleSequenceAlignment getAlignment() {
		return aln;
	}
	
	/**
	 * Returns a multiple sequence alignment of all valid nucleotide CDS sequences from the 
	 * UniprotHomologList plus the query's CDS by mapping the nucleotide sequences to the protein 
	 * sequences alignment. 
	 * @return
	 */
	public MultipleSequenceAlignment getNucleotideAlignment() {
		if (nucAln!=null) {
			return nucAln;
		}
		
		// first we gather all the aligned protein sequences and the unaligned nucleotide sequences
		List<String> alignedProtSeqs = new ArrayList<String>();
		List<ProteinToCDSMatch> protToCDSMatches = new ArrayList<ProteinToCDSMatch>();
		// query
		ProteinToCDSMatch queryMatching = ref.getRepresentativeCDS();
		if (queryMatching!=null) {
			alignedProtSeqs.add(aln.getAlignedSequence(ref.getUniprotSeq().getName()));
			protToCDSMatches.add(queryMatching);
		}
		// homologs
		for (UniprotHomolog hom:this) {
			ProteinToCDSMatch matching = hom.getUniprotEntry().getRepresentativeCDS();
			if (matching!=null){
				alignedProtSeqs.add(aln.getAlignedSequence(hom.getBlastHit().getSubjectId()));
				protToCDSMatches.add(matching);
			}
		}
		
		// put the gaps in and produce the final nucleotide alignment
		List<Sequence> allSeqs = new ArrayList<Sequence>();

		for (int seqIdx=0;seqIdx<alignedProtSeqs.size();seqIdx++) {
			String alignedProtSeq = alignedProtSeqs.get(seqIdx);
			ProteinToCDSMatch protToCDSMatch = protToCDSMatches.get(seqIdx);
			String unalignedNucSeq = protToCDSMatch.getNucleotideSeqForBestTranslation();
			
			StringBuffer nucSeqSB = new StringBuffer();
			int j = 0;
			for (int i=0;i<alignedProtSeq.length();i++) {
				char aa = alignedProtSeq.charAt(i);
				if (aa==MultipleSequenceAlignment.GAPCHARACTER) {
					nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
					nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
					nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
				} else {
					if (j+3<=unalignedNucSeq.length()) {
						nucSeqSB.append(unalignedNucSeq.substring(j, j+3));
					} else {
						nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
						nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
						nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
					}
					j+=3;
				}				
			}
			allSeqs.add(new Sequence(protToCDSMatch.getCDSName(),nucSeqSB.toString()));
		}
		
		try {
			nucAln = new MultipleSequenceAlignment(allSeqs);
		} catch(AlignmentConstructionError e) {
			LOGGER.fatal("Unexpected error while creating the nucleotides alignment");
			LOGGER.fatal(e.getMessage());
			System.exit(1);
		}
		
		return nucAln;
	}
	
	public void writeAlignmentToFile(File alnFile) throws FileNotFoundException {
		aln.writeFasta(new PrintStream(alnFile), 80, true);
	}
 
	public void writeNucleotideAlignmentToFile(File alnFile) throws FileNotFoundException {
		this.getNucleotideAlignment().writeFasta(new PrintStream(alnFile), 80, true);
	}
	
	/**
	 * Removes homologs below the given sequence identity value.
	 * @param idCutoff
	 */
	public void restrictToMinIdAndCoverage(double idCutoff, double queryCovCutoff) {
		this.idCutoff = idCutoff;

		Iterator<UniprotHomolog> it = list.iterator();
		while (it.hasNext()) {
			UniprotHomolog hom = it.next();
			if ((hom.getBlastHit().getPercentIdentity()/100.0)<=idCutoff || hom.getBlastHit().getQueryCoverage()<=queryCovCutoff) {
				it.remove();

			}
		}
		// finally we update the lookup table
		initialiseMap();
	}
	
	/**
	 * Removes the redundant sequences in this list of Homologs. The redundancy reduction procedes as follows:
	 * 1) It groups the sequences by taxonomy id and sequence identity
	 * 2) If any of the groups have more than one member then the pairwise identities within the group are calculated 
	 *    (all vs all Needleman-Wunsch). From the pairwise identity matrix sequences that are not 100% identity to all the others
	 *    are removed from the group, leaving groups that contain only identical sequences from same species. 
	 * 3) If after this second pruning any group has more than 1 member then a single member is chosen (first one with a
	 *    good matching corresponding CDS sequence or simply first one if no good CDS matchings exist) 
	 */
	public void removeRedundancy() {
		
		// 1) grouping by tax id and sequence identity
		Map<String,List<UniprotHomolog>> groups = new HashMap<String,List<UniprotHomolog>>();
		for (UniprotHomolog hom:this){
			double percentId = hom.getPercentIdentity();
			String taxId = hom.getUniprotEntry().getTaxId();
			String key = taxId+"_"+String.format("%5.1f",percentId);
			if (groups.containsKey(key)) {
				groups.get(key).add(hom);
			} else {
				List<UniprotHomolog> list = new ArrayList<UniprotHomolog>();
				list.add(hom);
				groups.put(key, list);
			}
		}
		LOGGER.debug("Number of protein sequence groups for redundancy elimination (based on same identity value and same tax id): "+groups.size());
		// 2) finding if group members are really identical (all vs all pairwise alignments)
		for (String key:groups.keySet()) {
			// all vs all pairwise alignments
			List<UniprotHomolog> list = groups.get(key);
			LOGGER.debug("Size of group "+key+": "+list.size());
			if (list.size()>1) {
				double[][] pairwiseIdMatrix = new double[list.size()][list.size()];
				for (int i=0;i<list.size();i++){
					for (int j=i+1;j<list.size();j++){
						try {
							PairwiseSequenceAlignment aln = new PairwiseSequenceAlignment(list.get(i).getUniprotSeq(), list.get(j).getUniprotSeq());
							pairwiseIdMatrix[i][j]=aln.getPercentIdentity();
						} catch (PairwiseSequenceAlignmentException e) {
							LOGGER.fatal("Unexpected error. Couldn't align sequences for redundancy removal procedure");
							LOGGER.fatal(e.getMessage());
							System.exit(1);
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
				Iterator<UniprotHomolog> it = list.iterator();
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
		//    in any case we first check that the one we want to eliminate doesn't have a better CDS representative
		List<UniprotHomolog> toRemove = new ArrayList<UniprotHomolog>();
		for (String key:groups.keySet()) {
			List<UniprotHomolog> list = groups.get(key);
			if (list.size()>1) {
				UniprotHomolog homToKeep = null;
				if (hasCDSData()) {
					for (UniprotHomolog hom:list){
						if (hom.getUniprotEntry().getRepresentativeCDS()!=null) {
							homToKeep = hom;
							break;
						}
					}
				}
				if (homToKeep==null) { // if there wasn't any good CDS homolog we remove all but first
					for (int i=1;i<list.size();i++) {
						toRemove.add(list.get(i));
					}
				} else { // if we found one good CDS homolog we remove all the others
					for (UniprotHomolog hom:list) {
						if (hom!=homToKeep) toRemove.add(hom);
					}
				}
			}
		}
		for (UniprotHomolog hom:toRemove){
			this.list.remove(hom);
			LOGGER.info("Homolog "+hom.getUniId()+" removed because it is redundant.");
		}
		LOGGER.info("Number of homologs after redundancy elimination: "+this.size());
		
		// finally we update the lookup table
		initialiseMap();
	}
	
	/**
	 * Reduces the number of homologs by skimming the list for homologs with same identities 
	 * until maxDesiredHomologs is reached.
	 * The procedure is as follows: sequences are grouped by identity to query and one of each group 
	 * eliminated (if possible the one without valid CDS matching) at each iteration until the 
	 * desired number of homologs is reached.
	 * This is needed for example to perform ka/ks calculations with selecton (too many sequences 
	 * are far too slow and a risk of ks saturation).
	 * @param maxDesiredHomologs
	 */
	public void skimList(int maxDesiredHomologs) {
		if (list.size()<=maxDesiredHomologs) {
			return;
		}
		LOGGER.info("List of homologs too long: "+list.size()+", skimming it.");
		// 1) grouping by sequence identity
		Map<Integer,List<UniprotHomolog>> groups = new HashMap<Integer,List<UniprotHomolog>>();
		for (UniprotHomolog hom:this){
			double percentId = hom.getPercentIdentity();
			int key =(int) Math.round(percentId);
			if (groups.containsKey(key)) {
				groups.get(key).add(hom);
			} else {
				List<UniprotHomolog> list = new ArrayList<UniprotHomolog>();
				list.add(hom);
				groups.put(key, list);
			}
		}
		// 2) skimming iteratively
		int countIterations = 0;
		outer:
		while (true) {
			countIterations++;
			for (int key:groups.keySet()) {
				List<UniprotHomolog> group = groups.get(key);
				if (group.size()>1) {
					// remove the last element of the group
					UniprotHomolog toRemove = null;
					if (hasCDSData()) {
						toRemove = getHomologNonValidCDS(group);
					}
					if (toRemove==null) {
						toRemove = group.get(group.size()-1);
					} 
					list.remove(toRemove);					
					group.remove(toRemove);
					
					LOGGER.info("Removed "+toRemove.getUniId());
					if (list.size()<=maxDesiredHomologs) break outer;
				}
			}
		}
		LOGGER.info("Size of homolog list after skimming: "+list.size()+" ("+countIterations+" iterations)");

		// and we reinitialise the lookup maps
		initialiseMap();
	}
	
	/**
	 * Given a list of homologs returns the first one that does not have a valid CDS match 
	 * or null if all homologs have valid matches. 
	 * @param list
	 * @return
	 */
	private static UniprotHomolog getHomologNonValidCDS(List<UniprotHomolog> list) {
		for (UniprotHomolog hom:list) {
			if (hom.getUniprotEntry().getRepresentativeCDS()==null) {
				return hom;
			}
		}
		return null;
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
	
	public int getNumHomologsWithCDS() {
		int count = 0;
		for (UniprotHomolog hom:this) {
			if (hom.getUniprotEntry().hasCDS()) {
				count++;
			}
		}				
		return count;
	}
	
	public int getNumHomologsWithValidCDS() {
		int count = 0;
		for (UniprotHomolog hom:this) {
			if (hom.getUniprotEntry().getRepresentativeCDS()!=null) {
				count++;
			}
		}
		return count;
	}

	/**
	 * Tells whether CDS data was retrieved from this list by calling {@link #retrieveEmblCdsSeqs(File)}
	 * @return
	 */
	public boolean hasCDSData() {
		return haveCDSData;
	}
	
	/**
	 * Checks whether this UniprotHomologList contains genes encoded with strictly 
	 * one genetic code type.
	 * @see GeneticCodeType
	 * @return
	 */
	public boolean isConsistentGeneticCodeType() {
		GeneticCodeType lastGct = null;
		for (UniprotHomolog hom:this) {
			GeneticCodeType gct = hom.getUniprotEntry().getGeneticCodeType();
			if (lastGct!=null && !lastGct.equals(gct)) {
				return false;
			}
			lastGct = gct;
		}
		return true;
	}
	
	/**
	 * Gets the genetic code type of this UniprotHomologList. It does it by simply returning
	 * the genetic code type of the first homolog in the list. It does NOT check for consistency
	 * within the list. For that use {@link #isConsistentGeneticCodeType()}
	 * @return
	 */
	public GeneticCodeType getGeneticCodeType() {
		return this.list.get(0).getUniprotEntry().getGeneticCodeType();
	}
	
	/**
	 * Compute the sequence entropies for all reference sequence (uniprot) positions
	 * @param reducedAlphabet
	 */
	public void computeEntropies(int reducedAlphabet) {
		this.reducedAlphabet = reducedAlphabet;
		this.entropies = new ArrayList<Double>(); 
		for (int i=0;i<ref.getUniprotSeq().getLength();i++){
			entropies.add(this.aln.getColumnEntropy(this.aln.seq2al(ref.getUniprotSeq().getName(),i+1), reducedAlphabet));
		}
	}
	
	/**
	 * Compute the sequence ka/ks ratios with selecton for all reference CDS sequence positions
	 * @param selectonBin
	 * @param resultsFile
	 * @param logFile
	 * @param treeFile
	 * @param globalResultsFile
	 * @param epsilon
	 * @throws IOException
	 */
	public void computeKaKsRatiosSelecton(File selectonBin, File resultsFile, File logFile, File treeFile, File globalResultsFile, double epsilon) 
	throws IOException {
		kaksRatios = new ArrayList<Double>();
		SelectonRunner sr = new SelectonRunner(selectonBin);
		if(!resultsFile.exists()) {
			File alnFile = File.createTempFile("selecton.", ".cds.aln");
			alnFile.deleteOnExit();
			this.nucAln.writeFasta(new PrintStream(alnFile), 80, true);
			sr.run(alnFile, resultsFile, logFile, treeFile, null, globalResultsFile, this.ref.getRepresentativeCDS().getCDSName(), this.getGeneticCodeType(),epsilon);
		} else {
			LOGGER.warn("Selecton output file "+resultsFile+" already exists. Using the file instead of running selecton.");
			try {
				sr.parseResultsFile(resultsFile, this.ref.getRepresentativeCDS().getCDSName());
			} catch (FileFormatError e) {
				LOGGER.warn("Cached output selecton file "+resultsFile+" does not seem to be in the righ format");
				LOGGER.warn(e.getMessage());
				LOGGER.warn("Running selecton and overwritting the file.");
				File alnFile = File.createTempFile("selecton.", ".cds.aln");
				alnFile.deleteOnExit();
				this.nucAln.writeFasta(new PrintStream(alnFile), 80, true);
				sr.run(alnFile, resultsFile, logFile, treeFile, null, globalResultsFile, this.ref.getRepresentativeCDS().getCDSName(), this.getGeneticCodeType(),epsilon);				
			}
		}
		kaksRatios = sr.getKaKsRatios();
		if (kaksRatios.size()!=this.ref.getUniprotSeq().getLength()) {
			LOGGER.info("Size of ka/ks ratio list ("+kaksRatios.size()+") is not the same as length of reference sequence ("+this.ref.getUniprotSeq().getLength()+")");
		}
	}

	public List<Double> getEntropies() {
		return entropies;
	}
	
	public List<Double> getKaksRatios() {
		return kaksRatios;
	}
	
	public int getReducedAlphabet() {
		return reducedAlphabet;
	}
	
	/**
	 * Tells whether a given position of the reference sequence (starting at 0)
	 * is reliable with respect to the CDS matching of the reference sequence and the CDS
	 * matchings of the homologs 
	 * (not reliable is considered any position where the CDS translation does not match exactly the 
	 * protein residue) 
	 * @param i
	 * @return
	 */
	public boolean isReferenceSeqPositionReliable(int i) {
		// check if matching CDS of ref sequence is reliable at position i
		if (!this.ref.isReliablePosition(i)) {
			return false;
		}
		// check if each of the matching CDS of the homologs is reliable at position i of the ref sequence
		for (UniprotHomolog hom:this){
			if (hom.getUniprotEntry().getRepresentativeCDS()==null) {
				continue; // the homolog may have no representative CDS, in that case we don't want to check the CDS alignment as it doesn't exist
			}
			int alnPos = aln.seq2al(this.ref.getUniprotSeq().getName(), i+1);
			int seqPos = aln.al2seq(hom.getLongSequenceTag(), alnPos);
			if (seqPos==-1) { // the position maps to a gap in the homolog sequence, there's no CDS alignment to check, we can continue to next homolog directly
				continue;
			}
			if (!hom.getUniprotEntry().isReliablePosition(seqPos-1)) {
				return false;
			}
		}
		return true;
	}
}
