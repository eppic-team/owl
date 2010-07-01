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

	
	public UniprotHomologList(UniprotEntry ref) {
		this.ref = ref;
		this.idCutoff = 0.0; // i.e. no filter
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
			// NOTE: we blast the reference uniprot sequence
			this.ref.getUniprotSeq().writeToFastaFile(inputSeqFile);
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
			if (!blastList.getQueryId().equals(this.ref.getUniprotSeq().getName())) {
				throw new IOException("Query id "+blastList.getQueryId()+" from cache file "+cacheFile+" does not match the id from the sequence: "+this.ref.getUniprotSeq().getName());
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
				list.add(new UniprotHomolog(hit,new UniprotEntry(uniId)));
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
				hom.getUniprotEntry().setUniprotSeq(new Sequence(hom.getUniId(),entry.getSequence().getValue()));
				
				List<NcbiTaxonomyId> ncbiTaxIds = entry.getNcbiTaxonomyIds();
				if (ncbiTaxIds.size()>1) {
					System.err.println("Warning! more than one taxonomy id for uniprot entry "+hom.getUniId());
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
			allIds.addAll(hom.getUniprotEntry().getEmblCdsIds());
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
			for (String emblCdsId:hom.getUniprotEntry().getEmblCdsIds()) {
				if (lookup.containsKey(emblCdsId)) {
					seqs.add(lookup.get(emblCdsId));
				} else { 
					// this will happen when the CDS sequence was not returned by embl dbfetch or the cache file does not have it (it's in the list of missing entries)
					// in either case we don't want the list of emblcs sequences to contain a null
					System.err.println("Warning! Sequence for EMBL CDS "+emblCdsId+" of uniprot entry "+hom.getUniId()+" could not be found. Not using it.");
					//TODO should we also remove from the list of embl cds ids this emblCdsId? Not sure if it can cause problems
					// to have an embl cds id without its corresponding sequence
				}
			}
			
			hom.getUniprotEntry().setEmblCdsSeqs(seqs);
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
	 * Write to the given file the query protein sequence and all the homolog protein 
	 * sequences (full uniprot sequences) in fasta format.
	 * @param outFile
	 * @throws FileNotFoundException
	 */
	public void writeToFasta(File outFile) throws FileNotFoundException {
		PrintWriter pw = new PrintWriter(outFile);
		
		int len = 80;

		pw.println(MultipleSequenceAlignment.FASTAHEADER_CHAR + this.ref.getUniprotSeq().getName());
		for(int i=0; i<this.ref.getLength(); i+=len) {
			pw.println(ref.getUniprotSeq().getSeq().substring(i, Math.min(i+len,ref.getLength())));
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
		this.writeToFasta(homologSeqsFile);
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
		List<Sequence> allSeqs = new ArrayList<Sequence>();
		
		// CDS of the query sequence
		StringBuffer nucSeqSB = new StringBuffer();
		String queryProtSeq = aln.getAlignedSequence(ref.getUniprotSeq().getName());
		ProteinToCDSMatch queryMatching = ref.getRepresentativeCDS();
		if (queryMatching!=null) {
			String bestNucSeq = queryMatching.getNucleotideSeqForBestTranslation();
			int j = 0;
			for (int i=0;i<queryProtSeq.length();i++) {
				char aa = queryProtSeq.charAt(i);
				if (aa==MultipleSequenceAlignment.GAPCHARACTER) {
					nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
					nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
					nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
				} else {
					if (j+3<bestNucSeq.length()) {
						nucSeqSB.append(bestNucSeq.substring(j, j+3));
					} else {
						nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
						nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
						nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
					}
					j+=3;
				}
			}
			allSeqs.add(new Sequence(queryMatching.getCDSName(),nucSeqSB.toString()));
		}
		
		for (UniprotHomolog hom:this) {
			nucSeqSB = new StringBuffer();
			String protSeq = aln.getAlignedSequence(hom.getBlastHit().getSubjectId());
			ProteinToCDSMatch matching = hom.getUniprotEntry().getRepresentativeCDS();
			if (matching!=null) {
				String bestNucSeq = matching.getNucleotideSeqForBestTranslation();
				int j = 0;
				for (int i=0;i<protSeq.length();i++) {
					char aa = protSeq.charAt(i);
					if (aa==MultipleSequenceAlignment.GAPCHARACTER) {
						nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
						nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
						nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
					} else {
						if (j+3<bestNucSeq.length()) {
							nucSeqSB.append(bestNucSeq.substring(j, j+3));
						} else {
							nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
							nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
							nucSeqSB.append(MultipleSequenceAlignment.GAPCHARACTER);
						}
						j+=3;
						
					}
				}
				allSeqs.add(new Sequence(matching.getCDSName(),nucSeqSB.toString()));
			}
		}
		try {
			nucAln = new MultipleSequenceAlignment(allSeqs);
		} catch(AlignmentConstructionError e) {
			System.err.println("Unexpected error while creating the nucleotides alignment");
			System.err.println(e.getMessage());
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
			System.out.println("Selecton output file "+resultsFile+" already exists. Using the file instead of running selecton.");
			try {
				sr.parseResultsFile(resultsFile, this.ref.getRepresentativeCDS().getCDSName());
			} catch (FileFormatError e) {
				System.err.println("Warning! cached output selecton file "+resultsFile+" does not seem to be in the righ format");
				System.err.println(e.getMessage());
				System.err.println("Running selecton and overwritting the file.");
				File alnFile = File.createTempFile("selecton.", ".cds.aln");
				alnFile.deleteOnExit();
				this.nucAln.writeFasta(new PrintStream(alnFile), 80, true);
				sr.run(alnFile, resultsFile, logFile, treeFile, null, globalResultsFile, this.ref.getRepresentativeCDS().getCDSName(), this.getGeneticCodeType(),epsilon);				
			}
		}
		kaksRatios = sr.getKaKsRatios();
		if (kaksRatios.size()!=this.ref.getUniprotSeq().getLength()) {
			System.err.println("Warning! Size of ka/ks ratio list ("+kaksRatios.size()+") is not the same as length of reference sequence ("+this.ref.getUniprotSeq().getLength()+")");
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
			int alnPos = aln.seq2al(this.ref.getUniprotSeq().getName(), i+1);
			int seqPos = aln.al2seq(hom.getLongSequenceTag(), alnPos);
			if (seqPos==-1) { // the position maps to a gap in the homolog sequence, there's no CDS alignment to check, we can continue to next homolog directly
				continue;
			}
			if (hom.getUniprotEntry().getRepresentativeCDS()==null) {
				continue; // the homolog may have no representative CDS, in that case we don't want to check the CDS alignment as it doesn't exist
			}
			if (!hom.getUniprotEntry().isReliablePosition(seqPos-1)) {
				return false;
			}
		}
		return true;
	}
}
