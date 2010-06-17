package owl.core.sequence;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import owl.core.connections.EmblWSDBfetchConnection;
import owl.core.connections.NoMatchFoundException;
import owl.core.connections.UniProtConnection;
import owl.core.runners.blast.BlastHit;
import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseType;
import uk.ac.ebi.kraken.interfaces.uniprot.GeneEncodingType;
import uk.ac.ebi.kraken.interfaces.uniprot.NcbiTaxonomyId;
import uk.ac.ebi.kraken.interfaces.uniprot.Organelle;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.dbx.embl.Embl;

/**
 * Class to encapsulate a Uniprot blast hit representing a homolog to a certain
 * sequence. It stores data from the blast hit, the uniprot identifiers as well 
 * as taxonomy ids and EMBL CDS (DNA coding sequence) ids and sequences.
 *  
 * @see UniprotHomologList
 * 
 * @author duarte_j
 *
 */
public class UniprotHomolog {
	
	
	private BlastHit blastHit;
	private String uniId;
	private Sequence uniproSeq;
	private List<String> taxIds; // many taxonomy ids for a uniprot entry should be a very pathological case. Should be safe to use the first one always
	private List<String> emblCdsIds;
	private List<Sequence> emblCdsSeqs; 
	private GeneEncodingType geneEncodingOrganelle; // the organelle where this gene is encoded (important for genetic code): if gene encoded in nucleus this is null
	
	
	public UniprotHomolog(BlastHit blastHit, String uniId) {
		this.blastHit = blastHit;
		this.uniId = uniId;
	}
	
	public UniprotHomolog(String uniId) {
		this.uniId = uniId;
		this.blastHit = null;
	}

	public BlastHit getBlastHit() {
		return blastHit;
	}

	public void setBlastHit(BlastHit blastHit) {
		this.blastHit = blastHit;
	}

	public String getUniId() {
		return uniId;
	}

	public void setUniId(String uniId) {
		this.uniId = uniId;
	}
	
	public Sequence getUniprotSeq() {
		return uniproSeq;
	}
	
	public void setUniprotSeq(Sequence seq) {
		this.uniproSeq = seq;
	}
	
	public List<String> getTaxIds() {
		return taxIds;
	}

	public void setTaxIds(List<String> taxIds) {
		this.taxIds = taxIds;
	}

	public List<String> getEmblCdsIds() {
		return emblCdsIds;
	}

	public void setEmblCdsIds(List<String> emblCdsIds) {
		this.emblCdsIds = emblCdsIds;
	}

	public List<Sequence> getEmblCdsSeqs() {
		return emblCdsSeqs;
	}

	public void setEmblCdsSeqs(List<Sequence> emblCdsSeqs) {
		this.emblCdsSeqs = emblCdsSeqs;
	}

	public double getPercentIdentity() {
		return this.blastHit.getPercentIdentity();
	}
	
	public GeneEncodingType getGeneEncodingOrganelle() {
		return this.geneEncodingOrganelle;
	}
	
	public void setGeneEncodingOrganelle(GeneEncodingType geneEncodingOrganelle){
		this.geneEncodingOrganelle = geneEncodingOrganelle;
	}
	
	/**
	 * Retrieves from UniprotKB the sequence, taxonomy and EMBL CDS ids data,
	 * by using the remote Uniprot API
	 */
	public void retrieveUniprotKBData() {
		this.taxIds = new ArrayList<String>();
		this.emblCdsIds = new ArrayList<String>();
		Set<String> tmpEmblCdsIdsSet = new TreeSet<String>();
		
		UniProtConnection uniprotConn = new UniProtConnection();
		UniProtEntry entry = null;
		try {
			entry = uniprotConn.getEntry(uniId);
		} catch (NoMatchFoundException e) {
			System.err.println("Warning: couldn't find uniprot id "+uniId+" through Uniprot JAPI");
			return;
		}
		
		this.setUniprotSeq(new Sequence(this.getUniId(),entry.getSequence().getValue()));
		
		for(NcbiTaxonomyId ncbiTaxId:entry.getNcbiTaxonomyIds()) {
			taxIds.add(ncbiTaxId.getValue());
		}
		
		Collection<Embl> emblrefs = entry.getDatabaseCrossReferences(DatabaseType.EMBL);
		for(Embl ref:emblrefs) {
			String emblCdsIdWithVer = ref.getEmblProteinId().getValue();
			if (!emblCdsIdWithVer.equals("-")) { // for non annotated genomic dna cds sequences the identifier is '-', we ignore them
				String emblCdsId = emblCdsIdWithVer.substring(0, emblCdsIdWithVer.lastIndexOf("."));
				//emblCdsIds.add(emblCdsId);
				// to ensure there are no duplicate ids (it happens sometimes) we use a set
				tmpEmblCdsIdsSet.add(emblCdsId);
			}
		}
		this.emblCdsIds.addAll(tmpEmblCdsIdsSet);
		List<Organelle> orglls = entry.getOrganelles();
		if (orglls.size()>0) {
			this.geneEncodingOrganelle = orglls.get(0).getType();
			if (orglls.size()>1) {
				for (Organelle orgll:orglls){ 
					if (!orgll.getType().equals(this.geneEncodingOrganelle)) {
						System.err.println("Warning! Different gene encoding organelles for Uniprot "+this.uniId);
					}
				}
			}
		}
	}
	
	/**
	 * Retrieves from EMBL DB fetch web service the EMBL CDS sequences
	 * @param cacheFile a FASTA file containing the sequences to retrieve. If present and if
	 * it contains the required sequence then it is read from cacheFile. If null or file
	 * does not exist or file older than {@link EmblWSDBfetchConnection.MAX_CACHE_AGE} then the sequences are 
	 * retrieved from EMBL DB fetch 
	 * @throws IOException
	 */
	public void retrieveEmblCdsSeqs(File cacheFile) throws IOException {
		this.emblCdsSeqs = new ArrayList<Sequence>();
		
		try {
			this.emblCdsSeqs = EmblWSDBfetchConnection.fetchEMBLCDS(this.emblCdsIds, cacheFile);
		} catch (NoMatchFoundException e) {
			// this is unlikely to happen here, that's why we don't write a better error message
			System.err.println("Couldn't retrieve EMBL CDS sequences for EMBL cds ids"); 
		}

	}

	public void checkEmblCDSMatching() {
		for (Sequence cds:emblCdsSeqs) {
			if (this.geneEncodingOrganelle!=null) {
				System.err.println("Warning! The cds sequence "+cds.getSecondaryAccession()+" is not encoded in nucleus!");
			}
			ProteinToCDSMatch matching = new ProteinToCDSMatch(this.getUniprotSeq(), cds, GeneticCodeType.STANDARD);
			matching.printSummary(99.9999f);
			
		}		
	}
}
