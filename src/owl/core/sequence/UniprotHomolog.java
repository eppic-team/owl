package owl.core.sequence;

import java.util.List;

import owl.core.runners.blast.BlastHit;

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
	private List<String> taxIds; // many taxonomy ids for a uniprot entry should be a very pathological case. Should be safe to use the first one always
	private List<String> emblCdsIds;
	private List<Sequence> emblCdsSeqs; 
	
	
	
	public UniprotHomolog(BlastHit blastHit, String uniId) {
		this.blastHit = blastHit;
		this.uniId = uniId;
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
	
}
