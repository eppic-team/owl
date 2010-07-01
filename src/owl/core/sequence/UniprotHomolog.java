package owl.core.sequence;

import owl.core.runners.blast.BlastHit;

/**
 * Class to encapsulate a Uniprot blast hit representing a homolog to a certain
 * sequence. It stores data from the blast hit and the uniprot entry 
 *  
 * @see UniprotHomologList
 * 
 * @author duarte_j
 *
 */
public class UniprotHomolog {
	
	
	private BlastHit blastHit;
	private UniprotEntry uniprotEntry;
	
	
	public UniprotHomolog(BlastHit blastHit, UniprotEntry uniprotEntry) {
		this.blastHit = blastHit;
		this.uniprotEntry = uniprotEntry;
	}
	
	public BlastHit getBlastHit() {
		return blastHit;
	}

	public void setBlastHit(BlastHit blastHit) {
		this.blastHit = blastHit;
	}

	public UniprotEntry getUniprotEntry() {
		return this.uniprotEntry;
	}
	
	public void setUniprotEntry(UniprotEntry uniprotEntry) {
		this.uniprotEntry = uniprotEntry;
	}
	
	public String getUniId() {
		return this.uniprotEntry.getUniId();
	}

	public Sequence getUniprotSeq() {
		return this.uniprotEntry.getUniprotSeq();
	}
	

	public double getPercentIdentity() {
		return this.blastHit.getPercentIdentity();
	}
	
	/**
	 * Returns the long sequence tag i.e. the original Uniprot long identifier (e.g. tr|Q6V461|Q6V461_AERHY)
	 * (without the comment field)
	 * @return
	 */
	public String getLongSequenceTag() {
		return this.blastHit.getSubjectId();
	}
}
