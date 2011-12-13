package owl.core.sequence;

import java.io.Serializable;

import owl.core.runners.blast.BlastHit;

/**
 * Class to encapsulate a blast hit representing a homolog to a certain
 * sequence. It stores data from the blast hit and the Uniprot entry or other 
 * data if not a Uniprot 
 *  
 * @see HomologList
 * 
 * @author duarte_j
 *
 */
public class Homolog implements Serializable {
	
	
	private static final long serialVersionUID = 1L;

	private BlastHit blastHit;
	private UniprotEntry uniprotEntry;
	private String id;
	private String sequence;
	
	/**
	 * Constructs a new Uniprot Homolog
	 * @param blastHit
	 * @param uniprotEntry
	 */
	public Homolog(BlastHit blastHit, UniprotEntry uniprotEntry) {
		this.blastHit = blastHit;
		this.uniprotEntry = uniprotEntry;
	}
	
	/**
	 * Constructs a new non-uniprot Homolog
	 * @param blastHit
	 * @param id
	 */
	public Homolog(BlastHit blastHit, String id) {
		this.blastHit = blastHit;
		this.id = id;
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
	
	public String getIdentifier() {
		if (uniprotEntry==null) return id;
		else return this.uniprotEntry.getUniId();
	}

	public Sequence getUniprotSeq() {
		return this.uniprotEntry.getUniprotSeq();
	}
	

	public double getPercentIdentity() {
		return this.blastHit.getTotalPercentIdentity();
	}
	
	/**
	 * Returns the long sequence tag i.e. the original Uniprot long identifier (e.g. tr|Q6V461|Q6V461_AERHY)
	 * (without the comment field)
	 * @return
	 */
	public String getLongSequenceTag() {
		return this.blastHit.getSubjectId();
	}
	
	public String getSequence() {
		if (uniprotEntry==null) return sequence;
		else return uniprotEntry.getUniprotSeq().getSeq();
	}
	
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
	
	/**
	 * Returns true if this Homolog is a Uniprot entry, false if a Uniparc
	 * @return
	 */
	public boolean isUniprot() {
		return uniprotEntry!=null;
	}
}
