package owl.core.sequence;

import java.io.Serializable;

import owl.core.runners.blast.BlastHit;

/**
 * Class to encapsulate a blast hit representing a homolog to a certain
 * sequence. 
 * It stores data from the blast hit and the UniRef entry (can be UniProt or UniParc)
 *  
 * @see HomologList
 * 
 * @author duarte_j
 *
 */
public class Homolog implements Serializable {
	
	
	private static final long serialVersionUID = 1L;

	private BlastHit blastHit;
	private UnirefEntry unirefEntry;
	
	/**
	 * Constructs a new UniRef Homolog
	 * @param blastHit
	 * @param uniprotEntry
	 */
	public Homolog(BlastHit blastHit, UnirefEntry unirefEntry) {
		this.blastHit = blastHit;
		this.unirefEntry = unirefEntry;
	}
	
	public BlastHit getBlastHit() {
		return blastHit;
	}

	public void setBlastHit(BlastHit blastHit) {
		this.blastHit = blastHit;
	}

	public UnirefEntry getUnirefEntry() {
		return this.unirefEntry;
	}
	
//	public UniprotEntry getUniprotEntry() {
//		return this.unirefEntry.getUniprotEntry();
//	}

	/**
	 * Returns the UniProt identifier if the UniRef entry corresponds to a UniProt or
	 * else it returns the UniParc identifier 
	 * @return
	 */
	public String getIdentifier() {
		if (unirefEntry.isUniprot()) return unirefEntry.getUniprotId();
		else return unirefEntry.getUniparcId();
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
		return unirefEntry.getSequence();
	}
	
	/**
	 * Returns true if this Homolog is a Uniprot entry, false if a Uniparc
	 * @return
	 */
	public boolean isUniprot() {
		return unirefEntry.isUniprot();
	}
}
