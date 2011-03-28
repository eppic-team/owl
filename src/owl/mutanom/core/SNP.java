package owl.mutanom.core;

import owl.core.structure.AminoAcid;

/**
 * A SNP in a protein sequence. That is, a germline mutation which is listed in the dbSNP database.
 * @author stehr
 */
public class SNP extends Mutation {

	/*--------------------------- member variables --------------------------*/
	int rsId;			// dbSNP identifier
	double freq;		// minor allele frequency (if known, otherwise NaN)
	
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Creates a new missense SNP.
	 */
	public SNP(AminoAcid before, AminoAcid after, int position, int rsId, double freq) {
		super(before,after, position);
		this.rsId = rsId;
		this.freq = freq;
		this.type = MutType.MISSENSE;
	}
	
	/**
	 * Creates a new SNP of the given type
	 */
	public SNP(AminoAcid before, AminoAcid after, int position, int rsId, double freq, MutType type) {
		super(before,after, position);
		this.rsId = rsId;
		this.freq = freq;
		this.type = type;
	}
	
	/*-------------------------- implemented methods ------------------------*/
	
	public String toString() {
		return "p." + before.getOneLetterCode() + position + after.getOneLetterCode();
	}
	
}
