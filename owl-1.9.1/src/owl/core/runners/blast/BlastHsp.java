package owl.core.runners.blast;

import java.io.Serializable;

import owl.core.sequence.alignment.AlignmentConstructionException;
import owl.core.sequence.alignment.MultipleSequenceAlignment;

public class BlastHsp implements Serializable {

	private static final long serialVersionUID = 1L;

	private BlastHit parent;
	private int aliLength;
	private int queryStart;
	private int queryEnd;
	private int subjectStart;
	private int subjectEnd;
	private double eValue;
	private double score;
	private double percentIdentity;

	private MultipleSequenceAlignment al; 		// if BlastHit create from tabular file parsing, this will be null 

	/**
	 * Constructor to be used when parsing from XML output, use the setters to 
	 * fill the values
	 */
	public BlastHsp(BlastHit hit) {
		this.parent = hit;
	}
	
	/**
	 * Constructor to be used when parsing from tabular output
	 * @param field0
	 * @param field1
	 * @param field2
	 * @param field3
	 * @param field4
	 * @param field5
	 * @param field6
	 * @param field7
	 * @param field8
	 * @param field9
	 * @param field10
	 * @param field11
	 */
	public BlastHsp(String field0, String field1, String field2, String field3, String field4, String field5, 
			String field6, String field7, String field8, String field9, String field10, String field11) {
		this.percentIdentity = Double.parseDouble(field2.trim());
		this.aliLength = Integer.parseInt(field3.trim());
		//this.mismatches = Integer.parseInt(field4.trim());
		//this.gapOpenings = Integer.parseInt(field5.trim());
		this.queryStart = Integer.parseInt(field6.trim());
		this.queryEnd = Integer.parseInt(field7.trim());
		this.subjectStart = Integer.parseInt(field8.trim());
		this.subjectEnd = Integer.parseInt(field9.trim());
		this.eValue = Double.parseDouble(field10.trim());
		this.score = Double.parseDouble(field11.trim());
		this.al = null;		
	}
	
	/**
	 * Returns the BlastHit parent of this hsp.
	 * @return
	 */
	public BlastHit getParent() {
		return parent;
	}
	
	/**
	 * Sets the BlastHit parent of this hsp.
	 * @param parent
	 */
	public void setParent(BlastHit parent) {
		this.parent = parent;
	}
	
	/**
	 * Returns the score for this hsp.
	 * @return
	 */
	public double getScore() {
		return score;
	}
	
	/**
	 * Sets the score for this hsp.
	 * @param score
	 */
	public void setScore(double score) {
		this.score = score;
	}
	
	/** 
	 * Returns the e-value for this hsp.
	 * @return
	 */
	public double getEValue() {
		return this.eValue;
	}
	
	/**
	 * Returns the alignment length of this hsp.
	 * @return
	 */
	public int getAliLength() {
		return aliLength;
	}

	/**
	 * Returns the subject start position of this hsp.
	 * @return
	 */
	public int getSubjectStart() {
		return subjectStart;
	}

	/**
	 * Returns the subject's end position of this hsp.
	 * @return
	 */
	public int getSubjectEnd() {
		return subjectEnd;
	}

	/**
	 * Returns the query's start position of this hsp.
	 * @return
	 */
	public int getQueryStart() {
		return queryStart;
	}

	/**
	 * Returns the query's end position of this hsp.
	 * @return
	 */
	public int getQueryEnd() {
		return queryEnd;
	}

	public void setAliLength(int aliLength) {
		this.aliLength = aliLength;
	}

	public void setQueryStart(int queryStart) {
		this.queryStart = queryStart;
	}

	public void setQueryEnd(int queryEnd) {
		this.queryEnd = queryEnd;
	}

	public void setSubjectStart(int subjectStart) {
		this.subjectStart = subjectStart;
	}

	public void setSubjectEnd(int subjectEnd) {
		this.subjectEnd = subjectEnd;
	}

	public void setEValue(double value) {
		eValue = value;
	}

	/**
	 * Returns the percent identity value of this hsp.
	 * @return
	 */
	public double getPercentIdentity() {
		return this.percentIdentity;
	}
	
	public void setPercentIdentity(double percentIdentity) {
		this.percentIdentity = percentIdentity;
	}
	
	public int getIdentities() {
		return (int)((percentIdentity/100.0)*aliLength);
	}
	
	public MultipleSequenceAlignment getAlignment() {
		return al;
	}
	
	/**
	 * Sets the alignment of this hsp given the 2 aligned sequences of query and subject
	 * @param querySeq
	 * @param subjectSeq
	 */
	public void setAlignment(String querySeq, String subjectSeq) {
		String[] tags = {parent.getQueryId(), parent.getSubjectId()};
		String[] seqs = {querySeq, subjectSeq};
		try {
			this.al = new MultipleSequenceAlignment(tags, seqs);
		} catch (AlignmentConstructionException e) {
			System.err.println("Error while constructing alignment from parsed blast output: "+e.getMessage());
		}
	}

	/**
	 * Returns an alignment result of transforming the hsp local alignment of this BlastHit 
	 * into one that contains the full sequences of query and subject given (padded 
	 * with gaps on the opposite sides). 
	 * e.g. If blast alignment is:
	 *  q:   ABC--DE--  (full q: bbABCDEeee) 
	 *  s:   -ABCD--EF  (full s: bbbABCDEFee)
	 * the new alignment will be:
	 *  q:  ---bbABC--DE--eee--
	 *  s:  bbb---ABCD--EF---ee
	 *      ^^^^^         ^^^^^
	 * The alignment of this BlastHit is unaffected. The returned alignment is a new object.
	 * @param fullQuerySeq
	 * @param fullSubjectSeq
	 * @return
	 */
	public MultipleSequenceAlignment getAlignmentFullSeqs(String fullQuerySeq, String fullSubjectSeq) {
		
		String querySeqNoGaps = this.al.getSequenceNoGaps(parent.getQueryId());
		String subjectSeqNoGaps = this.al.getSequenceNoGaps(parent.getSubjectId());

		if (!fullQuerySeq.contains(querySeqNoGaps)){
			throw new IllegalArgumentException("Given full query sequence is not a superstring of this BlastHit's alignment query sequence");
		}
		if (!fullSubjectSeq.contains(subjectSeqNoGaps)){
			throw new IllegalArgumentException("Given full subject sequence is not a superstring of this BlastHit's alignment subject sequence");
		}
		
		if (fullQuerySeq.length()==querySeqNoGaps.length() && fullSubjectSeq.length()==subjectSeqNoGaps.length()) {
			// the condition is equivalent to following, as a sanity check we also try it
			if (parent.getQueryLength()==fullQuerySeq.length() && parent.getSubjectLength()==fullSubjectSeq.length()) {
				// nothing to do, blast alignment is already spanning both full sequences of query and subject
				return this.al;				
			} else {
				System.err.println("Unexpected error: inconsistency between queryStart/End, subjectStart/End and sequences in stored alignment. Please report bug!");
				System.exit(1);
			}
		}
		
		String querySeq = this.al.getAlignedSequence(parent.getQueryId());
		String subjectSeq = this.al.getAlignedSequence(parent.getSubjectId());
		
		String newQuerySeq = getNGaps(this.subjectStart)+
							fullQuerySeq.substring(0, this.queryStart-1)+
							querySeq+
							fullQuerySeq.substring(this.queryEnd)+
							getNGaps(parent.getSubjectLength()-this.subjectEnd);
		String newSubjectSeq = fullSubjectSeq.substring(0, this.subjectStart-1)+
							getNGaps(this.queryStart)+
							subjectSeq+
							getNGaps(parent.getQueryLength()-this.queryEnd)+
							fullSubjectSeq.substring(this.subjectEnd);
		
		String[] tags = {parent.getQueryId(), parent.getSubjectId()};
		String[] seqs = {newQuerySeq, newSubjectSeq};
		MultipleSequenceAlignment newAln = null;
		try {
			newAln = new MultipleSequenceAlignment(tags, seqs);
		} catch (AlignmentConstructionException e) {
			System.err.println("Unexpected error: new alignment with full sequences from blast alignment couldn't be created. Please report the bug! Error: "+e.getMessage());
			System.exit(1);
		}
		return newAln;
	}
	
	/**
	 * Produces a string of n gap characters
	 * @param n
	 * @return
	 */
	private String getNGaps(int n) {
		StringBuffer buf = new StringBuffer();
		for (int i=0; i<n; i++)
		      buf.append (MultipleSequenceAlignment.GAPCHARACTER);
		return buf.toString();	
	}
	
	/**
	 * Return the hsp alignment of this BlastHit with the subjectId tag replaced by templateId 
	 * (pdbCode+pdbChaincode). If subjectId doesn't match regex {@value #ID_REGEX} then 
	 * alignment with normal tags is returned. 
	 * is returned.
	 * @see {@link #getTemplateId()}
	 * @param fullQuerySeq
	 * @param fullSubjectSeq
	 * @return
	 */
	public MultipleSequenceAlignment getAlignmentFullSeqsWithPDBTag(String fullQuerySeq, String fullSubjectSeq) {
		if (parent.getTemplateId()!=null) {
			MultipleSequenceAlignment aln = this.getAlignmentFullSeqs(fullQuerySeq, fullSubjectSeq);
			aln.resetTag(parent.getSubjectId(), parent.getTemplateId());
			return aln;
		} else {
			return getAlignmentFullSeqs(fullQuerySeq, fullSubjectSeq);
		}
	}

	/**
	 * Returns the query coverage for this hsp, i.e. 
	 * the ratio of aligned residues of the query compared to its length
	 * @return
	 */
	public double getQueryCoverage() {
		return  ((double)(getQueryEnd()-getQueryStart())/parent.getQueryLength());
	}
	
	/**
	 * Prints this BlastHsp in tabular format replicating blast's own 
	 * tabular output (blast's command line option -m 8). 
	 * The only field that we can't reproduce (because we don't parse it) is 
	 * the gap openings (column 6) where we print always a 0 instead.
	 */
	public void printTabular() {
		String scoreFormat = "%4.0f";
		if (score<100) {
			scoreFormat = "%4.1f";
		}
		System.out.printf("%s\t%s\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.0e\t"+scoreFormat+"\n",
				parent.getQueryId(),parent.getSubjectId(),
				percentIdentity,aliLength,aliLength-getIdentities(),0,queryStart,queryEnd,subjectStart,subjectEnd,eValue,score);
	}
	
}
