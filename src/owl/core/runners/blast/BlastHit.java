package owl.core.runners.blast;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.sequence.alignment.AlignmentConstructionError;
import owl.core.sequence.alignment.MultipleSequenceAlignment;


/**
 * A blast output record.
 * If a subject matches multiple times in different regions we consider that separate blast records, 
 * i.e. a BlastHit correspond to a <Hsp> record in Blast XML output.   
 * 
 * 
 */
public class BlastHit {
	
	public static final int OUTPUT_LENGTH = 80;
	private static final String ID_REGEX = "pdb\\|(\\d\\w\\w\\w)\\|(\\w)";
	
	private String queryId; // queryId and queryLength are redundant here (they belong in BlastHitList) but anyway useful to have copies here
	private int queryLength;
	private String subjectId;
	private String subjectDef;
	private double percentIdentity;
	private int aliLength;
	private int queryStart;
	private int queryEnd;
	private int subjectStart;
	private int subjectEnd;
	private int subjectLength;
	private double eValue;
	private double score;
	private MultipleSequenceAlignment al; 		// if BlastHit create from tabular file parsing, this will be null 

	public BlastHit(String field0, String field1, String field2, String field3, String field4, String field5, 
			String field6, String field7, String field8, String field9, String field10, String field11) {
		this.queryId = field0.trim();
		this.subjectId = field1.trim();
		this.percentIdentity = Double.parseDouble(field2.trim());
		this.aliLength = Integer.parseInt(field3.trim());
		//this.missmatches = Integer.parseInt(field4.trim());
		//this.gapOpenings = Integer.parseInt(field5.trim());
		this.queryStart = Integer.parseInt(field6.trim());
		this.queryEnd = Integer.parseInt(field7.trim());
		this.subjectStart = Integer.parseInt(field8.trim());
		this.subjectEnd = Integer.parseInt(field9.trim());
		this.eValue = Double.parseDouble(field10.trim());
		this.score = Double.parseDouble(field11.trim());
		this.al = null;
	}

	public BlastHit() {
		
	}
	
	/**
	 * Prints a few selected fields for this blast hit 
	 */
	public void print() {
		System.out.println(queryId+"\t"+subjectId+"\t"+percentIdentity+"\t"+eValue+"\t"+score);
	}
	
	/**
	 * Prints a few selected fields for this blast hit plus a graphical representation of the match.
	 * The match is scaled by the given scale factor and rounded to screen columns. 
	 * @param scaleFactor
	 */
	public void printWithOverview(double scaleFactor) {
		System.out.printf("%"+queryId.length()+"s\t%10s\t%5.1f\t%8.1e\t%4.0f ", queryId, subjectId, percentIdentity, eValue, score);
		int beg = (int) Math.floor(scaleFactor * this.queryStart);
		int end = (int) Math.ceil(scaleFactor * this.queryEnd);
		printOverviewLine(beg, end);
	}
	
	/**
	 * Print the column headers corresponding to the printWithOverview() method.
	 * Additionally prints a graphical overview of the query (queryLength scaled by scaleFactor).
	 * @param queryLength
	 * @param scaleFactor
	 */
	public static void printHeaderWithOverview(int queryLength, double scaleFactor, int queryIDlength) {
		System.out.printf("%"+queryIDlength+"s\t%10s\t%5s\t%8s\t%4s ", "query", "subject", "id%", "e-val", "sc");
		int beg = 1;
		int end = (int) Math.ceil(scaleFactor * queryLength);
		printOverviewLine(beg, end);		
	}
	
	/**
	 * Print one line of the match overview.
	 * @param beg the beginning of the match in screen columns
	 * @param end the end of the match in screen columns
	 */
	private static void printOverviewLine(int beg, int end) {
		for (int i = 1; i < beg; i++) {
			System.out.print(" ");
		}
		for (int i = beg; i <= end; i++) {
			System.out.print("-");
		}
		System.out.println();
	}
	
	/** 
	 * Returns the e-value of this hit.
	 * @return
	 */
	public double getEValue() {
		return this.eValue;
	}
	
	/**
	 * Returns the percent identity value of this hit.
	 * @return
	 */
	public double getPercentIdentity() {
		return this.percentIdentity;
	}
	
	/**
	 * Returns the alignment length of this hit.
	 * @return
	 */
	public int getAliLength() {
		return aliLength;
	}

	/**
	 * Returns the subject start position of this hit.
	 * @return
	 */
	public int getSubjectStart() {
		return subjectStart;
	}

	/**
	 * Returns the subject's end position of this hit.
	 * @return
	 */
	public int getSubjectEnd() {
		return subjectEnd;
	}

	/**
	 * Returns the query's start position of this hit.
	 * @return
	 */
	public int getQueryStart() {
		return queryStart;
	}

	/**
	 * Returns the query's end position of this hit.
	 * @return
	 */
	public int getQueryEnd() {
		return queryEnd;
	}

	public String getQueryId() {
		return queryId;
	}

	public void setQueryId(String queryId) {
		this.queryId = queryId;
	}

	public int getQueryLength() {
		return this.queryLength;
	}
	
	public void setQueryLength(int queryLength) {
		this.queryLength = queryLength;
	}
	
	public String getSubjectId() {
		return subjectId;
	}

	public void setSubjectId(String subjectId) {
		this.subjectId = subjectId;
	}
	
	public String getSubjectDef() {
		return subjectDef;
	}
	
	public void setSubjectDef(String subjectDef) {
		this.subjectDef = subjectDef;
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}

	public void setPercentIdentity(double percentIdentity) {
		this.percentIdentity = percentIdentity;
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

	public void setSubjectLength(int subjectLength) {
		this.subjectLength = subjectLength;
	}

	public int getSubjectLength() {
		return subjectLength;
	}
	
	/**
	 * Sets the alignment of this hit given the 2 aligned sequences of query and subject
	 * @param querySeq
	 * @param subjectSeq
	 */
	public void setAlignment(String querySeq, String subjectSeq) {
		String[] tags = {queryId, subjectId};
		String[] seqs = {querySeq, subjectSeq};
		try {
			this.al = new MultipleSequenceAlignment(tags, seqs);
		} catch (AlignmentConstructionError e) {
			System.err.println("Error while constructing alignment from parsed blast output: "+e.getMessage());
		}
	}

	/**
	 * Return the hsp alignment of this BlastHit with tags queryId and subjectId
	 * @return
	 */
	public MultipleSequenceAlignment getAlignment() {
		return this.al;	
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
		
		String querySeqNoGaps = this.al.getSequenceNoGaps(queryId);
		String subjectSeqNoGaps = this.al.getSequenceNoGaps(subjectId);

		if (!fullQuerySeq.contains(querySeqNoGaps)){
			throw new IllegalArgumentException("Given full query sequence is not a superstring of this BlastHit's alignment query sequence");
		}
		if (!fullSubjectSeq.contains(subjectSeqNoGaps)){
			throw new IllegalArgumentException("Given full subject sequence is not a superstring of this BlastHit's alignment subject sequence");
		}
		
		if (fullQuerySeq.length()==querySeqNoGaps.length() && fullSubjectSeq.length()==subjectSeqNoGaps.length()) {
			// the condition is equivalent to following, as a sanity check we also try it
			if (this.queryLength==fullQuerySeq.length() && this.subjectLength==fullSubjectSeq.length()) {
				// nothing to do, blast alignment is already spanning both full sequences of query and subject
				return this.al;				
			} else {
				System.err.println("Unexpected error: inconsistency between queryStart/End, subjectStart/End and sequences in stored alignment. Please report bug!");
				System.exit(1);
			}
		}
		
		String querySeq = this.al.getAlignedSequence(queryId);
		String subjectSeq = this.al.getAlignedSequence(subjectId);
		
		String newQuerySeq = getNGaps(this.subjectStart)+
							fullQuerySeq.substring(0, this.queryStart-1)+
							querySeq+
							fullQuerySeq.substring(this.queryEnd)+
							getNGaps(this.subjectLength-this.subjectEnd);
		String newSubjectSeq = fullSubjectSeq.substring(0, this.subjectStart-1)+
							getNGaps(this.queryStart)+
							subjectSeq+
							getNGaps(this.queryLength-this.queryEnd)+
							fullSubjectSeq.substring(this.subjectEnd);
		
		String[] tags = {this.queryId, this.subjectId};
		String[] seqs = {newQuerySeq, newSubjectSeq};
		MultipleSequenceAlignment newAln = null;
		try {
			newAln = new MultipleSequenceAlignment(tags, seqs);
		} catch (AlignmentConstructionError e) {
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
		if (this.getTemplateId()!=null) {
			MultipleSequenceAlignment aln = this.getAlignmentFullSeqs(fullQuerySeq, fullSubjectSeq);
			aln.resetTag(this.subjectId, this.getTemplateId());
			return aln;
		} else {
			return getAlignmentFullSeqs(fullQuerySeq, fullSubjectSeq);
		}
	}
	
	/**
	 * Returns the query coverage for this hit's alignment, i.e. 
	 * the ratio of aligned residues of the query compared to its length
	 * @return
	 */
	public double getQueryCoverage() {
		return  ((double)(this.getQueryEnd()-this.getQueryStart())/this.getQueryLength());
	}
	
	/**
	 * Returns the template id as concatenated pdb code + chain code e.g. 1abcA
	 * @return the template id or null if queryId is not in the right format
	 */
	public String getTemplateId() {
		Pattern p = Pattern.compile(ID_REGEX);
		Matcher m = p.matcher(subjectId);
		if (m.matches()) {
			return m.group(1).toLowerCase()+m.group(2).toUpperCase();
		}
		return null;
	}
	
	// statics
	
	/**
	 * From a template id in the form pdbCode+pdbChainCode (e.g. 1abcA) returns a subjectId
	 * as it is in the blastable pdb sequences file: pdb|1ABC|A
	 * This is exactly the reverse of {@link #getTemplateId()}
	 * See also {@link #ID_REGEX} 
	 * @param templateId
	 * @return 
	 */
	public static String templateIdToSubjectId(String templateId) {
		Pattern p = Pattern.compile("(\\d\\w\\w\\w)(\\w)");
		Matcher m = p.matcher(templateId);
		if (m.matches()) {
			return "pdb|"+m.group(1).toUpperCase()+"|"+m.group(2);
		} else {
			throw new IllegalArgumentException("Given templateId must be of the form pdbCode+pdbChainCode, e.g. 1abcA");
		}
	}
}
