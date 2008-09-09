package sequence;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A blast output record  
 * 
 */
public class BlastHit {
	
	public static final int OUTPUT_LENGTH = 80;
	private static final String ID_REGEX = "pdb\\|(\\d\\w\\w\\w)\\|(\\w)";
	
	private String queryId;
	private String subjectId;
	private double percentIdentity;
	private int aliLength;
	private int missmatches;
	private int gapOpenings;
	private int queryStart;
	private int queryEnd;
	private int subjectStart;
	private int subjectEnd;
	private double eValue;
	private double score;

	public BlastHit(String queryId, String subjectId, double percentIdentity, int aliLength, int missmatches, int gapOpenings, 
			int queryStart, int queryEnd, int subjectStart, int subjectEnd, double eValue, double score) {
		this.queryId = queryId;
		this.subjectId = subjectId;
		this.percentIdentity = percentIdentity;
		this.aliLength = aliLength;
		this.missmatches = missmatches;
		this.gapOpenings = gapOpenings;
		this.queryStart = queryStart;
		this.queryEnd = queryEnd;
		this.subjectStart = subjectStart;
		this.subjectEnd = subjectEnd;
		this.eValue = eValue;
		this.score = score;
		
	}
	
	public BlastHit(String field0, String field1, String field2, String field3, String field4, String field5, 
			String field6, String field7, String field8, String field9, String field10, String field11) {
		this.queryId = field0.trim();
		this.subjectId = field1.trim();
		this.percentIdentity = Double.parseDouble(field2.trim());
		this.aliLength = Integer.parseInt(field3.trim());
		this.missmatches = Integer.parseInt(field4.trim());
		this.gapOpenings = Integer.parseInt(field5.trim());
		this.queryStart = Integer.parseInt(field6.trim());
		this.queryEnd = Integer.parseInt(field7.trim());
		this.subjectStart = Integer.parseInt(field8.trim());
		this.subjectEnd = Integer.parseInt(field9.trim());
		this.eValue = Double.parseDouble(field10.trim());
		this.score = Double.parseDouble(field11.trim());
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
		System.out.printf("%5s\t%10s\t%5.1f\t%8.1e\t%4.0f ", queryId, subjectId, percentIdentity, eValue, score);
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
	public static void printHeaderWithOverview(int queryLength, double scaleFactor) {
		System.out.printf("%5s\t%10s\t%5s\t%8s\t%4s ", "query", "subject", "id%", "e-val", "sc");
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
	 * Returns the number of missmatches for this hit.
	 * @return
	 */
	public int getMissmatches() {
		return missmatches;
	}

	/**
	 * Returns the number of gap openings for this hit.
	 * @return
	 */
	public int getGapOpenings() {
		return gapOpenings;
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
}
