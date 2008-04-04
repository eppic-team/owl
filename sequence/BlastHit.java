package sequence;

/**
 * A blast output record  
 * 
 */
public class BlastHit {
	
	String queryId;
	String subjectId;
	double percentIdentity;
	int aliLength;
	int missmatches;
	int gapOpenings;
	int queryStart;
	int queryEnd;
	int subjectStart;
	int subjectEnd;
	double eValue;
	double score;

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
}
