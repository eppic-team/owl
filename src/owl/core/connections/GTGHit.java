package owl.core.connections;

import java.sql.SQLException;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbCodeNotFoundException;
import owl.core.structure.PdbLoadException;
import owl.core.structure.PdbaseParser;
import owl.core.util.MySQLConnection;


/**
 * A GTG output hit
 *
 */
public class GTGHit {
	
	private static final String NULL_PDBCHAINCODE = "A";
	
	private static final String ID_REGEX = "(\\d\\w\\w\\w)(\\w)";

	String queryId;
	String subjectId;

	int queryStart;
	int queryEnd;
	int subjectStart;
	int subjectEnd;
	
	int identities;
	//int missmatches; // do we want this?
	int alnLength;
	// alignment // do we want to store the full alignment or is it enough with the above?

	int totalScore;
	int consistencyScore;
	int motifScore;

	String subjectSequence;
	

	
	public GTGHit(String queryId, String subjectId, int totalScore, int consistencyScore, int motifScore) {
		this.queryId = queryId;
		this.subjectId = subjectId;
		
		// We asume here that the subjectId is a PDB id in the form pdbCode+pdbChainCode
		// It seems that GTG uses pre-remediated PDB data: we have to add an "A" when length of subjectId is only 4, i.e. missing (NULL) chain code 
		if (this.subjectId.length()==4) {
			this.subjectId += NULL_PDBCHAINCODE; 
		}
		
		this.totalScore = totalScore;
		this.consistencyScore = consistencyScore;
		this.motifScore = motifScore;
		
	}

	public int getTotalScore() {
		return this.totalScore;
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
	
	public void setSubjectSequence(String subjectSequence) {
		this.subjectSequence = subjectSequence;
	}

	public String getSubjectSequence() {
		return subjectSequence;
	}
	
	public void setIdentities(int identities) {
		this.identities = identities;
	}
	
	public void setAlnLength(int alnLength) {
		this.alnLength = alnLength;
	}
	
	public double getPercentIdentity() {
		return ((double) identities/alnLength)*100;
	}
	
	public void print() {
		System.out.printf("%8s %8s %6d %6d %6d %6d %6d %6d %5.1f %6d %6d %6d \n", 
				queryId, subjectId, queryStart, queryEnd, subjectStart, subjectEnd, identities, alnLength, getPercentIdentity(),totalScore, consistencyScore, motifScore);
	}
	
	/**
	 * Checks the subject sequence and compares it to the pdbase sequence, 
	 * provided that the subjectId is a PDB identifier in the form pdbCode+pdbChainCode
	 * @param conn
	 * @param pdbaseDb
	 * @return false if sequences don't match, true if they do
	 * @throws PdbCodeNotFoundException
	 * @throws SQLException
	 * @throws PdbLoadException
	 */
	public boolean checkSubjectSequence(MySQLConnection conn, String pdbaseDb) throws PdbCodeNotFoundException, SQLException, PdbLoadException {
		if (subjectSequence.equals("")) {
			System.err.println("Subject sequence is empty");
			return false;
		}
		String pdbCode = subjectId.substring(0, 4);
		String pdbChainCode = subjectId.substring(4);
		PdbAsymUnit fullpdb = new PdbAsymUnit(pdbCode,conn,pdbaseDb);
		PdbChain pdb = fullpdb.getChain(pdbChainCode);
		TreeMap<Integer, Integer> mapObsSerials2resser = PdbaseParser.getObservedResMapping(pdbCode,pdbChainCode,conn,pdbaseDb);
		String pdbaseSequence = pdb.getSequence().getSeq();
		
		for (int i=0; i<subjectSequence.length();i++) {
			int indexInFullSeq = mapObsSerials2resser.get(i+1)-1;
			if (subjectSequence.charAt(i)!=pdbaseSequence.charAt(indexInFullSeq)) {
				System.err.println("GTG output subject's sequence doesn't match sequence from pdbase at position "+(i+1)+"(internal res. serial: "+(indexInFullSeq+1)+"). subject: "+subjectSequence.charAt(i)+", pdbase: "+pdbaseSequence.charAt(indexInFullSeq));
				return false;
			}
		}

// 		the following would have worked as well. But our getObservedSequence() was getting observed sequences without non-standard aas, so the matching of sequences was failing a lot 
//		String pdbaseSequence = pdb.getObservedSequence();
//		
//		// 1st check: subject sequence must not be longer than full sequence from pdbase
//		// this also prevents for getting an out of bound index error in following loop
//		if (subjectSequence.length()>pdbaseSequence.length()) {
//			System.err.println("GTG output subject's sequence is longer than sequence from pdbase.");
//			return false;
//		}
//
//		for (int i=0; i<subjectSequence.length();i++) {					
//			// 2nd check: all residues match
//			if (subjectSequence.charAt(i)!=pdbaseSequence.charAt(i)) {
//				System.err.println("GTG output subject's sequence doesn't match sequence from pdbase at position "+(i+1)+". subject: "+subjectSequence.charAt(i)+", pdbase: "+pdbaseSequence.charAt(i));
//				return false;
//			}
//		}
		return true;
	}
	
	/**
	 * Reassigns the serials of subjectStart and subjectEnd positions based on the 
	 * pdbase full SEQRES sequence, i.e. converts GTG serials to our internal residue
	 * serials (cif serials)
	 * NOTE: there are some data compatibility issues with the output of GTG. Sometimes 
	 * our mapping of observed sequence serials (numbering used by GTG) and pdbase full 
	 * sequence is not perfect, so it might be off a few residues at most. Anyway this 
	 * happens in extremely rare cases where a residue is considered observed in pdbase 
	 * and not by GTG (examples are: 1fdpB residue 136, 1tonA residue 79)
	 * @param conn
	 * @param pdbaseDb
	 * @throws SQLException
	 * @throws PdbCodeNotFoundException
	 * @throws PdbLoadException
	 */
	public void reassignSubjectSerials(MySQLConnection conn, String pdbaseDb) throws SQLException, PdbCodeNotFoundException, PdbLoadException {
		if (subjectStart==0 && subjectEnd ==0 ) return; // subject start and end not known: nothing to do
		
		String pdbCode = subjectId.substring(0, 4);
		String pdbChainCode = subjectId.substring(4);
		TreeMap<Integer, Integer> mapObsSerials2resser = PdbaseParser.getObservedResMapping(pdbCode,pdbChainCode,conn,pdbaseDb);
		
		if (!mapObsSerials2resser.containsKey(subjectStart) || !mapObsSerials2resser.containsKey(subjectEnd)) {
			System.err.println("Cannot map observed residue sequence serials to internal residue serials for "+subjectId+": either subject start or end serial are out of range");
			return; // can't map: return 
		}
		subjectStart = mapObsSerials2resser.get(subjectStart);
		subjectEnd = mapObsSerials2resser.get(subjectEnd);
	}
	
	/**
	 * Prints a few selected fields for this GTG hit plus a graphical representation of the match.
	 * The match is scaled by the given scale factor and rounded to screen columns. 
	 * @param scaleFactor
	 */
	public void printWithOverview(double scaleFactor) {
		System.out.printf("%5s\t%5s\t%5.1f\t%7d ", queryId, subjectId, getPercentIdentity(), totalScore);
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
		System.out.printf("%5s\t%5s\t%5s\t%7s ", "query", "sbjct", "id%", "sc");
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
	 * Returns the template id as concatenated pdbCode+pdbChainCode e.g. 1abcA
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
