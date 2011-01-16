package owl.core.runners.blast;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import owl.core.util.FileFormatError;


/**
 * A simple parser for blast tabular output files 
 *
 */
public class BlastTabularParser {

	private File blastOutFile;
	private BlastHitList hits;
	
	/**
	 * Constructs the parser object for given file and performs parsing
	 * To get the parsed values use {@link #getHits()}
	 * @param blastOutFile
	 * @throws IOException if I/O error while reading the file
	 * @throws FileFormatError if given file doesn't follow blast's tabular 
	 * output format
	 */
	public BlastTabularParser(File blastOutFile) throws IOException, FileFormatError {
		this.blastOutFile = blastOutFile;
		this.parse();
	}
	
	private void parse() throws IOException, FileFormatError {
		hits = new BlastHitList();
		String lastSubjectId = null;
		BlastHit currentHit = null;
		BufferedReader br = new BufferedReader(new FileReader(blastOutFile));
		String line;
		int lineCount = 0;
		while ((line=br.readLine())!=null) {
			lineCount++;
			if (line.startsWith("#")) continue;
			
			String[] fields = line.split("\t");
			if (fields.length!=12) 
				throw new FileFormatError("Blast output file "+blastOutFile+" doesn't seem to have the right format at line "+lineCount);
			
			BlastHsp hsp = new BlastHsp(fields[0],fields[1],fields[2],fields[3],fields[4],fields[5],fields[6],fields[7],fields[8],fields[9],fields[10],fields[11]);
			String subjectId = fields[1].trim();
			if (lastSubjectId!=null && lastSubjectId.equals(subjectId)) {
				currentHit.addHsp(hsp);				
			} else {
				currentHit = new BlastHit();
				currentHit.addHsp(hsp);
				currentHit.setQueryId(fields[0].trim());
				currentHit.setQueryLength(-1);
				currentHit.setSubjectId(subjectId);
				currentHit.setSubjectDef(null);
				currentHit.setSubjectLength(-1);
				hits.add(currentHit);
			}
			lastSubjectId = subjectId;
		}
		br.close();
		
	}
	
	/**
	 * Gets the list of blast hits parsed from the blast tabular output file
	 * @return
	 */
	public BlastHitList getHits() {
		return this.hits;
	}
}
