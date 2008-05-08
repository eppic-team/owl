package sequence;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;

import tools.MySQLConnection;

/**
 * A GTG (global trace graph) output parser
 * See http://ekhidna.biocenter.helsinki.fi/gtg/start and
 * http://ekhidna.biocenter.helsinki.fi/casp8/
 *
 */
public class GTGParser {
	
	private File gtgOutputFile;
	private GTGHitList hits;

	/**
	 * Constructs a GTGParser and parses the given file.
	 * To get the results call {@link #getHits()}
	 * @param gtgOutputFile
	 * @throws IOException
	 */
	public GTGParser(File gtgOutputFile) throws IOException {
		this.gtgOutputFile = gtgOutputFile;
		this.parse();
	}
	
	private void parse() throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(gtgOutputFile));
		String line;
		int hitCounter = 0;
		int lastNon0QuerySerial = 0, lastNon0SubjectSerial = 0;
		hits = new GTGHitList();
		GTGHit hit = null;
		String subjectSequence = "";
		int identities = 0;
		int alnLength = 0;
		boolean onAlnRegion = false; // flag to indicate whether we are within the aligned region or out of it
		while ((line=br.readLine())!=null) {
			if (line.equals("")) continue;
			if (line.startsWith("#")) {
				if (hit!=null) { //i.e. we are not in first hit
					// we assign the last seen non-zero query/subject serial as the end query/subject serials for the last hit
					hit.setSubjectEnd(lastNon0SubjectSerial);
					hit.setQueryEnd(lastNon0QuerySerial);
					hit.setSubjectSequence(subjectSequence);
					hit.setIdentities(identities);
					hit.setAlnLength(alnLength);
					lastNon0SubjectSerial = 0;
					lastNon0QuerySerial = 0;
					subjectSequence = "";
					onAlnRegion = false;
					alnLength = 0;
					identities = 0;
				}
				hitCounter++;
				Pattern p = Pattern.compile("^#\\s+target=(\\w+)\\s+score1=(\\d+)\\s+score2=(\\d+)\\s+score3=(\\d+)\\s+snid=(\\w+)\\s+pdbid=(\\w+)");
				Matcher m = p.matcher(line);
				if (m.matches()) {
					String queryId = m.group(1);
					int totalScore = Integer.parseInt(m.group(2));
					int consistencyScore = Integer.parseInt(m.group(3));
					int motifScore = Integer.parseInt(m.group(4));
					// String snid = m.group(5);// we don't need this for now, what is it anyway?
					String subjectId = m.group(6);
					hit = new GTGHit(queryId, subjectId, totalScore, consistencyScore, motifScore);
					hits.add(hit);
				}
			} else {
				String[] fields = line.split("\\t");
				int querySerial = Integer.parseInt(fields[2]);
				int subjectSerial = Integer.parseInt(fields[3]);

				if (querySerial!=0 && subjectSerial!=0) {
					onAlnRegion = true;
					// if this hit doesn't have a start assigned yet we assign it, i.e. firt non-zero query/subject serial is assigned as the start  
					if (hit.queryStart==0 && hit.subjectStart==0) {
						hit.setQueryStart(querySerial);
						hit.setSubjectStart(subjectSerial);
					}
					
					// we keep the last found query/subject serials so that we can assign the end query/subject positions the next time we find a "#" line
					lastNon0QuerySerial = querySerial;
					lastNon0SubjectSerial = subjectSerial;
				}
				if (onAlnRegion) alnLength++;
				
				String alnField = fields[4];
				String[] tokens = alnField.split("\\s");
				if (!tokens[2].equals("-"))
					subjectSequence += tokens[2];
				if (tokens[1].equals("=="))
					identities++;
			}
		}
		
		// after we finish reading the last line we still miss to add the alignment data for the last hit, we do it now:
		hit.setSubjectEnd(lastNon0SubjectSerial);
		hit.setQueryEnd(lastNon0QuerySerial);
		hit.setSubjectSequence(subjectSequence);
		hit.setIdentities(identities);
		hit.setAlnLength(alnLength);
			
		br.close();
	}
	
	/**
	 * Gets the GTGHitList result of parsing the file given in constructor
	 * @return
	 */
	public GTGHitList getHits() {
		return this.hits;
	}
	
	/**
	 * Main method to test the class
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		File file = new File("/project/StruPPi/CASP8/dryrun/casp7_gtg/T0345.gtg");
		file = new File("/project/StruPPi/CASP8/gtg/T0387.gtg");
		GTGParser gtgparser = new GTGParser(file);
		GTGHitList hits = gtgparser.getHits();
		
		System.out.println("Number of hits: "+hits.size());
		hits.printTabular();
		
		MySQLConnection conn = new MySQLConnection("white","duarte","nieve");
		String pdbaseDb ="pdbase";


		for (GTGHit hit:hits) {
			System.out.println(hit.subjectId);
			try {
				hit.checkSubjectSequence(conn, pdbaseDb);
				hit.reassignSubjectSerials(conn, pdbaseDb);
			} catch (PdbLoadError e) {
				System.err.println(e.getMessage());
			} catch (PdbCodeNotFoundError e) {
				System.err.println(e.getMessage());
			}
		}
		
		conn.close();

		System.out.println("After reassignment of subject serials");
		hits.print();
	}
}
