package owl.core.runners;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.sequence.alignment.AlignmentConstructionException;
import owl.core.sequence.alignment.MultipleSequenceAlignment;
import owl.core.structure.Pdb;
import owl.core.structure.PdbLoadError;
import owl.core.util.FileFormatError;

/**
 * A class for performing DALI Structural alignments via a locally installed DALI executable
 * @author Matthias Winkelmann
 */

public class DaliRunner {	
	
	private static final String DALI_REGEX = "(Query|Sbjct)\\s+([a-zA-Z\\.]+)";
	private Pdb first;
	private Pdb second;
	private File workdir;

	/**
	 * Using the provided query and subject Pdb files, DALI is executed in a temporary directory and 
	 * the resulting alignment converted to CLUSTAL format.  
	 * @param query
	 * @param subj
	 * @param dali_executable
	 * @param tempdir
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws AlignmentConstructionException
	 */
	
	public DaliRunner(Pdb query, Pdb subj,String dali_executable, String tempdir) throws IOException, InterruptedException, AlignmentConstructionException {
		
		first = query;
		second = subj;
		workdir = createTempDirectory(tempdir);
		first.writeToPDBFile(workdir.getAbsolutePath()+"/mod1.pdb");
		second.writeToPDBFile(workdir.getAbsolutePath()+"/mod2.pdb");
		
		String daliOutputFilename = workdir.getAbsolutePath()+"/"+first.getPdbCode()+first.getPdbChainCode()+
									"-"+second.getPdbCode()+second.getPdbChainCode()+".html";
		try {
		Process dali = Runtime.getRuntime().exec(new String[]{dali_executable,
				"-pairwise",
				workdir.getAbsolutePath()+"/mod1.pdb",
				workdir.getAbsolutePath()+"/mod2.pdb"}, 
				new String[]{"",""}, workdir);
		dali.waitFor();
		} catch (IOException e) {
			throw new AlignmentConstructionException("Could not run DALI. Check if "+dali_executable+" exists");
		}
		new File(workdir.getAbsolutePath()+"/aln.html").renameTo(new File(daliOutputFilename));
		correctFileDALIFormat(daliOutputFilename);
	}
	
	public String getClustalFile() {
		return workdir.getAbsolutePath()+"/alignment.clustal";
	}
	
	/**
	 * DALI output is html, which is converted to a CLUSTAL-formatted alignment. PDB and Chain identifier 
	 * are also restored. For PDB files with unobserved residues, these are reinserted into the alignment
	 * with corresponding gaps in the other sequence.  
	 * @param fileName a DALI html-formatted output file containing the CLUSTAL-like structure alignment
	 * @throws IOException
	 * @throws FileFormatError
	 */
	
  private void correctFileDALIFormat(String fileName) throws IOException, AlignmentConstructionException {

	String nextLine = "";
	String subj = "";
	String query = "";
	// open file

	BufferedReader fileIn = new BufferedReader(new FileReader(fileName));

	// We parse the provided filename to recover the Pdb and 
	// chain identifiers which are unfortunately not preserved by DALI
	
	Pattern p = Pattern.compile("(\\d[a-z0-9]{3}[A-Z])\\-(\\d[a-z0-9]{3}[A-Z])");
	Matcher m = p.matcher(fileName);
	m.find();
	String firstPdbID = m.group(1);
	String secondPdbID = m.group(2);
	
	p = Pattern.compile(DALI_REGEX);

	// read sequences
	try {
	while ((nextLine = fileIn.readLine()) != null) {
		if (nextLine.startsWith("Query")) {
			m = p.matcher(nextLine);
			m.find();
			query += m.group(2);
		} else if (nextLine.startsWith("Sbjct")) {
			m = p.matcher(nextLine);
			m.find();
			subj += m.group(2);
		}
	}
	} catch (IllegalStateException e) {
		throw new AlignmentConstructionException("Could not read DALI alignment. Check "+fileName+" for errors");
	}

	// We convert the dot used by DALI to whatever we are using
	
	query = query.toUpperCase().replace('.',MultipleSequenceAlignment.GAPCHARACTER);
	subj = subj.toUpperCase().replace('.',MultipleSequenceAlignment.GAPCHARACTER);
	fileIn.close();
	
	
	// DALI removes unobserved residues from the alignment, so we have to reinsert those
	
	Entry<Integer, Character> missing; // the next missing residue that we have to insert
	int entryKey;
	int entryPos;
	if (first.countUnobserved() > 0) {
	
		TreeMap<Integer, Character> unobserved1 = first.getUnobservedResidues();
		entryKey = 0;
		for (int i = 0; i < unobserved1.size();i++) {
			missing = unobserved1.higherEntry(entryKey);
			entryKey = missing.getKey();
			entryPos = indexWithoutGaps2IndexWithGaps(entryKey-1,query);
			query = query.substring(0,entryPos)+missing.getValue()+query.substring(entryPos);
			subj = subj.substring(0,entryPos)+MultipleSequenceAlignment.GAPCHARACTER+subj.substring(entryPos);
		}
		
	}
	
	if (second.countUnobserved() > 0) {
		entryKey = 0;
		TreeMap<Integer, Character> unobserved2 = second.getUnobservedResidues();
		System.out.println(second.countUnobserved());
		for (int i = 0; i < unobserved2.size();i++) {
			missing = unobserved2.higherEntry(entryKey);
			entryKey = missing.getKey();
			entryPos = indexWithoutGaps2IndexWithGaps(entryKey-1,subj);
			query = query.substring(0,entryPos)+MultipleSequenceAlignment.GAPCHARACTER+query.substring(entryPos);
			subj = subj.substring(0,entryPos)+missing.getValue()+subj.substring(entryPos);
		}
		
	}
	 
	
	// write sequences in clustal format
	
	Writer fileOut = new FileWriter(workdir.getAbsolutePath()+"/"+"alignment.clustal");
	fileOut.write("CLUSTAL\n");
	fileOut.append(firstPdbID+"   "+query+"\n");
	fileOut.append(secondPdbID+"   "+subj+"\n");
	fileOut.close();

  }
  
  
  	private int indexWithoutGaps2IndexWithGaps(int index,String s) {
  		
  		int nongap = 0;
  		int i = 0;
  		for (i = 0;i<s.length() && nongap < index;i++) {
  			if (s.charAt(i) != MultipleSequenceAlignment.GAPCHARACTER) {
  				nongap++;
  			}
  		}
  		return i;
  			
  	}
	
	private static File createTempDirectory(String tempdir)
    throws IOException
{
    final File temp;

    temp = new File(tempdir+"/cmviewDALI"+Long.toString(System.nanoTime()));

   
    if(!(temp.mkdir()))
    {
        throw new IOException("Could not create temp directory: " + temp.getAbsolutePath());
    }
    temp.deleteOnExit();
    return (temp);
}

	/**
	 * Main (Test) method
	 * @throws IOException 
	 * @throws PdbLoadError 
	 * @throws InterruptedException 
	 * @throws FileFormatError 
	 */
	
	public static void main(String[] args) throws IOException, PdbLoadError, InterruptedException {
		

	}

}
