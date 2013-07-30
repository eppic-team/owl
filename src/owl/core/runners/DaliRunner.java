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
import owl.core.structure.PdbChain;
import owl.core.util.FileFormatException;

/**
 * A class for performing DALI Structural alignments via a locally installed DALI executable
 * @author Matthias Winkelmann
 */
public class DaliRunner {	

	private static final String DALI_REGEX = "(Query|Sbjct)\\s+([a-zA-Z\\.\\-]+)";
	private PdbChain first;
	private PdbChain second;
	private File workdir;

	/**
	 * Using the provided query and subject Pdb files, DALI is executed in a temporary directory and 
	 * the resulting alignment converted to CLUSTAL format.  
	 * @param query
	 * @param subj
	 * @param queryTag
	 * @param subjTag
	 * @param dali_executable
	 * @param tempdir
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws AlignmentConstructionException
	 */
	public DaliRunner(PdbChain query, PdbChain subj, String queryTag, String subjTag, String dali_executable, String tempdir) throws IOException, InterruptedException, AlignmentConstructionException {

		first = query;
		second = subj;
		workdir = createTempDirectory(tempdir);
		first.writeToPDBFile(new File(workdir,"mod1.pdb"));
		second.writeToPDBFile(new File(workdir,"mod2.pdb"));

		Process dali = Runtime.getRuntime().exec(new String[]{dali_executable,
				"-pairwise",
				workdir.getAbsolutePath()+"/mod1.pdb",
				workdir.getAbsolutePath()+"/mod2.pdb"}, 
				new String[]{"",""}, workdir);
		dali.waitFor();

		// output can be aln.html or index.html depending on dali version (I think!)
		File daliActualOutputFile = new File(workdir.getAbsolutePath(),"aln.html");
		if (!daliActualOutputFile.exists()) {
			daliActualOutputFile = new File(workdir.getAbsolutePath(),"index.html");
			if (!daliActualOutputFile.exists()) {
				throw new IOException("Could not find DALI output file aln.html or index.html in DALI working dir "+workdir);
			}
		}

		correctFileDALIFormat(daliActualOutputFile, queryTag, subjTag);
	}

	public String getClustalFile() {
		return workdir.getAbsolutePath()+"/alignment.clustal";
	}

	/**
	 * DALI output is html, which is converted to a CLUSTAL-formatted alignment. 
	 * For PDB files with unobserved residues, these are reinserted into the alignment
	 * with corresponding gaps in the other sequence.  
	 * @param file a DALI html-formatted output file containing the CLUSTAL-like structure alignment
	 * @param queryTag
	 * @param subjTag
	 * @throws IOException
	 * @throws FileFormatException
	 */
	private void correctFileDALIFormat(File file, String queryTag, String subjTag) throws IOException, AlignmentConstructionException {

		String nextLine = "";
		String subj = "";
		String query = "";
		// open file

		BufferedReader fileIn = new BufferedReader(new FileReader(file));


		// with the regex we capture any letter, . (dot) or - (hyphen). Apparently some version of DALI use '.', some '-' for the gap character
		Pattern p = Pattern.compile(DALI_REGEX);

		// read sequences
		try {
			while ((nextLine = fileIn.readLine()) != null) {
				if (nextLine.startsWith("Query")) {
					Matcher m = p.matcher(nextLine);
					if (m.find()) {
						query += m.group(2);
					}
				} else if (nextLine.startsWith("Sbjct")) {
					Matcher m = p.matcher(nextLine);
					if (m.find()) {
						subj += m.group(2);
					}
				}
			}
		} catch (IllegalStateException e) {
			fileIn.close();
			throw new AlignmentConstructionException("Could not read DALI alignment. Check "+file+" for errors");
		}

		// We convert the dot used by DALI to whatever we are using (note newer versions of DALI use hyphens directly)

		query = query.toUpperCase().replace('.',MultipleSequenceAlignment.GAPCHARACTER);
		subj = subj.toUpperCase().replace('.',MultipleSequenceAlignment.GAPCHARACTER);
		fileIn.close();


		// DALI removes unobserved residues from the alignment, so we have to reinsert those

		Entry<Integer, Character> missing; // the next missing residue that we have to insert
		int entryKey;
		int entryPos;
		if ((first.getFullLength()-first.getStdAaObsLength()) > 0) {

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

		if ((second.getFullLength()-second.getStdAaObsLength()) > 0) {
			entryKey = 0;
			TreeMap<Integer, Character> unobserved2 = second.getUnobservedResidues();
			System.out.println((second.getFullLength()-second.getStdAaObsLength()));
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
		fileOut.append(queryTag+"   "+query+"\n");
		fileOut.append(subjTag+"   "+subj+"\n");
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

	private static File createTempDirectory(String tempdir) throws IOException {
		
		final File temp;

		temp = new File(tempdir+"/cmviewDALI"+Long.toString(System.nanoTime()));


		if(!(temp.mkdir()))
		{
			throw new IOException("Could not create temp directory: " + temp.getAbsolutePath());
		}
		temp.deleteOnExit();
		return (temp);
	}


}
