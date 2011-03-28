package owl.core.runners;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.sequence.GeneticCodeType;
import owl.core.structure.AminoAcid;
import owl.core.util.FileFormatException;

/**
 * Class to run the selecton program for identification of site-specific positive selection
 * and purifying selection.
 * 
 * See http://selecton.tau.ac.il/overview.html
 * 
 * @author duarte_j
 *
 */
public class SelectonRunner {

	private File selectonBin;
	
	private List<Double> kaksRatios;
	
	public SelectonRunner(File selectonBin) {
		this.selectonBin = selectonBin;
	}
	
	/**
	 * Runs selecton parsing the output file. Use {@link #getKaKsRatios()} to get the results.
	 * @param inputAlnFile input nucleotide alignment file in FASTA format, must contain the refSequenceName
	 * @param resultsFile selecton results output file
	 * @param logFile selecton log output file, if null it will be written to a temp file and removed 
	 * @param treeFile selecton tree output file, if null it will be written to a temp file and removed
	 * @param colorBinFile selecton color bin output file, if null it will be written to a temp file and removed
	 * @param globalResultsFile selecton global output file, if null it will be written to a temp file and removed
	 * @param refSequenceName the name (FASTA tag) of the reference sequence for which we want ka/ks results
	 * @param gcType the GeneticCodeType to be used by selecton
	 * @param epsilon the epsilon value for selecton's likelihood optimization. The smaller the higher precision
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void run(File inputAlnFile, File resultsFile, File logFile, File treeFile, File colorBinFile, File globalResultsFile, String refSequenceName, GeneticCodeType gcType, double epsilon) 
	throws IOException, InterruptedException {
		String prefix = "selecton";
		if (logFile==null) {
			logFile = File.createTempFile(prefix, ".log");
			logFile.deleteOnExit();
		}
		if (treeFile==null) {
			treeFile = File.createTempFile(prefix, ".tree");
			treeFile.deleteOnExit();
		}
		if (colorBinFile==null) {
			colorBinFile = File.createTempFile(prefix, ".colorBins");
			colorBinFile.deleteOnExit();
		}
		if (globalResultsFile==null) {
			globalResultsFile = File.createTempFile(prefix, ".global");
			globalResultsFile.deleteOnExit();
		}
		String[] cmdLine = {selectonBin.toString(), 
		"-q",refSequenceName,
		"-i",inputAlnFile.toString(),
		"-e",String.format("%4.2f",epsilon),
		"-g",String.valueOf(gcType.getSelectonId()),
		"-r",resultsFile.toString(),
		"-l",logFile.toString(),
		"-t",treeFile.toString(),
		"-o",globalResultsFile.toString(),
		"-c",colorBinFile.toString()};
		Process selectonProc = Runtime.getRuntime().exec(cmdLine);

		int exitValue = selectonProc.waitFor();
		if (exitValue>0) {
			throw new IOException("Selecton exited with error value " + exitValue);
		}
		try {
			parseResultsFile(resultsFile, refSequenceName);
		} catch(FileFormatException e) {
			throw new IOException("Selecton output file "+resultsFile+" is not in the expected format: "+e.getMessage());
		}
	}
	
	/**
	 * Parses selecton results output file (the one obtained with the -r option). Get the 
	 * ka/ks values by calling {@link #getKaKsRatios()}
	 * @param resultsFile
	 * @param sequenceName the name (FASTA tag) of the reference sequence passed with the
	 * -q option to selecton
	 * @throws IOException
	 * @throws FileFormatException if sequence serials not sequential, unknown aminoacid code 
	 * or sequence reference name given not matching the one in the file
	 */
	public void parseResultsFile(File resultsFile, String sequenceName) throws IOException, FileFormatException {
		kaksRatios = new ArrayList<Double>();
		BufferedReader br = new BufferedReader(new FileReader(resultsFile));
		String line;
		Pattern seqNamePat = Pattern.compile("^Displayed\\son\\ssequence\\s(.*)$");
		String seqName = null;
		int lineCount = 0;
		int lastSerial = 0;
		while ((line=br.readLine())!=null) {
			lineCount++;
			if (lineCount<=6) {
				Matcher m = seqNamePat.matcher(line);
				if (m.matches()) {
					seqName = m.group(1).trim();
					if (!seqName.equals(sequenceName)) {
						throw new FileFormatException("The reference sequence name in output selecton file "+resultsFile+" ("+seqName+") does not coincide with the expected sequence name "+sequenceName);
					}
				}
				continue;
			}
			String[] tokens = line.split("\\s+");
			int serial = Integer.parseInt(tokens[0]);
			if (lastSerial+1!=serial) {
				throw new FileFormatException("Aminoacids in the selecton results file "+resultsFile+" are not sequentially numbered at line "+lineCount);
			}
			AminoAcid aa = AminoAcid.getByOneLetterCode(tokens[1].toCharArray()[0]);
			if (aa==null) {
				throw new FileFormatException("Unknown amino acid one letter code in selecton results file "+resultsFile+" at line "+lineCount);
			}
			kaksRatios.add(Double.parseDouble(tokens[2]));
			lastSerial = serial;
		}
		br.close();
		if (seqName==null) {
			throw new FileFormatException("Could not find a sequence name in selecton results file "+resultsFile);
		}
	}
	
	public List<Double> getKaKsRatios() {
		return kaksRatios;
	}
}
