package tinker;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.channels.FileChannel;
import java.util.Formatter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class TinkerRunner {
	
	private static final String PROTEIN_PROG = "protein";
	private static final String DISTGEOM_PROG = "distgeom";
	private static final String PDBXYZ_PROG = "pdbxyz";
	private static final String XYZPDB_PROG = "xyzpdb";
	private static final String CYCLISE_PROTEIN_STR = "N";
	private static final String DGEOM_PARAMS = "Y N Y Y N N A";
	private static final String TINKER_ERROR_STR = " TINKER is Unable to Continue";
	
	private String tinkerBinDir;
	
	private String forceFieldFileName;
	
	private String proteinProg;
	private String distgeomProg;
	private String pdbxyzProg;
	private String xyzpdbProg;
	
	private File logFile;
	private PrintWriter log;
	
	// arrays for storing distgeom output data
	private int[] numUpperBoundViol;
	private int[] numLowerBoundViol;
	private double[] maxUpperBoundViol;
	private double[] maxLowerBoundViol;
	private double[] rmsBoundViol;
	private int[] numUpperViol;
	private int[] numLowerViol;
	private double[] maxUpperViol;
	private double[] maxLowerViol;
	private double[] rmsRestViol;

	
	/**
	 * Constructs a TinkerRunner object by passing initial parameters
	 * @param tinkerBinDir The directory where the tinker executables are
	 * @param forceFieldFileName The force field file
	 * @param logFile File where all tinker output will be logged to
	 * @throws FileNotFoundException If logFile can't be written
	 */
	public TinkerRunner(String tinkerBinDir, String forceFieldFileName, File logFile) throws FileNotFoundException {
		this.tinkerBinDir = tinkerBinDir;
		this.forceFieldFileName = forceFieldFileName;
		this.proteinProg = new File(this.tinkerBinDir,PROTEIN_PROG).getAbsolutePath();
		this.distgeomProg = new File(this.tinkerBinDir,DISTGEOM_PROG).getAbsolutePath();
		this.pdbxyzProg = new File(this.tinkerBinDir,PDBXYZ_PROG).getAbsolutePath();
		this.xyzpdbProg = new File(this.tinkerBinDir,XYZPDB_PROG).getAbsolutePath();
		
		this.logFile = logFile;
		this.log = new PrintWriter(new FileOutputStream(logFile));
	}
	
	/**
	 * To get the expected File that a tinker program will output given an input file and an extension for the output files
	 * The directory where the input file is will be scanned to see if it contains files of the form basename.ext, basename.ext_2, basename.ext_3 etc.
	 * @param file
	 * @param ext
	 * @return
	 */
	private File getTinkerOutputFileName(File file, String ext){
		String basename = file.getName();
		basename = basename.substring(0, basename.lastIndexOf("."));
		String dirname = file.getParent();
		
		String tinkerOutFileName = basename + "." + ext;
		
		if (new File(dirname,tinkerOutFileName).exists()) {
			int i = 2;
			tinkerOutFileName = basename + "." + ext + "_" + i;
			while (new File(dirname,tinkerOutFileName).exists()) {
				i++;
				tinkerOutFileName = basename + "." + ext + "_" + i;
			}
		}
		return new File(dirname,tinkerOutFileName);
	}
	
	/**
	 * Runs tinker's protein program to generate an elongated protein structure given a sequence
	 * @param sequence
	 * @param outPath The directory where output files will be written
	 * @param outBasename The base name for the output files
	 * @throws IOException
	 * @throws TinkerError
	 */
	public void runProtein(String sequence, String outPath, String outBasename) throws IOException, TinkerError {
		boolean tinkerError = false; // to store the exit state of the tinker program
		
		if (!new File(outPath).exists()) {
			throw new FileNotFoundException("Specified directory "+outPath+" does not exist");
		}
		File tinkerxyzout = getTinkerOutputFileName(new File(outPath,outBasename+".xyz"),"xyz");
		File tinkerintout = getTinkerOutputFileName(new File(outPath,outBasename+".int"),"int");
		tinkerintout.deleteOnExit();
		File tinkerseqout = getTinkerOutputFileName(new File(outPath,outBasename+".seq"),"seq");
		
		// running protein program in outPath dir (so that output files are written to outPath)
		Process protProc = Runtime.getRuntime().exec(proteinProg, null, new File(outPath));
		// piping input
		PrintWriter protInput = new PrintWriter(protProc.getOutputStream());
		protInput.println(outBasename);
		protInput.println("Unfolded chain created by tinker's protein program");
		protInput.println(forceFieldFileName);
		for (int i=0;i<sequence.length();i++) {
			// we've got to use 3 letter code for CYS, otherwise tinker takes the default CYX (which means cystein with disulfide bridge)
			if (sequence.charAt(i)=='C') {
				protInput.println("CYS");
			} else {
				protInput.println(sequence.charAt(i));
			}
		}
		protInput.println();
		protInput.println(CYCLISE_PROTEIN_STR);
		protInput.close();
		
		// logging output
		BufferedReader protOutput = new BufferedReader(new InputStreamReader(protProc.getInputStream()));
		String line;
		while((line = protOutput.readLine()) != null) {
			log.println(line);
			if (line.startsWith(TINKER_ERROR_STR)) {
				tinkerError = true;
			}
		}
		
		tinkerxyzout.renameTo(new File(outPath,outBasename+".xyz"));
		tinkerseqout.renameTo(new File(outPath,outBasename+".seq"));
		
		if (tinkerError) {
			log.close();
			throw new TinkerError("Tinker error, revise log file "+logFile.getAbsolutePath());
		}
	}

	/**
	 * Runs tinker's distgeom program capturing output with restrain violation statistics into member variable arrays
	 * that can be retrieved using the getters: getMaxLowerBoundViol, getMaxUpperBoundViol, getMaxLowerViol etc... 
	 * Two files are needed as input for distgeom: an xyz file and a key file, the latter is not passed but instead implicitely 
	 * defined by xyzFile: must be in same directory and must have same basename with extension .key 
	 * @param xyzFile
	 * @param outPath Directory where output files will be written
	 * @param outBasename Base name of the output files
	 * @param n Number of models that we want distgeom to produce
	 * @throws TinkerError If an error seen in tinker's output
	 * @throws IOException
	 */
	public void runDistgeom(File xyzFile, String outPath, String outBasename, int n) throws TinkerError, IOException {
		boolean tinkerError = false; // to store the exit state of the tinker program
		if (!new File(outPath).exists()) {
			throw new FileNotFoundException("Specified directory "+outPath+" does not exist");
		}
		if (!xyzFile.exists()){
			throw new FileNotFoundException("Specified xyz file "+xyzFile.getAbsolutePath()+" does not exist");
		}
		String basename = xyzFile.getName();
		basename = basename.substring(0, basename.lastIndexOf("."));
		File keyFile = new File(xyzFile.getParent(),basename+".key");
		if (! keyFile.exists()) {
			throw new FileNotFoundException("Key file "+keyFile.getAbsolutePath()+" not present in input directory "+xyzFile.getParent());
		}
		// getting names of tinker output files
		File[] tinkerout = new File[n+1];
		for (int i=1;i<=n;i++) {
			String ext = new Formatter().format("%03d", i).toString();
			tinkerout[i] = getTinkerOutputFileName(xyzFile, ext);
		}
		// initialising arrays were we store captured output data
		numUpperBoundViol = new int[n+1];
		numLowerBoundViol = new int[n+1];
		maxUpperBoundViol = new double[n+1];
		maxLowerBoundViol = new double[n+1];
		rmsBoundViol = new double[n+1];
		numUpperViol = new int[n+1];
		numLowerViol = new int[n+1];
		maxUpperViol = new double[n+1];
		maxLowerViol = new double[n+1];
		rmsRestViol = new double[n+1];
		// running distgeom program
		Process dgeomProc = Runtime.getRuntime().exec(distgeomProg+" "+xyzFile.getAbsolutePath()+" "+n+" "+DGEOM_PARAMS);
		// logging and capturing output
		BufferedReader dgeomOutput = new BufferedReader(new InputStreamReader(dgeomProc.getInputStream()));
		String line;
		int i=1;
		while((line = dgeomOutput.readLine()) != null) {
			log.println(line);
			if (line.startsWith(TINKER_ERROR_STR)) {
				tinkerError = true;
			}
			Pattern p = Pattern.compile("^ Num Upper Bound Violations :\\s+(\\d+)");
			Matcher m = p.matcher(line);
			if (m.find()) {
				numUpperBoundViol[i]=Integer.parseInt(m.group(1));
			}
			p = Pattern.compile("^ Num Lower Bound Violations :\\s+(\\d+)");
			m = p.matcher(line);
			if (m.find()) {
				numLowerBoundViol[i]=Integer.parseInt(m.group(1));
			}
			p = Pattern.compile("^ Max Upper Bound Violation :\\s+(\\d+\\.\\d\\d\\d\\d)");
			m = p.matcher(line);
			if (m.find()) {
				maxUpperBoundViol[i]=Double.parseDouble(m.group(1));
			}
			p = Pattern.compile("^ Max Lower Bound Violation :\\s+(\\d+\\.\\d\\d\\d\\d)");
			m = p.matcher(line);
			if (m.find()) {
				maxLowerBoundViol[i]=Double.parseDouble(m.group(1));
			}
			p = Pattern.compile("^ RMS Deviation from Bounds :\\s+(\\d+\\.\\d\\d\\d\\d)");
			m = p.matcher(line);
			if (m.find()) {
				rmsBoundViol[i]=Double.parseDouble(m.group(1));
			}
			p = Pattern.compile("^ Num Upper Restraint Violations :\\s+(\\d+)");
			m = p.matcher(line);
			if (m.find()) {
				numUpperViol[i]=Integer.parseInt(m.group(1));
			}
			p = Pattern.compile("^ Num Lower Restraint Violations :\\s+(\\d+)");
			m = p.matcher(line);
			if (m.find()) {
				numLowerViol[i]=Integer.parseInt(m.group(1));
			}
			p = Pattern.compile("^ Max Upper Restraint Violation :\\s+(\\d+\\.\\d\\d\\d\\d)");
			m = p.matcher(line);
			if (m.find()) {
				maxUpperViol[i]=Double.parseDouble(m.group(1));
			}
			p = Pattern.compile("^ Max Lower Restraint Violation :\\s+(\\d+\\.\\d\\d\\d\\d)");
			m = p.matcher(line);
			if (m.find()) {
				maxLowerViol[i]=Double.parseDouble(m.group(1));
			}
			p = Pattern.compile("^ RMS Restraint Dist Violation :\\s+(\\d+\\.\\d\\d\\d\\d)");
			m = p.matcher(line);
			if (m.find()) {
				rmsRestViol[i]=Double.parseDouble(m.group(1));
				i++;
			}
		}
		//renaming files to our chosen outBasename+ext
		for (i=1;i<=n;i++) {
			String ext = new Formatter().format("%03d", i).toString();
			tinkerout[i].renameTo(new File(outPath,outBasename+"."+ext));
		}
		// throwing exception if error string was caught in output
		if (tinkerError) {
			log.close();
			throw new TinkerError("Tinker error, revise log file "+logFile.getAbsolutePath());
		}
		
	}
	
	/**
	 * Runs tinker's xyzpdb program to convert a given xyzFile (needing also a seqFile) to a pdbFile
	 * @param xyzFile
	 * @param seqFile
	 * @param pdbFile
	 * @throws IOException 
	 * @throws TinkerError If an error seen in tinker's output
	 */
	public void runXyzpdb(File xyzFile, File seqFile, File pdbFile) throws IOException, TinkerError {
		boolean tinkerError = false; // to store the exit state of the tinker program
		if (!xyzFile.exists()){
			throw new FileNotFoundException("Specified xyz file "+xyzFile.getAbsolutePath()+" does not exist");
		}
		if (!seqFile.exists()){
			throw new FileNotFoundException("Specified seq file "+seqFile.getAbsolutePath()+" does not exist");
		}
		
		String basename = xyzFile.getName();
		basename = basename.substring(0, basename.lastIndexOf("."));
		File tmpSeqFile = new File(seqFile.getParent(),basename+".seq");
				
		// if seqFile doesn't follow the naming convention (basename of xyzFile+seq extension) that tinker expects, we copy it to tmpSeqFile (which has right name) 
		if (!tmpSeqFile.equals(seqFile)) {
			FileChannel srcChannel = new FileInputStream(seqFile).getChannel();
			FileChannel dstChannel = new FileOutputStream(tmpSeqFile).getChannel();
			dstChannel.transferFrom(srcChannel, 0, srcChannel.size());
			srcChannel.close();
			dstChannel.close();
			// if we copied then that means tmpSeqFile is different from seqFile and thus we want to delete the tmp file on exit
			tmpSeqFile.deleteOnExit(); 
		}
        
		File tinkerpdbout = getTinkerOutputFileName(xyzFile, "pdb");
		
		// running tinker's xyzpdb
		// beware: it takes as a silent input the seq file seqFile (or tmpSeqFile if the input seqFile didn't have the right name) 
		Process xyzpdbProc = Runtime.getRuntime().exec(xyzpdbProg+" "+xyzFile.getAbsolutePath()+" "+forceFieldFileName);

		// logging output
		BufferedReader xyzpdbOutput = new BufferedReader(new InputStreamReader(xyzpdbProc.getInputStream()));
		String line;
		while((line = xyzpdbOutput.readLine()) != null) {
			log.println(line);
			if (line.startsWith(TINKER_ERROR_STR)) {
				tinkerError = true;
			}
		}
		
		tinkerpdbout.renameTo(pdbFile);
		
		if (tinkerError) {
			log.close();
			throw new TinkerError("Tinker error, revise log file "+logFile.getAbsolutePath());
		}
	}
	
	/**
	 * Runs tinker's pdbxyz program to convert a pdbFile to a xyzFile
	 * @param pdbFile
	 * @param xyzFile
	 * @throws IOException
	 * @throws TinkerError If an error seen in tinker's output
	 */
	public void runPdbxyz(File pdbFile, File xyzFile) throws IOException, TinkerError{
		boolean tinkerError = false; // to store the exit state of the tinker program
		if (!pdbFile.exists()){
			throw new FileNotFoundException("Specified pdb file "+pdbFile.getAbsolutePath()+" does not exist");
		}
		File tinkerxyzout = getTinkerOutputFileName(pdbFile, "xyz");
		// running tinker's pdbxyz
		Process pdbxyzProc = Runtime.getRuntime().exec(pdbxyzProg+" "+pdbFile.getAbsolutePath()+" "+forceFieldFileName);

		// logging output
		BufferedReader pdbxyzOutput = new BufferedReader(new InputStreamReader(pdbxyzProc.getInputStream()));
		String line;
		while((line = pdbxyzOutput.readLine()) != null) {
			log.println(line);
			if (line.startsWith(TINKER_ERROR_STR)) {
				tinkerError = true;
			}
		}
		tinkerxyzout.renameTo(xyzFile);
		if (tinkerError) {
			log.close();
			throw new TinkerError("Tinker error, revise log file "+logFile.getAbsolutePath());
		}
	}
	
	/**
	 * Closes log stream, must be called after no other tinker program will be run with this TinkerRunner object
	 * (otherwise log is not flushed to file)
	 */
	public void closeLog() {
		log.close();
	} 
	
	public double[] getMaxLowerBoundViol() {
		return maxLowerBoundViol;
	}

	public double[] getMaxLowerViol() {
		return maxLowerViol;
	}

	public double[] getMaxUpperBoundViol() {
		return maxUpperBoundViol;
	}

	public double[] getMaxUpperViol() {
		return maxUpperViol;
	}

	public int[] getNumLowerBoundViol() {
		return numLowerBoundViol;
	}

	public int[] getNumLowerViol() {
		return numLowerViol;
	}

	public int[] getNumUpperBoundViol() {
		return numUpperBoundViol;
	}

	public int[] getNumUpperViol() {
		return numUpperViol;
	}

	public double[] getRmsBoundViol() {
		return rmsBoundViol;
	}

	public double[] getRmsRestViol() {
		return rmsRestViol;
	}
}
