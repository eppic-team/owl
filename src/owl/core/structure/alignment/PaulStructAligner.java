package owl.core.structure.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;

import owl.core.sequence.Sequence;
import owl.core.sequence.alignment.AlignmentConstructionException;
import owl.core.sequence.alignment.MultipleSequenceAlignment;
import owl.core.structure.Pdb;
import owl.core.structure.Template;
import owl.core.structure.TemplateList;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.FileFormatException;
import owl.core.util.MySQLConnection;


/**
 * A structural aligner with paul: performs multiple structure alignment
 * using the paul program. See https://www.mi.fu-berlin.de/w/LiSA/Paul
 * @author duarte
 *
 */
public class PaulStructAligner implements StructAligner {
	
	private static final boolean DEBUG = false;
	
	public static final String  PAULMODE_VERYFAST = "veryfast";
	public static final String  PAULMODE_FAST = "fast";
	public static final String  PAULMODE_SLOW = "slow";
	public static final String  PAULMODE_VERYSLOW = "veryslow";
	
	private static final String DEFAULT_CONTACT_TYPE = "Cb";
	private static final double DEFAULT_CUTOFF = 8;
	
	private File paulProg;
	private File logFile;
	private PrintWriter paulLog;
	private String tempBaseName;
	private String paulMode;
	
	private String contactType;
	private double cutoff;
	
	/**
	 * Constructs a new PaulStructAligner
	 * @param paulProg
	 * @param paulMode one of the constants {@link #PAULMODE_VERYFAST}, {@link #PAULMODE_FAST}, 
	 * {@link #PAULMODE_SLOW} or {@link #PAULMODE_VERYSLOW}
	 * @param logFile where paul's output/erro will be logged
	 * @throws FileNotFoundException if logFile can't be found
	 */
	public PaulStructAligner(File paulProg, String paulMode, File logFile) throws FileNotFoundException {
		this.paulProg = paulProg;
		this.logFile = logFile;
		this.paulLog = new PrintWriter(this.logFile);
		this.tempBaseName = "temp"+System.currentTimeMillis();
		this.paulMode = paulMode;
		
		this.contactType = DEFAULT_CONTACT_TYPE;
		this.cutoff = DEFAULT_CUTOFF;
	}

	public MultipleSequenceAlignment alignStructures(TemplateList templates) throws StructAlignmentException, IOException {
		if (!templates.isGraphDataLoaded()) {
			throw new IllegalArgumentException("Given TemplateList does not have graph data loaded");
		}
		RIGraph[] graphs = new RIGraph[templates.size()];
		String[] tags = new String[templates.size()];
		int i = 0;
		for (Template template: templates) {
			graphs[i] = template.getRIGraph();
			tags[i] = template.getId();
			i++;
		}
		return alignStructures(graphs, tags);
	}

	/**
	 * Structurally aligns the given pdbs
	 * Uses the default contact type and cutoff ({@link #DEFAULT_CONTACT_TYPE}, {@link #DEFAULT_CUTOFF}) 
	 * to get the contact maps for the paul structural alignment.  
	 * @param pdbs
	 * @param tags
	 * @return
	 * @throws StructAlignmentException if a problem occurs while running the structural alignment
	 * @throws IOException
	 */
	public MultipleSequenceAlignment alignStructures(Pdb[] pdbs, String[] tags) throws StructAlignmentException, IOException {
		if (pdbs.length!=tags.length) {
			throw new IllegalArgumentException("Given pdbs and tags have different sizes (pdbs "+pdbs.length+", tags "+tags.length+").");
		}
		RIGraph[] graphs = new RIGraph[pdbs.length];
		for (int i=0;i<pdbs.length;i++) {
			graphs[i]=pdbs[i].getRIGraph(contactType, cutoff);
		}
		return alignStructures(graphs,tags);
	}
	
	public MultipleSequenceAlignment alignStructures(RIGraph[] graphs, String[] tags) throws StructAlignmentException, IOException {
		if (graphs.length!=tags.length) {
			throw new IllegalArgumentException("Given graphs and tags have different sizes (graphs "+graphs.length+", tags "+tags.length+").");
		}
		// we have to write all temp files to current dir: paul just doesn't work otherwise
		File cmListFile = new File(tempBaseName+".paul"); 
		if (!DEBUG) cmListFile.deleteOnExit();
		File seqListFile = new File(tempBaseName+".paul.seqs");
		if (!DEBUG) seqListFile.deleteOnExit();
		PrintWriter pwc = new PrintWriter(cmListFile);
		PrintWriter pws = new PrintWriter(seqListFile);
		for (int i=0;i<graphs.length;i++) {
			File cmFile = new File(tags[i]+"."+ tempBaseName+".cm");
			if (!DEBUG) cmFile.deleteOnExit();
			pwc.println(cmFile.getName());
			File seqFile = new File(tags[i]+"."+ tempBaseName+ ".cm.seq");
			pws.println(seqFile.getName());
			if (!DEBUG) seqFile.deleteOnExit();
			String[] seq = {graphs[i].getSequence()};
			String[] tag4seq = {tags[i]};
			Sequence.writeSeqs(seqFile, seq, tag4seq);
			graphs[i].writeToPaulFile(cmFile.getAbsolutePath());
		}
		pwc.close();
		pws.close();
		MultipleSequenceAlignment al = runPaul(cmListFile, seqListFile);
		// resetting tags from paul's output (uses the file name of the cm file as tags) to our original tags
		for (int i=0;i<tags.length;i++) {
			al.resetTag(tags[i]+"."+ tempBaseName, tags[i]);
		}
		return al; 
	}
	
	private MultipleSequenceAlignment runPaul(File cmListFile, File seqListFile) throws StructAlignmentException, IOException {
		String cmdLine = paulProg+" -i "+cmListFile+" -s "+seqListFile+" -m "+paulMode;
				
		Process proc = Runtime.getRuntime().exec(cmdLine);
		// logging and capturing stdout/stderr
		BufferedReader stdout = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		BufferedReader stderr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));

		String line;

		paulLog.println("#cmd: "+cmdLine);
		paulLog.println("#################");
		paulLog.println("# paul stdout ");
		paulLog.println("#################");
		while((line = stdout.readLine()) != null) {
			paulLog.println(line);
		}
		paulLog.println("#################");
		paulLog.println("# paul stderr ");
		paulLog.println("#################");
		while((line = stderr.readLine()) != null) {
			paulLog.println(line);
		}

		try {
			int exitValue = proc.waitFor();
			// throwing exception if exit state is not 0 
			if (exitValue!=0) {
				paulLog.flush();
				throw new StructAlignmentException("paul exited with value "+exitValue+". Revise log file "+logFile);
			}
		} catch (InterruptedException e) {
			System.err.println("Unexpected error while waiting for paul to exit. Error: "+e.getMessage());
			System.exit(1);
		}
		String base = cmListFile.getName();
		File outAliClustal = new File(base+".prog.aln");
		MultipleSequenceAlignment al = null;
		try {
			al = new MultipleSequenceAlignment(outAliClustal.getAbsolutePath(),MultipleSequenceAlignment.CLUSTALFORMAT);
		} catch (AlignmentConstructionException e) {
			this.cleanUp(base);
			this.paulLog.flush();
			throw new StructAlignmentException(e);
		}
		catch (FileFormatException e) {
			this.cleanUp(base);
			this.paulLog.flush();
			throw new StructAlignmentException(e);
		}
		
		this.cleanUp(base);
		this.paulLog.close();
		
		return al;
	}
	
	private void cleanUp(String base) {
		if (!DEBUG) {
			new File(base+".prog.aln").delete();
			new File(base+".prog.dnd").delete();
			new File(base+".prog.lib").delete();
			new File(base+".prog.mistral").delete();
			new File(base+".prog.results").delete();
			new File(base+".pw.aln").delete();
			new File(base+".pw.dnd").delete();
			new File(base+".pw.lib").delete();
			new File(base+".pw.mistral").delete();
			new File(base+".pw.results").delete();
		}
	}
	
	/**
	 * testing
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		File paulProg = new File("/project/StruPPi/bin/paul");
		File templatesFile = new File("/project/StruPPi/OncoGenesHL/modeling/results/FBXW7_WT/P2/FBXW7_WT.blast.templates");
		File logFile = new File("/scratch/local/temp/tmp2/paul.log");
		
		MySQLConnection conn = new MySQLConnection("talyn", "pdbase");
		
		TemplateList templates = new TemplateList(templatesFile);
		templates.loadPdbData(conn, "pdbase");
		templates.loadRIGraphs("Cb", 8.0);
		PaulStructAligner pr = new PaulStructAligner(paulProg, PaulStructAligner.PAULMODE_VERYFAST, logFile);
		MultipleSequenceAlignment al = pr.alignStructures(templates);
		al.printFasta();
		
	}
	
}
