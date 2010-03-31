package casp.tools;

import gnu.getopt.Getopt;

import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.structure.Alignment;
import owl.core.util.MySQLConnection;


public class DisplaySecStructMatching {
	
	private static final String PROGRAM_NAME = "displaySSMatching";
	
	private static final String PDBASE_DB = "pdbase";
	
	private static final String DSSP_EXE = "/project/StruPPi/bin/dssp";
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {

		File alnFile = null;
		File psipredFile = null;
		
		String help = "Usage: \n" +
		PROGRAM_NAME+"\n" +
		"   -a :  file with multiple sequence alignment of target and templates\n"+
		"   -s :  file with PsiPred prediction (horizontal format) of secondary structure for target\n"+
		"   -n :  no sequence, show only secondary structure assignments and supress sequence output\n\n";

		Getopt g = new Getopt(PROGRAM_NAME, args, "a:s:nh?");
		int c;
		boolean showSequence = true;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'a':
				alnFile = new File(g.getOptarg());
				break;
			case 's':
				psipredFile = new File(g.getOptarg());
				break;
			case 'n':
				showSequence = false;
				break;
			case 'h':
			case '?':
				System.out.println(help);
				System.exit(0);
				break; // getopt() already printed an error
			}
		}

		if (alnFile==null || psipredFile==null) {
			System.err.println("Missing option, you must specify an alignment file and a psipred horizontal file");
			System.out.println(help);
			System.exit(1);
		}
		
		MySQLConnection conn = new MySQLConnection();
		
		Alignment aln = new Alignment(alnFile.getAbsolutePath(),"FASTA");

		String targetTag = "";
		Pattern p = Pattern.compile("(\\w+)\\..*"); // i.e. we take the basename of the file (everything before the first dot)
		Matcher m = p.matcher(psipredFile.getName());
		if (m.matches()) {
			targetTag = m.group(1);
		} 
		
		if (!aln.hasTag(targetTag)) {
			System.err.println("Given alignment file "+alnFile+" doesn't contain the tag '"+targetTag+"' (parsed from the psipred horizontal file name "+psipredFile+")");
			System.exit(1);
		}
		
		aln.addSecStructAnnotation(conn, PDBASE_DB, DSSP_EXE);
		aln.writeWithSecStruct(System.out, targetTag, psipredFile, showSequence);
	}

}
