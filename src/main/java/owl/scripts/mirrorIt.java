package owl.scripts;
import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.structure.PdbAsymUnit;



public class mirrorIt {

	private static final String PROG_NAME = "mirrorIt";
	
	private static final String CIFREPODIR = "/path/to/mmCIF/gz/all/repo/dir";


	public static void main(String[] args) throws Exception {

		String help = "\nMirrors the given pdb structure or file\n" +
					"Usage:\n" +
					PROG_NAME+" <pdbCode/pdb file name> <out_file>\n\n" +
							"example: "+ PROG_NAME +" 1abc mirrored.pdb\n\n";

		if (args.length<2) {
			System.err.println("Missing argument");
			System.err.println(help);
			System.exit(1);
		}
		String pdbId = args[0];
		File outFile = new File(args[1]);
		
		String pdbCode = "";
		boolean pdbFromFile = false;
		
		// pdb code and chain code
		Pattern p = Pattern.compile("(\\d\\w\\w\\w)");
		Matcher m = p.matcher(pdbId);
		if (m.matches()) {
			pdbCode = m.group(1);
			pdbFromFile = false;
		} 
		else if (new File(pdbId).exists()) {
			pdbFromFile = true;
		}
		else {				
			System.err.println("Either PDB id given not in correct format or it is not an existing PDB file: "+pdbId);
			System.exit(1);
		}		

		PdbAsymUnit pdb = null;
		
		if (pdbFromFile) {
			pdb = new PdbAsymUnit(new File(pdbId));
		} else {
			
			File cifFile = new File(System.getProperty("java.io.tmpdir"),pdbCode+".cif");
			cifFile.deleteOnExit();
			PdbAsymUnit.grabCifFile(CIFREPODIR, null, pdbCode, cifFile, false);				
			pdb = new PdbAsymUnit(cifFile);
		}
		
		pdb.mirror();
		pdb.writeToPdbFile(outFile);
	}

}
