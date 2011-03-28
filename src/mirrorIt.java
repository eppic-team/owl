import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.structure.Pdb;
import owl.core.structure.PdbasePdb;
import owl.core.structure.PdbfilePdb;



public class mirrorIt {

	private static final String PROG_NAME = "mirrorIt";
	
	public static void main(String[] args) throws Exception {

		String help = "\nMirrors the given pdb structure or file\n" +
					"Usage:\n" +
					PROG_NAME+" <pdbCode+pdbChainCode/pdb file name> <out_file>\n\n" +
							"example: "+ PROG_NAME +" 1abcA mirrored.pdb\n\n";

		if (args.length<2) {
			System.err.println("Missing argument");
			System.err.println(help);
			System.exit(1);
		}
		String pdbId = args[0];
		File outFile = new File(args[1]);
		
		String pdbCode = "";
		String pdbChainCode = "";
		boolean pdbFromFile = false;
		
		// pdb code and chain code
		Pattern p = Pattern.compile("(\\d\\w\\w\\w)(\\w)");
		Matcher m = p.matcher(pdbId);
		if (m.matches()) {
			pdbCode = m.group(1);
			pdbChainCode = m.group(2);
			pdbFromFile = false;
		} 
		else if (new File(pdbId).exists()) {
			pdbFromFile = true;
		}
		else {				
			System.err.println("Either PDB id given not in correct format or it is not an existing PDB file: "+pdbId);
			System.exit(1);
		}		

		Pdb pdb = null;
		
		if (pdbFromFile) {
			pdb = new PdbfilePdb(pdbId);
			pdbChainCode = pdb.getChains()[0];
			pdb.load(pdbChainCode);
		} else {
			pdb = new PdbasePdb(pdbCode);
			pdb.load(pdbChainCode);
		}
		
		pdb.mirror();
		pdb.writeToPDBFile(outFile.getAbsolutePath());
	}

}
