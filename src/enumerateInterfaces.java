import gnu.getopt.Getopt;

import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.structure.ChainInterface;
import owl.core.structure.ChainInterfaceList;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.SpaceGroup;
import owl.core.util.MySQLConnection;


public class enumerateInterfaces {

	private static final File NACCESS_EXE = new File("/home/duarte_j/bin/naccess");

	private static final Pattern  PDBCODE_PATTERN = Pattern.compile("^\\d\\w\\w\\w$");
	
	// 6.0 seems to be PISA's cutoff, found for structure 1pmm where with 5.5 there is one interface (tiny, 1 atom contacting) missing
	// 5.0  gives 25 for 1pmo (right number) but 6.0 gives one more (26) for which NACCESS measures a negative area...
	// what's the cutoff then? I'm trying a value in between but it seems strange to choose such fractional values
	// 5.75 gives 25 for 1pmo (right) but 26 for 1pmm (instead of 27)
	// 5.90 gives 25 for 1pmo (right)  and 27 for 1pmm (right)  
	private static final double CUTOFF = 5.9; 
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {

		
		String help = 
			"Usage: enumerateInterfaces -i <pdb code> [-w <out dir for pdb files>]\n" +
			"If -w specified PDB files for each interface will be written to given out dir\n\n";
		
		String pdbStr = null;
		File writeDir = null;

		Getopt g = new Getopt("enumerateInterfaces", args, "i:w:h?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'i':
				pdbStr = g.getOptarg();
				break;
			case 'w':
				writeDir = new File(g.getOptarg());
				break;
			case 'h':
			case '?':
				System.out.println(help);
				System.exit(0);
				break; // getopt() already printed an error
			}
		}

		if (pdbStr == null) {
			System.err.println("Missing input PDB code/file (-i)");
			System.err.println(help);
			System.exit(1);
		}
		
		File inputFile = new File(pdbStr);
		Matcher m = PDBCODE_PATTERN.matcher(pdbStr);
		if (m.matches()) {
			inputFile = null;
		}
		
		if (inputFile!=null && !inputFile.exists()){
			System.err.println("Given file "+inputFile+" does not exist!");
			System.exit(1);
		}

		String outBaseName = pdbStr;
		
		PdbAsymUnit pdb = null;
		if (inputFile==null) {
			pdb = new PdbAsymUnit(pdbStr, new MySQLConnection(), "pdbase");
		} else {
			pdb = new PdbAsymUnit(inputFile);
			outBaseName = inputFile.getName().substring(0, inputFile.getName().lastIndexOf("."));
		}

		System.out.println(pdb.getSpaceGroup().getShortSymbol()+" ("+pdb.getSpaceGroup().getId()+")");
		
		System.out.println("Calculating possible interfaces...");
		long start = System.currentTimeMillis();
		ChainInterfaceList interfaces = pdb.getAllInterfaces(CUTOFF, NACCESS_EXE);
		long end = System.currentTimeMillis();
		System.out.println("Done. Time "+(end-start)/1000+"s");
		
		System.out.println("Total number of interfaces found: "+interfaces.size());

					
		for (int i=0;i<interfaces.size();i++) {
			ChainInterface interf = interfaces.get(i);
			System.out.println("\n##Interface "+(i+1));
			System.out.println("Transf1: "+SpaceGroup.getAlgebraicFromMatrix(interf.getFirstTransf())+
					". Transf2: "+SpaceGroup.getAlgebraicFromMatrix(interf.getSecondTransf()));
			System.out.println(interf.getFirstMolecule().getPdbChainCode()+" - "+interf.getSecondMolecule().getPdbChainCode());
			System.out.println("Number of contacts: "+interf.getNumContacts());
			System.out.println("Number of contacting atoms (from both molecules): "+interf.getNumAtomsInContact());
			System.out.printf("Interface area: %8.2f\n",interf.getInterfaceArea());
		
			if (writeDir!=null) {
				interf.writeToPdbFile(new File(writeDir,outBaseName+"."+(i+1)+".interface.pdb"));
			}
		}
	}

}
