import gnu.getopt.Getopt;

import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.runners.PymolRunner;
import owl.core.structure.Asa;
import owl.core.structure.ChainInterface;
import owl.core.structure.ChainInterfaceList;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.SpaceGroup;


public class enumerateInterfaces {

	private static final String LOCAL_CIF_DIR = "/nfs/data/dbs/pdb/data/structures/all/mmCIF";
	private static final String BASENAME = "interf_enum";
	private static final String TMPDIR = System.getProperty("java.io.tmpdir");
	private static final File   PYMOL_EXE = new File("/usr/bin/pymol");
	private static final int[] HEIGHTS = {300};
	private static final int[] WIDTHS = {300};
	
	private static final Pattern  PDBCODE_PATTERN = Pattern.compile("^\\d\\w\\w\\w$");
	
	// 6.0 seems to be PISA's cutoff, found for structure 1pmm where with 5.5 there is one interface (tiny, 1 atom contacting) missing
	// 5.0  gives 25 for 1pmo (right number) but 6.0 gives one more (26) for which NACCESS measures a negative area...
	// what's the cutoff then? I'm trying a value in between but it seems strange to choose such fractional values
	// 5.75 gives 25 for 1pmo (right) but 26 for 1pmm (instead of 27)
	// 5.90 gives 25 for 1pmo (right)  and 27 for 1pmm (right)  
	private static final double CUTOFF = 5.9; 
	
	private static final double CLASH_DISTANCE = 1.5;
	
	private static final int NTHREADS = Runtime.getRuntime().availableProcessors();
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {

		
		String help = 
			"Usage: enumerateInterfaces -i <pdb code> [-w <out dir for pdb files>]\n" +
			"If -w specified PDB files for each interface will be written to given out dir\n" +
			"If -l specified cartoon PNG images of each interface will be written to given \n" +
			"out dir (must use -w also)\n\n";
		
		String pdbStr = null;
		File writeDir = null;
		int nThreads = NTHREADS;
		boolean generatePngs = false;

		Getopt g = new Getopt("enumerateInterfaces", args, "i:w:t:lh?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'i':
				pdbStr = g.getOptarg();
				break;
			case 'w':
				writeDir = new File(g.getOptarg());
				break;
			case 't':
				nThreads = Integer.parseInt(g.getOptarg());
				break;
			case 'l':
				generatePngs = true;
				break;
			case 'h':
				System.out.println(help);
				System.exit(0);
				break;
			case '?':
				System.err.println(help);
				System.exit(1);
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
		
		if (generatePngs==true && writeDir==null) {
			System.err.println("Can't generate images if a write directory not specified (use -w)");
			System.exit(1);
		}

		String outBaseName = pdbStr;
		
		PdbAsymUnit pdb = null;
		if (inputFile==null) {
			File cifFile = new File(TMPDIR,BASENAME+"_"+pdbStr+".cif");
			PdbAsymUnit.grabCifFile(LOCAL_CIF_DIR, null, pdbStr, cifFile, false);
			pdb = new PdbAsymUnit(cifFile);
		} else {
			pdb = new PdbAsymUnit(inputFile);
			outBaseName = inputFile.getName().substring(0, inputFile.getName().lastIndexOf("."));
		}

		System.out.println(pdb.getSpaceGroup().getShortSymbol()+" ("+pdb.getSpaceGroup().getId()+")");
		
		System.out.println("Calculating possible interfaces... (using "+nThreads+" CPUs for ASA calculation)");
		long start = System.currentTimeMillis();
		ChainInterfaceList interfaces = pdb.getAllInterfaces(CUTOFF, null, Asa.DEFAULT_N_SPHERE_POINTS, nThreads);
		long end = System.currentTimeMillis();
		System.out.println("Done. Time "+(end-start)/1000+"s");
		
		System.out.println("Total number of interfaces found: "+interfaces.size());

		PymolRunner pr = new PymolRunner(PYMOL_EXE);
					
		for (int i=0;i<interfaces.size();i++) {
			ChainInterface interf = interfaces.get(i);
			System.out.println("\n##Interface "+(i+1));
			if (interf.hasClashes(CLASH_DISTANCE)) System.out.println("CLASHES!!!");
			System.out.println("Transf1: "+SpaceGroup.getAlgebraicFromMatrix(interf.getFirstTransf())+
					". Transf2: "+SpaceGroup.getAlgebraicFromMatrix(interf.getSecondTransf()));
			System.out.println(interf.getFirstMolecule().getPdbChainCode()+" - "+interf.getSecondMolecule().getPdbChainCode());
			System.out.println("Number of contacts: "+interf.getNumContacts());
			System.out.println("Number of contacting atoms (from both molecules): "+interf.getNumAtomsInContact());
			System.out.printf("Interface area: %8.2f\n",interf.getInterfaceArea());
		
			if (writeDir!=null) {
				File pdbFile = new File(writeDir,outBaseName+"."+(i+1)+".interface.pdb");
				File[] outPngFiles = {new File(writeDir,outBaseName+"."+(i+1)+".interface.png")};
				interf.writeToPdbFile(pdbFile);
				if (generatePngs) {
					pr.generatePng(pdbFile, outPngFiles, "cartoon", "white", HEIGHTS, WIDTHS);
				}
			}
		}
	}

}
