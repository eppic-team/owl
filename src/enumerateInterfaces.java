import gnu.getopt.Getopt;

import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.runners.PymolRunner;
import owl.core.structure.Asa;
import owl.core.structure.ChainInterface;
import owl.core.structure.ChainInterfaceList;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbChain;
import owl.core.structure.SpaceGroup;
import owl.core.util.Goodies;


public class enumerateInterfaces {

	private static final String LOCAL_CIF_DIR = "/nfs/data/dbs/pdb/data/structures/all/mmCIF";
	private static final String BASENAME = "interf_enum";
	private static final String TMPDIR = System.getProperty("java.io.tmpdir");
	private static final File   PYMOL_EXE = new File("/usr/bin/pymol");
	
	private static final double BSATOASA_CUTOFF = 0.95;

	private static final Pattern  PDBCODE_PATTERN = Pattern.compile("^\\d\\w\\w\\w$");
	
	// 6.0 seems to be PISA's cutoff, found for structure 1pmm where with 5.5 there is one interface (tiny, 1 atom contacting) missing
	// 5.0  gives 25 for 1pmo (right number) but 6.0 gives one more (26) for which NACCESS measures a negative area...
	// what's the cutoff then? I'm trying a value in between but it seems strange to choose such fractional values
	// 5.75 gives 25 for 1pmo (right) but 26 for 1pmm (instead of 27)
	// 5.90 gives 25 for 1pmo (right)  and 27 for 1pmm (right)  
	private static final double CUTOFF = 5.9; 
	
	private static final int NTHREADS = Runtime.getRuntime().availableProcessors();
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {

		
		String help = 
			"Usage: \n" +
			"enumerateInterfaces \n" +
			" -i <string> : input pdb code\n" +
			" [-t <int>]  : number of threads for calculating ASAs. Default: "+NTHREADS + "\n"+
			" [-w <dir>]  : output dir to write PDB files for each interface \n" +
			" [-l]        : cartoon PNG images of each interface will be written to\n" +
			"               output dir given in -w \n" +
			" [-s]        : write a serialized interfaces.dat file to output dir given\n" +
			"               in -w\n" +
			" [-n]        : non-polymer chains will also be considered\n" +
			" [-d]        : more verbose output for debugging\n\n";
		
		String pdbStr = null;
		File writeDir = null;
		int nThreads = NTHREADS;
		boolean generatePngs = false;
		boolean debug = false;
		boolean serialize = false;
		boolean nonPoly = false;

		Getopt g = new Getopt("enumerateInterfaces", args, "i:w:t:lsndh?");
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
			case 's':
				serialize = true;
				break;
			case 'n':
				nonPoly = true;
				break;
			case 'd':
				debug = true;
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

		// we remove H atoms
		pdb.removeHatoms();
		
		System.out.println(pdb.getPdbCode()+" - "+pdb.getNumPolyChains()+" polymer chains ("+pdb.getAllRepChains().size()+" sequence unique), " +
				pdb.getNumNonPolyChains()+" non-polymer chains.");

		for (String repChain:pdb.getAllRepChains()) {
			System.out.println(pdb.getSeqIdenticalGroupString(repChain));
		}
		System.out.println("Polymer chains: ");
		for (PdbChain chain:pdb.getPolyChains()) {
			System.out.println(chain.getChainCode()+"("+chain.getPdbChainCode()+")");
		}
		System.out.println("Non-polymer chains: ");
		for (PdbChain chain:pdb.getNonPolyChains()) {
			System.out.println(chain.getChainCode()+"("+chain.getPdbChainCode()+") "+" residues: "+chain.getObsLength()+
					" ("+chain.getFirstResidue().getLongCode()+"-"+chain.getFirstResidue().getSerial()+")");
		}
		
		System.out.println(pdb.getSpaceGroup().getShortSymbol()+" ("+pdb.getSpaceGroup().getId()+")");
		if (debug) System.out.println("Symmetry operators: "+pdb.getSpaceGroup().getNumOperators());
		
		System.out.println("Calculating possible interfaces... (using "+nThreads+" CPUs for ASA calculation)");
		long start = System.currentTimeMillis();
		ChainInterfaceList interfaces = pdb.getAllInterfaces(CUTOFF, null, Asa.DEFAULT_N_SPHERE_POINTS, nThreads, true, nonPoly, debug);
		long end = System.currentTimeMillis();
		System.out.println("Total time for interface calculation: "+(end-start)/1000+"s");
		
		System.out.println("Total number of interfaces found: "+interfaces.size());

		PymolRunner pr = new PymolRunner(PYMOL_EXE);
					
		for (int i=0;i<interfaces.size();i++) {
			ChainInterface interf = interfaces.get(i+1);
			interf.calcRimAndCore(BSATOASA_CUTOFF);
			String parallel = "";
			if (interf.isParallel()) parallel = " -- PARALLEL interface";
			System.out.println("\n##Interface "+(i+1)+" "+
					interf.getFirstSubunitId()+"-"+
					interf.getSecondSubunitId()+parallel);
			if (interf.hasClashes()) System.out.println("CLASHES!!!");
			System.out.println("Transf1: "+SpaceGroup.getAlgebraicFromMatrix(interf.getFirstTransf())+
					". Transf2: "+SpaceGroup.getAlgebraicFromMatrix(interf.getSecondTransf()));
			System.out.println(interf.getFirstMolecule().getChainCode()+" - "+interf.getSecondMolecule().getChainCode());
			System.out.printf("Connection vector: (%5.2f %5.2f %5.2f)\n",interf.getConnectionVector().x,interf.getConnectionVector().y,interf.getConnectionVector().z);
			System.out.println("Number of contacts: "+interf.getNumContacts());
			System.out.println("Number of contacting atoms (from both molecules): "+interf.getNumAtomsInContact());
			System.out.println("Number of core residues at "+String.format("%4.2f", BSATOASA_CUTOFF)+
					" bsa to asa cutoff: "+interf.getFirstRimCore().getCoreSize()+" "+interf.getSecondRimCore().getCoreSize());
			System.out.printf("Interface area: %8.2f\n",interf.getInterfaceArea());
		
			if (writeDir!=null) {
				File pdbFile = new File(writeDir,outBaseName+"."+(i+1)+".interface.pdb");
				interf.writeToPdbFile(pdbFile);
				if (generatePngs) {
					pr.generateInterfPngPsePml(interf, pdbFile, 
							new File(writeDir,outBaseName+"."+(i+1)+".pse"), 
							new File(writeDir,outBaseName+"."+(i+1)+".pml"), outBaseName);
				}
			}
		}
		
		if (serialize && writeDir!=null) {
			Goodies.serialize(new File(writeDir,outBaseName+".interfaces.dat"), interfaces);
		}
		
	}

}
