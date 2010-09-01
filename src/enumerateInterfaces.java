import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
//import java.util.Set;

import owl.core.runners.NaccessRunner;
import owl.core.structure.ChainInterface;
import owl.core.structure.Pdb;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.SpaceGroup;
import owl.core.util.MySQLConnection;


public class enumerateInterfaces {

	private static final File NACCESS_EXE = new File("/home/duarte_j/bin/naccess");

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

		if (args.length!=1) {
			System.err.println("Usage: enumerateInterfaces <pdb code>");
			System.exit(1);
		}
		String pdbCode = args[0];
		
		PdbAsymUnit pdb = new PdbAsymUnit(pdbCode, new MySQLConnection(), "pdbase");

		System.out.println(pdb.getSpaceGroup().getShortSymbol()+" ("+pdb.getSpaceGroup().getId()+")");
		
		System.out.println("Calculating possible interfaces...");
		long start = System.currentTimeMillis();
		Set<ChainInterface> interfaces = pdb.getAllInterfaces(CUTOFF);
		long end = System.currentTimeMillis();
		System.out.println("Done. Time "+(end-start)/1000+"s");

		
		System.out.println("Total number of interfaces found: "+interfaces.size());
		System.out.println("Calculating BSAs with NACCESS");
		
		for (Pdb chain:pdb.getAllChains()) {
			NaccessRunner nar = new NaccessRunner(NACCESS_EXE, "");
			nar.runNaccess(chain);
		}
		HashMap<String, HashMap<Integer,Double>> asas = pdb.getAbsSurfaceAccessibilities();
		
		for (ChainInterface interf:interfaces) {
			System.out.print(".");
			interf.calcBSAnaccess(NACCESS_EXE,asas);
		}
		System.out.println();
		
		List<ChainInterface> sortedInterfaces = new ArrayList<ChainInterface>();
		sortedInterfaces.addAll(interfaces);
		Collections.sort(sortedInterfaces);
		
		for (int i=sortedInterfaces.size()-1;i>=0;i--) {
			ChainInterface interf = sortedInterfaces.get(i);
			int j= sortedInterfaces.size()-i;
			System.out.println("\n##Interface "+j);
			System.out.println("Transf1: "+SpaceGroup.getAlgebraicFromMatrix(interf.getFirstTransf())+
					". Transf2: "+SpaceGroup.getAlgebraicFromMatrix(interf.getSecondTransf()));
			System.out.println(interf.getFirstMolecule().getPdbChainCode()+" - "+interf.getSecondMolecule().getPdbChainCode());
			System.out.println("Number of contacts: "+interf.getNumContacts());
			System.out.println("Number of contacting atoms (from both molecules): "+interf.getNumAtomsInContact());
			System.out.printf("Interface area: %8.2f (%8.2f)\n",interf.getInterfaceArea(),interf.getInterfaceArea()/2.0);
			
		}
	}

}
