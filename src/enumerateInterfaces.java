import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
//import java.util.Set;

import owl.core.structure.ChainInterface;
import owl.core.structure.PdbAsymUnit;
import owl.core.util.MySQLConnection;


public class enumerateInterfaces {

	private static final File naccessExecutable = new File("/home/duarte_j/bin/naccess");
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
		
		
		//Set<ChainInterface> interfaces = pdb.getAllInterfaces(5.0);
		List<ChainInterface> interfaces = pdb.getAllInterfaces(5.0);
		System.out.println("Total number of interfaces found: "+interfaces.size());
		System.out.println("Calculating BSAs with NACCESS");
		for (ChainInterface interf:interfaces) {
			System.out.print(".");
			interf.calcBSAnaccess(naccessExecutable);
		}
		System.out.println();
		
		List<ChainInterface> sortedInterfaces = new ArrayList<ChainInterface>();
		sortedInterfaces.addAll(interfaces);
		Collections.sort(sortedInterfaces);
		
		for (int i=sortedInterfaces.size()-1;i>=0;i--) {
			ChainInterface interf = sortedInterfaces.get(i);
			int j= sortedInterfaces.size()-i;
			System.out.println("\n##Interface "+j);
			System.out.println("Transf1: "+interf.getFirstTransf()+". Transf2: "+interf.getSecondTransf());
			System.out.println(interf.getFirstMolecule().getPdbChainCode()+" - "+interf.getSecondMolecule().getPdbChainCode());
			System.out.println("Contacting atoms: "+interf.getAICGraph().getEdgeCount());
			System.out.printf("Interface area: %8.2f\n",interf.getInterfaceArea());
			
		}
	}

}
