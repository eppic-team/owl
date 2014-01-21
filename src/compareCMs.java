import java.io.File;
import java.util.HashMap;

import owl.core.structure.*;
import owl.core.structure.graphs.RIGraph;


public class compareCMs {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		
		String onlyIn1File="onlyin1.graph";
		String onlyIn2File="onlyin2.graph";
		String commonFile="common.graph";
		
		
		String pdbcode1="1ses";
		String chaincode1="A";
		String pdbcode2="1set";
		String chaincode2="A";

		
		System.out.println("loading structures from pdbase");
		PdbAsymUnit fullpdb1 = new PdbAsymUnit(new File(pdbcode1+".cif"));
		PdbChain pdb1 = fullpdb1.getChain(chaincode1);
		PdbAsymUnit fullpdb2 = new PdbAsymUnit(new File(pdbcode2+".cif"));
		PdbChain pdb2 = fullpdb2.getChain(chaincode2);

		System.out.println("getting graphs");
		RIGraph graph1 = pdb1.getRIGraph("ALL", 4.2);
		RIGraph graph2 = pdb2.getRIGraph("ALL", 4.2);

		HashMap<String,RIGraph> comparison = graph1.compare(graph2);
		
		RIGraph onlyIn1 = comparison.get("onlythis");
		RIGraph onlyIn2 = comparison.get("onlyother");
		RIGraph common = comparison.get("common");
		
		System.out.println("writing output files "+onlyIn1File+", "+onlyIn2File+", "+commonFile);
		onlyIn1.writeToFile(onlyIn1File);
		onlyIn2.writeToFile(onlyIn2File);
		common.writeToFile(commonFile);
		
		

	}

}
