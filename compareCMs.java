import java.io.IOException;

import proteinstructure.*;
import java.util.HashMap;


public class compareCMs {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception, IOException, PdbaseInconsistencyError, PdbaseAcCodeNotFoundError, MsdsdAcCodeNotFoundError, MsdsdInconsistentResidueNumbersError{
		
		String onlyIn1File="onlyin1.graph";
		String onlyIn2File="onlyin2.graph";
		String commonFile="common.graph";
		
		
		String pdbcode1="1ses";
		String chaincode1="A";
		String pdbcode2="1set";
		String chaincode2="A";

		
		System.out.println("loading structures from pdbase");
		Pdb pdb1 = new PdbasePdb(pdbcode1,chaincode1); 
		Pdb pdb2 = new PdbasePdb(pdbcode2,chaincode2);

		System.out.println("getting graphs");
		Graph graph1 = pdb1.get_graph("ALL", 4.2);
		Graph graph2 = pdb2.get_graph("ALL", 4.2);

		HashMap<String,Graph> comparison = graph1.compare(graph2);
		
		Graph onlyIn1 = comparison.get("onlythis");
		Graph onlyIn2 = comparison.get("onlyother");
		Graph common = comparison.get("common");
		
		System.out.println("writing output files "+onlyIn1File+", "+onlyIn2File+", "+commonFile);
		onlyIn1.write_graph_to_file(onlyIn1File);
		onlyIn2.write_graph_to_file(onlyIn2File);
		common.write_graph_to_file(commonFile);
		
		

	}

}
