import proteinstructure.*;
import java.util.HashMap;


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
		Pdb pdb1 = new PdbasePdb(pdbcode1);
		pdb1.load(chaincode1);
		Pdb pdb2 = new PdbasePdb(pdbcode2);
		pdb2.load(chaincode2);

		System.out.println("getting graphs");
		RIGraph graph1 = pdb1.getRIGraph("ALL", 4.2);
		RIGraph graph2 = pdb2.getRIGraph("ALL", 4.2);

		HashMap<String,RIGraph> comparison = graph1.compare(graph2);
		
		RIGraph onlyIn1 = comparison.get("onlythis");
		RIGraph onlyIn2 = comparison.get("onlyother");
		RIGraph common = comparison.get("common");
		
		System.out.println("writing output files "+onlyIn1File+", "+onlyIn2File+", "+commonFile);
		onlyIn1.write_graph_to_file(onlyIn1File);
		onlyIn2.write_graph_to_file(onlyIn2File);
		common.write_graph_to_file(commonFile);
		
		

	}

}
