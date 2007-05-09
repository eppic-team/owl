import java.io.IOException;

import proteinstructure.Graph;
import proteinstructure.Pdb;
import proteinstructure.PdbaseAcCodeNotFoundError;
import proteinstructure.PdbaseInconsistencyError;
public class testPdb {

	/**
	 * Test class for classes in aglappe.proteinstructure package (Pdb, Graph, AA, Contact) 
	 * Loads structure data from db and calculates graph dumping contacts to file
	 */
	
	public static void main(String[] args) throws IOException, PdbaseInconsistencyError, PdbaseAcCodeNotFoundError{
		System.out.println("loading structure");
		Pdb pdb = new Pdb("1bxy","A");
		System.out.println("getting graph");
		Graph graph = pdb.get_graph("Ca", 8.0);
		System.out.println("dumping contacts to file");
		graph.write_contacts_to_file("test.txt");
		

	}

}
