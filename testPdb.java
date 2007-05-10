import java.io.IOException;
import java.util.ArrayList;

import proteinstructure.*;

import java.util.TreeMap;

public class testPdb {

	/**
	 * Test class for classes in aglappe.proteinstructure package (Pdb, Graph, AA, Contact) 
	 * Loads structure data from pdbase, msdsd and pdb file and calculates graph dumping contacts to file
	 * Also gets contact map object and some common neighbours from it
	 */
	
	public static void main(String[] args) throws IOException, PdbaseInconsistencyError, PdbaseAcCodeNotFoundError, MsdsdAcCodeNotFoundError, MsdsdInconsistentResidueNumbersError{
		
		// data from pdbase
		System.out.println("loading structure from pdbase");
		Pdb pdb = new Pdb("1bxy","A"); // pdbase is default source for constructor, a name for a pdbase database can also be passed
		System.out.println("dumping structure to pdb file");
		pdb.dump2pdbfile("testpdb.txt");
		System.out.println("getting graph");
		Graph graph = pdb.get_graph("ALL", 4.1);
		System.out.println("dumping contacts to file");
		graph.write_contacts_to_file("test.txt");
		System.out.println("getting contact map object from graph");
		ContactMap cm = graph.getCM();
		ArrayList<Contact> testconts = new ArrayList<Contact>();
		for (int i=1;i<3;i++){
			for (int j=1;j<11;j++){
				testconts.add(new Contact(i,j));
			}
		}
		System.out.println("printing some common neighbours");
		for (Contact con:testconts){
			int i=con.i;
			int j=con.j;
			if (i!=j){
				TreeMap<Integer,String> cns =cm.getComNbs(i,j);
				System.out.println("\n"+i+","+j);
				for (int k:cns.keySet()){
					System.out.print(k+": "+cns.get(k)+", ");
				}
			}
		}

		// data from msdsd
		System.out.println("loading structure from msdsd");
		Pdb pdb2 = new Pdb("1bxy","A","msdsd_00_07_a");
		System.out.println("dumping structure to pdb file");
		pdb2.dump2pdbfile("testpdb2.txt");
		System.out.println("getting graph");
		Graph graph2 = pdb2.get_graph("ALL", 4.1);
		System.out.println("dumping contacts to file");
		graph2.write_contacts_to_file("test2.txt");

		// data from pdb
		System.out.println("reading from dumped pdb file");
		Pdb pdb3 = new Pdb("testpdb.txt"); // we read from the file dumped from pdbase
		System.out.println("getting graph");
		Graph graph3 = pdb3.get_graph("ALL", 4.1);
		System.out.println("dumping contacts to file");
		graph3.write_contacts_to_file("test3.txt");
	}

}
