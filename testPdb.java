import java.io.IOException;
import java.sql.SQLException;
//import java.util.ArrayList;

import proteinstructure.*;

//import java.util.TreeMap;

public class testPdb {

	/**
	 * Test class for classes in aglappe.proteinstructure package (Pdb, Graph, AA, Contact) 
	 * Loads structure data from pdbase, msdsd and pdb file and calculates graph dumping contacts to file
	 * Also gets contact map object and some common neighbours from it
	 * @throws SQLException 
	 */
	
	public static void main(String[] args) throws IOException, PdbaseInconsistencyError, PdbaseAcCodeNotFoundError, MsdsdAcCodeNotFoundError, MsdsdInconsistentResidueNumbersError, SQLException{
		String accode="1bxy";
		String chaincode="A";
		// data from pdbase
		System.out.println("loading structure from pdbase");
		Pdb pdb = new PdbasePdb(accode,chaincode); // by default it goes to database "pdbase", a name for another pdbase database + a MySQLConnection can also be passed
		System.out.println("dumping structure to pdb file");
		pdb.dump2pdbfile("testpdb.txt");
		String chain = pdb.getChainCode();
		System.out.println("getting graph");
		Graph graph = pdb.get_graph("ALL", 4.1);
		System.out.println("dumping contacts to file");
		graph.write_contacts_to_file("test.txt");
		System.out.println("getting start of contact map matrix from graph");
		int[][] mat = graph.getIntMatrix();
		for (int i=0;i<10;i++){
			for (int j=0;j<10;j++){
				System.out.print(mat[i][j]);
			}
			System.out.println();
		}

		// data from msdsd
		System.out.println("loading structure from msdsd");
		Pdb pdb2 = new MsdsdPdb(accode,chaincode); // by default it goes to database "msdsd_00_07_a", a name for another msdsd database + a MySQLConnection can also be passed 
		System.out.println("dumping structure to pdb file");
		pdb2.dump2pdbfile("testpdb2.txt");
		System.out.println("getting graph");
		Graph graph2 = pdb2.get_graph("ALL", 4.1);
		System.out.println("dumping contacts to file");
		graph2.write_contacts_to_file("test2.txt");

		// data from pdb
		System.out.println("reading from dumped pdb file");
		Pdb pdb3 = new PdbfilePdb("testpdb.txt",chain); // we read from the file dumped from pdbase, note that dumped file contains internal chain identifier
		System.out.println("getting graph");
		Graph graph3 = pdb3.get_graph("ALL", 4.1);
		System.out.println("dumping contacts to file");
		graph3.write_contacts_to_file("test3.txt");
	}

}
