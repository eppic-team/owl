import java.io.IOException;
import java.sql.SQLException;

import edu.uci.ics.jung.graph.util.Pair;

import proteinstructure.*;

/**
 * An example program for some of the classes in the proteinstructure package (Pdb, RIGraph, ...) 
 * Loads structure data from pdbase and pdb file and calculates graph dumping contacts to file 
 * 
 */
public class testPdb {

	public static void main(String[] args) throws SQLException, PdbCodeNotFoundError, PdbLoadError, IOException {
		
		String pdbCode="1bxy";
		String pdbChainCode="A";
		
		// data from pdbase database
		System.out.println("loading structure from pdbase");
		Pdb pdbFromPdbase = new PdbasePdb(pdbCode); // by default it goes to database "pdbase", a name for another pdbase database + a MySQLConnection can also be passed
		pdbFromPdbase.load(pdbChainCode);
		System.out.println("dumping structure to pdb file");
		pdbFromPdbase.writeToPDBFile("test_dump_from_pdbase.pdb");
		// note that the chainCode is not necessarily the same as the pdbChainCode
		String chainCode = pdbFromPdbase.getChainCode();
		System.out.println("getting graph");
		RIGraph graph = pdbFromPdbase.getRIGraph("ALL", 4.1);
		System.out.println("writing graph to file");
		graph.write_graph_to_file("test_from_pdbase.cm");
		
		// getting edges and some information from the end nodes
		System.out.println("i\ti_res\tj\tj_res\tdistance");
		for (RIGEdge edge:graph.getEdges()) {
			Pair<RIGNode> pair = graph.getEndpoints(edge);
			RIGNode nodei = pair.getFirst();
			RIGNode nodej = pair.getSecond();
			System.out.println(nodei.getResidueSerial()+"\t"+nodei.getResidueType()+"\t"+nodej.getResidueSerial()+"\t"+nodej.getResidueType()+"\t"+edge.getDistance());
		}

		// data from pdb file
		System.out.println("reading from dumped pdb file");
		Pdb pdbFromFile = new PdbfilePdb("test_dump_from_pdbase.pdb"); // we read from the file dumped from pdbase, note that dumped file contains internal chain identifier
		pdbFromFile.load(chainCode);
		System.out.println("getting graph");
		RIGraph graph2 = pdbFromFile.getRIGraph("ALL", 4.1);
		System.out.println("writing graph to file");
		graph2.write_graph_to_file("test_pdb_from_file.cm");
	}

}
