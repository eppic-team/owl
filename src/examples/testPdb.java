package examples;
import java.io.File;
import java.io.IOException;

import owl.core.structure.*;
import owl.core.structure.graphs.RIGEdge;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.FileFormatException;

import edu.uci.ics.jung.graph.util.Pair;


/**
 * An example program for some of the classes in the owl.core.structure package (PdbChain, RIGraph, ...) 
 * Loads structure data from pdbase and pdb file and calculates graph, writing contacts to file 
 * 
 */
public class testPdb {

	public static void main(String[] args) throws PdbLoadException, IOException, FileFormatException {
		
		String pdbCode="1bxy";
		String pdbChainCode="A";
		
		// data from pdbase database
		System.out.println("loading structure from pdbase");
		File cifFile = new File(System.getProperty("java.io.tmpdir"),pdbCode+".cif");
		cifFile.deleteOnExit();
		PdbAsymUnit.grabCifFile("/path/to/mmCIF/gz/all/repo/dir", null, pdbCode, cifFile, false);				
		PdbAsymUnit fullpdb = new PdbAsymUnit(cifFile);
				
		PdbChain pdbFromPdbase = fullpdb.getChain(pdbChainCode);
		System.out.println("dumping structure to pdb file");
		pdbFromPdbase.writeToPDBFile(new File("test_dump_from_pdbase.pdb"));
		// note that the chainCode is not necessarily the same as the pdbChainCode
		String chainCode = pdbFromPdbase.getChainCode();
		System.out.println("getting graph");
		RIGraph graph = pdbFromPdbase.getRIGraph("ALL", 4.1);
		System.out.println("writing graph to file");
		graph.writeToFile("test_from_pdbase.cm");
		
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
		// we read from the file dumped from pdbase, note that dumped file contains internal (CIF) chain identifier
		PdbAsymUnit pdbFromFile = new PdbAsymUnit(new File("test_dump_from_pdbase.pdb"));
		PdbChain chain = pdbFromFile.getChain(chainCode);
		System.out.println("getting graph");
		RIGraph graph2 = chain.getRIGraph("ALL", 4.1);
		System.out.println("writing graph to file");
		graph2.writeToFile("test_pdb_from_file.cm");
	}

}
