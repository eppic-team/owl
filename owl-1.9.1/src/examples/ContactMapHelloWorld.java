package examples;

import java.io.IOException;

import edu.uci.ics.jung.graph.util.Pair;

import owl.core.structure.graphs.FileRIGraph;
import owl.core.structure.graphs.RIGEdge;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.FileFormatError;

public class ContactMapHelloWorld {
	
	/**
	 * A simple example script that takes to contact maps and
	 * reports on all contacts in a if they exist in b.
	 * @param args requires path and filenames of two contact maps as command line parameters
	 * @throws IOException
	 * @throws FileFormatError
	 * @author matthias winkelmann
	 */
	public static void main(String[] args) throws IOException, FileFormatError {
		
		// loading the contact maps into graph objects
		
		RIGraph graph, graph2;
		graph = new FileRIGraph(args[0]);
		graph2 = new FileRIGraph(args[1]);
		
		// we get all edges (contacts) from the first graph and loop through them
		
		for (RIGEdge edge:graph.getEdges()) {
			
			// the edge doesn't actually contain information about the residues it connects.
			// we get that information from the graph
			
			Pair<RIGNode> endpoints = graph.getEndpoints(edge);
			int i = endpoints.getFirst().getResidueSerial();
			int j = endpoints.getSecond().getResidueSerial();
			
			// check if the edge between residues i and j also exists in the second graph 
			// and write results to output
			
			if (graph2.containsEdgeIJ(i,j)) {
				System.out.println("Edge "+i+" "+j+" is in both graphs");
			} else {
				System.out.println("Edge "+i+" "+j+" is in first but not second graph");
			}
			
		}
		
		
			
		
		
		
		
	}
}
