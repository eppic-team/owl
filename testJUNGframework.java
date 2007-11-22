import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;


import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

import proteinstructure.*;


public class testJUNGframework {

	public static void main(String[] args) throws FileNotFoundException, IOException, GraphFileFormatError, PdbaseInconsistencyError, PdbCodeNotFoundError, SQLException, PdbChainCodeNotFoundError, GraphIdNotFoundError {
		Pdb pdb = new PdbasePdb("7adh", "A");
		
		RIGraph graph = pdb.get_graph("Ca", 8.0);
		RIGraph graph2 = pdb.get_graph("Ca/Cb", 8.0);
		//RIGraph graph = new FileRIGraph("/scratch/local/temp.graph");
		//RIGraph graph = new DbRIGraph("casp_decoys",9);
		//System.out.println(graph);
		System.out.println("#pdb code:"+graph.getPdbCode());
		System.out.println("#chain code:"+graph.getChainCode());
		System.out.println("#pdb chain code:"+graph.getPdbChainCode());
		System.out.println("#full length:"+graph.getFullLength());
		System.out.println("#obs length:"+graph.getObsLength());
		System.out.println("#ct:"+graph.getContactType());
		System.out.println("#cutoff:"+graph.getCutoff());
		System.out.println("#edge count:"+graph.getEdgeCount());
		System.out.println("#vertex count:"+graph.getVertexCount()); 

		// TEST getting atom weights and weights
		for (RIGEdge e:graph.getEdges()) {
			Pair<RIGNode> pair = graph.getEndpoints(e);
			System.out.println(pair.getFirst()+" "+pair.getSecond()+" atom w: "+e.getAtomWeight()+" w: "+e.getWeight() );
		}
		
		
		// TEST printing node information
		System.out.println("TEST printing node information");
		for (int ind:graph.getSerials()){
			RIGNode n = graph.getNodeFromSerial(ind);
			SecStrucElement sselem = n.getSecStrucElement();
			char sstype = sselem==null?0:sselem.getType();
			System.out.println("V"+n.getResidueSerial()+" "+n.getResidueType()+" "+sstype);
			
		}		
		System.out.println();
		
		
		// TEST of neighbourhood methods
		System.out.println("TEST neighbor methods");
		HashMap<Pair<RIGNode>,Integer> cmNbhSizes = graph.getAllCommonNbhSizes();
		int i = 109;
		for (int j=170;j<=180;j++) {
			RIGNbhood nbh = graph.getNbhood(graph.getNodeFromSerial(j));
			System.out.println("nbhood of "+j+": "+ nbh);
			for (RIGNode nb:nbh.values()) {
				System.out.print(" "+nb);
			}
			System.out.println();
			int size = cmNbhSizes.get(new Pair<RIGNode>(graph.getNodeFromSerial(i),graph.getNodeFromSerial(j)));
			System.out.println("common neighbs for edge "+i+", "+j+": "+size);
			RIGCommonNbhood cmNbh = graph.getCommonNbhood(graph.getNodeFromSerial(i), graph.getNodeFromSerial(j));
			for (RIGNode nb:cmNbh.values()){
				System.out.print(" "+nb);
			}
			System.out.println();
			System.out.println(cmNbh);
		}

		// TEST of evaluate prediction
		System.out.println("TEST eval prediction");
		PredEval pev = graph.evaluatePrediction(graph2);
		pev.print();

		
		// TEST of range restrict methods
		System.out.println();
		System.out.println("TEST restrict methods");
		System.out.println("num contacts full range: "+graph.getEdgeCount());
		graph.restrictContactsToMaxRange(10);
		System.out.println("num contacts max range 10: "+graph.getEdgeCount());
		
		
		
		// TEST of isDirected()
		System.out.println("TEST isDirected");
		if (graph.isDirected()) {
			System.out.println("directed");
		}
		
		// TEST of Pair<Integer>, does it respect equals
		System.out.println("TEST Pair<Integer>");
		Pair<Integer> p1 = new Pair<Integer>(1,2);
		Pair<Integer> p2 = new Pair<Integer>(1,2);
		Set<Pair<Integer>> set = new HashSet<Pair<Integer>>();
		set.add(p1);
		set.add(p2);
		if (p1.equals(p2)) {
			System.out.println("2 pairs are equal");
		}
		System.out.println(set);
		
		// TEST of getSuccessorCount, getPredecessorCount, getNeighborCount
		System.out.println("TEST degree methods: getNeighborCount, getSuccessorCount, getPredecessorCount");
		// graph is undirected, graph2 is directed
		System.out.println("undirected graph: ");
		for (int serial=1;serial<=10;serial++){
			RIGNode node = graph.getNodeFromSerial(serial);
			System.out.println("degree for node "+serial+": "+graph.getNeighborCount(node)+", indegree: "+graph.getPredecessorCount(node)+", outdegree: "+graph.getSuccessorCount(node));
		}
		System.out.println("directed graph: ");
		for (int serial=1;serial<=10;serial++){
			RIGNode node = graph2.getNodeFromSerial(serial);
			System.out.println("degree for node "+serial+": "+graph2.getNeighborCount(node)+", indegree: "+graph2.getPredecessorCount(node)+", outdegree: "+graph2.getSuccessorCount(node));
		}

		// TEST adding and getting undirected edges
		System.out.println("TEST adding and getting undirected edges: ");
		RIGraph g = new RIGraph("ABCD");
		EdgeType edgeType = EdgeType.UNDIRECTED;
		RIGNode n1 = g.getNodeFromSerial(1);
		RIGNode n2 = g.getNodeFromSerial(2);
		RIGNode n3 = g.getNodeFromSerial(3);
		g.addEdge(new RIGEdge(1), n1, n2, edgeType);
		g.addEdge(new RIGEdge(2), n2, n1, edgeType);
		g.addEdge(new RIGEdge(3), n1, n2, edgeType);
		g.addEdge(new RIGEdge(4), n1, n3, edgeType);
		for (RIGEdge e:g.getEdges()) {
			System.out.println(g.getEndpoints(e).getFirst()+" "+g.getEndpoints(e).getSecond()+", w "+e.getAtomWeight());
		}
		System.out.println("vertex count: "+g.getVertexCount());
		System.out.println("edge count: "+g.getEdgeCount());
		System.out.println("edge w "+ g.findEdge(n3, n1).getAtomWeight());
	}
	

}
