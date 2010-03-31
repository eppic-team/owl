import java.sql.SQLException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import owl.core.structure.*;


import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;



public class testJUNGframework {

	public static void main(String[] args) throws FileFormatError, PdbLoadError, SQLException, PdbCodeNotFoundError {
		Pdb pdb = new PdbasePdb("7adh");
		pdb.load("A");
		
		RIGraph graph = pdb.getRIGraph("Ca", 8.0);
		RIGraph graph2 = pdb.getRIGraph("Ca/Cg", 8.0);
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
		HashMap<Pair<Integer>,Integer> cmNbhSizes = graph.getAllCommonNbhSizes();

		for (Pair<Integer> pair:cmNbhSizes.keySet()) {
			//for (int j=170;j<=180;j++) {
			if (pair.getFirst()==109 && pair.getSecond()>=170 && pair.getSecond()<=180) {
				RIGNbhood nbh = graph.getNbhood(graph.getNodeFromSerial(pair.getSecond()));
				System.out.println("nbhood of "+pair.getSecond()+": "+ nbh);
				for (RIGNode nb:nbh.values()) {
					System.out.print(" "+nb);
				}
				System.out.println();

				int size = cmNbhSizes.get(pair);
				System.out.println("common neighbs for edge "+pair.getFirst()+", "+pair.getSecond()+": "+size);
				RIGCommonNbhood cmNbh = graph.getCommonNbhood(graph.getNodeFromSerial(pair.getFirst()), graph.getNodeFromSerial(pair.getSecond()));
				for (RIGNode nb:cmNbh.values()){
					System.out.print(" "+nb);
				}
				System.out.println();
				System.out.println(cmNbh);
			}
		}

		// TEST of evaluate prediction
		System.out.println("TEST eval prediction");
		PredEval pev = graph.evaluatePrediction(graph2);
		pev.print();
		
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
		
		// TEST of degree methods and getSuccessorCount, getPredecessorCount, getNeighborCount
		// degree and neighbor counts are different when the graph has parallel edges (our graphs shouldn't have them anyway)
		System.out.println("TEST degree and neighbor count methods");
		// graph is undirected, graph2 is directed
		System.out.println("undirected graph: ");
		for (int serial=1;serial<=10;serial++){
			RIGNode node = graph.getNodeFromSerial(serial);
			System.out.println("node "+serial+": nbr count "+graph.getNeighborCount(node)+", predec count: "+graph.getPredecessorCount(node)+", successor count: "+graph.getSuccessorCount(node));
			System.out.println("node "+serial+": degree "+graph.degree(node)+", in degree: "+graph.inDegree(node)+", out degree: "+graph.outDegree(node));
		}
		System.out.println("directed graph: ");
		for (int serial=1;serial<=10;serial++){
			RIGNode node = graph2.getNodeFromSerial(serial);
			System.out.println("node "+serial+": nbr count "+graph2.getNeighborCount(node)+", predec count: "+graph2.getPredecessorCount(node)+", successor count: "+graph2.getSuccessorCount(node));
			System.out.println("node "+serial+": degree "+graph2.degree(node)+", in degree: "+graph2.inDegree(node)+", out degree: "+graph2.outDegree(node));
		}

		// TEST adding and getting directed/undirected edges
		System.out.println("TEST adding and getting directed/undirected edges: ");
		RIGraph g = new RIGraph("AMCD");
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
	
		if (g.containsEdgeIJ(3, 1)){
			System.out.println("Edge 3,1 found");
		}
		if (g.containsEdgeIJ(1, 3)){
			System.out.println("Edge 1,3 found");
		}
		
		
		// TEST copying of RIGraph objects and restrictRange methods
		System.out.println("TEST copying of RIGraph objects");
		Pdb pdb1 = new PdbasePdb("1bxy");
		pdb1.load("A");
		Pdb pdb2 = new PdbasePdb("1sha");
		pdb2.load("A");
		RIGraph g1 = pdb1.getRIGraph("Ca", 8);
		RIGraph g2 = pdb2.getRIGraph("Ca", 8);
		System.out.println("g1: "+g1.getPdbCode()+" nodes: "+g1.getVertexCount()+" edges: "+g1.getEdgeCount());
		System.out.println("g2: "+g2.getPdbCode()+" nodes: "+g2.getVertexCount()+" edges: "+g2.getEdgeCount());
		RIGraph newgraph = g1.copy();
		System.out.println("new g: "+newgraph.getPdbCode()+" nodes: "+newgraph.getVertexCount()+" edges: "+newgraph.getEdgeCount());
		newgraph = g2;
		System.out.println("new g reassigned to g2: "+newgraph.getPdbCode()+" nodes: "+newgraph.getVertexCount()+" edges: "+newgraph.getEdgeCount());
		newgraph.restrictContactsToMaxRange(10);
		System.out.println("new g reassigned to g2, max range 10: "+newgraph.getPdbCode()+" nodes: "+newgraph.getVertexCount()+" edges: "+newgraph.getEdgeCount());
		System.out.println("g1: "+g1.getPdbCode()+" nodes: "+g1.getVertexCount()+" edges: "+g1.getEdgeCount());
		
	}
	

}
