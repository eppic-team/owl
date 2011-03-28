package examples;
import java.io.IOException;
import java.sql.SQLException;

import edu.uci.ics.jung.graph.util.Pair;
import owl.core.structure.Pdb;
import owl.core.structure.PdbCodeNotFoundException;
import owl.core.structure.PdbLoadException;
import owl.core.structure.PdbfilePdb;
import owl.core.structure.graphs.*;





public class HelloAGLappe {
	/* hello World program for using AG Lappe PDB* classes 

* load PDF file (given as first argument to main / chain args1) to Object 
* generate random subset of size x% 
* extract nbhoodstrings and score them 
* make a copy and simulate edge movements 	

	/**
	 * @param args
	 */
	public static void main(String[] args) throws PdbLoadException, IOException, SQLException, PdbCodeNotFoundException {
		String cType = "Ca"; // contact type like "Ca"/"Cb"/"ALL"
		double cutoff = 8.0; 
		String myCode = args[0]; 
		String myChain = args[1]; 
		String outfile = myCode+"_"+myChain+"_"+cType+"_"+cutoff+".cm"; 
		// under UNIX : String infile = "/home/lappe/"+myCode+".pdb";
		// on Mac 
		String infile = "/Users/lappe/_PDBs/"+myCode+".pdb";
		String nbString=""; 
		RIGNbhood myNbs;
		RIGNode iNode, jNode, myNode; 
		
		// alternative A : load from file 
		System.out.println("Loading PDB "+myCode+" chain "+myChain+" from file "+infile+"."); 
		Pdb original = new PdbfilePdb( infile);
		
		// alternative B : load from database  
		// System.out.println("Loading PDB "+myCode+" chain "+myChain+" from database."); 
		// Pdb original = new PdbasePdb( myCode); 
		
		System.out.println("original Pdb object created.");
		original.load( myChain); 
		System.out.println("chain loaded.");
		RIGraph origraph = original.getRIGraph(cType, cutoff); 
		System.out.println("Converted to RIG ("+cType+","+cutoff+")"); 
		// for test write out the contacts as a file 
		System.out.println("Writing contact map to "+outfile); 
		origraph.writeToFile( outfile);
		
		// example of looping through the collection of nodes, ordered by serial  
		System.out.println("\n\nExtracting nbstrings for each node (ordered):"); 
		for ( int i:origraph.getSerials()) {
			myNode=origraph.getNodeFromSerial( i); 
			System.out.print("\n> "+myNode.getResidueSerial()+":"+myNode.getResidueType());
			myNbs = origraph.getNbhood(myNode); 
			nbString = myNbs.getNbString(); 
			System.out.print("\t nbs: "+nbString);	
		} // next node 	
		
		// example of looping through the collection of nodes, unordered from collection 
		System.out.println("\n\nExtracting nbstrings for each node (collection, unordered):"); 
		for (RIGNode cNode:origraph.getVertices()) {
			System.out.print("\n> "+cNode.getResidueSerial()+":"+cNode.getResidueType());
			myNbs = origraph.getNbhood(cNode); 
			nbString = myNbs.getNbString(); 
			System.out.print("\t nbs: "+nbString);	
		} // next node 
		
		// example of looping through collection of edges 
		System.out.println("\n\nExtracting edgelist :");
		for (RIGEdge myEdge:origraph.getEdges()) {
			Pair<RIGNode> twoNodes = origraph.getEndpoints(myEdge);
			iNode = twoNodes.getFirst(); 
			jNode = twoNodes.getSecond(); 
			System.out.println( iNode.getResidueType()+" <==> "+jNode.getResidueType()); 
		} // next edge 
		
		
	}
}

/* for (int i=origraph.getFirstResidueSerial(); i<=origraph.getLastResidueSerial(); i++) {
System.out.print("\n> "+i+"\t"+origraph.getSequence().substring(i-1,i)+"\t"+origraph.getNeighborhoodString(i)); 

} // next residue serial
*/ 
/*
String myfile = args[0]; 
		Pdb original = new PdbfilePdb( myfile);  
		original.load( args[1]); 
		RIGraph origraph = original.getRIGraph("Ca", 8.0); // contact type like "Ca"/"Cb"/"ALL"... cutoff double 
		// for test write out the contacts as a file 
		origraph.write_graph_to_file(myfile+"_Ca_8.cm"); 
*/ 