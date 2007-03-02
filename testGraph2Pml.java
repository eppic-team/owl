import java.sql.*;
import java.io.*;
import tools.*;

/**
 * Package:		tools
 * Class: 		testGraph2Pml
 * Author:		Ioannis Filippis, filippis@molgen.mpg.de
 * Date:		21/03/2006
 * 
 * Simple test class for visualizing contact graphs using Graph2Pml class plus
 * PyMol class. Run PyMol with the -R option to run the server.
 * 
 * Notes:
 * - You can find the results from all examples here in 
 *   /project/StruPPi/ioannis/workspace/aglappe/vis_examples/
 * 
 * Changelog:
 * 28/03/06 modified by IF (more examples added)
 * 27/03/06 modified by IF (functional with comments - added cgo examples)
 * 21/03/06 first created by IF (non functional)
 */

public class testGraph2Pml {

	public static void main(String[] args) {
		
		String connFile = "/project/StruPPi/ioannis/cnfs/msdgraph.my.cnf", pdbFileName = "", molObjName = "";
		Graph2Pml graphPml = null;
		PyMol pml = null;
		
		mySQLConnect SQLC =  new mySQLConnect();
		SQLC.readConnectionFile(connFile);
		Connection conn = SQLC.openConnection();
	
		//Create a new output stream to the pymol server running at the localhost
		//all commands will be send there
		PrintWriter serverOutPw = new PrintWriter(new PymolServerOutputStream("http://"+Machine.getClient()+":9123"), true);
		pml = new PyMol(serverOutPw);

		pml.openLog("test", "/home/filippis/Desktop");
		
		//Edges as lines
		//export the pdb file directly from msdsd and load it in pymol
		pdbFileName =  "1rx4_20717_52567";
		try {
			Msdsd2Pdb.export2File("1rx4", 20717, 52567, "/project/StruPPi/ioannis/tmp/"+pdbFileName, "filippis");
		} catch (Exception e) {
			e.printStackTrace();
		}
		molObjName = pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		System.out.println(molObjName);
		
		//create a pymol contact graph object for the graph 33729
		//use only SC_SC edges with SC_SC edge weight and nodes that have SC_SC edges
		//no further filtering in the edges based on the contact range is applied
		//no user based graph filtering is applied
		//the graph comes from msdsd
		//normal edges will be utilised
		//the graph will not be visualised as "directed" 
		graphPml = new Graph2Pml(serverOutPw, molObjName, 33729, "SC_SC", "SC_SC_in+SC_SC_out", "SC_SC", "true", "true", true, false, false, conn);
		//draw nodes, edges, special residues but not surface
		//the node size will be calculated based on the SC_SC_out field on asc order
		//the node color will be calculated based on the BB_BB_out field not discretised
		//all ALA nodes will be transparent
		//all residues that do not have any SC interactions will be the special ones colored green
		//all edges will be lines of width determined by the SC_SC field in asc order
		//all edges lines will have color calculated based on the BB_BB_out field not discretised
		//all edges lines with adjacent ALA's will have gaps
		graphPml.draw(true, true, true, false);		
		graphPml.setNodeSizeMethod("SC_SC_out", false);
		graphPml.setNodeColorMethod("BB_BB_out", false);
		graphPml.setNodeTranspCondition("(res = \"ALA\")");
		graphPml.setSpecialRes("green", "SELECT cid, num FROM newmsdgraph.nodes WHERE (graph_id = 33729) AND (sc_sc_in+sc_sc_out = 0);");
		graphPml.setEdgeSizeMethod("SC_SC", false);
		graphPml.setEdgeColorMethod("BB_BB", false);
		graphPml.setEdgeGapCondition("(i_res = \"ALA\" OR j_res = \"ALA\")");
		//create the contact graph
		graphPml.outputGraph();
		//save the image not ray-traced
		pml.saveImage("test1", "/home/filippis/Desktop", false);
		//if you want later to read the view from a log file
		//pml.openLog("test1_log","/home/filippis/Desktop");
		//get the view
		pml.getView("test1_view");
		//pml.closeLog();

		// delete all and reload the structure
		pml.delete("all", false);
		pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		
		//to check reversibility in setNodeSizeMethod, setEdgeSizeMethod 
		//and discretised coloring in setNodeColorMethod, setEdgeColorMethod
		graphPml.setNodeSizeMethod("SC_SC_out", true);
		graphPml.setNodeColorMethod("BB_BB_out", true);
		graphPml.setEdgeSizeMethod("SC_SC", true);
		graphPml.setEdgeColorMethod("BB_BB", true);
		graphPml.outputGraph();
		//set view by reading a log file
		//pml.getFileView("test2_view", "/home/filippis/Desktop/test1_log.pml");
		//pml.setView("test2_view");
		//set view based on the view of the first structure so that you can compare
		//different graphs for the same structure
		pml.setView("test1_view");
		pml.saveImage("test2","/home/filippis/Desktop", false);
		
		// delete all and reload the structure
		pml.delete("all", false);
		pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		
		//just check user defined set/discretised colors
		//all green nodes must become red
		graphPml.setNodeColors(new String[] {"red", "magenta", "orange", "chocolate"});
		//switched red-green, yellow-blue
		graphPml.setEdgeColors(new String[] {"red", "green", "yellow", "blue", "violet", "cyan", "salmon", "lime", "pink", "slate", "magenta", "orange", "marine", "olive", "purple", "teal", "forest", "firebrick", "chocolate", "wheat", "white", "grey"});
	    graphPml.outputGraph();
		pml.setView("test1_view");
		pml.saveImage("test2a","/home/filippis/Desktop", false);
		
		pml.delete("all", false);
		pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		
		//to check node method for edge coloring when nodeColor method not discretised
		graphPml.setUniformNodeSize(0.6);
		graphPml.setNodeColorMethod("BB_BB_out", false);
		graphPml.setUniformEdgeSize(4);
		graphPml.setEdgeColorMethod("node", true);
		graphPml.outputGraph();
		pml.setView("test1_view");
		pml.saveImage("test3","/home/filippis/Desktop", false);
		
		pml.delete("all", false);
		pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
				
		//to check node method for edge coloring when nodeColor method is discretised
		graphPml.setUniformNodeSize(0.6);
		graphPml.setNodeColorMethod("BB_BB_out", true);
		graphPml.setUniformEdgeSize(4);
		graphPml.setEdgeColorMethod("node", true);
		graphPml.outputGraph();
		pml.setView("test1_view");
		pml.saveImage("test4","/home/filippis/Desktop", false);
		
		pml.delete("all", false);
		pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
				
		//to check node method for edge coloring when nodeColor method is uniform
		graphPml.setUniformNodeSize(0.6);
		graphPml.setUniformNodeColor("blue");
		graphPml.setUniformEdgeSize(4);
		graphPml.setEdgeColorMethod("node", true);
		graphPml.outputGraph();
		pml.setView("test1_view");
		pml.saveImage("test5","/home/filippis/Desktop", false);		
		
		pml.delete("all", false);
		pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		
		//to check the trick of directionality with gapped lines for "half" edges
		//for this purpose a SC_BB graph model is used that doesn't ensure bidirectionality
		graphPml = new Graph2Pml(serverOutPw, molObjName, 33729, "SC_BB", "SC_BB_in+SC_BB_out", "SC_BB", "true", "true", true, false, true, conn);
		graphPml.outputGraph();
		pml.setView("test1_view");
		pml.saveImage("test6","/home/filippis/Desktop", false);		
		
		
		pml.delete("all", false);
		
		//test multi-chain macromolecule with RNA also
		pdbFileName =  "1a0a_1107_2560";
		try {
			Msdsd2Pdb.export2File("1a0a", 1107, 2560, "/project/StruPPi/ioannis/tmp/"+pdbFileName, "filippis");
		} catch (Exception e) {
			e.printStackTrace();
		}
		molObjName = pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		System.out.println(molObjName);
		
		graphPml = new Graph2Pml(serverOutPw, molObjName, 95, "SC_SC", "SC_SC_in+SC_SC_out", "SC_SC", "true", "true", true, false, false, conn);
		graphPml.draw(true, true, false, false);
		//color nodes based on the chains they belong to
		graphPml.setChainNodeColor();
		graphPml.setNodeTransp(0.2);		
		graphPml.setUniformEdgeColor("red");
		graphPml.setUniformEdgeSize(4.5);
		graphPml.setEdgeGap(0.5);
		graphPml.outputGraph();
		pml.saveImage("test7","/home/filippis/Desktop", false);
		pml.getView("test2_view");
		
		pml.delete("all", false);
		pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		
		//color edges based on the chains they belong to
		graphPml.draw(true, true, false, false);		
		graphPml.setUniformNodeColor("red");
		graphPml.setChainEdgeColor();
		graphPml.outputGraph();
		pml.setView("test2_view");
		pml.saveImage("test8","/home/filippis/Desktop", false);		
		
		pml.delete("all", false);
		pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		
		//to check node method for edge coloring when nodeColor method is based on chain coloring
		//also check the user defined chain colors for nodes
		graphPml.setChainNodeColor(new String[] {"blue"});
		graphPml.setEdgeColorMethod("node", true);
		graphPml.outputGraph();
		pml.setView("test2_view");
		pml.saveImage("test9","/home/filippis/Desktop", false);
		
		pml.delete("all", false);
		pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		
		//draw also surface - you can not see the graph anymore
		graphPml.draw(true, true, false, true);
		graphPml.setSurfTransp(0.7);
		graphPml.outputGraph();
		pml.setView("test2_view");
		pml.saveImage("test10","/home/filippis/Desktop", false);
		
		try { Thread.sleep(10000); } catch (Exception e) {}
		
		//hide the surface, show the RNA also
		//change the color of the nodes for each chain using the pymol defined list of nodes for each chain 
		pml.hideWhat("surface", "graphMol", false);
		pml.showWhat("lines", "restMol", false);
		graphPml.setNodesColorUniform("nodes_cidA","red");
		graphPml.setNodesColorUniform("nodes_cidA","yellow");
		
		try { Thread.sleep(10000); } catch (Exception e) {}
		
		//if you believe that the objects contained in the list won't exceed the pymol command length limit,
		//you can concatenate all objects as a single string selection and handle them all together 
		pml.concatList("nodes", "+", "nodesStr");
		pml.setNodeColor("blue", "nodesStr", true);
		pml.selPml("nodes", "nodesStr", true, true);
				
		pml.delete("all", false);
		pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		
		//check a graph model with only inter-secondary-structure-elements SC_SC edges
		//remember that all the nodes with SC_SC will be selected
		//use default behaviour
		graphPml = new Graph2Pml(serverOutPw, molObjName, 95, "SC_SC", "SC_SC_in+SC_SC_out", "SC_SC", "(!((i_ssid = j_ssid) AND (i_sstype = j_sstype)))", "true", true, false, false, conn);
		graphPml.draw(true, true, false, false);
		graphPml.outputGraph();
		pml.setView("test2_view");
		pml.saveImage("test11","/home/filippis/Desktop", false);
		
		pml.delete("all", false);
		
		//CGO EDGES SECTION
		pdbFileName =  "1kj0_21698_55132";
		try {
			Msdsd2Pdb.export2File("1kj0", 21698, 55132, "/project/StruPPi/ioannis/tmp/"+pdbFileName, "filippis");
		} catch (Exception e) {
			e.printStackTrace();
		}
		molObjName = pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		
		//Cgo makes sense to be used with directed graph
		graphPml = new Graph2Pml(serverOutPw, molObjName, 35666, "SC_SC", "SC_SC_in+SC_SC_out", "SC_SC", "true", "true", true, true, true, conn);
		graphPml.draw(true, true, false, false);		
		graphPml.outputGraph();
		pml.saveImage("test12","/home/filippis/Desktop", false);
		pml.getView("test3_view");

		// delete all edges
		pml.iterateList("edges", "edge");
		pml.delete("edge", true);
		
		graphPml.draw(false, true, false, false);
		graphPml.setUniformCgoEdgeSize(0.2);
		graphPml.setEdgeColorMethod("node", true);
		graphPml.outputGraph();
		pml.setView("test3_view");
		pml.saveImage("test13","/home/filippis/Desktop", false);
		
		// delete all edges
		pml.iterateList("edges", "edge");
		pml.delete("edge", true);
				
		graphPml.setEdgeSizeMethod("SC_SC", false);
		graphPml.setEdgeColorMethod("BB_BB", false);
		graphPml.outputGraph();
		pml.setView("test3_view");
		pml.saveImage("test14","/home/filippis/Desktop", false);
		
		pml.iterateList("edges", "edge");
		pml.delete("edge", true);
		
		graphPml.setEdgeSizeMethod("SC_SC", true);
		graphPml.setEdgeColorMethod("BB_BB", true);
		graphPml.outputGraph();
		pml.setView("test3_view");
		pml.saveImage("test15","/home/filippis/Desktop", false);
				
		pml.delete("all", false);
		pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
				
		//to check node method for edge coloring when nodeColor method not discretised
		graphPml.draw(true, true, false, false);
		graphPml.setUniformNodeSize(0.6);
		graphPml.setNodeColorMethod("BB_BB_out", false);
		graphPml.setUniformEdgeSize(0.4);
		graphPml.setEdgeColorMethod("node", true);
		graphPml.outputGraph();
		pml.setView("test3_view");
		pml.saveImage("test16","/home/filippis/Desktop", false);
		
		pml.delete("all", false);
		pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		
		//to check node method for edge coloring when nodeColor method not discretised
		graphPml.setUniformNodeSize(0.6);
		graphPml.setNodeColorMethod("BB_BB_out", true);
		graphPml.setUniformEdgeSize(0.4);
		graphPml.setEdgeColorMethod("node", true);
		graphPml.outputGraph();
		pml.setView("test3_view");
		pml.saveImage("test17","/home/filippis/Desktop", false);
		
		pml.delete("all", false);
		pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		
		//see "real" half-edges
		graphPml = new Graph2Pml(serverOutPw, molObjName, 35666, "SC_BB", "SC_BB_in+SC_BB_out", "SC_BB", "true", "true", true, true, true, conn);
		graphPml.outputGraph();
		pml.setView("test3_view");
		pml.saveImage("test18","/home/filippis/Desktop", false);		
		
		pml.delete("all", false);
		
		//test multi-chain macromolecule with RNA also
		pdbFileName = "1a0a_1107_2560.pdb";
		molObjName = pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		System.out.println(molObjName);
		
		graphPml = new Graph2Pml(serverOutPw, molObjName, 95, "SC_SC", "SC_SC_in+SC_SC_out", "SC_SC", "true", "true", true, true, true, conn);
		graphPml.draw(true, true, false, false);
		graphPml.setNodeColorMethod("SC_SC_out", false);
		graphPml.setNodeTransp(0);		
		graphPml.setEdgeSizeMethod("SC_SC", false);
		graphPml.setChainEdgeColor(new String[] {"yellow","red"});
		graphPml.outputGraph();
		pml.saveImage("test19","/home/filippis/Desktop", false);
		
		pml.delete("all", false);
		
		pdbFileName = "1rx4_20717_52567.pdb";
		pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		
		graphPml = new Graph2Pml(serverOutPw, molObjName, 33729, "SC_SC", "SC_SC_in+SC_SC_out", "SC_SC", "true", "true", true, true, true, conn);
		graphPml.draw(true, true, false, false);
		graphPml.setNodeSizeMethod("SC_SC_out", false);
		graphPml.setNodeColorMethod("SC_SC_out", false);
		graphPml.setNodeTransp(0.6);		
		graphPml.setEdgeSizeMethod("SC_SC", false);
		graphPml.setEdgeColorMethod("SC_SC", false);
		graphPml.outputGraph();
		pml.saveImage("test20","/home/filippis/Desktop", false);
		
		pml.delete("all", false);
		pdbFileName = "1kj0_21698_55132.pdb";
		pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		
		//check the setNodeSizeRange method
		graphPml = new Graph2Pml(serverOutPw, molObjName, 35666, "SC_SC", "SC_SC_in+SC_SC_out", "SC_SC", "true", "true", true, true, true, conn);
		graphPml.draw(true, false, false, false);
		graphPml.setChainColors(new String[] {"red"});
		graphPml.setNodeSizeMethod("SC_SC_out", false);
		graphPml.setNodeSizeRange(new double[] {0.6, 1.2});
		graphPml.outputGraph();
		pml.saveImage("test21","/home/filippis/Desktop", false);
		
		pml.delete("all", false);
		pdbFileName = "1kj0_21698_55132.pdb";
		pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		
		//check the setDefaults method
		graphPml.setDefaults(new String[] {"nodeSizeRange", "chainColors"});
		graphPml.outputGraph();
		pml.saveImage("test22","/home/filippis/Desktop", false);
		
		pml.closeLog();
		
		SQLC.closeConnection(conn);	

    }

} // end of class testGraph2Pml
