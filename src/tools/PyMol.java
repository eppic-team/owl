package tools;

import java.sql.*;
import java.io.*;
import java.util.*;
import java.text.*;

/**
 * Package:		tools
 * Class: 		PyMol
 * Author:		Ioannis Filippis, filippis@molgen.mpg.de
 * Date:		20/03/2006
 *
 * This class serves as a simple PyMol API. Most methods send PyMol commands 
 * to the PrintWriter taken by the constructor, while others return pymol 
 * selection strings. There is also the possibility to send the atomic
 * coordinates of a specific model of a macromolecule directly from msdsd
 * and load them on-the-fly.
 *  
 * Notes: 
 * - variable: this boolean parameter appears in many methods. The idea is that if
 * 	 the objectName provided is really a PyMol object or selection, then it must be
 * 	 quoted. However, it may be just a string variable, so it must be evaluated.In 
 * 	 that case parameter variable should be true.
 * - This is not a full PyMol API. Existing methods wrap basic pymol commands
 * 	 in a simple, rather stupid way to facilitate contact graph visualization.
 * - If you want to use cgo edges, create a .pymolrc file in your home directory
 *   with the command "run /project/StruPPi/PyMolAll/pymol/scripts/ioannis/graph.py"
 * 
 * Changelog:
 * 27/03/06 modified by IF (refresh method plus refresh added in saveImage)
 * 23/03/06 modified by IF (python string methods added plus removeList,concatList methods)
 * 22/03/06 modified by IF (objectNameQuotes class variable replaced by 
 * 			method parameter variable) 
 * 20/03/06 first created by IF
 */
public class PyMol {

	// Out: the PrintWriter where PyMol commands are sent to
	// attrs: a HashMap with keys the PyMol state variables allowed to be set
	//		  by the set class. Values are integers (could be used in a "switch")
	// DF: used to format decimal numbers so commands issued are readable
    private PrintWriter Out = null;
    private HashMap<String, Integer> attrs = new HashMap<String, Integer>();
    private DecimalFormat DF = new DecimalFormat("#0.000");

 
    
    // constructor
    public PyMol(PrintWriter out) {
    	
    	this.Out = out;
    	
    	attrs.put("sphere_transparency", new Integer(1));
    	attrs.put("sphere_scale", new Integer(2));
    	attrs.put("dash_gap", new Integer(3));
    	attrs.put("dash_width", new Integer(4));
    	attrs.put("transparency", new Integer(5));
    	
    }

    /**
     * sends the atom lines of a model (modelId) of a biological unit (assemblyId) of a protein
     * (accessionCode) directly from msdsd. In this way, the structure is loaded wihout temporary
     * pdb files. The structure object is named using the pattern accessionCode_assemblyId_modelId.
     * 
     * Notes:
     * - A connection file is needed to connect to msdsd
     * - The chain_pdb_code is used in the chainID field in the atom line, while the chain_code is used 
     * 	 in the segID field (due to its length). Therefore, "segi" and not "chain" must be used in pymol
     * 	 selections (all methods take care of that based on the boolean parameter msdsd)
     * - There are two versions of sendAtomLines. One that takes the atomic coordinates from the 
     * 	 partial atom_data tables (needs the table number e.g. 1 for atom_data_1, but is faster), 
     *	 while the other uses the merged table (really slow - should be avoided)
     * - In general, the method is slow. The use of temporary files should be preferred. Have a look
     *	 at Msdsd2Pdb class.
     */ 
    public void sendAtomLines(String accessionCode, int assemblyId, int modelId, int atomDataTblNum, String connFile) {

		String query = "SELECT CONCAT("+
		    "RPAD(\"ATOM\", 6, \" \"), "+
		    "LPAD(serial, 5, \" \"), "+
		    "\" \", "+
		    "LPAD(chem_atom_name, 4, \" \"), "+
		    "IF(alt_code IS NULL, \" \", alt_code), "+
		    "code_3_letter, "+
		    "\" \", "+
		    "IF(chain_pdb_code IS NULL, \" \", chain_pdb_code), "+
		    "LPAD(residue_serial, 4, \" \"), "+
		    "IF(residue_pdb_insert_code IS NULL, \" \", residue_pdb_insert_code), "+
		    "REPEAT(\" \", 3), "+
		    "LPAD(x, 8, \" \"), "+
		    "LPAD(y, 8, \" \"), "+
		    "LPAD(z, 8, \" \"), "+
		    "LPAD(occupancy, 6, \" \"), "+
		    "REPEAT(\" \", 6), "+
		    "REPEAT(\" \", 6), "+
		    "RPAD(chain_code, 4, \" \") "+
		    ") AS atom_lines FROM msdsd.atom_data_"+atomDataTblNum+" WHERE "+
		    "(assembly_id = "+assemblyId+") AND "+
		    "(model_id = "+modelId+") AND "+
		    "((alt_code = \"A\") OR (alt_code IS NULL)) AND "+
		    "(pdb_group = \"A\") "+
		    "ORDER BY chain_code, residue_serial, serial;";
	
		mySQLConnect SQLC = new mySQLConnect();
		SQLC.readConnectionFile(connFile);
		Connection conn = SQLC.openConnection();
	
		Statement S;
		ResultSet R;
	
		Out.print("cmd.read_pdbstr(\"\"\"");
		
		try { 
		    S = conn.createStatement();	
		    R = S.executeQuery(query);
		    while (R.next()) {
		    	Out.println(R.getString(1)+"\\");
		    }
		    R.close();
		    S.close();   
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch 
		
		SQLC.closeConnection(conn);
	
		Out.println("END\"\"\", \""+accessionCode+"_"+assemblyId+"_"+modelId+"\")");

    }
    
    /**
     * sends the atom lines of a model (modelId) of a biological unit (assemblyId) of a protein
     * (accessionCode) directly from msdsd. In this way, the structure is loaded wihout temporary
     * pdb files. The structure object is named using the pattern accessionCode_assemblyId_modelId.
     * 
     * Notes:
     * - A connection file is needed to connect to msdsd
     * - The chain_pdb_code is used in the chainID field in the atom line, while the chain_code is used 
     * 	 in the segID field (due to its length). Therefore, "segi" and not "chain" must be used in pymol
     * 	 selections (all methods take care of that based on the boolean parameter msdsd)
     * - There are two versions of sendAtomLines. One that takes the atomic coordinates from the 
     * 	 partial atom_data tables (needs the table number e.g. 1 for atom_data_1, but is faster), 
     *	 while the other uses the merged table (really slow - should be avoided)
     * - In general, the method is slow. The use of temporary files should be preferred. Have a look
     *	 at Msdsd2Pdb class.
     */ 
    public void sendAtomLines(String accessionCode, int assemblyId, int modelId, String connFile) {

		String query = "SELECT CONCAT("+
		    "RPAD(\"ATOM\", 6, \" \"), "+
		    "LPAD(serial, 5, \" \"), "+
		    "\" \", "+
		    "LPAD(chem_atom_name, 4, \" \"), "+
		    "IF(alt_code IS NULL, \" \", alt_code), "+
		    "code_3_letter, "+
		    "\" \", "+
		    "IF(chain_pdb_code IS NULL, \" \", chain_pdb_code), "+
		    "LPAD(residue_serial, 4, \" \"), "+
		    "IF(residue_pdb_insert_code IS NULL, \" \", residue_pdb_insert_code), "+
		    "REPEAT(\" \", 3), "+
		    "LPAD(x, 8, \" \"), "+
		    "LPAD(y, 8, \" \"), "+
		    "LPAD(z, 8, \" \"), "+
		    "LPAD(occupancy, 6, \" \"), "+
		    "REPEAT(\" \", 6), "+
		    "REPEAT(\" \", 6), "+
		    "RPAD(chain_code, 4, \" \") "+
		    ") AS atom_lines FROM msdsd.atom_data WHERE "+
		    "(assembly_id = "+assemblyId+") AND "+
		    "(model_id = "+modelId+") AND "+
		    "((alt_code = \"A\") OR (alt_code IS NULL)) AND "+
		    "(pdb_group = \"A\") "+
		    "ORDER BY chain_code, residue_serial, serial;";
	
		mySQLConnect SQLC = new mySQLConnect();
		SQLC.readConnectionFile(connFile);
		Connection conn = SQLC.openConnection();
	
		Statement S;
		ResultSet R;
	
		Out.print("cmd.read_pdbstr(\"\"\"");
		
		try { 
		    S = conn.createStatement();	
		    R = S.executeQuery(query);
		    while (R.next()) {
		    	Out.println(R.getString(1)+"\\");
		    }
		    R.close();
		    S.close();   
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch 
		
		SQLC.closeConnection(conn);
	
		Out.println("END\"\"\", \""+accessionCode+"_"+assemblyId+"_"+modelId+"\")");

    }

    /**
     * loads a pdb file
     * 
     * Notes:
     * - There are two versions. One of them loads the pdb file specified by its name and the directory path, while
     *   the other just uses the filename.
     * - Pdb file is expected to have the extension ".pdb".
     * - The structure object is named using the accession code (extension ".pdb" is trimmed)
     */ 
    public String loadPDB(String pdbFileName, String pdbDir) {
       	
		Out.println("cmd.load(\""+pdbDir+"/"+pdbFileName+"\",  \""+pdbFileName.substring(0, pdbFileName.length()-4)+"\", 1)");
	
		return (pdbFileName.substring(0, pdbFileName.length()-4));

    }
    
    /**
     * loads a pdb file
     * 
     * Notes:
     * - There are two versions. One of them loads the pdb file specified by its name and the directory path, while
     *   the other just uses the filename.
     * - Pdb file is expected to have the extension ".pdb".
     * - The structure object is named using the accession code (extension ".pdb" is trimmed)
     */ 
    public String loadPDB(String pdbFileName) {
       	
		Out.println("cmd.load(\""+pdbFileName+"\",  \""+pdbFileName.substring(0, pdbFileName.length()-4)+"\", 1)");
	
		return (pdbFileName.substring(0, pdbFileName.length()-4));

    }

    /**
     * adds a node as a sphere centered on the Ca of a residue
     * 
     * Notes:
     * - The node object is named n.cid.num (e.g. n.A.15 is the 15th residue-node in chain A)
     * - The msdsd boolean is used to denote whether msdsd is the source of the structure. In that case,
     * 	 segi is used instead of chain
     */ 
    public String addNode(String cid, int num, boolean msdsd) {
	
		String nodeName = "n."+cid+"."+num;
		String nodeSel = selectNode(cid, num, msdsd, true);
		
		Out.println("cmd.create(\""+nodeName+"\", \""+nodeSel+"\")");
		Out.println("cmd.show(\"sphere\", \""+nodeName+"\")");
	
		return nodeName;

    }

    /**
     * adds an edge as a distance object between the Ca's of 2 residues
     * 
     * Notes:
     * - The edge object is named e.i_cid.i_num.j_cid.j_num (e.g. e.A.1.B.10 is the edge
     *   between the 1st residue-node in chain A and the 10th residue in chain B)
     * - The msdsd boolean is used to denote whether msdsd is the source of the structure. In that case,
     * 	 segi is used instead of chain
     */    
    public String addEdge(String i_cid, int i_num, String j_cid, int j_num, boolean msdsd) {
    	
  

		String edgeName = "e."+i_cid+"."+i_num+"."+j_cid+"."+j_num;
		String iNodeSel = selectNode(i_cid, i_num, msdsd, true);
		String jNodeSel = selectNode(j_cid, j_num, msdsd, true);
	
		Out.println("cmd.distance(\""+edgeName+"\", \""+iNodeSel+"\", \""+jNodeSel+"\")");
		Out.println("cmd.hide(\"labels\")");
	
		return edgeName;

    }
    
    /** Creates an edge between the C-alpha atoms of the given residues in the given chain. 
     *  The selection in pymol will be names pdbFileName+"Sel"+selNum 
     */
    public void setDistance(int resi1, int resi2, String pdbFilename, int selNum, String chain_pdb_code){   	
    	Out.println("distance "+ pdbFilename+"Sel"+selNum+" , chain "+chain_pdb_code+" and resi " + resi1 + " and name ca, chain "+chain_pdb_code+" and resi " + resi2 + " and name ca;");
    }
    
    /** Creates an edge between the C-alpha atoms of the given residues.
     *  Use this variant if there is only one unnamed chain in the current structure.
     *  The selection in pymol will be names pdbFileName+"Sel"+selNum 
     */
    public void setDistance(int resi1, int resi2, String pdbFilename, int selNum){  	
    	Out.println("distance "+ pdbFilename+"Sel"+selNum+" , resi " + resi1 + " and name ca, resi " + resi2 + " and name ca;");
    }
    
    /**
     * adds an edge as a sausage Compiled Graphic Object (cgo) between the Ca's of 2 residues
     * 
     * Notes:
     * - The edge object is named e.i_cid.i_num.j_cid.j_num (e.g. e.A.1.B.10 is the edge
     *   between the 1st residue-node in chain A and the 10th residue in chain B)
     * - If the graph is directed, then only half of the user-formatted sausage will be
     * 	 drawn towards the target node and the rest will be drawn as a thin sausage. 
     * - rgb(color) and half_rgb(half_color) define the color of the two parts of the edge (if directed).
     * 	 There are two versions of addCgoEdge, one with string (color) and one with array of doubles
     * 	 (rgb) parameters.
     * - The msdsd boolean is used to denote whether msdsd is the source of the structure. In that case,
     * 	 segi is used instead of chain.
     */    
    public String addCgoEdge(String i_cid, int i_num, String j_cid, int j_num, double[] rgb, double dashGap, double dashLength, double dashRadius, boolean directed, double[] half_rgb, boolean msdsd) {

		String edgeName = "e."+i_cid+"."+i_num+"."+j_cid+"."+j_num;
		String iNodeSel = selectNode(i_cid, i_num, msdsd, true);
		String jNodeSel = selectNode(j_cid, j_num, msdsd, true);
	
		Out.println("edge name="+edgeName+", i_node="+iNodeSel+", j_node="+jNodeSel+", r="+DF.format(rgb[0])+", g="+DF.format(rgb[1])+", b="+DF.format(rgb[2])+", dg="+DF.format(dashGap)+", dl="+DF.format(dashLength)+", dr="+DF.format(dashRadius)+", dir="+(directed?1:0)+", dir_r="+DF.format(half_rgb[0])+", dir_g="+DF.format(half_rgb[1])+", dir_b="+DF.format(half_rgb[2])+"");
	
		return edgeName;

    }
    
    /**
     * adds an edge as a sausage Compiled Graphic Object (cgo) between the Ca's of 2 residues
     * 
     * Notes:
     * - The edge object is named e.i_cid.i_num.j_cid.j_num (e.g. e.A.1.B.10 is the edge
     *   between the 1st residue-node in chain A and the 10th residue in chain B)
     * - If the graph is directed, then only half of the user-formatted sausage will be
     * 	 drawn towards the target node and the rest will be drawn as a thin sausage. 
     * - rgb(color) and half_rgb(half_color) define the color of the two parts of the edge (if directed).
     * 	 There are two versions of addCgoEdge, one with string (color) and one with array of doubles
     * 	 (rgb) parameters.
     * - The msdsd boolean is used to denote whether msdsd is the source of the structure. In that case,
     * 	 segi is used instead of chain.
     * - Have a look at the /project/StruPPi/PyMolAll/pymol/scripts/ioannis/graph.py with the 
     * 	 implementation of the edge command
     */    
    public String addCgoEdge(String i_cid, int i_num, String j_cid, int j_num, String color, double dashGap, double dashLength, double dashRadius, boolean directed, String half_color, boolean msdsd) {

		String edgeName = "e."+i_cid+"."+i_num+"."+j_cid+"."+j_num;
		String iNodeSel = selectNode(i_cid, i_num, msdsd, true);
		String jNodeSel = selectNode(j_cid, j_num, msdsd, true);
	
		Out.println("edge name="+edgeName+", i_node="+iNodeSel+", j_node="+jNodeSel+", color="+color+", dg="+DF.format(dashGap)+", dl="+DF.format(dashLength)+", dr="+DF.format(dashRadius)+", dir="+(directed?1:0)+", dir_color="+color+"");
	
		return edgeName;

    }
 
    /**
     * deletes an object or selection
     */    
    public void delete(String objectName, boolean variable) {
	
		Out.println("cmd.delete("+((!variable)?"\""+objectName+"\"":objectName)+")");

    }
    
    /**
     * deletes a node
     */  
    public void delNode(String cid, int num) {
	
		String nodeName = "\"n."+cid+"."+num+"\"";
		Out.println("cmd.delete("+nodeName+")");

    }
    
    /**
     * deletes an edge
     */ 
    public void delEdge(String i_cid, int i_num, String j_cid, int j_num) {

		String edgeName = "\"e."+i_cid+"."+i_num+"."+j_cid+"."+j_num+"\"";
		Out.println("cmd.delete("+edgeName+")");

    }

    
    /**
     * turns on atom/bond representation for all bonds for an object or selection
     */ 
    public void myShow(String objectName) {
	
		Out.println("cmd.show('"+objectName+"')");

    }
    
    
    
    /**
     * turns on atom/bond representation specified by what for an object or selection
     */ 
    public void showWhat(String what, String objectName, boolean variable) {

		Out.println("cmd.show(\""+what+"\", "+((!variable)?"\""+objectName+"\"":objectName)+")");
		
    }

    /**
     * turns on atom/bond representation for all bonds for an object or selection
     */ 
    public void show(String objectName, boolean variable) {
	
		Out.println("cmd.show(\"everything\", "+((!variable)?"\""+objectName+"\"":objectName)+")");

    }


    /**
     * turns on node representation
     */     
    public void showNode(String cid, int num) {
	
		String nodeName = "n."+cid+"."+num;
		Out.println("cmd.show(\"sphere\", \""+nodeName+"\")");

    }
    
    /**
     * turns on edge representation
     */     
    public void showEdge(String i_cid, int i_num, String j_cid, int j_num) {

		String edgeName = "e."+i_cid+"."+i_num+"."+j_cid+"."+j_num;
		Out.println("cmd.show(\"everything\", \""+edgeName+"\"");

    }
    
    /**
     * turns off atom/bond representation specified by what for an object or selection
     */ 
    public void hideWhat(String what, String objectName, boolean variable) {

		Out.println("cmd.hide(\""+what+"\", "+((!variable)?"\""+objectName+"\"":objectName)+")");
		
    }
    
    /**
     * turns off atom/bond representation for all bonds for an object or selection
     */ 
    public void hide(String objectName, boolean variable) {
	
		Out.println("cmd.hide(\"everything\", "+((!variable)?"\""+objectName+"\"":objectName)+")");

    }

    
    /**
     * turns off atom/bond representation for all bonds for an object or selection
     */ 
    public void myHide(String objectName) {
	
		Out.println("cmd.hide('"+objectName+"')");

    }
    
    
    /**
     * turns off node representation
     */   
    public void hideNode(String cid, int num) {
	
		String nodeName = "n."+cid+"."+num;
		Out.println("cmd.hide(\"everything\", \""+nodeName+"\")");

    }
    
    /**
     * turns off edge representation
     */   
    public void hideEdge(String i_cid, int i_num, String j_cid, int j_num) {

		String edgeName = "e."+i_cid+"."+i_num+"."+j_cid+"."+j_num;
		Out.println("cmd.hide(\"everything\", \""+edgeName+"\")");

    }

    /**
     * creates a named neighbourhood selection using all atoms within distCutOff Angstrom from 
     * all atoms of the specified residue
     * 
     * Notes:
     * - The neighbourhood selection is named n.cid.num.N (e.g. n.A.1.N is the neighbourhood of
     * 	 the 1st residue-node in chain A)
     * - The msdsd boolean is used to denote whether msdsd is the source of the structure. In that case,
     * 	 segi is used instead of chain.
     */        
    public void selNbrPml(String cid, int num, double distCutOff, boolean msdsd) {

		String nodeNbrName = "n."+cid+"."+num+".N";
		String nodeSel = selectNode(cid, num, msdsd, false);

		Out.println("cmd.select(\""+nodeNbrName+"\", \""+nodeSel + " around "+DF.format(distCutOff)+"\")");

    }

    /**
     * creates a selection called name using the atom selection "encoded" in the objectName
     *  
     * Notes:
     * - If pink is false, the selection display is disabled (pink stuff go away!!!!) 
     */   
    public void selPml(String name, String objectName, boolean pink, boolean variable) {

		Out.println("cmd.select(\""+name+"\", "+((!variable)?"\""+objectName+"\"":objectName)+")");
		if (!pink) { Out.println("cmd.disable(\""+name+"\")"); };

    }
    
    /**
     * returns an atom selection string for a node
     *  
     * Notes:
     * - The CA boolean is used to denote whether the Ca atom should be only selected
     * - The msdsd boolean is used to denote whether msdsd is the source of the structure. In that case,
     * 	 segi is used instead of chain.
     */   
    public String selectNode(String cid, int num, boolean msdsd, boolean CA) {

    	return "("+(msdsd?"segi ":"chain ")+(cid.equals("")?"\"\"":cid)+" and resi "+num+(CA?" and name CA":"")+")";
	
    }
    
    /**
     * returns an atom selection string for a chain
     *  
     * Notes:
     * - The msdsd boolean is used to denote whether msdsd is the source of the structure. In that case,
     * 	 segi is used instead of chain.
     */   
    public String selectChain(String cid, boolean msdsd) {

    	return "("+(msdsd?"segi ":"chain ")+(cid.equals("")?"\"\"":cid)+")";

    }
    
    public void select(String name , String residue_nr){
    	
    	Out.println("select "+ name+", resi " + residue_nr );
    }
    
    
    
    public int set(String objectName, double value, String object) {
    	if (object   == ""){
		
		Out.println("set "+objectName+", " + value);}
    	
    	else {
    		Out.println("set "+objectName+", " + value + " in object " + object);
    	}
		//Out.println("cmd.set(\""+attribute+"\", "+DF.format(value)+", "+((!variable)?"\""+objectName+"\"":objectName)+")");
	
		return 0;
		
    }
    
    
    
    
    /**
     * changes one of the PyMol state variable (attribute) for a specific object or selection
     * and sets it equal to the provided value.
     *  
     * Notes:
     * - If the state variable is not defined in the attrs hashMap, then set returns -1 else 0.
     */   
    public int set(String attribute, double value, String objectName, boolean variable) {

		Integer key = null;
		key = attrs.get(attribute);
	
		if (key == null) { return -1; }
	
		Out.println("cmd.set(\""+attribute+"\", "+DF.format(value)+", "+((!variable)?"\""+objectName+"\"":objectName)+")");
	
		return 0;
		
    }
    
    /**
     * sets the color for a node
     */  
    public void setNodeColor(String color, String nodeName, boolean variable) {

		Out.println("cmd.set(\"sphere_color\", \""+color+"\", "+((!variable)?"\""+nodeName+"\"":nodeName)+")");

    }
    
    /**
     * sets the color for an object or selection
     */ 
    public void setColor(String color, String objectName, boolean variable) {

		Out.println("cmd.color(\""+color+"\", "+((!variable)?"\""+objectName+"\"":objectName)+")");

    }

    /**
     * define a new color or rather color name based on an existing color
     */ 
    public void createColor(String colorName, String color) {

		Out.println("cmd.set_color(\""+colorName+"\", cmd.get_color_tuple(cmd.get_color_index(\""+color+"\")))");

    }

    /**
     * define a new color providing the RGB array of doubles
     */     
    public void createColor(String colorName, double[] rgb) {
	
		Out.println("cmd.set_color(\""+colorName+"\", ["+DF.format(rgb[0])+", "+DF.format(rgb[1])+", "+DF.format(rgb[2])+"])");

    }
    
    /**
     * create/set a string
     */  
    public void initString(String stringName, String value) {

    	Out.println(stringName+" = \""+value+"\"");
	
    }
    
    /**
     * append to String
     */  
    public void appendString(String stringName, String value) {

    	Out.println(stringName+" += \""+value+"\"");
	
    }
    
    /**
     * rstring a String
     */  
    public void rstripString(String stringName, String value) {

    	Out.println("string.rstrip("+stringName+", \""+value+"\")");
	
    }    
    
    /**
     * replace a specific substring within a string with another one
     */  
    public void replaceString(String stringName, String oldValue, String newValue) {

    	Out.println("string.replace("+stringName+", \""+oldValue+"\", \""+newValue+"\")");
	
    }        
    
    /**
     * create/initialize a list
     */  
    public void initList(String listName) {

    	Out.println(listName+" = []");
	
    }
    
    /**
     * append objects to a list
     */  
    public void appendList(String listName, String value) {

    	Out.println(listName+".append(\""+value+"\")");
	
    }
    
    /**
     * remove object from a list
     */  
    public void removeList(String listName, String value) {

    	Out.println(listName+".remove(\""+value+"\")");
	
    }
    
    /**
     * concatenate all members of a list of strings into a string using specific seperator
     */  
    public void concatList(String listName, String sep, String stringName) {

    	Out.println(stringName+" = \""+sep+"\".join("+listName+")");
	
    }
    
    /**
     * iterate through list objects
     * 
     * Notes:
     * - item is the name of the variable containing the object
     * - This is the only method that "print" instead of "println". 
     * 	 In this way, a method can be executed at each iteration e.g. 
     *		PyMol pml = new PyMol(out);
     *		...
     *		pml.iterateList("edges", "edge");
     *		pml.setColor("red", edge, true);
     * - In the previous example, the "edge" parameter in the setColor command
     *   is the name of tha variable holding the actual object. Therefore, it must
     *   be evaluated and so not quoted. This is why setColor is called with true
     *   value for the argument variable
     */
    public void iterateList(String listName, String item) {

		Out.print("for "+item+" in "+listName+":");
	
    }
    
    /**
     * opens a log file
     * 
     * Notes:
     * - Be careful. All log files are saved with extension .pml. Remember to include it in the fileName!!!!!
     * - There are two versions. One of them opens the log file specified by its name and the directory path, while
     *   the other just uses the filename.
     */ 
    public void openLog(String fileName, String dir) {
       	
		Out.println("log_open "+dir+"/"+fileName+".pml");

    }
    
    /**
     * opens a log file
     * 
     * Notes:
     * - Be careful. All log files are saved with extension .pml. Remember to include it in the fileName!!!!!
     * - There are two versions. One of them opens the log file specified by its name and the directory path, while
     *   the other just uses the filename.
     */ 
    public void openLog(String fileName) {
       	
		Out.println("log_open "+fileName+".pml");

    }
    
    /**
     * closes the last opened log file
     */ 
    public void closeLog() {
       	
		Out.println("log_close");

    }
    
    /**
     * get the current view into a string
     */ 
    public void getView(String viewName) {
       	
		Out.println(viewName+" = cmd.get_view()");

    }
    
    /**
     * get the current view into a string from a log file
     */ 
    public void getFileView(String viewName, String fileName) {
       	
    	String viewStr = "";
    	String line = "", viewLine = "";
    	boolean viewFinished = false;
    	
    	try {
    		//read the first view in the file
    		BufferedReader fileIn = new BufferedReader(new FileReader(new File(fileName)));
    		while ( ((line = fileIn.readLine()) != null) && (!viewFinished) ){
    			if (line.equals("_ set_view (\\")) {
    				while ( ((viewLine = fileIn.readLine()) != null) && (!viewFinished) ) {
    					if (viewLine.startsWith("_")) {
    						viewStr = viewStr + viewLine.substring(1,viewLine.length()-1);
    					} else {
    						viewFinished = true;
    					}
    				}
    			}
    		}
    		if (fileIn != null) { fileIn.close(); }
    		viewStr = "("+viewStr+")";
    	}
    	catch (Exception e) { System.out.println(e); }
		Out.println(viewName+" = "+viewStr);

    }
    
    /**
     * set the view based on a string's value
     */ 
    public void setView(String view) {
       	
		Out.println("cmd.set_view("+view+")");

    }
    
    /**
     * writes a png format image file of the current frame
     * 
     * Notes:
     * - rayTraced defines whether the image will be ray-traced or not
     * - There are two versions. One of them saves the image file as specified by its name and the directory path, 
     *   while the other just uses the filename.
     */ 
    public void saveImage(String imageName, String dir, boolean rayTraced) {
	
    	refresh();
		if (rayTraced) { Out.println("ray"); }
		Out.println("png "+dir+"/"+imageName+".png");

    }
    
    /**
     * writes a png format image file of the current frame
     * 
     * Notes:
     * - rayTraced defines whether the image will be ray-traced or not
     * - There are two versions. One of them saves the image file as specified by its name and the directory path, 
     *   while the other just uses the filename.
     */ 
    public void saveImage(String imageName, boolean rayTraced) {
	
    	refresh();
		if (rayTraced) { Out.println("ray"); }
		Out.println("png "+imageName+".png");

    }
    
    /**
     * sources a pymol command script
     * 
     * Notes:
     * - There are two versions. One of them runs the file specified by its name and the directory path, 
     *   while the other just uses the filename.
     */ 
    public void runScript(String fileName, String dir) {

		Out.println("@"+dir+"/"+fileName+".pml");
	
    }
    
    /**
     * sources a pymol command script
     * 
     * Notes:
     * - There are two versions. One of them runs the file specified by its name and the directory path, 
     *   while the other just uses the filename.
     */ 
    public void runScript(String fileName) {

		Out.println("@"+fileName+".pml");
	
    }

    /**
     * zooms on an object or selection
     */     
    public void zoom(String objectName, boolean variable) {
	
		Out.println("cmd.zoom("+((!variable)?"\""+objectName+"\"":objectName)+")");

    }
    
    /**
     * sets the background color
     */    
    public void background(String color) {

		Out.println("cmd.bg_color(\""+color+"\")");
	
    }
    
    /**
     * refresh the scene
     */       
    public void refresh() {

		Out.println("cmd.refresh()");
	
    }
    
    /**
     * some initialization commands for structure visualization
     * 
     * Notes:
     * - These are defined according to Ioannis' preferences
     */    
    public void init() {

		Out.println("cmd.set(\"depth_cue\", 0)");
		Out.println("cmd.set(\"ray_trace_fog\", 0)");
		Out.println("cmd.hide()");
		Out.println("cmd.show(\"cartoon\")");
		Out.println("cmd.cartoon(\"automatic\")");
		Out.println("cmd.set(\"cartoon_flat_sheets\", 0)");
		Out.println("cmd.set(\"cartoon_fancy_helices\", 1)");
		Out.println("cmd.hide(\"spheres\", \"hetatm\")");

    }

}
