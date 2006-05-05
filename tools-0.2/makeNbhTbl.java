/****************************************************************************
 * Make neighbourhood table                                             
 * Author:      Henning Stehr
 * Date:        21/12/2005
 * Usage:       makenbh [start_index] [end_index]
 * Description: Simple command line program to generate the 'neighbourhoods'
 *              table in the 'newmsdgraph' database. The program will
 *              generate rows in the neighbourhood table for each node in
 *              each graph with start_index <= graph_id <= end_index.
 * Contact:     stehr@molgen.mpg.de
 *****************************************************************************/

import java.sql.*;

public class makeNbhTbl {

/* -------------------------------- constants ------------------------------------------*/
	
	// constants for database connection
	static final String    dbDriver = "org.gjt.mm.mysql.Driver";
	static final String    dbServer = "jdbc:mysql://fliederlila.molgen.mpg.de:3036/newmsdgraph";
	static final String    dbUserName = "henning";
	static final String    dbPassword = "das4ke";

	// constants for graph database
	//static final String    listTable = "henning.myfamily_list";
	//static final String    nodeTable = "henning.myfamily_nodes";
	//static final String    edgeTable = "henning.myfamily_edges";
	static final String    listTable = "list";
	static final String    nodeTable = "nodes";
	static final String    edgeTable = "edges";
	
	
	// constants for neighbourhood table
	//static final String    nbhTable = "henning.myfamily_nbhs";
	static final String    nbhTable = "neighbourhoods_sc_sc";
	static final String    edgeType = "SC_SC";
	static final String    tempNodes = "temp_nodes";
	static final String    tempEdges = "temp_edges";
	static final String    nbhTableDefinition = "(`graph_id` int(10) unsigned default NULL," +
												" `cid` varchar(6) default NULL," +
												" `num` int(5) unsigned default NULL," +
												" `res` char(3) default NULL," +
												" `entry_id` int(5) unsigned default NULL," +
												" `size` int(3) unsigned default NULL, " +
												" `hash` char(20) default NULL, " +
												" `ALA` int(3) unsigned default NULL," +
												" `ARG` int(3) unsigned default NULL," +
												" `ASN` int(3) unsigned default NULL," +
												" `ASP` int(3) unsigned default NULL," +
												" `CYS` int(3) unsigned default NULL," +
												" `GLN` int(3) unsigned default NULL," +
												" `GLU` int(3) unsigned default NULL," +
												" `GLY` int(3) unsigned default NULL," +
												" `HIS` int(3) unsigned default NULL," +
												" `ILE` int(3) unsigned default NULL," +
												" `LEU` int(3) unsigned default NULL," +
												" `LYS` int(3) unsigned default NULL," +
												" `MET` int(3) unsigned default NULL," +
												" `PHE` int(3) unsigned default NULL," +
												" `PRO` int(3) unsigned default NULL," +
												" `SER` int(3) unsigned default NULL," +
												" `THR` int(3) unsigned default NULL," +
												" `TRP` int(3) unsigned default NULL," +
												" `TYR` int(3) unsigned default NULL," +
												" `VAL` int(3) unsigned default NULL," +
												" CONSTRAINT PRIMARY KEY (graph_id, cid, num));";
	
	// This array is used to convert an amino acid number to its three letter code.
	// The numbering is from 1 to 20 to ensure compatibility with function aa_three2num().
	static final String[]    aa_num2three = {"???", "ALA", "ARG", "ASN", "ASP", "CYS",
		                                            "GLN", "GLU", "GLY", "HIS", "ILE",
		                                            "LEU", "LYS", "MET", "PHE", "PRO",
		                                            "SER", "THR", "TRP", "TYR", "VAL"};
	
/* ------------------------------ class variables -------------------------------------*/	

	static Connection conn; // global database connection
	
/* -------------------------------- functions -----------------------------------------*/	
	
	// return the number of an amino acid by its three letter code
	public static int aa_three2num(String threeLetterCode) {
		if(threeLetterCode.equalsIgnoreCase("ALA")) return 1;
		if(threeLetterCode.equalsIgnoreCase("ARG")) return 2;
		if(threeLetterCode.equalsIgnoreCase("ASN")) return 3;
		if(threeLetterCode.equalsIgnoreCase("ASP")) return 4;
		if(threeLetterCode.equalsIgnoreCase("CYS")) return 5;
		if(threeLetterCode.equalsIgnoreCase("GLN")) return 6;
		if(threeLetterCode.equalsIgnoreCase("GLU")) return 7;
		if(threeLetterCode.equalsIgnoreCase("GLY")) return 8;
		if(threeLetterCode.equalsIgnoreCase("HIS")) return 9;
		if(threeLetterCode.equalsIgnoreCase("ILE")) return 10;
		if(threeLetterCode.equalsIgnoreCase("LEU")) return 11;
		if(threeLetterCode.equalsIgnoreCase("LYS")) return 12;
		if(threeLetterCode.equalsIgnoreCase("MET")) return 13;
		if(threeLetterCode.equalsIgnoreCase("PHE")) return 14;
		if(threeLetterCode.equalsIgnoreCase("PRO")) return 15;
		if(threeLetterCode.equalsIgnoreCase("SER")) return 16;
		if(threeLetterCode.equalsIgnoreCase("THR")) return 17;
		if(threeLetterCode.equalsIgnoreCase("TRP")) return 18;
		if(threeLetterCode.equalsIgnoreCase("TYR")) return 19;
		if(threeLetterCode.equalsIgnoreCase("VAL")) return 20;
		return -1;
	}
	
	public static String nbhGetHashValue(int[] nbh) {
		
		char		 currentChar;
		StringBuffer resultBuffer = new StringBuffer("????????????????????");
		
		// nbh[i] is defined for i=1..20
		for(int i = 1; i <= 20; i++) {
			currentChar = '?';
			if(nbh[i] < 0) 
				System.err.println("Error in function nbhGetHashValue: " +
					               "array entry nbh[" + i + "] = " + nbh[i] +" is negative.");
			//else if(nbh[i] <= 9) currentChar = '0' + nbh[i];
			//else if(nbh[i] <= 9 + 26) currentChar = 'A' + nbh[i];
			else if(nbh[i] > 26 + 9)
				System.err.println("Error in function nbhGetHashValue: " +
						           "array entry nbh[" + i + "] = " + nbh[i] + " is too big.");
			else currentChar = Character.forDigit(nbh[i], 26 + 10);
			resultBuffer.setCharAt(i-1, currentChar);
		}
		return resultBuffer.toString();
	}
	
	
	// open connection to database
	public static void openDBConnection() {
		
		try { // try to load driver
		    Class.forName(dbDriver).newInstance(); 
		    try {
				conn = DriverManager.getConnection(dbServer, dbUserName, dbPassword);
		    } catch (SQLException E) {
		    	System.err.println("Error: Unable to open database connection.");
				System.err.println("SQLException: " + E.getMessage());
				System.err.println("SQLState:     " + E.getSQLState());
				System.err.println("VendorError:  " + E.getErrorCode());
		    } // end try/catch connection 
		} // end try load Driver 
		catch (Exception E) {
		    System.err.println("Error: Unable to load database driver.");
		    E.printStackTrace();
		} // end catch 
		
		return;
	}
	
	// close database connection
	public static void closeDBConnection() {
		
	try {
		conn.close();
	} catch (SQLException E) {
		System.out.println("Error: Unable to close connection.");
		System.out.println("SQLException: " + E.getMessage());
		System.out.println("SQLState:     " + E.getSQLState());
		System.out.println("VendorError:  " + E.getErrorCode());
		System.exit(1);
	} // end try/catch connection 
		
		return;
	}	
	
	// create new neighbourhood table if not exists
	public static void createNeighbourhoodTableIfNotExists() {
		Statement stmt;
		String query;
		try {
			query = "CREATE TABLE IF NOT EXISTS " + nbhTable + " " + nbhTableDefinition + ";";
		    stmt = conn.createStatement();	
		    stmt.execute(query);
		} catch (SQLException E) {
			System.err.println("Error: Unable to create neighbourhood table.");
		    System.err.println("SQLException: " + E.getMessage());
		    System.err.println("SQLState:     " + E.getSQLState());
		    System.err.println("VendorError:  " + E.getErrorCode());
		} // end catch		
		
		return;
	}
	
	// insert neighbourhood entries for each graph with from <= graph_id <= to
	public static void insertNeighbourhoods(int fromIdx, int toIdx) {
		Statement   graphsStmt,
		            nodesStmt,
		            edgesStmt;
		ResultSet   graphsRs,
		            nodesRs,
					edgesRs;
		String      query,
		            chainId = "?",
		            resType,
		            nbhResType,
		            hash;
		int         graphId = 0,
		            entryId,
		            resNum = 0,
		            numGraphs = 0,
		            numNodes = 0,
//		            existingNodes,
		            resTypeNum,
		            nbhResCount,
		            nbhResTypeNum,
		            nbhSize;
		int[]       nbh;
		
		try {
			// get number of graphs
			graphsStmt = conn.createStatement();
			query = "SELECT COUNT(*) FROM " + listTable + 
			        " WHERE graph_id >= " + fromIdx + " AND graph_id <= " + toIdx + ";";
		    graphsRs = graphsStmt.executeQuery(query);
		    if(graphsRs.next()) {
		    	numGraphs = graphsRs.getInt(1);
		    }
		    graphsRs.close();		    
		    graphsStmt.close();
		    System.out.println("Processing " + numGraphs + (numGraphs==1?" graph...":" graphs..."));
			
			// iterate over all graphs
			graphsStmt = conn.createStatement();
			query = "SELECT graph_id, entry_id FROM " + listTable + 
			        " WHERE graph_id >= " + fromIdx + " AND graph_id <= " + toIdx + ";";
		    graphsRs = graphsStmt.executeQuery(query);
		    numGraphs = 0;
		    while(graphsRs.next()) {
		    	graphId = graphsRs.getInt(1);
		    	entryId = graphsRs.getInt(2);
		    	
		    	// check whether neighbourhoods for this graph already exist
//		    	nodesStmt = conn.createStatement();
//		    	query = "SELECT COUNT(*) FROM " + nbhTable + " WHERE graph_id = " + graphId + ";";
//		    	nodesRs = nodesStmt.executeQuery(query);
//		    	if(nodesRs.next()) {
//		    		existingNodes = nodesRs.getInt(1);
//		    		if(existingNodes > 0) {
//		    			System.out.println("Warning: Neighbourhoods for the graph with graph_id = "
//		    				           		+ graphId + " already exist.");
//		    		}
//		    	}
//		    	nodesRs.close();
//		    	nodesStmt.close();
		    	 	
		    	// extract all nodes and edges for this graph into temporary tables
		    	nodesStmt = conn.createStatement();
		    	nodesStmt.addBatch("DROP TEMPORARY TABLE IF EXISTS " + tempNodes + ";");
		    	nodesStmt.addBatch("DROP TEMPORARY TABLE IF EXISTS " + tempEdges + ";");
		    	nodesStmt.addBatch("CREATE TEMPORARY TABLE " + tempNodes + 
		    			           " SELECT * FROM " + nodeTable + " WHERE graph_id = " + graphId + ";");
		    	nodesStmt.addBatch("CREATE TEMPORARY TABLE " + tempEdges + 
 			                       " SELECT * FROM " + edgeTable + " WHERE graph_id = " + graphId + 
 			                       " AND " + edgeType + " > 0;"); // take only relevant edges
		    	// if necessary create indices to speed up queries
		    	nodesStmt.executeBatch();
		    	nodesStmt.close();
		    	
		    	// iterate over all nodes for this graph
		    	nodesStmt = conn.createStatement();		    	
		    	query = "SELECT cid, num, res FROM "+ tempNodes + ";"; // graph_id matches automatically
		    	nodesRs = nodesStmt.executeQuery(query);
		    	while(nodesRs.next()) {
		    		chainId = nodesRs.getString(1);
		    		resNum = nodesRs.getInt(2);
		    		resType = nodesRs.getString(3);
		    		resTypeNum = aa_three2num(resType);
		    		if(resTypeNum < 1 || resTypeNum > 20) {
		    			System.err.println("Warning: Non standard amino acid " + resType + 
		    					          " at position " + resNum + " in chain " + chainId + 
		    					          " of graph " + graphId + ".");
		    		}
		    				    		
		    		// get residue counts from edges table
		    		edgesStmt = conn.createStatement();
		    		query = "SELECT j_res AS res, COUNT(*) AS count FROM "+ tempEdges + 
		    		        " WHERE i_cid = '" + chainId + "' AND i_num = " + resNum + 
		    		        " GROUP BY j_res;";
		    		//System.out.println(query);
		    		edgesRs = edgesStmt.executeQuery(query);
		    		
		    		nbh = new int[21]; // set all residue counts to zero
		    		nbhSize = 0;
		    		while(edgesRs.next()) {
		    			nbhResType = edgesRs.getString(1);
		    			nbhResCount = edgesRs.getInt(2);
		    			nbhResTypeNum = aa_three2num(nbhResType);
		    			if(nbhResTypeNum >= 1 && nbhResTypeNum <= 20) {
		    				nbh[nbhResTypeNum] = nbhResCount;
		    				nbhSize += nbhResCount;
		    			}
		    		}
		    		edgesRs.close();
		    		edgesStmt.close();
		    		
		    		// update nbh table
		    		hash = nbhGetHashValue(nbh);
		    		query = "INSERT INTO " + nbhTable + " (graph_id, cid, num, res, entry_id, size, hash, "
		    		                       + "ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, "
		    		                       + "LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL) " + 
		    		        "VALUES (" + graphId + ",'" + chainId + "'," + resNum + ",'" 
		    		                   + resType + "'," + entryId + "," + nbhSize + ",'" + hash + "',"
		    		                   + nbh[1]  + "," + nbh[2]  + "," + nbh[3]  + "," + nbh[4]  + "," + nbh[5]  + ","
		    		                   + nbh[6]  + "," + nbh[7]  + "," + nbh[8]  + "," + nbh[9]  + "," + nbh[10] + ","
		    		                   + nbh[11] + "," + nbh[12] + "," + nbh[13] + "," + nbh[14] + "," + nbh[15] + ","
		    		                   + nbh[16] + "," + nbh[17] + "," + nbh[18] + "," + nbh[19] + "," + nbh[20] + ");";
		    		edgesStmt = conn.createStatement();
		    		if(edgesStmt.executeUpdate(query) != 1) {
		    			System.err.println("Error: Insert into table " + nbhTable + " failed.");
		    			System.exit(1);
		    		}
		    		edgesStmt.close();
		    		numNodes++;
		    	} // end of while over nodes
		    	nodesRs.close();
		    	nodesStmt.close();
		       	numGraphs++;
		       	if(numGraphs % 1000 == 0) {
		       		System.out.print(numGraphs + " ");
		       		System.out.flush();
		       	}
		       	System.gc(); // clean up
		    } // end while over all graphs
		    graphsRs.close();		    
		    graphsStmt.close();
		} catch (SQLException E) {
			System.err.println("Error in function insertNeighbourhoods:");
			System.err.println("graphId = "+ graphId +", chainId = "+ chainId +", resNum = "+ resNum);
		    System.err.println("SQLException: " + E.getMessage());
		    System.err.println("SQLState:     " + E.getSQLState());
		    System.err.println("VendorError:  " + E.getErrorCode());
		} // end catch
		
		System.out.println();
		System.out.println("Extracted "+ numNodes + " neighbourhood" + (numNodes==1?"":"s") 
				        + " from " + numGraphs + " graph" + (numGraphs==1?"":"s") + ".");			
		return;
	}
	
	public static void createIndicesOnNeighbourhoodTable() {
		// index on graph_id
		// index on graph_id, cid, num
		// index on res
		// index on size
		// index on hash
	}
	
/* ----------------------------------- main -------------------------------------------*/	
	
	public static void main(String[] args) {
		
		int from = 0,
		    to = 0;
		
		// read command line parameters		
		if ( args.length != 2) {
		    System.err.println("Usage: makeNbhTbl [start_index] [end_index]");
		    System.exit(1);
		} else {
		
		// read parameter 1
			
			try {
				from = Integer.parseInt(args[0]); 
				if(from <= 0) { 
					System.err.println("Parameter 1 has to be positive");
					System.exit(1);
				}
			}
			catch(NumberFormatException e) {
				System.err.println("Parameter 1 has to be an integer");
				System.exit(1);
			}

		// read parameter 2
			
			try { 
				to = Integer.parseInt(args[1]);
				if(to <= 0) { 
					System.err.println("Parameter 2 has to be positive");
					System.exit(1);
				} else
				if(from > to) {
					System.err.println("Start index has to be smaller or equal than end index");
					System.exit(1);
				}
			}
			catch(NumberFormatException e) { 
				System.err.println("Parameter 2 has to be an integer");
				System.exit(1);
			}
			
		} 

		// execute		
		openDBConnection();
		createNeighbourhoodTableIfNotExists();
		insertNeighbourhoods(from, to);
		createIndicesOnNeighbourhoodTable();
		closeDBConnection();
		
	} // end of main

} // end of class makeNbhTbl
