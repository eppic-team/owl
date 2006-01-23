package tools; 

import java.sql.*;
import java.io.*;
import java.util.*;

/* find the maximum Spanning Tree structure data 
    author: Michael Lappe, 2005-04-15  
*/ 

public class calcMST  
{

    static int graphIDX, reportLevel=0, pre_num=0, act_num=0, root_num=0; 
    static String pre_cid="", act_cid="", root_cid=""; 
    static Connection C = null; // The Connection to the SQL-database 
    
    public static void setConnection( Connection xC) { C = xC; }
   
    public static void setRootCID( String CID) { root_cid = CID; }
    public static String getRootCID() { return root_cid; }

    public static void setRootNUM( int NUM) { root_num = NUM; }
    public static int getRootNUM() { return root_num; }

    public static void setGraphIDX( int GIDX) { graphIDX = GIDX; }
    public static int getGraphIDX() { return graphIDX; }


    /** Sets the reportLevel, so that only messages below or equal the level given are displayed.
	The following levels are available: (0 is default - no messages) 
	3 - talkative mode 
	2 - normal mode 
	1 - basic mode 
	0 - silent mode
	This means by a higher level you get more detailed information, while on a reportLevel 
	of 0 only the most basic results and errors are displayed */
    public static void setReportLevel( int newLevel) { reportLevel = newLevel; }
    /** gets the currently set reportLevel */ 
    static int getReportLevel() { return reportLevel; }

    static void report(String text, int level) {
	if (level<=reportLevel) {
	    System.out.print( text);
	} // end if level of the curent message is above the current level 
    } // end of report

    static void reportln(String text, int level) {
	if (level<=reportLevel) {
	     System.out.println( text);
	} // end if level of the curent message is above the current level 
    } // end of report


    public static void initSums() {
	try {
	  
	    reportln("Total Sums for Graph IDX#"+getGraphIDX()+" are initialised.", 2); 
	    Statement S;
	    S = C.createStatement();
	    S.executeQuery( "update spath_nodes set used_total=0, dist_total=0, reached=0 where IDX="+getGraphIDX()+";");
	    S.close();
	    
	    S = C.createStatement();
	    S.executeQuery( "update spath_edges set used_total = 0 where IDX="+getGraphIDX()+";");
	    S.close();
	} catch (SQLException E) {
	    System.out.println("SQLException:\t" + E.getMessage());
	    System.out.println("SQLState:\t" + E.getSQLState());
	    System.out.println("VendorError:\t" + E.getErrorCode());
	} 
    } // end initStartTable()


    /** Initialising the calculation by first emptying all status and then marking adjacent edges 
	of the start-node as Q 
	back-edges to the start node are excluded by marking them as B 
    */
    public static void initCounts() {
	reportln("Counts for Graph IDX#"+getGraphIDX()+" are initialised.", 2); 
	try {
	    // empty the edge table  
	    Statement S = C.createStatement();
	    S.executeQuery( "UPDATE spath_edges set i_stat=\'R\', j_stat=\'R\', used=0 where IDX="+getGraphIDX()+";");
	    // and Initialise it with the start node
	    S = C.createStatement();
	    S.executeQuery( "UPDATE spath_nodes set stat=\'R\', pre_cid=\"\", pre_num=0, used=0, dist=0 where IDX="+getGraphIDX()+";");
	   
    	} catch (SQLException E) {
	    System.out.println("SQLException:\t" + E.getMessage());
	    System.out.println("SQLState:\t" + E.getSQLState());
	    System.out.println("VendorError: \t" + E.getErrorCode());
	} 
    } // end 


    public static void addSums() {
	try {
	    Statement S;
	    
	    S = C.createStatement();
	    S.executeQuery( "update spath_nodes set used_total=used_total+used, dist_total=dist_total+dist where IDX="+getGraphIDX()+";");
	    S.close();
	    
	    S = C.createStatement();
	    S.executeQuery( "update spath_edges set used_total = used_total+used where IDX="+getGraphIDX()+";");
	    S.close();
	} catch (SQLException E) {
	    System.out.println("SQLException:\t" + E.getMessage());
	    System.out.println("SQLState:\t" + E.getSQLState());
	    System.out.println("VendorError:\t" + E.getErrorCode());
	} 
    } // end initStartTable()


    /** Performs the maximum Spanning Tree Calculation from the root 
     */
    public static void performMST() {
	int rootdist = 1; 
	initCounts(); 
	act_cid = getRootCID();
	act_num = getRootNUM();
	pre_cid = act_cid; 
	pre_num = act_num; 
	reportln("Performing MaximumSpanningTree (MST) calculation from root node ["+act_cid+":"+act_num+"]", 2);

	while (act_num > 0) { // while we get a NodeID out of Q 	  
	    report(".", 0); 
	    reportln( "\n----------------> ["+act_cid+":"+act_num+"]", 2);
	    addNode2P( rootdist);
	    backtrack(); 
	    followMaxEdgeP2Q(); // setting new (act & pre) cid:num 
	    reportln("new node is ["+act_cid+":"+act_num+"]", 3);
	    rootdist = 0; // root to root distance correction only on first run
	} // end while there are nodes left in Q 
    } // end of performMST

    /*
      get the distance of the predescessing node 
    */ 
    public static int getPreDistance() {
	String Query;
	int distance = 0; 
	try {
	    if (! ( pre_cid.equals(root_cid) && (pre_num==root_num) ) )  {

		// getting the distance of the predecessor
		Statement S = C.createStatement();
		Query = "select dist from spath_nodes where cid=\'"+pre_cid+"\' and num="+pre_num+" and IDX="+getGraphIDX()+" limit 1;";
		// reportln( ">2P> "+Query, 3); 
		ResultSet R = S.executeQuery( Query);
		if ( R.next()) distance = R.getInt( 1);  
		R.close(); 
		S.close(); 
	    }
    	} catch (SQLException E) {
	    System.out.println("SQLException:\t" + E.getMessage());
	    System.out.println("SQLState:\t" + E.getSQLState());
	    System.out.println("VendorError: \t" + E.getErrorCode());
	} 
	return distance; 
    } // end 


    /* 
       Moving actual Node from Q to P 
    */
    public static void addNode2P( int root) {
	String Query; 
	reportln("adding node "+act_cid+":"+act_num+" to P", 2); 
	try {
	    // getting the distance of the predecessor
	    int dist = getPreDistance()-root; 
	    // updating the actual Node 
	    Statement S = C.createStatement();
	    Query = "UPDATE spath_nodes set stat=\'P\', pre_cid=\'"+pre_cid+"\', pre_num="+pre_num+", reached=reached+1, dist="+(dist+1)+" where IDX="+getGraphIDX()+" and cid=\'"+act_cid+"\' and num="+act_num+";"; 
	    // reportln( ">2P> "+Query, 3); 
	    S.executeQuery( Query);
	    // mark outgoing edges 
	    S = C.createStatement();
	    S.executeQuery( "UPDATE spath_edges set i_stat=\'P\' where i_cid=\'"+act_cid+"\' and i_num="+act_num+" and IDX="+getGraphIDX()+";");
	    S.close(); 
	     // mark incoming edges 
	    S = C.createStatement();
	    S.executeQuery( "UPDATE spath_edges set j_stat=\'P\' where j_cid=\'"+act_cid+"\' and j_num="+act_num+" and IDX="+getGraphIDX()+";");
	    S.close(); 
    	} catch (SQLException E) {
	    System.out.println("SQLException:\t" + E.getMessage());
	    System.out.println("SQLState:\t" + E.getSQLState());
	    System.out.println("VendorError: \t" + E.getErrorCode());
	} 
    } // end 


    static void followMaxEdgeP2Q() {
	
	reportln("\nRetrieving the max.edge P -> Q", 2);
	try {
	    Statement S = C.createStatement();

	    ResultSet R = S.executeQuery("select i_cid, i_num, SC_SC, j_cid, j_num from spath_edges where IDX="+getGraphIDX()+" and i_stat=\"P\" and j_stat=\"R\" and SC_SC>0 order by SC_SC DESC limit 1;");

	    if ( R.next() ) {
		pre_cid = R.getString( 1); 
		pre_num = R.getInt( 2);
		act_cid = R.getString( 4); 
		act_num = R.getInt( 5);
		reportln("nextEdge ["+pre_cid+":"+pre_num+"] --("+R.getInt(3)+")--> ["+act_cid+":"+act_num+"]", 2);
	    } else {
		pre_cid =""; 
		pre_num = 0; 
		act_cid =""; 
		act_num = 0;
	    }
	    R.close(); 
	    S.close();	  
	} catch (SQLException E) {
	    System.out.println("SQLException:\t" + E.getMessage());
	    System.out.println("SQLState:\t" + E.getSQLState());
	    System.out.println("VendorError:\t" + E.getErrorCode());
	} 
	    
    } // end followMaxEdgeP2Q()


    /** backtracking from actual node to root 
     */ 

    static void backtrack() {
	int s_num=0, t_num=act_num, weight=0, dista, counter=0;
	String s_cid="", t_cid=act_cid, Query; 
	Statement S, Stmt; 
	ResultSet R, RS; 
	reportln("backtracking path from ["+act_cid+":"+act_num+"] to root ["+root_cid+":"+root_num+"]", 3);
	reportln("#n target -(dist)-> source", 3); 
	// report the path that we have found 
	try {	
	    while ( !( s_cid.equals(root_cid) && (s_num==root_num) && t_cid.equals(root_cid) && (t_num==root_num)) ) {  // while not back to root yet  
		counter ++; 
		Stmt = C.createStatement();
		RS = Stmt.executeQuery( "select pre_cid, pre_num, dist from spath_nodes where IDX="+getGraphIDX()+" and cid=\'"+t_cid+"\' and num="+t_num+";");
		if ( RS.next()) {  // get the predecessor 
		    s_cid = RS.getString(1); // previous node cid 
		    s_num = RS.getInt(2); // previous node num  
		    dista = RS.getInt(3); // distance 
		    reportln("#"+counter+" ["+t_cid+":"+t_num+"] -("+dista+")-> ["+s_cid+":"+s_num+"]", 2); 
		    // increase usage source node 
		    S = C.createStatement();
		    S.executeQuery( "update spath_nodes set used = used+1 where IDX="+getGraphIDX()+" and cid =\'"+s_cid+"\' and num="+s_num+";");
		    S.close();

		    // increase usage edge 
		    S = C.createStatement();
		    S.executeQuery( "update spath_edges set used = used+1 where IDX="+getGraphIDX()+" and i_cid =\'"+s_cid+"\' and i_num="+s_num+" and j_cid =\'"+t_cid+"\' and j_num="+t_num+";");
		    S.close();
		    t_cid = s_cid; 
		    t_num = s_num; 
		} else {
		    reportln("\n",2); 
		    break;
		} // end if R.next()
	    } // end while we haven't reached back to the source 
	    
	} catch (SQLException E) {
	    System.out.println("SQLException:\t" + E.getMessage());
	    System.out.println("SQLState:\t" + E.getSQLState());
	    System.out.println("VendorError:\t" + E.getErrorCode());
	} 
    } // end of reportShortestPath()

} // end class calcPathway 
