import tools.MySQLConnection;

import java.sql.SQLException;
import java.sql.Statement;
import java.sql.ResultSet;

public class writeTargetMoves {

	/**
	 * writes all the moves for a given residue for later scoring. 
	 * @author lappe
	 */
	static int maxRank = 21; // value to replace for non-existence of central redue in the resultvector (rank=0) 
	// higher values should penalize non-existence more
	static int VL=1; // Verbosity Level 
	static String user = "lappe"	; // change user name!!
	static MySQLConnection conn;
	static String backgrndDB = "CASP_decoys"; 
	static String targetNodes = "single_model_node";
	static String targetEdges = "single_model_edge";
	static String targetScore = "single_model_score";
	static int cn1[], cn2[]; 

	public static void main(String[] args) {
		int graph_id=0, num=0, node_id=0; 
		String cid="", res="", sstype=""; 
		String sql;  
		Statement mstmt;  
		ResultSet mrsst; 
		
		try {
			conn = new MySQLConnection("white",user,"nieve", backgrndDB); // the UPPERCASE DB!  
			System.out.println("Writing Target neighborhoods v.0.2. (SC&BB)"); 
			// retrieve node_id   | cid  | res  | sstype defined by graph_id, num  
			// sql = "select node_id, cid, num, res, sstype from "+targetNodes+" where graph_id="+graph_id+" limit 3;"; 
			sql = "select graph_id, node_id, cid, num, res, sstype from "+targetNodes+";"; //  +" limit 1;"; 
			mstmt = conn.createStatement();
			mrsst = mstmt.executeQuery(sql); 
			while (mrsst.next()) {
			    // this is the central node -> get type and secondary structure
				graph_id = mrsst.getInt( 1); 
				node_id = mrsst.getInt( 2); 
				cid = mrsst.getString( 3);
				num = mrsst.getInt( 4); 
				res = mrsst.getString( 5).toUpperCase();
				sstype= mrsst.getString( 6).toUpperCase();
				System.out.println("GraphID "+graph_id+" Chain "+cid+", Central residue ("+node_id+") "+res+":"+num+":"+sstype); 
				writeMoves( graph_id, node_id, cid, num, res, sstype); 		
			} // end if central residue
			mrsst.close();
			mstmt.close(); 			
		} catch (SQLException e) {
			e.printStackTrace();
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
		} // end try/catch 
		System.out.println("fin."); 
	}	// end main 


	
	public static void writeMoves( int cgraph_id, int cnode_id, String ccid, int cnum, String cres, String csstype) {
		int n1=0, n2=0, ni=0, nj=0, i=0, j=0, j_num=0, j_shell=0, j_cnsize=0, j_bb, j_sc;
		int minus, mcn, plus, pcn; 
		String sql, j_res, j_sec, mres, mss, pres, pss, nn, nb; 
		boolean overx = false;
		Statement stmt, mst;  
		ResultSet rsst;
		
		try {
			
			// retrieve the original nbhood into orig_shell
			System.out.println("retrieving first shell ... "); 
			stmt = conn.createStatement();
			stmt.executeUpdate("drop table if exists temp_shell;"); 
			stmt.close(); 
	
			stmt = conn.createStatement();
			// stmt.executeUpdate("create table temp_shell as select i_num, i_res, j_num, j_res, j_sstype, 1 as shell from "+targetEdges+" where graph_id="+cgraph_id+" and i_num="+cnum+";");
			stmt.executeUpdate("create table temp_shell as select i_num, i_res, j_num, j_res, j_sstype, 1 as shell, weight_bb as BB, weight_SC as SC " +
					"from "+targetEdges+" where graph_id="+cgraph_id+" and i_cid='"+ccid+"' and i_num="+cnum+";"); 
			stmt.close(); 
	
			System.out.println("adding the 2nd shell ...");  
			sql = "select j_num from temp_shell where shell=1;";
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			i=0; 
			while (rsst.next()) {
				i++; 
				j_num = rsst.getInt(1);
				System.out.println(i+":"+j_num); 
				mst = conn.createStatement();
				sql = "insert into temp_shell select i_num, i_res, j_num, j_res, j_sstype, 2 as shell, weight_bb as BB, weight_SC as SC " +
						"from "+targetEdges+" where graph_id="+cgraph_id+" and i_cid='"+ccid+"' and i_num="+j_num+";";
				// System.out.println(">"+sql); 
				mst.executeUpdate( sql); 
				mst.close(); 
			} // end while 
			rsst.close(); 
			stmt.close(); 

			System.out.println("retrieving the entire nbhood (1st and 2nd shell)");  
			sql = "select j_num, j_res, j_sstype, min(shell) as shell, count(*) as cn, max(BB) as BB, max(SC) as SC " +
					"from temp_shell where j_num!="+cnum+" group by j_num order by j_num;";
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			n1=0; 
			n2=0; 		
			while (rsst.next()) {
				if ( rsst.getInt( 4)==1) { // count 1st shell entry 
					n1++;
					System.out.print("1#"+n1);
				} // end if 2st shell 
				if ( rsst.getInt( 4)==2) { // count 2nd shell entry 
					n2++;
					System.out.print("2#"+n2);
				} // end if 2nd shell 
				System.out.println(" :\t"+rsst.getInt( 1)+"\t"+rsst.getString( 2)+"\t"+rsst.getString( 3)+"\t"+rsst.getInt( 4)+"\t"+rsst.getInt( 5)+"\t"+rsst.getInt( 6)+"\t"+rsst.getInt( 7)); 
			} // end while
			
			System.out.println("|1st shell|="+n1+" \tx\t |2nd shell|="+n2);
			mst = conn.createStatement();
			sql = "delete from target_score where graph_id="+cgraph_id+" and node_id="+cnode_id+" and cid='"+ccid+"' and num="+cnum+";"; 
			System.out.println(">"+sql); 
			mst.executeUpdate( sql); 
			mst.close(); 
			
			cn1 = new int[n1+1]; 
			cn2 = new int[n2+1]; 
			// minus | mres | mss  | mcn | plus | pres | pss  | pcn | nn
			for (j=0; j<=n2; j++) { // outer loop through all indirect contacts
				for (i=0; i<=n1; i++) { // inner loop through all direct contacts
					overx = false;
					if (VL>=1) {
						System.out.print("("+i+","+j+")\t");
					}
					minus = 0; 
					mres  = ""; 
					mss   = ""; 
					mcn   = n1;  
					plus  = 0; 
					pres  = ""; 
					pss   = ""; 
					pcn   = n1;
					nn="%";
					nb="%";
					ni=0; nj=0; 
					rsst.beforeFirst();
					
					while (rsst.next()) {
						j_num = rsst.getInt(1); 
						j_res = rsst.getString(2);
						j_sec = rsst.getString(3); 
						j_shell = rsst.getInt(4);
						j_cnsize = rsst.getInt(5);
						j_bb = rsst.getInt(6);
						j_sc = rsst.getInt(7);
						if (VL>=2) System.out.print("\n"+rsst.getInt( 1)+"\t"+rsst.getString( 2)+"\t"+rsst.getString( 3)+"\t"+rsst.getInt( 4)+"\t"+rsst.getInt( 5)+"\t"+rsst.getInt( 6)+"\t"+rsst.getInt( 7));

						if (j_num>cnum) { // we are over central residue  
							if (!overx) {
								nn+="x%";
								nb+="x%";
								overx=true; 
							} // end if over x 
						} // END IF J > X
						
						if (j_shell==1) { // a direct 1st shell neighbour 
							ni++;
							if (ni!=i) {// if this is NOT the one direct nb 2B dropped
								if (j_sc>=j_bb) { // SC dominated 
									nn+=j_res.toUpperCase()+"%"; // it is included
									nb+=j_res.toUpperCase()+"%";
								} else { // BB dominated 
									nn+=xlateSS( j_sec)+"%"; // it is included
									nb+=xlateSS( j_sec)+"%";
								} // end if SC/BB domin. 
							} else { // this one IS dropped
								minus = j_num; 
								mres  = j_res; 
								mss   = j_sec; 
								mcn   = j_cnsize;  
							} // end if ni!=i 
							
						} else { // 2nd shell neighbour 
							nj++; 
							if (nj==j) { // this is the 2nd shell nb 2B included
								if (j_sc>=j_bb) nn+=j_res.toUpperCase()+"%"; // SC dominated 
								else nn+=xlateSS(j_sec)+"%";
								plus  = j_num; 
								pres  = j_res; 
								pss   = j_sec; 
								pcn   = j_cnsize;  
							} // end if 
			
						} // end if 1st/2nd shell
						
					} // end while through the entire nbhood
					if (!overx) { // in case x is the very last we haven't seen it yet  
						nn+="x%"; // add it in the end
						nb+="x%"; 
						overx=true; 
					} // end if over x
					if (VL>=1) {
						System.out.print("("+nn+")\t("+nb+")\t");
					}
					
					// Store the resulting moves (nn for +SC, nb for +BB contact) 
					// insert into target_score values ( 1, 2, 'C', 0, 'A', 'H', 1, 2, 112, 'V', 'H', 3, 123, 'P', 'O', 2, 'TESTIT', 283, 12, 0, 0.00);
					// SC move into resulting table 
					mst = conn.createStatement();
					sql = "insert into target_score values ( "+cgraph_id+", "+cnode_id+", '"+ccid+"', "+cnum+"," +
							" '"+cres+"', '"+csstype+"', "+i+", "+j+", "+
							minus+", '"+mres+"', '"+mss+"', "+mcn+", "+
							plus +", '"+pres+"', '"+pss+"', "+pcn+", "+
							"'"+nn+"', 0, 0, 0, 0.00);";
							
					// System.out.println(">"+sql); 
					mst.executeUpdate( sql); 
					mst.close();
					
				} // close inner loop (i)
				if (VL>=1) {
				   System.out.println(".");
				} else {
					System.out.print(".");
				} 
			} // next outerloop (j)
			
			rsst.close(); 
			stmt.close(); 
			
		} catch (SQLException e) {
			e.printStackTrace();
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
		} // end try/catch 

	}	// end writeMoves 
	
	public static String xlateSS( String sstype) { // xlate SecSTructure-Symbol 
		String SSymbol = "o";
		sstype=sstype.toUpperCase(); 
		if (sstype.equals("H")) SSymbol="z"; // Helix 
		if (sstype.equals("S")) SSymbol="b"; // beta strand 
		return SSymbol; 
	} // end xlateSS

} // end class 


/*
drop table target_score;
CREATE TABLE `target_score` (
  `graph_id` int(10) unsigned NOT NULL,
  `node_id` int(10) unsigned DEFAULT NULL,
  `cid` varchar(6) COLLATE latin1_general_cs NOT NULL,
  `num` int(10) unsigned NOT NULL,
  `res` char(1) COLLATE latin1_general_cs DEFAULT NULL,
  `sstype` char(1) COLLATE latin1_general_cs DEFAULT NULL,
  `i` int(2) NOT NULL DEFAULT '0',
  `j` int(2) NOT NULL DEFAULT '0',
  `minus` int(2) NOT NULL DEFAULT '0',
  `mres` char(1) COLLATE latin1_general_cs DEFAULT NULL,
  `mss` char(1) COLLATE latin1_general_cs DEFAULT NULL,
  `mcn` int(2) NOT NULL DEFAULT '0',
  `plus` int(2) NOT NULL DEFAULT '0',
  `pres` char(1) COLLATE latin1_general_cs DEFAULT NULL,
  `pss` char(1) COLLATE latin1_general_cs DEFAULT NULL,
  `pcn` int(2) NOT NULL DEFAULT '0',
  `nn` varchar(50) COLLATE latin1_general_cs,
  `total` int(8) NOT NULL DEFAULT '0',
  `rank` int(3) NOT NULL DEFAULT '0',
  `deltarank` int(3) NOT NULL DEFAULT '0',
  `score` decimal(6,2) NOT NULL DEFAULT '0.00'
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_general_cs;


*/ 