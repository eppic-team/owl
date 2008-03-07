import tools.MySQLConnection;

import java.sql.SQLException;
import java.sql.Statement;
import java.sql.ResultSet;

public class scoreTargetMoves {

	/**
	 * 
	 * calculate count, score etc. for all targetMoves 
	 * 
	 * @author lappe
	 */
	private static int maxRank = 21; // value to replace for non-existence of central redue in the resultvector (rank=0) 
	// higher values should penalize non-existence more
	private static int VL=1; // Verbosity Level 
	
	private static final String dbserver = "white";
	private static final String user = MySQLConnection.getUserName();
	private static final String pwd = "nieve";
	private static MySQLConnection conn;
	 
	private static String backgrndDB = "cullpdb_90";   // change!!
	private static String hashDB = "nbhashing"; 
	
	static int sTotal, sRank; 
	
	public static void main(String[] args) throws SQLException{
		
		if (args.length<2) {
			System.err.println("Must pass 2 parameters: <graph_id>, <target_db.target_moves_table>");
			System.exit(1);
		} 
		int graph_id = Integer.parseInt(args[0]);
		String targetMovesTable = args[1];
		
		int num, i, j, total, rank, deltaRank=0, counter=0, nullrank=maxRank, minus, mcn, plus, pcn;
		String sql, cid, res, sstype, nn, pred="", mres, mss, pres, pss;  
		Statement mstmt;  
		ResultSet mrsst;
				
		conn = new MySQLConnection(dbserver,user,pwd, backgrndDB); // connection to the the background DB!  
		//System.out.println("Scoring Target neighborhoods "+prgID); 
 

		// preparing the result db 
		//sql = "drop table IF EXISTS "+scoreTableName+";";
		//mstmt = conn.createStatement();
		//mstmt.executeUpdate(sql); 
		//mstmt.close(); 
		//sql = "create table "+scoreTableName+" select * from "+targetScore;  
		//mstmt = conn.createStatement();
		//mstmt.executeUpdate(sql); 
		//mstmt.close(); 

		sql = "SELECT graph_id, cid, num, res, sstype, i, j, minus, mres, mss, mcn, plus, pres, pss, pcn, nn" +
		" FROM "+targetMovesTable+
		" WHERE graph_id="+graph_id+
		" ORDER BY graph_id, cid, num, i, j ;";

		mstmt = conn.createStatement();
		mrsst = mstmt.executeQuery(sql); 
		counter=0; 
		while (mrsst.next()) {
			counter++;
			// System.out.print("\n"+counter+":\t"); 

			//graph_id = mrsst.getInt( 1);
			cid      = mrsst.getString(2);
			num      = mrsst.getInt(3);
			res      = mrsst.getString(4);
			sstype   = mrsst.getString(5);
			i        = mrsst.getInt(6);
			j        = mrsst.getInt(7);
			minus    = mrsst.getInt(8);
			mres     = mrsst.getString(9);
			mss      = mrsst.getString(10);
			mcn      = mrsst.getInt(11);
			plus     = mrsst.getInt(12);
			pres     = mrsst.getString(13);
			pss      = mrsst.getString(14);
			pcn      = mrsst.getInt(15);
			nn       = mrsst.getString(16);

//			graph_id | cid | num | res  | sstype | i | j  | minus | mres | mss  | mcn | plus | pres | pss  | pcn | nn | total | rank | deltarank | score
			System.out.print(graph_id+"\t"+cid+"\t"+num+"\t"+res+"\t"+sstype+"\t"+i+"\t"+j+"\t"); 
			System.out.print(minus+"\t"+mres+"\t"+mss+"\t"+mcn+"\t"+plus+"\t"+pres+"\t"+pss+"\t"+pcn+"\t"+nn+"\t");

			if (j==0) { // top entry / column of movematrix 
				pred=nn; 
			}
			if (VL>=2) System.out.println("\n["+i+","+j+"]");
			getCountRank( nn, res, pred); 
			total = sTotal; 
			rank = sRank;
			if (i==0 && j==0) nullrank = rank; 
			deltaRank = rank-nullrank; 
			System.out.println(total+"\t"+rank+"\t"+deltaRank+"\t"+(pcn*mcn*deltaRank)); 
			// graph_id | node_id | cid | num | res  | sstype | i | j  | minus | mres | mss  | mcn | plus | pres | pss  | pcn | nn | total | rank | deltarank | score

			/*sql = "update "+scoreTableName+" set total="+sTotal+", rank="+sRank+", deltarank="+deltaRank
				+" where graph_id="+graph_id+" and node_id="+node_id+" and cid='"+cid+"' and num="+num
				+" and i="+i+" and j="+j+";";
				nstmt = conn.createStatement();
				nstmt.executeUpdate(sql); 
				nstmt.close(); 
			 */ 
		} // end while 

		// Cleanup ... 
		mrsst.close(); 
		mstmt.close(); 

		//System.out.println("fin."); 
	}	// end main 


	private static void getCountRank( String nbs, String centRes, String predec) throws SQLException {
		String sql, res; 
		Statement stmt;  
		ResultSet rsst;
		int counter=0, c=0, lastc=0, rank = 0; 
		boolean seenCentRes = false; 
		
		if (VL>=2) {
			System.out.println("getCountRank for ");
			System.out.println("nbs    : "+nbs); 
			System.out.println("centRes: "+centRes);
		}
		String hashKey = ""; 
		String hashTable = "";
		if (nbs.replace("%", "").length()>2) { // case 1: neighbourhood at least length 2
			hashKey = nbs.replace("%","").substring(0, 2);
			hashTable = hashDB+"."+backgrndDB+"_"+hashKey;
			sql = "select sum(c) from "+hashTable+" where n like '"+nbs+"';";
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			if (rsst.next()) sTotal = rsst.getInt( 1);
			rsst.close(); 
			stmt.close(); 
		} else { // case 2: length on neighbourhood is 1: we just have 'x', we query whole table
			sql = "select sum(c) from "+backgrndDB+".nb_equals;";
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			if (rsst.next()) sTotal = rsst.getInt( 1);
			rsst.close(); 
			stmt.close(); 				
		}

		if (VL>=2) System.out.println( "Total : "+sTotal); 

		if (sTotal>0) {
			if (!hashTable.equals("")) { // i.e. we set hashTable (case 1 above) 
				sql = "select res, sum(c) as t from "+hashTable+" where n like '"+nbs+"' group by res order by t DESC;";
			} else { // i.e. we didn't set hashTable (case 2 above)
				sql = "select res, sum(c) as t from "+backgrndDB+".nb_equals group by res order by t DESC;";
			}
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			if (VL>=2) System.out.println("###\tres\ttotal t");
			sRank = maxRank; 
			while (rsst.next()) {	
				counter++; 
				res = rsst.getString(1); // 1st column -- res
				c = rsst.getInt( 2); // 2nds column : count/residue
				if (VL>=2) System.out.print(counter+"\t"+res+"\t"+c);

				if ((c == lastc) && (lastc>0) && seenCentRes) { // tie 
					if (VL>=2) System.out.print(" <-- TIE!");
					rank = counter; 		
				} // end if 
				if (res.equals(centRes)) { 
					if (VL>=2) System.out.print(" <== " + centRes);
					seenCentRes = true;
					rank = counter; 
				}
				if (VL>=2) System.out.println(".");
			}
			if (VL>=2) System.out.println("=> rank "+rank); 
			rsst.close(); 
			stmt.close(); 
		} // endif sTotal > 0 
		if (rank==0) sRank = maxRank;
		else sRank = rank; 
	} // end of getCountRank 

} // end class 
