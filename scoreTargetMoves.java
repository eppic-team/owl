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
	static int maxRank = 21; // value to replace for non-existence of central redue in the resultvector (rank=0) 
	// higher values should penalize non-existence more
	static int VL=1; // Verbosity Level 
	static String user = "lappe"	; // change user name!!
	static MySQLConnection conn;
	static String prgID = "V03"; 
	static String backgrndDB = "cullpdb_90";  
	static String targetDB = "CASP_decoys"; 
	static String targetNodes = "target_node";
	static String targetEdges = "target_edge";
	static int sTotal, sRank; 
	
	public static void main(String[] args) {
		
		int graph_id, node_id, num, i, j, total, rank, deltaRank=0, counter=0, nullrank=maxRank, minus, mcn, plus, pcn;
		String sql, scoreTableName, cid, res, sstype, nn, pred="", mres, mss, pres, pss;  
		Statement mstmt;  
		ResultSet mrsst;
				
		try {
			conn = new MySQLConnection("white",user,"nieve", backgrndDB); // connection to the the background DB!  
			System.out.println("Scoring Target neighborhoods v.0.4. "); 
			scoreTableName = targetDB+".score_"+prgID+"_"+backgrndDB; 
			System.out.println("scoreTableName:"+scoreTableName);
			
			// preparing the result db 
			sql = "drop table IF EXISTS "+scoreTableName+";";
			mstmt = conn.createStatement();
			mstmt.executeUpdate(sql); 
			mstmt.close(); 
			//sql = "create table "+scoreTableName+" select * from "+targetDB+".target_score where i=0 and j=0;";
			//sql = "create table "+scoreTableName+" select * from "+targetDB+".target_score where i=0 or j=0;";
			sql = "create table "+scoreTableName+" select * from "+targetDB+".target_score;";
				//	"order by graph_id, node_id, cid, num, i, j limit 23;"; // where pss='';"; to exclude ss inserts  
			mstmt = conn.createStatement();
			mstmt.executeUpdate(sql); 
			mstmt.close(); 

			sql = "select graph_id, node_id, cid, num, res, sstype, i, j, minus, mres, mss, mcn, plus, pres, pss, pcn, nn" +
					" from "+scoreTableName+" order by graph_id, node_id, cid, num, i, j;"; 
			mstmt = conn.createStatement();
			mrsst = mstmt.executeQuery(sql); 
			counter=0; 
			while (mrsst.next()) {
				counter++;
				// System.out.print("\n"+counter+":\t"); 
				
				graph_id = mrsst.getInt( 1);
				node_id  = mrsst.getInt( 2);
				cid      = mrsst.getString( 3);
				num      = mrsst.getInt( 4);
				res      = mrsst.getString( 5);
				sstype   = mrsst.getString( 6);
				i        = mrsst.getInt( 7);
				j        = mrsst.getInt( 8);
				minus    = mrsst.getInt( 9);
				mres     = mrsst.getString( 10);
				mss      = mrsst.getString( 11);
				mcn      = mrsst.getInt( 12);
				plus     = mrsst.getInt( 13);
				pres     = mrsst.getString( 14);
				pss      = mrsst.getString( 15);
				pcn      = mrsst.getInt( 16);
				
				nn       = mrsst.getString( 17);

//				graph_id | node_id | cid | num | res  | sstype | i | j  | minus | mres | mss  | mcn | plus | pres | pss  | pcn | nn | total | rank | deltarank | score
				System.out.print(graph_id+"\t"+node_id+"\t"+cid+"\t"+num+"\t"+res+"\t"+sstype+"\t"+i+"\t"+j+"\t"); 
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
				System.out.println("\t"+total+"\t"+rank+"\t"+deltaRank+"\t"+0.0); 
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

		} catch (SQLException e) {
			e.printStackTrace();
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
		} // end try/catch 
		System.out.println("fin."); 
	}	// end main 


	public static void getCountRank( String nbs, String centRes, String predec) {
		String sql, res; 
		Statement stmt;  
		ResultSet rsst;
		int counter=0, c=0, lastc=0, rank = 0; 
		boolean seenCentRes = false; 
		
		try {
			// Hashing first row tables comes first
//			 Hashing first row tables comes first
			if (VL>=2) {
				System.out.println("getCountRank for ");
				System.out.println("nbs    : "+nbs); 
				System.out.println("centRes: "+centRes);
				System.out.println("predec : "+predec);
			}
			String this_n = nbs.replace("%",""); 
			String prec_n = predec.replace("%","");
			// Check that this is really a superstring of predec ?! 
			if (VL>=2) System.out.println("this_n: ["+this_n+"]"); 
			if (VL>=2) System.out.println("prec_n: ["+prec_n+"]");
			if (prec_n.equals(this_n)) {
				if (VL>=2) System.out.println("have to create db for this "+prec_n);
				sql = "create table IF NOT EXISTS nbhashtables."+prec_n+" as select res, n, k from "+backgrndDB+".nbstrings where n like '"+nbs+"';"; 
				if (VL>=2) System.out.println(" >> "+sql); 
				stmt = conn.createStatement();
				stmt.executeUpdate( sql); 
				stmt.close(); 
			} else if (VL>=2) System.out.println("using preceding db of "+prec_n); 
			
			if (VL>=2) {
				System.out.println("getCountRank for ");
				System.out.println("nbs    : "+nbs); 
				System.out.println("centRes: "+centRes);
			}
			
			// sql = "select count(*) from "+backgrndDB+".single_model_node where n like '"+nbs+"';";
			sql = "select count(*) from nbhashtables."+prec_n+" where n like '"+nbs+"';";
			if (VL>=2) System.out.println( sql); 
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			if (rsst.next()) sTotal = rsst.getInt( 1);
			rsst.close(); 
			stmt.close(); 
			if (VL>=2) System.out.println( "Total : "+sTotal); 
			
			if (sTotal>0) {
				// sql = "select res, count(*) as t from "+backgrndDB+".single_model_node where n like '"+nbs+"' group by res order by t DESC;";
				sql = "select res, count(*) as t from nbhashtables."+prec_n+" where n like '"+nbs+"' group by res order by t DESC;";
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
		} catch (SQLException e) {
			e.printStackTrace();
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
		}
	} // end of getCountRank 

} // end class 
