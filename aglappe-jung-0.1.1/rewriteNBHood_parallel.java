import tools.MySQLConnection;

import java.sql.SQLException;
import java.sql.Statement;
import java.sql.ResultSet;

public class rewriteNBHood_parallel {

	/**
	 * reriting the nbhoodstring such that BB dominated nbs are represented by SSType 
	 * 
	 * @author lappe
	 */

	static String user = "lappe"	; // change user name!!
	static MySQLConnection conn;
	static int k=0, k_sc=0, k_bb=0; 
	static String process,  cid, res, sectype, fullnbs;
	static int graph_id=0, node_id=0, resnr=0;


	public static void main(String[] args) throws SQLException {
		String nbhood = "", isql; 
		Statement istmt; 

		System.out.println("rewriteNBHood");
		if (args.length!=1){
			System.err.println("The processname has to be given .... i.e. tla01");
			System.exit(1);
		}		
		process = args[0];
		System.out.println("processID="+process);

		try {

			conn = new MySQLConnection("white",user,"nieve","pdb_reps_graph_4_2"); 
			node_id = getNextNode(); 
			while (node_id>0) {
				System.out.print(graph_id); // graph
				System.out.print("\t"+node_id); // node 
				System.out.print("\t"+cid); // cid
				System.out.print("\t"+resnr); // num
				System.out.print("\t"+res); // res
				System.out.print("\t"+sectype); // sstype 
				fullnbs = "";
				nbhood= getSSNbHoodString( graph_id, node_id, resnr); 

				System.out.println("\t"+nbhood+"\t"+fullnbs);
				istmt = conn.createStatement();
				isql = "update single_model_cofull set k="+k+", k_bb="+k_bb+", k_sc="+k_sc+", n=\'"+nbhood+"\', nf=\'"+fullnbs+"\'" +
						" where ( graph_id="+graph_id+" and node_id="+node_id+");";
				// System.out.println("i>"+isql); 
				istmt.executeUpdate( isql);

				node_id = getNextNode(); 
			} // end while we have another node 

		} catch (SQLException e) {
			e.printStackTrace();
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
		}

		System.out.println("fin."); 
	}	

	public static int getNextNode() {
		String sql; 
		Statement stmt;  
		ResultSet rsst;

		graph_id=0;  
		node_id =0; 
		try {
			sql = "update single_model_cofull set n=\'>"+process+"\' where length(n)=0 limit 1;";
			stmt = conn.createStatement();
			stmt.executeUpdate(sql);
			sql = "select graph_id, node_id, cid, num, res, sstype from single_model_cofull where n=\'>"+process+"\';";
			// System.out.println("node> "+ sql); 
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);

			if (rsst.next()) {	
				// System.out.print(">"); 
				graph_id=rsst.getInt(1);  
				node_id =rsst.getInt(2); 
				cid     =rsst.getString(3); 
				resnr   =rsst.getInt(4);
				res     =rsst.getString(5);
				sectype =rsst.getString(6); 
			}
		} catch (SQLException e) {
			e.printStackTrace();
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
		}
		return node_id; 
	} // end of getnbhoodstring 


	public static String getSSNbHoodString( int gid, int nid, int num ) {
		int j_num, j_bb, j_sc;
		boolean overx=false; 
		String sql, nbs="", j_res, j_ss, nextnb; 
		Statement stmt;  
		ResultSet rsst;

		try {
			k=0; k_sc=0; k_bb=0; 
			sql = "select j_num, j_res, j_sstype, BB, SC from chain_edge where graph_id="+gid+" and i_node_id="+nid+" order by j_num;";
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			fullnbs=""; 
			while (rsst.next()) {	
				j_num = rsst.getInt(1); 
				j_res = rsst.getString(2); 				
				j_ss = rsst.getString(3); 
				j_bb = rsst.getInt(4); 
				j_sc = rsst.getInt(5);
				nextnb="?"; 
				if (!overx && j_num>num) { 
					overx=true; 
					nbs+="x";
					fullnbs+="x"; 
				}
				if ( (j_bb>j_sc) ) { // bb dominated -> only symbol for sstype
					k_bb++; 
					nextnb="o"; // for "other" 
					if (j_ss.equals("H")) nextnb="z"; // Helix 
					if (j_ss.equals("S")) nextnb="b"; // beta strand 
				} else { // sc dominated contact
					k_sc++; 
					nextnb=j_res; 
				} // end if bb dominated
				nbs+= nextnb;
				fullnbs+= j_res; 
				k++; 
				// System.out.print( nextnb);
			}
			rsst.close(); 
			stmt.close();
			if (!overx) {
				nbs+="x"; // add if last in string == not seen yet
				fullnbs+="x"; 
			}

		} catch (SQLException e) {
			e.printStackTrace();
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
		}
		return nbs; 
	} // end of getnbhoodstring 
}
