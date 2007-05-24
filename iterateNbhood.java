import tools.MySQLConnection;

import java.sql.SQLException;
import java.sql.Statement;
import java.sql.ResultSet;

public class iterateNbhood {

	/**
	 *  
	 * exchange of one neighbor at a time by a common neighbor 
	 * @author lappe
	 */
	
	static String user = "lappe"	; // change user name!!
	static MySQLConnection conn;
	
	
	public static void main(String[] args) {
		
	    if (args.length<2){
			System.err.println("The graph_id and residue-nr. needs to be given .... i.e. 9 28");
			System.exit(1);
		}		
		int graphid = Integer.parseInt( args[0]);
		int resnr = Integer.parseInt( args[1]); 
		int n1=0, n2=0, ni=0, nj=0, j_num=0, lastNum=0, i, j; 
		conn = new MySQLConnection("white",user,"nieve","pdb_reps_graph_4_2_a"); 
		String sql, j_res; 
		Statement stmt, jst;  
		ResultSet rsst, jrs;
		
		try {
			System.out.println("getting direct neighborhood ... "); 
			stmt = conn.createStatement();
			stmt.executeUpdate("drop table if exists temp_shell;"); 
			stmt.close(); 
	
			stmt = conn.createStatement();
			stmt.executeUpdate("create table temp_shell as select i_num, i_res, j_num, j_res, j_sstype, 1 as shell from single_model_edge where graph_id="+graphid+" and i_num="+resnr+";");
			stmt.close(); 
	
			System.out.println("building the 2nd shell");  
			sql = "select j_num, j_res, j_sstype from temp_shell where shell=1;";
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			while (rsst.next()) {
				n1++; 
				j_num = rsst.getInt(1);
				System.out.println(n1+":"+j_num); 
				jst = conn.createStatement();
				sql = "insert into temp_shell select i_num, i_res, j_num, j_res, j_sstype, 2 as shell from single_model_edge where graph_id="+graphid+" and i_num="+j_num+";";
				// System.out.println(">"+sql); 
				jst.executeUpdate( sql); 
				jst.close(); 
			} // end while 
			rsst.close(); 
			stmt.close(); 
			
			System.out.println("geting the entire 2st and 2nd shell");  
			sql = "select j_num, j_res, j_sstype, min(shell) as shell, count(*) as cn from temp_shell group by j_num;";
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			// counting shell2
			n2=0; 
			while (rsst.next()) {
				if ( rsst.getInt( 4)==2) n2++;
				System.out.println(n2+":"+rsst.getInt( 1)+"\t"+rsst.getString( 2)+"\t"+rsst.getString( 3)+"\t"+rsst.getInt( 4)+"\t"+rsst.getInt( 5)); 
			} // end while 
			System.out.println("SIZE 1st shell "+n1); 
			System.out.println("SIZE 2nd shell "+n2);
			
			for (i=0; i<=n1; i++) { // 1st loop through all direct contacts
				System.out.print(i+" - "); 
				for (j=1; j<=n2; j++) { // 2nd loop through all indirect contacts
					System.out.print(j+":");
					
					// if j_num == resnr -> then x 
				/*
				lastNum=0;
				nj = 0; 
				rsst.beforeFirst();
				while (rsst.next()) {
					nj++; 
					j_num = rsst.getInt( 1); 
					j_res = rsst.getString(2);
					if (nj!=i) { // this is the neighbor 2 B dropped 
						if (lastNum<resnr && resnr<j_num) {
							System.out.print("%"+nj+"x"); 
						}
						System.out.print("%"+j_res);
					} 
					lastNum = j_num; 
				} // end while
				*/  
					
			 
				} // close for loop 2 (j)
				System.out.println(".");
			} // next loop 1 (i)
			rsst.close(); 
			stmt.close(); 
			
			// Summarize 
			
			
			
		} catch (SQLException e) {
			e.printStackTrace();
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
		} // end try/catch 
		System.out.println("fin."); 
	}	

		

}
