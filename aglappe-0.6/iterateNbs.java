import tools.MySQLConnection;

import java.sql.SQLException;
import java.sql.Statement;
import java.sql.ResultSet;

public class iterateNbs {

	/**
	 *  
	 * exchange of one neighbor at a time by a common neighbor 
	 * @author lappe
	 */
	
	static String user = "lappe"	; // change user name!!
	static MySQLConnection conn;
	static double lastEntropy=0.0, lastFreq, lastAUC, lastavgk, lastdevk; 
	static double  orgEntropy=0.0, orgFreq, orgAUC, orgavgk, orgdevk; 
	static int lastRank, lastTotal;  
	static int orgRank, orgTotal;
	
	public static void main(String[] args) {
		
	    if (args.length<2){
			System.err.println("The graph_id and residue-nr. needs to be given .... i.e. 9 28");
			System.exit(1);
		}		
		int graphid = Integer.parseInt( args[0]);
		int resnr = Integer.parseInt( args[1]); 
		int n1=0, n2=0, ni=0, nj=0, j_num=0, j_shell, j_cnsize, i, j; 
		conn = new MySQLConnection("white",user,"nieve","pdb_reps_graph_4_2_a"); 
		String sql, j_res, j_sec, restype="?", ressec="?", nbs_lo, nbs_up; 
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
			
			System.out.println("retrieving the entire 1st and 2nd shell");  
			sql = "select j_num, j_res, j_sstype, min(shell) as shell, count(*) as cn from temp_shell group by j_num;";
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			// counting shell2
			n2=0; 
			while (rsst.next()) {
				if ( rsst.getInt( 4)==2) { // count 2nd shell entry 
					n2++;
					if ( rsst.getInt( 1)==resnr) { // this is the central node -> get type and secondary structure 
						restype =  rsst.getString( 2).toUpperCase(); 
						ressec  =  rsst.getString( 3).toUpperCase();
					} // end if central residue 
				} // end if 2nds shell 
				System.out.println(n2+":"+rsst.getInt( 1)+"\t"+rsst.getString( 2)+"\t"+rsst.getString( 3)+"\t"+rsst.getInt( 4)+"\t"+rsst.getInt( 5)); 
			} // end while 
			System.out.println("SIZE 1st shell "+n1); 
			System.out.println("SIZE 2nd shell "+n2);
			System.out.println("GraphID "+graphid+" Central residue is "+restype+":"+resnr+":"+ressec); 
			
			for (i=0; i<=n1; i++) { // 1st loop through all direct contacts
				// System.out.print(i+" - "); 
				for (j=0; j<=n2; j++) { // 2nd loop through all indirect contacts
					System.out.print("("+i+","+j+")");
					ni = 0; 
					nj = 0; 
					nbs_lo="%";
					nbs_up="%";
					rsst.beforeFirst();
					while (rsst.next()) {
						j_num = rsst.getInt( 1); 
						j_res = rsst.getString(2);
						j_sec = rsst.getString(3);
						j_shell = rsst.getInt( 4);
						j_cnsize = rsst.getInt( 5);
					  
						if (j_shell==1) { // a direct 1st shell neighbour 
							ni++;
							if (ni!=i) {// if this is NOT the one direct nb 2B dropped 
								nbs_lo+=j_res+"%";
								nbs_up+=j_res+"%";								
							} else { // this one IS dropped 
								System.out.print("(-"+j_res+":"+j_num+":"+j_sec+"/"+j_cnsize+")"); 
							} // end if ni!=i 
						} else { // 2nd shell neighbour 
							nj++; 
							if (j_num==resnr) { // the central residue is part if the 2nd shell
								if (nj!=j) { // drop x if marked for inclusion 
									nbs_lo +="x%"; 
									nbs_up +="x%"; 
								} else {
									System.out.print(" no x ..."); 
								}
							} else { // this is not x 
								if (nj==j) { // this is the 2nd shell nb 2B included
									nbs_lo+=j_res.toLowerCase()+"%";
									nbs_up+=j_res.toUpperCase()+"%";
									System.out.print("(+"+j_res+":"+j_num+":"+j_sec+"/"+j_cnsize+")"); 
								} // end if 
							} // end if this is central residue x 
						} // end if 1st/2nd shell
					} // end while through the entire nbhood 
					
 					getEntropy( nbs_lo, restype); 
					if (i==0 && j==0) { // original nbhoodstring without any insertions/deletions 
						orgEntropy = lastEntropy;
						orgFreq = lastFreq;
						orgAUC  = lastAUC;
						orgavgk = lastavgk;
						orgdevk = lastdevk; 
						orgRank = lastRank;
						orgTotal= lastTotal; 
					} // end if 0/0 for defining org* 
					System.out.println(" t="+orgTotal+" \tentropy = "+String.format("%.5f", orgEntropy)+" \trank#"+String.format("%2d",orgRank)+" \tp("+restype+") = "+String.format("%.5f",orgFreq)+" \t\tAUC = "+String.format("%.5f",orgAUC)+" \t\tavg(k)="+String.format("%.2f",orgavgk)+"\tstddev="+String.format("%.2f",orgdevk));
					System.out.print("\tlo "+nbs_lo);
					printValues(); 				
					System.out.print("\tup "+nbs_up);
					getEntropy(nbs_up, restype); 
					printValues(); 
					
					// System.out.println(".");
				} // close for loop 2 (j)
			} // next loop 1 (i)
			rsst.close(); 
			stmt.close(); 
			
		} catch (SQLException e) {
			e.printStackTrace();
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
		} // end try/catch 
		System.out.println("fin."); 
	}	// end main 
	
	public static void printValues() {
		System.out.print( "\t"+lastTotal+"("+(lastTotal-orgTotal)+")"); 
		System.out.print( "\t"+String.format("%.5f", lastEntropy) +" ("+String.format("%.5f", lastEntropy-orgEntropy)+")");
		System.out.print( "\t#"+String.format("%2d",lastRank) +" ("+String.format("%2d",(lastRank-orgRank))+")");
		System.out.print( "\t"+String.format("%.5f", lastFreq) +" ("+String.format("%.5f", (lastFreq-orgFreq))+")");
		System.out.print( "\t"+String.format("%.5f", lastAUC) +" ("+String.format("%.5f", (lastAUC-orgAUC))+")");
		System.out.print( "\t"+String.format("%.2f", lastavgk) +" ("+String.format("%.2f", (lastavgk-orgavgk))+")");
		System.out.print( "\t"+String.format("%.2f", lastdevk) +" ("+String.format("%.2f", (lastdevk-orgdevk))+")");
		System.out.println("");
	}
	
	public static void getEntropy( String nbs, String centRes) {
		String sql, res; 
		Statement stmt;  
		ResultSet rsst;
		double p, psum=0.0, logp, plogp, plogpsum=0.0; 
		int num; 
		try {
			sql = "select count(*) from single_model_node where n like '"+nbs+"';";
			// System.out.println( sql); 
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			if (rsst.next()) lastTotal = rsst.getInt( 1); 
			rsst.close(); 
			stmt.close(); 
			
			sql = "select res, count(*) as t, count(*)/"+lastTotal+" as p, avg( k), stddev( k) from single_model_node where n like '"+nbs+"' group by res order by p DESC;";
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			// System.out.println("rank : res : total t : fraction p : log2(p) : -p*log2(p)");
			int rank = 0; 
			boolean seenCentRes = false;  
			lastAUC = 0.0; 
			lastRank = 0; 
			lastFreq = 0.0;
			lastavgk = 0.0;
			lastdevk = 0.0;
			while (rsst.next()) {	
				rank ++;
				res = rsst.getString(1); // 1st column -- res
				num = rsst.getInt(2); // 2nd column -- num
				p = rsst.getDouble(3); // 3rd: fraction p 
				// System.out.print(rank+ " : " + res+"   : "+num+ " : " + p);
				logp = Math.log(p)/Math.log(2.0); // to basis 2 for info in bits 
				// System.out.print(" : " + logp); 
				plogp = -1.0 * p * logp;  
				// System.out.print(" : " + plogp);
				plogpsum += plogp;  
				psum += p; 
				
				if (res.equals(centRes)) { 
					// System.out.print(" <==" + centRes);
					seenCentRes = true;
					lastFreq = p;
					lastRank = rank; 
					lastavgk = rsst.getDouble(4);
					lastdevk = rsst.getDouble(5);
				}
				if (seenCentRes) lastAUC += p; 
				// System.out.println("");
			}
			// System.out.println("Sum :"+lastTotal+"       : "+psum+"       : "+plogpsum);
			rsst.close(); 
			stmt.close(); 
			lastEntropy = plogpsum; 
			
		} catch (SQLException e) {
			e.printStackTrace();
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
		}
		
	} // end of getEntropy 

} // end class 
