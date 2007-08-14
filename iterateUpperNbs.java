import tools.MySQLConnection;

import java.sql.SQLException;
import java.sql.Statement;
import java.sql.ResultSet;

public class iterateUpperNbs {

	/**
	 *  
	 * exchange of one neighbor at a time by a common neighbor 
	 * @author lappe
	 */
	static int maxRank = 31; // value to replace for non-existence of central redue in the resultvector (rank=0) 
	// higher values should penalize non-existence more
	static String user = "lappe"	; // change user name!!
	static MySQLConnection conn;
	static double lastEntropy=0.0, lastFreq, lastAUC, lastavgk, lastdevk; 
	static double  orgEntropy=0.0, orgFreq, orgAUC, orgavgk, orgdevk; 
	static int lastRank, lastTotal;  
	static int orgRank, orgTotal;
	
	public static void main(String[] args) throws SQLException {
		
	    if (args.length<2){
			System.err.println("The graph_id and residue-nr. needs to be given .... i.e. 9 28");
			System.exit(1);
		}		
		int graphid = Integer.parseInt( args[0]);
		int resnr = Integer.parseInt( args[1]); 
		int n1=0, n2=0, ni=0, nj=0, j_num=0, j_shell, j_cnsize, i, j, sumdelta, j_num12=0, j_num21=0; 
		conn = new MySQLConnection("white",user,"nieve","pdb_reps_graph_4_2"); // the UPPERCASE DB!  
		String sql, j_res, j_sec, restype="?", ressec="?", nbs, mymove; 
		Statement stmt, jst;  
		ResultSet rsst;
		
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
			
			// remove from here for unperturbed version 
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
			rsst.close(); 
			stmt.close(); 
			
			// create matrices accordingly 
			
			// percolate the environment  
			System.out.println("removing a contact from shell1 -> shell2");  
			sql = "select j_num, j_res, j_sstype, shell from temp_shell where i_num="+resnr+" and shell = 1 order by rand() limit 1;";
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			if (rsst.next()) {
				j_num12 = rsst.getInt( 1); 
				System.out.println(rsst.getInt( 4)+" -> 2 : ("+rsst.getString( 2)+":"+j_num12+":"+rsst.getString( 3)+")");
			}
			rsst.close();
			stmt.close(); 
			stmt = conn.createStatement();
			stmt.executeUpdate("delete from temp_shell where j_num="+j_num12+";"); 
			stmt.close(); 
			
			
	        // and move an indirect nbor from shell2 -> shell 1; 
			System.out.println("adding a contact from shell2 -> shell1");  
			sql = "select j_num, j_res, j_sstype, shell from temp_shell where i_num!="+resnr+" and j_num!="+resnr+" and i_num!="+j_num12+" and j_num!="+j_num12+" and shell = 2 order by rand() limit 1;";
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			if (rsst.next()) {
				j_num21 = rsst.getInt( 1); 
				System.out.println(rsst.getInt( 4)+" -> 1 : ("+rsst.getString( 2)+":"+j_num21+":"+rsst.getString( 3)+")");
			}
			rsst.close();
			stmt.close(); 
			stmt = conn.createStatement();
			stmt.executeUpdate("update temp_shell set i_num="+resnr+", i_res='"+restype+"', j_sstype='"+ressec+"', shell=1 where j_num="+j_num21+";");
			stmt.close(); 
			stmt = conn.createStatement();
			stmt.executeUpdate("delete from temp_shell where shell>1;");
			stmt.close(); 
			n1 = 0; 
			System.out.println("re-building the 2nd shell");  
			sql = "select distinct j_num from temp_shell where shell=1 order by j_num;";
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
			// end of remove for unperturbed version 
			
			System.out.println("retrieving the entire 1st and 2nd shell");  
			sql = "select j_num, j_res, j_sstype, min(shell) as shell, count(*) as cn from temp_shell group by j_num order by j_num;";
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
			 
			
			for (j=0; j<=n2; j++) { // outer loop through all indirect contacts
				// System.out.print(i+" - "); 
				sumdelta=0; 
				for (i=0; i<=n1; i++) { // inner loop through all direct contacts
					mymove = "("+i+","+j+")";
					ni = 0; 
					nj = 0; 
					nbs="%";
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
								nbs+=j_res.toUpperCase()+"%";								
							} else { // this one IS dropped 
								mymove += "(-"+j_res+":"+j_num+":"+j_sec+"/"+j_cnsize+")"; 
							} // end if ni!=i 
						} else { // 2nd shell neighbour 
							nj++; 
							if (j_num==resnr) { // the central residue is part if the 2nd shell
								if (nj!=j) { // drop x if marked for inclusion 
									nbs+="x%"; 
								} else {
									mymove += " no x ..."; 
								}
							} else { // this is not x 
								if (nj==j) { // this is the 2nd shell nb 2B included
									nbs+=j_res.toUpperCase()+"%";
									mymove += "(+"+j_res+":"+j_num+":"+j_sec+"/"+j_cnsize+")"; 
								} // end if 
							} // end if this is central residue x 
						} // end if 1st/2nd shell
					} // end while through the entire nbhood 
					
 					getEntropy( nbs, restype); 
					if (i==0 && j==0) { // original nbhoodstring without any insertions/deletions 
						orgEntropy = lastEntropy;
						orgFreq = lastFreq;
						orgAUC  = lastAUC;
						orgavgk = lastavgk;
						orgdevk = lastdevk; 
						orgRank = lastRank;
						orgTotal= lastTotal; 
					} // end if 0/0 for defining org*
					if (lastRank > 0) { 
						sumdelta += (lastRank-orgRank);
					} else {
						sumdelta += (maxRank-orgRank);
					}
					
					// System.out.println(" t="+orgTotal+" \tentropy = "+String.format("%.5f", orgEntropy)+" \trank#"+String.format("%2d",orgRank)+" \tp("+restype+") = "+String.format("%.5f",orgFreq)+" \t\tAUC = "+String.format("%.5f",orgAUC)+" \t\tavg(k)="+String.format("%.2f",orgavgk)+"\tstddev="+String.format("%.2f",orgdevk));
					//if ((lastRank>0 && lastRank<orgRank) || (i==0 && j==0)) {
						System.out.print(mymove+"\t"+nbs);
						printValues(); 	
					//}
								
					
					// System.out.println(".");
				} // close inner loop (i)
				System.out.println( "Summ "+j+":"+sumdelta);
			} // next outerloop (j)
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
			if (lastRank==0) lastRank = maxRank;  
		} catch (SQLException e) {
			e.printStackTrace();
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
		}
		
	} // end of getEntropy 

} // end class 
