import tools.MySQLConnection;

import java.sql.SQLException;
import java.sql.Statement;
import java.sql.ResultSet;

public class listInfoLoss {

	/**
	 * info loss / entropy calculations given a nbhoodstring 
	 * dropping one neighbor at a time 
	 * @author lappe
	 */
	
	static String user = "lappe"	; // change user name!!
	static MySQLConnection conn;
	static double orgFreq, lastFreq, orgAUC, lastAUC; 
	static int orgRank, lastRank, orgTotal, lastTotal; 
	
	public static void main(String[] args) throws SQLException {
		double entropy = 0.0, newentropy=0.0, gain; 
		String nbhood = "", front, middle, tail, newhood="", central ="I";
		if (args.length<2){
			System.err.println("The starting neighborhood-string and central residue type needs to be given .... i.e. %K%D%L%I%x%D%C% I");
			System.exit(1);
		}		
		nbhood = args[0];
		central = args[1]; 
		
		conn = new MySQLConnection("white",user,"nieve","pdb_reps_graph_4_2_a"); 
		int l = (int)((nbhood.length()-1)/2); 
		
			System.out.println("ListInfoLoss");
			System.out.print("0 - ("+nbhood+")("+l+") "); 
			entropy = getEntropy( nbhood, central); 
			
			// System.out.println("Symbols in nbhood :"+l); 
			orgFreq = lastFreq;
			orgAUC = lastAUC; 
			orgRank = lastRank; 
			orgTotal = lastTotal; 
			System.out.println("central "+central+". \tt="+orgTotal+" \tentropy = "+String.format("%.5f", entropy)+" \trank#"+lastRank+" \tp("+central+") = "+String.format("%.5f",lastFreq)+" \t\tAUC = "+String.format("%.5f",lastAUC));
//			 System.out.println("Symbols in nbhood :"+l);
			
			for (int i = 0; i<l; i++) {
				front = nbhood.substring(0,i*2); 
				middle = nbhood.substring(i*2, i*2+2); 
				tail = nbhood.substring(i*2+2);
				System.out.print((i+1)+" - "+front+"("+middle+")"+tail);
				
				newhood = front+tail; 
				
				System.out.print( " -> " + newhood); 
				newentropy = getEntropy( newhood, central); 
				gain = newentropy-entropy; 
				System.out.print( "\t"+lastTotal+"("+(lastTotal-orgTotal)+")"); 
				System.out.print( "\t"+String.format("%.5f", newentropy) +" ("+String.format("%.5f", gain)+")");
				System.out.print( "\t#"+lastRank +" ("+(lastRank-orgRank)+")");
				System.out.print( "\t"+String.format("%.5f", lastFreq) +" ("+String.format("%.5f", (lastFreq-orgFreq))+")");
				System.out.print( "\t"+String.format("%.5f", lastAUC) +" ("+String.format("%.5f", (lastAUC-orgAUC))+")");
				
				System.out.println( ".");  
			} // next symbol in nbhood;  
		
		System.out.println("fin."); 
	}	

		
	
	public static double getEntropy( String nbs, String centRes) { 
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
			
			sql = "select res, count(*) as t, count(*)/"+lastTotal+" as p from single_model_node where n like '"+nbs+"' group by res order by p DESC;";
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			// System.out.println("rank : res : total t : fraction p : log2(p) : -p*log2(p)");
			int rank = 0; 
			boolean seenCentRes = false;  
			lastAUC= 0; 
			while (rsst.next()) {	
				rank ++;
				res = rsst.getString(1); // 1st column -- res
				//num = rsst.getInt(2); // 2nd column -- num
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
				}
				if (seenCentRes) lastAUC += p; 
				// System.out.println("");
			}
			// System.out.println("Sum :"+lastTotal+"       : "+psum+"       : "+plogpsum);
			rsst.close(); 
			stmt.close(); 
			
		} catch (SQLException e) {
			e.printStackTrace();
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
		}
		return plogpsum; 
	} // end of getEntropy 

}
