import tools.MySQLConnection;

import java.sql.SQLException;
import java.sql.Statement;
import java.sql.ResultSet;

public class listInfoGain {

	/**
	 * "Hello World" for entropy calculations given a nbhoodstring 
	 * 
	 * @author lappe
	 */
	
	static String user = "lappe"	; // change user name!!
	static MySQLConnection conn;
	static double orgFreq, lastFreq, orgAUC, lastAUC; 
	static int orgRank, lastRank; 
	
	public static void main(String[] args) {
		double entropy = 0.0, newentropy=0.0, gain; 
		String nbhood = "", front, middle, tail, newhood="", central ="I";
		if (args.length<2){
			System.err.println("The starting neighborhood-string and central residue type needs to be given .... i.e. %K%D%L%I%x%D%C% I");
			System.exit(1);
		}		
		nbhood = args[0];
		central = args[1]; 
		
		conn = new MySQLConnection("white",user,"nieve","pdb_reps_graph_4_2"); 
		int l = (int)((nbhood.length()-1)/2); 
		int N = 1; 
			System.out.println("ListInfoGain");
			System.out.print("0 - (%x%)("+l+") "); 
			entropy = getEntropy( "%x%", central); 
			System.out.print( entropy + " bits."); 
			// System.out.println("Symbols in nbhood :"+l); 
			orgFreq = lastFreq;
			orgAUC = lastAUC; 
			orgRank = lastRank; 
			System.out.println("central Residue "+central+ " rank#"+lastRank+" p="+lastFreq+"  AUC="+String.format("%.5f",lastAUC));
//			 System.out.println("Symbols in nbhood :"+l);
			
			for (int i = 0; i<l; i++) {
				front = nbhood.substring(0,i*2); 
				middle = nbhood.substring(i*2, i*2+2); 
				tail = nbhood.substring(i*2+2);
				System.out.print((i+1)+" - "+front+"("+middle+")"+tail);
				
				if (middle.equals("%x")) { // switch from N to C of X  
					N = -1; 
					newhood = "%x%"; 
				} else  {
				   if (N < 0) // we are in the C-terminal section 
				       newhood = "%x"+middle+"%"; 
				   else // N terminal (before X) 
					   newhood = middle+"%x%";
				} // end if    
				System.out.print( " -> " + newhood); 
				newentropy = getEntropy( newhood, central); 
				gain = newentropy-entropy; 
				System.out.print( " : "+String.format("%.5f", newentropy) +"bits ("+String.format("%.5f", gain)+")");
				System.out.print( " : #"+lastRank +" ("+(lastRank-orgRank)+")");
				System.out.print( " : "+String.format("%.5f", lastFreq) +"("+String.format("%.5f", (lastFreq-orgFreq))+")");
				System.out.print( " : "+String.format("%.5f", lastAUC) +"("+String.format("%.5f", (lastAUC-orgAUC))+")");
				
				System.out.println( ".");  
			} // next symbol in nbhood;  
		
		System.out.println("fin."); 
	}	

		
	
	public static double getEntropy( String nbs, String centRes) {
		int total = 0, num; 
		String sql, res; 
		Statement stmt;  
		ResultSet rsst;
		double p, psum=0.0, logp, plogp, plogpsum=0.0; 
		
		try {
			sql = "select count(*) from single_model_node where n like '"+nbs+"';";
			// System.out.println( sql); 
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			if (rsst.next()) total = rsst.getInt( 1); 
			rsst.close(); 
			stmt.close(); 
			
			sql = "select res, count(*) as t, count(*)/"+total+" as p from single_model_node where n like '"+nbs+"' group by res order by p DESC;";
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			// System.out.println("rank : res : total t : fraction p : log2(p) : -p*log2(p)");
			int rank = 0; 
			boolean seenCentRes = false;  
			lastAUC= 0; 
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
				}
				if (seenCentRes) lastAUC += p; 
				// System.out.println("");
			}
			// System.out.println("Sum :"+total+"       : "+psum+"       : "+plogpsum);
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
