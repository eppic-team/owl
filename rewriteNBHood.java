import tools.MySQLConnection;

import java.sql.SQLException;
import java.sql.Statement;
import java.sql.ResultSet;

public class rewriteNBHood {

	/**
	 * rewriting the nbhoodstring such that BB dominated nbs are represented by SSType
	 * only single gapsymbol is used for nbhood and fullstring _With _Gaps (nwg and fwg)  
	 * output to stdout such that it can be read by
	 * mysql> load data infile '/scratch/local/cullpdb/cullpdb_20.nbhoods' into table nbstrings;
	 * see bottom for table definition 
	 * @author lappe
	 */

	static String user = "lappe", fullnbs="", gapnbs="", nums="", gaphood=""; // change user name!!
	static MySQLConnection conn;
	static int k=0, k_sc=0, k_bb=0;
	static int minDmeansSC=10; // above this threshold contacts are assumed to be SC dominated 
	static String backgrndDB = "cullpdb_40"; 

	public static void main(String[] args) throws SQLException {
		String nbhood = "", msql, cid, res, sectype;
		int graph_id=0, node_id=0, resnr=0; 
		// System.out.println("rewriteNBHood v.0.2");
		conn = new MySQLConnection("white",user,"nieve", backgrndDB);
		Statement mstmt; 
		ResultSet mrsst;

		try {
			msql = "select graph_id, node_id, cid, num, res, sstype from single_model_node where sstype is not null;";
			// System.out.println("nodes> "+ msql);
			mstmt = conn.createStatement();
			mrsst = mstmt.executeQuery(msql);

			while (mrsst.next()) {
				// System.out.print(">");
				graph_id=mrsst.getInt(1);
				node_id =mrsst.getInt(2);
				cid     =mrsst.getString(3);
				resnr   =mrsst.getInt(4);
				res     =mrsst.getString(5);
				sectype =mrsst.getString(6);
				System.out.print(graph_id); // graph
				System.out.print("\t"+node_id); // node
				System.out.print("\t"+cid); // cid
				System.out.print("\t"+resnr); // num
				System.out.print("\t"+res); // res
				System.out.print("\t"+sectype); // sstype

				nbhood= getSSNbHoodString( graph_id, cid, resnr);
				System.out.print("\t"+k);
				System.out.print("\t"+k_bb);
				System.out.print("\t"+k_sc);  
				System.out.print("\t"+nbhood);
				System.out.print("\t"+gaphood);
				System.out.print("\t"+fullnbs);
				System.out.print("\t"+gapnbs);
				// System.out.print("\t"+nums);
				System.out.println(); 

			}
			mrsst.close();
			mstmt.close();

		} catch (SQLException e) {
			e.printStackTrace();
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
		}

		// System.out.println("fin.");
	}

	public static String getSSNbHoodString( int gid, String cid, int num ) {
		int j_num, j_bb, j_sc, lastnum=0;
		boolean overx=false; 
		String sql, nbs="", j_res, j_ss, nextnb, gapsym="", xgapsym=""; 
		Statement stmt;  
		ResultSet rsst;

		try {
			k=0; k_sc=0; k_bb=0; // initialising the degreecounts with 0  
			sql = "select j_num, j_res, j_sstype, weight_BB, weight_SC from single_model_edge where graph_id="+gid+" and i_cid='"+cid+"' and i_num="+num+" order by j_num;";
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			fullnbs=""; 
			gapnbs="";
			gaphood="";
			nums=""; 
			while (rsst.next()) {	
				gapsym=""; 
				xgapsym=""; 
				j_num = rsst.getInt(1); 
				j_res = rsst.getString(2); 				
				j_ss = rsst.getString(3); 
				j_bb = rsst.getInt(4); 
				j_sc = rsst.getInt(5);
				// System.out.print("\n-> "+j_res+" "+j_num+" "+j_ss+" "+"("+j_sc+"/"+j_bb+")");
				nextnb="?"; 
				if (lastnum>0 && ((j_num-lastnum)>=2)) gapsym="%";   
				if (!overx && j_num>num) { 
					overx=true; 
					nbs+="x";
					fullnbs+="x";
					nums+=(" ("+num+") "); 
					if (lastnum>0 && ((num-lastnum)>2)) xgapsym="%";   
					gapnbs+=xgapsym+"x"; 
					gaphood+=xgapsym+"x"; 
					gapsym=""; 
					if ((j_num-num)>=2) gapsym="%";  
				}
				// conversion of BB dominated contacts
				if (j_sc==0 && j_bb==0) { // both values are 0 -> decide by sequence sparation
					// System.out.print("SC&BB=0"); 
					if ( Math.abs( num-j_num)>=minDmeansSC) {
						nextnb=j_res;
						k_sc++; 
						// System.out.print(" -> SC ");
					} else {
						nextnb=xlateSS( j_ss);
						k_bb++; 
						// System.out.print(" -> BB ");
					} // end if Math.abs( num-j_num)>=minDmeansSC
				} else { // at least 1 value is non zero 
					// System.out.print("SC>0 or BB>0");
					if ( (j_bb>j_sc) ) { // bb dominated -> only symbol for sstype	
						k_bb++; 
						nextnb=xlateSS( j_ss);
						// System.out.print(" -> BB ");
					} else { // sc dominated contact
						k_sc++; 
						nextnb=j_res;
						// System.out.print(" -> SC ");
					} // end if bb dominated
				} // end if we have non-zero entries
				// System.out.print(" : "+nextnb);
				nbs+= nextnb;
				nums+=" "+j_num; 
				fullnbs+= j_res; 
				k++;

				gapnbs+=gapsym+j_res;
				gaphood+=gapsym+nextnb; 
				lastnum=j_num; 
				// System.out.print( nextnb);
			}
			rsst.close(); 
			stmt.close();
			if (!overx) {
				nbs+="x"; // add if last in string == not seen yet
				fullnbs+="x"; 
				if (lastnum>0 && ((num-lastnum)==2)) xgapsym="_";
				if (lastnum>0 && ((num-lastnum)>2)) xgapsym="%";
				gapnbs+=xgapsym+"x"; 
				gaphood+=xgapsym+"x"; 
			}

		} catch (SQLException e) {
			e.printStackTrace();
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
		}
		return nbs; 
	} // end of getnbhoodstring

	public static String xlateSS( String sstype) { // xlate SecSTructure-Symbol 
		String SSymbol = "o";
		sstype=sstype.toUpperCase(); 
		if (sstype.equals("H")) SSymbol="z"; // Helix 
		if (sstype.equals("S")) SSymbol="b"; // beta strand 
		return SSymbol; 
	} // end xlateSS

} // end 

/* table definition for the results

CREATE TABLE `nbstrings` (
`graph_id` int(10) unsigned DEFAULT NULL,
`node_id` int(10) unsigned DEFAULT NULL,
`cid` varchar(6) COLLATE latin1_general_cs DEFAULT NULL,
`num` int(10) unsigned DEFAULT NULL,
`res` char(1) COLLATE latin1_general_cs DEFAULT NULL,
`sstype` char(1) COLLATE latin1_general_cs DEFAULT NULL,
`k` int(10) unsigned DEFAULT NULL,
`k_BB` int(10) unsigned DEFAULT NULL,
`k_SC` int(10) unsigned DEFAULT NULL,
`n` varchar(30) COLLATE latin1_general_cs DEFAULT NULL,
`nwg` varchar(40) COLLATE latin1_general_cs DEFAULT NULL,
`f` varchar(30) COLLATE latin1_general_cs DEFAULT NULL,
`fwg` varchar(40) COLLATE latin1_general_cs DEFAULT NULL
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_general_cs; 

*/ 