import tools.MySQLConnection;

import java.sql.SQLException;
import java.sql.Statement;
import java.sql.ResultSet;

public class rewriteNBHood {

        /**
         * reriting the nbhoodstring such that BB dominated nbs are represented by SSType
         *
         * @author lappe
         */

        static String user = "lappe", fullnbs="", gapnbs="", nums="", gaphood=""; // change user name!!
        static MySQLConnection conn;
        static int k=0, k_sc=0, k_bb=0;

        public static void main(String[] args) throws SQLException {
                String nbhood = "", msql, cid, res, sectype;
                int graph_id=0, node_id=0, resnr=0; 
                // System.out.println("rewriteNBHood");
                conn = new MySQLConnection("white",user,"nieve","pdb_reps_graph_4_2");
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

                                nbhood= getSSNbHoodString( graph_id, node_id, resnr);
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

        public static String getSSNbHoodString( int gid, int nid, int num ) {
        	int j_num, j_bb, j_sc, lastnum=0;
        	boolean overx=false; 
        	String sql, nbs="", j_res, j_ss, nextnb, gapsym="", xgapsym=""; 
        	Statement stmt;  
        	ResultSet rsst;

        	try {
        		k=0; k_sc=0; k_bb=0; 
        		sql = "select j_num, j_res, j_sstype, BB, SC from chain_edge where graph_id="+gid+" and i_node_id="+nid+" order by j_num;";
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
        			nextnb="?"; 
        			if (lastnum>0 && ((j_num-lastnum)==2)) gapsym="_";
    				if (lastnum>0 && ((j_num-lastnum)>2)) gapsym="%";   
        			if (!overx && j_num>num) { 
        				overx=true; 
        				nbs+="x";
        				fullnbs+="x";
        				nums+=(" ("+num+") "); 
        				if (lastnum>0 && ((num-lastnum)==2)) xgapsym="_";
        				if (lastnum>0 && ((num-lastnum)>2)) xgapsym="%";   
        				gapnbs+=xgapsym+"x"; 
        				gaphood+=xgapsym+"x"; 
        				gapsym=""; 
        				if ((j_num-num)==2) gapsym="_";
        				if ((j_num-num)>2) gapsym="%";  
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

}
