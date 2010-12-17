package owl.gmbp;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Vector;

import owl.core.util.MySQLConnection;

public class CMPdb_nbhString_traces {
	
    private MySQLConnection conn;
	
    // variables for queries
	private String host; // = "localhost" "talyn";
	private String username; // = "root" "vehlow";
	private String password; // = "nieve";
	private String db; // = "bagler_all13p0_alledges";
	
	private String jatom = "CA";
	private String nbhs = "%P%R%T%x%W%";
	private char sstype = 'S';
	private boolean diffSSType = true;
	private int maxNumLines = 100;
	
	// outputs (results)
	private Vector<float[]> nbhsNodes;
	private int[] numNodesPerLine;
	private int numLines;
	
	private static final char[] aas = new char[]{'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'}; // possible Residues
	private final String AAStr = new String(aas); 
	private static final char[] sstypes = new char[]{'H','S','O'};
	public final String SSTStr = new String(sstypes);
	public static final char AnySStype = 'A';
	
	public CMPdb_nbhString_traces(String nbhs, String jatom, String db){ // throws SQLException {
		this.db = db;
		this.jatom = jatom;
		this.nbhs = nbhs;
//		conn = new MySQLConnection();
//		conn = new MySQLConnection(this.host,this.username,this.password,this.db);
	}
	
	public CMPdb_nbhString_traces(String nbhs, String jatom){ // throws SQLException {
		this.jatom = jatom;
		this.nbhs = nbhs;
//		conn = new MySQLConnection();
//		conn = new MySQLConnection(this.host,this.username,this.password,this.db);
	}
	
	public void run() throws SQLException {
		String query;
		Statement stmt;
		ResultSet result_nbhs, result_nodes;
		
		int gID=0, num=0;
		int j_num=0;
		float theta=0, phi=0;
		char jRes='A', jSSType='H', iSSType='H';
		int resID=0, ssTypeID=0;
		
		this.nbhsNodes = new Vector<float[]>();
		
		stmt = conn.createStatement();	

//		System.out.println("extracting nbhstrings of certain type");
		query = "SELECT count(*) from "+db+".nbhstrings where nbhstring like '"+this.nbhs+"';";
		result_nbhs = stmt.executeQuery(query);
		int countAll =0;
		if (result_nbhs.next())
			countAll = result_nbhs.getInt(1); // extract raw count 
//		query = "SELECT graph_id, num from "+db+".nbhstrings where nbhstring like '"+this.nbhs+"' limit "+this.maxNumLines+";";
		query = "SELECT graph_id, num from "+db+".nbhstrings where nbhstring like '"+this.nbhs+"' ;";
		result_nbhs = stmt.executeQuery(query);
		Vector<int[]> nodes = new Vector<int[]>();
		while (result_nbhs.next()) {
			int[] val = {result_nbhs.getInt(1), result_nbhs.getInt(2)}; //new int[2];
			nodes.add(val);			
		}
		System.out.print("Orig num Traces = "+countAll+" -> number of nbhsTraces: "+this.maxNumLines+"\t");
		result_nbhs.close();
		stmt.close();
		this.numNodesPerLine = new int[nodes.size()];
		this.numLines = nodes.size();
		
		int[] numNodesSSType = new int[sstypes.length];
		int[][] numNodesSSTypeLine = new int[sstypes.length][nodes.size()];
		
//		System.out.println("graphi_id + '\t' + i_num + '\t' + j_num + '\t' + theta + '\t' + phi");
		int foundTraces = 0;
		for (int i=0; i<nodes.size(); i++){
			if (foundTraces<this.maxNumLines){

				stmt = conn.createStatement();
				int[] val = (int[]) nodes.get(i);
				gID = val[0];
				num = val[1];


				query = "SELECT j_num as j, theta, phi, j_res, j_sstype, i_sstype from "+db+".edges where graph_id="+gID+" and i_num="+num
                 +" and j_atom='"+this.jatom+"' order by j_num;";
//				if (this.diffSSType)
//					query = "SELECT j_num as j, theta, phi, j_res, j_sstype, i_sstype from "+db+".edges where graph_id="+gID+" and i_num="+num
//			         +" and i_sstype='"+this.sstype+"' and j_atom='"+this.jatom+"' order by j_num;";
				
				result_nodes = stmt.executeQuery(query);
				int cntNodes = 0;
				while(result_nodes.next()) {
					j_num = result_nodes.getInt(1); //("j_num");
					theta = result_nodes.getFloat(2);
					phi = result_nodes.getFloat(3);
					jRes = result_nodes.getString(4).charAt(0);
					jSSType = result_nodes.getString(5).charAt(0);
					iSSType = result_nodes.getString(6).charAt(0);
					resID = this.AAStr.indexOf(jRes);
					ssTypeID = this.SSTStr.indexOf(jSSType);
					
					int index = SSTStr.indexOf(iSSType);
					if (index>=0 && index<sstypes.length){
						numNodesSSType[index] += 1;
						numNodesSSTypeLine[index][i] += 1;						
					}
					
					if (iSSType == this.sstype || this.sstype == AnySStype || !this.diffSSType ){
						cntNodes++;
						this.numNodesPerLine[i] += 1;
						
//						if (i==0)
//						System.out.println(gID+" , "+num+" , "+j_num+" , "+theta+" , "+phi+" , "+jRes+"="+resID+" , "+jSSType+"="+ssTypeID);
						
						float[] node = {gID, num, j_num, theta, phi, resID, ssTypeID};
						this.nbhsNodes.add(node);					
					}
				}
//				if (cntNodes==0)
//					System.out.println("No result for query: "+query);
				if (cntNodes>0)
					foundTraces++;
							
				result_nodes.close();
				stmt.close();
				
			}
	    }
		System.out.print("foundTraces: "+foundTraces+"\n");
		
		// ------ToDo
		/* Count number of nodes for each sstype
		 * evaluate if correct sstype was handed over (significant more modes than for set sstype)
		 * change sstype in contactView
		 * */		
		System.out.print("Histogram Nodes for SSType"+"\t");
		int id=this.SSTStr.indexOf(String.valueOf(this.sstype));
		int mostOccSSType = 0;
		if (this.sstype!=AnySStype)
			mostOccSSType = numNodesSSType[id];
		else
			id = 0;
//		id = this.SSTStr.indexOf("S"); // this.sstypes.toString().indexOf("H");
//		// check
//		if (this.sstype != this.sstypes[id])
//			System.out.println("failed");
		for (int ssT=0; ssT<sstypes.length; ssT++){
			char ssType = sstypes[ssT];
			System.out.print(String.valueOf(ssType)+"\t");
//			for (int i=0; i<nodes.size(); i++){
//				System.out.print(String.valueOf(numNodesSSTypeLine[ssT][i])+"\t");
//		    }
			System.out.print("sum= "+String.valueOf(numNodesSSType[ssT])+"\t");
//			System.out.println();
			if (numNodesSSType[ssT]> 1.3*mostOccSSType){
				id = ssT;
				mostOccSSType = numNodesSSType[ssT];
			}
		}
		System.out.println();
		if (sstypes[id]!=this.sstype && this.sstype!=AnySStype){
			System.out.println("Mayority of traces for different ssType "+sstypes[id]+" instead of "+this.sstype);
//			this.sstype = sstypes[id]; 
//			this.run();
		}		
	}
	
	public void writeNbhsTracesOutput(String filename){
		CSVhandler csv = new CSVhandler();
		csv.generateCSVFile(this.nbhsNodes, filename);
	}
	
	public void setDb(String db) {
		this.db = db;
	}
	public String getDb() {
		return this.db;
	}
	
	public void setJAtom(String s){
		this.jatom = s;
	}
	public String getJAtom(){
		return this.jatom;
	}
	
	public void setNBHS(String s){
		this.nbhs = s;
	}
	public String getNBHS(){
		return this.nbhs;
	}
	
	public void setSSType(char c){
		this.sstype = c;
	}
	public char getSSType(){
		return this.sstype;
	}
	
	public void setDiffSSType(boolean b){
		this.diffSSType = b;
	}
	public boolean getDiffSSType(){
		return this.diffSSType;
	}
	public Vector<float[]> getNBHSnodes(){
		return this.nbhsNodes;
	}
	public int[] getNumNodesPerLine(){
		return this.numNodesPerLine;
	}
	public int getNumLines(){
		return this.numLines;
	}
	
	public int getMaxNumLines() {
		return maxNumLines;
	}
	public void setMaxNumLines(int maxNumLines) {
		this.maxNumLines = maxNumLines;
	}

	public void setDBaccess(String dbUSER, String dbPWD, String dbHOST, String dbNAME) throws SQLException {
		this.host = dbHOST;
		this.username = dbUSER;
		this.password = dbPWD;
		this.db = dbNAME;
		conn = new MySQLConnection(this.host,this.username,this.password,this.db);
	}

}
