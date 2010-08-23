package owl.gmbp;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.Vector;

import javax.vecmath.Vector3d;

import edu.uci.ics.jung.graph.util.Pair;

import owl.core.structure.AAinfo;
import owl.core.structure.Pdb;
import owl.core.structure.PdbLoadError;
import owl.core.structure.PdbfilePdb;
import owl.core.structure.Residue;
import owl.core.structure.graphs.RIGEdge;
import owl.core.structure.graphs.RIGGeometry;
import owl.core.structure.graphs.RIGNbhood;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.MySQLConnection;

public class CMPdb_sphoxel_greedy extends CMPdb_sphoxel{
	
	// --- neighbourhood defining variables
	private RIGNbhood nbhood;
	private String nbhString;

	private int[] nbSerials;
	private char[] nbhsRes;
	private Vector3d[] coordNBs;
	
	// --- contact defining parameters
	private Residue iResidue;
	private Residue jResidue;
	
	// --- position defining parameters
	private Vector3d contactPos;
	private static int factDelta = 2; 
	private double deltaR = defaultDeltaR; //0.5;
	private double deltaT = factDelta*defaultDeltaTheta;
	private double deltaP = factDelta*defaultDeltaPhi;
	
	private PrintWriter logOut;
	private Vector<String> tempTables;
	
	double expProb;

	/*
	 * Constructors
	 */
	public CMPdb_sphoxel_greedy() throws SQLException{
		this.iResidue = null;
		this.jResidue = null;
		this.nbhood = null;
		
		this.expProb = 0;
		this.tempTables = new Vector<String>();
		this.logOut =null;
		
		this.conn = new MySQLConnection(this.host,this.username,this.password,this.db);
	}
	
	public CMPdb_sphoxel_greedy(Residue iResidue, Residue jResidue, RIGNbhood nbhood) throws SQLException {
		this.iResidue = iResidue;
		this.jResidue = jResidue;
		initContactType();
		this.nbhood = nbhood;
		initNBHood();
		initPosParameters();
		
		this.logOut =null;

		this.tempTables = new Vector<String>();
		
		this.conn = new MySQLConnection(this.host,this.username,this.password,this.db);
	}
	
	public void initLogFile(String logFileFN) throws IOException{
		this.logOut = new PrintWriter(new FileWriter(logFileFN));		
	}
	
	/*
	 * private methods for initialisation
	 */
	private void initNBHood(){
		if (this.nbhood != null){
			this.nbhString = nbhood.getNbString();
			this.nbSerials = new int[nbhood.getSize()];
			this.nbhsRes = new char[this.nbhString.length()];
			int cnt=0;
			for (RIGNode node:nbhood.getNeighbors()){
				int resSer = node.getResidueSerial();
				this.nbSerials[cnt] = resSer;
				this.nbhsRes[cnt] = this.nbhString.charAt(cnt);
				cnt++;
				if (logOut!=null)
					logOut.print(resSer+"_"+node.getResidueType()+"\t");
				else
					System.out.print(resSer+"_"+node.getResidueType()+"\t");
//				System.out.print(resSer+"_"+this.mod.getNodeFromSerial(resSer).getResidueType()+"\t");
			}
			if (logOut!=null)
				logOut.println();
			else
				System.out.println();
		}
	}
	private void initContactType(){
		if (iResidue!=null && jResidue!=null){
			this.iRes = this.iResidue.getAaType().getOneLetterCode();
			this.jRes = this.jResidue.getAaType().getOneLetterCode();
			if (this.iResidue.getSsElem()!=null){
				this.diffSSType = true;
				this.issType = iResidue.getSsElem().getType();
			}
			else
				this.diffSSType = false;
			if (!diffSSType)
				if (logOut!=null)
					logOut.println("No Secondary Structure for Contact between "+this.iRes+this.iResidue.getSerial()
							+" and "+this.jRes+this.jResidue.getSerial());
				else
					System.out.println("No Secondary Structure for Contact between "+this.iRes+this.iResidue.getSerial()
						+" and "+this.jRes+this.jResidue.getSerial());
		}
	}
	private void initPosParameters(){
		this.contactPos = new Vector3d(new double[]{0,0,0});
		this.deltaR = defaultDeltaR; //0.5;
		this.deltaT = defaultDeltaTheta;
		this.deltaP = defaultDeltaTheta;
	}
	
	/*
	 * LOS computation --> counting methos 
	 */
		
	/**
	 * returns number of contacts between iRes and jRes
	 * @return count: number of contacts
	 * @throws SQLException
	 */
	private int countContactsOfType(String tName) throws SQLException{
		int count = 0;
		
		String query;
		Statement stmt;
		ResultSet result;
		
		stmt = conn.createStatement();				
		if (diffSSType)
			query = "SELECT count(*) from "+this.db+"."+tName+" where i_res='"+this.iRes+"' and i_sstype='"
				+this.issType+"' and j_res='"+jRes+"';";					
		else 
			query = "SELECT count(*) from "+this.db+"."+tName+" where i_res='"+this.iRes+"' and j_res='"+jRes+"';";	
//		System.out.println("countContactsOfType(): "+query);
		result = stmt.executeQuery(query);
		if(result.next()) 
			count = result.getInt(1); // extract raw count 
		result.close();
		stmt.close();
		
		return count;
	}
	
	/**
	 * returns number of contacts between iRes and jRes at specific position
	 * @param position
	 * @return count: number of contacts
	 * @throws SQLException
	 */
	private int countContactsOfTypeAt(Vector3d position, String tName) throws SQLException{
		int count = 0;
		double r = position.x;
		double theta = position.y;
		double phi = position.z;
		
		String query;
		Statement stmt;
		ResultSet result;
		
		stmt = conn.createStatement();	
		if (diffSSType)
			query = "SELECT count(*) from "+this.db+"."+tName+" where i_res='"+this.iRes+"' and i_sstype='"
				+this.issType+"' and j_res='"+jRes+"' and r >= "+(r-this.deltaR)+" and r < "+(r+this.deltaR)
				+" and theta>"+(theta-this.deltaT)+" and theta<"+(theta+this.deltaT)
				+" and phi>"+(phi-this.deltaP)+" and phi<"+(phi+this.deltaP)+";";					
		else 
			query = "SELECT count(*) from "+this.db+"."+tName+" where i_res='"+this.iRes+"' and j_res='"+jRes+
				"' and r >= "+(r-this.deltaR)+" and r < "+(r+this.deltaR)
				+" and theta>"+(theta-this.deltaT)+" and theta<"+(theta+this.deltaT)
				+" and phi>"+(phi-this.deltaP)+" and phi<"+(phi+this.deltaP)+";";	
//		System.out.println("countContactsOfTypeAt(): "+query);
		result = stmt.executeQuery(query);
		if(result.next()) 
			count = result.getInt(1); // extract raw count 
		result.close();
		stmt.close();
		
		return count;
	}
	
	/**
	 * Create new table of all edges starting in iRes for those graphs which contain an edge starting in iRes and ending in kRes
	 * @param kResidues (kRes that is element of neighbourhood of iRes)
	 * @param tName (Name for temporal table)
	 * @return
	 * @throws SQLException 
	 */
	private boolean createEdgeTable(char kRes, String baseTable, String tName) throws SQLException{
		boolean successful = false;
		
		Vector<Integer> graphIDs = new Vector<Integer>();
		
		String query="";
		Statement stmt;
		ResultSet result;
		
		stmt = conn.createStatement();
		
		String kCondition = " and j_res='"+kRes+"'";
		
		if(diffSSType)
			query = "SELECT graph_id from "+db+"."+baseTable+" where i_res='"+this.iRes+"' and i_sstype='"+this.issType+"'"+kCondition+";";
		else
			query = "SELECT graph_id from "+db+"."+baseTable+" where i_res='"+this.iRes+"' "+kCondition+";";
//		System.out.println("kConditions: "+query);
		result = stmt.executeQuery(query);
		while (result.next()){
			int id = result.getInt(1);
			if (!graphIDs.contains(id))
				graphIDs.add(id);
		}
		String idCond = "";
		int cnt = 0;
		if (graphIDs.size()>0){
			for (int id:graphIDs){
				if (cnt==0)
					idCond += ("and (graph_id="+id);
				else
					idCond += (" or graph_id="+id);
//				System.out.print(id+", ");
				cnt++;
			}
			idCond += ")";
//			System.out.println();
//			System.out.println(idCond);			
		}
		
		query = "DROP TABLE IF EXISTS "+db+"."+tName+";";
		stmt.execute(query);
		
		if (diffSSType)
			query = "CREATE table "+this.db+"."+tName+" (SELECT * from "+this.db+"."+baseTable+" where i_res='"+this.iRes+"' and i_sstype='"
				+this.issType+"' "+idCond+");";
		else
			query = "CREATE table "+this.db+"."+tName+" (SELECT * from "+this.db+"."+baseTable+" where i_res='"+this.iRes+"' "+idCond+");";
//		System.out.println("idCond: "+query);
		stmt.execute(query);
		this.tempTables.add(tName);
				
		result.close();
		stmt.close();
		
		return successful;
	}
	
	/**
	 * Create new table of all edges starting in iRes for those graphs which contain an edge starting in iRes and ending in kRes
	 * @param kResidues (char-Array that contains all residues of neighbourhood of iRes)
	 * @param tName (Name for temporal table)
	 * @return
	 * @throws SQLException 
	 */
	@SuppressWarnings("unused")
	private boolean createEdgeTable(char[] kResidues, String tName) throws SQLException{
		boolean successful = false;
		
		Vector<Integer> graphIDs = new Vector<Integer>();
		
		String query="";
		Statement stmt;
		ResultSet result;
		
		stmt = conn.createStatement();
		
		String kCondition = "";
		if (kResidues.length>0){
			kCondition += " and (j_res='"+kResidues[0]+"'";
			if (kResidues.length>1){
				for (int i=1; i<kResidues.length; i++){
					kCondition += " or j_res='"+kResidues[i]+"'";
				}				
			}
			kCondition += ")";
		}		
		if(diffSSType)
			query = "SELECT graph_id from "+db+".edges where i_res='"+this.iRes+"' and i_sstype='"+this.issType+"'"+kCondition+";";
		else
			query = "SELECT graph_id from "+db+".edges where i_res='"+this.iRes+"' "+kCondition+";";
		if (logOut!=null)
			logOut.println("kConditions: "+query);
		else
			System.out.println("kConditions: "+query);
		result = stmt.executeQuery(query);
		while (result.next()){
			int id = result.getInt(1);
			if (!graphIDs.contains(id))
				graphIDs.add(id);
		}
		String idCond = "";
		int cnt = 0;
		if (graphIDs.size()>0){
			for (int id:graphIDs){
				if (cnt==0)
					idCond += ("and (graph_id="+id);
				else
					idCond += (" or graph_id="+id);
//				System.out.print(id+", ");
				cnt++;
			}
			idCond += ")";
//			System.out.println();
//			System.out.println(idCond);			
		}
		
		query = "DROP TABLE IF EXISTS "+db+"."+tName+";";
		stmt.execute(query);
		
		if (diffSSType)
			query = "CREATE table "+this.db+"."+tName+" (SELECT * from "+this.db+".edges where i_res='"+this.iRes+"' and i_sstype='"
				+this.issType+"' "+idCond+");";
		else
			query = "CREATE table "+this.db+"."+tName+" (SELECT * from "+this.db+".edges where i_res='"+this.iRes+"' "+idCond+");";
		if (logOut!=null)
			logOut.println("idCond: "+query);
		else
			System.out.println("idCond: "+query);
		stmt.execute(query);	
		this.tempTables.add(tName);
				
		result.close();
		stmt.close();
		
		return successful;
	}
	
	/**
	 * @throws SQLException 
	 * 
	 */
	public void greedy4ObsScore() throws SQLException{
		
		String baseTName = CMPdb_sphoxel.defaultEdgeTable; // "edges";
		Vector<Integer> kResIDs = new Vector<Integer>();
		kResIDs = greedy4ObsScore(baseTName, kResIDs);
//		if (kResIDs.size()!=this.nbhsRes.length)
		if (logOut!=null)
			logOut.print(this.iRes+"_"+this.jRes);
		else
			System.out.print(this.iRes+"_"+this.jRes);
		for (int kID:kResIDs){
			if (logOut!=null)
				logOut.print("\t-->"+this.nbhsRes[kID]);
			else
				System.out.print("\t-->"+this.nbhsRes[kID]);
		}
		if (logOut!=null)
			logOut.print("\n");
		else
			System.out.print("\n");
		
//		for (int k=0; k<this.nbhsRes.length; k++){
//			String tName = "edges_kRes_"+nbhsRes[k];
//			createEdgeTable(nbhsRes[k], tName);
//			int cntTijk = countContactsOfType(tName) +1;
//			int cntTijkPos = countContactsOfTypeAt(contactPos, tName) +1;
//			double obs = (double)cntTijkPos/(double)cntTijk;
//			System.out.println("Obs for k="+nbhsRes[k]+": "+cntTijkPos+"/"+cntTijk+" = "+obs);
//			if (obs>bestObs){
//				bestObs = obs;
//				bestK = k;
//				bestTN = tName;
//			}
//		}
//		System.out.println("Best k="+bestK+" with obs="+bestObs+" of table="+bestTN);
	}
	
	/**
	 * @throws SQLException 
	 * 
	 */
	public Vector<Integer> greedy4ObsScore(String baseTName, Vector<Integer> kResIDs) throws SQLException{
		int bestK = 0;
		double bestObs = 0;
		String bestTN = "";
		for (int k=0; k<this.nbhsRes.length; k++){
			if (!kResIDs.contains(k)){
				String tName = baseTName+"_"+nbhsRes[k];
				createEdgeTable(nbhsRes[k], baseTName, tName);
				int cntTijk = countContactsOfType(tName) +1;
				int cntTijkPos = countContactsOfTypeAt(contactPos, tName) +1;
				double obs = (double)cntTijkPos/(double)cntTijk;
				if (logOut!=null)
					logOut.println("tName="+tName+"Obs for k="+nbhsRes[k]+": "+cntTijkPos+"/"+cntTijk+" = "
							+obs+" obs/exp="+(obs/this.expProb)+" log="+Math.log(obs/this.expProb));
				else
					System.out.println("tName="+tName+"Obs for k="+nbhsRes[k]+": "+cntTijkPos+"/"+cntTijk+" = "
						+obs+" obs/exp="+(obs/this.expProb)+" log="+Math.log(obs/this.expProb));
				if (obs>bestObs){
					bestObs = obs;
					bestK = k;
					bestTN = tName;
				}				
			}
		}
		kResIDs.add(bestK);
		if (logOut!=null)
			logOut.println("Best k="+bestK+" with obs="+bestObs+" of table="+bestTN);
		else
			System.out.println("Best k="+bestK+" with obs="+bestObs+" of table="+bestTN);
		if (kResIDs.size()!=this.nbhsRes.length)
			kResIDs = greedy4ObsScore(bestTN, kResIDs);
		
//		for (int k=0; k<this.nbhsRes.length; k++){
//			String tName = "edges_kRes_"+nbhsRes[k];
//					
//			String query = "DROP TABLE IF EXISTS "+db+"."+tName+";";
//			Statement stmt;			
//			stmt = conn.createStatement();
//			stmt.execute(query);
//			stmt.close();
//		}
		
		return kResIDs;
	}
	
	/**
	 * @throws SQLException 
	 * 
	 */
	public void createEnvTable() throws SQLException{
		String baseTName = CMPdb_sphoxel.defaultEdgeTable;
		String templateT = "edges_resenv";
		String tName = templateT+"_"+this.iRes+"_"+String.valueOf(nbhsRes);
		
		if (diffSSType)
			tName += "_diffSST";

		logOut.println("createEnvTable for iRes="+this.iRes+" jRes="+this.jRes+"  within env:"+String.valueOf(this.nbhsRes));

//		double r = this.contactPos.x;
//		double theta = this.contactPos.y;
//		double phi = this.contactPos.z;
		
		String query="";
		Statement stmt;
//		ResultSet result;
				
		stmt = conn.createStatement();

		query = "DROP TABLE IF EXISTS "+db+"."+tName+";";
		stmt.execute(query);
		query = "CREATE table "+this.db+"."+tName+" (SELECT * from "+this.db+"."+templateT+" limit 0);";
		stmt.execute(query);
		logOut.println(query);
//		
		char kRes = this.jRes;
//		query = "INSERT INTO "+tName+" SELECT * from "+this.db+"."+baseTName+" where i_res='"+this.iRes+"' and i_sstype='"
//			+this.issType+"' and j_res='"+kRes+"' and r >= "+(r-this.deltaR)+" and r < "+(r+this.deltaR)
//			+" and theta>"+(theta-this.deltaT)+" and theta<"+(theta+this.deltaT)
//			+" and phi>"+(phi-this.deltaP)+" and phi<"+(phi+this.deltaP)+";";
//		logOut.println(query);
//		stmt.execute(query);	
		if (this.nbhsRes.length!=this.coordNBs.length)
			System.out.println("Warning: input arays of different length!");
		for (int k=0; k<this.nbhsRes.length; k++){
			Vector3d pos = this.coordNBs[k];
			double r = pos.x;
			double theta = pos.y;
			double phi = pos.z;
			kRes = nbhsRes[k];
			int jNum = this.nbSerials[k];
			if (diffSSType)
				query = "INSERT INTO "+tName+" SELECT * from "+this.db+"."+baseTName+" where i_res='"+this.iRes+"' and i_sstype='"
				    +this.issType+"' and j_res='"+kRes+"' and r >= "+(r-this.deltaR)+" and r < "+(r+this.deltaR)
				    +" and theta>"+(theta-this.deltaT)+" and theta<"+(theta+this.deltaT)
				    +" and phi>"+(phi-this.deltaP)+" and phi<"+(phi+this.deltaP)+";";
			else
				query = "INSERT INTO "+tName+" SELECT * from "+this.db+"."+baseTName+" where i_res='"+this.iRes
					+"' and j_res='"+kRes+"' and r >= "+(r-this.deltaR)+" and r < "+(r+this.deltaR)
					+" and theta>"+(theta-this.deltaT)+" and theta<"+(theta+this.deltaT)
					+" and phi>"+(phi-this.deltaP)+" and phi<"+(phi+this.deltaP)+";";
			logOut.println("ResNr:"+jNum+"  "+query);
			stmt.execute(query);	
		}
				
//		result.close();
		stmt.close();	
	}
		
	/**
	 * 
	 * @param contactCoord
	 * @return
	 * @throws SQLException
	 */
	public String createEnvTable(HashMap<Pair<Integer>,Vector3d> contactCoord, String pdbName) throws SQLException{
		String baseTName = CMPdb_sphoxel.defaultEdgeTable;
		String templateT = "edges_resenv";
		char kRes = this.jRes;
		int iNum = this.iResidue.getSerial();
		int jNum = this.jResidue.getSerial();
		
		String tName = templateT+"_"+pdbName+"_"+iNum+this.iRes+"_"+jNum+"_"+String.valueOf(nbhsRes);
		
		if (diffSSType)
			tName += "_diffSST";
		
		logOut.println("E"+iNum+this.iRes+"-SSType="+issType+"  nbhString="+this.nbhString); //+"_"+this.jResidue.getSerial()+this.jRes);
		for (int i=0; i<this.nbSerials.length; i++)
			logOut.print(this.nbSerials[i]+"_"+this.nbhsRes[i]+"\t");
		logOut.print("\n");
		logOut.println("createEnvTable "+tName+" for iRes="+this.iRes+" jRes="+this.jRes+"  within env:"+String.valueOf(this.nbhsRes));

		String query="";
		Statement stmt;
				
		stmt = conn.createStatement();

		query = "DROP TABLE IF EXISTS "+db+"."+tName+";";
		stmt.execute(query);
		query = "CREATE table "+this.db+"."+tName+" (SELECT * from "+this.db+"."+templateT+" limit 0);";
		stmt.execute(query);
		logOut.println(query);
//		
//		if (this.nbhsRes.length!=this.coordNBs.length)
//			System.out.println("Warning: input arays of different length!");
		for (int k=0; k<this.nbhsRes.length; k++){
			kRes = nbhsRes[k];
			jNum = this.nbSerials[k];
			Pair<Integer> key = new Pair<Integer>(iNum, jNum);
			Vector3d pos = contactCoord.get(key);
			double r = pos.x;
			double theta = pos.y;
			double phi = pos.z;
			
			if (diffSSType)
				query = "INSERT INTO "+tName+" SELECT * from "+this.db+"."+baseTName+" where i_res='"+this.iRes+"' and i_sstype='"
				    +this.issType+"' and j_res='"+kRes+"' and r >= "+(r-this.deltaR)+" and r < "+(r+this.deltaR)
				    +" and theta>"+(theta-this.deltaT)+" and theta<"+(theta+this.deltaT)
				    +" and phi>"+(phi-this.deltaP)+" and phi<"+(phi+this.deltaP)+";";
			else
				query = "INSERT INTO "+tName+" SELECT * from "+this.db+"."+baseTName+" where i_res='"+this.iRes
					+"' and j_res='"+kRes+"' and r >= "+(r-this.deltaR)+" and r < "+(r+this.deltaR)
					+" and theta>"+(theta-this.deltaT)+" and theta<"+(theta+this.deltaT)
					+" and phi>"+(phi-this.deltaP)+" and phi<"+(phi+this.deltaP)+";";
			logOut.println("ResNr:"+jNum+"  "+query);
			stmt.execute(query);	
		}
		
		stmt.close();	
		
		return tName;
	}
	

	public void evaluateEnvTable(String tName) throws SQLException{
//		Select accession_code, cid, i_num,count(*) as c from edges_resenv_t_vpftgriia group by accession_code, cid, i_num order by c desc;
		String query="";
		Statement stmt;
		ResultSet result;
		
		String cntTName = "nbh_count";
				
		stmt = conn.createStatement();
		
		query = "DROP TABLE IF EXISTS "+db+"."+cntTName+";";
		stmt.execute(query);
		query = "CREATE table "+this.db+"."+cntTName+" (accession_code char(4), cid varchar(6), i_num int(11), c int(2));";
		stmt.execute(query);
		
		query = "INSERT INTO "+cntTName+" Select accession_code, cid, i_num,count(*) as c from "+this.db+"."+tName
			+" group by accession_code, cid, i_num order by c desc;";
		stmt.execute(query);
		
		query = "Select c,count(*) as cntC from "+this.db+"."+cntTName+" group by c order by cntC desc;";
		result = stmt.executeQuery(query);
		while (result.next()){
			// count all occurances of all counts: e.g. 1x cnt=9, 2x cnt=2, 70x cnt=1
			int c = result.getInt(1);
			int cntC = result.getInt(2);
			
			logOut.print(c+"->"+cntC+"\t");
		}
		logOut.print("\n");
		
		query = "DROP TABLE IF EXISTS "+db+"."+cntTName+";";
		stmt.execute(query);
		
		stmt.close();			
	}
	
	
	/**
	 * @throws SQLException 
	 * 
	 */
	public void computeExpScore() throws SQLException{
		if (logOut!=null){
			logOut.println("Greedy for Contact of type: "+this.iRes+"->"+this.jRes+" "+this.issType
					+" at pos: "+contactPos.toString()+" kRes="+String.valueOf(nbhsRes));
			logOut.println("Range of r = +-"+this.deltaR+" and theta= +-"+this.deltaT+" and phi= +-"+this.deltaP);
		}
		else {
			System.out.println("Greedy for Contact of type: "+this.iRes+"->"+this.jRes+" "+this.issType
					+" at pos: "+contactPos.toString()+" kRes="+String.valueOf(nbhsRes));
			System.out.println("Range of r = +-"+this.deltaR+" and theta= +-"+this.deltaT+" and phi= +-"+this.deltaP);			
		}
//		System.out.println("SELECT count(*) from "+this.db+".edges where i_res='"+this.iRes+"' and i_sstype='"
//				+this.issType+"' and j_res='"+jRes+"';");
//		System.out.println("SELECT count(*) from "+this.db+".edges where i_res='"+this.iRes+"' and i_sstype='"
//				+this.issType+"' and j_res='"+jRes+"' and r >= (r"+(-this.deltaR)+" and r < r"+(+this.deltaR)
//				+" and theta> theta"+(-this.deltaT)+" and theta< theta"+(+this.deltaT)
//				+" and phi> phi"+(-this.deltaP)+" and phi< phi"+(+this.deltaP)+";");
		// compute expected value
		int cntTij = countContactsOfType(defaultEdgeTable) +1;
		int cntTijPos = countContactsOfTypeAt(this.contactPos, defaultEdgeTable) +1;
		double exp = (double)cntTijPos/(double)cntTij;
		
		this.expProb = exp;
		
		if (logOut!=null)
			logOut.println("Exp: "+cntTijPos+"/"+cntTij+" = "+exp);
		else
			System.out.println("Exp: "+cntTijPos+"/"+cntTij+" = "+exp);
	}
	
	public void dropTemporalTables() throws SQLException{
		if (this.tempTables.size()>0){
			Statement stmt;	
			stmt = conn.createStatement();
			for (String tName:this.tempTables){
				String query = "DROP TABLE IF EXISTS "+db+"."+tName+";";
				stmt.execute(query);
			}			
			stmt.close();
			
			this.tempTables = new Vector<String>();
		}
	}
	
	public void closeLogFile() {
		if (logOut!=null)
			this.logOut.close();
	}
	
	/* ----- GETTERS and SETTERS ------------	 */	
	public char[] getNbhsRes() {
		return nbhsRes;
	}
	public void setNbhsRes(char[] nbhsRes) {
		this.nbhsRes = nbhsRes;
	}
	
	public void setCoordNBs(Vector3d[] coordNBs){
		this.coordNBs = coordNBs;
	}
	
	public void setNbhsSer(int[] serials){
		this.nbSerials = serials;
	}
	
	public Residue getiResidue() {
		return iResidue;
	}

	public void setiResidue(Residue iResidue) {
		this.iResidue = iResidue;
	}

	public Residue getjResidue() {
		return jResidue;
	}

	public void setjResidue(Residue jResidue) {
		this.jResidue = jResidue;
	}

	public String getNbhString() {
		return nbhString;
	}
	public void setNbhString(String nbhString) {
		this.nbhString = nbhString;
	}





	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

//		RIGNbhood nbhood = this.mod.getGraph().getNbhood(nodeI);
		
//		double r=resDist;
//		theta = graphGeom.getRotCoordOfContacts().get(new Pair<Integer>(iNum, j)).y;
//		phi = graphGeom.getRotCoordOfContacts().get(new Pair<Integer>(iNum, j)).z;
		
		CMPdb_sphoxel_greedy greedy = null;
		String proteinName = "1bxy";
		String proteinFN = "/Users/vehlow/Documents/workspace/PDBs/1bxy.pdb";
		double cutOff=8.0;
		String contactType="SC";	
		RIGraph graph = null;
		RIGGeometry graphGeom = null;
		
		Pdb pdb = new PdbfilePdb(proteinFN);
		int modNr = 1;
		try {
			pdb.load("A",modNr);
		} catch (PdbLoadError e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if (pdb!=null && pdb.getRIGraph(contactType, cutOff)!=null) {
			graph = pdb.getRIGraph(contactType, cutOff);
			graphGeom = new RIGGeometry(graph, pdb.getResidues());
			
			try {
				greedy = new CMPdb_sphoxel_greedy();
				String logFileFN = "/Users/vehlow/Documents/workspace/outputFiles/evaluateNBH_1bxy.txt";
				greedy.initLogFile(logFileFN);
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			TreeMap<Integer,Residue> residues = pdb.getResidues();
			
			System.out.println(proteinName+" has "+graph.getEdgeCount()+" edges");
			
			for (RIGEdge edge:graph.getEdges()) {
				// get Neighbourhood
				Pair<RIGNode> nodes = graph.getEndpoints(edge);
				RIGNode iNode = nodes.getFirst();
				RIGNode jNode = nodes.getSecond();
				int iNum = iNode.getResidueSerial();
				int jNum = jNode.getResidueSerial();
				String iResType = iNode.getResidueType();
				String jResType = jNode.getResidueType();
				char issType = 0;
				boolean diffSST = false;
				if (iNode.getSecStrucElement() != null){
					diffSST = true;
					issType = iNode.getSecStrucElement().getType();
				}
				
				Residue iRes = residues.get(iNum);
				Residue jRes = residues.get(jNum);
				
				System.out.println("E"+iNum+iResType+"-"+issType+"_"+jNum+jResType);
				
				RIGNbhood nbhood = graph.getNbhood(iNode);
				String nbhString = nbhood.getNbString();
				int[] nbSerials = new int[nbhood.getSize()];
				char[] nbhsRes = new char[nbhood.getSize()];
				int cnt=0;
				int i = 0;
				for (RIGNode node:nbhood.getNeighbors()){
					if (nbhString.charAt(cnt)=='x'){
						cnt++;
					}
					{
						int resSer = node.getResidueSerial();
						nbSerials[i] = resSer;
						nbhsRes[i] = nbhString.charAt(cnt);
						i++;
						System.out.print(resSer+"_"+node.getResidueType()+"\t");
//						System.out.print(resSer+"_"+this.mod.getNodeFromSerial(resSer).getResidueType()+"\t");						
					}
					cnt++;
				}
				System.out.println();
				// create Table
				greedy.setIRes(AAinfo.threeletter2oneletter(iResType).charAt(0)); //('T');
				greedy.setJRes(AAinfo.threeletter2oneletter(jResType).charAt(0)); //('A');			
				greedy.setNbhsRes(nbhsRes); //.nbhsRes = new char[]{'V','P','F','T','G','R','I','I','A'};
				greedy.setNbhsSer(nbSerials); // .nbSerials = new int[]{62,64,65,85,86,270,271,272,273};
				greedy.setDiffSSType(diffSST);
				greedy.setISSType(issType);
				greedy.setiResidue(iRes);
				greedy.setNbhString(nbhString);
				greedy.setjResidue(jRes);
				
				String tName = null;
				try {
					tName = greedy.createEnvTable(graphGeom.getRotCoordOfContacts(),proteinName);
				} catch (SQLException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				// count accession code appearances
				if (tName!=null){
					try {
						greedy.evaluateEnvTable(tName);
					} catch (SQLException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
//				greedy.closeLogFile();
			}
			greedy.closeLogFile();
		}
		
		//Select * from edges_resenv_t_vpftgri order by graph_id, accession_code;
		//Select accession_code, cid, i_num,count(*) as c from edges_resenv_t_vpftgriia group by accession_code, cid, i_num;
	    // Select accession_code, cid, i_num,count(*) as c from edges_resenv_t_vpftgriia group by accession_code, cid, i_num order by c desc;
		// Select * from edges_resenv_a_dritfatltsf order by graph_id, accession_code;
		// Select accession_code, cid, i_num,count(*) as c from edges_resenv_t_vpftgriia group by accession_code, cid, i_num order by c desc;
		// SELECT * from bagler_all13p0_alledges.edges where graph_id=1 and i_num=354 order by j_num;
		
		
//		try {
//			/*
//			 * Selected contact changed to: 63T_S  273A
//			 * Filename= sphoxelBG_T-s_A-a_rSR.csv
//			 * Edge type: Cb
//			 * 62_VAL	64_PRO	65_PHE	85_THR	86_GLY	270_ARG	271_ILE	272_ILE	273_ALA
//			 * 62_VAL: (3.790542309485543, 0.41603364996679265, -2.411607228178422)
//			 * 64_PRO: (3.8022236914731895, 1.7471926158531828, -0.27430555619890346)
//			 * 65_PHE: (6.546326068872526, 1.911898086059309, 0.2548014843008439)
//			 * 85_THR: (7.229295263025298, 2.6384044929128865, -0.38285092655450975)
//			 * 86_GLY: (7.422997642462245, 2.441628539832542, 0.4942550998255656)
//			 * 270_ARG: (8.360183311387377, 1.0832183461851157, 2.169960157364696)
//			 * 271_ILE: (6.474047651971675, 0.6548182513501589, 1.9608615320692695)
//			 * 272_ILE: (4.517631016362445, 0.8384230696283604, 1.0719900380035916)
//			 * 273_ALA: (6.255787959961556, 0.9652239495130668, 0.24672550715132224)
//			 * NBHString: VxPFTGRIIA-->%V%x%P%F%T%G%R%I%I%A%
//			 * r=6.26 theta=0.97 phi=0.25-->3.39
//			 */
//			/*
//			 * coordinates within edges table: r>0, theta[0:Pi], phi[-Pi:=Pi]
//			 */
//			greedy = new CMPdb_sphoxel_greedy();
//			String logFileFN = "/Users/vehlow/Documents/workspace/outputFiles/greedyTest2.txt";
//			greedy.initLogFile(logFileFN);
//			greedy.setIRes('T');
//			greedy.setJRes('A');			
//			greedy.nbhsRes = new char[]{'V','P','F','T','G','R','I','I','A'};
//			greedy.nbSerials = new int[]{62,64,65,85,86,270,271,272,273};
//			greedy.setDiffSSType(false);
//			greedy.setISSType('S');
//			greedy.contactPos = new Vector3d(new double[]{6.26, 0.97, 0.25});
//			
//			Vector3d[] coordNBs = new Vector3d[greedy.nbhsRes.length];
//			coordNBs[0] = new Vector3d(new double[]{3.790542309485543, 0.41603364996679265, -2.411607228178422});
//			coordNBs[1] = new Vector3d(new double[]{3.8022236914731895, 1.7471926158531828, -0.27430555619890346});
//			coordNBs[2] = new Vector3d(new double[]{6.546326068872526, 1.911898086059309, 0.2548014843008439});
//			coordNBs[3] = new Vector3d(new double[]{7.229295263025298, 2.6384044929128865, -0.38285092655450975});
//			coordNBs[4] = new Vector3d(new double[]{7.422997642462245, 2.441628539832542, 0.4942550998255656});
//			coordNBs[5] = new Vector3d(new double[]{8.360183311387377, 1.0832183461851157, 2.169960157364696});
//			coordNBs[6] = new Vector3d(new double[]{6.474047651971675, 0.6548182513501589, 1.9608615320692695});
//			coordNBs[7] = new Vector3d(new double[]{4.517631016362445, 0.8384230696283604, 1.0719900380035916});
//			coordNBs[8] = new Vector3d(new double[]{6.255787959961556, 0.9652239495130668, 0.24672550715132224});
//			greedy.setCoordNBs(coordNBs);
//			
//			TreeMap<Integer, Vector3d> nbhResCoord = new TreeMap<Integer, Vector3d>();
//			nbhResCoord.put(0, new Vector3d(new double[]{3.790542309485543, 0.41603364996679265, -2.411607228178422}));
//			nbhResCoord.put(1, new Vector3d(new double[]{3.8022236914731895, 1.7471926158531828, -0.27430555619890346}));
//			nbhResCoord.put(2, new Vector3d(new double[]{6.546326068872526, 1.911898086059309, 0.2548014843008439}));
//			nbhResCoord.put(3, new Vector3d(new double[]{7.229295263025298, 2.6384044929128865, -0.38285092655450975}));
//			nbhResCoord.put(4, new Vector3d(new double[]{7.422997642462245, 2.441628539832542, 0.4942550998255656}));
//			nbhResCoord.put(5, new Vector3d(new double[]{8.360183311387377, 1.0832183461851157, 2.169960157364696}));
//			nbhResCoord.put(6, new Vector3d(new double[]{6.474047651971675, 0.6548182513501589, 1.9608615320692695}));
//			nbhResCoord.put(7, new Vector3d(new double[]{4.517631016362445, 0.8384230696283604, 1.0719900380035916}));
//			nbhResCoord.put(8, new Vector3d(new double[]{6.255787959961556, 0.9652239495130668, 0.24672550715132224}));
//			
//			
////			System.out.println("iRes="+greedy.iRes+" iSS="+greedy.issType+"   jRes="+greedy.jRes+"    contactPos={3.5, 1.1, 0.5}    greedy.nbhsRes={'P','L','F'}");
//			
////			greedy.computeExpScore();
////			greedy.greedy4ObsScore();
////			greedy.dropTemporalTables();
//			
//			greedy.createEnvTable();
//			
//			//Select * from edges_resenv_t_vpftgri order by graph_id, accession_code;
//			//Select accession_code, cid, i_num,count(*) as c from edges_resenv_t_vpftgriia group by accession_code, cid, i_num;
		    // Select accession_code, cid, i_num,count(*) as c from edges_resenv_t_vpftgriia group by accession_code, cid, i_num order by c desc;
//			
//			greedy.closeLogFile();
//			
//		} catch (SQLException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		
		
	}

}
