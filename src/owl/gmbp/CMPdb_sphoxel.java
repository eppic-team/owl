package owl.gmbp;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import edu.uci.ics.jung.graph.util.Pair;

import owl.core.structure.Residue;
import owl.core.structure.graphs.RIGGeometry;
import owl.core.util.MySQLConnection;

public class CMPdb_sphoxel {
	
	public static final int defaultNumSteps = 4;
	public static final float defaultResol = 180/(float)defaultNumSteps; // in degrees
	protected static final String defaultDB = "bg";
	public static final String[] radiusRanges = {"rSR","rMR","rLR"};
	public static final char AnySStype = 'A';
	
	/*--------------------------- member variables --------------------------*/		
	private MySQLConnection conn;
	
	private String host = "talyn";
	private String username = "vehlow";
	private String password = "nieve";
	private String db = "bagler_all13p0_alledges";
	private String tableNameRes = "edgesCertainRes";
	private String tableNameResR = "edgesCertainResRRange";
	private String tableNameResRT = "edgesCertainResRTRange";
	private String tableNameR = "edgesRrange";
	private String tableNameRT = "edgesRTrange";
	private final String tn = "edges";
	private String tnRes = "";
	
    private boolean diffSSType = false;
	
	private int numSteps = defaultNumSteps;
	private int numRatiosX = numSteps;
	private int numRatiosY = 2*numSteps;
	private double [][] ratios;
	private double [][][] bayesRatios;
	private double minRatio = 0;
	private double maxRatio = 0;
	
	private float minr=0.0f, maxr=9.2f;
	private String radiusPrefix = "rSR"; // "rMR" or "rLR"
	private final double mintheta=0.0, minphi=-Math.PI;
	private double svoxelsize=(Math.PI)/numSteps; //, deltar=1.0 ;
	private char iRes='A', jRes='A';
	private char issType='H', jssType='H';
//	private final char[] sstype = {'H', 'S', 'O', 'A'};
	
	public CMPdb_sphoxel(char iRes, char jRes, String db) throws SQLException {
		this.db = db;
		this.iRes = iRes;
		this.jRes = jRes;
//		conn = new MySQLConnection();
		conn = new MySQLConnection(this.host,this.username,this.password,this.db);
	}
	
	public CMPdb_sphoxel() throws SQLException {
		conn = new MySQLConnection(this.host,this.username,this.password,this.db);
	}
	
	public double getLogOddsScore(){
		double score = 0;
		
		return score;
	}
	
	public double getLogOddsScore(Residue iRes, Residue jRes, double resDist, RIGGeometry graphGeom, double dAngle, double dRad) throws SQLException{
		String query;
		Statement stmt;
		ResultSet result_angle, result_type, result_angle_type, result_all;
		double theta=0.0, phi=0.0;
		double countAll = 0, countAngle = 0, countType = 0, countAngleType = 0;
		double countExp = 0;
		double countObs = 0;
		double lOS=0;
		
		// extract parameters
		if (iRes.getSsElem()!=null)
			diffSSType=true;
		else
			diffSSType=false;
		
		int iNum=iRes.getSerial();
		int j=jRes.getSerial();
		
		double r=resDist;
		theta = graphGeom.getRotCoordOfContacts().get(new Pair<Integer>(iNum, j)).y;
		phi = graphGeom.getRotCoordOfContacts().get(new Pair<Integer>(iNum, j)).z;
		
		this.db = defaultDB;
//		if (diffSSType)
//			tnRes = "edges_"+this.iRes+"_"+String.valueOf(this.issType).toLowerCase()+"_"+this.jRes+"_"+String.valueOf(AnySStype).toLowerCase();
//		else
//			tnRes = "edges_"+this.iRes+"_"+String.valueOf(AnySStype).toLowerCase()+"_"+this.jRes+"_"+String.valueOf(AnySStype).toLowerCase();	
		
		stmt = conn.createStatement();		
		// ---- count all
		query = "SELECT count(*) from "+this.db+".edges;";
		result_all = stmt.executeQuery(query);
		if(result_all.next()) 
			countAll = result_all.getInt(1); // extract raw count 
		result_all.close();
		
		// ---- count where type	
		if (diffSSType)
			query = "SELECT count(*) from "+this.db+".edges where i_res='"+iRes.getAaType().getOneLetterCode()+"' and i_sstype='"
				+iRes.getSsElem().getType()+"' and j_res='"+jRes.getAaType().getOneLetterCode()+"';";					// System.out.println(query);
		else 
			query = "SELECT count(*) from "+this.db+".edges where i_res='"+iRes.getAaType().getOneLetterCode()+"' and j_res='"
			+jRes.getAaType().getOneLetterCode()+"';";	
		result_type = stmt.executeQuery(query);		
		if(result_type.next()) 
			countType = result_type.getInt(1); // extract raw count 
		result_type.close();
			
		// ---- count where angle
		query = "SELECT count(*) from "+db+".edges where r >= "+(r-dRad)+" and r < "+(r+dRad)
		   	+" and theta>="+(theta-dAngle)+" and theta<"+(theta+dAngle)+" and phi>="+(phi-dAngle)+" and phi<"
		   	+(phi+dAngle)+" ;";
		result_angle = stmt.executeQuery(query);
		if(result_angle.next())
			countAngle = result_angle.getInt(1); // extract raw count 
		result_angle.close();
				
		// ---- count where angle and type	
		if (diffSSType)
			query = "SELECT count(*) from "+db+".edges where r >= "+(r-dRad)+" and r < "+(r+dRad)
		    	+" and i_res='"+iRes.getAaType().getOneLetterCode()+"' and i_sstype='"+iRes.getSsElem().getType()+"' and j_res='"+jRes.getAaType().getOneLetterCode()+"' and theta>"
		    	+(theta-dAngle)+" and theta<"+(theta+dAngle)+" and phi>"+(phi-dAngle)+" and phi<"+(phi+dAngle)+" ;";					// System.out.println(query);
		else 
			query = "SELECT count(*) from "+db+".edges where r >= "+(r-dRad)+" and r < "+(r+dRad)
			   +" and i_res='"+iRes.getAaType().getOneLetterCode()+"' and j_res='"+jRes.getAaType().getOneLetterCode()+"' and theta>"+(theta-dAngle)+" and theta<"+(theta+dAngle)+" and phi>"+(phi-dAngle)
			   +" and phi<"+(phi+dAngle)+" ;";					
		result_angle_type = stmt.executeQuery(query);		
		if(result_angle_type.next()) 
			countAngleType = result_angle_type.getInt(1); // extract raw count 
		result_angle_type.close();	
				
		countAll++; countAngle++; countType++; countAngleType++;
				
		countObs = countAngleType;
		countExp = (countAngle*countType)/countAll;
		// -- +factor to avoid division by 0 --> equals countObs++ and countExp++
		double ratio = countObs/countExp;
		double ratio1 = Math.log (ratio);
		lOS=ratio1;
//		System.out.print(ratio+"_"+ratio1+"\t");
		stmt.close();
		
		return lOS;
	}
	
	public double getLogOddsScore(Residue iRes, Residue jRes, double resDist, RIGGeometry graphGeom) throws SQLException{

		double deltar=0.5;
		double deltatheta=Math.PI/72;
		deltatheta /= 2;
//		double deltaphi=Math.PI/72;
		
		double lOS = getLogOddsScore(iRes, jRes, resDist, graphGeom, deltatheta, deltar);
		
		return lOS;
	}
	
	// actual fastest version to compute LOS's
	public void runBayesPreComp() throws SQLException {
		String query;
		Statement stmt;
		ResultSet result_angle, result_type, result_angle_type, result_all;
		double theta=0.0, phi=0.0;
		double countAll = 0, countAngle = 0, countType = 0, countAngleType = 0;
		double countExp = 0;
		double countObs = 0;

		this.db = defaultDB;
		this.tableNameR = this.tableNameR+this.issType+this.iRes+this.jRes;
//		this.tableNameRes = this.tableNameRes+this.issType+this.iRes+this.jRes;
		this.tableNameResR = this.tableNameResR+this.issType+this.iRes+this.jRes;
		this.tableNameResRT = this.tableNameResRT+this.issType+this.iRes+this.jRes;
		this.tableNameRT = this.tableNameRT+this.issType+this.iRes+this.jRes;
						
		this.ratios = new double [this.numRatiosX][this.numRatiosY];
		this.bayesRatios = new double [this.numRatiosX][this.numRatiosY][3];
		
		System.out.println("BayesRatios______extracting voxel density for: "+this.iRes+"_"+this.jRes+" sstype:"+this.issType+ " "+this.diffSSType);		
		
		if (diffSSType)
			tnRes = "edges_"+this.iRes+"_"+String.valueOf(this.issType).toLowerCase()+"_"+this.jRes+"_"+String.valueOf(AnySStype).toLowerCase();
		else
			tnRes = "edges_"+this.iRes+"_"+String.valueOf(AnySStype).toLowerCase()+"_"+this.jRes+"_"+String.valueOf(AnySStype).toLowerCase();	
		
		stmt = conn.createStatement();
		
		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameResRT+";";
		stmt.execute(query);
		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameRT+";";
		stmt.execute(query);
		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameResR+";";
		stmt.execute(query);
		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameR+";";
		stmt.execute(query);
		
		query = "CREATE table "+this.db+"."+this.tableNameR+" (SELECT * from "+this.db+"."+this.tn+" where r >= "+this.minr+" and r < "+this.maxr+");";
		stmt.execute(query);				
		query = "ALTER TABLE "+this.db+"."+this.tableNameR+" add index t(theta);";
		stmt.execute(query);	
//		System.out.println(this.tableNameR+" created");
		
		query = "CREATE table "+this.db+"."+this.tableNameResR+" (SELECT * from "+this.db+"."+this.tnRes+" where r >= "+this.minr+" and r < "+this.maxr+");";
		stmt.execute(query);		
		query = "ALTER TABLE "+this.db+"."+this.tableNameResR+" add index t(theta);";
		stmt.execute(query);		
//		System.out.println(this.tableNameResR+" created");
		
		// ---- count all
		query = "SELECT count(*) from "+this.db+"."+tn+";";
		result_all = stmt.executeQuery(query);
		if(result_all.next()) {
			countAll = result_all.getInt(1); // extract raw count 
		}// and while next results			
		result_all.close();
		
		// ---- count where type
		query = "SELECT count(*) from "+this.db+"."+tnRes+";";
		result_type = stmt.executeQuery(query);
		if(result_type.next()) {
			countType = result_type.getInt(1); // extract raw count 
		}// and while next results	
		result_type.close();
				
		for (int i=0; i<this.ratios.length; i++){
		    theta = this.mintheta + i*this.svoxelsize;
		    
			query = "CREATE table "+this.db+"."+this.tableNameRT+" (SELECT * from "+this.db+"."+this.tableNameR+" where theta>="+theta+" and theta<"+(theta+this.svoxelsize)+");";
			stmt.execute(query);		
			query = "ALTER TABLE "+this.db+"."+this.tableNameRT+" add index p(phi);";
			stmt.execute(query);
			
			query = "CREATE table "+this.db+"."+this.tableNameResRT+" (SELECT * from "+this.db+"."+this.tableNameResR+" where theta>="+theta+" and theta<"+(theta+this.svoxelsize)+");";
			stmt.execute(query);			
			query = "ALTER TABLE "+this.db+"."+this.tableNameResRT+" add index p(phi);";
			stmt.execute(query);

		   	for (int j=0; j<this.ratios[i].length; j++){
	    		phi = this.minphi + j*this.svoxelsize;
				
				// count where angle
				query = "SELECT count(*) from "+db+"."+this.tableNameRT+" where phi>="+phi+" and phi<"+(phi+this.svoxelsize)+" ;";
				result_angle = stmt.executeQuery(query);
				if(result_angle.next()) {
					countAngle = result_angle.getInt(1); // extract raw count 
				}// and while next results	
				result_angle.close();
				
				// count where angle and type
				query = "SELECT count(*) from "+db+"."+this.tableNameResRT+" where phi>="+phi+" and phi<"+(phi+this.svoxelsize)+" ;";
				result_angle_type = stmt.executeQuery(query);
				if(result_angle_type.next()) {
					countAngleType = result_angle_type.getInt(1); // extract raw count 
				}// and while next results	
				result_angle_type.close();	
				
				countAll++; countAngle++; countType++; countAngleType++;
				
				countObs = countAngleType;
				countExp = (countAngle*countType)/countAll;
				// -- +factor to avoid division by 0 --> equals countObs++ and countExp++
				double ratio = countObs/countExp;
				double ratio1 = Math.log (ratio);
								
				//---replace SoP with array entry
				this.ratios[i][j] = ratio1;
				this.bayesRatios[i][j][0] = ratio1;
				this.bayesRatios[i][j][1] = countObs;
				this.bayesRatios[i][j][2] = countExp;
				if (ratio1<this.minRatio)
					this.minRatio = ratio1;
				if (ratio1>this.maxRatio)
					this.maxRatio = ratio1;
					
			} //--end for phi
			
			query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameRT+";";
			stmt.execute(query);
			query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameResRT+";";
			stmt.execute(query);
		} //--end for theta

		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameResR+";";
		stmt.execute(query);
		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameR+";";
		stmt.execute(query);
		stmt.close();
		
	}
	
	public void runBayes() throws SQLException {
		String query;
		Statement stmt;
		ResultSet result_angle, result_type, result_angle_type, result_all;
		double theta=0.0, phi=0.0;
		double countAll = 0, countAngle = 0, countType = 0, countAngleType = 0;
		double countExp = 0;
		double countObs = 0;
		
		this.db = defaultDB;
		String tn = "edges";
		String tnRes = "";
//		String tnResR = "";
//		String tnR = "";
		// initialize table names
//		tnR = tn+"_"+String.valueOf(this.minr)+"_"+String.valueOf(this.maxr);
		if (diffSSType)
			tnRes = "edges_"+this.iRes+"_"+String.valueOf(this.issType).toLowerCase()+"_"+this.jRes+"_"+String.valueOf(this.jssType).toLowerCase();
		else
			tnRes = "edges_"+this.iRes+"_"+String.valueOf(AnySStype).toLowerCase()+"_"+this.jRes+"_"+String.valueOf(AnySStype).toLowerCase();
//		tnResR = tnRes + "_" +String.valueOf(this.minr)+"_"+String.valueOf(this.maxr);
		
		this.ratios = new double [this.numRatiosX][this.numRatiosY];
		this.bayesRatios = new double [this.numRatiosX][this.numRatiosY][3];
		
		System.out.println("BayesRatios______extracting voxel density for: "+this.iRes+"_"+this.jRes+" sstype:"+this.issType+ " "+this.diffSSType);		
		
		stmt = conn.createStatement();
			
		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameResRT+";";
		stmt.execute(query);
		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameRT+";";
		stmt.execute(query);
		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameResR+";";
		stmt.execute(query);
		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameR+";";
		stmt.execute(query);
		
		query = "CREATE table "+this.db+"."+this.tableNameR+" (SELECT * from "+this.db+"."+tn+" where r >= "+this.minr+" and r < "+this.maxr+");";
		stmt.execute(query);		
		query = "ALTER TABLE "+this.db+"."+this.tableNameR+" add index t(theta);";
		stmt.execute(query);
		
		query = "CREATE table "+this.db+"."+this.tableNameResR+" (SELECT * from "+this.db+"."+tnRes+" where r >= "+this.minr+" and r < "+this.maxr+");";
		stmt.execute(query);		
		query = "ALTER TABLE "+this.db+"."+this.tableNameResR+" add index t(theta);";
		stmt.execute(query);	
		
		// ---- count all
		query = "SELECT count(*) from "+this.db+"."+tn+";";
		result_all = stmt.executeQuery(query);
		if(result_all.next()) {
			countAll = result_all.getInt(1); // extract raw count 
		}// and while next results			
		result_all.close();
		
		// ---- count where type
		query = "SELECT count(*) from "+this.db+"."+tnRes+";";
		result_type = stmt.executeQuery(query);
		if(result_type.next()) {
			countType = result_type.getInt(1); // extract raw count 
		}// and while next results	
		result_type.close();
		
		for (int i=0; i<this.ratios.length; i++){
		    theta = this.mintheta + i*this.svoxelsize;
		    
			query = "CREATE table "+this.db+"."+this.tableNameRT+" (SELECT * from "+this.db+"."+this.tableNameR+" where theta>="+theta+" and theta<"+(theta+this.svoxelsize)+");";
			stmt.execute(query);		
			query = "ALTER TABLE "+this.db+"."+this.tableNameRT+" add index p(phi);";
			stmt.execute(query);
			
			query = "CREATE table "+this.db+"."+this.tableNameResRT+" (SELECT * from "+this.db+"."+this.tableNameResR+" where theta>="+theta+" and theta<"+(theta+this.svoxelsize)+");";
			stmt.execute(query);			
			query = "ALTER TABLE "+this.db+"."+this.tableNameResRT+" add index p(phi);";
			stmt.execute(query);

		   	for (int j=0; j<this.ratios[i].length; j++){
	    		phi = this.minphi + j*this.svoxelsize;
				
				// count where angle
				query = "SELECT count(*) from "+db+"."+this.tableNameRT+" where phi>="+phi+" and phi<"+(phi+this.svoxelsize)+" ;";
				result_angle = stmt.executeQuery(query);
				if(result_angle.next()) {
					countAngle = result_angle.getInt(1); // extract raw count 
				}// and while next results	
				result_angle.close();
				
				// count where angle and type
				query = "SELECT count(*) from "+db+"."+this.tableNameResRT+" where phi>="+phi+" and phi<"+(phi+this.svoxelsize)+" ;";
				result_angle_type = stmt.executeQuery(query);
				if(result_angle_type.next()) {
					countAngleType = result_angle_type.getInt(1); // extract raw count 
				}// and while next results	
				result_angle_type.close();	
				
				countAll++; countAngle++; countType++; countAngleType++;
				
				countObs = countAngleType;
				countExp = (countAngle*countType)/countAll;
				// -- +factor to avoid division by 0 --> equals countObs++ and countExp++
				double ratio = countObs/countExp;
				double ratio1 = Math.log (ratio);
								
				//---replace SoP with array entry
				this.ratios[i][j] = ratio1;
				this.bayesRatios[i][j][0] = ratio1;
				this.bayesRatios[i][j][1] = countObs;
				this.bayesRatios[i][j][2] = countExp;
				if (ratio1<this.minRatio)
					this.minRatio = ratio1;
				if (ratio1>this.maxRatio)
					this.maxRatio = ratio1;
					
			} //--end for phi
			
			query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameRT+";";
			stmt.execute(query);
			query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameResRT+";";
			stmt.execute(query);
		} //--end for theta
		stmt.close();
	}
	
	public void runBayesByCreateTables() throws SQLException {
		String query;
		Statement stmt;
		ResultSet result_angle, result_type, result_angle_type, result_all;
		double theta=0.0, phi=0.0;
		double countAll = 0, countAngle = 0, countType = 0, countAngleType = 0;
		double countExp = 0;
		double countObs = 0;
		
		this.db = defaultDB;
		
		this.ratios = new double [this.numRatiosX][this.numRatiosY];
		this.bayesRatios = new double [this.numRatiosX][this.numRatiosY][3];
		
		System.out.println("BayesRatios______extracting voxel density for: "+this.iRes+"_"+this.jRes+" sstype:"+this.issType+ " "+this.diffSSType);		
		
		stmt = conn.createStatement();
		
		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameRes+";";
		stmt.execute(query);	
		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameResRT+";";
		stmt.execute(query);
		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameRT+";";
		stmt.execute(query);
		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameResR+";";
		stmt.execute(query);
		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameR+";";
		stmt.execute(query);
				
		if (diffSSType){
			query = "CREATE table "+this.db+"."+this.tableNameRes+" (SELECT * from "+this.db+"."+this.tn+" where i_res='"+this.iRes+"' and i_sstype='"
				+issType+"' and j_res='"+this.jRes+"');";
		}
		else{
			query = "CREATE table "+this.db+"."+this.tableNameRes+" (SELECT * from "+this.db+"."+this.tn+" where i_res='"+this.iRes+"' and j_res='"
				+this.jRes+"');";
		}		
		stmt.execute(query);		
		query = "ALTER TABLE "+this.db+"."+this.tableNameRes+" add index r(r);";
		stmt.execute(query);
		
		query = "CREATE table "+this.db+"."+this.tableNameR+" (SELECT * from "+this.db+"."+this.tn+" where r >= "+this.minr+" and r < "+this.maxr+");";
		stmt.execute(query);		
		query = "ALTER TABLE "+this.db+"."+this.tableNameR+" add index t(theta);";
		stmt.execute(query);
		
		query = "CREATE table "+this.db+"."+this.tableNameResR+" (SELECT * from "+this.db+"."+this.tableNameRes+" where r >= "+this.minr+" and r < "+this.maxr+");";
		stmt.execute(query);		
		query = "ALTER TABLE "+this.db+"."+this.tableNameResR+" add index t(theta);";
		stmt.execute(query);	
		
		// ---- count all
		query = "SELECT count(*) from "+this.db+"."+this.tn+";";
		result_all = stmt.executeQuery(query);
		if(result_all.next()) {
			countAll = result_all.getInt(1); // extract raw count 
		}// and while next results			
		result_all.close();
		
		// ---- count where type
		query = "SELECT count(*) from "+this.db+"."+this.tableNameRes+";";
		result_type = stmt.executeQuery(query);
		if(result_type.next()) {
			countType = result_type.getInt(1); // extract raw count 
		}// and while next results	
		result_type.close();
		
		for (int i=0; i<this.ratios.length; i++){
		    theta = this.mintheta + i*this.svoxelsize;
		    
			query = "CREATE table "+this.db+"."+this.tableNameRT+" (SELECT * from "+this.db+"."+this.tableNameR+" where theta>="+theta+" and theta<"+(theta+this.svoxelsize)+");";
			stmt.execute(query);		
			query = "ALTER TABLE "+this.db+"."+this.tableNameRT+" add index p(phi);";
			stmt.execute(query);
			
			query = "CREATE table "+this.db+"."+this.tableNameResRT+" (SELECT * from "+this.db+"."+this.tableNameResR+" where theta>="+theta+" and theta<"+(theta+this.svoxelsize)+");";
			stmt.execute(query);			
			query = "ALTER TABLE "+this.db+"."+this.tableNameResRT+" add index p(phi);";
			stmt.execute(query);

		   	for (int j=0; j<this.ratios[i].length; j++){
	    		phi = this.minphi + j*this.svoxelsize;
				
				// count where angle
				query = "SELECT count(*) from "+db+"."+this.tableNameRT+" where phi>="+phi+" and phi<"+(phi+this.svoxelsize)+" ;";
				result_angle = stmt.executeQuery(query);
				if(result_angle.next()) {
					countAngle = result_angle.getInt(1); // extract raw count 
				}// and while next results	
				result_angle.close();
				
				// count where angle and type
				query = "SELECT count(*) from "+db+"."+this.tableNameResRT+" where phi>="+phi+" and phi<"+(phi+this.svoxelsize)+" ;";
				result_angle_type = stmt.executeQuery(query);
				if(result_angle_type.next()) {
					countAngleType = result_angle_type.getInt(1); // extract raw count 
				}// and while next results	
				result_angle_type.close();	
				
				countAll++; countAngle++; countType++; countAngleType++;
				
				countObs = countAngleType;
				countExp = (countAngle*countType)/countAll;
				// -- +factor to avoid division by 0 --> equals countObs++ and countExp++
				double ratio = countObs/countExp;
				double ratio1 = Math.log (ratio);
								
				//---replace SoP with array entry
				this.ratios[i][j] = ratio1;
				this.bayesRatios[i][j][0] = ratio1;
				this.bayesRatios[i][j][1] = countObs;
				this.bayesRatios[i][j][2] = countExp;
				if (ratio1<this.minRatio)
					this.minRatio = ratio1;
				if (ratio1>this.maxRatio)
					this.maxRatio = ratio1;
					
			} //--end for phi
			
			query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameRT+";";
			stmt.execute(query);
			query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameResRT+";";
			stmt.execute(query);
		} //--end for theta

		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameRes+";";
		stmt.execute(query);	
		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameResR+";";
		stmt.execute(query);
		query = "DROP TABLE IF EXISTS "+db+"."+this.tableNameR+";";
		stmt.execute(query);
		stmt.close();
	}
	
	public void runBayesSlow() throws SQLException {
		String query;
		Statement stmt;
		ResultSet result_angle, result_type, result_angle_type, result_all;
		double theta=0.0, phi=0.0;
		double countAll = 0, countAngle = 0, countType = 0, countAngleType = 0;
		double countExp = 0;
		double countObs = 0;
		
		this.ratios = new double [this.numRatiosX][this.numRatiosY];
		this.bayesRatios = new double [this.numRatiosX][this.numRatiosY][3];
		
		System.out.println("BayesRatios______extracting voxel density for: "+this.iRes+"_"+this.jRes+" sstype:"+this.issType+ " "+this.diffSSType);	
//		System.out.println("countAll: countAngle: countType: countAngleType:");	
		
		// ---- count all
		stmt = conn.createStatement();
		query = "SELECT count(*) from "+db+".edges;";
		result_all = stmt.executeQuery(query);
		if(result_all.next()) {
			countAll = result_all.getInt(1); // extract raw count 
		}// and while next results				
		stmt.close();	
		result_all.close();
		// ---- count where type
		stmt = conn.createStatement();		
		if (diffSSType){
			query = "SELECT count(*) from "+db+".edges where i_res='"+this.iRes+"' and i_sstype='"
			+issType+"' and j_res='"+this.jRes+"';";					// System.out.println(query);
		}
		else {
			query = "SELECT count(*) from "+db+".edges where i_res='"+this.iRes+"' and j_res='"
			+this.jRes+"';";								
		}
		result_type = stmt.executeQuery(query);
		if(result_type.next()) {
			countType = result_type.getInt(1); // extract raw count 
		}// and while next results	
		stmt.close();
		result_type.close();
		
		for (int i=0; i<this.ratios.length; i++){
		    theta = this.mintheta + i*this.svoxelsize;
		   	for (int j=0; j<this.ratios[i].length; j++){
	    		phi = this.minphi + j*this.svoxelsize;
				
				// count where angle
				stmt = conn.createStatement();
				query = "SELECT count(*) from "+db+".edges where r >= "+this.minr+" and r < "+this.maxr
			    +" and theta>="+theta+" and theta<"+(theta+this.svoxelsize)+" and phi>="+phi+" and phi<"
			    +(phi+this.svoxelsize)+" ;";
				result_angle = stmt.executeQuery(query);
				if(result_angle.next()) {
					countAngle = result_angle.getInt(1); // extract raw count 
				}// and while next results				
				stmt.close();	
				result_angle.close();
				
				// count where angle and type
				stmt = conn.createStatement();		
				if (diffSSType){
					query = "SELECT count(*) from "+db+".edges where r >= "+this.minr+" and r < "+this.maxr
				    +" and i_res='"+this.iRes+"' and i_sstype='"+issType+"' and j_res='"+this.jRes+"' and theta>"
				    +theta+" and theta<"+(theta+this.svoxelsize)+" and phi>"+phi+" and phi<"+(phi+this.svoxelsize)+" ;";					// System.out.println(query);
				}
				else {
					query = "SELECT count(*) from "+db+".edges where r >= "+this.minr+" and r < "+this.maxr
				    +" and i_res='"+this.iRes+"' and j_res='"+this.jRes+"' and theta>"+theta+" and theta<"+(theta+this.svoxelsize)+" and phi>"+phi
				    +" and phi<"+(phi+this.svoxelsize)+" ;";					
				}
				result_angle_type = stmt.executeQuery(query);
				if(result_angle_type.next()) {
					countAngleType = result_angle_type.getInt(1); // extract raw count 
				}// and while next results				
				stmt.close();
				result_angle_type.close();	
				
				countAll++; countAngle++; countType++; countAngleType++;
				
				countObs = countAngleType;
				countExp = (countAngle*countType)/countAll;
				// -- +factor to avoid division by 0 --> equals countObs++ and countExp++
				double ratio = countObs/countExp;
				double ratio1 = Math.log (ratio);
			
				this.ratios[i][j] = ratio1;
				this.bayesRatios[i][j][0] = ratio1;
				this.bayesRatios[i][j][1] = countObs;
				this.bayesRatios[i][j][2] = countExp;
				if (ratio1<this.minRatio)
					this.minRatio = ratio1;
				if (ratio1>this.maxRatio)
					this.maxRatio = ratio1;
			} //--end for phi
		} //--end for theta
	}		
	
	public void writeSphoxelOutput(String filename){
		String[] names = {"log(ratio)","Obs","Exp"};
		CSVhandler csv = new CSVhandler();
		csv.generateCsvFile(this.bayesRatios, names, filename);
	}

	public void setDb(String db) {
		this.db = db;
	}
	public String getDb() {
		return this.db;
	}
	
	public void setMinr(float minr){
		this.minr = minr;
	}
	public double getMinr(){
		return this.minr;
	}
	public void setMaxr(float maxr){
		this.maxr = maxr;
	}
	public double getMaxr(){
		return this.maxr;
	}
	
	public void setIRes(char c){
		this.iRes = c;
		System.out.println("iRes: "+c);
	}
	public char getIRes(){
		return this.iRes;
	}
	
	public void setJRes(char c){
		this.jRes = c;
	}
	public char getJRes(){
		return this.jRes;
	}	

	public void setISSType(char c){
		this.issType = c;
	}
	public char getISSType(){
		return this.issType;
	}
	public void setJSSType(char c){
		this.jssType = c;
	}
	public char getJSSType(){
		return this.jssType;
	}
	
	public void setDiffSSType(boolean b){
		this.diffSSType = b;
	}
	public boolean getDiffSSType(){
		return this.diffSSType;
	}
	
	public void setNumSteps(int num){
		this.numSteps = num;
		this.numRatiosX = this.numSteps;
		this.numRatiosY = 2*this.numSteps;
		this.svoxelsize = Math.PI/this.numSteps;
	}
	public int getNumSteps(){
		return this.numSteps;
	}
	public int getNumRatiosX(){
		return this.numRatiosX;
	}
	public int getNumRatiosY(){
		return this.numRatiosY;
	}
	public double getSVoxelsize(){
		return this.svoxelsize;
	}
	public double getMinRatio(){
		return this.minRatio;
	}
	public double getMaxRatio(){
		return this.maxRatio;
	}
	public double[][] getRatios(){
		return this.ratios;
	}
	public double[][][] getBayesRatios(){
		return this.bayesRatios;
	}
	public int getDefaultNumSteps(){
		return defaultNumSteps;
	}
	public float getDefaultResol(){
		return defaultResol;
	}
	public String getRadiusPrefix() {
		return radiusPrefix;
	}
	public void setRadiusPrefix(String radiusPrefix) {
		this.radiusPrefix = radiusPrefix;
	}


}
