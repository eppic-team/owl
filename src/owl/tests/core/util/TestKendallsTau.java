package owl.tests.core.util;

import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Vector;

import owl.core.util.KendallsTau;
import owl.core.util.MySQLConnection;
import owl.gmbp.CSVhandler;



public class TestKendallsTau {	
	

	private static String host = "talyn";
	private static String username = "vehlow";
	private static String password = "nieve";
	private static String db = "scoring";
	
	private static MySQLConnection conn;
	
	private static double thresGDT = 0.0;

	/**
	 * @param args
	 * @throws SQLException 
	 */
	public static void main(String[] args) {

//		TestKendallsTau kdTest = new TestKendallsTau();
		
//		try {
////			runByDB1();
//			combineTables();
//		} catch (SQLException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}	
		
//		thresGDT = 80.0;
		evaluateScoreList();

	}
	
	private static void evaluateScoreList(){
		CSVhandler csv = new CSVhandler();
		String inputDir = "/Volumes/StruPPi/CASP8/server_models/";
		inputDir = "/Users/vehlow/Documents/workspace/outputFiles/";
		String filename = csv.openCSVFile(inputDir);
		System.out.println(filename);
		Vector<String[]> vec=null;
		try {
			vec = csv.readCsvFile(filename);
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if (vec!=null){
			int num = vec.size()-1;
			int[] col1 = new int[num]; 
			String[] col0 = new String[num]; 
			double[] col2 = new double[num];
			double[] col3 = new double[num];
			
			int cnt=0;
			for (String[] row:vec){
				if(vec.indexOf(row)==0 && row[0].equals("target")){	}
				else {
					double rank = Double.valueOf(row[2]);
					int r = (int) rank;
					if (r>0){
						if (!(row[4].equals("NaN")) && Double.valueOf(row[3])>=thresGDT){
//						if (!(row[4].equals("NaN") || row[4].equals("0")) && Double.valueOf(row[3])>=thresGDT){
							col0[cnt] = row[1];
							col1[cnt] = r;
							col2[cnt] = Double.valueOf(row[3]);
							col3[cnt] = Double.valueOf(row[4]);
							cnt++;
//							System.out.println(row[0]+" "+r);												
						}
					}
				}
			}
			int[] rank = new int[cnt]; 
			double[] ranking = new double[cnt];
			String[] server = new String[cnt]; 
			double[] gdts = new double[cnt];
			double[] geom = new double[cnt];
			if (cnt<num){
				for (int i=0; i<cnt; i++){
					server[i] = col0[i];
					rank[i] = col1[i];
					ranking[i] = (double)col1[i];
					gdts[i] = Double.valueOf(col2[i]);
					geom[i] = Double.valueOf(col3[i]);					
				}				
			}
			else {
				server = col0;
				rank = col1;
				for (int i=0; i<cnt; i++)
					ranking[i] = col1[i];
				gdts = col2;
				geom = col3;
			}
//			for (int i=0; i<cnt; i++){
//				System.out.print(server[i]+"\t"+rank[i]+"\t"+gdts[i]+"\t"+geom[i]+"\n");
//			}
			
			// RCC's:
			System.out.println("TauB's");
			System.out.println("for all scores with GDT>="+thresGDT);
			double tB=0;
			tB = computeKendallsTau(gdts, ranking, rank);
			System.out.println("gdts vs ranking,"+tB); 
			tB = computeKendallsTau(gdts, geom, rank);
			System.out.println("gdts vs geom,"+tB); 
		}
	}
	
	@SuppressWarnings("unused")
	private static void combineTables() throws SQLException{
		conn = new MySQLConnection(host,username,password,db);	
		ResultSet result;
		String query;
		Statement stmt;
		
		int num =0;
				
		stmt = conn.createStatement();	
		
		query = "Select count(*) from scoring.t471_gdts;";
		result = stmt.executeQuery(query);
		while (result.next()){
			num = result.getInt(1);
//			System.out.println("Num= "+num); 
		}
		
		int[] rank = new int[num]; 
		double[] ranking = new double[num];
		double[] antiRanking = new double[num];
		String[] server = new String[num]; 
		double[] gdts = new double[num]; 
		double[] dopescore = new double[num]; 
		double[] roscore = new double[num]; 
		double[] fascore = new double[num]; 		
		
		query = "Select * from scoring.t471_gdts group by rank;";
//		System.out.println(query);
		result = stmt.executeQuery(query);
		
		int cnt = 0;
		while (result.next()){
			server[cnt] = result.getString(1);
			rank[cnt] = result.getInt(2);
			gdts[cnt] = result.getDouble(3);
			
			cnt++;
		}
		
		for (int i=0; i<server.length; i++){
			query = "Select dopescore from scoring.t471_dope where server like '"+server[i]+"';";
			result = stmt.executeQuery(query);
			while (result.next()){
				dopescore[i] = result.getDouble(1);
			}
			query = "Select roscore, fascore from scoring.t471_rosetta where server like '"+server[i]+"';";
			result = stmt.executeQuery(query);
			while (result.next()){
				roscore[i] = result.getDouble(1);
				fascore[i] = result.getDouble(2);
			}
			ranking[i] = (double)rank[i];
			antiRanking[i] = (double)rank[rank.length-i-1];
		}
		
		result.close();
		stmt.close();

		System.out.println("server"+","+"rank"+","+"gdts"+","+"dopeS"+","+"rosS"+","+"faS");
		for (int i=0; i<server.length; i++)
			System.out.println(server[i]+","+rank[i]+","+gdts[i]+","+dopescore[i]+","+roscore[i]+","+fascore[i]);
		System.out.println();
		
		// RCC's:
		System.out.println("TauB's");
		double tB=0;
		tB = computeKendallsTau(dopescore, roscore, rank);
		System.out.println("dopeS vs rosetta,"+tB); 
		tB = computeKendallsTau(gdts, ranking, rank);
		System.out.println("gdts vs ranking,"+tB);
		tB = computeKendallsTau(gdts, antiRanking, rank);
		System.out.println("gdts vs antiRanking,"+tB);  
		tB = computeKendallsTau(fascore, roscore, rank);
		System.out.println("fascore vs roscore,"+tB);  
		tB = computeKendallsTau(gdts, roscore, rank);
		System.out.println("gdts vs roscore,"+tB);
		tB = computeKendallsTau(gdts, dopescore, rank);
		System.out.println("gdts vs dopescore,"+tB);
		
	}
	
	@SuppressWarnings("unused")
	private static void runByDB1() throws SQLException{		
		conn = new MySQLConnection(host,username,password,db);	
		ResultSet result;
		String query;
		Statement stmt;
		
		int num =0;
		boolean reverse = true;

		stmt = conn.createStatement();	
		query = "Select count(*) from scoring.t531_exscores;";
		result = stmt.executeQuery(query);
		while (result.next()){
			num = result.getInt(1);
			System.out.println("Num= "+num); 
		}
//		query = "Select * from scoring.t531_exscores group by dopescore;";
//		query = "Select * from scoring.t531_exscores group by roscore;";
		query = "Select * from scoring.t531_exscores group by fascore;";
		System.out.println(query+"    "+String.valueOf(reverse));
		result = stmt.executeQuery(query);
		
		int ranking[] = new int[num]; 
		double rank[] = new double[num];
		for (int i=0; i<ranking.length; i++){
			if (reverse){
				ranking[i] = ranking.length-i;
				rank[i] = ranking.length-i;
			}
			else {
				ranking[i] = i+1;
				rank[i] = i+1;				
			}
		}
		double dopeS[] = new double[num];
		double rosS[] = new double[num];
		double faS[] = new double[num];
		String keys[] = new String[num];
		
		int cnt = 0;
		System.out.print("ranking"+"\t"+"keys"+"\t"+"dopeS"+"\t"+"rosS"+"\t"+"faS"+"\n");
		while (result.next()){
			keys[cnt] = result.getString(1);
			dopeS[cnt] = result.getDouble(2);
			rosS[cnt] = result.getDouble(3);
			faS[cnt] = result.getDouble(4);
			System.out.print(ranking[cnt]+"\t"+keys[cnt]+"\t"+dopeS[cnt]+"\t"+rosS[cnt]+"\t"+faS[cnt]+"\n");
			cnt++;
		}
		
		stmt.close();
		result.close();		
		
		double tB = 0;
		
		tB = computeKendallsTau(dopeS, rank, ranking);
		System.out.println("dopeS vs rank    tau B : "+tB); 
		tB = computeKendallsTau(rosS, rank, ranking);
		System.out.println("rosS vs rank     tau B : "+tB); 
		tB = computeKendallsTau(faS, rank, ranking);
		System.out.println("faS vs rank      tau B : "+tB); 
		tB = computeKendallsTau(dopeS, rosS, ranking);
		System.out.println("dopeS vs rosS    tau B : "+tB); 
		tB = computeKendallsTau(dopeS, faS, ranking);
		System.out.println("dopeS vs faS     tau B : "+tB); 
		tB = computeKendallsTau(rosS, faS, ranking);
		System.out.println("rosS vs faS      tau B : "+tB); 
			
	}
	
	@SuppressWarnings("unused")
	private double computeKendallsTau(KendallsTau kt1, KendallsTau kt2){		
		return kt1.kendallsTauB( kt2);
	}
	
	private static double computeKendallsTau(double[] values1, double[] values2, int[] ranking){
		KendallsTau kt1 = new KendallsTau(values1, ranking);
		KendallsTau kt2 = new KendallsTau(values2, ranking);
//		for (int i=0; i<values1.length; i++)
//			System.out.println(kt1.getRanking()[i]+" "+kt1.getValues()[i]+"   "+kt2.getRanking()[i]+" "+kt2.getValues()[i]);
//		int xChanges = kt1.bubbleSort();
//		System.out.println("nr of Exchanges total kt1: "+xChanges);
//		xChanges = kt2.bubbleSort();
//		System.out.println("nr of Exchanges total kt2: "+xChanges);
//		for (int i=0; i<values1.length; i++)
//			System.out.println(kt1.getRanking()[i]+" "+kt1.getValues()[i]+"   "+kt2.getRanking()[i]+" "+kt2.getValues()[i]);
		
//		// compare the two RVs 
//		double ktVal = kt1.kendallsTau(kt2);
//		System.out.println("normalized Kendalls tau = (x*2)/((n-1)*n) = "+ktVal); 
//		for (int i=0; i<values1.length; i++)
//			System.out.println(kt1.getRanking()[i]+" "+kt1.getValues()[i]+"   "+kt2.getRanking()[i]+" "+kt2.getValues()[i]);
				
//		double tB= kt1.kendallsTauB( kt2);
//		System.out.println("calculating Kendalls tau B : "+tB); 
//		for (int i=0; i<values1.length; i++)
//			System.out.println(kt1.getRanking()[i]+" "+kt1.getValues()[i]+"   "+kt2.getRanking()[i]+" "+kt2.getValues()[i]);
//		return tB;
		
		return kt1.kendallsTauB( kt2);
	}

}
