import java.sql.*;

import tools.ClusterConnection;

public class testClusterConnection {

	/**
	 * @param args
	 */
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int asuId=30;
		//ClusterConnection conn=new ClusterConnection();
		ClusterConnection conn=new ClusterConnection("pdbgraph","asu","duarte","nieve");
		System.out.println(conn.db);
		try {
			//System.out.println(conn.getHost4Idx(asuId));		
			Statement S;
			ResultSet R;
			//S = conn.createNConStatement("asu",asuId);
			S=conn.createStatement(asuId);
			System.out.println(conn.host);
			R=S.executeQuery("SELECT asu_id,accession_code FROM asu_list WHERE asu_id="+asuId+";");
			System.out.println("asu_id\taccession_code");
			System.out.println("Fetch size is="+R.getFetchSize());
			while (R.next()){
				System.out.println(R.getInt(1)+"\t"+R.getString(2));
			}
			R.close();
			S.close();
		}
		catch (SQLException e) {
			e.printStackTrace();
		}
		conn.close();
		/*
		try {
		Runtime.getRuntime().exec("touch jose");
		}
		catch(IOException e){
			e.printStackTrace();
		}
		*/
	}

}
