package tools;

import java.sql.SQLException;

public class QueryThread extends Thread {

	String query;
	String host;
	String user;
	String pwd;
	String db;
	
	public QueryThread(String query, String host, String user, String pwd, String db){
		this.query = query;
		this.host = host;
		this.user = user;
		this.pwd = pwd;
		this.db = db;
		start();
	}
	
	public QueryThread(String query, MySQLConnection conn){
		this.query = query;
		this.host = conn.getHost();
		this.user = conn.getUser();
		this.pwd = conn.getPassword();
		this.db = conn.getDbname();
		start();		
	}
	
	public void run(){
		try {
			long start = System.currentTimeMillis();
			MySQLConnection conn = new MySQLConnection(host,user,pwd,db);
			conn.executeSql(query);
			conn.close();
			long end = System.currentTimeMillis();
			System.out.println("Executed QUERY="+query+" DB="+db+", HOST="+host+". Time was: "+(end-start)/1000+" seconds.");
		}
		catch(SQLException e){
			e.printStackTrace();
			System.err.println("Couldn't execute QUERY="+query+" DB="+db+", HOST="+host);
		}
	}
	
}
