import java.sql.*;
import java.io.*;
import tools.*;

public class testGraph2Pml {

	public static void main(String[] args) {

		String connFile = "/project/StruPPi/ioannis/cnfs/msdgraph.my.cnf";
	
		mySQLConnect SQLC =  new mySQLConnect();
		SQLC.readConnectionFile(connFile);
		Connection conn = SQLC.openConnection();
	
		PrintWriter serverOutPw = new PrintWriter(new PymolServerOutputStream("http://blau:9123"), true);
		PyMol pml = new PyMol(serverOutPw);
	
		String pdbFileName = Msdsd2Pdb.export2File("1rx4", 20717, 52567, 9, "/project/StruPPi/ioannis/tmp");
		String molObjName = pml.loadPDB(pdbFileName, "/project/StruPPi/ioannis/tmp");
		System.out.println(molObjName);
	
		Graph2Pml graphPml = new Graph2Pml(serverOutPw, molObjName, 33729, "SC_SC", "SC_SC_in+SC_SC_out", "SC_SC", "true", "true", true, false, false, conn);
		graphPml.exportPML();
	
		SQLC.closeConnection(conn);	

    }

} // end of class Graph2Pml
