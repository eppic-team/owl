package tools;

import tools.*;
import java.sql.*;
import java.io.*;

public class Report {

    private Connection myConnection;
    private mySQLConnect SQLC;
    private PrintStream Log;
    private int reportLevel;
    private String hostName;
    private String clientName;
    
    public Report(String connFile, String logFile, int reportLevel) {

	SQLC = new mySQLConnect();
	SQLC.readConnectionFile(connFile);
	myConnection = SQLC.openConnection();
	
	try {
	    Log = new PrintStream(new FileOutputStream(logFile));
	} catch (Exception e) {
	    System.out.println(e);
	}

	this.reportLevel = reportLevel;
	clientName = Machine.getClient();
	hostName = SQLC.getHost();
	
    }

    /** Sets the reportLevel, so that only messages below or equal the level given are displayed.
	The following levels are available: (0 is default - no messages) 
	3 - talkative mode 
	2 - normal mode 
	1 - basic mode 
	0 - silent mode
	This means by a higher level you get more detailed information, while on a reportLevel 
	of 0 only the most basic results and errors are displayed */
    public void setReportLevel(int newLevel) { reportLevel = newLevel; }

    public int getReportLevel() { return reportLevel; }

    public void print(String text, int level) {
	if (level<=reportLevel) {
	    System.out.print(text);
	} // end if level of the curent message is above the current level 
    } // end of report

    public void println(String text, int level) {
	if (level<=reportLevel) {
	     System.out.println(text);
	} // end if level of the curent message is above the current level 
    } // end of report

    public void printLog(String text, int level) {
	if (level<=reportLevel) {
	    Log.print(text);
	} // end if level of the curent message is above the current level 
    } // end of report

    public void printlnLog(String text, int level) {
	if (level<=reportLevel) {
	     Log.println(text);
	} // end if level of the curent message is above the current level 
    } // end of report

    public void newProc(String description, int level) {
	
	if (level<=reportLevel) {
	    try {
		Statement S = myConnection.createStatement();
		S.executeUpdate("INSERT INTO report (Client, Host, Description, Start) VALUES (\""+clientName+"\", \""+hostName+"\", \""+description+"\", NOW());");
		S.close();
	    } catch (SQLException E) {
		System.out.println("SQLException:\t" + E.getMessage());
		System.out.println("SQLState:\t" + E.getSQLState());
		System.out.println("VendorError: \t" + E.getErrorCode());
	    }
	}
	    
    }

    public void updateProc(int level) {

	if (level<=reportLevel) {
	    try {
		Statement S = myConnection.createStatement();
		S.executeUpdate("UPDATE report SET End = NOW(), Total = TIMEDIFF(End, Start) WHERE ID = LAST_INSERT_ID();");
		S.close();
	    } catch (SQLException E) {
		System.out.println("SQLException:\t" + E.getMessage());
		System.out.println("SQLState:\t" + E.getSQLState());
		System.out.println("VendorError: \t" + E.getErrorCode());
	    }
	}
 
    }

    public void updateProc(String description, int level) {
	
	if (level<=reportLevel) {
	    try {
		Statement S = myConnection.createStatement();
		S.executeUpdate("UPDATE report SET End = NOW(), Total = TIMEDIFF(End, Start) WHERE (Client = \""+clientName+"\") AND (Host = \""+hostName+"\") AND (Description = \""+description+"\");");
		S.close();
	    } catch (SQLException E) {
		System.out.println("SQLException:\t" + E.getMessage());
		System.out.println("SQLState:\t" + E.getSQLState());
		System.out.println("VendorError: \t" + E.getErrorCode());
	    }
	}

    }

    public void close() {
	
	SQLC.closeConnection(myConnection);
	Log.close();

    }

}
