package tools; 

import tools.*;
import java.io.*;
import java.util.*;
import org.rosuda.JRclient.*;
import Jama.*;

public class R {

    private static Rconnection c = null;
    private static String mySqlConfString = null;

    public R(String username, String password, String connFile) {

	try {
	    c = new Rconnection();
	    c.login(username, password);
	    mySqlConfString = readConnectionFile(connFile);
	    // treat warnings as errors
	    c.voidEval("options(warn=0)");
	    c.voidEval("library(RMySQL);drv<-dbDriver(\"MySQL\")");
	} catch(RSrvException rse) {
            System.out.println("Rserve exception: "+rse.getMessage());
	    cleanUp();
        } catch(Exception e) {
            System.out.println("Something went wrong, but it's not the Rserve: "+e.getMessage());
            e.printStackTrace();
	    cleanUp();
        }
	
    }

    public double[][] simpleLinRegr(String predictorQuery, String responseQuery) {
	
	//predictor = intercept + slope*response
	
	double[][] lr = null;
	boolean singlePredictor = true;

	try {

	    c.voidEval("con<-dbConnect(drv, "+mySqlConfString+")");
	    c.voidEval("x<-dbGetQuery(con,  \""+predictorQuery+"\")");
	    singlePredictor = c.eval("ncol(x)==1").asBool().isTRUE();

	    if (!singlePredictor) {
		c.voidEval("discon<-dbDisconnect(con)");
		return lr;
	    }

	    c.voidEval("x<-x[,1]");
	    c.voidEval("rs<-dbGetQuery(con, \""+responseQuery+"\")");
	    c.voidEval("lr.m<-matrix(0, ncol(rs), 4)");

	    c.voidEval("for (j in 1:ncol(rs)) {"+
		       "   y<-rs[,j];"+
		       "   model<-lm(y~x);"+
		       "   info<-summary(model);"+
		       "   lr.m[j,1:2]<-model$coefficients;"+ //intercept - slope
		       "   lr.m[j,3]<-info$coefficients[2,4];"+ //p-value
		       "   lr.m[j,4]<-sqrt(info$r.squared);"+ //r
		       "}");

	    lr = c.eval("lr.m").asDoubleMatrix();

	    c.voidEval("discon<-dbDisconnect(con)");

	} catch(RSrvException rse) {
            System.out.println("Rserve exception: "+rse.getMessage());
	    cleanUp();
        } catch(Exception e) {
            System.out.println("Something went wrong, but it's not the Rserve: "+e.getMessage());
            e.printStackTrace();
	    cleanUp();
        }
	
	return lr;

    }

    private static void cleanUp() {

	try {
	    c.voidEval("try(discon<-dbDisconnect(con), silent = TRUE)");
	    c.voidEval("dbUnloadDriver(drv)");
	} catch(RSrvException rse) {
            System.out.println("Rserve exception: "+rse.getMessage());
        } finally {
	    if (c != null) {
		c.close();	    
	    }
	}
    }

    public void close() {

	try {
	    c.voidEval("try(discon<-dbDisconnect(con), silent = TRUE)");
	    c.voidEval("dbUnloadDriver(drv)");
	} catch(RSrvException rse) {
            System.out.println("Rserve exception: "+rse.getMessage());
        } finally {
	    if (c != null) {
		c.close();	    
	    }
	}
    }

    private static String readConnectionFile(String connFile) {

	FileReader theFile = null;
	BufferedReader fileIn = null;
	StringTokenizer str;
	String item, oneLine, dummy; 
	String host = null, port = null, user = null, password = null, dbname = null;
	
	// list the entries in the file and decompose them 
	try {
	    File inputFile = new File(connFile);
	    theFile = new FileReader(inputFile); // open the File
	    fileIn = new BufferedReader( theFile); // open BufferedReader 
	    while ((oneLine = fileIn.readLine() ) != null ) {
		// Write the line at hand to stdout, just for testing purposes 
		// System.out.println("["+oneLine+"]");
		// Construct a stringTokenizer for the line that we read with : delimited
		str = new StringTokenizer( oneLine, " ="); // true sets returnDelimiters flag 
		while ( str.hasMoreTokens()) {
		    item = str.nextToken();
		    // System.out.println("item:"+item);
		    if( item.equals("host")) {
			host=str.nextToken();
			// System.out.println("host:"+host);
			break; 
		    } // end if host 
		    if( item.equals("port")) {
			port=str.nextToken();
			// System.out.println("port:"+port);
			break; 
		    } // end if port
		    if( item.equals("user")) {
			user=str.nextToken();
			// System.out.println("user:"+user);
			break; 
		    } // end if password 
		    if( item.equals("password")) {
			password=str.nextToken();
			// System.out.println("password:"+password);
			break; 
		    } // end if password
		    if( item.equals("database")) {
			dbname=str.nextToken();
			// System.out.println("database:"+dbname);
			break; 
		    } // end if password 
		    
		} // next token in this line 
	    } // next line in the file 
	} // end try opening the file 
	catch ( Exception e ) { System.out.println( e); }  
	
	try { // closing the file
	    if( fileIn != null) fileIn.close();
	    if( theFile != null) theFile.close();
	} catch ( Exception e ) { System.out.println( e); }
	
	return ("username = \""+user+"\", password = \""+password+"\", dbname = \""+dbname+"\", host = \""+host+"\", port = "+port);

    } // end class readConnectionFile

} // end class mySQLconnect 

