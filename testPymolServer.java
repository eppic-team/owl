import org.apache.xmlrpc.*;
import java.net.*;
import java.util.*;
import java.io.*;

/**
 * Package:		tools
 * Class: 		testPymolServer
 * Author:		Henning Stehr, stehr@molgen.mpg.de
 * Date:		2/Feb/2006
 * 
 * Simple test class for XML-RPC calls to a remote PyMol server. Run PyMol
 * with the -R option to run the server.
 * 
 * Changelog:
 * 06/02/02 first created by HS
 */
public class testPymolServer {

	/**
	 * Send the first command line parameter to the pymol server.
	 */
	public static void main(String[] args) {
		
		String command = "load /project/StruPPi/PDBs/mainPDB/1RX4.pdb, hello";
		if(args.length > 0) {
			command = args[0];
		}
		
	    try {
	    String myURL = "http://gelb:9123";
	    XmlRpcClient client = new XmlRpcClient(myURL);
	    Vector<String> myvector = new Vector<String>();
	    myvector.add(command);
	    try {
	    	client.execute("do",myvector);	    
	    }
	    catch (IOException e) {
	    	e.printStackTrace();
	    }
	    catch( XmlRpcException e) {
	    	e.printStackTrace();
	    }
	    }
	    catch(MalformedURLException e){
	    	e.printStackTrace();
	    }
	    
	    
	}

}
