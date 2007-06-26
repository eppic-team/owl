import org.apache.xmlrpc.*;
import org.apache.xmlrpc.client.XmlRpcClient;
import org.apache.xmlrpc.client.XmlRpcClientConfigImpl;

import java.net.*;
import java.util.*;

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
		
		String command = "load /project/StruPPi/Databases/pdb/rx/1rx4.pdb, hello";
		if(args.length > 0) {
			command = args[0];
		}
		
	    try {
	    String myURL = "http://anthrazit:9123";
	    XmlRpcClient client = new XmlRpcClient();
	    XmlRpcClientConfigImpl config = new XmlRpcClientConfigImpl();
	    config.setServerURL(new URL(myURL));
	    client.setConfig(config);

	    Vector<String> myvector = new Vector<String>();
	    myvector.add(command);
	    try {
	    	client.execute("do",myvector);	    
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
