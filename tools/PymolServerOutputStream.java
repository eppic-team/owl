package tools;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Vector;
//import org.apache.xmlrpc.*;
import org.apache.xmlrpc.client.*;
//import org.apache.xmlrpc.common.*;
import org.apache.xmlrpc.XmlRpcException;
/**
 * Package:		tools
 * Class: 		PymolServerOutput
 * Author:		Henning Stehr, stehr@molgen.mpg.de
 * Date:		2/Feb/2006
 * 
 * This class provides an output stream to send commands to a
 * pymol server. It wraps the XML-RPC connection to the server
 * using the url taken by the constructor. For convenience the
 * OutputStream should be wrapped into a PrintStream.
 * 
 * Code example:
 * url = "http://localhost:9123";
 * PrintWriter serverOut = new PrintWriter(new PymolServerOutputStream(url), true);
 * serverOut.println("load /project/StruPPi/PDBs/mainPDB/1RX4.pdb, hello");
 * 
 * Notes: 
 * - When constructing the PrintWriter, make sure to turn on automatic flushing
 *   by giving "true" as the second parameter, otherwise flush() has to be invoked
 *   manually to send the command to the server.
 * - Automatic flushing is only invoked when using 'println' rather than 'print'.
 * - For some reason, this does not work with a PrintStream only with a PrintWriter.
 * - '\n' in command strings causes problems
 * 
 * Changelog:
 * 06/02/02 first created by HS
 */
public class PymolServerOutputStream extends OutputStream {

	static final String 	DEFAULTXMLRPCCOMMAND	= "do";
	static final int		INITIALBUFFERCAPACITY   = 255;
	static final int		BUFFERCAPACITYINCREASE  = 255;
	public static final int        PYMOLCOMMANDLENGTHLIMIT = 1060;
	
	byte[]			commandBuffer;
	int				commandBufferPtr,  // current load
					commandBufferCap;  // capacity
	XmlRpcClient    client;
	
	/**
	 * Create a new output stream to the server running at the given url
	 */ 
	public PymolServerOutputStream(String url) {
		
		// create empty command buffer
		commandBuffer = new byte[INITIALBUFFERCAPACITY];
		commandBufferCap = INITIALBUFFERCAPACITY;
		commandBufferPtr = 0;
		
		try {
		// initialize connection to server
	    this.client = new XmlRpcClient();
	    XmlRpcClientConfigImpl config = new XmlRpcClientConfigImpl();
	    config.setServerURL(new URL(url));
	    client.setConfig(config);
		} catch(MalformedURLException e) {
			System.err.println("Error: Malformed URL in constructor of PymolServerOutputStream");
			e.printStackTrace();
		}
	}


	/**
	 * Store bytes in command buffer
	 */ 
	@Override
	public void write(int arg0) throws IOException {
		if(commandBufferPtr >= commandBufferCap) {
			// get more space
			byte[] newBuffer = new byte[commandBufferCap + BUFFERCAPACITYINCREASE];
			for(int i = 0; i < commandBufferPtr; i++) {
				newBuffer[i] = commandBuffer[i];
			}
			commandBuffer = newBuffer;
			commandBufferCap += BUFFERCAPACITYINCREASE;
		}
		commandBuffer[commandBufferPtr++] = (byte) arg0;
	}
	

	/**
	 * Send command to server whenever flush is invoked
	 */ 
	@Override
	public void flush() throws IOException {
		// prevent overlong commands from being sent to server
		if(commandBufferPtr > PYMOLCOMMANDLENGTHLIMIT + 1) {
			throw new IOException();
		}
		// prepare parameter vector
		Vector<String> myvector = new Vector<String>();
	    String commandString = new String(this.commandBuffer, 0, commandBufferPtr).trim();
	    myvector.add(commandString);
	    // reset command buffer
	    this.commandBuffer = new byte[INITIALBUFFERCAPACITY];
	    this.commandBufferCap = INITIALBUFFERCAPACITY;
	    this.commandBufferPtr = 0;
	    // submit command
	    try {
	    	this.client.execute(DEFAULTXMLRPCCOMMAND, myvector);
	    	
	    } catch(XmlRpcException e) {
	    	// System.out.println("XMP-RPC exception occured.");
	    	// send IOException instead, so use of XML-RPC is transparent
	    	throw new IOException(e.getMessage()); 
	    }
	}
	
	// test function for current class
	public static void main(String[] args) {

		String	serverUrl = 	"http://gelb:9123";
		String  command0 =      "delete hello";
		String	command1 =		"load /project/StruPPi/PDBs/mainPDB/1RX4.pdb, hello";
		String	command2 =		"show cartoon";
		String	command3 =		"select hello2, resi 12";
		String	command4 =		"show sticks, hello2";
		String	command5 =		"zoom hello2";
		
		// try to connect directly
    	System.out.println("Trying to connect directly by XML-RPC...");		
	    try {
		    
		    XmlRpcClient client = new XmlRpcClient();
		    XmlRpcClientConfigImpl config = new XmlRpcClientConfigImpl();
		    config.setServerURL(new URL(serverUrl));
		    client.setConfig(config);

		    
		    Vector<String> myvector = new Vector<String>();
		    myvector.add(command1);
		    try {
		    	client.execute(DEFAULTXMLRPCCOMMAND,myvector);
		    	System.out.println("done.");
		    }
		    catch( XmlRpcException e) {
		    	e.printStackTrace();
		    }
	    }
	    catch(MalformedURLException e){
	    	e.printStackTrace();
	    }		
		
		// try using OutputStream
    	System.out.println("Trying to connect via PymolServerOutputStream...");	
		try {
			// create new output stream
			OutputStream serverOut = new PymolServerOutputStream(serverUrl);
			
			// send data to output stream
			serverOut.write(command0.getBytes());
			// execute command
			serverOut.flush();
			
	    	System.out.println("done.");
			
		} catch(MalformedURLException e) {
			System.out.println("Error: The server URL is wrong.");
			e.printStackTrace();
			System.exit(1);
		} catch(IOException e) {
			System.out.println("Error writing to server url.");
			e.printStackTrace();
			System.exit(1);			
		}
			
		// now try using OutputStream wrapped into PrintStream
    	System.out.println("Trying to connect via PrintWriter...");		
		PrintWriter serverOutPw = new PrintWriter(new PymolServerOutputStream(serverUrl), true);
		serverOutPw.println(command0);
		serverOutPw.println(command1);	
		serverOutPw.println(command2);
		serverOutPw.println(command3);
		serverOutPw.println(command4);	
		serverOutPw.println(command5);
		System.out.println("done.");
	}
	
}
