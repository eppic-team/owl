package tools.trunk;

import tools.trunk.*;
import java.io.*;
import java.net.*;
 
public class Machine {

    public static String getClient() {
	
	String client = "";
 
	try {
	    client = InetAddress.getLocalHost().getHostName();
	} catch (Exception e) {
	    System.out.println(e);
	}
	
	return client;

    }

    public static String getOSArch() {

	return (System.getProperty("os.arch"));

    }
    
    public static String getOSName() {

	return (System.getProperty("os.name"));
	
    }

    public static String getOSver() {

	return (System.getProperty("os.version"));
	
    }

}
 
