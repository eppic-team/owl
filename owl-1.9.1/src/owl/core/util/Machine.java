package owl.core.util;

import java.net.*;
import java.util.Properties;
 
public class Machine {

    public static String getClient() {
		String client = "";	 
		try {
		    client = InetAddress.getLocalHost().getHostName();
		} catch (Exception e) {
		    System.err.println(e);
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

    public static String getUserName() {
    	return (System.getProperty("user.name"));
    }
    
    public static String getTempDir() {
    	return(System.getProperty("java.io.tmpdir"));
    }
    
    public static Properties getAllProperties() {
    	return System.getProperties();
    }
    
    public static long getTotalMemory() {
    	return Runtime.getRuntime().totalMemory(); 
    }
    
    public static long getFreeMemory() {
    	return Runtime.getRuntime().freeMemory();
    }
    
    public static void main(String[] args) {
    	System.out.println("Hi, I am " + getClient() + " and these are my properties:");
    	getAllProperties().list(System.out);
    	System.out.printf("Total memory: %,d bytes\n", getTotalMemory());
    	System.out.printf("Free memory : %,d bytes\n", getFreeMemory());
    }
}
 
