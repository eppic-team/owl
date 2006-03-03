package tools;

import java.io.*;
 
public class SystemCmd {

    public static String exec(String[] cmd) {
	
	String output = "";

	try {

	    String line = "";
	    Process p = Runtime.getRuntime().exec(cmd);
	    BufferedReader input = new BufferedReader(new InputStreamReader(p.getInputStream()));
	    while ((line = input.readLine()) != null) {
		output = output + line + "\n";
	    }
	    input.close();
	
	} catch (Exception err) {
	    err.printStackTrace();
	}
	
	return output;

    }

}
 
