package tools;

import java.io.*;

/**
 * Package:		tools
 * Class: 		PyMol
 * Author:		Real Gagnon (http://www.rgagnon.com/javadetails/java-0014.html)
 * 				Copyright by Ioannis Filippis, filippis@molgen.mpg.de
 * Date:		03/03/2006
 * 
 * SystemCmd's static exec method will spawn an external process to execute an external 
 * program or to launch a unix script or just execute a simple linux command. 
 * The output is captured and returned as string.
 *  
 * Notes:
 * - There are two versions of the exec method, one taking a String array and one just a String
 * - If you want to launch a Unix script and especially use the < > | piping options,
 * 	 you must invoke the command processor
 * 		String[] cmd = {"/bin/sh", "-c", "ls > hello"};
 * 		SystemCmd.exec(cmd);
 * - If you need to pass arguments, it's safer to a String array especially if they contain spaces.
 * 		String[] cmd = { "myProgram.exe", "-o=This is an option" };
 * 		SystemCmd.exec(cmd);
 * 
 * Changelog:
 * 14/03/06 Added exec method to take only a string as argument by JD
 * 03/03/06 first created/copied by IF
 */

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
	
	public static String exec(String cmd) {
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

