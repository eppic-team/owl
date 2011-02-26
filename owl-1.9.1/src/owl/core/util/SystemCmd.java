package owl.core.util;

import java.io.*;

/**
 * Package:		tools
 * Class: 		SystemCmd
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
 * 19/09/06 Added functionality to run commands in threads by JD
 * 14/03/06 Added exec method to take only a string as argument by JD
 * 03/03/06 first created/copied by IF
 */

public class SystemCmd extends Thread {
	
	String command;
	int exitVal;
	
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

	/**
	 * Constructor. To be used if we want to run a command as a new thread.
	 * When running the command in threads the only thing we get back from them is the exit value (see run method)
	 * @param cmd
	 */
	public SystemCmd(String cmd) {
		this.command=cmd;
		this.exitVal=-1;
	}

	@Override
	public void run(){
		try {
			Process p = Runtime.getRuntime().exec(command);
			p.waitFor();
			exitVal = p.exitValue();
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Couldn't execute command: "+command);
		} catch (InterruptedException e) {			
			e.printStackTrace();
		}		
	}
	
	public int getExitVal(){
		return exitVal;
	}
}

