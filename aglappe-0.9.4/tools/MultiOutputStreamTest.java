package tools;

import junit.framework.TestCase;

import java.io.OutputStream;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

public class MultiOutputStreamTest extends TestCase {

	/*
	 * Test method for 'tools.MultiOutputStream'
	 */
	public void testMultiOutputStream() throws FileNotFoundException {
		
		 final String url1 = "http://gelb:9123";
		 final String url2 = "http://gelb:9124";		 
		 final String url3 = "http://gelb:9125";
		 final String file1 = "/home/stehr/test.log";
		 final String command1 = "load /project/StruPPi/PDBs/mainPDB/1RX4.pdb, hello";
		 final String command2 = "rock";
		 
		 OutputStream os1 = new PymolServerOutputStream(url1);
		 OutputStream os2 = new PymolServerOutputStream(url2);
		 OutputStream os3 = new PymolServerOutputStream(url3);
		 OutputStream os4 = new FileOutputStream(file1);
		 OutputStream multi = new MultiOutputStream(os1, os2, os3, os4);
		 PrintWriter multiOut = new PrintWriter(multi, true);
		 multiOut.println(command1);
		 multiOut.println(command2);		 
	}
}
