package owl.tests.core.util;

import java.io.File;
import java.io.OutputStream;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.net.InetAddress;
import java.net.UnknownHostException;

import org.junit.Test;

import owl.core.util.MultiOutputStream;
import owl.core.util.PymolServerOutputStream;


public class MultiOutputStreamTest {

	/*
	 * Test method for 'tools.MultiOutputStream'
	 */
	@Test
	public void testMultiOutputStream() throws FileNotFoundException, UnknownHostException {
		
		final String tmpDir = System.getProperty("java.io.tmpdir"); 
		final String host = InetAddress.getLocalHost().getHostName();
		
		final String url1 = "http://"+host+":9123";
		final String url2 = "http://"+host+":9124";		 
		final String url3 = "http://"+host+":9125";
		final File file1 = new File(tmpDir,"test.log");
		file1.deleteOnExit();
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
		multiOut.close();
	}
}
