package tests.tools;

import java.io.PrintWriter;
import java.io.OutputStream;
import java.net.InetAddress;
import java.net.UnknownHostException;

//import org.junit.Test;

import tools.PymolServerOutputStream;

import static org.junit.Assert.*;

public class PymolServerOutputStreamTest {

	/**
	 * Test method for 'tools.PymolServerOutputStream'
	 * Note that PyMol doesn't eat commands longer than 1060
	 * characters. It even crashes with a segfault.
	 * @throws UnknownHostException 
	 */
	//@Test // commented out because this test needs the pymol server to be running and that's not implemented yet
	public void testPymolServerOutputStream() throws UnknownHostException {
		
		final String host = InetAddress.getLocalHost().getHostName();

		final String	serverUrl = 	"http://"+host+":9123";
		final int		from	   =	   10,
		to		   =	 1060,
		inc		   = 	   10;
		final char		dummyChar  =      'a';

		// open connection
		OutputStream serverOutOs = new PymolServerOutputStream(serverUrl);
		assertNotNull(serverOutOs);

		PrintWriter serverOutPw = new PrintWriter(serverOutOs, true);
		assertNotNull(serverOutPw);

		// send commands of different lenghths
		for(int l = from; l <= to; l += inc) {
			// make dummy command
			StringBuffer cmd = new StringBuffer(l);
			for(int i = 0; i < l; i++) {
				cmd.append(dummyChar);
			}
			System.out.println("Sending string of length " + l + " to PyMol server");
			serverOutPw.println(cmd.toString());
			assertFalse(serverOutPw.checkError());
		}
	}
}
