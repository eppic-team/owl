package tests.tools;

import java.io.PrintWriter;
import java.io.OutputStream;

import tools.PymolServerOutputStream;

import junit.framework.TestCase;

public class PymolServerOutputStreamTest extends TestCase {

	/*
	 * Test method for 'tools.PymolServerOutputStream'
	 * Note that PyMol doesn't eat commands longer than 1060
	 * characters. It even crashes with a segfault.
	 */
	public void testPymolServerOutputStream() {
		
			final String	serverUrl = 	"http://gelb:9123";
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
