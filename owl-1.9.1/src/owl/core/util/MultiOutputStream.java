package owl.core.util;

import java.io.IOException;
import java.io.OutputStream;

/**
 * Package:		tools
 * Class: 		MultiOutputStream
 * Author:		Henning Stehr, stehr@molgen.mpg.de
 * Date:		9/Mar/2006
 * 
 * An OutputStream which allows to write to multiple output streams.
 * Invoking a method of MultiOutputStream invokes the same method
 * with the same parameters in each of the specified output streams.
 * 
 * Code example:
 * url1 = "http://gelb:9123";
 * url2 = "http://blau:9123";
 * file1 = "test.log";
 * OutputStream os1 = new PymolServerOutputStream(url1);
 * OutputStream os2 = new PymolServerOutputStream(url2);
 * OutputStream os3 = new FileOutputStream(file1);
 * OutputStream multi = new MultiOutputStream(os1, os2, os3);
 * PrintWriter multiOut = new PrintWriter(multi, true);
 * multiOut.println("load /project/StruPPi/PDBs/mainPDB/1RX4.pdb, hello");
 * 
 * Changelog:
 * 2006/03/09 first created by HS
 */
public class MultiOutputStream extends OutputStream {

	OutputStream[] streams;
	
	/**
	 * Creates a new MultiOutputStream.
	 * @param outputStreams The output streams which should be contained
	 * in this MultiOutputStream.
	 */
	public MultiOutputStream(OutputStream...outputStreams) {
			this.streams = outputStreams;
	}

	/**
	 * Writes a byte to all output streams contained in this
	 * MultiOutputStream.
	 * @param The byte to be written to the output streams
	 */
	@Override
	public void write(int arg0) throws IOException {
		for(OutputStream stream : streams) {
			stream.write(arg0);
		}
	}
	
	/**
	 * Flushes all output streams contained in this
	 * MultiOutputStream.
	 */ 
	@Override
	public void flush() throws IOException {
		for(OutputStream stream : streams) {
			stream.flush();
		}	
	}
	
	/**
	 * Closes all output streams contained in this
	 * MultiOutputStream.
	 */ 
	@Override
	public void close() throws IOException {
		for(OutputStream stream : streams) {
			stream.close();
		}
	}	

}
