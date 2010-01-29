import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.OutputStream;
import java.io.FileOutputStream;

import tools.PyMol;
import tools.PymolServerOutputStream;
import tools.MultiOutputStream;

public class testPyMol {

	/**
	 * Test class for the PyMol class (our java to PyMol API) 
	 * @param file where PyMol script will be written
	 * @author duarte
	 */
	public static void main(String[] args) {
		String file = "";
		boolean server=false;
		if (args.length<1){
			System.err.println("Give at least one file argument");
			System.exit(1);
		}		
		file = args[0];
		if (args.length>1){ // two arguments: output both file and server
			server=true;
		}
		
		
		String url = "http://anthrazit:9123";
		
		String pdbFileName = "/project/StruPPi/jose/tinker/benchmarking/1bxy_A.pdb";

		PrintWriter Out = null;		
		
		// to output only to server we would need the following PrintWriter
		//Out = new PrintWriter(new PymolServerOutputStream(url),true);
		
		if (server){
			OutputStream os1 = new PymolServerOutputStream(url);
			OutputStream os2 = null;
			try {
				os2 = new FileOutputStream(file);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
				System.exit(2);
			}
			OutputStream multi = new MultiOutputStream(os1, os2);
			Out = new PrintWriter(multi, true);
		} else {
			try {
				Out = new PrintWriter(new FileOutputStream(file),true);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		PyMol mypymol = new PyMol(Out);
		
		mypymol.loadPDB(pdbFileName);
		// more pymol commands

	}

}
