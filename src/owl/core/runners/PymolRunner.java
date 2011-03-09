package owl.core.runners;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;

import owl.core.structure.Pdb;
import owl.core.structure.PdbLoadException;
import owl.core.structure.PdbfilePdb;
import owl.core.util.StreamGobbler;

public class PymolRunner {
	
	private static final String[] CHAIN_COLORS = {"green","cyan","yellow","white","lightblue","magenta"};
	
	private File pymolExec;
	
	public PymolRunner(File pymolExec) {
		this.pymolExec = pymolExec;
	}

	/**
	 * Generates png images of the desired heights and widths with the specified style and 
	 * coloring each chain with a color of {@link #CHAIN_COLORS}
	 * @param pdbFile
	 * @param outPngFiles output png file names
	 * @param style can be cartoon, surface, spheres
	 * @param bgColor the background color for the image: black, white, gray
	 * @param heights
	 * @param widths 
	 * @throws IOException 
	 * @throws InterruptedException 
	 * @throws PdbLoadException 
	 * @throws IllegalArgumentException if heights length differs from widhts length
	 */
	public void generatePng(File pdbFile, File[] outPngFiles, String style, String bgColor, int[] heights, int[] widths) throws IOException, InterruptedException, PdbLoadException {
		if (heights.length!=widths.length || heights.length!=outPngFiles.length) 
			throw new IllegalArgumentException("The number of heights is different from the number of widths or the number of output png files");
		String molecName = pdbFile.getName().substring(0, pdbFile.getName().lastIndexOf('.'));
		Pdb pdb = new PdbfilePdb(pdbFile.getAbsolutePath());
		String[] chains = pdb.getChains();

		Process pymolProcess = Runtime.getRuntime().exec(pymolExec+" -q -c -p");
		
		new StreamGobbler("pymol_stdout", pymolProcess.getInputStream()).start();
		new StreamGobbler("pymol_stderr", pymolProcess.getErrorStream()).start();

		PrintWriter pymolIn = new PrintWriter(new BufferedOutputStream(pymolProcess.getOutputStream()));
		
		
		InputStreamReader pymolOut = new InputStreamReader(new BufferedInputStream(pymolProcess.getInputStream()));
		while (pymolOut.ready()) { 			
			System.out.write(pymolOut.read());
		}
		
		pymolIn.println("load "+pdbFile.getAbsolutePath());
		pymolIn.println("bg "+bgColor);
		pymolIn.println("remove solvent");
		pymolIn.println("as "+style);
		for (int c=0;c<chains.length;c++) {
			pymolIn.println("color "+CHAIN_COLORS[c%CHAIN_COLORS.length]+", "+molecName+" and chain "+chains[c]);
		}

		for (int i=0;i<heights.length;i++) {
			pymolIn.println("viewport "+heights[i]+","+widths[i]);
			pymolIn.println("ray");
			pymolIn.println("png "+outPngFiles[i].getAbsolutePath());
		}
		
		pymolIn.println("quit");
		pymolIn.flush();
		pymolProcess.waitFor();
	}
	

}
