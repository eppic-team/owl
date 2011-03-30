package owl.core.runners;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import owl.core.structure.PdbLoadException;
import owl.core.structure.PdbfileParser;
import owl.core.util.StreamGobbler;

public class PymolRunner {
	
	/**
	 * We use 26 colors corresponding to chain letters A to Z (second 13 are repeated from first 13)
	 */
	private static final String[] CHAIN_COLORS = 
	{"green","cyan","yellow","white","lightblue","magenta","red","orange","wheat","limon","salmon","palegreen","lightorange",
	 "green","cyan","yellow","white","lightblue","magenta","red","orange","wheat","limon","salmon","palegreen","lightorange",};
	
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
	public void generatePng(File pdbFile, File[] outPngFiles, String style, String bgColor, int[] heights, int[] widths) 
	throws IOException, InterruptedException, PdbLoadException {
		
		if (heights.length!=widths.length || heights.length!=outPngFiles.length) 
			throw new IllegalArgumentException("The number of heights is different from the number of widths or the number of output png files");
		String molecName = pdbFile.getName().substring(0, pdbFile.getName().lastIndexOf('.'));
		PdbfileParser parser = new PdbfileParser(pdbFile.getAbsolutePath());
		String[] chains = parser.getChains();

		Process pymolProcess = Runtime.getRuntime().exec(pymolExec+" -q -c -p");
		
		new StreamGobbler("pymol_stdout", pymolProcess.getInputStream()).start();
		new StreamGobbler("pymol_stderr", pymolProcess.getErrorStream()).start();

		PrintWriter pymolIn = new PrintWriter(new BufferedOutputStream(pymolProcess.getOutputStream()));
				
		pymolIn.println("load "+pdbFile.getAbsolutePath());
		pymolIn.println("bg "+bgColor);
		pymolIn.println("remove solvent");
		pymolIn.println("as "+style);
		for (int c=0;c<chains.length;c++) {
			char letter = chains[c].charAt(0);
			String color = null;
			if (letter<'A' || letter>'Z') {
				// if out of the range A-Z then we assign simply a color based on the chain index
				color = CHAIN_COLORS[c%CHAIN_COLORS.length];
			} else {
				// A-Z correspond to ASCII codes 65 to 90. The letter ascii code modulo 65 gives an indexing of 0 (A) to 25 (Z)
				// a given letter will always get the same color assigned
				color = CHAIN_COLORS[letter%65];	
			}
			pymolIn.println("color "+color+", "+molecName+" and chain "+letter);
		}

		for (int i=0;i<heights.length;i++) {
			pymolIn.println("viewport "+heights[i]+","+widths[i]);
			pymolIn.println("ray");
			pymolIn.println("png "+outPngFiles[i].getAbsolutePath());
		}
		
		pymolIn.println("quit");
		pymolIn.flush();
		int exit = pymolProcess.waitFor();
		if (exit!=0) {
			throw new IOException("Pymol exited with error status "+exit);
		}
	}
	

}
