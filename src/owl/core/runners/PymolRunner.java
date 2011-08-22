package owl.core.runners;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

import owl.core.structure.ChainInterface;
import owl.core.structure.Residue;

public class PymolRunner {
	
	/**
	 * We use 26 colors corresponding to chain letters A to Z (second 13 are repeated from first 13)
	 */
	private static final String[] DEF_CHAIN_COLORS = 
	{"green","cyan","yellow","white","lightblue","magenta","red","orange","wheat","limon","salmon","palegreen","lightorange",
	 "green","cyan","yellow","white","lightblue","magenta","red","orange","wheat","limon","salmon","palegreen","lightorange",};
	
	private static final String DEF_SYM_RELATED_CHAIN_COLOR = "grey";
	
	private static final String DEF_TN_STYLE = "cartoon";
	private static final String DEF_TN_BG_COLOR = "white";
	private static final int[] DEF_TN_HEIGHTS = {75};
	private static final int[] DEF_TN_WIDTHS = {75};
	
	private File pymolExec;
	private String[] chainColors;
	private String symRelatedColor;
	
	public PymolRunner(File pymolExec) {
		this.pymolExec = pymolExec;
		chainColors = DEF_CHAIN_COLORS;
		symRelatedColor = DEF_SYM_RELATED_CHAIN_COLOR;
	}
	
	public void setColors(String[] chainColors, String symRelatedColor) {
		this.chainColors = chainColors;
		this.symRelatedColor = symRelatedColor;
	}
	
	
	/**
	 * Generate thumbnail files for this interface with PyMol for given pdbFile. 
	 * Output files will be written to same dir as pdbFile using given base name.
	 * @param pymolExe
	 * @param pdbFile
	 * @param base
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void generateThumbnails(ChainInterface interf, File pdbFile, String base) throws IOException, InterruptedException {
		
		File[] pngFiles = new File[DEF_TN_HEIGHTS.length];
		for (int i=0;i<DEF_TN_HEIGHTS.length;i++) {
			pngFiles[i] = new File(pdbFile.getParent(),base+"."+DEF_TN_WIDTHS[i]+"x"+DEF_TN_HEIGHTS[i]+".png");
		}
		String[] chains = new String[2];
		chains[0] = interf.getFirstMolecule().getPdbChainCode();
		chains[1] = interf.getSecondPdbChainCodeForOutput();
		generatePng(pdbFile, chains, interf.isSymRelated(), pngFiles, DEF_TN_STYLE, DEF_TN_BG_COLOR, DEF_TN_HEIGHTS, DEF_TN_WIDTHS);
	}
	
	/**
	 * Generates png images of the desired heights and widths with the specified style and 
	 * coloring each chain with a color of {@link #CHAIN_COLORS}
	 * @param pdbFile
	 * @param chains the chains present in the file
	 * @param isSymRelated whether the PDB file contains crystal-symmetry-related chains (originally they had the same chain code) 
	 * @param outPngFiles output png file names
	 * @param style can be cartoon, surface, spheres
	 * @param bgColor the background color for the image: black, white, gray
	 * @param heights
	 * @param widths 
	 * @throws IOException 
	 * @throws InterruptedException 
	 * @throws IllegalArgumentException if heights length differs from widhts length
	 */
	public void generatePng(File pdbFile, String[] chains, boolean isSymRelated, File[] outPngFiles, String style, String bgColor, int[] heights, int[] widths) 
	throws IOException, InterruptedException {
		
		if (heights.length!=widths.length || heights.length!=outPngFiles.length) 
			throw new IllegalArgumentException("The number of heights is different from the number of widths or the number of output png files");
		String molecName = pdbFile.getName().substring(0, pdbFile.getName().lastIndexOf('.'));
		
		List<String> command = new ArrayList<String>();
		command.add(pymolExec.getAbsolutePath());
		command.add("-q");
		command.add("-c");

		StringBuffer pymolScriptBuilder = new StringBuffer();
		
		pymolScriptBuilder.append("load "+pdbFile.getAbsolutePath() + ";");
		
		pymolScriptBuilder.append("bg "+bgColor + ";");
		
		pymolScriptBuilder.append("orient;");
		
		pymolScriptBuilder.append("remove solvent;");
		
		pymolScriptBuilder.append("as "+style + ";");
		
		
		for (int c=0;c<chains.length;c++) {
			char letter = chains[c].charAt(0);
			String color = getChainColor(letter, c, isSymRelated);
			pymolScriptBuilder.append("color "+color+", "+molecName+" and chain "+letter + ";");
		}

		for (int i=0;i<heights.length;i++) {
			pymolScriptBuilder.append("viewport "+heights[i]+","+widths[i] + ";");
			
			pymolScriptBuilder.append("ray;");
			
			pymolScriptBuilder.append("png "+outPngFiles[i].getAbsolutePath() + ";");
		}
		
		pymolScriptBuilder.append("quit;");
		
		command.add("-d");
		command.add(pymolScriptBuilder.toString());
		
		Process pymolProcess = new ProcessBuilder(command).start();
		int exit = pymolProcess.waitFor();
		if (exit!=0) {
			throw new IOException("Pymol exited with error status "+exit);
		}
	}
	
	/**
	 * Generates pymol pse file and pml script for given interface producing a 
	 * mixed cartoon/surface representation of interface with selections 
	 * coloring each chain with a color of {@link #CHAIN_COLORS}
	 * @param pymolExec
	 * @param interf
	 * @param pdbFile
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public void generateInterfPse(ChainInterface interf, File pdbFile, File pseFile, File pmlFile) 
	throws IOException, InterruptedException {
		
		String molecName = pdbFile.getName().substring(0, pdbFile.getName().lastIndexOf('.'));

		List<String> command = new ArrayList<String>();
		command.add(pymolExec.getAbsolutePath());
		command.add("-q");
		command.add("-c");


		// NOTE we used to pass all commands in one string after -d (with the pymolScriptBuilder StringBuffer.
		//      But pymol 1.3 and 1.4 seem to have problem with very long strings (causing segfaults)
		//      Because of that now we write most commands to pml file (which we were doing anyway so that users can 
		//      use the pml scripts if they want) and then load the pmls with pymol "@" command
		
		
		StringBuffer pymolScriptBuilder = new StringBuffer();
		PrintStream pml = new PrintStream(pmlFile);
		
		pymolScriptBuilder.append("load "+pdbFile.getAbsolutePath()+";");

		String cmd;

		cmd = "orient";
		writeCommand(cmd, pml);
		
		cmd = "remove solvent";
		writeCommand(cmd, pml);
		
		cmd = "as cartoon";
		writeCommand(cmd, pml);
		
		char chain1 = interf.getFirstMolecule().getPdbChainCode().charAt(0);
		char chain2 = interf.getSecondPdbChainCodeForOutput().charAt(0);
		
		String color1 = getChainColor(chain1, 0, interf.isSymRelated());
		String color2 = getChainColor(chain2, 1, interf.isSymRelated());
		
		cmd = "color "+color1+", "+molecName+" and chain "+chain1;
		writeCommand(cmd, pml);
		cmd = "color "+color2+", "+molecName+" and chain "+chain2;
		writeCommand(cmd, pml);

		cmd = getSelString("core", chain1, interf.getFirstRimCore().getCoreResidues());
		writeCommand(cmd, pml);
		cmd = getSelString("core", chain2, interf.getSecondRimCore().getCoreResidues());
		writeCommand(cmd, pml);
		cmd = getSelString("rim", chain1, interf.getFirstRimCore().getRimResidues());
		writeCommand(cmd, pml);
		cmd = getSelString("rim", chain2, interf.getSecondRimCore().getRimResidues());
		writeCommand(cmd, pml);
		
		cmd = "select interface"+chain1+", core"+chain1+" or rim"+chain1;
		writeCommand(cmd, pml);
		cmd = "select interface"+chain2+", core"+chain2+" or rim"+chain2;
		writeCommand(cmd, pml);
		cmd = "select bothinterf , interface"+chain1+" or interface"+chain2;
		writeCommand(cmd, pml);
		cmd = "show surface, chain "+chain1;
		writeCommand(cmd, pml);
		cmd = "show surface, chain "+chain2;
		writeCommand(cmd, pml);
		//pymolScriptBuilder.append("color blue, core"+chains[0]+";");
		//pymolScriptBuilder.append("color red, rim"+chains[0]+";");
		cmd = "color red, interface"+chain1;
		writeCommand(cmd, pml);
		//pymolScriptBuilder.append("color slate, core"+chains[1]+";");
		//pymolScriptBuilder.append("color raspberry, rim"+chains[1]+";");
		cmd = "color raspberry, interface"+chain2;
		writeCommand(cmd, pml);
		cmd = "show sticks, bothinterf";
		writeCommand(cmd, pml);
		cmd = "set transparency, 0.35";
		writeCommand(cmd, pml);
		//pymolScriptBuilder.append("zoom bothinterf"+";");
		cmd = "select resi 0";// so that the last selection is deactivated
		writeCommand(cmd, pml);
		
		pml.close();
		
		pymolScriptBuilder.append("@ "+pmlFile+";");
		
		pymolScriptBuilder.append("save "+pseFile+";");
		
		pymolScriptBuilder.append("quit;");
		
		command.add("-d");
		
		//System.out.println(pymolScriptBuilder.toString());
		
		command.add(pymolScriptBuilder.toString());

		
		Process pymolProcess = new ProcessBuilder(command).start();
		int exit = pymolProcess.waitFor();
		if (exit!=0) {
			throw new IOException("Pymol exited with error status "+exit);
		}
	}
	
	private String getChainColor(char letter, int index, boolean isSymRelated) {
		String color = null;
		if (isSymRelated && index!=0) {
			color = symRelatedColor;
		} else {
			if (letter<'A' || letter>'Z') {
				// if out of the range A-Z then we assign simply a color based on the chain index
				color = chainColors[index%chainColors.length];
			} else {
				// A-Z correspond to ASCII codes 65 to 90. The letter ascii code modulo 65 gives an indexing of 0 (A) to 25 (Z)
				// a given letter will always get the same color assigned
				color = chainColors[letter%65];	
			}
		}
		return color;
	}
	
	private String getResiSelString(List<Residue> list) {
		if (list.isEmpty()) return "0";
		StringBuffer sb = new StringBuffer();
		for (int i=0;i<list.size();i++) {
			sb.append(list.get(i).getSerial());
			if (i!=list.size()-1) sb.append("+");
		}
		return sb.toString();
	}

	private String getSelString(String namePrefix, char chainName, List<Residue> list) {
		return "select "+namePrefix+chainName+", chain "+chainName+" and resi "+getResiSelString(list);
	}
	
	private void writeCommand(String cmd, PrintStream ps) {
		//if (sb!=null) {
		//	sb.append(cmd+";");
		//}
		if (ps!=null) {
			ps.println(cmd);
		}
		
	}
	
	/**
	 * Reads from properties file the chain colors: 26 colors, one per alphabet letter
	 * and a color for the sym related chain
	 * @param is
	 * @throws IOException
	 */
	public void readColorsFromPropertiesFile(InputStream is) throws IOException {
		
		Properties p = new Properties();
		p.load(is);

		chainColors = new String[26];
		char letter = 'A';
		for (int i=0;i<26;i++) {
			chainColors[i] = p.getProperty(Character.toString(letter)); 
			letter++;
		}
		symRelatedColor = p.getProperty("SYMCHAIN");
	}	

}