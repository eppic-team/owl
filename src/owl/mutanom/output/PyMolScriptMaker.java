package owl.mutanom.output;

import java.io.*;

import owl.core.util.Interval;

/**
 * Generates Pymol scripts for visualizing mutations and other features in a given protein.
 * @author stehr
 *
 */
public class PyMolScriptMaker {
	
	/*--------------------------- type definitions --------------------------*/
	public enum Color {RED, GREEN, BLUE, CYAN, MAGENTA, YELLOW, BLACK, WHITE, ORANGE, WHEAT, PALEGREEN, LIGHTBLUE, LIGHTPINK, PALEYELLOW, PALECYAN, LIGHTORANGE, BLUEWHITE};
	
	/*------------------------------ constants ------------------------------*/
	public static final String PYMOL_EXECUTABLE = "/project/StruPPi/bin/pymol";
	
	/*--------------------------- member variables --------------------------*/
	private File scriptFile;
	private PrintStream script;
	
	/*----------------------------- constructors ----------------------------*/
	public PyMolScriptMaker(boolean keepScript) throws IOException {
		scriptFile = File.createTempFile("pymol", ".pml");
		if(keepScript) {
			System.out.println("Writing script file to " + scriptFile);
		} else {
			scriptFile.deleteOnExit();
		}
		script = new PrintStream(scriptFile);
	}
	
	/*---------------------------- public methods ---------------------------*/
	/**
	 * Load a number of global visualization settings we commonly use.
	 */
	public void loadDefaultSettings() {
		script.println("bg white");
		script.println("set sphere_scale, 1.2"); // 5.0
		script.println("set sphere_transparency, 0.0"); // 0.7
		script.println("set transparency, 0.9");
		script.println("set surface_color, white");
		script.println("set side_chain_helper, 1");
	}	
	
	/**
	 * Load a file in PyMol. File can be a PDB or pymol session file.
	 */
	public void load(File file, String objName) {
		script.println("load " + file + ", " + objName);
		script.println("color gray");
		script.println("remove solvent");
		script.println("orient");
		script.println("as cartoon");		
		script.println("show surface");
	}
	
	/**
	 * Highlight a residue in the currently loaded object. Creates a selection with the given name,
	 * turns on stick representation for the selection, assigns the given color and draws a sphere around the Ca.
	 * @param objName the object in which a residue is to be highlighted
	 * @param pos the residue number in the given chain to be highlighted
	 * @param chain the chain in which a residue is to be highlighted
	 * @param color the color to assign to the residue
	 * @param selName the name of the selection to be created
	 * @param grpName TODO
	 */
	public void highlightResidue(String objName, int pos, char chain, Color color, String selName, String grpName) {
		//script.printf("pseudoatom %s, %s and chain %s and resi %d and name ca\n", selName, objName, chain, pos);
		script.printf("select %s, %s and chain %s and resi %d and name ca\n", selName, objName, chain, pos);
		script.println("color " + color + ", " + selName);
		//script.println("show sticks, " + selName);
		script.println("show spheres, " + selName);
		if(grpName != null) {
			script.println("group " + grpName + ", " + selName);
		}
	}
	
	/**
	 * Highlight an interval in the currently loaded object. Creates a selection with the given name,
	 * and assigns the given color to it. Multiple calls with the same selName should add to the selection.
	 * @param objName the object in which the interval is to be highlighted
	 * @param intv the interval in the given chain to be highlighted
	 * @param chain the chain in which a residue is to be highlighted (has to match the loaded pdb file)
	 * @param color the color to assign to the residue (Color is a local subclass)
	 * @param selName the name of the selection to be created
	 */
	public void highlightInterval(String objName, Interval intv, char chain, Color color, String selName) {
		//script.printf("pseudoatom %s, %s and chain %s and resi %d and name ca\n", selName, objName, chain, pos);
		script.printf("select %s, %s and chain %s and resi %d-%d\n", selName, objName, chain, intv.beg, intv.end);
		script.println("color " + color + ", " + selName);
		//script.println("show sticks, " + selName);
		//script.println("show spheres, " + selName);
	}
	
	/**
	 * Mutates the given residue to the given amino acid type.
	 * TODO: Not yet implemented.
	 * @param objName
	 * @param pos
	 * @param targetAminoAcid
	 */
	public void mutateResidue(String objName, int pos, String targetAminoAcid) {
		// TODO: Not yet implemented
	}
	
	/**
	 * Creates a duplicate of the given object with the givem new name.
	 * @param objName the name of the original object
	 * @param newObjName the name of the new (duplicate) object
	 */
	public void copyObject(String objName, String newObjName) {
		script.println("copy " + newObjName + ", " + objName);
	}
	
	/**
	 * Create an object from a selection
	 * @param objName the target object name
	 * @param selName the source selection name
	 */
	public void createObject(String objName, String selName) {
		script.println("create " + objName + ", " + selName);
		script.println("zoom");
	}
	
	/**
	 * Saves the current view as a raytraced png image.
	 * @param pngFile the image file to be written
	 */
	public void writePng(File pngFile) {
		script.println("draw");
		script.println("png " + pngFile);
	}
	
	/**
	 * Save the current PyMol session to a file.
	 * @param outFile the name of the session file to be written
	 */
	public void saveSession(File outFile) {
		script.println("save " + outFile);
	}
	
	/**
	 * Close output stream and run script
	 * @throws IOException 
	 */
	public void executeAndClose() throws IOException {
		script.println("quit");
		script.close();
		String cmd = PYMOL_EXECUTABLE + " " + scriptFile;
		System.out.println(cmd);
		Runtime.getRuntime().exec(cmd);
	}
	
}
