package sequence;

import java.io.File;
import java.io.IOException;

/**
 * A collection of little tools related to Blast and Processing Blast output.
 * @author stehr
 *
 */
public class BlastUtils {

	/**
	 * Uses a perl script to render a PNG image of the given blast output file.
	 * @param blastOutputFile
	 * @param imgFile
	 */
	public static void renderBlast(File blastOutputFile, File imgFile) throws IOException {
		String renderScript = "/project/StruPPi/CASP8/scripts/render_blast.pl";
		String cmdLine = String.format("%s %s > %s", renderScript, blastOutputFile.getAbsolutePath(), imgFile.getAbsolutePath());
		Runtime.getRuntime().exec(cmdLine);
	}
	
	/**
	 * Testing some of the method in this class.
	 * @param args
	 */
	public static void main(String[] args) {
		File blastOutput = new File("");
		File imgFile = new File("");
		try {
			renderBlast(blastOutput, imgFile);
		} catch(IOException e) {
			System.out.println("RenderBlast failed: " + e.getMessage());
		}
	}
	
}
