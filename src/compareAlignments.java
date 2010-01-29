import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import proteinstructure.Alignment;
import proteinstructure.AlignmentConstructionError;
import proteinstructure.FileFormatError;
import proteinstructure.PairwiseSequenceAlignment;
import proteinstructure.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;


public class compareAlignments {

	public static int doCompare(File al1File, File al2File, PrintStream out) {
		Alignment al1 = null;
		Alignment al2 = null;
		try {
			al1 = new Alignment(al1File.getAbsolutePath(), Alignment.FASTAFORMAT);
			al2 = new Alignment(al2File.getAbsolutePath(), Alignment.FASTAFORMAT);
			
			for(String tag:al1.getTags()) {
				//out.println(tag);
				String s1 = al1.getSequenceNoGaps(tag);
				String s2 = al2.getSequenceNoGaps(tag);
				if(!s1.equals(s2)) {
					out.println(">" + tag);
					PairwiseSequenceAlignment pw = new PairwiseSequenceAlignment(s1,s2,al1File.getName(), al2File.getName());
					pw.printAlignment();
				} else {
					//out.println(">" + tag);
				}
			}	
			
			return 0; // all right
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (FileFormatError e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (AlignmentConstructionError e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (PairwiseSequenceAlignmentException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return 1; // error
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if(args.length < 2) {
			System.out.println("Usage: java compareAlignments <alignment1> <alignment2>");
			System.exit(1);
		}
		
		String al1FileName = args[0];
		String al2FileName = args[1];
		PrintStream out = System.out;
		
		File al1File = new File(al1FileName);
		File al2File = new File(al2FileName);
		int ret = doCompare(al1File, al2File, out);
		System.exit(ret);

	}

}
