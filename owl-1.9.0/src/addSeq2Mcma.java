import java.io.*;

import owl.core.sequence.alignment.AlignmentConstructionError;
import owl.core.sequence.alignment.MultipleSequenceAlignment;
import owl.core.util.FileFormatError;


/**
 * A script to add sequences to a raw MCMA alignment (which contains only the aligned positions
 * @author stehr
 */
public class addSeq2Mcma {

	public static int addSequence(File alFile, File seqFile, PrintStream out) {
		MultipleSequenceAlignment al = null;
		MultipleSequenceAlignment seqs = null;
		try {
			al = new MultipleSequenceAlignment(alFile.getAbsolutePath(), MultipleSequenceAlignment.FASTAFORMAT);
			seqs = new MultipleSequenceAlignment(seqFile.getAbsolutePath(), MultipleSequenceAlignment.FASTAFORMAT);
			
			for(String tag:al.getTags()) {
				out.println(">" + tag);
				String alSeq = al.getAlignedSequence(tag);
				String seq = seqs.getSequenceNoGaps(tag);
				int c = 0; // counter in sequence
				for (int i = 0; i < alSeq.length(); i++) {
					if(alSeq.charAt(i) != '-') {
						out.print(seq.charAt(c));
						c++;
					} else {
						out.print('-');
					}
				}
				out.println();
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
		}
		return 1; // error
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length < 2) {
			System.out.println("Usage: java addSeq3Mcma <alignment> <sequences> [<outfile>]");
			System.exit(1);
		}
		
		String alFileName = args[0];
		String seqFileName = args[1];
		PrintStream out = System.out;
		if(args.length > 2) {
			out = new PrintStream(new FileOutputStream(args[2]));
		}
		
		File alFile = new File(alFileName);
		File seqFile = new File(seqFileName);
		int ret = addSequence(alFile, seqFile, out);
		System.exit(ret);
	}

}
