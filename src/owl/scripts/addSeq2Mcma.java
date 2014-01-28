package owl.scripts;
import java.io.*;
import java.util.*;

import owl.core.sequence.alignment.AlignmentConstructionException;
import owl.core.sequence.alignment.MultipleSequenceAlignment;
import owl.core.sequence.Sequence;
import owl.core.util.FileFormatException;


/**
 * A script to add sequences to a raw MCMA alignment (which contains only the aligned positions
 * @author stehr
 */
public class addSeq2Mcma {

	/**
	 * @param alFile
	 * @param seqFile
	 * @param out
	 * @param aligned TRUE - sequences in seqFile are aligned; FALSE - sequences in seqFile are not aligned
	 * @return
	 */
	public static int addSequence(File alFile, File seqFile, PrintStream out, boolean aligned) {
		MultipleSequenceAlignment al = null;
		HashMap<String, String> seqs = new HashMap<String, String>();
		
		try {
			al = new MultipleSequenceAlignment(alFile.getAbsolutePath(), MultipleSequenceAlignment.FASTAFORMAT);
			if (aligned){
				MultipleSequenceAlignment seqsGap = new MultipleSequenceAlignment(seqFile.getAbsolutePath(), MultipleSequenceAlignment.FASTAFORMAT);
				for (int i = 0; i < seqsGap.getNumberOfSequences(); i++)
					seqs.put(seqsGap.getTagFromIndex(i), seqsGap.getSequenceNoGaps(seqsGap.getTagFromIndex(i)));
			}
			else{
				List<Sequence> seqsNoGap = Sequence.readSeqs(seqFile, null);
				for (int i=0; i < seqsNoGap.size(); i++)
					seqs.put(seqsNoGap.get(i).getName(), seqsNoGap.get(i).getSeq());
			}
			
			for(String tag:al.getTags()) {
				out.println(">" + tag);
				String alSeq = al.getAlignedSequence(tag);
				String seq = seqs.get(tag);
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
		} catch (FileFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (AlignmentConstructionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return 1; // error
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length < 3) {
			System.out.println("Usage: java addSeq2Mcma <alignment> <sequences> <aligned> [<outfile>]");
			System.exit(1);
		}
		
		String alFileName = args[0];
		String seqFileName = args[1];
		boolean aligned = Boolean.parseBoolean(args[2]);
		PrintStream out = System.out;
		if(args.length > 3) {
			out = new PrintStream(new FileOutputStream(args[3]));
		}
		
		File alFile = new File(alFileName);
		File seqFile = new File(seqFileName);
		int ret = addSequence(alFile, seqFile, out, aligned);
		System.out.println("addSeq2Mcma: Done");
		System.exit(ret);
	}

}
