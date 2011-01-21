package owl.core.sequence.alignment;

import java.util.Iterator;
import java.util.TreeSet;

import owl.core.util.IntPairComparator;


import edu.uci.ics.jung.graph.util.Pair;

/**
 * Instances of this class construct pairwise alignments based on a matching
 * between position in the first sequence onto positions in the second
 * sequence.
 * <p>
 * Execute method {@link main(String)} to get an impression of the functionality of this class. 
 * 
 * @author Lars Petzold
 * @see owl.sadp.SADP
 * */
public class PairwiseAlignmentConverter {

	private static char    GAP       	= '-';
	private static char    UNDEF     	= 'X';
	private MultipleSequenceAlignment      ali       	= null;
	private Iterator<Pair<Integer>> it  = null;
	private int[]          lengths   	= new int[2];
	private String[]       sequences 	= new String[2];
	private String[]       tags     	= new String[2];
	private int            offset    	= 0;

	/**
	 * Constructs sequence alignment for the given edge set.
	 * The tags are used as sequence identifiers. The length parameters as the
	 * length of the pseudo-sequence which consists only of 'X' characters
	 * denoting undefined amino acids. Use constructors taking the sequences
	 * itself to avoid pseudo-sequences.
	 * @param it   edge iterator
	 * @param len1 length of the first pseudo-sequence (assumed length of the original sequence)
	 * @param len2 length of the second pseudo-sequence (assumed length of the original sequence)
	 * @param tag1 sequence name of the first sequence
	 * @param tag2 sequence name of the second sequence
	 * */
	public PairwiseAlignmentConverter( Iterator<Pair<Integer>> it, int len1, int len2, String tag1, String tag2, int fi ) 
	throws AlignmentConstructionException {
		this(it,len1,len2,null,null,tag1,tag2,fi);
	}

	/**
	 * The actual constructor. Each public constructor calls this constructor.
	 * @param it   edge iterator
	 * @param len1 length of the first sequence
	 * @param len2 length of the second sequence
	 * @param seq1 first sequence
	 * @param seq2 second sequence
	 * @param tag1 name of the first sequence
	 * @param tag2 name of the second sequence
	 * */
	private PairwiseAlignmentConverter( Iterator<Pair<Integer>> it, int len1, int len2, String seq1, String seq2, String tag1, String tag2, int fi) 
	throws AlignmentConstructionException {

		sequences[0] = seq1;
		sequences[1] = seq2;
		tags[0]      = tag1;
		tags[1]      = tag2;
		lengths[0]   = len1;
		lengths[1]   = len2;
		this.it      = it;
		offset       = -fi; // flip the sign to achieve a index conversion where 0 is the first index

		buildAlignment();
	}

	public PairwiseAlignmentConverter( Iterator<Pair<Integer>> it, int len1, String seq2, String tag1, String tag2, int fi)
	throws AlignmentConstructionException {
		this(it,len1,seq2.length(),null,seq2,tag1,tag2,fi);
	}

	public PairwiseAlignmentConverter( Iterator<Pair<Integer>> it, String seq1, int len2, String tag1, String tag2, int fi)
	throws AlignmentConstructionException {
		this(it,seq1.length(),len2,seq1,null,tag1,tag2,fi);
	}

	public PairwiseAlignmentConverter( Iterator<Pair<Integer>> it, String seq1, String seq2, String tag1, String tag2, int fi ) 
	throws AlignmentConstructionException {
		this(it,seq1.length(),seq2.length(),seq1,seq2,tag1,tag2,fi);
	}

	/**
	 * Builds the alignment based on the preset members.
	 */
	private void buildAlignment( ) throws AlignmentConstructionException {

		// gapped sequences for the respective contact maps
		StringBuffer s1 = new StringBuffer(Math.max(lengths[0],lengths[1]));
		StringBuffer s2 = new StringBuffer(s1.capacity());

		// sequence positions
		int  prev1 = -1;
		int  prev2 = -1;
		int  cur1  =  0;
		int  cur2  =  0;
		Pair<Integer> e     =  null;

		while( it.hasNext() ) {

			e    = it.next();
			cur1 = e.getFirst()+offset;
			cur2 = e.getSecond()+offset;

			// puts columns ("non-trusted" MATCHes and MISMATCHes) between two
			// consecutive matching edges
			putColumnsBetweenMatches(prev1+1,cur1-prev1-1,prev2+1,cur2-prev2-1,s1,s2);

			// add recently found MATCH (cur1,cur2)
			putMatch(cur1,cur2,s1,s2);

			// ... and update position counters
			prev1 = cur1;
			prev2 = cur2;			
		}

		// insert trailing GAPs and probably some "non-trusted" MATCHes
		putColumnsBetweenMatches(cur1+1, lengths[0]-cur1-1, cur2+1, lengths[1]-cur2-1, s1, s2);

		// compose the name to sequence mapping 
		String[] seqs = {s1.toString(), s2.toString()};
		ali = new MultipleSequenceAlignment(tags,seqs);
	}

	/**
	 * Retrieves alignment.
	 * 
	 * @return pairwise alignment based on the set matching edges as being
	 *         passed to any constructor. 
	 */
	public MultipleSequenceAlignment getAlignment() {

		return ali;
	}

	/**
	 * Puts a character into the given gapped sequence.
	 * <p>
	 * If the actual sequence is not provided an X is inserted instead denoting
	 * a non standard sequence item.
	 * 
	 * @param pos       index of character in the original sequence to be set
	 * @param s         sequence already containing GAP-characters
	 * @param whichOne  toggles the sequences, 1 -> first sequence, 2 -> second
	 *                  sequence
	 */
	private void putCharacter( int pos, StringBuffer s, int whichOne ) {

		if( sequences[whichOne-1] != null ) {
			s.append(sequences[whichOne-1].charAt(pos));
		} else {
			s.append(UNDEF);
		}
	}

	/**
	 * In one or both sequences serveral positions might have been skipped for
	 * some reasons. This method therefore applies a greedy strategy to first
	 * introduce "non-trusted" MATCHes and secondly to align the remaining
	 * positions -- in at most one sequence -- to GAPs.
	 */
	private void putColumnsBetweenMatches( int beg1, int len1, int beg2, int len2, StringBuffer s1, StringBuffer s2 ) {

		int pos = 0;

		// we do model preceding gaps first
		if( len1 > 0 ) {
			if( len2 > 0 ) {
				// model non-trusted MATCHes first
				for( pos = 0; pos < Math.min(len1,len2); ++pos ) {
					putMatch(beg1+pos,beg2+pos,s1,s2);
				}
				// secondly, introduce GAPs in one of the sequences
				if( len1 > len2 ) {
					for( ; pos < len1; ++pos ) {
						putGapInSecond(beg1+pos,s1,s2);
					}
				} else {
					for( ; pos < len2; ++pos ) {
						putGapInFirst(beg2+pos,s1,s2);
					}				
				}
			} else {
				// introduce GAPS in the second sequence
				for( pos = 0; pos < len1; ++pos ) {
					putGapInSecond(beg1+pos,s1,s2);
				}
			}
		} else if (len2 > 0) {
			// introduce GAPs in the first sequence
			for( pos = 0; pos < len2; ++pos ) {
				putGapInFirst(beg2+pos,s1,s2);
			}				
		}
	}

	/**
	 * Puts a GAP in the first sequence and a character in the second sequence.
	 * 
	 * @param pos2  character in the second sequence to be set
	 * @param s1    first sequence already containing GAP-characters
	 * @param s2    second sequence already containing GAP-characters
	 */
	private void putGapInFirst( int pos2, StringBuffer s1, StringBuffer s2 ) {

		s1.append(GAP);
		putCharacter(pos2,s2,2);
	}

	/**
	 * Puts a character in the first sequence and a GAP in the second sequence.
	 * 
	 * @param pos1  index of character in the second sequence to be set
	 * @param s1    first sequence already containing GAP-characters
	 * @param s2    second sequence already containing GAP-characters
	 */
	private void putGapInSecond( int pos1, StringBuffer s1, StringBuffer s2 ) {

		putCharacter(pos1,s1,1);
		s2.append(GAP);
	}

	/**
	 * Puts a MATCH which means that pos1'th character is to be inserted in s1 and the pos2'th in s2.
	 * 
	 * @param pos1  current position in the first sequence
	 * @param pos2  current position in the second sequence
	 * @param s1    first sequence already containing GAP-characters
	 * @param s2    second sequence already containing GAP-characters 
	 */
	private void putMatch( int pos1, int pos2, StringBuffer s1, StringBuffer s2 ) {

		putCharacter(pos1,s1,1);
		putCharacter(pos2,s2,2);
	}

	/**
	 * Tests PairwiseAlignmentConverter and gives Examples how to use it.
	 * */
	public static void main( String args[]) {

		// MATCHING:
		//     0   1 2 3     4 5 6 7    
		// - - o - o o o - - o o o o 
		//     |   | | :     | :   
		// o o o o o o o o o o o - -
		// 0 1 2 3 4 5 6 7 8 9 1
		//                     0 
		//
		// LEGEND: 
		//  '-' -> GAP
		//  'o' -> position
		//  '|' -> MATCH
		//  ':' -> "non trusted" MATCH (does not have a reference in the edge set)

		TreeSet<Pair<Integer>> edges = new TreeSet<Pair<Integer>>(new IntPairComparator());
		edges.add(new Pair<Integer>(0,2));
		edges.add(new Pair<Integer>(1,4));
		edges.add(new Pair<Integer>(2,5));
		edges.add(new Pair<Integer>(4,9));

		String seq1 = "ABCDEFGH";
		String seq2 = "ABCDEFGHIJK";

		String tag1 = "seq1";
		String tag2 = "seq2";

		PairwiseAlignmentConverter[] pab = new PairwiseAlignmentConverter[4];

		try {       
			pab[0] = new PairwiseAlignmentConverter( edges.iterator(), seq1,          seq2,          tag1, tag2, 0 );
			pab[1] = new PairwiseAlignmentConverter( edges.iterator(), seq1.length(), seq2,          tag1, tag2, 0 );
			pab[2] = new PairwiseAlignmentConverter( edges.iterator(), seq1,          seq2.length(), tag1, tag2, 0 );
			pab[3] = new PairwiseAlignmentConverter( edges.iterator(), seq1.length(), seq2.length(), tag1, tag2, 0 );
		} catch(Exception e) {
			System.err.println(e.getMessage());
			System.exit(-1);
		}

		for( int i=0; i<pab.length; ++i ) {

			MultipleSequenceAlignment a = pab[i].getAlignment();

			System.out.println("\n" + i + "\n``````````````````");

			for( String seqTag:a.getTags() ) {
				System.out.println(seqTag + ": " + a.getAlignedSequence(seqTag));
			}
		}        
	}

}
