package owl.core.sequence;

import owl.core.structure.AminoAcid;

/**
 * A codon. 
 * 
 * Note that by convention we always use lower case for nucleotides.
 * 
 * @author duarte_j
 *
 */
public class Codon {

	// codes for ambiguous nucleotides (from wikipedia/iupac http://en.wikipedia.org/wiki/Nucleic_acid_sequence)
	private static final char AMB_ACGT = 'n';
	private static final char AMB_AC   = 'm';
	private static final char AMB_AT   = 'w';
	private static final char AMB_TC   = 'y';
	
	
	private char[] codon;
	
	private int hash; // the cached hash
	
	public Codon(char first, char second, char third){
		if (!isLowerCase(first,second,third)){
			throw new IllegalArgumentException("Codon must be lower case! Codon: "
					+Character.toString(first)+Character.toString(second)+Character.toString(third));
		}
		codon = new char[3];
		codon[0]=first;
		codon[1]=second;
		codon[2]=third;
	}
	
	public Codon(String codon){
		if (codon.length()!=3) throw new IllegalArgumentException("Codon length is not 3! Codon: "+codon);
		this.codon = new char[3];
		this.codon[0] = codon.charAt(0);
		this.codon[1] = codon.charAt(1);
		this.codon[2] = codon.charAt(2);
		if (!isLowerCase(this.codon[0],this.codon[1],this.codon[2])){
			throw new IllegalArgumentException("Codon must be lower case! Codon: "
					+codon);
		}		
	}
	
	private boolean isLowerCase(char first, char second, char third) {
		return (Character.isLowerCase(first) && 
				Character.isLowerCase(second) && 
				Character.isLowerCase(third));
	}
	
	public int hashCode() {
		// taken from java's String hashCode implementation
		int h = hash;
		if (h == 0) {
			for (int i=0; i<3; i++)	{
				h = 31 * h + codon[i];
			}
			hash = h;
		}
		return h;
	}
	
	public boolean equals(Object o) {
		if (o==this) return true;
		if (! (o instanceof Codon)) {
			return false;
		}
		Codon other = (Codon) o;
		if (other.codon[0]==this.codon[0] && 
				other.codon[1]==this.codon[1] && 
				other.codon[2]==this.codon[2]) {
			return true;
		}
		return false;
	}
	
	public String toString() { 
		return new String(codon);
	}
	
	public boolean isAmbiguous() {
		for (char c:codon) {
			if (c==AMB_ACGT || c==AMB_AC || c== AMB_AT || c==AMB_TC) return true;
		}
		return false;
	}
	
	public AminoAcid getTranslation(GeneticCodeType gct) throws TranslationException {
		return Translator.translate(gct, this);
	}
	
	// testing
	public static void main(String[] args) {
		Codon c1 = new Codon('a','c','t');
		Codon c2 = new Codon('a','c','t');
		Codon c3 = new Codon('a','c','g');
		Codon c4 = new Codon('c','a','t');
		Codon c5 = new Codon('a','t','c');
		printComparison(c1, c1);
		printComparison(c1, c2);
		printComparison(c1, c3);
		printComparison(c1, c4);
		printComparison(c1, c5);
	}
	
	private static void printComparison(Codon c1, Codon c2) {
		if (c1==c2) System.out.println("same object");
		System.out.println(c1.toString()+" "+c2.toString());
		System.out.println("Equals: "+c1.equals(c2));
		System.out.println("Hash:   "+c1.hashCode()+"\t"+c2.hashCode());
		System.out.println();
	}
	
}
