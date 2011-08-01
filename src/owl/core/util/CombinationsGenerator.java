package owl.core.util;

import java.math.BigInteger;

/**
 * Class to generate combinations of n elements in groups of size r
 * 
 * Taken from http://www.merriampark.com/comb.htm
 *
 */
public class CombinationsGenerator {


	private int[] a;
	private int n;
	private int r;
	private BigInteger numLeft;
	private BigInteger total;


	public CombinationsGenerator (int n, int r) {
		if (r > n) {
			throw new IllegalArgumentException ();
		}
		if (n < 1) {
			throw new IllegalArgumentException ();
		}
		this.n = n;
		this.r = r;
		a = new int[r];
		BigInteger nFact = getFactorial (n);
		BigInteger rFact = getFactorial (r);
		BigInteger nminusrFact = getFactorial (n - r);
		total = nFact.divide (rFact.multiply (nminusrFact));
		reset ();
	}


	public void reset () {
		for (int i = 0; i < a.length; i++) {
			a[i] = i;
		}
		numLeft = new BigInteger (total.toString ());
	}


	/** 
	 * Return number of combinations not yet generated
	 */
	public BigInteger getNumLeft () {
		return numLeft;
	}

	/**
	 * Are there more combinations?
	 * @return
	 */
	public boolean hasMore () {
		return numLeft.compareTo (BigInteger.ZERO) == 1;
	}

	/**
	 * Return total number of combinations
	 * @return
	 */
	public BigInteger getTotal () {
		return total;
	}

	/**
	 * Compute factorial
	 * @param n
	 * @return
	 */
	private static BigInteger getFactorial (int n) {
		BigInteger fact = BigInteger.ONE;
		for (int i = n; i > 1; i--) {
			fact = fact.multiply (new BigInteger (Integer.toString (i)));
		}
		return fact;
	}

	/**
	 * Generate next combination (algorithm from Rosen p. 286)
	 * @return
	 */
	public int[] getNext () {

		if (numLeft.equals (total)) {
			numLeft = numLeft.subtract (BigInteger.ONE);
			return a;
		}

		int i = r - 1;
		while (a[i] == n - r + i) {
			i--;
		}
		a[i] = a[i] + 1;
		for (int j = i + 1; j < r; j++) {
			a[j] = a[i] + j - i;
		}

		numLeft = numLeft.subtract (BigInteger.ONE);
		return a;

	}

	/**
	 * Example
	 * @param args
	 */
	public static void main(String[] args) {
		String[] elements = {"a", "b", "c", "d", "e", "f", "g"};
		int[] indices;
		CombinationsGenerator x = new CombinationsGenerator (elements.length, 3);
		StringBuffer combination;
		while (x.hasMore ()) {
			combination = new StringBuffer ();
			indices = x.getNext ();
			for (int i = 0; i < indices.length; i++) {
				combination.append (elements[indices[i]]);
			}
			System.out.println (combination);
		}
	}

}
