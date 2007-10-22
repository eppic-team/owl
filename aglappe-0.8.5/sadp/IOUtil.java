package sadp;
import java.io.*;

/**
 * This class parses contacts maps from files. Each contact map must be stored
 * in a file in the following format:
 * 
 * The first line contains the number of nodes. The following lines contain four
 * values n m v w, where the first pair (n,m) are nodes that form a contact
 * between the respective nodes and the second pair (u, v) are weights of the
 * contacts. In the current form, weights are read but ignored.
 */
public class IOUtil {

    private static final String ERR_MSG = "ERROR: Corrupted file: ";


    /**
     * Reads a contact maps given in the specified file of the specified format
     * and returns a contact map.
     */
    public static ContactMap read(String f) {

	try {
	    FileReader fr = new FileReader(f);
	    StreamTokenizer st = new StreamTokenizer(fr);
	    st.parseNumbers();

	    // read number of nodes
	    st.nextToken();
	    if (st.ttype == StreamTokenizer.TT_EOF) {
		return null;
	    }
	    if (st.ttype != StreamTokenizer.TT_NUMBER) {
		System.err.println(ERR_MSG + f);
		System.exit(0);
	    }
	    int n = (int) st.nval;

	    // initialize contact matrix
	    boolean[][] A = new boolean[n][n];

	    // read contacts
	    int i, j;
	    while (st.nextToken() != StreamTokenizer.TT_EOF) {

		// first node
		if (st.ttype != StreamTokenizer.TT_NUMBER) {
		    System.err.println(ERR_MSG + f);
		    System.exit(0);
		}
		i = (int) st.nval;

		// second node
		st.nextToken();
		if (st.ttype != StreamTokenizer.TT_NUMBER) {
		    System.err.println(ERR_MSG + f);
		    System.exit(0);
		}
		j = (int) st.nval;

		A[i][j] = true;
		A[j][i] = true;

		// ignore weights
		for (int k = 0; k < 2; k++) {
		    st.nextToken();
		    if (st.ttype != StreamTokenizer.TT_NUMBER) {
			System.err.println(ERR_MSG + f);
			System.exit(0);
		    }
		}
	    }
	    ContactMap cm = new ContactMap(A);
	    cm.setFileName(f);
	    return cm;
	} catch (IOException e) {
	    System.err.print("Exception while reading from file ");
	    System.err.println(f + ". " + e.getMessage());
	}
	return null;
    }

    public static void main(String[] args) {

	// String source = "data/family16X/data/bS/4.1/";
	// String target = "data/family16X_L/data/bS/4.1/";
	// IOUtil.convert2LanciaFormat(source, target);
    }
}
