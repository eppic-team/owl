import proteinstructure.NodeSet;

/**
 * Program to test parsing of a residue selection string (e.g. 1-3,5,7-8)
 */
public class testResidueSelectionString {
	
	public static void main(String[] args) {

		String[] tp = {"1","10","10,1","10-9","10-10,5","1,6-7,8-9"};
		String[] tn = {"", "-",",","1,1,","1,-2","10-9-7",",4","4--5","1.5","1,,2"," "};
		int fp = 0, fn = 0;

		System.out.println("Positives:");
		for(String selStr:tp) {
			System.out.print(selStr + "\t\t");
			if(NodeSet.isValidSelectionString(selStr)) {
				System.out.println("ok");
			} else {
				System.out.println("error");
				fn++;
			}
		}
		System.out.println("Negatives:");
		for(String selStr:tn) {
			System.out.print(selStr + "\t\t");
			if(!NodeSet.isValidSelectionString(selStr)) {
				System.out.println("ok");
			} else {
				System.out.println("error");
				fn++;
			}
		}
		assert(fp == 0);
		assert(fn == 0);
		
		System.out.println("Parsing:");
		for(String selStr:tp) {
			System.out.print(selStr + "\t\t");
			NodeSet nodeSet = NodeSet.parseSelectionString(selStr);
			System.out.println(nodeSet);
		}
	}

}
