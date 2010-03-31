package casp.benchmarking;

import java.io.File;

public class benchmarkOne {
	
	/*--------------------------------- main --------------------------------*/
	
	public static void main(String[] args) {
		
		if(args.length < 1) {
			System.out.println("Usage: Benchmarking <target>");
			System.exit(0);
		}
		
		String target = args[0];
		//String target = "T0399";
		
		File predDir = new File("/project/StruPPi/CASP8/submitted");
		String groupSuffix = "TS183";
		boolean eval3D = true;
		
		Benchmarking bm = new Benchmarking(null);
		int ret = bm.getTargetResult(predDir, groupSuffix, target, eval3D);
		if(ret == Benchmarking.NO_ERROR) {
			bm.printLastResultTable(System.out);
			System.out.println();
			bm.printLastResultSummary(System.out);
		}
	}
}
