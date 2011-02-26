package owl.casp.benchmarking;

import java.io.File;

public class benchmarkOne {
	
	/*--------------------------------- main --------------------------------*/
	
	public static void main(String[] args) {
		
		if(args.length < 1) {
			System.out.println("Usage: benchmarkOne <base_dir> <group_suffix> <target>");
			System.out.println("e.g. benchmarkOne /project/StruPPi/CASP8/ TS183 T0318");
			System.exit(0);
		}
		
		String baseDir = args[0];
		String groupSuffix = args[1];
		String target = args[2];
		//String target = "T0399";
		
		File predDir = new File(new File(baseDir), "submitted");
		//String groupSuffix = "TS183";
		boolean eval3D = true;
		
		Benchmarking bm = new Benchmarking(null, new File(baseDir));
		int ret = bm.getTargetResult(predDir, groupSuffix, target, eval3D);
		if(ret == Benchmarking.NO_ERROR) {
			bm.printLastResultTable(System.out);
			System.out.println();
			bm.printLastResultSummary(System.out);
		}
	}
}
