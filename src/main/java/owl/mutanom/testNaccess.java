package owl.mutanom;

import java.io.File;
import java.io.IOException;

import owl.core.runners.NaccessRunner;
import owl.core.structure.PdbChain;


public class testNaccess {

	private static final String NACCESS_EXECUTABLE = "/project/StruPPi/bin/naccess";
	private static final String NACCESS_PARAMETERS = "";
	
	public static void main(String[] args) throws IOException {
		
		PdbChain p1 = PdbChain.readStructureOrExit("3gftA");
		PdbChain p2 = PdbChain.readStructureOrExit("2vukA");
		
		NaccessRunner nar = new NaccessRunner(new File(NACCESS_EXECUTABLE), NACCESS_PARAMETERS);
		nar.runNaccess(p1);
		nar.runNaccess(p2);
		System.out.println("done.");
	}

}
