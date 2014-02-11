package owl.scripts;
import java.io.File;

import owl.core.runners.CalcSurfVolRunner;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.TemplateList;


/**
 * Executable class to calculate volumes and surfaces for a list of pdb codes 
 * using the calc-volume and calc-surface programs
 * @author duarte
 *
 */
public class calcVolsAndSurfs {

	private static final String CIFREPODIR = "/path/to/mmCIF/gz/all/repo/dir";


	private static final String LIST = "/project/StruPPi/jose/optimal_reconstruction/model_pdbs.txt";
	private static final String CALCSURF_EXE = "/project/StruPPi/Software/libproteingeometry-2.3.1/bin/calc-surface";
	private static final String CALCVOL_EXE = "/project/StruPPi/Software/libproteingeometry-2.3.1/bin/calc-volume";
	
	public static void main(String[] args) throws Exception {
				
		String[] pdbIds = TemplateList.readIdsListFile(new File(LIST));

		for (String pdbId:pdbIds) {
			System.out.print(pdbId+"\t");
			String pdbCode = pdbId.substring(0,4);
			String pdbChainCode = pdbId.substring(4,5);
			
			

			File cifFile = new File(System.getProperty("java.io.tmpdir"),pdbCode+".cif");
			cifFile.deleteOnExit();
			PdbAsymUnit.grabCifFile(CIFREPODIR, null, pdbCode, cifFile, false);				
			PdbAsymUnit fullpdb = new PdbAsymUnit(cifFile);

			PdbChain pdb = fullpdb.getChain(pdbChainCode);
			System.out.printf("%10.3f\t",CalcSurfVolRunner.calcVolume(pdb,CALCVOL_EXE, ""));
			System.out.printf("%10.3f\n",CalcSurfVolRunner.calcSurface(pdb,CALCSURF_EXE, ""));
		}
			

	}

}
