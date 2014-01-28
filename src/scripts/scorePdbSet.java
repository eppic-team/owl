package scripts;
import java.io.File;
import java.io.IOException;

import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbLoadException;
import owl.core.structure.TemplateList;
import owl.core.util.FileFormatException;
import owl.decoyScoring.AtomTypeScorer;
import owl.decoyScoring.ResTypeScorer;
import owl.decoyScoring.Scorer;


/**
 * Executable class to score a given list of PDB codes in a list file (first argument passed)
 * Useful to find outliers in the scoring of native structures.
 * @author duarte
 *
 */
public class scorePdbSet {

	private static final String CIFREPODIR = "/path/to/mmCIF/gz/all/repo/dir";

	private static final File atomScMatFile = new File("/project/StruPPi/jose/emp_potential/scoremat.atom.cullpdb20");
	private static final File resScMatFile = new File("/project/StruPPi/jose/emp_potential/scoremat.res.cullpdb20");


	public static void main(String[] args) throws IOException, FileFormatException  {

		File listFile = new File(args[0]);
			
		String[] pdbIds = TemplateList.readIdsListFile(listFile);

		AtomTypeScorer atomScorer = new AtomTypeScorer(atomScMatFile);
		ResTypeScorer resScorer = new ResTypeScorer(resScMatFile);
		
		for (String pdbId:pdbIds) {
			String pdbCode = pdbId.substring(0,4);
			String pdbChainCode = pdbId.substring(4,5);
		
			try {
				
				File cifFile = new File(System.getProperty("java.io.tmpdir"),pdbCode+".cif");
				cifFile.deleteOnExit();
				PdbAsymUnit.grabCifFile(CIFREPODIR, null, pdbCode, cifFile, false);				
				PdbAsymUnit fullpdb = new PdbAsymUnit(cifFile);

				PdbChain pdb = fullpdb.getChain(pdbChainCode);
				if (!Scorer.isValidPdb(pdb)) {
					continue;
				}
				System.out.printf("%5s\t%4d\t%7.2f\t%7.2f\n",pdbId,pdb.getStdAaObsLength(),resScorer.scoreIt(pdb),atomScorer.scoreIt(pdb));
			} catch (PdbLoadException e) {
				System.err.println("Couldn't load "+pdbId);
				continue;
			} 
			
		}
		
		
	}

}
