import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import owl.core.structure.PdbCodeNotFoundError;
import owl.core.structure.PdbLoadError;
import owl.core.structure.PdbasePdb;
import owl.core.structure.TemplateList;
import owl.core.util.FileFormatError;
import owl.core.util.MySQLConnection;
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

	private static final String PDBASEDB = "pdbase";
	private static final File atomScMatFile = new File("/project/StruPPi/jose/emp_potential/scoremat.atom.cullpdb20");
	private static final File resScMatFile = new File("/project/StruPPi/jose/emp_potential/scoremat.res.cullpdb20");


	public static void main(String[] args) throws SQLException, IOException, FileFormatError  {

		File listFile = new File(args[0]);
		
		MySQLConnection conn = new MySQLConnection();
		
		String[] pdbIds = TemplateList.readIdsListFile(listFile);

		AtomTypeScorer atomScorer = new AtomTypeScorer(atomScMatFile);
		ResTypeScorer resScorer = new ResTypeScorer(resScMatFile);
		
		for (String pdbId:pdbIds) {
			String pdbCode = pdbId.substring(0,4);
			String pdbChainCode = pdbId.substring(4,5);
		
			try {
				PdbasePdb pdb = new PdbasePdb(pdbCode,PDBASEDB,conn);
				pdb.load(pdbChainCode);
				if (!Scorer.isValidPdb(pdb)) {
					continue;
				}
				System.out.printf("%5s\t%4d\t%7.2f\t%7.2f\n",pdbId,pdb.getObsLength(),resScorer.scoreIt(pdb),atomScorer.scoreIt(pdb));
			} catch (PdbLoadError e) {
				System.err.println("Couldn't load "+pdbId);
				continue;
			} catch (PdbCodeNotFoundError e) {
				System.err.println("Couldn't find "+pdbId);
				continue;				
			}
			
		}
		
		
	}

}
