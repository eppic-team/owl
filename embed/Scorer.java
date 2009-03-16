package embed;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Arrays;

import proteinstructure.FileRIGraph;
import proteinstructure.GraphFileFormatError;
import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbasePdb;
import proteinstructure.RIGraph;

import tools.MySQLConnection;
import tools.RegexFileFilter;

/**
 * Class containing static methods to score contact map subsets.
 * 
 * @author duarte
 *
 */
public class Scorer {

	private static File dir = new File("/project/LitNet/Essence/Data/Contactmaps"); 


	
	/*------------------------- privates  -------------------------------*/
	
	/**
	 * Infer bounds for all pairs from a spare bounds matrix via the triangle inequality
	 * @param sparseBounds
	 * @return
	 */
	private static Bound[][] inferAllBounds(Bound[][] sparseBounds) {
		BoundsSmoother bs = new BoundsSmoother(sparseBounds);
		return bs.getInitialBoundsAllPairs();
	}
	
	/*------------------------- publics  -------------------------------*/
	
	/**
	 * Evaluates the error value for the given sparseBounds matrix (a subset of the contact map) 
	 * by first inferring bounds for all pairs through the triangle inequality and then 
	 * measuring how well the all pairs matrix fits the contact map.  
	 * Thus a high error value means the given sparseBounds matrix has low information content
	 * about the rest of the contact map and low error value that the matrix has high 
	 * information content.
	 * The error measure consists of the sum of the deviation of the upper bounds 
	 * to the ones in the contact map: sum(max(0,(u'i-ui))), i.e. if upper bound in 
	 * given matrix is below the one in the contact map there is no penalty.
	 * @param sparseBounds
	 * @param fullContactMap
	 * @return
	 */
	public static double getCMError(Bound[][] sparseBounds, RIGraph fullContactMap) {
		Bound[][] cmapBounds = Reconstructer.convertRIGraphToBoundsMatrix(fullContactMap);
		// infer bounds for all pairs through triangle inequality
		Bound[][] bounds = inferAllBounds(sparseBounds);
		double sumDev = 0;
		for (int i=0;i<bounds.length;i++) {
			for (int j=i+1;j<bounds.length;j++) {
				if (cmapBounds[i][j]!=null) {
					sumDev += Math.max(0,bounds[i][j].upper-cmapBounds[i][j].upper);
				}
			}
		}
		return sumDev/(double) fullContactMap.getFullLength(); 
	}
	
	/**
	 * Evaluates the error value for the given contact map subset 
	 * by first inferring bounds for all pairs through the triangle inequality and then 
	 * measuring how well the all pairs matrix fits the contact map.  
	 * Thus a high error value means the given sparseBounds matrix has low information content
	 * about the rest of the contact map and low error value that the matrix has high 
	 * information content  
	 * @param subset
	 * @param fullContactMap
	 * @return
	 */
	public static double getCMError(RIGraph subset, RIGraph fullContactMap) {
		Bound[][] subsetBounds = Reconstructer.convertRIGraphToBoundsMatrix(subset);
		return getCMError(subsetBounds, fullContactMap);
	}

	/**
	 * Evaluates the error value for the given contact map file containing a subset. 
	 * It will be compared to the contact map taken from the pdbCode/contactType/cutoff
	 * present in the contact map file.
	 * @param rigFile
	 * @param pdbaseDb
	 * @param conn
	 * @return
	 * @throws GraphFileFormatError
	 * @throws IOException
	 * @throws PdbCodeNotFoundError
	 * @throws SQLException
	 * @throws PdbLoadError
	 */
	public static double getCMError(File rigFile, String pdbaseDb, MySQLConnection conn) throws GraphFileFormatError, IOException, PdbCodeNotFoundError, SQLException, PdbLoadError {
		RIGraph graph = new FileRIGraph(rigFile.getAbsolutePath());
		int numContSubset = graph.getEdgeCount();
		String pdbCode = graph.getPdbCode();
		String pdbChainCode = graph.getPdbChainCode();
		String ct = graph.getContactType();
		double cutoff = graph.getCutoff();
		
		Pdb pdb = new PdbasePdb(pdbCode, pdbaseDb, conn);
		pdb.load(pdbChainCode);
		RIGraph fullContactMap = pdb.get_graph(ct, cutoff);
		
		int numTotalCont = fullContactMap.getEdgeCount();
		System.out.printf("size: %4.2f\n",(double)numContSubset/(double)numTotalCont);
		
		return getCMError(graph, fullContactMap);
		
	}
	
	/*-------------------------- main  -------------------------------*/

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		
		File[] selSubsets = dir.listFiles(new RegexFileFilter(".*_all\\.graph"));
		File[] rndSubsets = dir.listFiles(new RegexFileFilter(".*_rnd\\.graph"));
		
		Arrays.sort(selSubsets);
		Arrays.sort(rndSubsets);
		
		MySQLConnection conn = new MySQLConnection();
		
		for (int i=0;i<selSubsets.length;i++) {
			File selSubset = selSubsets[i];
			File rndSubset = rndSubsets[i];
			System.out.printf("%s\t"+Distiller.SCORE_PRINT_FORMAT+"\n",selSubset.getName(),getCMError(selSubset, "pdbase", conn));
			System.out.printf("%s\t"+Distiller.SCORE_PRINT_FORMAT+"\n",rndSubset.getName(),getCMError(rndSubset, "pdbase", conn));
			System.out.println();
		}
		
	}

	
}
