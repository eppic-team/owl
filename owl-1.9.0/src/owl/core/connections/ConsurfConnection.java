package owl.core.connections;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.structure.Pdb;
import owl.core.structure.Residue;


/**
 * Class to get conservation data from the consurf server or from local copies 
 * of the consurf files
 *   
 * @author duarte_j
 *
 */
public class ConsurfConnection {
	
	// NOTE: the URL doesn't seem to work anymore - Jose 21.07.2009
	private static final String CONSURF_URL_PREFIX = "http://consurf.tau.ac.il/results/";
	private static final String CONSURF_DIR = "/project/StruPPi/Databases/ConSurf-HSSP/ConservationGrades";


	/**
	 * Constructs a new ConsurfConnection. Use subsequently {@link #getConsurfDataFromWeb(Pdb, String)} or
	 * {@link #getConsurfDataLocal(Pdb, String)} depending where one wants the data to be read from.
	 */
	public ConsurfConnection() {
		
	}
	
	/**
	 * Parses the consurf data from the consurf web site, updating the residues of the
	 * given Pdb object with the consurf data 
	 * @param pdb
	 * @param consurfURLPrefix a consurf URL prefix or null (default will be taken)
	 * @return
	 * @throws IOException
	 */
	public int getConsurfDataFromWeb(Pdb pdb, String consurfURLPrefix) throws IOException {
		if (consurfURLPrefix==null) consurfURLPrefix = CONSURF_URL_PREFIX;
		
		String pdbCode = pdb.getPdbCode();
		String pdbChainCode = pdb.getPdbChainCode();

		// TODO: Check if url exists and if not do the same as for the offline case 
		URL consurfhssp = new URL(consurfURLPrefix+"HSSP_ML_"+pdbCode+(pdbChainCode.equals(Pdb.NULL_CHAIN_CODE)?"_":pdbChainCode)+"/pdb"+pdbCode+".gradesPE");
		URLConnection ch = consurfhssp.openConnection();
		return parseConsurfData(pdb, new BufferedReader(new InputStreamReader(ch.getInputStream())));
	}
	
	/**
	 * Parses the consurf data from a consurf local directory, updating the residues of the
	 * given Pdb object with the consurf data 
	 * @param pdb
	 * @param consurfDir a local directory with files with consurf data or null (default will be taken)
	 * @return
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public int getConsurfDataLocal(Pdb pdb, String consurfDir) throws FileNotFoundException, IOException {
		if (consurfDir==null) consurfDir = CONSURF_DIR;
		
		String pdbCode = pdb.getPdbCode();
		String pdbChainCode = pdb.getPdbChainCode();

		File consurfhssp = new File(consurfDir,pdbCode+(pdbChainCode.equals(Pdb.NULL_CHAIN_CODE)?"_":pdbChainCode)+".grades");
		if (!consurfhssp.exists() && pdbChainCode.equals("A")) {
			System.out.println("consurf");
			consurfhssp = new File(consurfDir,pdbCode+"_.grades");
		}
		return parseConsurfData(pdb,new BufferedReader(new FileReader(consurfhssp)));		
	}
	
	/**
	 * Parses consurf data from the given BufferedReader updating the residues of the given
	 * Pdb object with the consurf data
	 * @return the number of mistakes found while parsing consurf results
	 * @throws IOException when something goes wrong with parsing the BufferedReader
	 */
	private int parseConsurfData(Pdb pdb, BufferedReader in) throws IOException {
		
		String inputLine;
		Pattern p = Pattern.compile("^\\s+\\d+");
		int lineCount = 0;
		int consurfHsspMistakes = 0;
			
		Integer[] ressers = new Integer[pdb.getObsLength()];
		pdb.getAllSortedResSerials().toArray(ressers);

		while ((inputLine = in.readLine()) != null) { 
			Matcher m = p.matcher(inputLine);
			if (m.find()) {
				lineCount++;
				int resser = ressers[lineCount-1];
				String[] fields = inputLine.split("\\s+");
				String pdbresser = fields[3].equals("-")?"-":fields[3].substring(3, fields[3].indexOf(':'));
				if (fields[2].equals(String.valueOf(pdb.getResidue(resser).getAaType().getOneLetterCode())) &
						(pdbresser.equals("-") | pdbresser.equals(pdb.getPdbResSerFromResSer(resser)))) {
					if (pdb.containsResidue(resser)) {
						Residue residue = pdb.getResidue(resser);
						residue.setConsurfScore(Double.valueOf(fields[4]));
						residue.setConsurfColor(Integer.valueOf(fields[5]));
					} else {
						consurfHsspMistakes++;
					}
				} else {
					consurfHsspMistakes++;
				}
			}
		}
		in.close();

		// checking how many were actually assigned
		int assigned = 0;
		for (int resser:pdb.getAllSortedResSerials()) {
			Residue residue = pdb.getResidue(resser);
			if (residue.getConsurfScore()!=null)	assigned++;
		}
        consurfHsspMistakes += Math.abs(pdb.getObsLength() - assigned);
        if (consurfHsspMistakes > 0) {
    		for (int resser:pdb.getAllSortedResSerials()) {
    			Residue residue = pdb.getResidue(resser);
        		residue.setConsurfScore(null);
        		residue.setConsurfColor(null);
        	}
		}

		return consurfHsspMistakes;

	}

}
