package owl.core.connections;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.structure.Pdb;
import owl.core.structure.features.EC;
import owl.core.structure.features.ECRegion;

/**
 * 
 * Connection class to parse the EC annotation data from UCL's Andrew C.R. Martin group.
 * It uses their own PDB to Uniprot mapping, see http://www.bioinf.org.uk/pdbsprotec/ 
 * and http://www.bioinf.org.uk/pdbsws/
 * 
 * NOTE: the SIFTS mapping should be used in preference to this one.
 * The SIFTS initiative is supposed to be the gold standard for PDB to Uniprot 
 * mapping. See http://www.ebi.ac.uk/pdbe/docs/sifts/
 * 
 * @see The {@link SiftsConnection} class.
 *
 */
public class ECConnection {

	private static final String PDB2EC_MAPPING_URL = "http://www.bioinf.org.uk/pdbsprotec/mapping.txt";
	private static final String PDB2EC_MAPPING_FILE = "/project/StruPPi/Databases/PDBSProtEC/mapping.txt";

	
	/**
	 * Parses local/remote PDB to EC mapping for the PDB code of the given Pdb object 
	 * and sets its EC annotation object 
	 * @param pdb
	 * @param online
	 * @throws IOException
	 */
	public static void parseEC(Pdb pdb, boolean online) throws IOException {

		EC ec = new EC();	
		ECRegion er = null;
		String startPdbResSer = "", endPdbResSer = "";
		BufferedReader in;
		String inputLine;
		Pattern p = Pattern.compile("^ \\d");

		if (online) {
			URL pdb2ecMapping = new URL(PDB2EC_MAPPING_URL);
			URLConnection p2e= pdb2ecMapping.openConnection();
			in = new BufferedReader(new InputStreamReader(p2e.getInputStream()));
		} else {
			File pdb2ecMapping = new File(PDB2EC_MAPPING_FILE);
			in = new BufferedReader(new FileReader(pdb2ecMapping));
		}

		while ((inputLine = in.readLine()) != null) { 
			Matcher m = p.matcher(inputLine);
			if (m.find()) {
				String curPdbCode = inputLine.substring(0,9).trim();
				String curPdbChainCode = (inputLine.charAt(11) == ' ')?Pdb.NULL_CHAIN_CODE:String.valueOf(inputLine.charAt(11));
				if (curPdbCode.equals(pdb.getPdbCode()) && curPdbChainCode.equals(pdb.getPdbChainCode())) {
					startPdbResSer = inputLine.substring(20,26).trim();
					endPdbResSer = inputLine.substring(27,33).trim();
					String id = inputLine.substring(43).trim();
					//System.out.println(curPdbCode+":"+curPdbChainCode+":"+startPdbResSer+"-"+endPdbResSer+":"+ec);
					er = new ECRegion(id, startPdbResSer, endPdbResSer, pdb.getResSerFromPdbResSer(startPdbResSer), pdb.getResSerFromPdbResSer(endPdbResSer));
					ec.add(er);
				}
			}
		}

		in.close();
		
		pdb.setEC(ec);

	}
}
