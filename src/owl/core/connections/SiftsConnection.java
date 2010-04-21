package owl.core.connections;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.features.SiftsFeature;



/**
 * Connection class to get data from EBI's SIFTS resource mapping PDB identifiers to 
 * other databases: Uniprot, EC, taxonomy, GO, pubmed...
 * This seems to be at the moment the gold standard for PDB to Uniprot mapping
 * See {@link http://www.ebi.ac.uk/msd/sifts}
 * 
 * @author duarte
 *
 */
public class SiftsConnection {
	
	public static final String PDB2UNIPROT_URL = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/text/pdb_chain_uniprot.lst";

	private static final Pattern URL_PATTERN = Pattern.compile("^\\w+://.*"); 
	
	private HashMap<String,ArrayList<SiftsFeature>> chain2uniprot;
	
	/**
	 * Constructs a SiftsConnection parsing the SIFTS data and storing it.
	 * To access the SIFTS data use {@link #getMappings(String, String)}
	 * @param pdb2uniprotURL a URL pointing to the SIFTS pdb to uniprot mapping file or simply 
	 * a path to a local file
	 * @throws IOException
	 */
	public SiftsConnection(String pdb2uniprotURL) throws IOException{
		chain2uniprot = new HashMap<String, ArrayList<SiftsFeature>>();
		parsePdb2Uniprot(pdb2uniprotURL);
	}


	/**
	 * Parses the SIFTS pdb to uniprot mapping file and stores the resu.
	 * @param fileURL a URL pointing to the SIFTS pdb to uniprot mapping file or simply 
	 * a path to a local file
	 * @throws IOException
	 */
	private void parsePdb2Uniprot(String fileURL) throws IOException {
		Reader reader = null;
		 
		Matcher m = URL_PATTERN.matcher(fileURL);
		if (m.matches()) {
			// it is a URL
			URL pdb2enzymeURL = new URL(fileURL);
			URLConnection urlConn = pdb2enzymeURL.openConnection();
			reader = new InputStreamReader(urlConn.getInputStream());
		} else {
			// it is a file
			reader = new FileReader(new File(fileURL));
		}
		BufferedReader br = new BufferedReader(reader);
		String line;
		while ((line=br.readLine())!=null) {
			if (line.startsWith("PDB")) continue;
			String[] fields = line.split("\\s+");
			String pdbCode = fields[0];
			String pdbChainCode = fields[1];
			String id = pdbCode+pdbChainCode;
			String uniprotId = fields[2];
			int cifBeg = Integer.parseInt(fields[3]);
			int cifEnd = Integer.parseInt(fields[4]);
			int uniBeg = Integer.parseInt(fields[7]);
			int uniEnd = Integer.parseInt(fields[8]);
			

			SiftsFeature siftsMapping = new SiftsFeature(pdbCode, pdbChainCode, uniprotId, cifBeg, cifEnd, uniBeg, uniEnd);
			if (chain2uniprot.containsKey(id)) {
				chain2uniprot.get(id).add(siftsMapping);
			} else {
				ArrayList<SiftsFeature> ups = new ArrayList<SiftsFeature>();
				ups.add(siftsMapping);
				chain2uniprot.put(id, ups);				
			}
		}
		
		br.close();
		

	}

	/**
	 * Gets the Collection of SiftsFeatures for the given PDB chain 
	 * @param pdbCode
	 * @param pdbChainCode
	 * @return
	 * @throws NoMatchFoundException
	 */
	public Collection<SiftsFeature> getMappings(String pdbCode, String pdbChainCode) throws NoMatchFoundException{
		if (!chain2uniprot.containsKey(pdbCode+pdbChainCode)) 
			throw new NoMatchFoundException("No SIFTS mapping for PDB "+pdbCode+", chain "+pdbChainCode);
		return chain2uniprot.get(pdbCode+pdbChainCode);
	}
	
	/**
	 * Gets the Collection of all mappings for all PDB chains found in the SIFTS repository
	 * @return
	 */
	public Collection<ArrayList<SiftsFeature>> getAllMappings() {
		return chain2uniprot.values();
	}
}
