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
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.features.SiftsFeature;



/**
 * Connection class to get data from EBI's SIFTS resource mapping PDB chains to Uniprot entries.
 * This seems to be at the moment the gold standard for PDB to Uniprot mapping. The class loads
 * the results from a given URL or file pointer and caches the results so that subsequent queries
 * are done in O(1) time (for the price of memory consumption).
 * See {@link http://www.ebi.ac.uk/msd/sifts}
 * 
 * @author duarte, stehr
 */
public class SiftsConnection {
	
	public static final String PDB2UNIPROT_URL = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/text/pdb_chain_uniprot.lst";

	private static final Pattern URL_PATTERN = Pattern.compile("^\\w+://.*"); 
	
	private HashMap<String,ArrayList<SiftsFeature>> chain2uniprot;
	private HashMap<String,ArrayList<SiftsFeature>> uniprot2chain;	
	
	/**
	 * Constructs a SiftsConnection using the default online URL, parsing the SIFTS data and storing it.
	 * @throws IOException
	 */
	public SiftsConnection() throws IOException {
		this(PDB2UNIPROT_URL);
	}
	
	/**
	 * Constructs a SiftsConnection parsing the SIFTS data and storing it.
	 * To access the SIFTS data use {@link #getMappings(String, String)}
	 * @param pdb2uniprotURL a URL pointing to the SIFTS pdb to uniprot mapping file or simply 
	 * a path to a local file
	 * @throws IOException
	 */
	public SiftsConnection(String pdb2uniprotURL) throws IOException{
		chain2uniprot = new HashMap<String, ArrayList<SiftsFeature>>();
		uniprot2chain = new HashMap<String, ArrayList<SiftsFeature>>();		
		parsePdb2Uniprot(pdb2uniprotURL);
	}

	/**
	 * Parses the SIFTS table and stores pdb2uniprot and uniprot2pdb maps.
	 * @param fileURL a URL pointing to the SIFTS pdb to uniprot mapping file or 
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
			String[] fields = line.split("\t");
			String pdbCode = fields[0];
			String pdbChainCode = fields[1];
			String id = pdbCode+pdbChainCode;
			String uniprotId = fields[2];
			int cifBeg = Integer.parseInt(fields[4]);
			int cifEnd = Integer.parseInt(fields[5]);
			String pdbBeg = fields[6];
			String pdbEnd = fields[7];
			int uniBeg = Integer.parseInt(fields[8]);
			int uniEnd = Integer.parseInt(fields[9]);
			
			SiftsFeature siftsMapping = new SiftsFeature(pdbCode, pdbChainCode, uniprotId, cifBeg, cifEnd, pdbBeg, pdbEnd, uniBeg, uniEnd);
			
			// store pdb2uniprot record
			if (chain2uniprot.containsKey(id)) {
				chain2uniprot.get(id).add(siftsMapping);
			} else {
				ArrayList<SiftsFeature> ups = new ArrayList<SiftsFeature>();
				ups.add(siftsMapping);
				chain2uniprot.put(id, ups);				
			}
			// store uniprot2pdb record
			if (uniprot2chain.containsKey(uniprotId)) {
				uniprot2chain.get(uniprotId).add(siftsMapping);
			} else {
				ArrayList<SiftsFeature> ups = new ArrayList<SiftsFeature>();
				ups.add(siftsMapping);
				uniprot2chain.put(uniprotId, ups);				
			}			
		}
		br.close();
	}

	/**
	 * Gets a list of SiftsFeatures for the given PDB chain 
	 * @param pdbCode
	 * @param pdbChainCode
	 * @return
	 * @throws NoMatchFoundException if no matching Uniprot entry is found
	 */
	public List<SiftsFeature> getMappings(String pdbCode, String pdbChainCode) throws NoMatchFoundException{
		if (!chain2uniprot.containsKey(pdbCode+pdbChainCode)) 
			throw new NoMatchFoundException("No SIFTS mapping for PDB "+pdbCode+", chain "+pdbChainCode);
		return chain2uniprot.get(pdbCode+pdbChainCode);
	}
	
	/**
	 * Returns a list of SiftsFeatures for the given Uniprot ID.
	 * Warning: If this is not the primary ID, no results may be found.
	 * @param uniprotId
	 * @return the SiftsFeatures containing information about PDB chains associated with this Uniprot entry
	 * @throws NoMatchFoundException if no matching PDB chains are found
	 */
	public List<SiftsFeature> getUniprot2PdbMappings(String uniprotId) throws NoMatchFoundException {
		if (!uniprot2chain.containsKey(uniprotId)) 
			throw new NoMatchFoundException("No SIFTS mapping for UniprotID "+uniprotId);
		return uniprot2chain.get(uniprotId);
	}
	
	/**
	 * Gets the Collection of all mappings for all PDB chains found in the SIFTS repository
	 * @return
	 */
	public Collection<ArrayList<SiftsFeature>> getAllMappings() {
		return chain2uniprot.values();
	}
	
	/**
	 * Gets the total number of available pdb to uniprot mappings available in the SIFTS 
	 * database.
	 * @return
	 */
	public int getMappingsCount() {
		return chain2uniprot.size();
	}
}
