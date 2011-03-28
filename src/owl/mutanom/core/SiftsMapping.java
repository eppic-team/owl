package owl.mutanom.core;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;


/**
 * Queries the SIFTS uniprot to pdb mapping and stores the results. The mapping can be local or the default FTP location
 * TODO: Merge with owl.core.connections.SiftsConnection and owl.core.features.SiftsFeature
 * @author stehr
 */
public class SiftsMapping {
	
	/*------------------------------ constants ------------------------------*/
	public static final String SIFTS_URL = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/text/pdb_chain_uniprot.lst"; // TODO: move to properties file
	
	/*--------------------------- member variables --------------------------*/
	String pdbCode;
	String chainCode;
	String uniprotId;
	int cifBeg;			// beginning residue in mmCIF sequence
	int cifEnd;			// end residue in mmCIF sequence
	int pdbBeg;			// beginning residue in pdb sequence
	int pdbEnd;			// end residue in pdb sequence
	int uniBeg;			// beginning residue in Uniprot sequence
	int uniEnd;			// end residue in Uniprot sequence
	
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Queries the SIFTS mapping for the given Pdb and Chain code and stores the results in the member variables.
	 * If localSiftsFile is null, queries the default FTP location.
	 * @throws IOException if the SIFTS source was not found or an error occured while reading from it
	 * @throws SiftsMappingNotFoundException if the SIFTS source was found but contained no mapping for the given pdb and chain code
	 */
	public SiftsMapping(String pdbCode, String chainCode, File localSiftsFile) throws SiftsMappingNotFoundException, IOException {		
		BufferedReader in = null;
		try {
			if(localSiftsFile == null) {
				URL url = new URL(SIFTS_URL);
				in = new BufferedReader(new InputStreamReader(url.openStream()));
			} else {
				in = new BufferedReader(new FileReader(localSiftsFile));
			}
			String line;
			while((line = in.readLine()) != null) {
				if(line.startsWith(pdbCode)) {
					//DEBUG: System.out.println(line);
					String[] fields = line.split("\t");
					//DEBUG: System.out.println(fields);
					String pdb = fields[0];
					String chain = fields[1];
					String uniprotId = fields[2];
					int cifBeg = Integer.parseInt(fields[4]);
					int cifEnd = Integer.parseInt(fields[5]);				
					int pdbBeg = Integer.parseInt(fields[6]);
					int pdbEnd = Integer.parseInt(fields[7]);				
					int uniBeg = Integer.parseInt(fields[8]);
					int uniEnd = Integer.parseInt(fields[9]);
					if(pdb.equals(pdbCode) && chain.equals(chainCode)) {
						// entry found
						this.pdbCode = pdb;
						this.chainCode = chain;
						this.uniprotId = uniprotId;
						this.cifBeg = cifBeg;
						this.cifEnd = cifEnd;
						this.pdbBeg = pdbBeg;
						this.pdbEnd = pdbEnd;
						this.uniBeg = uniBeg;
						this.uniEnd = uniEnd;
						in.close();
						return;
					}
				}
			}
			in.close();
		} catch (MalformedURLException e) {
			e.printStackTrace();
		} catch (NumberFormatException e) {
			throw new IOException(e.getMessage());
		} finally {
			if(in != null) in.close();
		}
		// entry not found;
		throw new SiftsMappingNotFoundException(pdbCode+chainCode + " not found in document at " + SIFTS_URL);
	}
	
	/**
	 * Queries the SIFTS mapping at the default FTP location for the given Pdb and Chain code and stores the results in the member variables.
	 * @throws IOException if the SIFTS source was not found or an error occured while reading from it
	 * @throws SiftsMappingNotFoundException if the SIFTS source was found but contained no mapping for the given pdb and chain code
	 */	
	public SiftsMapping(String pdbCode, String chainCode) throws SiftsMappingNotFoundException, IOException {
		this(pdbCode, chainCode, null);
	}
}
