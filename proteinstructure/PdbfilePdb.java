package proteinstructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PdbfilePdb extends Pdb {
	
	private static final String UNKNOWN_STRING ="Unknown";
	private static final String NULL_chainCode = "A";
	
	private String pdbfile;

	/**
	 * Constructs Pdb object reading from pdb file given pdb chain code. Model will be DEFAULT_MODEL
	 * @param pdbfile
	 * @param pdbChainCode
	 * @throws FileNotFoundException
	 * @throws IOException
	 * @throws PdbfileFormatError
	 */
	public PdbfilePdb (String pdbfile, String pdbChainCode) throws FileNotFoundException, IOException, PdbfileFormatError {
		this(pdbfile,pdbChainCode,DEFAULT_MODEL);
	}
	
	/**
	 * Constructs Pdb object reading from pdb file given pdb chain code and model serial
	 * @param pdbfile
	 * @param pdbChainCode
	 * @param model_serial
	 * @throws FileNotFoundException
	 * @throws IOException
	 * @throws PdbfileFormatError
	 */
	public PdbfilePdb (String pdbfile, String pdbChainCode, int model_serial) throws FileNotFoundException, IOException, PdbfileFormatError {
		this.pdbfile = pdbfile;
		this.model=model_serial;
		this.pdbCode=UNKNOWN_STRING; // we initialise to unknown in case we don't find it in pdb file 
		this.pdbChainCode=pdbChainCode;
		// we set chainCode to pdbChainCode except for case NULL where we use "A"
		this.chainCode=pdbChainCode;
		if (pdbChainCode.equals("NULL")) this.chainCode=NULL_chainCode;
		
		this.sequence=""; // we initialise it to empty string, then is set inread_pdb_data_from_file 
		
		read_pdb_data_from_file();
		
		// when reading from pdb file we have no information of residue numbers or author's (original) pdb residue number, so we fill the mapping with the residue numbers we know
		//TODO eventually we could assign our own internal residue numbers when reading from pdb and thus this map would be used
		this.resser2pdbresser = new HashMap<Integer, String>();
		this.pdbresser2resser = new HashMap<String, Integer>();
		for (int resser:resser2restype.keySet()){
			resser2pdbresser.put(resser, String.valueOf(resser));
			pdbresser2resser.put(String.valueOf(resser), resser);
		}
	}

	
	
	/**
	 * To read the pdb data (atom coordinates, residue serials, atom serials) from file
	 * chainCode gets set to internal identifier: if input chain code NULL then chainCode will be ' '
	 * pdbCode gets set to the one parsed in HEADER or to 'Unknown' if not found
	 * sequence gets set to the sequence read from ATOM lines (i.e. observed resdiues only)
	 * @param pdbfile
	 * @throws FileNotFoundException
	 * @throws IOException
	 * @throws PdbfileFormatError
	 */
	private void read_pdb_data_from_file() throws FileNotFoundException, IOException, PdbfileFormatError{
		resser_atom2atomserial = new HashMap<String,Integer>();
		resser2restype = new HashMap<Integer,String>();
		atomser2coord = new HashMap<Integer,Double[]>();
		atomser2resser = new HashMap<Integer,Integer>();
		boolean empty = true; // controls whether we don't find any atom line for given pdbChainCode and model
		// we set chainCodeStr (for regex) to pdbChainCode except for case NULL where we use " " (NULL is a blank chain code in pdb files)
		String chainCodeStr=pdbChainCode;
		if (pdbChainCode.equals("NULL")) chainCodeStr=" ";
		
		int thismodel=DEFAULT_MODEL; // we initialise to DEFAULT_MODEL, in case file doesn't have MODEL lines 
		BufferedReader fpdb = new BufferedReader(new FileReader(new File(pdbfile)));
		int linecount=0;
		String line;
		while ((line = fpdb.readLine() ) != null ) {
			linecount++;
			Pattern p = Pattern.compile("^HEADER");
			Matcher m = p.matcher(line);
			if (m.find()){
				Pattern ph = Pattern.compile("^HEADER.{56}(\\d\\w{3})");
				Matcher mh = ph.matcher(line);
				if (mh.find()) {
					pdbCode=mh.group(1).toLowerCase();
				}
			} else if (linecount==1) { // header not found and we are in line 1
				// a HEADER is the minimum we ask at the moment for a pdb file to have, if we don't find it in line 1 we throw an exception
				throw new PdbfileFormatError("The pdb file "+pdbfile+" doesn't seem to have the right format");
			}
			p = Pattern.compile("^MODEL\\s+(\\d+)");
			m = p.matcher(line);
			if (m.find()){
				thismodel=Integer.parseInt(m.group(1));
			}
			if (thismodel!=model) continue; // we skip reading of atom lines if we are not in the desired model
			p = Pattern.compile("^ATOM");
			m = p.matcher(line);
			if (m.find()){
				//                                 serial    atom   res_type      chain 	   res_ser     x     y     z
				Pattern pl = Pattern.compile(".{6}(.....).{2}(...).{1}(...).{1}"+chainCodeStr+"(.{4}).{4}(.{8})(.{8})(.{8})",Pattern.CASE_INSENSITIVE);
				Matcher ml = pl.matcher(line);
				if (ml.find()) {
					empty=false;
					int atomserial=Integer.parseInt(ml.group(1).trim());
					String atom = ml.group(2).trim();
					String res_type = ml.group(3).trim();
					int res_serial = Integer.parseInt(ml.group(4).trim());
					double x = Double.parseDouble(ml.group(5).trim());
					double y = Double.parseDouble(ml.group(6).trim());
					double z = Double.parseDouble(ml.group(7).trim());
					Double[] coords = {x, y, z};
					ArrayList<String> aalist=AA.aas();
					if (aalist.contains(res_type)) {
						atomser2coord.put(atomserial, coords);
						atomser2resser.put(atomserial, res_serial);
						resser2restype.put(res_serial, res_type);
						ArrayList<String> atomlist = aas2atoms.get(res_type);
						if (atomlist.contains(atom)){
							resser_atom2atomserial.put(res_serial+"_"+atom, atomserial);
						}
					}					
				}
			}
		}
		fpdb.close();
		if (empty) System.err.println("Couldn't find any atom line for given pdbChainCode: "+pdbChainCode+", model: "+model);
		// now we read the sequence from the resser2restype HashMap
		// NOTE: we must make sure elsewhere that there are no unobserved residues, we can't check that here!
		ArrayList<Integer> ressers = new ArrayList<Integer>();
		for (int resser:resser2restype.keySet()) {
			ressers.add(resser);
		}
		Collections.sort(ressers);
		for (int resser:ressers){
			String oneletter = AA.threeletter2oneletter(resser2restype.get(resser));
			sequence += oneletter;
		}
	}

}
