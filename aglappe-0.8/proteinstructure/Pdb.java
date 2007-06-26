package proteinstructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.ArrayList;
import java.util.Collections;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Pdb {
	
	public final static int DEFAULT_MODEL=1;
	
	HashMap<String,Integer> resser_atom2atomserial;
	HashMap<Integer,String> resser2restype;
	HashMap<Integer,Double[]> atomser2coord;
	HashMap<Integer,Integer> atomser2resser;
	HashMap<String,Integer> pdbresser2resser; // only used when reading from database
	HashMap<Integer,String> resser2pdbresser; // only used when reading from database
	
	HashMap<String,ArrayList<String>> aas2atoms = AA.getaas2atoms();
	public String sequence="";
	public String accode="";
    // given "external" pdb chain code, i.e. the classic pdb code ("NULL" if it is blank in original pdb file)	
	public String chaincode="";
	public int model=DEFAULT_MODEL;
	String db;
	
    // Our internal chain identifier (taken from dbs or pdb):
    // - in reading from pdbase or from msdsd it will be set to the internal chain id (asym_id field for pdbase, pchain_id for msdsd)
    // - in reading from pdb file it gets set to whatever parsed from the pdb file (i.e. can be also ' ')
	public String chain;
	
	/**
	 * Constructs Pdb object given accession code and pdb chain code. Model will be DEFAULT_MODEL
	 * Takes data from default pdbase-like database: PdbaseInfo.pdbaseDB
	 * @param accode
	 * @param chaincode
	 * @throws PdbaseInconsistencyError
	 * @throws PdbaseAcCodeNotFoundError
	 * @throws MsdsdAcCodeNotFoundError
	 * @throws MsdsdInconsistentResidueNumbersError
	 * @throws SQLException 
	 */
	public Pdb (String accode, String chaincode) throws PdbaseInconsistencyError, PdbaseAcCodeNotFoundError, MsdsdAcCodeNotFoundError, MsdsdInconsistentResidueNumbersError, SQLException {
		this(accode,chaincode,DEFAULT_MODEL,PdbaseInfo.pdbaseDB);		
	}

	/**
	 * Constructs Pdb object given accession code, pdb chain code and a model serial
	 * Takes data from default pdbase-like database: PdbaseInfo.pdbaseDB
	 * @param accode
	 * @param chaincode
	 * @param model_serial
	 * @throws PdbaseInconsistencyError
	 * @throws PdbaseAcCodeNotFoundError
	 * @throws MsdsdAcCodeNotFoundError
	 * @throws MsdsdInconsistentResidueNumbersError
	 * @throws SQLException 
	 */
	public Pdb (String accode, String chaincode, int model_serial) throws PdbaseInconsistencyError, PdbaseAcCodeNotFoundError, MsdsdAcCodeNotFoundError, MsdsdInconsistentResidueNumbersError, SQLException {
		this(accode,chaincode,model_serial,PdbaseInfo.pdbaseDB);		
	}
	
	/**
	 * Constructs Pdb object given accession code, pdb chain code and a source db for the data. Model will be DEFAULT_MODEL
	 * The db can be pdbase-like or msdsd-like
	 * @param accode
	 * @param chaincode
	 * @param db
	 * @throws PdbaseInconsistencyError
	 * @throws PdbaseAcCodeNotFoundError
	 * @throws MsdsdAcCodeNotFoundError
	 * @throws MsdsdInconsistentResidueNumbersError
	 * @throws SQLException 
	 */
	public Pdb (String accode, String chaincode, String db) throws PdbaseInconsistencyError, PdbaseAcCodeNotFoundError, MsdsdAcCodeNotFoundError, MsdsdInconsistentResidueNumbersError, SQLException {
		this(accode,chaincode,DEFAULT_MODEL,db);		
	}

	/**
	 * Constructs Pdb object given accession code, pdb chain code, model serial and a source db for the data.
	 * The db can be pdbase-like or msdsd-like 
	 * @param accode
	 * @param chaincode
	 * @param model_serial
	 * @param db
	 * @throws PdbaseInconsistencyError
	 * @throws PdbaseAcCodeNotFoundError
	 * @throws MsdsdAcCodeNotFoundError
	 * @throws MsdsdInconsistentResidueNumbersError
	 * @throws SQLException 
	 */
	public Pdb (String accode, String chaincode, int model_serial, String db) throws PdbaseInconsistencyError, PdbaseAcCodeNotFoundError, MsdsdAcCodeNotFoundError, MsdsdInconsistentResidueNumbersError, SQLException {
		this.accode=accode;
		this.chaincode=chaincode;
		this.model=model_serial;
		this.db=db;
		// we initialise it to chaincode, in read_pdb_data_from_pdbase gets reset to the right internal chain id. If reading msdsd it stays as chaincode
		this.chain=chaincode; 
		if (db.contains("msdsd")){
			read_pdb_data_from_msdsd(db);
		} else if (db.contains("pdbase")){
			read_pdb_data_from_pdbase(db);
		}
		// we initialise resser2pdbresser from the pdbresser2resser HashMap
		resser2pdbresser = new HashMap<Integer, String>();
		for (String pdbresser:pdbresser2resser.keySet()){
			resser2pdbresser.put(pdbresser2resser.get(pdbresser), pdbresser);
		}
	}
	
	/**
	 * Constructs Pdb object reading from pdb file given pdb chain code. Model will be DEFAULT_MODEL
	 * @param pdbfile
	 * @param chaincode
	 * @param dummy dummy parameter to distinguish the method from other with the same signature
	 * TODO dummy parameter is a dirty hack, must solve it in other way: subclassing
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public Pdb (String pdbfile, String chaincode, boolean dummy) throws FileNotFoundException, IOException {
		this(pdbfile,chaincode,DEFAULT_MODEL,dummy);
	}
	
	/**
	 * Constructs Pdb object reading from pdb file given pdb chain code and model serial
	 * @param pdbfile
	 * @param chaincode
	 * @param model_serial
	 * @param dummy dummy parameter to distinguish the method from other with the same signature
	 * TODO dummy parameter is a dirty hack, must solve it in other way: subclassing
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public Pdb (String pdbfile, String chaincode, int model_serial, boolean dummy) throws FileNotFoundException, IOException {
		this.model=model_serial;
		this.chaincode=chaincode;
		read_pdb_data_from_file(pdbfile);
	}

	public void read_pdb_data_from_pdbase(String db) throws PdbaseInconsistencyError, PdbaseAcCodeNotFoundError, SQLException{
		resser_atom2atomserial = new HashMap<String,Integer>();
		resser2restype = new HashMap<Integer,String>();
		atomser2coord = new HashMap<Integer,Double[]>();
		atomser2resser = new HashMap<Integer,Integer>();

		PdbaseInfo mypdbaseinfo = new PdbaseInfo(accode,chaincode,model,db);
		ArrayList<ArrayList> resultset = mypdbaseinfo.read_atomData();
		chain = mypdbaseinfo.chain;
		sequence = mypdbaseinfo.read_seq();
		pdbresser2resser = mypdbaseinfo.get_ressers_mapping();
		mypdbaseinfo.close();
		
		for (ArrayList result:resultset){
			int atomserial = (Integer) result.get(0);
			String atom = (String) result.get(1);
			String res_type = (String) result.get(2);
			int res_serial = (Integer) result.get(3);
			double x = (Double) result.get(4);
			double y = (Double) result.get(5);
			double z = (Double) result.get(6);
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
	
	public void read_pdb_data_from_msdsd(String db) throws MsdsdAcCodeNotFoundError, MsdsdInconsistentResidueNumbersError, SQLException {
		resser_atom2atomserial = new HashMap<String,Integer>();
		resser2restype = new HashMap<Integer,String>();
		atomser2coord = new HashMap<Integer,Double[]>();
		atomser2resser = new HashMap<Integer,Integer>();
		
		MsdsdInfo mymsdsdinfo = new MsdsdInfo(accode,chaincode,model,db);
		ArrayList<ArrayList> resultset = mymsdsdinfo.read_atomData();
		sequence = mymsdsdinfo.read_seq();
		chain = mymsdsdinfo.chain;
		pdbresser2resser = mymsdsdinfo.get_ressers_mapping();
		mymsdsdinfo.close();

		for (ArrayList result:resultset){
			int atomserial = (Integer) result.get(0);
			String atom = (String) result.get(1);
			String res_type = (String) result.get(2);
			int res_serial = (Integer) result.get(3);
			double x = (Double) result.get(4);
			double y = (Double) result.get(5);
			double z = (Double) result.get(6);
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
	
	/**
	 * To read the pdb data (atom coordinates, residue serials, atom serials) from file
	 * chain gets set to internal identifier: if input chain code NULL then chain will be ' '
	 * accode gets set to the one parsed in HEADER or to 'Unknown' if not found
	 * sequence gets set to the sequence read from ATOM lines (i.e. observed resdiues only)
	 * @param pdbfile
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public void read_pdb_data_from_file(String pdbfile) throws FileNotFoundException, IOException{
		resser_atom2atomserial = new HashMap<String,Integer>();
		resser2restype = new HashMap<Integer,String>();
		atomser2coord = new HashMap<Integer,Double[]>();
		atomser2resser = new HashMap<Integer,Integer>();
		boolean empty = true; // controls whether we don't find any atom line for given chaincode and model
		// we set chain to chaincode except for case NULL where we use " " (NULL is a blank chain code in pdb files)
		chain=chaincode;
		if (chaincode.equals("NULL")) chain=" ";
		int thismodel=DEFAULT_MODEL; // we initialise to DEFAULT_MODEL, in case file doesn't have MODEL lines 
		BufferedReader fpdb = new BufferedReader(new FileReader(new File(pdbfile)));
		String line;
		while ((line = fpdb.readLine() ) != null ) {
			Pattern p = Pattern.compile("^HEADER");
			Matcher m = p.matcher(line);
			if (m.find()){
				Pattern ph = Pattern.compile("^HEADER.{56}(\\d\\w{3})");
				Matcher mh = ph.matcher(line);
				if (mh.find()) {
					accode=mh.group(1).toLowerCase();
				} else {
					accode="Unknown";
				}
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
				//                                 serial    atom   res_type      chain res_ser     x     y     z
				Pattern pl = Pattern.compile(".{6}(.....).{2}(...).{1}(...).{1}"+chain+"(.{4}).{4}(.{8})(.{8})(.{8})",Pattern.CASE_INSENSITIVE);
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
		if (empty) System.err.println("Couldn't find any atom line for given chaincode: "+chaincode+", model: "+model);
		// now we read the sequence from the resser2restype HashMap
		// NOTE: we must make sure elsewhere that there are no unobserved residues, we can't check that here!
		ArrayList<Integer> ressers = new ArrayList<Integer>();
		for (int resser:resser2restype.keySet()) {
			ressers.add(resser);
		}
		Collections.sort(ressers);
		sequence="";
		for (int resser:ressers){
			String oneletter = AA.threeletter2oneletter(resser2restype.get(resser));
			sequence += oneletter;
		}
	}

	/**
	 * Dumps coordinate data into a file in pdb format (ATOM lines only)
	 * The chain dumped is the value of the chain field, i.e. our internal chain identifier for Pdb objects
	 * @param outfile
	 * @throws IOException
	 */
	public void dump2pdbfile(String outfile) throws IOException {
		TreeMap<Integer,Object[]> lines = new TreeMap<Integer,Object[]>();
		PrintStream Out = new PrintStream(new FileOutputStream(outfile));
		Out.println("HEADER  Dumped from "+db+". pdb accession code="+accode+", chain='"+chain+"'");
		for (String resser_atom:resser_atom2atomserial.keySet()){
			int atomserial = resser_atom2atomserial.get(resser_atom);
			int res_serial = Integer.parseInt(resser_atom.split("_")[0]);
			String atom = resser_atom.split("_")[1];
			String res_type = resser2restype.get(res_serial);
			Double[] coords = atomser2coord.get(atomserial);
			Object[] fields = {atomserial, atom, res_type, chain, res_serial, coords[0], coords[1], coords[2]};
			lines.put(atomserial, fields);
		}
		for (int atomserial:lines.keySet()){
			Out.printf("ATOM  %5d  %3s %3s %1s%4d    %8.3f%8.3f%8.3f\n",lines.get(atomserial));
		}
		Out.println("END");
		Out.close();
	}
	
	public void dumpseq(String seqfile) throws IOException {
		PrintStream Out = new PrintStream(new FileOutputStream(seqfile));
		Out.println(">"+accode+"_"+chaincode);
		Out.println(sequence);
		Out.close();
	}
	
	public int get_length(){
		return resser2restype.size();
	}
	
	public HashMap<Integer,Double[]> get_coords_for_ct(String ct) {
		HashMap<Integer,Double[]> coords = new HashMap<Integer,Double[]>(); 
		HashMap<String,String[]> restype2atoms = AA.ct2atoms(ct);
		for (int resser:resser2restype.keySet()){
			String[] atoms = restype2atoms.get(resser2restype.get(resser));
			for (String atom:atoms){
				if (resser_atom2atomserial.containsKey(resser+"_"+atom)){
					int atomser = resser_atom2atomserial.get(resser+"_"+atom);
					Double[] coord = atomser2coord.get(atomser);
					coords.put(atomser, coord);
				}
				else if (atom.equals("O") && resser_atom2atomserial.containsKey(resser+"_"+"OXT")){
					int atomser = resser_atom2atomserial.get(resser+"_"+"OXT");
					Double[] coord = atomser2coord.get(atomser);
					coords.put(atomser, coord);
				}
				else {
					System.err.println("Couldn't find "+atom+" atom for resser="+resser+". Continuing without that atom for this resser.");
				}
			}
		}
		return coords;
	}
	
	public TreeMap<Contact, Double> calculate_dist_matrix(String ct){
		TreeMap<Contact,Double> dist_matrix = new TreeMap<Contact,Double>();
		if (!ct.contains("/")){
			HashMap<Integer,Double[]> coords = get_coords_for_ct(ct);
			for (int i_atomser:coords.keySet()){
				for (int j_atomser:coords.keySet()){
					if (j_atomser>i_atomser) {
						Contact pair = new Contact(i_atomser,j_atomser);
						dist_matrix.put(pair, distance(coords.get(i_atomser),coords.get(j_atomser)));
					}
				}
			}
		} else {
			String i_ct = ct.split("/")[0];
			String j_ct = ct.split("/")[1];
			HashMap<Integer,Double[]> i_coords = get_coords_for_ct(i_ct);
			HashMap<Integer,Double[]> j_coords = get_coords_for_ct(j_ct);
			for (int i_atomser:i_coords.keySet()){
				for (int j_atomser:j_coords.keySet()){
					if (j_atomser!=i_atomser){
						Contact pair = new Contact(i_atomser,j_atomser);
						dist_matrix.put(pair, distance(i_coords.get(i_atomser),j_coords.get(j_atomser)));
					}
				}
			}
		}
		return dist_matrix;
	}
	
	/**
	 * Get the contacts for given contact type and cutoff for this Pdb object.
	 * Returns a Graph object with the contacts
	 * The graph is always directional, i.e. in crossed cases xx/yy: i corresponds to xx and j to yy
	 * ct here can be single contact type (e.g. Ca, BB) or crossed (e.g. BB/SC or Ca/Cb)
	 * @param ct
	 * @param cutoff
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public Graph get_graph(String ct, double cutoff){
		TreeMap<Contact,Double> dist_matrix = calculate_dist_matrix(ct);
		ContactList contacts = new ContactList();
        // we loop here over all indices of dist_matrix, 
        // we took care already that in symmetric cases (i.e. single contact type, not crossed) we have only one side of the matrix and 
        // in asymmetrical cases (i.e. crossed cases) we have both sides of the matrix
		for (Contact pair:dist_matrix.keySet()){
			int i_atomser=pair.i;
			int j_atomser=pair.j;
			if (dist_matrix.get(pair)<=cutoff){
				int i_resser = atomser2resser.get(i_atomser);
				int j_resser = atomser2resser.get(j_atomser);
				Contact resser_pair = new Contact(i_resser,j_resser);
                // for multi-atom models (BB, SC, ALL or BB/SC) we need to make sure that we don't have contacts from residue to itself or that we don't have duplicates				
				if (i_resser!=j_resser && (! contacts.contains(resser_pair))){
					contacts.add(resser_pair);
				}
			}
		}
		Collections.sort(contacts);
		TreeMap<Integer,String> nodes = new TreeMap<Integer,String>();
		for (int resser:resser2restype.keySet()){
			nodes.put(resser,resser2restype.get(resser));
		}
		Graph graph = new Graph (contacts,nodes,sequence,cutoff,ct,accode,chain,chaincode);
		return graph;
	}
	
	public Double distance(Double[] coords1, Double[] coords2){
		return Math.sqrt(Math.pow((coords1[0]-coords2[0]),2)+Math.pow((coords1[1]-coords2[1]),2)+Math.pow((coords1[2]-coords2[2]),2));
	}
	
	public int get_resser_from_pdbresser (String pdbresser){
		return pdbresser2resser.get(pdbresser);
	}
	
	public String get_pdbresser_from_resser (int resser){
		return resser2pdbresser.get(resser);
	}

	public int get_resser_from_atomser(int atomser){
		return atomser2resser.get(atomser);
	}
}
