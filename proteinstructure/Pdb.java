package proteinstructure;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Locale;
import java.util.TreeMap;
import java.util.ArrayList;
import java.util.Collections;

/**
 * A single chain pdb protein structure
 * 
 * @author		Jose Duarte
 * Class:		Pdb
 * Package:		proteinstructure
 */
public abstract class Pdb {
	
	protected final static int DEFAULT_MODEL=1;				// default model serial (NMR structures)
	public final static String NONSTANDARD_AA_LETTER="X";   // letter we assign to nonstandard aas to use in sequence
	
	protected HashMap<String,Integer> resser_atom2atomserial; // residue serial+atom name (separated by underscore) to atom serials
	protected HashMap<Integer,String> resser2restype;   // residue serial to 3 letter residue type 
	protected HashMap<Integer,Double[]> atomser2coord;  // atom serials to 3D coordinates
	protected HashMap<Integer,Integer> atomser2resser;  // atom serials to residue serials
	protected HashMap<String,Integer> pdbresser2resser; // pdb (author) residue serials (can include insetion codes so they are strings) to internal residue serials
	protected HashMap<Integer,String> resser2pdbresser; // internal residue serials to pdb (author) residue serials (can include insertion codes so they are strings)
	
	protected HashMap<Integer,String> resser2secstruct; // residue serials to secondary structure
	
	protected HashMap<String,ArrayList<String>> aas2atoms = AA.getaas2atoms(); // contains atom names for each aminoacid
	
	protected String sequence; 		// full sequence as it appears in SEQRES field
	protected String pdbCode;
    // given "external" pdb chain code, i.e. the classic (author's) pdb code ("NULL" if it is blank in original pdb file)	
	protected String pdbChainCode;
    // Our internal chain identifier:
    // - in reading from pdbase or from msdsd it will be set to the internal chain id (asym_id field for pdbase, pchain_id for msdsd)
    // - in reading from pdb file it coincides with pdbChainCode except for "NULL" where we use "A"
	protected String chainCode;
	protected int model;  			// the model serial for NMR structures
	
	protected String db;			// the db from which we have taken the data (if applies)
	
	
	

	/**
	 * Dumps coordinate data into a file in pdb format (ATOM lines only)
	 * The chain dumped is the value of the chainCode variable, i.e. our internal chain identifier for Pdb objects
	 * @param outfile
	 * @throws IOException
	 */
	public void dump2pdbfile(String outfile) throws IOException {
		TreeMap<Integer,Object[]> lines = new TreeMap<Integer,Object[]>();
		PrintStream Out = new PrintStream(new FileOutputStream(outfile));
		Out.println("HEADER  Dumped from "+db+". pdb code="+pdbCode+", chain='"+chainCode+"'");
		for (String resser_atom:resser_atom2atomserial.keySet()){
			int atomserial = resser_atom2atomserial.get(resser_atom);
			int res_serial = Integer.parseInt(resser_atom.split("_")[0]);
			String atom = resser_atom.split("_")[1];
			String res_type = resser2restype.get(res_serial);
			Double[] coords = atomser2coord.get(atomserial);
			Object[] fields = {atomserial, atom, res_type, chainCode, res_serial, coords[0], coords[1], coords[2]};
			lines.put(atomserial, fields);
		}
		for (int atomserial:lines.keySet()){
			// Local.US is necessary, otherwise java prints the doubles locale-dependant (i.e. with ',' for some locales)
			Out.printf(Locale.US,"ATOM  %5d  %3s %3s %1s%4d    %8.3f%8.3f%8.3f\n",lines.get(atomserial));
		}
		Out.println("END");
		Out.close();
	}
	
	public void dumpseq(String seqfile) throws IOException {
		PrintStream Out = new PrintStream(new FileOutputStream(seqfile));
		Out.println(">"+pdbCode+"_"+pdbChainCode);
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
		Graph graph = new Graph (contacts,nodes,sequence,cutoff,ct,pdbCode,chainCode,pdbChainCode);
		return graph;
	}
	
	private Double distance(Double[] coords1, Double[] coords2){
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
	
	public String getChainCode(){
		return this.chainCode;
	}
	
	public String getPdbChainCode(){
		return this.pdbChainCode;
	}
	
	public String getSecStructure(int resser){
		return this.resser2secstruct.get(resser);
	}
}
