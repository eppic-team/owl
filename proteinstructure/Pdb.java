package proteinstructure;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Locale;
import java.util.TreeMap;
import java.util.ArrayList;

import javax.vecmath.Point3d;
import javax.vecmath.Point3i;

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
	protected HashMap<Integer,Point3d> atomser2coord;  // atom serials to 3D coordinates
	protected HashMap<Integer,Integer> atomser2resser;  // atom serials to residue serials
	protected HashMap<String,Integer> pdbresser2resser; // pdb (author) residue serials (can include insetion codes so they are strings) to internal residue serials
	protected HashMap<Integer,String> resser2pdbresser; // internal residue serials to pdb (author) residue serials (can include insertion codes so they are strings)
	
	protected HashMap<Integer,String> resser2secstruct;   // residue serials to secondary structure
	protected TreeMap<String,Interval> secstruct2resinterval;// secondary structure element to residue serial intervals
	
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
			Point3d coords = atomser2coord.get(atomserial);
			Object[] fields = {atomserial, atom, res_type, chainCode, res_serial, coords.x, coords.y, coords.z};
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
	
	public TreeMap<Integer,Point3d> get_coords_for_ct(String ct) {
		TreeMap<Integer,Point3d> coords = new TreeMap<Integer,Point3d>(); 
		HashMap<String,String[]> restype2atoms = AA.ct2atoms(ct);
		for (int resser:resser2restype.keySet()){
			String[] atoms = restype2atoms.get(resser2restype.get(resser));
			for (String atom:atoms){
				if (resser_atom2atomserial.containsKey(resser+"_"+atom)){
					int atomser = resser_atom2atomserial.get(resser+"_"+atom);
					Point3d coord = atomser2coord.get(atomser);
					coords.put(atomser, coord);
				}
				else if (atom.equals("O") && resser_atom2atomserial.containsKey(resser+"_"+"OXT")){
					int atomser = resser_atom2atomserial.get(resser+"_"+"OXT");
					Point3d coord = atomser2coord.get(atomser);
					coords.put(atomser, coord);
				}
				else {
					System.err.println("Couldn't find "+atom+" atom for resser="+resser+". Continuing without that atom for this resser.");
				}
			}
		}
		return coords;
	}

	/**
	 * Returns the distance matrix as a HashMap with Contacts (residue serial pairs) as keys
	 * It doesn't make sense to call this method for multi atom contact 
	 * types (for each residue serial pair there's more than 1 distance)
	 * Thus before calling this we should check AA.isValidSingleAtomCT(ct)
	 * @param ct the contact type
	 * @return
	 */
	public HashMap<Edge, Double> calculate_dist_matrix(String ct){
		HashMap<Edge,Double> distMatrixAtoms = new HashMap<Edge,Double>();
		if (!ct.contains("/")){
			TreeMap<Integer,Point3d> coords = get_coords_for_ct(ct);
			for (int i_atomser:coords.keySet()){
				for (int j_atomser:coords.keySet()){
					if (j_atomser>i_atomser) {
						Edge pair = new Edge(i_atomser,j_atomser);
						distMatrixAtoms.put(pair, coords.get(i_atomser).distance(coords.get(j_atomser)));
					}
				}
			}
		} else {
			String i_ct = ct.split("/")[0];
			String j_ct = ct.split("/")[1];
			TreeMap<Integer,Point3d> i_coords = get_coords_for_ct(i_ct);
			TreeMap<Integer,Point3d> j_coords = get_coords_for_ct(j_ct);
			for (int i_atomser:i_coords.keySet()){
				for (int j_atomser:j_coords.keySet()){
					if (j_atomser!=i_atomser){
						Edge pair = new Edge(i_atomser,j_atomser);
						distMatrixAtoms.put(pair, i_coords.get(i_atomser).distance(j_coords.get(j_atomser)));
					}
				}
			}
		}

		HashMap<Edge,Double> distMatrixRes = new HashMap<Edge, Double>();
		for (Edge cont: distMatrixAtoms.keySet()){
			int i_resser = get_resser_from_atomser(cont.i);
			int j_resser = get_resser_from_atomser(cont.j);
			distMatrixRes.put(new Edge(i_resser,j_resser), distMatrixAtoms.get(cont));
		}

		return distMatrixRes;
	}
	
	/**
	 * Get the graph for given contact type and cutoff for this Pdb object.
	 * Returns a Graph object with the contacts
	 * We do geometric hashing for fast contact computation (without needing to calculate full distance matrix)
	 * @param ct
	 * @param cutoff
	 * @return
	 */
	public Graph get_graph(String ct, double cutoff){ 
		TreeMap<Integer,Point3d> i_coords = null;
		TreeMap<Integer,Point3d> j_coords = null;
		boolean directed = false;
		if (!ct.contains("/")){
			i_coords = get_coords_for_ct(ct);
			directed = false;
		} else {
			String i_ct = ct.split("/")[0];
			String j_ct = ct.split("/")[1];
			i_coords = get_coords_for_ct(i_ct);
			j_coords = get_coords_for_ct(j_ct);
			directed = true;
		}
		int[] i_atomserials = new  int[i_coords.size()]; // map from matrix indices to atomserials
		int[] j_atomserials = null;
		
		int SCALE=100; // i.e. we use units of hundredths of Amstrongs (thus cutoffs can be specified with a maximum precission of 0.01A)
		
		int boxSize = (int) Math.floor(cutoff*SCALE);
		
		HashMap<Point3i,Box> boxes = new HashMap<Point3i,Box>();
		int i=0;
		for (int i_atomser:i_coords.keySet()){
			//coordinates for atom serial atomser, we will use i as its identifier below
			Point3d coord = i_coords.get(i_atomser);
			int floorX = boxSize*((int)Math.floor(coord.x*SCALE/boxSize));
			int floorY = boxSize*((int)Math.floor(coord.y*SCALE/boxSize));
			int floorZ = boxSize*((int)Math.floor(coord.z*SCALE/boxSize));
			Point3i floor = new Point3i(floorX,floorY,floorZ);
			if (boxes.containsKey(floor)){
				// we put the coords for atom i in its corresponding box (identified by floor)
				boxes.get(floor).put_i_Point(i, coord);
				if (!directed){
					boxes.get(floor).put_j_Point(i, coord);
				}
			} else {
				Box box = new Box(floor);
				box.put_i_Point(i, coord);
				if (!directed){
					box.put_j_Point(i, coord);
				}
				boxes.put(floor,box);
			}
			i_atomserials[i]=i_atomser; //as atomserials in coords were ordered (TreeMap) the new indexing will still be ordered
			i++;
		}
		int j=0;
		if (directed) {
			j_atomserials = new  int[j_coords.size()];
			for (int j_atomser:j_coords.keySet()){
				//coordinates for atom serial atomser, we will use j as its identifier below
				Point3d coord = j_coords.get(j_atomser);
				int floorX = boxSize*((int)Math.floor(coord.x*SCALE/boxSize));
				int floorY = boxSize*((int)Math.floor(coord.y*SCALE/boxSize));
				int floorZ = boxSize*((int)Math.floor(coord.z*SCALE/boxSize));
				Point3i floor = new Point3i(floorX,floorY,floorZ);
				if (boxes.containsKey(floor)){
					// we put the coords for atom j in its corresponding box (identified by floor)
					boxes.get(floor).put_j_Point(j, coord);
				} else {
					Box box = new Box(floor);
					box.put_j_Point(j, coord);
					boxes.put(floor,box);
				}
				j_atomserials[j]=j_atomser; //as atomserials in coords were ordered (TreeMap) the new indexing will still be ordered
				j++;
			}
		} else {
			j_atomserials = i_atomserials;
		}

		
		double[][]distMatrix = new double[i_atomserials.length][j_atomserials.length];
		
		for (Point3i floor:boxes.keySet()){ // for each box
			// distances of points within this box
			boxes.get(floor).getDistancesWithinBox(distMatrix,directed);

			//TODO should iterate only through half of the neighbours here 
			// distances of points from this box to all neighbouring boxes: 26 iterations (26 neighbouring boxes)
			for (int x=floor.x-boxSize;x<=floor.x+boxSize;x+=boxSize){
				for (int y=floor.y-boxSize;y<=floor.y+boxSize;y+=boxSize){
					for (int z=floor.z-boxSize;z<=floor.z+boxSize;z+=boxSize){
						if (!((x==floor.x)&&(y==floor.y)&&(z==floor.z))) { // skip this box
							Point3i neighbor = new Point3i(x,y,z);
							if (boxes.containsKey(neighbor)){
								boxes.get(floor).getDistancesToNeighborBox(boxes.get(neighbor),distMatrix,directed);
							}
						}
					}
				}
			} 
		} 
		
		// getting the contacts (in residue serials) from the atom serials (partial) distance matrix 
		EdgeSet contacts = new EdgeSet();
		for (i=0;i<distMatrix.length;i++){
			for (j=0;j<distMatrix[i].length;j++){
				// the condition distMatrix[i][j]!=0.0 takes care of skipping several things: 
				// - diagonal of the matrix in case of undirected
				// - lower half of matrix in case of undirected
				// - cells for which we didn't calculate a distance as the 2 points were not in same or neighbouring boxes (i.e. too far apart)
				if (distMatrix[i][j]!=0.0 && distMatrix[i][j]<=cutoff){
					int i_resser = atomser2resser.get(i_atomserials[i]);
					int j_resser = atomser2resser.get(j_atomserials[j]);
					Edge resser_pair = new Edge(i_resser,j_resser);
					// for multi-atom models (BB, SC, ALL or BB/SC) we need to make sure that we don't have contacts from residue to itself or that we don't have duplicates				
					if (i_resser!=j_resser){ // duplicates are automatically taken care by the EdgeSet class wich is a TreeSet and doesn't allow duplicates 
						contacts.add(resser_pair);
					}
				}

			}
		}

		// creating and returning the graph object
		TreeMap<Integer,String> nodes = new TreeMap<Integer,String>();
		for (int resser:resser2restype.keySet()){
			nodes.put(resser,resser2restype.get(resser));
		}
		Graph graph = new Graph (contacts,nodes,sequence,cutoff,ct,pdbCode,chainCode,pdbChainCode);

		return graph;
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
	
	public TreeMap<String,Interval> getAllSecStructElements(){
		return secstruct2resinterval;
	}
}
