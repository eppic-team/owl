package embed.contactmaps;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.*;

import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbasePdb;
import proteinstructure.RIGraph;

import tools.MySQLConnection;
import edu.uci.ics.jung.graph.util.Pair;
import embed.*;

public class ContactMap {
	
	private static RIGraph rig;
	
	private static SparseMatrix full;
	
	private SparseMatrix mat;
	
	private SparseMatrix square;
	
	private String pdb;
	
	private int native_counter;
	
	private int non_natives;
	
	public ContactMap (){};
	
	public ContactMap (String pdb_code) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		setContactMap(pdb_code);
	}
	
	public ContactMap (ContactMap cm){
		setContactMaps(cm);
	}
	
	
	public void setContactMap (String pdb_code) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		MySQLConnection conn = new MySQLConnection ();
		Pdb pdb = new PdbasePdb(pdb_code, "pdbase_20090728", conn);
		pdb.load("A");
		rig = pdb.getRIGraph("Ca", 9);
		int coverage  = (int) (((double) rig.getEdgeCount())*0.05);
		Bound[][] sparse = Individuals.randomSet(rig, conn, coverage);
		Bound[][] fullcm = Reconstructer.convertRIGraphToBoundsMatrix(rig);
		mat = convertToSparseMatrix(sparse);
		setPdbCode(pdb_code);
		if(full == null) setFullContactMap(convertToSparseMatrix(fullcm));
		setContactMap2();
		setNativeCounter();
	}
	
	public void setContactMaps (ContactMap cm){
		if(full == null) setFullContactMap(cm.getFullMap());
		setContactMap(cm);
		setPdbCode(cm.getPdbCode());
		setContactMap2();
		native_counter = cm.getNativeCounter();
		non_natives    = cm.getNonNatives();
	}
	
	public void setContactMap (HashMap<Pair<Integer>,Double> map, int dim){
		mat = new SparseMatrix (map,dim,dim);
	}
	
	public void setContactMap (ContactMap map){
		mat = new SparseMatrix(map.getMap());
	}
	
	public void setFullContactMap (SparseMatrix full_cm){
		full = new SparseMatrix (full_cm);
	}
	
	public void setContactMap2 (){
		if(mat != null){
			square = squareContactMap();
		}
	}
	
	public void setNativeCounter (){
		native_counter = countCommonContactsInSquareMap();
		non_natives    = square.getMatrix().size() - native_counter;
	}
	
	public void setPdbCode (String pdb_code){
		pdb = new String (pdb_code);
	}
	
	public SparseMatrix convertToSparseMatrix(Bound[][] bounds){
		if(bounds.length == bounds[0].length){
			int length = bounds.length;
			HashMap<Pair<Integer>,Double> map = new HashMap<Pair<Integer>,Double> ();
			for(int i = 0; i < length-1; i++){
				for(int j = i + 1; j < length; j++){
					if(bounds[i][j] != null && i != j - 1){
						Pair<Integer> pair1 = new Pair<Integer> (new Integer (i), new Integer (j));
						Pair<Integer> pair2 = new Pair<Integer> (new Integer (j), new Integer (i));
						map.put(pair1, new Double(1.0));
						map.put(pair2, new Double(1.0));
					}
				}
			}
			return new SparseMatrix (map,length,length);
		}
		else throw new IllegalArgumentException ("Bounds must always have matching dimensions!");
	}
	
	public int countCommonContactsInSquareMap (){
		Set<Pair<Integer>> key = getSquareMap().getIndexPairs();
		key.retainAll(getFullMap().getIndexPairs());
		return key.size();
	}
	
	public ContactMap merge (ContactMap map){
		if(getPdbCode().matches(map.getPdbCode())){
			int dim = map.getMap().getColumnDimension();
			HashMap<Pair<Integer>,Double> m = new HashMap<Pair<Integer>,Double>(); 
			Set<Pair<Integer>> key1  = getMap().getIndexPairs();
			Set<Pair<Integer>> key2  = map.getMap().getIndexPairs();
			Set<Pair<Integer>> key1a = getMap().getIndexPairs();
			Set<Pair<Integer>> key2a = map.getMap().getIndexPairs();
			int size = key1.size();
			key1.retainAll(key2);
			key1a.removeAll(key1);
			key2a.removeAll(key1);
			Iterator<Pair<Integer>> it = key1.iterator();
			while(it.hasNext()){
				Pair<Integer> pair = it.next();
				m.put(pair, new Double(1.0));
			}
			while(m.size() < size){
				Random rand = new Random ();
				Pair<Integer> pair1 = null;
				Pair<Integer> pair2 = null;
				int ind = rand.nextInt(2);
				if(ind == 0){
					Iterator<Pair<Integer>> it1 = key1a.iterator();
					if(it1.hasNext()){
						pair1 = it1.next();
						pair2 = new Pair<Integer> (pair1.getSecond(),pair1.getFirst());
					}
				}
				if(ind == 1){
					Iterator<Pair<Integer>> it1 = key2a.iterator();
					if(it1.hasNext()){
						pair1 = it1.next();
						pair2 = new Pair<Integer> (pair1.getSecond(),pair1.getFirst());
					}
				}
				if(m.put(pair1, new Double(1.0)) != null) key1a.remove(pair1); key2a.remove(pair1);  
				if(m.put(pair2, new Double(1.0)) != null) key1a.remove(pair2); key2a.remove(pair2);
			}
			ContactMap cm = new ContactMap();
			cm.setContactMap(m, dim);
			cm.setPdbCode(getPdbCode());
			if(full == null) cm.setFullContactMap(getFullMap());
			cm.setContactMap2();
			cm.setNativeCounter();
			return cm;
		}
		else throw new IllegalArgumentException ("Merging contact maps of different pdb codes does not make sense!");
	}
	
	public SparseMatrix getMap (){
		if(mat != null) return new SparseMatrix (mat);
		else throw new NullPointerException ("The contact map field must be initialized before calling this method.");
	}
	
	public SparseMatrix getFullMap (){
		if(full != null) return new SparseMatrix (full);
		else throw new NullPointerException ("The full contact map field was not initialized before calling this method.");
	}
	
	public String getPdbCode (){
		if(pdb != null) return new String (pdb);
		else throw new NullPointerException ("The pdb code field was not initialized before calling this method.");
	}
	
	public SparseMatrix getSquareMap (){
		return new SparseMatrix (square);
	}
	
	public int getNativeCounter (){
		return native_counter;
	}
	
	public int getNonNatives (){
		return non_natives;
	}
	
	public SparseMatrix squareContactMap (){
		if(getMap() != null){
			return (getMap().multiply(getMap()));
		}
		else throw new NullPointerException ("The contact map field must be initialized before calling this method.");
	}
	
	public Individuals convertToIndividuals () throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		if(mat != null && square != null && rig != null){
			HashSet<Pair<Integer>> pairs = convertToHashSet(mat.getIndexPairs());
			Individuals in = new Individuals ();
			in.setName(pdb);
			in.setChainCode(rig.getChainCode());
			Individuals.setFullContactMap(rig);
			in.setSequence(rig.getSequence());
			in.storer(pairs);
			in.setEntries(pairs);
			in.setNumOfContacts(pairs);
			in.setFullContact(rig.getEdgeCount());
			in.setErrorValues();
			return in;
		}
		else throw new NullPointerException ("Fields must be initialized before calling this method!");
	}
	
	public String toString (){
		String str = "";
		int[][] index_array = SortIntArray.converter(getMap().getIndexPairs());
		int length = index_array[0].length;
		for(int i = 0; i < length; i++){
			int fi = index_array[0][i], si = index_array[1][i];
			str += fi + "\t" + si + "\t" + getMap().getMatrixEntry(fi,si) + "\n"; 
		}
		return str;
	}
	
	public static HashSet<Pair<Integer>> convertToHashSet (Set<Pair<Integer>> oldset){
		HashSet<Pair<Integer>> newset = new HashSet<Pair<Integer>> (oldset.size());
		Iterator<Pair<Integer>> it = oldset.iterator();
		while(it.hasNext()){
			Pair<Integer> pair  = it.next();
			int f_in = pair.getFirst().intValue() + 1, s_in = pair.getSecond().intValue() + 1;
			Pair<Integer> npair = null;
			if(f_in < s_in) npair = new Pair<Integer> (new Integer(f_in), new Integer (s_in));
			else npair = new Pair<Integer> (new Integer(s_in), new Integer (f_in));
			newset.add(npair);
		}
		return newset;
	}
	
	public static void main (String[] args) throws SQLException, PdbCodeNotFoundError, PdbLoadError, FileNotFoundException, IOException{
		String dir      = "/project/StruPPi/gabriel/Arbeiten/";
		String addname1 = "cmEvolver/"; 
		ContactMap cm = new ContactMap ("1bkr");
		System.out.println(cm.getNativeCounter());
		Individuals in = cm.convertToIndividuals();
		System.out.println(in.toString());
		in.printToFile(dir, addname1, "test01"+in.getName());
	}

}
