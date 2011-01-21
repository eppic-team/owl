package owl.embed.contactmaps;

import java.io.*;
import java.sql.SQLException;
import java.util.*;

import owl.core.structure.Pdb;
import owl.core.structure.PdbCodeNotFoundException;
import owl.core.structure.PdbLoadError;
import owl.core.structure.PdbasePdb;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.MySQLConnection;
import owl.core.util.RegexFileFilter;
import owl.embed.*;


import edu.uci.ics.jung.graph.util.Pair;

public class ContactMap {
	
	private static SparseMatrix backbone;
	
	private static RIGraph rig;
	
	private static SparseMatrix full;
	
	private static SparseMatrix full2;
	
	private static SparseMatrix full3;
	
	private static int max_potence;
	
	private static double contact_coverage;
	
	private SparseMatrix mat;
	
	private SparseMatrix square;
	
	private HashMap<Integer,SparseMatrix> potence_series;
	
	private HashMap<Integer,SparseMatrix> potence_normalized;
	
	private String pdb;
	
	private int native_counter2;
	
	private int non_natives2;
	
	private int native_counter3;
	
	private int non_natives3;
	
	private double sensitivity2;
	
	private double sensitivity3;
	
	public ContactMap (){};
	
	public ContactMap (String pdb_code, int max_pot) throws SQLException, PdbCodeNotFoundException, PdbLoadError{
		contact_coverage = 0.1;
		max_potence = max_pot;
		setContactMap(pdb_code);
	}
	
	public ContactMap (String pdb_code, int max_pot, double contact) throws SQLException, PdbCodeNotFoundException, PdbLoadError, ContactMatrixException{
		if(contact > 0.0 && contact < 1.0) contact_coverage = contact;
		else throw new ContactMatrixException("Parameter 'contact' must be element in the open unit intervall, but was set to :"+contact+"!");
		max_potence = max_pot;
		setContactMap(pdb_code);
	}
	
	public ContactMap (ContactMap cm){
		setContactMaps(cm);
	}
	
	public ContactMap (int max_pot, String path) throws IOException, SQLException, PdbCodeNotFoundException, PdbLoadError{
		setContactMap(max_pot,path,0.1);
	}
	
	public ContactMap (int max_pot, String path, double contact) throws IOException, SQLException, PdbCodeNotFoundException, PdbLoadError{
		setContactMap(max_pot,path, contact);
	}

	public ContactMap (Individuals in, int max_pot){
		max_potence = max_pot;
		contact_coverage = 0.1;
		setContactMaps(in);
	}
	
	public ContactMap (Individuals in, int max_pot, double contact) throws ContactMatrixException {
		max_potence = max_pot;
		if(contact > 0.0 && contact < 1.0) contact_coverage = contact;
		else throw new ContactMatrixException("Parameter 'contact' must be element in the open unit intervall, but was set to :"+contact+"!");
		setContactMaps(in);
	}
	
	public void setContactMap (String pdb_code) throws SQLException, PdbCodeNotFoundException, PdbLoadError{
		MySQLConnection conn = new MySQLConnection ();
		Pdb pdb = new PdbasePdb(pdb_code, "pdbase_20090728", conn);
		pdb.load("A");
		rig = pdb.getRIGraph("Ca", 9);
		int coverage  = (int) (((double) rig.getEdgeCount())*contact_coverage);
		Bound[][] sparse = Individuals.randomSet(rig, conn, coverage);
		Bound[][] fullcm = Reconstructer.convertRIGraphToBoundsMatrix(rig);
		mat = convertToSparseMatrix(sparse);
		setBackbone();
		setPdbCode(pdb_code);
		if(full == null) setFullContactMap(convertToSparseMatrix(fullcm).add(backbone));
		setContactMap2();
		setPotenceSeries();
		setNativeCounter();
		setNativeCounter3();
		setSensitivity();
		conn.close();
	}
	
	public void setContactMaps (ContactMap cm){
		if(full == null) setFullContactMap(cm.getFullMap().add(backbone));
		setContactMap(cm);
		setBackbone();
		setPdbCode(cm.getPdbCode());
		setContactMap2(cm.getSquareMap());
		setPotenceSeries(cm.getPotenceSeries(),cm.getPotenceNormalized());
		native_counter2 = cm.getNativeCounter();
		non_natives2    = cm.getNonNatives();
		native_counter3 = cm.getNativeCounter3();
		non_natives3    = cm.getNonNatives3();
		setSensitivity();
	}
	
	public void setContactMap (int max_pot, String path, double contact) throws IOException, SQLException, PdbCodeNotFoundException, PdbLoadError, ContactMatrixException {
		if(contact > 0.0 && contact < 1.0) contact_coverage = contact;
		else throw new ContactMatrixException("Parameter 'contact' must be element in the open unit intervall, but was set to :"+contact+"!");
		max_potence = max_pot;
		Individuals in = new Individuals(path);
		setContactMaps(in);
	}
	
	public void setContactMaps (Individuals in){
		if(in != null){
			HashMap<Pair<Integer>,Double> map = new HashMap<Pair<Integer>,Double>();
			HashSet<Pair<Integer>> set = in.getHashSet();
			Iterator<Pair<Integer>> it = set.iterator();
			while(it.hasNext()){
				Pair<Integer> pair = it.next();
				Pair<Integer> npair1 = new Pair<Integer>(pair.getFirst().intValue()-1,pair.getSecond().intValue()-1);
				Pair<Integer> npair2 = new Pair<Integer>(pair.getSecond().intValue()-1,pair.getFirst().intValue()-1);
				map.put(npair1, new Double(1.0));
				map.put(npair2, new Double(1.0));
			}
			setContactMap(map,in.getSequence().length());
			setBackbone();
			if(rig == null) rig = Individuals.fullcontactmap;
			if(full == null) setFullContactMap(convertToSparseMatrix(Reconstructer.convertRIGraphToBoundsMatrix(rig)).add(backbone));
			setPdbCode(in.getPdbCode());
			setContactMap2();
			setPotenceSeries();
			setNativeCounter();
			setNativeCounter3();
			setSensitivity();
		}
		else throw new NullPointerException ("Parameter 'in' not initialized!");
	}
	
	public void setContactMap (HashMap<Pair<Integer>,Double> map, int dim){
		mat = new SparseMatrix (map,dim,dim);
	}
	
	public void setContactMap (ContactMap map){
		mat = map.getMap();
	}
	
	public void setFullContactMap (SparseMatrix full_cm){
		full  = new SparseMatrix (full_cm);
		full2 = full_cm.multiply(full_cm);
		full3 = full2.multiply(full_cm); 
	}
	
	public void setContactMap2 (){
		if(mat != null){
			square = squareContactMap();
		}
	}
	
	public void setContactMap2 (SparseMatrix sqr){
		square = sqr;
	}
	
	public void setNativeCounter (){
		native_counter2 = countCommonContactsInSquareMap();
		non_natives2    = square.getMatrix().size() - native_counter2;
	}
	
	public void setNativeCounter3 (){
		if(potence_normalized != null){
			Set<Pair<Integer>> ind = getPotenceSeries(3).getIndexPairs();
			Set<Pair<Integer>> ful = full.getIndexPairs();
			ful.retainAll(ind);
			native_counter3 = ful.size();
			non_natives3 = ind.size() - native_counter3;
		}
	}
	
	public void setPdbCode (String pdb_code){
		pdb = new String (pdb_code);
	}
	
	public void setPotenceSeries (){
		if(mat != null && square != null){
			potence_series = new HashMap<Integer,SparseMatrix>();
			potence_series.put(new Integer (1), mat);
			SparseMatrix id = new SparseMatrix (square);
			potence_series.put(new Integer(2), id);
			for(int i = 3; i < max_potence; i++){
				id = id.multiply(mat);
				potence_series.put(new Integer (i), id);
			}
			setNormalizedSeries();
		}
	}
	
	public void setPotenceSeries (HashMap<Integer,SparseMatrix> pot, HashMap<Integer,SparseMatrix> normal){
		potence_series     = pot;
		potence_normalized = normal;
	}
	
	@SuppressWarnings("deprecation")
	public void setNormalizedSeries (){
		if(potence_series != null){
			potence_normalized = new HashMap<Integer,SparseMatrix> ();
			for(int i = 1; i < max_potence; i++){
				Integer index = new Integer (i);
				SparseMatrix pot = potence_series.get(index);
				double greatest  = pot.getLargestEntrie();
				potence_normalized.put(index, pot.multiply(1.0/greatest));
			}
		}
	}
	
	public void setSensitivity (){
		if(potence_series != null){
			Set<Pair<Integer>> set2 = full2.getIndexPairs();
			Set<Pair<Integer>> set3 = full3.getIndexPairs();
			Set<Pair<Integer>> fl2  = getFullMap().getIndexPairs();
			Set<Pair<Integer>> fl3  = getFullMap().getIndexPairs();
			fl2.retainAll(set2);
			fl3.retainAll(set3);
			int size2 = fl2.size();
			int size3 = fl3.size();
			sensitivity2 = ((double) native_counter2)/((double) size2 + native_counter2);
			sensitivity3 = ((double) native_counter3)/((double) size3 + native_counter3);
		}
	}
	
	public void setBackbone (){
		if(backbone == null) backbone = addBackbone(mat.getColumnDimension());
		mat = mat.add(backbone);		
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
				else{
					Iterator<Pair<Integer>> it1 = key2a.iterator();
					if(it1.hasNext()){
						pair1 = it1.next();
						pair2 = new Pair<Integer> (pair1.getSecond(),pair1.getFirst());
					}
				}
				if(pair1 != null && pair2 != null){
					if(m.put(pair1, new Double(1.0)) != null) key1a.remove(pair1); key2a.remove(pair1);  
					if(m.put(pair2, new Double(1.0)) != null) key1a.remove(pair2); key2a.remove(pair2);
				}
			}
			ContactMap cm = new ContactMap();
			cm.setContactMap(m, dim);
			cm.setPdbCode(getPdbCode());
			if(full == null) cm.setFullContactMap(getFullMap().add(backbone));
			cm.setContactMap2();
			cm.setPotenceSeries();
			cm.setNativeCounter();
			cm.setNativeCounter3();
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
		return native_counter2;
	}
	
	public int getNonNatives (){
		return non_natives2;
	}
	
	public int getNativeCounter3 (){
		return native_counter3;
	}
	
	public int getNonNatives3 (){
		return non_natives3;
	}
	
	public double getSensitivity2(){
		return sensitivity2;
	}
	
	public double getSensitivity3(){
		return sensitivity3;
	}
	
	public HashMap<Integer,SparseMatrix> getPotenceSeries (){
		if(potence_series != null) return new HashMap<Integer,SparseMatrix>(potence_series);
		else throw new NullPointerException ("Potence series field not initialized!");
	}
	
	public SparseMatrix getPotenceSeries (int i) throws ContactMatrixException {
		Integer index = new Integer (i);
		if(getPotenceSeries().containsKey(index)){
			return getPotenceSeries().get(index);
		}
		else throw new ContactMatrixException ("No such entry!");
	}
	
	public HashMap<Integer,SparseMatrix> getPotenceNormalized (){
		if(potence_normalized != null) return new HashMap<Integer,SparseMatrix>(potence_normalized);
		else throw new NullPointerException ("Potence series field not initialized!");
	}
	
	public SparseMatrix getPotenceNormalized (int i) throws ContactMatrixException {
		Integer index = new Integer (i);
		if(getPotenceNormalized().containsKey(index)){
			return getPotenceNormalized().get(index);
		}
		else throw new ContactMatrixException ("No such entry!");
	}
	
	public SparseMatrix squareContactMap (){
		if(getMap() != null){
			return (getMap().multiply(getMap()));
		}
		else throw new NullPointerException ("The contact map field must be initialized before calling this method.");
	}
	
	public Individuals convertToIndividuals () throws SQLException, PdbCodeNotFoundException, PdbLoadError{
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
	
	public Individuals convertToIndividuals (SparseMatrix mat) throws SQLException, PdbCodeNotFoundException, PdbLoadError{
		if(mat != null&& rig != null){
			Set<Pair<Integer>> pairs = mat.getIndexPairs();
			HashSet<Pair<Integer>> hash = new HashSet<Pair<Integer>>();
			Iterator<Pair<Integer>> it = pairs.iterator();
			while(it.hasNext()){
				Pair<Integer> pair = it.next();
				int f_val = pair.getFirst().intValue(), s_val = pair.getSecond().intValue();
				if(f_val < s_val) hash.add(pair);
				else hash.add(new Pair<Integer>(new Integer(s_val),new Integer(f_val)));
			}
			Individuals in = new Individuals ();
			in.setName(pdb);
			in.setChainCode(rig.getChainCode());
			Individuals.setFullContactMap(rig);
			in.setSequence(rig.getSequence());
			in.storer(hash);
			in.setEntries(hash);
			in.setNumOfContacts(hash);
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
	
	public String toString (String dummy){
		return new String ("native2: "+native_counter2+"\tnative3: "+native_counter3+"\tsensitivity2: "+sensitivity2+"\tsensitivity3: "+sensitivity3);
	}
	
	public void writeToFile (String dir, String file_name) throws FileNotFoundException, IOException {
		if(potence_series != null && rig != null){
			File file = new File(dir);
			if(!file.exists()) file.mkdirs();
			HashMap<Integer,SparseMatrix> map = getPotenceSeries();
			for(int i = 1; i < max_potence; i++){
				String str = "#CMVIEW GRAPH FILE ver: 1.0\n#SEQUENCE: "+rig.getSequence()+"\n"+
				"#PDB: "+rig.getPdbCode()+ "\n#PDB CHAIN CODE: "+rig.getChainCode()+"\n#CT: "+rig.getContactType()+ "\n#CUTOFF: "+rig.getCutoff() + "\n";
				SparseMatrix mat = map.get(new Integer(i));
				Set<Pair<Integer>> indipair = mat.getIndexPairs();
				Iterator<Pair<Integer>> it = indipair.iterator();
				while(it.hasNext()){
					Pair<Integer> pair = it.next();
					int f_ind = pair.getFirst().intValue() + 1, s_ind = pair.getSecond().intValue() + 1;
					if(f_ind < s_ind){
						if(rig.containsEdgeIJ(f_ind, s_ind)) str += f_ind + "\t" + s_ind + "\t" + 1.0 + "\n";
						else str += f_ind + "\t" + s_ind + "\t" + 0.5 + "\n";
					}
				}
				FileOutputStream output = new FileOutputStream(dir+file_name+i+".cm");
				PrintStream printa      = new PrintStream (output);
				printa.print(str);
				printa.close();
				output.close();
				System.out.println("Potence series written to file...");
			}
		}
	}
	
	public void writeToFile2 (String dir, String file_name) throws FileNotFoundException, IOException {
		if(potence_normalized != null && rig != null){
			File file = new File(dir);
			if(!file.exists()) file.mkdirs();
			HashMap<Integer,SparseMatrix> map = getPotenceNormalized();
			for(int i = 1; i < max_potence; i++){
				String str = "#CMVIEW GRAPH FILE ver: 1.0\n#SEQUENCE: "+rig.getSequence()+"\n"+
				"#PDB: "+rig.getPdbCode()+ "\n#PDB CHAIN CODE: "+rig.getChainCode()+"\n#CT: "+rig.getContactType()+ "\n#CUTOFF: "+rig.getCutoff() + "\n";
				SparseMatrix mat = map.get(new Integer(i));
				Set<Pair<Integer>> indipair = mat.getIndexPairs();
				int[][] ind_array = SortIntArray.converter(indipair);
				int length = ind_array[0].length;
				for(int j = 0; j < length; j++){
					int f_ind = ind_array[0][j], s_ind = ind_array[1][j];
					if(f_ind <= s_ind) str += (f_ind+1) + "\t" + (s_ind+1) + "\t" + mat.getMatrixEntry(f_ind,s_ind) + "\n";
				}
				FileOutputStream output = new FileOutputStream(dir+file_name+i+".cm");
				PrintStream printa      = new PrintStream (output);
				printa.print(str);
				printa.close();
				output.close();
				System.out.println("Normalized potence series written to file...");
			}
		}
	}
	
	public static SparseMatrix addBackbone (int dim){
		HashMap<Pair<Integer>,Double> backbn = new HashMap<Pair<Integer>,Double>();
		int counter = 0;//, dim = mat.getColumnDimension();
		while(counter < dim-1){
			Integer f_val = new Integer (counter+1), s_val = new Integer (counter);
			Pair<Integer> pair1 = new Pair<Integer>(f_val,s_val), pair2 = new Pair<Integer>(s_val,f_val);
			backbn.put(pair1, new Double(1.0)); backbn.put(pair2, new Double(1.0));
			counter++;
		}
		return new SparseMatrix (backbn,dim,dim);		
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
	
	public static void main (String[] args) throws SQLException, PdbCodeNotFoundException, PdbLoadError, FileNotFoundException, IOException{
		/*String dir      = "/project/StruPPi/gabriel/Arbeiten/";
		String addname1 = "cmEvolver/"; 
		ContactMap cm = new ContactMap ("1bkr",11);
		System.out.println(cm.getNativeCounter());
		Individuals in = cm.convertToIndividuals();
		System.out.println(in.toString());
		in.printToFile(dir, addname1, "test01"+in.getName());
		cm.writeToFile(dir+addname1, "potence");
		cm.writeToFile2(dir+addname1, "potence_norm");*/
		String str = "/project/StruPPi/gabriel/Arbeiten/evolution/run_031209/1st_run/1bkr/Starter/";
		File file   = new File(str);
		File[] list = file.listFiles(new RegexFileFilter(".*.indi"));
		int length = list.length;
		Individuals[] in  = new Individuals[length];
		ContactMap[]  map = new ContactMap[length];
		for(int i = 0; i < length;i++){
			in[i] = new Individuals(list[i].getAbsolutePath());
			map[i] = new ContactMap(in[i],4);
			System.out.println("index="+i+"\t"+map[i].toString("dummy")+"\n");
		}
		String str1 = "/project/StruPPi/gabriel/Arbeiten/evolution/run_031209/1st_run/1bkr/";
		File dir = new File(str1);
		File[] list1 = dir.listFiles(new RegexFileFilter("deme.*"));
		int olength = list1.length;
		for(int i = 0; i < olength; i++){
			File hlp = new File(list1[i].getAbsolutePath()+"/temp2/");
			File[] list2 = hlp.listFiles(new RegexFileFilter(".*.indi"));
			int ilength = list2.length;
			Individuals[] in1 = new Individuals[ilength];
			ContactMap[] mat = new ContactMap[ilength];
			for(int j = 0; j < ilength; j++){
				in1[j] = new Individuals(list2[j].getAbsolutePath());
				mat[j] = new ContactMap(in1[j],4);
				System.out.println("index="+(ilength*i+j)+"\t"+mat[j].toString("dummy"));
			}
		}
	}
	
	class ContactMatrixException extends RuntimeException {
		private static final long serialVersionUID = 1L;
		public ContactMatrixException (String message){
			super(message);
		}
	}

}
