package embed.contactmaps;

import proteinstructure.*;
import tools.MySQLConnection;
import edu.uci.ics.jung.graph.util.*;
import embed.Bound;
import embed.Reconstructer;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;
import java.util.*;

/**
 * A class stripping contact maps, starting from a full contact map. All contacts are removed, ranked according to
 * the <code>{@link #getCMError(Bound[],RIGraph)}</code> and (!) <code>{@link #Scorer.getDMError(Bound[],RIGraph)}</code> defined in the <code>{@link #Scorer}</code> class.
 * If a contact removed yields a better CMError and DMError, the contact map (or <code>{@link #Individuals}</code> will be kept. The run continues until a final
 * number of contacts is reached. The contact maps ranked best will be written to a contact map file having the extension ".ust".
 * @author gmueller
 *
 */
public class Stripper {
	
	private static final String database = "pdbase_20090728";
	
	private int num_of_cycles;
	
	private Pdb protein;
	
	private RIGraph prot_graph;
	
	private HashMap<Integer,Individuals> unstrip_array;
	
	public Stripper (){
		protein = new Pdb();
		prot_graph = new RIGraph ();
		unstrip_array = new HashMap<Integer,Individuals> ();
	}
	
	public Stripper (String pdb_code) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		setUnStripper(pdb_code, 100);
	}
	
	public Stripper (String pdb_code, int num_of_runs) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		setUnStripper(pdb_code, num_of_runs);
	}
	
	/*-----------------------------------------void-------------------------------------------*/
	
	
	public void setUnStripper (String protein_PDB_code, int num_of_runs) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		MySQLConnection conn = new MySQLConnection ();
		protein = new PdbasePdb (protein_PDB_code, database, conn);
		protein.load("A");
		prot_graph = new RIGraph ();
		prot_graph = protein.getRIGraph("Ca", 9.0);
		unstrip_array = new HashMap<Integer,Individuals> (num_of_runs * 2);
		Individuals.setFullContactMap(prot_graph);
		Individuals ind = new Individuals (prot_graph);
		Integer index1 = new Integer (0);
		unstrip_array.put(index1,ind);
		num_of_cycles = num_of_runs;
		deselectContacts();
		conn.close();
	}
	
	public void deselectContacts () throws PdbCodeNotFoundError, SQLException, PdbLoadError{
		if(unstrip_array != null){
			int cycles = num_of_cycles;
			for(int i = 1; i <= cycles; i++){
				System.out.println("Unstripping: cycle "+i+"...");
				Individuals in = getIndi(i-1);
				HashSet<Pair<Integer>> contactset1 = in.getHashSet();
				HashSet<Pair<Integer>> contactset2 = new HashSet<Pair<Integer>>(contactset1);
				Iterator<Pair<Integer>> it = contactset1.iterator();
				while(it.hasNext()){
					Pair<Integer> pair = it.next();
					contactset2.remove(pair);
					Individuals newindi = new Individuals (in,contactset2);
					Integer first_index = new Integer(i);
					if(smallerError(newindi, first_index)){
						unstrip_array.put(first_index,newindi);
						//counter++;
					}
					contactset2 = new HashSet<Pair<Integer>> (contactset1);
				}
			}
		}
	}
	
	public void printToFile () throws IOException{
		if(this != null){
			String path = "/project/StruPPi/gabriel/Arbeiten/greedy/";
			File dir = new File (path);
			if(!dir.exists()) dir.mkdirs();
			File subdir = new File (dir.getAbsolutePath() + "/" + prot_graph.getPdbCode()+"/run_131009/");
			if(!subdir.exists()) subdir.mkdirs();
			int size = unstrip_array.size();
			for(int i = 0; i < size; i++){
				FileOutputStream out = new FileOutputStream (subdir.getAbsolutePath() +"/"+ i+".ust");
				PrintStream printa = new PrintStream (out);
				printa.print(toString(i));
				printa.close();
				out.close();
				System.out.println(i+"-th file written...");
			}
		}
	}
	
	/*---------------------------------------helper methods---------------------------------*/
	
	public HashMap<Integer,Individuals> getUnstripArray (){
		return new HashMap<Integer,Individuals> (unstrip_array);
	}
		
	public Individuals getIndi (int index){
		return new Individuals (unstrip_array.get(new Integer (index)));
	}
	
	public boolean haveSameContactMap (Integer index1, Integer index2){
		Individuals indi1 = getIndi(index1);
		Individuals indi2 = getIndi(index2);
		HashSet<Pair<Integer>> set1 = indi1.getHashSet(), set2 = indi2.getHashSet();
		set1.retainAll(set2);
		if(set1.size() < set2.size()){
			return false;
		}
		else{
			return true;
		}
	}
	
	public boolean smallerError (Individuals in, Integer index){
		if(unstrip_array.containsKey(index)){
			double CM = unstrip_array.get(index).getCM(), DM = unstrip_array.get(index).getDM();
			if(in.getCM() < CM && in.getDM() < DM){
				return true;
			}
			else{
				return false;
			}
		}
		else{
			return true;
		}
	}
	
	public String toString (){
		if(unstrip_array != null){
			String container = "";
			int size = unstrip_array.size();
			for(int i = 0; i < size; i++){
				container += getIndi(i).toString();
			}
			return container;
		}
		else{
			throw new NullPointerException ("The field 'unstrip_array' was not initialized before calling the 'toString()' method!");
		}
	}
	
	public String toString (int i){
		if(unstrip_array != null){
			String container = getIndi(i).toString();
			return container;
		}
		else{
			throw new NullPointerException ("The field 'unstrip_array' was not initialized before calling the 'toString()' method!");
		}
	}
	
	/*-----------------------------------------statics---------------------------------------*/
	
	
	public static Individuals reconstructIndi (Individuals in, HashSet<Pair<Integer>> set) throws PdbCodeNotFoundError, SQLException, PdbLoadError{
		Individuals newindi = new Individuals (in);
		newindi.storer(set);
		newindi.setEntries(set);
		newindi.setNumOfContact(in.getNumOfContacts());
		RIGraph rig = in.reconstructGraph ();
		Bound[][] bound = Reconstructer.convertRIGraphToBoundsMatrix(rig);
		newindi.setCMError(bound, Individuals.fullcontactmap);
		newindi.setDMError(bound,Individuals.fulldistancematrix);
		return newindi;
	}
	
	public static HashSet<Pair<Integer>> removeEntry (HashSet<Pair<Integer>> set, Pair<Integer> pair){
		HashSet<Pair<Integer>> copyset = new HashSet<Pair<Integer>> (set);
		if(copyset.contains(pair)){
			copyset.remove(pair);
			return new HashSet<Pair<Integer>> (copyset);
		}
		else{
			throw new IllegalArgumentException ("The pair ("+pair.getFirst()+","+pair.getSecond()+" is not contained in this set...");
		}
	}
	
	public static HashSet<Pair<Integer>> removeOneEntryAtRandom (HashSet<Pair<Integer>> set){
		Iterator<Pair<Integer>> it = set.iterator();
		Random rand = new Random ();
		int rand_val = rand.nextInt(set.size()), counter = 0;
		while(it.hasNext()){
			Pair<Integer> pair = it.next();
			if(counter == rand_val){
				set.remove(pair);
				break;
			}
			else{
				counter++;
			}
		}
		return new HashSet<Pair<Integer>> (set);
	}
	
	public static void main (String[] args) throws SQLException, PdbCodeNotFoundError, PdbLoadError, IOException{
		String[] prots = {"1bkr","1e0l","1e6k","1onl"};
		int[][] contacts = {{31,660},{4,171},{29,798},{39,816}};
		for(int i = 0; i < prots.length;i++){
			Stripper un = new Stripper (prots[i],contacts[i][1]-contacts[i][0]);
			un.printToFile();
		}
	}

}
