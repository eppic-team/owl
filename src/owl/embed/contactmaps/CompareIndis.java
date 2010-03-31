package owl.embed.contactmaps;

import java.sql.*;
import java.util.*;
import java.io.*;

import owl.core.structure.*;
import owl.core.util.MySQLConnection;
import owl.core.util.RegexFileFilter;

import edu.uci.ics.jung.graph.util.*;

public class CompareIndis extends Demes {
	
	//private static final String command1  = "select * from past_edges where generation=19 and stage=\"C\";";
	
	//private static final String trailer = "";
	
	private static final String pdb = "1e0l";
	
	private MySQLConnection conn;
	
	//private String path;
	
	private HashMap<Integer,HashSet<Pair<Integer>>> contact_sets;
	
	
	public CompareIndis (){
		super();
	}

	public CompareIndis (String dir) throws SQLException, PdbCodeNotFoundError, PdbLoadError, FileNotFoundException, IOException{
		super();
		conn = new MySQLConnection ();
		setContactSet(dir);
		generateIndis();
		conn.close();
	}
	
	public void setContactSet (String dir) throws SQLException, FileNotFoundException, IOException {
		File file = new File(dir);
		if(file.exists() && file.isDirectory()){
			File[] list = file.listFiles(new RegexFileFilter (".*.indi"));
			int length = list.length;
			if(length > 0){
				contact_sets = new HashMap<Integer,HashSet<Pair<Integer>>> ();
				for(int i = 0; i < length; i++){
					BufferedReader reader = new BufferedReader (new FileReader (list[i]));
					String linereader = null;
					HashSet<Pair<Integer>> set = new HashSet<Pair<Integer>> ();
					int counter = 0;
					while((linereader = reader.readLine()) != null){
						if(counter != 0){
							String[] str_ar = linereader.split("\t");
							Integer f_val = new Integer ((int) Double.parseDouble(str_ar[0]));
							Integer s_val = new Integer ((int) Double.parseDouble(str_ar[1]));
							Pair<Integer> pair = new Pair<Integer> (f_val,s_val);
							set.add(pair);
						}
						counter++;
					}
					Integer index = new Integer (i + 1);
					contact_sets.put(index, set);
				}
			}
		}
	}
	
	public void generateIndis () throws PdbCodeNotFoundError, SQLException, PdbLoadError{
		if(contact_sets !=null && contact_sets.size() > 0){
			HashMap<Integer,HashSet<Pair<Integer>>> copy = getContactSet();
			int size = copy.size();
			Individuals[] ar = new Individuals[size];
			Set<Integer> keys = copy.keySet();
			Iterator<Integer> it = keys.iterator();
			Pdb pr = new PdbasePdb(pdb,"pdbase_20090728",conn);
			pr.load("A");
			RIGraph rig = pr.getRIGraph("Ca", 9.0);
			Individuals.setFullContactMap(rig);
			while(it.hasNext()){
				Integer index = it.next();
				HashSet<Pair<Integer>> set = copy.get(index);
				int i = index.intValue() - 1;
				ar[i] = new Individuals();
				ar[i].storer(set);
				ar[i].setChainCode("A");
				ar[i].setName(pdb);
				ar[i].setSequence(rig.getSequence());
				ar[i].setNumOfContact(set.size());
				ar[i].setErrorValues();
				ar[i].setFullContact(rig.getEdgeCount());
			}
			super.setPop(ar);
			super.setCMErrorStats(ar);
			super.setDMErrorStats(ar);
			super.setGen(size);
			super.setSize(size);
			super.setWHasher();
			super.setMetrics();
		}
	}
	
	public HashMap<Integer,HashSet<Pair<Integer>>> getContactSet (){
		return new HashMap<Integer,HashSet<Pair<Integer>>> (contact_sets);
	}
	
	public String toString (){
		return new String (super.toString());
	}
	
	public static void main (String[] args) throws SQLException, PdbCodeNotFoundError, PdbLoadError, FileNotFoundException, IOException{
		String dir  = "/home/gmueller/shellscripts/1e0l/";
		String[] sub1 = {"Starter/","Final/"};
		for(int i = 0; i < 2; i++){
			CompareIndis ci = new CompareIndis (dir + sub1[i]);
			ci.printAllIndisToTempDir(dir + sub1[i], "");
		}
		
		//ci.printToFile(dir, "mich");
		
	}
}
