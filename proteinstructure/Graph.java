package proteinstructure;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;


public class Graph {

	ArrayList<Contact> contacts;
	TreeMap<Integer,String> nodes;
	String sequence;
	String accode;
	String chain;
	double cutoff;
	String ct;
	boolean directed;
	
	public Graph (ArrayList<Contact> contacts, TreeMap<Integer,String> nodes, String sequence, double cutoff,String ct, String accode, String chain) {
		this.contacts=contacts;
		this.cutoff=cutoff;
		this.nodes=nodes;
		this.sequence=sequence;
		this.accode=accode;
		this.chain=chain;
		this.ct=ct;
		directed=false;
		if (ct.contains("/")){
			directed=true;
		}
	}

	//TODO implement (from python) constructors for reading from file and db
	//TODO implement (from python) read_contacts_from_file
	//TODO implement (from python) read_contacts_from_db
	//TODO implement (from python) write_graph_to_db
	
	public void write_contacts_to_file (String outfile) throws IOException {
		PrintStream Out = new PrintStream(new FileOutputStream(outfile));
		for (Contact pair:contacts){
			int i_resser=pair.i;
			int j_resser=pair.j;
			Out.println(i_resser+"\t"+j_resser);
		}
		Out.close();		
	}
	
	public ContactMap getCM() {
		// residues is the map from residue nums to residue types used in ContactMap class, i.e. it is the same as Pdb.resser2restype or Graph.nodes
		TreeMap<Integer,String> residues = new TreeMap<Integer,String>();
		// we copy residues from nodes (deep copy) 
		for (int node:nodes.keySet()){
			residues.put(node, nodes.get(node));
		}
		// check if we are in directed or undirected case. If undirected we fill the opposite contacts to pass a full list of contacts to ContactMap (which contains full matrix)
		ArrayList<Contact> contacts2pass = new ArrayList<Contact>();
		if (directed){
			contacts2pass=contacts;
		} else {
			for (Contact cont:contacts){
				int i_resser = cont.i;
				int j_resser = cont.j;
				contacts2pass.add(new Contact(i_resser,j_resser));
				contacts2pass.add(new Contact(j_resser,i_resser));
			}
		}
		// construct the ContactMap object and return it
		ContactMap cm = new ContactMap(contacts2pass,residues,sequence);
		return cm;
		
	}
	
}

