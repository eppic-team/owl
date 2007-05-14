package proteinstructure;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
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
	
	/**
	 * Constructs Graph object by passing ArrayList with contacts and TreeMap with nodes (res serials and types)
	 * Must also pass contact type, cutoff, accession code and chain
	 * @param contacts
	 * @param nodes
	 * @param sequence
	 * @param cutoff
	 * @param ct
	 * @param accode
	 * @param chain
	 */
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

	/**
	 * Constructs Graph object by reading a file with contacts
	 * An object created with this constructor will be missing the fields sequence and nodes
	 * That means it's not possible to get a ContactMap from it using getCM because CM needs both sequence and nodes
	 * @param contactsfile
	 * @param cutoff
	 * @param ct
	 * @throws IOException
	 * @throws FileNotFoundException
	 */
	public Graph (String contactsfile, double cutoff,String ct) throws IOException, FileNotFoundException{
		this.cutoff=cutoff;
		this.ct=ct;
		directed=false;
		if (ct.contains("/")){
			directed=true;
		}
		read_contacts_from_file(contactsfile);
	}

	//TODO implement (from python) constructors for reading from db
	//TODO implement (from python) read_contacts_from_db
	//TODO implement (from python) write_graph_to_db, do we really need this here??
	
	public void read_contacts_from_file (String contactsfile) throws FileNotFoundException, IOException {
		contacts = new ArrayList<Contact>();
		System.out.println("Reading contacts from file "+contactsfile);
		BufferedReader fcont = new BufferedReader(new FileReader(new File(contactsfile)));
		String line;
		while ((line = fcont.readLine() ) != null ) {
			int i = Integer.parseInt(line.split("\\s+")[0]);
			int j = Integer.parseInt(line.split("\\s+")[1]);
			contacts.add(new Contact(i,j));
		}
		fcont.close();
	}
	
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

