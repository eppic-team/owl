package proteinstructure;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.TreeMap;
import tools.MySQLConnection;


public class Graph {

	public final static String MYSQLSERVER="white";
	public final static String MYSQLUSER=getUserName();
	public final static String MYSQLPWD="nieve";

	ArrayList<Contact> contacts;
	// nodes is a TreeMap of residue serials to residue types (3 letter code)
	TreeMap<Integer,String> nodes;
	String sequence;
	String accode;
	String chain;
	double cutoff;
	String ct;
	boolean directed=false;
	
	// these 2 fields only used when reading from db
	int graphid=0;
	int sm_id=0;
	
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
		if (ct.contains("/")){
			directed=true;
		}
	}
	
	/**
	 * Constructs Graph object from graph db
	 * ATTENTION!! chain is the internal database identifier, NOT! the pdb chain code
	 * TODO: we should also have a method to construct Graph from db using a pdb chain code 
	 * @param dbname
	 * @param accode
	 * @param chain
	 * @param cutoff
	 * @param ct
	 */
	public Graph(String dbname, String accode, String chain, double cutoff, String ct) throws GraphIdNotFoundError{
		this.cutoff=cutoff;
		this.accode=accode;
		this.chain=chain;
		this.ct=ct;
		//TODO graphs in db are never directed, so this doesn't really apply here. Must solve all this!
		if (ct.contains("/")){
			directed=true;
		}
		MySQLConnection conn = new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD,dbname);
		getgraphid(conn); // initialises graphid and sm_id
		read_graph_from_db(conn);
		conn.close();
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
		if (ct.contains("/")){
			directed=true;
		}
		read_contacts_from_file(contactsfile);
	}

	//TODO implement (from python) write_graph_to_db, do we really need it here??

	/** get user name from operating system (for use as database username) */
	private static String getUserName() {
		String user = null;
		user = System.getProperty("user.name");
		if(user == null) {
			System.err.println("Could not get user name from operating system. Exiting");
			System.exit(1);
		}
		return user;
	}

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
	
	/**
	 * Reads contacts and nodes from db.
	 * The db must be a graph db following our standard format, i.e. must have tables: 
	 * chain_graph, single_model_graph, single_model_node, single_model_edge
	 * We don't care here about the origin of the data (msdsd, pdbase, predicted) for the generation of the graph as long as it follows our data format
	 * We read both edges and nodes from single_model_edge and single_model_node.
	 * The sequence is taken from nodes, thus it won't have unobserved or non standard aas.
	 * @param conn
	 */
	public void read_graph_from_db(MySQLConnection conn){
		contacts = new ArrayList<Contact>();
		nodes = new TreeMap<Integer, String>();
		sequence = "";
		try {
			String sql="SELECT i_num,j_num FROM single_model_edge WHERE graph_id="+graphid+" AND j_num>i_num ORDER BY i_num,j_num ";
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			while (rsst.next()) {
				int i=rsst.getInt(1);
				int j=rsst.getInt(2);
				contacts.add(new Contact(i,j));
			}
			rsst.close();
			stmt.close();
			sql="SELECT num,res FROM single_model_node WHERE graph_id="+graphid+" ORDER BY num ";
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			while (rsst.next()){
				int num=rsst.getInt(1);
				String res=rsst.getString(2);
				nodes.put(num, AA.oneletter2threeletter(res));
				sequence+=res;
			}
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}

	}
	
	public void getgraphid (MySQLConnection conn) throws GraphIdNotFoundError{
        // NOTE: as chain we are using our internal identifier, which is the pchain_code in msdsd or the asym_id in pdbase
        // in the chain_graph table the internal chain identifier is called 'pchain_code'
		int pgraphid=0;
		try {
			String sql="SELECT graph_id FROM chain_graph WHERE accession_code='"+accode+"' AND pchain_code='"+chain+"' AND dist="+cutoff;
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			int check=0;
			while (rsst.next()) {
				check++;
				pgraphid=rsst.getInt(1);
			}
			if (check!=1){
				System.err.println("No pgraph_id match or more than 1 match for accession_code="+accode+", pchain_code="+chain+", dist="+cutoff);
			}
			rsst.close();
			stmt.close();
			// we set the ctstr to the same as ct except in ALL case, where it is BB+SC+BB/SC
			String ctstr=ct;
			if (ct.equals("ALL")){
				ctstr="BB+SC+BB/SC";
			}
			sql="SELECT graph_id,single_model_id FROM single_model_graph WHERE pgraph_id="+pgraphid+" AND CT='"+ctstr+"' AND dist="+cutoff+" AND CR='(true)' AND CW=1";
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			check=0;
			while (rsst.next()){
				check++;
				graphid=rsst.getInt(1);
				sm_id=rsst.getInt(2);
			}
			if (check!=1){
				System.err.println("No graph_id match or more than 1 match for pgraph_id="+pgraphid+", CT="+ctstr+" and cutoff="+cutoff);
				throw new GraphIdNotFoundError("No graph_id match or more than 1 match for pgraph_id="+pgraphid+", CT="+ctstr+" and cutoff="+cutoff);
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
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

