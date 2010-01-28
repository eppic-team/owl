package embed.contactmaps;

import proteinstructure.*;

import java.io.*;
import java.sql.SQLException;
import java.util.*;

//import edu.uci.ics.jung.*;
import edu.uci.ics.jung.graph.util.Pair;
import embed.Bound;
import embed.Distiller;
import embed.Reconstructer;
import embed.Scorer;

import Jama.Matrix;

/*import embed.BoundsSmoother.BoundsDigraphNode;
import embed.BoundsSmoother.SimpleEdge;*/
import tools.*;
//import org.jgap.Chromosome;

/**
 * This class is intended to deal with contact maps of a given protein. These contact maps contain information about the
 * sequence, pdb code and chain code. Additionally, it also contains a table of Integers, where the first value represents
 * the first amino acid <tt>i</tt> in the sequence and the second value represents the second amino acid <tt>j</tt>. Since the class
 * <code>{@link Scorer}</code> provides methods to compute error values to a given contact map, each instance of this class has
 * its own error, computed according to the <code>{@link Scorer#getDMError(RIGraph, double[][])</code> and <code>{@link Scorer#getCMError(Bound[][], RIGraph)</code> methods.
 * This class can also generate a random subset to a given number of contacts.
 * <p>
 * </p>
 * <p>
 * This class only considers <tt>C-alpha</tt> contacts, whose distance is equal to or less than 9.0 aangstroem. The default database
 * is <tt>pdbase_20090728</tt>. Additionally, this class contains methods to compare and breed individual contact maps, that is:
 * </p>
 * <p>
 * comparing the contact table, keep all common contacts and randomly filling up the missing contacts until the maximum number of
 * contacts is reached.
 * </p>
 * <p>
 * Methods, writing all important fields to a file and reading information from a contact map file, are also implemented. All file types,
 * accepted by this class, have either of mentioned extensions: '.indi', '.ust', '.cmap' or '.graph'. So prior to any call to such a method, one
 * must ensure, that such a file is present in the denoted directory. The standard output format is the '.indi' file type.   
 * </p>
 * <p>
 * This class can deal with multiple instances at the same time, if and only if all instances have the same pdb code. So in case of multiple
 * runs (like a loop), one has to explicitly call the method <code>{@link #clearFullCMandDM()}</code>, otherwise an Exception may occur.
 * </p>
 */
public class Individuals {//extends HashSet<Pair<Integer>> {
	
	/*----------------------------constants-----------------------------------------------------*/
	
	private static final long serialVersionUID = 1L;
	
	protected static final String path = "/project/StruPPi/jose/embed/essence/score_vs_rmsd";
	
	private static final String pdbaseDb = "pdbase_20090728";
	
	protected static final String[] file_extension = {"indi","ust","cmap","graph","cm"};
	
	/**
	 * contact type - default: C-alpha
	 */
	private static final String ct = "Ca";
	
	/**
	 * cutoff distance, set to default 9 angstroem
	 */
	private static final double di = 9.0;
	
	
	/*-----------------------------class variables----------------------------------------------*/
	
	/**
	 * the RIGraph with all contacts and information necessary to generate an instance of this class. It is important to notice,
	 * that in case of multiple runs with instances of this class having different pdb codes, one must call the method <code>{@link #clearFullCMandDM()}</code>
	 * before initializing any other instance of this class, in order to avoid Exception to be thrown.
	 */
	protected static RIGraph fullcontactmap;
	
	/**
	 * the full distance matrix, containing all contact pairs of the given protein below the threshold of 9 A. It is important to notice,
	 * that in case of multiple runs with instances of this class having different pdb codes, one must call the method <code>{@link #clearFullCMandDM()}</code>
	 * before initializing any other instance of this class, in order to avoid Exception to be thrown.
	 */
	protected static double[][] fulldistancematrix;
		
	/*-------------------------------------field members----------------------------------------*/

	/**
	 * field: array of integers representing contact pairs
	 */
	private int [][] entries;
	
	/**
	 * field
	 */
	protected boolean fifty;
	
	/**
	 * field: same as 'entries', in this case all contact pairs are stored in a HashSet, to avoid redundancy
	 */
	private HashSet<Pair<Integer>> store;
	
	/**
	 * field: number of contacts
	 */
	private int numOfContacts, fullContacts;
	
	/**
	 * String fields: 'name' = PDB code, 'sequence' = protein sequence and 'chainCode' = chain identifier
	 */
	private String name, sequence, chainCode;
	
	private double CMError, DMError;
	
	/*----------------------Constructors--------------------------------------------------*/
	/**
	 * zero parameter constructor, defines the default values of all non final fields
	 */
	public Individuals () {
		this.chainCode = " ";
		this.CMError = 0.0;
		this.DMError = 0.0;
		this.entries = new int[1][2];
		this.name = " ";
		this.sequence = " ";
		this.storer();
		this.numOfContacts = 0;
	}
	
	/**
	 * one parameter constructor, reads instance of Individuals as input parameter
	 * @param in
	 * @throws PdbLoadError 
	 * @throws PdbCodeNotFoundError 
	 * @throws SQLException 
	 */
	public Individuals (Individuals in) {
		this.setIndis(in);
		//setFullContactMap(in);
	}
	
	public Individuals (Individuals in, HashSet<Pair<Integer>> contactset) throws PdbCodeNotFoundError, SQLException, PdbLoadError{
		this.setIndis(in, contactset);
	}
		
	/**
	 * one parameter constructor, reads instance of RIGraph as input parameter
	 * @param rig
	 * @throws PdbLoadError 
	 * @throws PdbCodeNotFoundError 
	 * @throws SQLException 
	 * @throws Exception
	 */
	public Individuals (RIGraph rig) throws SQLException, PdbCodeNotFoundError, PdbLoadError {
		this.setIndis(rig);
	}
	
	/**
	 * Four parameter constructor, can either randomly sample a specified number of contacts or takes a specified number
	 * of contacts, the boolean parameter defines whether (<code> rand = true </code>) or not (false) random sampling is required.
	 * @param rig - an <code> RIGraph </code> instance that is converted to an <code> Individuals </code> instance
	 * @param conn - a connection to the database
	 * @param rand - parameter defining whether or not random sampling is wanted
	 * @param val - parameter specifying the number of contacts for the random sampling mode
	 * @throws ArrayIndexOutOfBoundsException
	 * @throws NullPointerException
	 * @throws PdbCodeNotFoundError
	 * @throws SQLException
	 * @throws PdbLoadError
	 */
	public Individuals (RIGraph rig, MySQLConnection conn, boolean rand, int val) throws ArrayIndexOutOfBoundsException, NullPointerException, PdbCodeNotFoundError, SQLException, PdbLoadError  {
		this.setIndis(rig, conn, rand, val);
		//conn.close();
	}
	
	/**
	 * four parameter constructor, can either randomly sample a specified number of contacts or takes a specified number
	 * of contacts, the boolean parameter defines whether (true) or not (false) random sampling is needed 
	 * @param rig - an RIGraph instance
	 * @param conn
	 * @param rand
	 * @param val
	 * @throws ArrayIndexOutOfBoundsException
	 * @throws NullPointerException
	 * @throws PdbCodeNotFoundError
	 * @throws SQLException
	 * @throws PdbLoadError
	 */
	public Individuals (RIGraph rig, MySQLConnection conn, boolean rand, double val) throws ArrayIndexOutOfBoundsException, NullPointerException, PdbCodeNotFoundError, SQLException, PdbLoadError  {
		this.setIndis(rig, conn, rand, val);
		//conn.close();
	}
	
	/**
	 * three parameter constructor, reads instances of Bound[][], String and a double array
	 * as input parameter
	 * @param bound
	 * @param i
	 * @param dm
	 */
	public Individuals (Bound[][] bound, String i, double[][] dm) {
		this.setIndis(bound, i, dm);
	}
	
	/**
	 * One-parameter constructor: reads all information needed to create an instance of <code>Individuals</code> from a
	 * contact map file. This constructor only accepts four types of contact maps with following file extensions: <tt>cmap</tt>, <tt>graph</tt>, <tt>indi</tt>
	 * or <tt>ust</tt>. If the denoted file does not exit or does not match the predefined extensions, this constructor exits with a referring message.
	 * 
	 * @param file
	 * @throws IOException
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public Individuals (String file) throws IOException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		File testfile = new File (file);
		if(testfile.exists()){
			this.setIndis(file);
		}
		else{
			System.err.println("The denoted file '"+file+"' does not exist.");
			System.exit(1);
		}
	}

	/*--------------------setters---------------------------------------------------------*/
	
	/**
	 * setter, copies all instance variables of input parameter 'in'.
	 */
	public void setIndis (Individuals in) {
		this.setName(in.getName());
		this.pdbCodeChecker(in.getName());
		//check, whether the two pdb codes match
		
		this.setIndis(in.getEntries());
		this.store = new HashSet<Pair<Integer>> (in.getHashSet());
		this.CMError = in.getCM();
		this.DMError = in.getDM();
		this.numOfContacts = in.getNumOfContacts();
		this.sequence = in.getSequence();
		this.chainCode = in.getChainCode();
		this.fullContacts = in.fullContacts;
	}
	
	/**
	 * setter, copies all instance variables of input parameter 'in', except for the field <code>{@link #store}</code>.
	 * Instead of the contact table of <tt>in</tt>, the parameter <tt>contactset</tt> is used as contact table.
	 * @param in
	 * @param contactset
	 * @throws PdbCodeNotFoundError
	 * @throws SQLException
	 * @throws PdbLoadError
	 */
	public void setIndis (Individuals in, HashSet<Pair<Integer>> contactset) throws PdbCodeNotFoundError, SQLException, PdbLoadError{
		this.setName(in.getName());
		this.pdbCodeChecker(in.getName());
		//check, whether the two pdb codes match
		
		this.setChainCode(in.getChainCode());
		this.setSequence(in.getSequence());
		this.store = new HashSet<Pair<Integer>> (contactset);
		this.setEntries(contactset);
		Bound[][] bound = Reconstructer.convertRIGraphToBoundsMatrix(this.getSubRIGraph());
		this.setCMError(bound, fullcontactmap);
		this.setDMError(bound,fulldistancematrix);
		this.setNumOfContact(contactset.size());
		this.setFullContact(fullcontactmap.getEdgeCount());
	}
	
	/**
	 * setter, reads an array of integers and uses them for the 'entries' field
	 * @param indices
	 */
	public void setIndis (int [][] indices) {
		int dim = indices.length;
		//int k = 0;
		//while(k < dim && (indices[k][0] != 0 || indices[k][1] != 0)){
			//k++;
		//}
		//TODO excluded since this code may cause an IllegalArgumentException, needs to be checked why this was put here
		
		this.entries = new int[dim][2];
		int counter1 = 0;
		int counter2 = 0;
		for(int i = 0; i < dim; i++){
			//for(int j = i + 1; j < dim; j++){
				if(indices[i][0] != indices[i][1] - 1){
					this.entries[counter1][0] = indices[i][0];
					this.entries[counter1][1] = indices[i][1];
					counter1++;
					counter2++;
				
			}
		}
	}
	
	/**
	 * Converts the <code> Bound[][] </code> instance to the field <code> entries </code> and checking
	 * whether or not the number of entries does not exceed the <code> numconts </code> 
	 * @param input - <code> Bound[][] </code> instance
	 * @param length - size of the <code> Bound[][] </code> array
	 * @param numconts - maximal number of entries
	 */
	public void setEntries(Bound[][] input, int length, int numconts){
		int counter = 0;
		this.entries = new int[numconts][2];
				//initializing the field entries, the first dimension containing the first contact partner
				//the second dimension containing the second contact partner
		
		for(int i = 0; i < length - 1; i++){						//loop for conversion
			
			boolean tester = false;									//if the number of contacts is exceeded
																	//this variable will be turned true
																	//and the loop stopped
			
			for(int j = i + 1; j < length; j++){
				if(input[i][j] != null && counter < numconts){		//check point to make sure no null contacts are
																	//taken and the number of contacts is not exceeded
					
					if(j - i > 1){									//check point to make sure no next neighbour contacts
																	//are taken
						
						this.entries[counter][0] = i;
						this.entries[counter][1] = j;
						counter++;
					}
				}
				else{												//if the number of contacts is exceeded a break option
																	//is executed
					if(counter >= numconts){
						tester = true;
						break;
					}
				}
			}
			if(tester){
				break;
			}
		}
	}
	
	/**
	 * Converts the <code> Bound[][] </code> instance to the field <code> entries </code> and checking
	 * whether or not the number of entries does not exceed the <code> numconts </code>.
	 * <p>
	 * The only difference to the method <code>{@link #setEntries(Bound[][], int, int)}</code> is that next neighbour contacts (i.e.
	 * are considered contacts.
	 * </p>
	 * @param length - size of the <code> Bound[][] </code> array
	 * @param input - <code> Bound[][] </code> instance
	 * @param numconts - maximal number of entries
	 */
	public void setEntries(int length, Bound[][] input, int numconts){
		int counter = 0;
		this.entries = new int[numconts][2];
		//initializing the field entries, the first dimension containing the first contact partner
		//the second dimension containing the second contact partner

		for(int i = 0; i < length - 1; i++){						//loop for conversion
			
			boolean tester = false;									//if the number of contacts is exceeded
																	//this variable will be turned true
																	//and the loop stopped
			
			for(int j = i + 1; j < length; j++){
				if(input[i][j] != null && counter < numconts){		//check point to make sure no null contacts are
																	//taken and the number of contacts is not exceeded
					//if(i != j - 1){
						this.entries[counter][0] = i;
						this.entries[counter][1] = j;
						counter++;
					//}
				}
				else{												//if the number of contacts is exceeded a break option
																	//is executed
					if(counter >= numconts){
						tester = true;
						break;
					}
				}
			}
			if(tester){
				break;
			}
		}
	}
	
	/**
	 * setter, converts a HashSet to an array of integers
	 * @param hash
	 */
	public void setEntries (HashSet<Pair<Integer>> hash) {
		int dim = hash.size();
		Iterator<Pair<Integer>> i = hash.iterator();
		this.entries = new int[dim][2];
		int j = 0;
		while(i.hasNext()) {
			Pair<Integer> pair = new Pair<Integer>(i.next());
			entries[j][0] = pair.getFirst() - 1;
			entries[j][1] = pair.getSecond() - 1;
			j++;
		}
	}
	
	/**
	 * setter, converts a RIGraph instance to a Individuals instance
	 * @param rig
	 * @throws PdbLoadError 
	 * @throws SQLException 
	 * @throws PdbCodeNotFoundError 
	 * @throws SQLException 
	 * @throws PdbCodeNotFoundError 
	 * @throws PdbLoadError 
	 */
	public void setIndis (RIGraph rig) throws PdbCodeNotFoundError, SQLException, PdbLoadError {
		this.areFullCMandDMinit(rig);

		this.setName(rig.getPdbCode());
		this.pdbCodeChecker(rig.getPdbCode());
		//PDB code is used as the 'name' of this instance
		
		Bound[][] bound = Reconstructer.convertRIGraphToBoundsMatrix(rig);
		int lengt = rig.getVertexCount();
		int numconts = rig.getEdgeCount();
		this.setEntries(lengt, bound, numconts);
		this.numOfContacts = rig.getEdgeCount();
		this.storer();
		this.setCMError(bound,fullcontactmap);
		this.setName(fullcontactmap.getPdbCode());
		this.setChainCode(fullcontactmap.getChainCode());
		this.setDMError(bound, fulldistancematrix);
		this.fullContacts = fullcontactmap.getEdgeCount();
		this.setSequence(rig.getSequence());
	}
	
	/**
	 * setter, does the same as setIndis (RIGraph rig) setter method, can do random sampling over a set of Individuals
	 * number of contacts specified by 'NumCont' parameter
	 * @param rig
	 * @param conn
	 * @param randomize
	 * @param NumCont
	 * @throws ArrayIndexOutOfBoundsException
	 * @throws NullPointerException
	 * @throws PdbCodeNotFoundError
	 * @throws SQLException
	 * @throws PdbLoadError
	 */
	public void setIndis (RIGraph rig, MySQLConnection conn, boolean randomize, int NumCont) throws ArrayIndexOutOfBoundsException, NullPointerException, PdbCodeNotFoundError, SQLException, PdbLoadError  {
		this.areFullCMandDMinit(rig);
		this.setName(rig.getPdbCode());
		this.pdbCodeChecker(rig.getPdbCode());
		//setFullContactMap(rig);
		Bound[][] bound = Reconstructer.convertRIGraphToBoundsMatrix(rig);
		int lengt = bound.length;
		int edgec = fullcontactmap.getEdgeCount();
		//int edgepercent = (int) ((double) NumCont/100.0*(double) edgec);
		if(randomize){
			bound = randomSet(rig,conn,NumCont);//edgepercent);
			this.setEntries(bound, lengt, NumCont);//edgepercent);
			//this.entries = new int[edgepercent][2];
		}
		else{
			this.setEntries(bound, lengt, rig.getEdgeCount());
			//this.entries = new int[rig.getEdgeCount()][2];
		}
		this.setCMError(bound, fullcontactmap);
		this.setDMError(bound, fulldistancematrix);
		/*for(int i = 0; i < lengt; i++){
			boolean tester = false;
			for(int j = i + 1; j < lengt; j++){
				if((bound[i][j] != null)&&(l < edgepercent)){
					if(i != j - 1){
						this.entries[l][0] = i;
						this.entries[l][1] = j;
						l++;
					}
				}
				else{
					tester = true;
					break;
				}
			}
			if(tester){
				break;
			}
		}*/
		this.fullContacts = edgec;
		this.storer();
		this.setNumOfContacts(this.getHashSet());
		this.setSequence(rig.getSequence());
		this.setChainCode(rig.getPdbChainCode());
		//conn.close();
	}
	
	/**
	 * setter, does the same as setIndis (RIGraph rig) setter method, can do random sampling over a set of Individuals
	 * number of contactes specified by 'NumCont' parameter
	 * @param rig
	 * @param conn
	 * @param randomize
	 * @param NumCont
	 * @throws ArrayIndexOutOfBoundsException
	 * @throws NullPointerException
	 * @throws PdbCodeNotFoundError
	 * @throws SQLException
	 * @throws PdbLoadError
	 */
	public void setIndis (RIGraph rig, MySQLConnection conn, boolean randomize, double NumCont) throws ArrayIndexOutOfBoundsException, NullPointerException, PdbCodeNotFoundError, SQLException, PdbLoadError  {
		this.areFullCMandDMinit(rig);
		this.setName(rig.getPdbCode());
		this.pdbCodeChecker(rig.getPdbCode());
		Bound[][] bound = Reconstructer.convertRIGraphToBoundsMatrix(rig);
		double cut = rig.getCutoff();
		String ct = rig.getContactType();
		int lengt = bound.length;
		Pdb n = getFullProt(rig, conn);
		RIGraph r = n.get_graph(ct, cut);
		int edgec = r.getEdgeCount();
		int edgepercent = (int) (NumCont/100.0*(double) edgec);
		double[][] dm = distMap(n);
		if(randomize){
			bound = randomSet(rig,conn,edgepercent);
			this.entries = new int[edgepercent][2];
		}
		else{
			this.entries = new int[rig.getEdgeCount()][2];
		}
		this.setCMError(bound, r);
		this.setDMError(bound, dm);
		this.setEntries(bound, lengt, edgepercent);
		/*for(int i = 0; i < lengt; i++){
			for(int j = i + 1; j < lengt; j++){
				if((bound[i][j] != null)&&(l < edgepercent)){
					if(i != j - 1){
						this.entries[l][0] = i;
						this.entries[l][1] = j;
						l++;
					}
				}
			}
		}*/
		this.fullContacts = edgec;
		this.storer();
		this.setNumOfContacts(this.getHashSet());
		this.setSequence(rig.getSequence());
		this.setChainCode(rig.getPdbChainCode());
		//conn.close();
	}
	
	/**
	 * setter, reads Bound instance, String instance an the distance map as a double array
	 * @param bound
	 * @param n
	 * @param dm
	 */
	public void setIndis (Bound[][] bound, String n, double[][] dm) {
		int lengt = bound.length, l = 0;
		this.entries = new int[lengt][2];
		for(int i = 0; i < lengt; i++){
			for(int j = i + 1; j < lengt; j++){
				if((bound[i][j] != null)&&(l < lengt)){
					if(i != j - 1){
						this.entries[l][0] = i;
						this.entries[l][1] = j;
						l++;
					}
				}
			}
		}
		this.setName(n);
		this.setDMError(bound, dm);
	}
	
	/**
	 * Setter: reads all important information from a contact map file or a graph file (files must have extensions specified by the constant
	 * String array <code>{@link #file_extension}</code>), initializing a <code>Individuals</code> instance. If either the PDB code or the contact table is missing
	 * this application will exit with an error message, specifying the cause of interruption. Make sure, the contact table is in the:
	 * <p>
	 * first contact as <code>(int)</code> \t second contact as <code>(int)</code>...
	 * </p>
	 * <p>
	 * If one of all the contact pairs exceeds the length of the amino acid sequence, the application will automatically exit, since this is
	 * a severe error in the contact map table part of the file.
	 * </p>
	 * <p>
	 * Note, that this method causes the thread to abruptly terminate, if either:
	 * </p>
	 * <p>
	 * - no file was found in the denoted directory,
	 * </p>
	 * <p>
	 * - no file in the predefined was found,
	 * </p>
	 * <p>
	 * - the file did not contain any pdb code (fatal error) or contact table (fatal error).
	 * </p>
	 * <p>
	 * In order to avoid this, one should make sure these requirements are met before calling this method.
	 * </p>
	 * @param file_path - String, denoting the absolute path and the file which will be read
	 * @throws IOException
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public void setIndis(String file_path) throws IOException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		File file = new File (file_path);
		if(file.exists()){
			//check, whether the file exists
			
			if(isReadableFileFormat(file)){
				//check, whether the string denoting file contains the correct extensions
				
				String filetype = getFileFormat(file);
				BufferedReader reader = new BufferedReader (new FileReader(file));
				String linereader = "";
				this.store = new HashSet<Pair<Integer>> ();
				while((linereader = reader.readLine()) != null){
					//loop, to initialize all fields of the Individuals class
					
					if(linereader.contains("#")){
						//the header starting with '#' contains important info as PDB code, chain code and sequence
						
						if(linereader.contains("PDB") && !linereader.contains("CHAIN")){
							this.name = linereader.split(": ")[1];
						}
						if(linereader.contains("PDB CHAIN")){
							this.chainCode = linereader.split(": ")[1];
						}
						if(linereader.contains("SEQUENCE")){
							this.sequence = linereader.split(": ")[1];
						}
						if(linereader.contains("CMError")){
							this.setCMError(linereader.split(": ")[1]);
						}
						if(linereader.contains("DMError")){
							this.setDMError(linereader.split(": ")[1]);
						}
						if(linereader.contains("NUMB. CONTACTS")){
							this.numOfContacts = (int) Double.parseDouble(linereader.split(": ")[1]);
						}
						if(linereader.contains("NUMB. OF ALL CONTACTS")){
							this.fullContacts = (int) Double.parseDouble(linereader.split(": ")[1]);
						}
					}
					else{
						if(linereader.contains("\t")){
							//the contact table is
							int fval = (int) Double.parseDouble(linereader.split("\t")[0]);
							int sval = (int) Double.parseDouble(linereader.split("\t")[1]);
							int seqlength = this.sequence.length();
							if(fval - 1 < seqlength && sval - 1 <= seqlength){
								//check, if each contact pair is within the range of the amino acid sequence length
								
								Integer first_val = new Integer (fval), sec_val = new Integer (sval);
								Pair<Integer> pair = new Pair<Integer> (first_val,sec_val);
								this.store.add(pair);
							}
							else{
								//exits, if at least one contact pair exceeds the sequence length, with the appropriate error message
								
								System.err.print("The pair '("+fval+","+sval+")' exceeded the length of this" +
										" protein "+this.name+", with "+seqlength+" amino acids!");
								System.exit(1);
							}
						}
					}
				}
				if(this.name == null){
					//exits, if the file did not contain any PDB code, causing a fatal error
					
					System.err.println("This '"+filetype+"' file in the directory '"+file_path+"' did not contain" +
							" any PDB code, causing an initialization problem.");
					System.exit(1);
				}
				if(this.store.size() == 0){
					//exits, if no contact was read from the file, causing a fatal error
					
					System.err.println("This '"+filetype+"' file in the directory '"+file_path+"' did not contain"+
							" any contacts as an array spaced by tabs.");
					System.exit(1);
				}
				if(this.sequence == null){
					//if no sequence was found
					
					MySQLConnection conn = new MySQLConnection ();
					Pdb pdb = new PdbasePdb (this.name, "pdbase_20090728", conn);
					this.setSequence(pdb.getSequence());
				}
				if(fullcontactmap == null && fulldistancematrix == null){
					MySQLConnection conn = new MySQLConnection ();
					Pdb pdb = new PdbasePdb(this.name, "pdbase_20090728", conn);
					pdb.load(this.chainCode);
					RIGraph fullCM = pdb.get_graph("Ca", 9);
					setFullContactMap(fullCM);
				}
				if(this.CMError == 0.0 || this.DMError == 0.0){					
					RIGraph subgraph = this.getSubRIGraph();
					Bound[][] sparse = Reconstructer.convertRIGraphToBoundsMatrix(subgraph);
					if(this.CMError == 0.0){
						this.setCMError(sparse,fullcontactmap);
					}
					if(this.DMError == 0.0){
						this.setDMError(sparse, fulldistancematrix);
					}
				}
				if(this.numOfContacts == 0){
					this.numOfContacts = this.store.size();
				}
				if(this.fullContacts == 0){
					this.fullContacts = fullcontactmap.getEdgeCount();
				}
				if(this.entries == null){
					this.setEntries(this.store);
				}
			}
			else{
				//exits, if the file is neither of the predefined file formats
				
				System.err.println("The denoted file '"+ file_path + "' is not in the accepted file format!");
				System.exit(1);
			}
		}
		else{
			//exits, if the file does not exist
			
			System.err.print("The file '"+file_path+"' does not exit!");
			System.exit(1);
		}
	}
	
	/**
	 * setter, sets the name of this instance using the String 'i'
	 * @param str
	 */
	public void setName (String str) {
		this.name = str;
	}
	
	/**
	 * setter, sets the instance variable CMError using Scorer.getCMError(Bound[][], RIGraph) method
	 * @param bound
	 * @param rig
	 * @throws PdbLoadError 
	 * @throws SQLException 
	 * @throws PdbCodeNotFoundError 
	 * @throws Exception
	 */
	public void setCMError (Bound[][] bound, RIGraph rig) throws PdbCodeNotFoundError, SQLException, PdbLoadError {
		this.CMError = getCMError(bound,rig);
	}
	
	/**
	 * a method, taking a String instance representing the 'CMError' and converts them to a double
	 * @param errorstring a String representing the 'CMError'
	 */
	public void setCMError (String errorstring){
		this.CMError = Double.parseDouble(errorstring);
	}
	
	/**
	 * setter, sets the instance variable DMError using Scorer.getDMError(Bound[][], double[][]) method
	 * @param bound
	 * @param dm
	 */
	public void setDMError (Bound[][] bound, double[][] dm) {
		this.DMError = getDMError(bound, dm);
	}
	
	/**
	 * a method, taking a String instance representing the 'CMError' and converts them to a double
	 * @param errorstring a String representing the 'DMError'
	 */
	public void setDMError (String errorstring){
		this.DMError = Double.parseDouble(errorstring);
	}
	
	/**
	 * a method converting a Individuals instance to RIGraph instance, uses this instance in
	 * order to compute the 'CMError' and 'DMError' 
	 * @throws PdbCodeNotFoundError
	 * @throws SQLException
	 * @throws PdbLoadError
	 */
	public void setErrorValues () throws PdbCodeNotFoundError, SQLException, PdbLoadError{
		RIGraph rig = reconstructGraph(this);
		Bound[][] bound = Reconstructer.convertRIGraphToBoundsMatrix(rig);
		this.setCMError(bound, fullcontactmap);
		this.setDMError(bound, fulldistancematrix);
	}
	
	/**
	 * setter, sets the instance variable numOfContacts using SparseGraph.getEdgeCount() method
	 * @param rig
	 */
	public void setNumOfContacts (RIGraph rig) {
		this.numOfContacts = rig.getEdgeCount();
	}
	
	/**
	 * sets the maximal number of contacts for this instance
	 * @param number the maximal number of contacts
	 */
	public void setNumOfContact (int number){
		this.numOfContacts = number;
	}
	
	/**
	 * setter, sets the instance variable numOfContacts using HashSet.size() method
	 * @param in
	 */
	public void setNumOfContacts (Individuals in) {
		this.numOfContacts = in.getHashSet().size();
	}

	/**
	 * setter, sets the instance variable numOfContacts using HashSet.size() method
	 * @param hash
	 */
	public void setNumOfContacts (HashSet<Pair<Integer>> hash) {
		this.numOfContacts = hash.size();
	}
	
	/**
	 * sets the field <code>{@link #fullContacts}</code>, which is the total number
	 * of contact for a given protein. Note, that this number is always constant for a given protein.
	 * @param contacts
	 */
	public void setFullContact (int contacts){
		this.fullContacts = contacts;
	}
	
	/**
	 * setter, stores all entries in 'entries' instance variable in a HashSet
	 */
	public void storer () {
		store = new HashSet<Pair<Integer>>();
		if(this.entries != null){
		int leng = (this.entries).length;
		//Integer[] inds = new Integer[leng];
		
		for(int i = 0; i < leng; i++){
			Integer in1 = new Integer(this.getNumbers(i, 0) + 1);
			Integer in2 = new Integer(this.getNumbers(i, 1) + 1);
				Pair<Integer> in = new Pair<Integer>(in1, in2);
				store.add(in);
		}
		store.hashCode();
		store.iterator();
		}
		else{
			throw new NullPointerException ("The field 'entries' must be initialized before calling this method!");
		}
	}
	
	/**
	 * takes a HashSet of Pairs of Integers and places them in the field <code>{@link #store}</code>.
	 * @param set
	 */
	public void storer (HashSet<Pair<Integer>> set){
		if(this.getSequence() != null){
			//if(set.size() == this.getSequence().length()){
				this.store = new HashSet<Pair<Integer>> (set);
				this.setEntries(set);
			/*}
			else{
				throw new IllegalArgumentException ("The number of entries in this HashSet exceeds the number of amino acids in this sequence.");
			}*/
		}
		else{
			this.store = new HashSet<Pair<Integer>> (set);
			this.setEntries(set);
		}
	}
	
	/**
	 * setter, sets this protein sequence
	 * @param seq
	 */
	public void setSequence(String seq) {
		this.sequence = seq;
	}
	
	/**
	 * setter, sets this protein chain code
	 * @param chain
	 */
	public void setChainCode (String chain) {
		this.chainCode = chain;
	}
	
	/**
	 * a sanity checker, any initialization of an instance of this class automatically calls this method, in 
	 * order to check, whether both Strings representing the pdb codes match. If they don't match an
	 * IllegalArgumentException is thrown.
	 * @param pdb_code a String representing the pdb code
	 * @throws IllegalArgumentException if both Strings dont match
	 */
	public void pdbCodeChecker (String pdb_code) throws IllegalArgumentException {
		if(fullcontactmap != null){
			if(!pdb_code.matches(this.getName())){
				throw new IllegalArgumentException ("Make sure, that only instances with matching pdb codes are treated at the same time.");
			}
		}
	}
	
	/**
	 * sanity checker: checks whether the constant <code>{@link #fullcontactmap}</code> null and in turn
	 * initializes this constant.
	 * @param rig an RIGraph instance
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public void areFullCMandDMinit (RIGraph rig) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		if(fullcontactmap == null){
			setFullContactMap(rig);
		}
	}
	
	/*-----------------------display methods-------------------------------------------*/
	
	/**
	 * display method, to display a HashSet
	 * @param hash
	 */
	public static void displayHash (HashSet<Pair<Integer>> hash) {
		System.out.println("Hash set: "+hash.toString());
	}
	
	/**
	 * display method, displays Individuals instances
	 */
	public void displayIndis (){
		System.out.print("{");
		int ind = (this.getEntries()).length - 1;
		System.out.print("# of contacts :"+this.getEntries().length+" and contacts of "+this.name+": ");
		for(int i = 0; i < ind; i++){
			System.out.print("{ "+this.getNumbers(i, 0)+", ");
			System.out.print(this.getNumbers(i, 1)+"},");//}
		}
		System.out.print("{ "+this.getNumbers(ind,0)+", "+this.getNumbers(ind,1)+"}");//}
		System.out.println("}, with CMError:"+this.getCM()+" and DMError: "+this.getDM());
	}
	
	/**
	 * Printing method: writes this 
	 * @throws IOException
	 * @throws FileNotFoundException
	 */
	public void printToFile() throws IOException, FileNotFoundException {
		FileOutputStream file = new FileOutputStream(path+name+"-"+this.getNumOfContacts()+file_extension);
		PrintStream printa = new PrintStream(file);
		printa.print(this.toString());
		printa.close();
		file.close();
		System.out.println(this.getName()+" written to file...");
		
	}
	
	/**
	 * additional printing method: writes an instance of this class to a specified file, using the <code>{@link #toString()}</code> method.
	 * The output exemplifies the standard output of this class. First, a header with all information needed by the program CMView, is written,
	 * followed by a table with all contact pairs in this instance. The default file extension is '.indi'. This is important since, this class
	 * only accepts only four types of contact map files.
	 * @param dir - a String denoting the absolute path of the output super directory
	 * @param addname1 - a sub directory of the output
	 * @param addname2 - a file name addition
	 * @throws IOException
	 * @throws FileNotFoundException
	 */
	public void printToFile(String dir, String addname1, String addname2) throws IOException, FileNotFoundException {
		File dirfile = new File(dir + addname1);
		if(!dirfile.exists()){
			dirfile.mkdirs();
		}
		FileOutputStream file = new FileOutputStream(dir + "/" + addname1 + "/" + this.name + addname2 + "-"+this.getNumOfContacts()+".indi");
		PrintStream printa = new PrintStream(file);
		printa.print(this.toString());
		printa.close();
		file.close();
		System.out.println(this.getName()+" written to file...");
	}
	
	/*----------------------getters----------------------------------------------*/
	
	/**
	 * getter, returns this chain code
	 * @return 'chainCode'
	 */
	public String getChainCode () {
		return new String (this.chainCode);
	}
	
	/**
	 * getter, returns this entries at position (i,j)
	 * @param i
	 * @param j
	 * @return 'entries[i][j]'
	 */
	public int getNumbers (int i, int j) {
		return this.entries[i][j];
	}

	/**
	 * getter, returns instance variable 'entries' as array of integers
	 * @return 'entries'
	 */
	public int[][] getEntries () {
		return this.entries;
	}
	
	/**
	 * getter, returns this name of the protein (using PDB code convention)
	 * @return 'name'
	 */
	public String getName () {
		return new String (this.name);
	}
	
	/**
	 * getter, returns this CMError
	 * @return 'CMError'
	 */
	public double getCM () {
		return this.CMError;
	}
	
	/**
	 * getter, returns this DMError
	 * @return 'DMError
	 */
	public double getDM () {
		return this.DMError;
	}
	
	/**
	 * auxiliary method: converts a distance matrix to a <code>{@link Bound}</code> matrix.
	 * @return
	 */
	public Bound[][] convertDistanceMatrixToBounds (){
		int protsize = this.sequence.length();
		Bound[][] bounds = new Bound[protsize][protsize];
		HashSet<Pair<Integer>> contactset = this.getHashSet();
		Iterator<Pair<Integer>> it = contactset.iterator();
		while(it.hasNext()){
			Pair<Integer> pair = it.next();
			int index1 = pair.getFirst().intValue()-1, index2 = pair.getSecond().intValue()-1;
			bounds[index1][index2] = new Bound(fulldistancematrix[index1][index2],fulldistancematrix[index1][index2]);
		}
		return bounds;
	}
	
	/**
	 * getter, returns this numOfContacts
	 * @return 'numOfContacts'
	 */
	public int getNumOfContacts () {
		return this.numOfContacts;
	}
	
	/**
	 * getter, returns this protein sequence
	 * @return 'sequence'
	 */
	public String getSequence () {
		return new String (this.sequence);
	}
	
	/**
	 * getter, returns this HashSet
	 * @return 'store'
	 */
	public HashSet<Pair<Integer>> getHashSet () {
		return new HashSet<Pair<Integer>> (this.store);
	}
		
	/**
	 * getter, returns boolean field 'fifty'
	 * @return 'fifty'
	 */
	public boolean getFifty() {
		return this.fifty;
	}
	
	/**
	 * compares this Individuals with in Individuals retaining only all contacts they have in common
	 * @param in
	 * @return 'hash1' - a Hashset where all entries are the same
	 */
	public HashSet<Pair<Integer>> comPareIndis(Individuals in) {
		HashSet<Pair<Integer>> hash1 = new HashSet<Pair<Integer>> (this.getHashSet());
		HashSet<Pair<Integer>> hash2 = new HashSet<Pair<Integer>> (in.getHashSet());
		hash1.retainAll(hash2);
		return new HashSet<Pair<Integer>> (hash1);
	}
	
	/**
	 * returns the field <code>{@link #fullContacts}</code>, the total number of contacts
	 * @return
	 */
	public int getFullContact() {
		return this.fullContacts;
	}
	
	/**
	 * Getter: converts this Individuals instance to a RIGraph instance, by copying the most
	 * relevant information to the new RIGraph instance, as
	 * <p>
	 * PDB code, sequence, chain code, cut off distance, contact type, the RIGNodes and the contact pairs
	 * </p>
	 * 
	 */
	public RIGraph getSubRIGraph (){
		RIGraph rig = new RIGraph ();
		//zero parameter constructor
		
		String seq = new String (this.sequence);
		rig.setSequence(seq);
		//setting the sequence field
		
		rig.setPdbCode(this.name);
		//setting the PDB code
		
		rig.setChainCode(this.chainCode);
		//setting the chain code
		
		rig.setContactType(ct);
		//setting the contact type
		
		rig.setCutoff(di);
		//setting the cut off distance
		for(int i = 0; i < seq.length(); i++){
			//loop to initialize all RIGNodes (i.e. the amino acids following the sequence)
			
			RIGNode node = new RIGNode (i,seq.substring(i,i+1));
			rig.addVertex(node);
			//adding each amino acid as a RIGNode
			
		}

		HashSet<Pair<Integer>> pairset = this.getHashSet();
		//copy of the contact map
		
		Iterator<Pair<Integer>> it = pairset.iterator();
		while(it.hasNext()){
			//loop over all contact pairs
			
			Pair<Integer> pair = it.next();
			int fval = pair.getFirst().intValue(), sval = pair.getSecond().intValue();
			rig.addEdgeIJ(fval,sval);
			//adding the edge of the i-th and j-th amino acids
		}
		return rig;
	}
	
	/**
	 * a method, converting all fields to a String representation. The standard output format is as follows:
	 * <p>
	 * first, header, in standard CMView Graph File format: protein sequence, pdb code, chain code, contact type (by default 'Ca'),
	 * cutoff distance (by default 9.0 aangstroem), both error values, the number of contacts and the total number of all contacts. 
	 * </p>
	 * <p>
	 * second, the actual contact table follows, which has the first value (corresponding to the lower value <tt>i</tt>, the second value
	 * <tt>j</tt>, by definition always greater than <tt>i</tt> and the frequency (for instances of this class, frequency is set to <tt>1.0</tt>.
	 */
	public String toString (){
		String str = "#CMVIEW GRAPH FILE ver: 1.0\n#SEQUENCE: "+this.getSequence()+"\n"+
		"#PDB: "+this.getName()+ "\n#PDB CHAIN CODE: "+this.getChainCode()+"\n#CT: "+getContactT()+ "\n#CUTOFF: "+getContactDist()+ "\n#CMError: "+ this.getCM() +  "\n#DMError: "+ this.getDM() + 
		 "\n#NUMB. CONTACTS: " + this.getNumOfContacts() +"\n" + "#NUMB. OF ALL CONTACTS: " + this.getFullContact() + "\n";
		int dim = 0;
		if(this.entries != null){
			dim = this.entries.length;
		}
		else{
			dim = this.store.size();
		}
		int[][] contact_array = SortIntArray.converter(this.store);
		for(int i = 0; i < dim; i++){
			str += + contact_array[0][i] + "\t" + contact_array[1][i] + "\t" + 1.0 + "\n";
		}
		return str;
	}
	
	/*------------------------------------------------statics---------------------------------------------------*/
	
	/**
	 * getter, returns contact type, by default: "CA" 
	 * @return 'ct' - this contact type
	 */
	public static String getContactT () {
		return ct;
	}
	
	/**
	 * getter, returns this cutoff distance
	 * @return 'di' - this cutoff distance
	 */
	public static double getContactDist () {
		return di;
	}
	
	/**
	 * getter, returns this CMError
	 * @param pop
	 * @return error
	 */
	public static double[] getCM (Individuals[] pop){
		int dim = pop.length;
		double[] error = new double[dim];
		for(int i = 0; i < dim; i++){
			error[i] = pop[i].getCM();
		}
		return error;
	}
	
	/**
	 * getter, returns this DMError
	 * @param pop
	 * @return error
	 */
	public static double[] getDM (Individuals[] pop){
		int dim = pop.length;
		double[] error = new double[dim];
		for(int i = 0; i < dim; i++){
			error[i] = pop[i].getDM();
		}
		return error;
	}
	
	/**
	 * compares two instances of Individuals by their HashSets retaining only those contacts both have in common,
	 * if contacts are missing the remaining gaps are filled randomly with non common contacts from both 'parents'.
	 * 
	 * @param in1
	 * @param in2
	 * @return new_indi
	 * @throws SQLException 
	 * @throws PdbCodeNotFoundError 
	 * @throws PdbLoadError
	 * 
	 */
	public static Individuals breedIndis (Individuals in1, Individuals in2) throws SQLException, PdbCodeNotFoundError, PdbLoadError, IllegalArgumentException {
		Individuals new_indi = new Individuals();
		//a new instance of Individuals class, with all fields set to null
		
		HashSet<Pair<Integer>> hash1 = new HashSet<Pair<Integer>>(in1.getHashSet());
		//contact set of the first parent
		
		HashSet<Pair<Integer>> hash2 = new HashSet<Pair<Integer>>(in2.getHashSet());
		//contact set of the second parent
		
		if(in1.getName().matches(in2.getName())){
			//check, whether both PDB codes match
			
			new_indi.setSequence(in1.getSequence());
			//initializing the field 'sequence'
			
			HashSet<Pair<Integer>> hash = in1.comPareIndis(in2);
			//set of common contacts
			
			new_indi.storer(hash);
			//common contacts kept for offspring
			
			new_indi.name = in1.getName();
			//PDB code
			
			new_indi.numOfContacts = hash.size();
			//number of common contacts
			
			new_indi.setEntries(hash);
			//initializing the field 'entries'
			
			new_indi.setChainCode(in1.getChainCode());
			//initializing the field 'chainCode'
			
			hash1 = removeSubSet(hash1,hash);
			hash2 = removeSubSet(hash2,hash);
			//removing all common contacts, so only non common contacts left in both HashSets
			
			Iterator<Pair<Integer>> i = hash1.iterator();
			Iterator<Pair<Integer>> j = hash2.iterator();
			while(new_indi.getNumOfContacts() != in1.getNumOfContacts()) {
				//loop to fill up the contact set to the number of contacts both parents have
				
				Random rand = new Random();
				//random number generator, to pick non common contacts at random
				
				if(rand.nextInt(2) % 2 == 0) {
					//if the random number is equivalent to 0 mod 2, pick contact from first HashSet (hash1)
					
					boolean test = new_indi.store.add(i.next());
					
					while(!test && i.hasNext()){
						test = new_indi.store.add(i.next());
						}
					i = hash1.iterator();
				}
				else {
					boolean test = new_indi.store.add(j.next()); 
					while(!test&&j.hasNext()){
						test = new_indi.store.add(j.next());
					}
					j = hash2.iterator();
				}
				new_indi.setNumOfContacts(new_indi.getHashSet());
			}
		new_indi.setEntries(new_indi.getHashSet());
		RIGraph rig = new RIGraph(new_indi.getSequence());
		rig = reconstructGraph(new_indi);		
		new_indi.setCMError(Reconstructer.convertRIGraphToBoundsMatrix(rig), fullcontactmap);
		new_indi.setDMError(Reconstructer.convertRIGraphToBoundsMatrix(rig), fulldistancematrix);
		new_indi.fullContacts = in1.fullContacts;
		}
		else{
			String name_mismatch = "PDB code of Individual 'in1' = " + in1.getName() + " does not match the PDB code of Individual 'in2' = " + in2.getName();
			throw new IllegalArgumentException("\n" + name_mismatch);
		}
		return new_indi;
	}
	
	/**
	 * reconstructer, reconstructs a RIGraph instance using an Individuals instance
	 * @param in an Individuals instance
	 * @return graph an RIGraph instance
	 */
	public static RIGraph reconstructGraph (Individuals in){
		String aa = in.getSequence();
		RIGraph graph = new RIGraph(aa);
		Iterator<Pair<Integer>> i = in.getHashSet().iterator();
		while(i.hasNext()) {
			Pair<Integer> pair = i.next();
			String aa1 = AAinfo.oneletter2threeletter(Character.toString(aa.charAt(pair.getFirst().intValue() - 1)));
			String aa2 = AAinfo.oneletter2threeletter(Character.toString(aa.charAt(pair.getSecond().intValue() - 1)));
			graph.addVertex(new RIGNode(pair.getFirst(), aa1));
			graph.addVertex(new RIGNode(pair.getSecond(), aa2));
			graph.addEdgeIJ(pair.getFirst(), pair.getSecond());
		}
		graph.setContactType(getContactT());
		graph.setCutoff(getContactDist());
		graph.setPdbChainCode(in.getChainCode());
		graph.setPdbCode(in.getName());
		graph.setChainCode(in.getChainCode());
		return graph;
	}
	
	/**
	 * setter: sets the constants <code>{@link #fullcontactmap}</code> and <code>{@link #fulldistancematrix}</code>, by getting all
	 * necessary information from a <code>{@link Pdb}</code> instance.
	 * @param in an Individuals instance
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public static final void setFullContactMap (Individuals in) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		String pdbCode = in.getName();

		//initializing String variable for chain code, only monomeric proteins excepted 
		String pdbChainCode = in.getChainCode();

		//initializing String variable for contact type
		String ct = Individuals.getContactT();

		//initializing double variable for cutoff length
		double cutoff = Individuals.getContactDist();
		
		MySQLConnection conn = new MySQLConnection();

		//initializing Pdb instance 'pdb' using subclass 'PdbasePdb' constructor 
		Pdb pdb = new PdbasePdb(pdbCode, pdbaseDb, conn);

		//loading PDB file of the same protein as specified by rigFile
		pdb.load(pdbChainCode);

		//initializing RIGraph instance with predefined parameters
		fullcontactmap = pdb.get_graph(ct, cutoff);
		
		Matrix fullDistanceMap = pdb.calculateDistMatrix(ct);
		
		fulldistancematrix =  fullDistanceMap.getArray();
		conn.close();
	}
	
	/**
	 * setter: sets the constants <code>{@link #fullcontactmap}</code> and <code>{@link #fulldistancematrix}</code>, by getting all
	 * necessary information from a <code>{@link Pdb}</code> instance.
	 * @param rig an RIGraph instance
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	protected static final void setFullContactMap (RIGraph rig) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		if(rig != null){
			String pdbCode = rig.getPdbCode();

			//initializing String variable for chain code, only monomeric proteins excepted 
			String pdbChainCode = rig.getChainCode();

			//initializing String variable for contact type
			String ct = Individuals.getContactT();

			//initializing double variable for cutoff length
			double cutoff = Individuals.getContactDist();

			MySQLConnection conn = new MySQLConnection();

			//initializing Pdb instance 'pdb' using subclass 'PdbasePdb' constructor 
			Pdb pdb = new PdbasePdb(pdbCode, pdbaseDb, conn);

			//loading PDB file of the same protein as specified by rigFile
			pdb.load(pdbChainCode);

			//initializing RIGraph instance with predefined parameters
			fullcontactmap = pdb.get_graph(ct, cutoff);

			Matrix fullDistanceMap = pdb.calculateDistMatrix(ct);

			fulldistancematrix =  fullDistanceMap.getArray();
			conn.close();
		}
		else{
			String rig_graph_error = "The residue interaction graph was not initialized.";
			throw new NullPointerException("\n" + rig_graph_error);
		}
	}	
	
	/*----------------------------statics-------------------------------------------*/
	
	/**
	 * Evaluates the error value for the given sparseBounds matrix (a subset of the contact map) 
	 * by first inferring bounds for all pairs through the triangle inequality and then 
	 * measuring how well the all pairs matrix fits the contact map.  
	 * Thus a high error value means the given sparseBounds matrix has low information content
	 * about the rest of the contact map and low error value that the matrix has high 
	 * information content.
	 * The error measure consists of the sum of the deviation of the upper bounds 
	 * to the ones in the contact map: sum(max(0,(u'i-ui))), i.e. if upper bound in 
	 * given matrix is below the one in the contact map there is no penalty.
	 * @param sparseBounds
	 * @param fullContactMap
	 * @return
	 * @throws SQLException 
	 * @throws PdbCodeNotFoundError 
	 * @throws PdbLoadError 
	 */
	public static double getCMError(Bound[][] sparseBounds, RIGraph fullContactMap) {
		Bound[][] cmapBounds = Reconstructer.convertRIGraphToBoundsMatrix(fullContactMap);
		// infer bounds for all pairs through triangle inequality
		Bound[][] bounds = Scorer.inferAllBounds(sparseBounds);
		double sumDev = 0;
		for (int i=0;i<bounds.length;i++) {
			for (int j=i+1;j<bounds.length;j++) {
				if (cmapBounds[i][j]!=null) {
					sumDev += Math.max(0,bounds[i][j].upper-cmapBounds[i][j].upper);
				}
			}
		}
		return sumDev/(double) fullContactMap.getFullLength(); 
	}
	
	/**
	 * Infer bounds for all pairs from a spare bounds matrix via the triangle inequality, evaluates 
	 * mean squared deviation for all given upper bounds per distance
	 * @param sub
	 * @param full
	 * @return
	 */
	public static double getDMError (Bound[][] sub, double[][] full) {
		Bound[][] sparse = Scorer.inferAllBounds(sub);
		//inferring all bounds via triangular inequality
		
		double error = 0;
		for(int index1 = 0; index1 < full.length; index1++){
			//looping over all bounds, first dimension
			
			for(int index2 = index1 +1; index2 < full.length; index2++){
				//looping over all bounds, second dimension
				
				if(full[index1][index2] != 0.0) {
					//only contacts of type (i,j), with j >= i + 1 are considered
					
					if(sparse[index1][index2].upper > full[index1][index2]){
					//only if the upper bound is greater than the exact distance value error is calculated
					
					error += Math.pow(sparse[index1][index2].upper - full[index1][index2],2.0 );
					//error function: square root of the sum of the square difference of upper bound and the exact distance,
					//divided by 1/2*n*(n - 1), where n is the number all contacts
					}
				}
			}
		}
		error = Math.pow(error, 0.5);
	return 2.0*error/((double) full.length*(full.length - 1));
	}
	/*-------------------------random generators------------------------------------*/
	
	/**
	 * random contact map generator, uses RIGraph MySQLConnection and a percentage input parameter
	 * @param input
	 * @param conn
	 * @param percent
	 * @return in
	 * @throws PdbCodeNotFoundError
	 * @throws NullPointerException
	 * @throws ArrayIndexOutOfBoundsException
	 * @throws SQLException
	 * @throws PdbLoadError
	 */
	public static Individuals randomSets (RIGraph input, MySQLConnection conn, double percent) throws NullPointerException, ArrayIndexOutOfBoundsException, PdbCodeNotFoundError, SQLException, PdbLoadError{
		Distiller dist = new Distiller(input);
		Bound[][] subset = dist.sampleSubset((int) (input.getEdgeCount()*percent));
		Individuals in = new Individuals(subset, input.getPdbCode(), distMap(getFullProt(input, conn)));
		//conn.close();
		return in;
	}
	
	/**
	 * method, that randomly selects contact of the full contact map
	 * @param input
	 * @param conn
	 * @return subset
	 * @throws Exception
	 * @throws SQLException
	 * @throws PdbLoadError
	 * @throws PdbCodeNotFoundError
	 * @throws NullPointerException
	 * @throws ArrayIndexOutOfBoundsException
	 */
	public static Bound[][] randomSet (RIGraph input, MySQLConnection conn) throws SQLException, PdbLoadError, PdbCodeNotFoundError, NullPointerException, ArrayIndexOutOfBoundsException {
		double cuto = input.getCutoff();
		String ct = input.getContactType();
		Pdb prot = getFullProt(input, conn);
		prot.load(input.getPdbChainCode());
		RIGraph full = prot.get_graph(ct, cuto);
		int numOfConts = input.getEdgeCount();
		Distiller dist = new Distiller(full);
		Bound[][] subset = dist.sampleSubset(numOfConts);
		//conn.close();
		return subset;
	}
	
	/**
	 * method that randomly selects contacts from the full contact map 
	 * @param input
	 * @param conn
	 * @param NumCont
	 * @return
	 * @throws SQLException
	 * @throws PdbLoadError
	 * @throws PdbCodeNotFoundError
	 * @throws NullPointerException
	 * @throws ArrayIndexOutOfBoundsException
	 */
	public static Bound[][] randomSet (RIGraph input, MySQLConnection conn, int NumCont) throws SQLException, PdbLoadError, PdbCodeNotFoundError, NullPointerException {
		double cuto = input.getCutoff();
		String ct = input.getContactType();
		Pdb prot = getFullProt(input, conn);
		prot.load(input.getPdbChainCode());
		RIGraph full = prot.get_graph(ct, cuto);
		int numOfConts = NumCont;
		Distiller dist = new Distiller(full);
		Bound[][] subset = dist.sampleSubset(numOfConts);
		//conn.close();
		return subset;
	}
	
	/**
	 * method, comparing the default file extensions of any output of this class with
	 * the absolute path of a File instance.
	 * @param file
	 * @return
	 */
	public static String getFileFormat (File file){
		String format = "";
		for(int i = 0; i < file_extension.length;i++){
			if(file.getName().contains(file_extension[i])){
				format = new String (file_extension[i]);
				break;
			}
		}
		return format;
	}
	/**
	 * method to recover a complete Pdb instance
	 * @param rig
	 * @param conn
	 * @return prot
	 * @throws SQLException 
	 * @throws PdbCodeNotFoundError 
	 * @throws PdbLoadError
	 */
	public static Pdb getFullProt (RIGraph rig, MySQLConnection conn) throws PdbCodeNotFoundError, SQLException, PdbLoadError {
		String pdbCode = rig.getPdbCode();
		String pdbaseDb = "pdbase_20090728";
		Pdb prot = new PdbasePdb(pdbCode, pdbaseDb, conn);
		prot.load(rig.getPdbChainCode());
		//conn.close();
		return prot;
	}
	
	/**
	 * method to calculate the distance map
	 * @param prot
	 * @return fullDistanceMap
	 */
	public static double[][] distMap (Pdb prot) {
		Matrix fullDistanceMap = prot.calculateDistMatrix("CA");
		return fullDistanceMap.getArray();
	}
	
	/**
	 * converts a Distance matrix to a <code>{@link Bound}</code> matrix.
	 * @param distmat
	 * @return
	 */			
	public static Bound[][] convertDistanceMatrixToBounds (double[][] distmat){
		int dim1 = distmat.length, dim2 = distmat[0].length;
		Bound[][] bounds = new Bound[dim1][dim2];
		for(int i = 0; i < dim1; i++){
			for(int j = 0; j < dim2; j++){
				if(distmat[i][j] <= 9.0){
					bounds[i][j] = new Bound(distmat[i][j],distmat[i][j]);
					bounds[j][i] = new Bound(distmat[i][j],distmat[i][j]);
				}
			}
		}
		return bounds;
	}
	
	/**
	 * checks, whether the denoted directory contains any standard output files of this class.
	 * If the parameter <tt>file</tt> is not a directory or it does not contain any files with
	 * predefined extensions, false will be returned 
	 * @param file a File instance representing a directory
	 * @return
	 */
	public static boolean containsReadableFiles(File file){
		if(file.exists() && file.isDirectory()){
			File[] list = file.listFiles();
			boolean test = false;
			for(int i = 0; i < list.length; i++){
				if(isReadableFileFormat(list[i])){
					test = true;
					break;
				}
			}
			return test;
		}
		else{
			return false;
		}
	}
	
	/**
	 * sanity checker, checks whether this File instance is compatible with the standard output format of
	 * this class.
	 * @param file a File instance representing a file
	 * @return true, if the File instance is compatible with the standard output of this class
	 */
	public static boolean isReadableFileFormat (File file){
		boolean test = false;
		for(int i = 0; i < file_extension.length;i++){
			if(file.getName().contains(file_extension[i])){
				test = true;
				break;
			}
		}
		return test;
	}
	
	/**
	 * 
	 * @param in
	 * @return
	 * @throws PdbCodeNotFoundError
	 * @throws SQLException
	 * @throws PdbLoadError
	 */
	public static HashSet<Individuals> localOptimizing (Individuals in) throws PdbCodeNotFoundError, SQLException, PdbLoadError{
		int[][] contacts = SortIntArray.converter(in.getHashSet());
		HashSet<Pair<Integer>> copy = in.getHashSet();
		int length = contacts.length, prot_length = in.getSequence().length();
		HashSet<Individuals> set = new HashSet<Individuals> ();
		double CM = in.getCM(), DM = in.getDM();
		for(int i = 0; i < length; i++){
			int f_val = contacts[0][i], s_val = contacts[1][i], counter = 0;
			Pair<Integer> pair = new Pair<Integer> (new Integer(f_val), new Integer (s_val));
			Individuals ini = new Individuals (in);
			ini.store.remove(pair);
			while(counter < prot_length - 1){
				int counterb = counter + 1;
				while(counterb < prot_length){
					Pair<Integer> newpair = new Pair<Integer> (new Integer(counter + 1), new Integer (counterb + 1));
					if(!copy.contains(newpair)){
						ini.store.add(newpair);
						ini.setErrorValues();
						double CM1 = ini.getCM(), DM1 = ini.getDM();
						if(CM1 < CM && DM1 < DM){
							set.add(ini);
						}
						else{
							ini.store.remove(newpair);
						}
					}
					counterb++;
				}
				counter++;
			}
		}
		return set;
	}
	
	/**
	 * converts a HashSet of Individuals instances to an array of Individuals instances
	 * @param set a HashSet of Individuals
	 * @return an array of Individuals
	 */
	public static Individuals[] convertToArray (HashSet<Individuals> set){
		Individuals[] array = new Individuals[set.size()];
		Iterator<Individuals> it = set.iterator();
		int counter = 0;
		while(it.hasNext()){
			array[counter] = new Individuals(it.next());
			counter++;
		}
		return array;
	}
	
	/**
	 * a method intended to locally optimize a given contact map, that is, remove one contact, add an other one and
	 * see, whether the error values decrease or not.
	 * @throws IOException
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public static void findLocalOptimum () throws IOException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		String dir = "/project/StruPPi/gabriel/Arbeiten/run_231009/run_two/1e0l/1e0l/deme0/temp2/1e0l0-4.indi";
		Individuals starter = new Individuals (dir);
		Individuals[] pop = convertToArray(localOptimizing(starter));
		Demes spec = new Demes (pop);
		spec.printAllIndisToTempDir(dir, "test");
	}
	
	/**
	 * method compares two HashSets and removes all common entries. If the HashSet <code>{@link #remove_set}</code> is not
	 * contained in <code>{@link #set1}</code>, an IllegalArgumentException is thrown. So the size of an HashSet received from
	 * 'retainAll(Collection<?>)' on both parameter must have the size of the remove_set parameter. 
	 * @param set1
	 * @param remove_set
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static HashSet<Pair<Integer>> removeSubSet (HashSet<Pair<Integer>> set1, HashSet<Pair<Integer>> remove_set) throws IllegalArgumentException {
		HashSet<Pair<Integer>> newset = new HashSet<Pair<Integer>> (set1);
		HashSet<Pair<Integer>> compset = comPareSets(set1,remove_set);
		Iterator<Pair<Integer>> it = remove_set.iterator();
		if(compset.size() == remove_set.size()){
			while(it.hasNext()){
				Pair<Integer> pair = it.next();
				newset.remove(pair);
			}
			return newset;
		}
		else{
			throw new IllegalArgumentException ("The HashSet 'remove_set' is not contained in the parameter 'set1'!");
		}
	}
	
	/**
	 * compares two HashSets of Pairs of Integers, retains all common entries and returns Set of all
	 * entries as a new HashSet.
	 * @param set1
	 * @param set2
	 * @return
	 */
	public static HashSet<Pair<Integer>> comPareSets(HashSet<Pair<Integer>> set1, HashSet<Pair<Integer>> set2){
		HashSet<Pair<Integer>> hash1 = new HashSet<Pair<Integer>> (set1);
		HashSet<Pair<Integer>> hash2 = new HashSet<Pair<Integer>> (set2);
		hash1.retainAll(hash2);
		return new HashSet<Pair<Integer>> (hash1);
	}
	
	/**
	 * an auxiliary method, that clears the constants: <code>{@link #fullcontactmap}</code> and <code>{@link #fulldistancematrix}</code>.
	 * After each run, with multiple instances of this class, one has to explicitly call this method, in order to clear the abovementioned
	 * constants, otherwise an Exception will occur.
	 */
	public static void clearFullCMandDM (){
		fullcontactmap = null;
		fulldistancematrix = null;
	}
	
	/*------------------------------main to test-----------------------------------------------*/
	/**
	 * Main: 
	 */
	public static void main (String[] args) throws Exception {
		String dir = "/project/StruPPi/gabriel/workspace_old/aglappe/embed/teststructures/";
		String[] prots = {"1bkr","1e0l","1e6k","1o8w","1odd","1onl","1pzc","1r9h","1sha","1ugm","2ci2"};
		HashMap<String,Integer> protmap = new HashMap<String,Integer> (11*2);
		File dirlist = new File (dir);
		File[] filelist = dirlist.listFiles(new RegexFileFilter (".*.graph"));
		Arrays.sort(filelist);
		int[] array = {2,3,7,8,9,10,11,12,13,14};
		int length = array.length;
		for(int i = 2; i < 3; i++){
			Integer protindex = new Integer (array[i]);
			protmap.put(prots[i],protindex);
		}
		//int length_arg = args.length;
		for(int i = 9; i < length; i++){
			//String arg = args[i];
			//int index = protmap.get(arg).intValue();
			int index = array[i];
			Individuals sathyasindis = new Individuals(filelist[index].getAbsolutePath());
			int numofconts = sathyasindis.getNumOfContacts();
			Species pop = new Species (index,numofconts,0.0003,true);
			Species.setPath("/project/StruPPi/gabriel/Arbeiten/"+sathyasindis.name+"/");
			//pop.evolve2(30);
			pop.copyParentalDemes();
			pop.evolve2(20);
			sathyasindis.printToFile(Species.path, "sathyas_set_"+sathyasindis.name+"/", "sathyas_set");
			pop.printStarterGenerationToFile(Species.path +  "/Starter/", "starter");
			clearFullCMandDM();
			
		}		
	}
		
	public static int compareTo (Pair<Double> tree, Pair<Double> input){
		double val1 = tree.getFirst().doubleValue();
		double val2 = input.getFirst().doubleValue();
		return (val1 < val2 ? -1 : (val1 == val2 ? 0 : 1));
	}
}