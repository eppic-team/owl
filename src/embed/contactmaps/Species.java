package embed.contactmaps;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import tools.RegexFileFilter;
import edu.uci.ics.jung.graph.util.*;

/**
 * This class allows stepwise multiple runs of evolution by creating so called demes. Each deme
 * consists of a matrix of <code> {@link #Species} </code> instances. The first entry in the array 
 * corresponds to the first or parental generation, having a fixed number of contacts
 * (default: 5%). The following entries correspond to the first, second ... filial
 * generation. All constructors of this class only use the four parameter constructor
 * of the <code> {@link #Species(int,int,int,int)} </code>.
 * class.
 * @author gmueller
 *
 */
public class Species {
	
	/*-----------------------------------fields-----------------------------------------*/
	
	/**
	 * field - deme, a matrix of {@link Demes} instances
	 */
	public Demes[][] deme;
	
	/**
	 * int fields: demenumb - number of demes; evosteps - number of evolution steps;
	 * gen - number of generation
	 */
	public int demenumb, evosteps, gen, pop_size, numofConts;
	
	/**
	 * field: a threshold value used to discriminate fast convergence
	 */
	public static double threshold;
	
	/**
	 * field: determines whether 'CMError' (evoType = true) or 'DMError' (evoType = false)
	 * are used for evolution
	 */
	public boolean evoType;
	
	public boolean homogenous;
	
	public HashMap<Integer,Integer> is_homogen;
	
	/**
	 * field: all contacts with a frequency higher than a given threshold are saved here
	 */
	public HashMap<Pair<Integer>,Double> contacts_hfreq;
	//public String character;
	
	public HashMap<Integer,HashMap<Integer,Demes>> deme2;
	
	/**
	 * static field: specifies the directory to which all files will be written
	 */
	public static String path = "/project/StruPPi/gabriel/Arbeiten/";
	
	public String temp_path;
	
	/**
	 * Zero-parameter constructor - default values: <code> demenumb = 5 </code>; <code> evosteps = 15</code>;
	 * <code> deme.length = 20</code>.
	 * @throws SQLException
	 * @throws IOException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public Species () throws SQLException, IOException, PdbCodeNotFoundError, PdbLoadError{
		this.setPop(0,5,15,20,0,true);
	}
	
	/**
	 * One-parameter constructor: the parameter specifies which protein will be evolved. The remaining
	 * parameters are set to:
	 * <p>
	 * <code> deme.length = 10</code>
	 *</p>
	 *<p>
	 * <code> evosteps = 20</code>
	 * </p>
	 * <p>
	 * <code> popsize = 20</code>
	 * </p>
	 * @param protentry - specifies which protein of Sathyas set is used. Must be chosen between 0 and 14 otherwise an IllegalArgumentException is thrown.
	 * @throws SQLException
	 * @throws IOException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 * @throws IllegalArgumentException
	 */
	public Species (int protentry) throws SQLException, IOException, PdbCodeNotFoundError, PdbLoadError, IllegalArgumentException {
		boolean tester1 = (protentry >= 0) && (protentry < 15);		//test if the parameter matches the requirements
		
		if(tester1){												//test
			
			this.setPop(0,1,20,20,protentry,true);						//initializing the Population instance
			
		}
		else{														//if the parameter does not match the requirements
																	//an IllegalArgumentException is thrown
			
			String proterror = "Parameter 'protentry' must be chosen between 0 and 14, but was set to 'protentry = " + protentry + "'.";
			throw new IllegalArgumentException(proterror);
		}
	}
	
	
	/**
	 * One-parameter constructor: this is a constructor copying the <code>spe</code> instance to this
	 * <code> Species </code> instance.
	 * @param spe
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public Species (Species spe) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		this.demenumb = spe.demenumb;								//copies the field demenumb
		
		this.evosteps = spe.evosteps;								//copies the field evosteps
		
		int dim1 = spe.deme.length, dim2 = spe.deme.length;			//the dimensions of the Species matrix are copied
		
		this.deme = new Demes[dim1][dim2];						//instantiate the matrix
		
		for(int i = 0; i < dim1; i++){								//loop to initilize all Species present in the spe.deme field
			
			for(int j = 0; j < dim2; j++){
				this.deme[i][j] = new Demes(spe.deme[i][j]);
			}
		}
	}

	/**
	 * Two parameter constructor: this constructor is intended to allow evolution runs without fixed
	 * number of evolution steps. Before each evolution step the standard deviation of the averaged
	 * error value is tested against a threshold value, which is set by default to 0.003.
	 * @param protentry - protein entry in Sathyas list - must be in {0,1,...,14}
	 * @param dummy - dummy parameter, to distinguish this constructor from <code> {@link #Population (int)} </code>
	 * @throws SQLException
	 * @throws IOException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public Species (int protentry, String dummy) throws SQLException, IOException, PdbCodeNotFoundError, PdbLoadError{
		threshold = 0.0003;
		this.setDeme2(5, protentry, 14, 5, dummy,true);
	}
	
	/**
	 * Two parameter constructor: this constructor is intended to allow evolution runs without fixed
	 * number of evolution steps. Before each evolution step the standard deviation of the averaged
	 * error value is tested against a threshold value, which is set by default to 0.003.
	 * @param protentry - protein entry in Sathyas list - must be in {0,1,...,14}
	 * @param threshold_val - the threshold value to which the mean standard deviation of a fixed generation will be compared
	 * @throws SQLException
	 * @throws IOException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public Species (int protentry, int numofconts, double threshold_val, boolean evo) throws SQLException, IOException, PdbCodeNotFoundError, PdbLoadError{
		threshold = threshold_val;
		this.evosteps = 30;
		this.setDeme2(5, protentry, 20, numofconts, "dummy",evo);
	}
	
	/**
	 * Three parameter constructor: allows to specifically pick one protein out of Sathyas set for evolution. Additionally, the number of
	 * contacts can be chosen to in order to deal with proteins of different lenght. The default values are the population size and the number
	 * of individuals in each deme: pop_size = 20; and deme_size = 10;
	 * @param protentry
	 * @param num_of_conts
	 * @param dummy
	 * @throws SQLException
	 * @throws IOException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public Species (int protentry, int num_of_conts, String dummy, boolean evo) throws SQLException, IOException, PdbCodeNotFoundError, PdbLoadError{
		threshold = 0.0003;
		this.setDeme2(5,protentry, 20, num_of_conts, dummy,evo);
	}
	
	/**
	 * Two-parameter constructor: sets this instance of this class to the specified parameters.
	 * @param protentry - must be chosen between 0 and 14 otherwise an IllegalArgumentException is thrown.
	 * @param demesize - must always be greater than 0 otherwise an IllegalArgumentException is thrown.
	 * @throws SQLException
	 * @throws IOException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 * @throws IllegalArgumentException
	 */
	public Species (int protentry, int demesize, boolean evo) throws SQLException, IOException, PdbCodeNotFoundError, PdbLoadError, IllegalArgumentException {
		boolean tester1 = (protentry >= 0) && (protentry < 15), tester2 = demesize > 0;		//test if the two parameters match the requirements
		
		if(tester1 && tester2){																//test
			
			this.setPop(0, demesize, 1, 500,protentry,evo);									//initializing the Population instance
			
		}
		else{																				//if at least one requirements are not fullfilled
																							//an IllegalArgumentException is thrown
			if(!tester1){
				String proterror = "Parameter 'protentry' must be chosen between 0 and 14, but was set to 'protentry = " + protentry + "'.";
				throw new IllegalArgumentException(proterror);
			}
			if(!tester2){
				String demeerror = "Parameter 'demesize' must never be less than zero, but was set to 'demesize = " + demesize + "'.";
				throw new IllegalArgumentException(demeerror);
			}
		}
	}
	
	/**
	 * Three-parameter constructor: sets this instance of this class to the specified parameters. Has 20
	 * evolution steps as default value.
	 * @param protentry - must be chosen between 0 and 14 otherwise an IllegalArgumentException is thrown.
	 * @param demesize - must always be greater than 0 otherwise an IllegalArgumentException is thrown.
	 * @param popsize - number of <code> {@link #Individuals} </code> in each Population - must always be
	 * greater than 0 otherwise an IllegalArgumentException is thrown.
	 * @throws SQLException
	 * @throws IOException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 * @throws IllegalArgumentException
	 */
	public Species (int protentry, int demesize, int popsize, boolean evo) throws SQLException, IOException, PdbCodeNotFoundError, PdbLoadError{
		boolean tester1 = (protentry >= 0) && (protentry < 15), tester2 = demesize > 0;		//test if the three parameters match the requirements
		
		boolean tester3 = popsize > 0;
		
		if(tester1 && tester2 && tester3){													//test
			
			this.setPop(0, demesize, 20, popsize,protentry,evo);								//initializing the Population instance
			
		}
		else{																				//if at least one requirements are not fullfilled
																							//an IllegalArgumentException is thrown
			if(!tester1){
				String proterror = "Parameter 'protentry' must be chosen between 0 and 14, but was set to 'protentry = " + protentry + "'.";
				throw new IllegalArgumentException(proterror);
			}
			if(!tester2){
				String demeerror = "Parameter 'demesize' must never be less than zero, but was set to 'demesize = " + demesize + "'.";
				throw new IllegalArgumentException(demeerror);
			}
			if(!tester3){
				String demeerror = "Parameter 'popsize' must never be less than zero, but was set to 'popsize = " + popsize + "'.";
				throw new IllegalArgumentException(demeerror);
			}
		}
	}
	
	public Species (String dir) throws FileNotFoundException, IOException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		this.readStarterFromFile(dir);
	}
	
	public Species (Demes[] spec, int evstep,String dir, String temp_dir, double thres_val){
		this.setPop(spec, evstep, dir, temp_dir, thres_val);
	}
	
	public Species (Individuals in){
		this.setPop(in,2,2);
	}
	
	/*------------------------------------------------------setter--------------------------------------------------------------*/
	
	/**
	 * setter method: this setter is used from all constructors, the only default value is the field <code> evoType = false</code>, meaning
	 * only the DMError is taken into account. 
	 * @param generation
	 * @param demes
	 * @param evsteps
	 * @param popsize
	 * @param protentry
	 * @throws SQLException
	 * @throws IOException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public void setPop(int generation, int demes, int evsteps, int popsize, int protentry, boolean evo) throws SQLException, IOException, PdbCodeNotFoundError, PdbLoadError{
		this.demenumb = demes;											//setting the number of demes in this Species
		
		this.evosteps = evsteps;										//setting the number of evolution steps
		
		this.gen = generation;											//setting the generation
		
		this.evoType = evo;											//setting the field 'evoType' to false (meaning only DMError used)
		
		this.deme = new Demes[demes][evsteps];						//initializzing the field 'deme'
		
		this.pop_size = popsize;
		
		for(int i = 0; i < this.demenumb; i++){
			this.deme[i][generation] = new Demes(popsize,protentry,5,this.gen, "deme" + i);
			System.out.println(this.deme[i][0].getName() + " initialized...");					//output indicating the initialization of the parental Species
			
		}
		this.numofConts = this.deme[0][0].getNumOfContacts();
	}
	
	/**
	 * @TODO: needs to be worked and checked for consistency
	 */
	public void setPop(Individuals in, int deme_number, int ratio){
		int[][] contacts = SortIntArray.converter(in.getHashSet());
		//HashMap<Integer,HashSet<Integer>> surveilleur = new HashMap<Integer,HashSet<Integer>> (2*deme_number);
		for(int i = 0; i<deme_number;i++){
			//Integer deme_index = new Integer (i);
			int counter = 0;
			while(counter < ratio){
				//Random rand = new Random ();
				HashSet<Pair<Integer>> contactset = new HashSet<Pair<Integer>> ();
				//Integer array_index = new Integer (rand.nextInt(contacts.length));
				Pair<Integer> pair = new Pair<Integer> (new Integer (contacts[i][0]), new Integer (contacts[i][1]));
				contactset.add(pair);
				contacts = removeEntrie(contacts,i);
				counter++;
			}
		}
		
	}
	
	/**
	 * setter: copying <code> Species </code> matrix <code> pop </code> to this <code> Species </code>
	 * starting from <code> start </code> to <code> end </code>.
	 * @param start
	 * @param end
	 * @param pop
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public void setDemes(int start, int end, Demes[][] pop) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim = pop[0].length;
		for(int i = start; i < end; i++){
			for(int j = 0; j < dim; j++){
				this.deme[i][j] = new Demes(pop[i][j]);
			}
		}
	}
	
	public void setDemes (int pos, Demes[] input) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim1 = input.length;
		for(int i = 0; i < dim1; i++){
			this.deme[pos][i] = new Demes(input[i]);
		}
	}
	
	public void setDeme (int deme_num, int gen_num, Demes input) throws SQLException, PdbCodeNotFoundError, PdbLoadError {
		this.deme[deme_num][gen_num] = new Demes(input);
	}
	
	public void setContactFrq (double threshold, int start) throws FileNotFoundException, IOException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim1 = this.demenumb, dim2 = this.evosteps;
		HashMap<Pair<Integer>,Double> subset = new HashMap<Pair<Integer>,Double> (); 
		for(int i = start; i < dim1 - start; i++){
			for(int j = 0; j < dim2; j++){
				subset.putAll(this.deme[i][j].contactsWithHighestFrequency(threshold));
			}
		}
		this.contacts_hfreq = new HashMap<Pair<Integer>,Double> (subset);
		this.printToFile(path, "1bkr");
	}
	
	public void setContactFrq (int start, double threshold) throws FileNotFoundException, IOException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim = this.getDemeSize();
		HashMap<Pair<Integer>,Double> subset = new HashMap<Pair<Integer>,Double> ();
		try{
			for(int j = 0; j < dim; j++){
				subset.putAll(this.deme[j][start].contactsWithHighestFrequency(threshold));
			}
			this.contacts_hfreq = new HashMap<Pair<Integer>,Double> (subset);
			this.printToFile(path, "1bkr");
		}
		catch (NullPointerException e){
			System.out.println("Field 'contact_hfreq' was not initialized.");
		}
	}
	
	/**
	 * This setter allows evolution runs with non fixed evolution steps, by using the field <code> {@link #deme2} </code>.
	 * @param demesize
	 * @param protentry
	 * @param popsize
	 * @throws SQLException
	 * @throws IOException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public void setDeme2 (int demesize, int protentry, int popsize, int num_of_conts, String dummy, boolean evo) throws SQLException, IOException, PdbCodeNotFoundError, PdbLoadError{
		if((0 <= protentry) && (protentry < 15) && (demesize > 0)){
			int deme_size_mult = demesize*2;
			this.deme2 = new HashMap<Integer,HashMap<Integer,Demes>> (deme_size_mult);
			Demes[] starter_pop = new Demes[demesize];
			for(int i = 0; i < demesize; i++){
				HashMap<Integer,Demes> submap = new HashMap<Integer,Demes> ();
				Integer deme_index = new Integer (i);
				Integer gen_index = new Integer (0);
				starter_pop[i] = new Demes(popsize, protentry, num_of_conts, 0, "deme" + i);
				//Pair<Integer> pair = new Pair<Integer> (deme_index, gen_index);
				submap.put(gen_index, starter_pop[i]);
				this.deme2.put(deme_index, submap);
				System.out.println(starter_pop[i].getPop(0).getName() + ", deme " + i + " initialized...");
			}
			this.gen = 0;
			this.pop_size = popsize;
			this.demenumb = demesize;
			this.numofConts = starter_pop[0].getNumOfContacts();
			this.evoType = evo;
		}
		else{
			if(protentry < 0){
				String arrayind = "Protein index = " + protentry + " must never be less than zero.";
				throw new ArrayIndexOutOfBoundsException(arrayind);
			}
			else{
				if(protentry >= 15){
					String arrayind = "Protein index = " + protentry + " must never be greater or equal to 15.";
					throw new ArrayIndexOutOfBoundsException(arrayind);
				}
				else{
					String arrayind = "Deme size = " + demesize + " must never be less than zero.";
					throw new ArrayIndexOutOfBoundsException(arrayind);
				}
			}
		}
	}
	
	/**
	 * Setter, initializing the field <code>{@link #deme2}</code>. Works in the same fashion as
	 * setter <code>{@link #setPop(int, int, int, int, int)}</code>.
	 * @param demesize
	 * @param protentry
	 * @param popsize
	 * @param evsteps
	 * @throws SQLException
	 * @throws IOException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public void setDeme2 (int demesize, int protentry, int popsize, int evsteps) throws SQLException, IOException, PdbCodeNotFoundError, PdbLoadError{
		this.demenumb = demesize;
		this.evosteps = evsteps;
		this.deme = new Demes[demesize][evsteps];
		this.deme2 = new HashMap<Integer,HashMap<Integer,Demes>> (2*demesize);
		for(int i = 0; i < demesize; i++){
			HashMap<Integer,Demes> submap = new HashMap<Integer,Demes> (2*evsteps);
			Integer gen_index  = new Integer (0);
			Integer deme_index = new Integer (i);
			this.deme[i][0] = new Demes (popsize, protentry, 5, 0, "deme" + i);
			submap.put(gen_index, this.deme[i][0]);
			this.deme2.put(deme_index, submap);
		}
		this.pop_size = popsize;
	}
	
	/**
	 * Method reading the default 'temp' files which are written during each run of evolution. To each run
	 * a temporary directory with all <code> {@link #Individuals} </code> instances as contact map files is made.
	 * In case of a severe error or exception occuring, each run can in principle be continued, just by
	 * calling the one-parameter constructor <code> {@link #Species (String)}</code>. 
	 * @param dir - absolut path of the temporary files are stored
	 * @throws IOException
	 * @throws FileNotFoundException
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 * @throws IllegalArgumentException - thrown if either the pdb code of two different Individuals are not the same or
	 * if the denoted path does not contain any contact map files.
	 * @throws NullPointerException the denoted path does not exist
	 */
	public void readStarterFromFile (String dir) throws IOException, FileNotFoundException, SQLException, PdbCodeNotFoundError, PdbLoadError , IllegalArgumentException , NullPointerException {
		File filedir = new File(dir);
		if(filedir.exists()){																					//tests whether this abstract path exists
			
			File[] filelist = filedir.listFiles(new RegexFileFilter("deme.*"));								//listing of all files in the denoted directory
			
			Arrays.sort(filelist);																				//sorting of the file list

			
			int length = filelist.length;
			this.evosteps = 20 - this.gen;															//number of evolution steps are set to default: 20 - gen (generation, where stopped)
			this.deme = new Demes[length][this.evosteps];
			if(length > 0){																						//tests whether this directory contains any cmap files
																												//if not: IllegalArgumentException thrown
				
				for(int k = 0; k < length; k++){
					
				
					File[] indiarray = filelist[k].listFiles(new RegexFileFilter (".*.cmap"));
					
					Arrays.sort(indiarray);
					
					int length1 = indiarray.length;
					
					//Species starter = new Species();																//new Species instance

					Individuals[] starterindis = new Individuals[length1];											//new array of Individuals instances

					for(int i = 0; i < length1; i++){																//looping over the files list

						//String dirstring = indiarray[i].getAbsolutePath().replace(dir,"");
						//if(!dirstring.contains("deme")){															//tests whether the file name contains keyword 'deme', if it does
							//contain it: skipping this file since Population contact files
						//}																							//only contain average info

						//if(dirstring.split("-").length == 3){

							BufferedReader reader = new BufferedReader(new FileReader(indiarray[i]));				//file reader

							HashSet<Pair<Integer>> contactset = new HashSet<Pair<Integer>> ();						//HashSet for the Individuals field 'store'

							String linereader = new String();														//a String instance reading each line of given file

							starterindis[i] = new Individuals();													//initializing the i-th Individuals instance

							while((linereader = reader.readLine()) != null){										//looping over the file

								if(linereader.contains("#")){														//all header lines of a contact map file start with '#'
									//representing all important fields of a single Individuals instance

									if(linereader.contains("SEQUENCE")){											//field: seq

										String seq = linereader.replace("#SEQUENCE: ", "");
										starterindis[i].setSequence(seq);
									}
									if(linereader.contains("PDB") && !linereader.contains("CHAIN")){				//field: name - pdb code

										starterindis[i].setName(linereader.replace("#PDB: ",""));
									}
									if(linereader.contains("PDB CHAIN")){											//field: chainCode

										starterindis[i].setChainCode(linereader.replace("#PDB CHAIN CODE: ", ""));
									}
									if(linereader.contains("CMError")){												//field: CMError

										starterindis[i].setCMError(linereader.replace("#CMError: ",""));
									}
									if(linereader.contains("DMError")){												//field: DMError

										starterindis[i].setDMError(linereader.replace("#DMError: ",""));
									}
									if(linereader.contains("NUMB. CONTACTS")){										//field: numOfContacts

										double nums = Double.parseDouble(linereader.replace("#NUMB. CONTACTS: ",""));
										int val = (int) nums;
										starterindis[i].setNumOfContact(val);
									}
									if(linereader.contains("NUMB. OF ALL CONTACTS")){								//field: fullContacts

										double nums = Double.parseDouble(linereader.replace("#NUMB. OF ALL CONTACTS: ",""));
										int val = (int) nums;
										starterindis[i].setFullContact(val);
									}
								}
								else{																				//the remaining content of a contact map file contains
									//the contacts in tabular form

									String[] array = linereader.split("\t");
									Integer contact1 = Integer.decode(array[0]);
									Integer contact2 = Integer.decode(array[1]);
									Pair<Integer> pair = new Pair<Integer> (contact1,contact2);
									contactset.add(pair);
								}
							}
							starterindis[i].storer(contactset);														//field: store initialized
							reader.close();
						//}
						/*else{																						//if the file name contains the keyword 'deme' the file
							//is checked for the entry 'GENERATION' to initialize
							//the field: gen of this Species instance

							BufferedReader reader = new BufferedReader(new FileReader(filelist[i]));
							String linereader = new String();
							while((linereader = reader.readLine()) != null){
								if(linereader.contains("GENERATION")){
									String num = linereader.replace("#GENERATION: ", "");
									double dnum = Double.parseDouble(num);
									this.gen = (int) dnum;															//field: gen - number of generation at which the evolution
									//was somehow interrupted
									break;
								}
							}*/
							//starterindis = Species.trim(starterindis);												//starter are trimmed in case of null entries in the array

							reader.close();
					}
					//length = starterindis.length;
					for(int j = 0; j < length1 - 1; j++){															//loop to check whether the i-th and i + 1st Individuals instances
						//have the same pdb code, otherwise an IllegalArgumentException is thrown

						int lastentry = starterindis[j + 1].getName().toCharArray().length;
						CharSequence s = starterindis[j + 1].getName().subSequence(0,lastentry);
						if(!(starterindis[j].getName().contentEquals(s))){
							String error = "The pdb code of Individuals " + j + " and Individuals" + (j + 1) + " do not match.\nMake sure only cmap.file with the same pdb code are place in the same directory.";
							//the Individuals.breed(Individuals,Individuals) throws as well an IllegalArgumentException if two pdb codes do not match

							throw new IllegalArgumentException(error);
						}
					}
					//starter = new Species(starterindis,0);

					this.deme[k][0] = new Demes(starterindis,0);															//starter are initialized

				}
				this.numofConts = this.deme[0][0].getPop(0).getNumOfContacts();
			}
			else{																								//if the denoted directory does not contain any contact map files

				String nullfiles = "The denoted path: " + dir + " does not contain any '.cmap' files.";
				throw new IllegalArgumentException(nullfiles);
			}
		}
		else{																									//if the denoted directory does not exist
			String nodir = "The path: " + dir + " does not exist.";
			throw new NullPointerException(nodir);
		}
	}
	
	/**
	 * method copying the parental generation of each deme from either the field <code>deme</code> to <code>deme2</code> or visa versa.
	 * This is helpful since both evolution methods <code>evolve()</code> and <code>evolve(String)</code> do not behave in the same manner.
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public void copyParentalDemes () {//throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		if(this.deme != null){
			this.deme2 = new HashMap<Integer,HashMap<Integer,Demes>> ();
			HashMap<Integer,Demes> submap = new HashMap<Integer,Demes> (2*this.evosteps);
			Integer gen_index = new Integer (0);
			for(int i = 0; i < this.demenumb; i++){
				Integer deme_index = new Integer(i);				
				Demes spec = new Demes (this.deme[i][0],0);
				submap.put(gen_index, spec);
				this.deme2.put(deme_index, submap);
			}
		}
		else{
			if(this.deme2 != null){
				HashMap<Integer,Demes> deme_map = this.getDemesOfSameGen(0);
				Set<Integer> keyset = deme_map.keySet();
				Iterator<Integer> it = keyset.iterator();
				this.deme = new Demes[this.demenumb][this.evosteps];
				while(it.hasNext()){
					Integer deme_index = it.next();
					int deme_val = deme_index.intValue();
					this.deme[deme_val][0] = new Demes(this.getDeme2(deme_val, 0),0);
				}
			}
		}
	}
	
	/**
	 * setter needed in order to create an <code>Population</code> instance by an array of <code>Species</code>.
	 * 
	 * @param spec
	 * @param evo_steps
	 * @param dir
	 * @param temp_dir
	 * @param thres_val
	 */
	public void setPop (Demes[] spec, int evo_steps, String dir, String temp_dir, double thres_val){
		threshold = thres_val;
		this.temp_path = new String (temp_dir); path = new String (dir);
		this.pop_size = spec[0].getSize();
		this.gen = 0;
		this.numofConts = spec[0].getNumOfContacts();
		this.demenumb = spec.length;
		this.evosteps = evo_steps;
		this.deme = new Demes [this.demenumb][evo_steps];
		this.deme2 = new HashMap<Integer,HashMap<Integer,Demes>> (this.demenumb *2);
		for(int i = 0; i < this.demenumb; i++){
			HashMap<Integer,Demes> submap = new HashMap<Integer,Demes> (evo_steps*2);
			this.deme[i][0] = new Demes (spec[i],0);
			Integer null_val = new Integer (0), deme_index = new Integer (i);
			submap.put(null_val, spec[i]);
			this.deme2.put(deme_index, submap);
		}
	}
	
	/*--------------------------------------------------evolution-----------------------------------------------------------*/
	
	/**
	 * Evolution method, each deme of this <code> Species </code> is evolved according to the 
	 * number of generations specified by the field <code> evosteps </code>. This method uses
	 * the method <code> evolve (int) </code> and copys each Population the the field <code> deme </code>
	 * using the setter <code> setDemes (int, Population[]) </code>. Since this method on acts on the field
	 * <code>{@link #deme}</code> it should only be used, when initialization of the specific field is
	 * guaranteed.
	 * @deprecated instead of calling this method, rather rely on <code>{@link #evolve(String)}</code>. 
	 * @throws FileNotFoundException
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 * @throws IOException
	 */
	@Deprecated
	public void evolve () throws FileNotFoundException, SQLException, PdbCodeNotFoundError, PdbLoadError, IOException{
		int dim1 = this.getDemeSize();				//number of demes in this Species
		
		for(int i = 0; i < dim1; i++){				//loop over all demes in this Species
			
			this.evolve(i);							//the actual evolution using the
		}
		this.setContactFrq(0.7,15);
	}

	/**
	 * Primary evolution method, each deme of this <code> Species </code> is evolved according to the 
	 * number of generations specified by the field <code> evosteps </code>. This method uses
	 * the method <code>{@link #evolve(int,String)}</code> and copies each Population the the field <code> deme </code>
	 * using the setter <code> setDemes (int, Population[]) </code>. Since this method on acts on the field
	 * <code>{@link #deme}</code> it should only be used, when initialization of the specific field is
	 * guaranteed.
	 * @param dummy
	 * @throws FileNotFoundException
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 * @throws IOException
	 */
	public void evolve (String dummy) throws FileNotFoundException, SQLException, PdbCodeNotFoundError, PdbLoadError, IOException{
		int dim1 = this.getDemeSize();				//number of demes in this Species

		int dim2 = this.getEvoSteps();				//number of evolution steps

		for(int j = 0; j < dim2-1; j++){

			if(j >= 10){	//after 5 generations all demes are checked for homogeneity
				
				boolean homo_tester = this.isGenXHomogenous(j);
				if(homo_tester){
					this.swapSpecies(j);
				}
			}
			for(int i = 0; i < dim1; i++){				//loop over all demes in this Species

				this.evolve(i,j);				//the actual evolution using the

			}
			this.gen = j;
		}
		System.out.println("Evolution of type 1 terminated..");
	}
	
	/**
	 * Evolution method, evolves the 'k'-th deme of this <code> {@link #Species} </code> at the given generation. Each <code> {@link #Population} </code>
	 * is written to a contact map file ('.cmap' extension). The evolution is performed on the field <code>{@link #deme}</code>.
	 * @param k - evolving the k-th deme
	 * @return
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	@Deprecated
	public Demes[] evolve (int k) throws SQLException, PdbCodeNotFoundError, PdbLoadError, FileNotFoundException, IOException{
		int generations = this.gen;
		int steps = this.evosteps;												//number of evolution steps
		
		Demes pop = new Demes(this.getDeme(k,generations));					//starter species at k-th generation
		
		//String name = pop.getName() + "gen" + pop.getGen();						//a string to specify the output file name
		String name = String.format("%sgen%02d", pop.getName(), pop.getGen());
		
		String printname = pop.getPop(0).getName();								//a string to specify the output file name
		
		Demes[] off = new Demes[steps];										// a new array of Species instances, with
																				//as many entries as specified by the field evosteps
		
		off[0] = new Demes(pop,this.gen);										//initializing the parental Species
		
		off[0].printToFile(path + printname + "/", name);						//writing the parental Species to the specified file
																				//as a contact map file
		
		off[0].printMetric(path + printname + "/", name,0);						//writing the parental Species to the specified file
																				//as a metric file, a metric file contains only the
																				//calculated pairwise distances between all Individuals
																				//of the given Species using the metric setter methods of the
																				//Species class
		
		File dirtest = new File(path + printname + "/tempfiles/");				//creating a temporary directory to save the current Species
																				//by writing each Individual to a contact map file
		
		if(!dirtest.exists()){													//test, if the specified path exists
			
			dirtest.mkdir();													//if the specified path does not exist, create the desired directory
			
		}
		for(int i = 1; i < steps; i++){											//the actual evolution loop
			
			String name1 = String.format("%sgen%02d", pop.getName(), i);;							//a helper string for the contact map files
			
			off[i] = new Demes(this.getDeme(k,i).evolve(false,"dummy"),i);	//the actual evolution step using the specific
																										//evolution method of the Species class

									
			off[i].printToFile(path + printname + "/", name1);					//writing the contact map file of the i-th generation
			
			off[i].printMetric(path + printname + "/", off[i].getName(),k);		//writing the metric file of the i-th generation
			
			this.setDeme(k, i, off[i]);
			this.printDemesToTempFiles(k,i);
		}
		return off;																//returning the Species array of this evolution run
	}
	
	/**
	 * Evolution method, evolves the 'k'-th deme of this <code> {@link #Species} </code> at the given generation. Each <code> {@link #Population} </code>
	 * is written to a contact map file ('.cmap' extension). The evolution is performed on the field <code>{@link #deme}</code>.
	 * @param k - evolving the k-th deme
	 * @param genindex - the generation index
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public void evolve (int k, int genindex) throws SQLException, PdbCodeNotFoundError, PdbLoadError, FileNotFoundException, IOException{
		int generations = genindex;
		//int steps = this.evosteps;												//number of evolution steps
		
		Demes pop = new Demes(this.getDeme(k,generations));					//starter species at k-th generation
		
		String name = pop.getPop(0).getName();						//a string to specify the output file name
		
		String printname = name+"/deme"+k+"/evo1/";								//a string to specify the output file name
		
		//Species off = new Species();											// a new array of Species instances, with
																				//as many entries as specified by the field evosteps
		
		//off = new Species(pop,this.gen);										//initializing the parental Species
		
		//off.printToFile(path + printname + "/", name);						//writing the parental Species to the specified file
																				//as a contact map file
		
		//off.printMetric(path + printname + "/", name,0);						//writing the parental Species to the specified file
																				//as a metric file, a metric file contains only the
																				//calculated pairwise distances between all Individuals
																				//of the given Species using the metric setter methods of the
																				//Species class
		
		File dirtest = new File(path + printname + "/tempfiles/");				//creating a temporary directory to save the current Species
																				//by writing each Individual to a contact map file
		
		if(!dirtest.exists()){													//test, if the specified path exists
			
			dirtest.mkdir();													//if the specified path does not exist, create the desired directory
			
		}			
		String name1 = String.format(name+"deme%sgen%02d",k,generations+1);					//a helper string for the contact map files
		
		boolean evo = this.evoType;
		
		Demes offi = new Demes(this.getDeme(k,genindex).evolve(evo, "dummy"),genindex+1);
		//the actual evolution step using the specific evolution method of the Species class
		
			this.setDeme(k, genindex+1, offi);
			offi.printAllIndisToTempDir(path,name+"/deme"+k+"/temp1/");
			
			
				//this.setContactFrq(genindex,0.5);
									
			offi.printToFile(path + printname + "/", name1);					//writing the contact map file of the i-th generation
			
			offi.printMetric(path + printname + "/", offi.getName(),k+1);		//writing the metric file of the i-th generation
		
		//return off;																//returning the Species array of this evolution run
	}
	
	/**
	 * Evolution method: the evolution is performed on the field <code>{@link #deme2}</code> - resulting in
	 * a slightly different approach, since this field is a HashMap instance. This method is called by the
	 * method <code>{@link #evolve2(int)}</code>, so it should never be used as 'true' evolutionary method. Instead,
	 * use the abovementioned method.
	 * @throws FileNotFoundException
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 * @throws IOException
	 */
	public void evolve2 () throws FileNotFoundException, SQLException, PdbCodeNotFoundError, PdbLoadError, IOException{
		HashMap<Integer,HashMap<Integer,Demes>> map = this.getDeme2();
		Set<Integer> demeset = map.keySet();
		Iterator<Integer> it = demeset.iterator();
		if(threshold == 0.0){
			threshold = 0.003;
		}
		int gennum = this.getGen();
		//if(gennum >= 10){
			boolean isless = this.isGenXHomogenous(gennum);//isLessThanThreshold(gennum, false, threshold);
			if(isless){
				this.swapSpecies(gennum);
				map = this.getDeme2();
				demeset = map.keySet();
				it = demeset.iterator();
			}
		//}
		while(it.hasNext()){
			Integer deme_index = it.next();
			this.evolve2(deme_index.intValue(),gennum);
		}
		this.gen += 1;
	}
	
	/**
	 * Primary evolution method: the evolution is performed on the field <code>{@link #deme2}</code> - resulting in
	 * a slightly different approach, since this field is a HashMap instance. This method calls the
	 * method <code>{@link #evolve2()}</code>, so should never be used as 'true' evolutionary method. Instead,
	 * use the abovementioned method. 
	 * @param final_generation - the last generation in the evolution (inclusive)
	 * @throws FileNotFoundException
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 * @throws IOException
	 */
	public void evolve2 (int final_generation) throws FileNotFoundException, SQLException, PdbCodeNotFoundError, PdbLoadError, IOException{
		int currentgen = this.gen;
		this.evosteps = final_generation + 1;
		//boolean homogenity = false;
		//if(currentgen >= 20){
			//homogenity = this.isGenXHomogenous(currentgen);
		//}
		boolean homogeneity = this.isGenXHomogenous(currentgen);
		if(currentgen >= 10){
			while(currentgen < final_generation && !homogeneity){
				homogeneity = this.isGenXHomogenous(currentgen);
				this.evolve2();
				currentgen = this.gen;
			}
		}
		else{
			if(!homogeneity){
				if(currentgen == 0){
					for(int i = 0; i < this.demenumb; i++){
						Demes spec = this.getDeme2(i, 0);
						String name = spec.getName().replace("deme" + i, "");
						File dirs = new File (path + name + "/deme" + i + "/evo2/");
						if(!dirs.exists()){
							dirs.mkdirs();
						}
						String helpstr = String.format("deme%sgen%02d", i, 0);
						spec.printToFile(path, name + "/deme" + i + "/evo2/" + name + helpstr);
						spec.printMetric(path, name + "/deme" + i + "/evo2/" + name + helpstr, 0);
					}
				}
				if(currentgen >= 10){
					System.out.println("Help line");
				}
				while(currentgen < final_generation){
					this.checkForHomogeneity(currentgen,0.8);
					if(this.homogenous){
						break;
					}
					this.evolve2();
					currentgen = this.gen;
				}
				if(this.homogenous){
					System.out.println("Evolution was stopped due to the lack of variability within the population at generation " + currentgen + ".");
				}
				else{
					System.out.println("Evolution of type 2 terminated...");
				}
			}
			else{
				System.out.println("Evolution was stopped due to the lack of variability within the population at generation " + currentgen + ".");
			}
		}
	}
	
	/**
	 * Evolution method: the evolution is performed on the field <code>{@link #deme2}</code> - resulting in
	 * a slightly different approach, since this field is a HashMap instance. This method is called by the
	 * method <code>{@link #evolve2()}</code>, so should never be used as 'true' evolutionary method. Instead,
	 * use the <code>{@link #evolve2(int)}</code> method.
	 * @param i
	 * @param current_gen
	 * @throws FileNotFoundException
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 * @throws IOException
	 */
	public void evolve2 (int i, int current_gen) throws FileNotFoundException, SQLException, PdbCodeNotFoundError, PdbLoadError, IOException{
		Integer deme_index    = new Integer (i);
		//current generation index as an Integer instance		
		
		Integer gen_plus_one  = new Integer(current_gen + 1);
		HashMap<Integer,Demes> submap    = this.getDeme2(deme_index.intValue());
		if(!submap.containsKey(gen_plus_one)){
			boolean evo = this.evoType;
			Demes off = new Demes (this.getDeme2(i,current_gen).evolve(evo, "dummy"), gen_plus_one);
			String name = off.getName();
			String dir_name = path + "/" + name;
			File dir = new File (dir_name + "/deme"+ i + "/evo2/");
			if(!dir.exists()){
				dir.mkdirs();
			}
			File dir1 = new File (dir_name + "/deme" + i + "/temp/");
			if(!dir1.exists()){
				dir1.mkdirs();
			}
			String helpstr = String.format("deme%sgen%02d", i, (current_gen + 1));
			off.printAllIndisToTempDir(path, name + "/deme" + i + "/temp2/");
			off.printMetric(path, name + "/deme" + i + "/evo2/" + name + helpstr, i);
			off.printToFile(path, name + "/deme" + i + "/evo2/" + name + helpstr);
			//Pair<Integer> pair = new Pair<Integer> (deme_index, gen_plus_one);
			submap.put(gen_plus_one, off);
			this.deme2.put(deme_index, submap);
			//this.gen = current_gen + 1;
		}
	}
	
	/**
	 * This auxiliary method chekcs whether all contact are already present in the Population and therefore
	 * any evolution will not yield any better result than the current one.
	 * @param gen_index
	 */
	public void checkForHomogeneity(int gen_index, double threshold_val){
		if(deme2 != null){
			HashMap<Integer,Demes> submap = this.getDemesOfSameGen(gen_index);
			Integer zero = new Integer (0);
			HashSet<Pair<Integer>> common_contacts = new HashSet<Pair<Integer>> (submap.get(zero).getSubHash());
			submap.remove(zero);
			Set<Integer> keyset = submap.keySet();
			Iterator<Integer> it = keyset.iterator();
			while(it.hasNext()){
				Integer deme_index = it.next();
				common_contacts.retainAll(submap.get(deme_index).getSubHash());
			}
			int frequency = (int) (((double) common_contacts.size())/threshold_val);
			if(frequency >= this.numofConts){
				this.homogenous = true;
			}
		}
		else{
			if(this.deme != null){
				HashSet<Pair<Integer>> common_contacts = new HashSet<Pair<Integer>> (this.deme[0][gen_index].getSubHash());
				for(int i = 1; i< this.demenumb;i++){
					common_contacts.retainAll(this.deme[i][gen_index].getSubHash());
				}
				int frequency = (int) (((double) common_contacts.size())/threshold_val);
				if(frequency >= this.numofConts){
					this.homogenous = true;
				}
			}
		}
	}
	
	public void printStarterGenerationToFile (String dir, String addname) throws FileNotFoundException, IOException{
		if(this.deme2 != null){
			HashMap<Integer,Demes> starterset = this.getDemesOfSameGen(0);
			Set<Integer> keyset = starterset.keySet();
			Iterator<Integer> it = keyset.iterator();
			int counter = 0;
			while(it.hasNext()){
				Integer index = it.next();
				Demes spec = starterset.get(index);
				spec.printAllIndisToTempDir(dir, addname + counter);
				counter++;
			}
		}
		else{
			if(this.deme != null){
				for(int i = 0; i < this.getDemeSize(); i++){
					Demes spec = new Demes (this.deme[i][0]);
					spec.printAllIndisToTempDir(dir, addname + 1);
				}
			}
		}
	}
	
	/*---------------------------------------------------display------------------------------------------------------------*/
	
	
	/**
	 * Printing method: creates 'cmap' files by extracting are important informations of each of the demes
	 * and write it to a specified file. This method is called by each evolution methods, defined in this
	 * class. Each <code>{@link #Species}</code> present in either one of the fields <code>{@link #deme}</code>
	 * or <code>{@link #deme2}</code> is written to a 'cmap' file.
	 * <p>
	 * the standard output is as follows:
	 * </p>
	 * <p>
	 * the header: with sequence, pdb code, chain code, generation, error values and deviations,
	 * </p>
	 * <p>
	 * an the actual contact map table, where the first and second column correspond to the first and second contact pair
	 * and the third column is the frequency.
	 * </p>
	 * @param path a String denoting the output directory
	 * @param name a String specifying the output file name
	 */
	public void printToFile(String path, String name) throws IOException, FileNotFoundException, SQLException, PdbCodeNotFoundError, PdbLoadError {
		Individuals ne = new Individuals(this.getDeme(0,0).getPop(0));
		String cont = "#CMVIEW GRAPH FILE ver: 1.0\n#SEQUENCE: "+ne.getSequence()+"\n"+
		"#PDB: "+ne.getName()+ "\n#PDB CHAIN CODE: "+ne.getChainCode()+"\n#CT: "+Individuals.getContactT()+ "\n#CUTOFF: "+Individuals.getContactDist()+"\n"+ 
		"#GENERATION: "+0 + "\n#Species SIZE: " + this.getDeme(0,0).getSize() + "\n#CMError: "+ 0 + "\n#CMstDev: " + 0 + "\n#DMError: "+ 0 + 
		"\n#DMstDev: " + 0 + "\n#NUMB. CONTACTS: " + ne.getNumOfContacts() +"\n" + "#NUMB. OF ALL CONTACTS: " + ne.getFullContact() + "\n";
		HashSet<Pair<Integer>> keyset = new HashSet<Pair<Integer>> (this.contacts_hfreq.keySet());
		int[][] index_array = SortIntArray.converter(keyset);
		int length = index_array[0].length;
		//TreeSet<CompPair> trset = new TreeSet<CompPair> (trmap.keySet());
		HashMap<Pair<Integer>,Double> hashmap = new HashMap<Pair<Integer>,Double> (this.contacts_hfreq);
		//HashSet<Pair<Integer>> hashset = this.getKey();
		//Iterator<CompPair> it = trset.iterator();
		cont = cont + "\n";
		//while(it.hasNext()){
		for(int i = 0; i < length; i++){
			int index1 = index_array[0][i], index2 = index_array[1][i];
			cont = cont + index1 + "\t" + index2 + "\t" + hashmap.get(new Pair<Integer>(new Integer(index1),new Integer(index2))) + "\n"; 
		}
		FileOutputStream file = new FileOutputStream(path+name+"-cons"+this.getDeme(0,0).getNumOfContacts()+".cmap");
		PrintStream printa = new PrintStream(file);
		printa.print(cont);
		printa.close();
		file.close();
		System.out.println(ne.getName()+" at generation "+this.gen+" written to file...");
	}
	
	public Demes[] copyPops(Demes[] input) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim1 = input.length;
		Demes[] pop = new Demes[dim1];
		for(int i = 0; i < dim1; i++){
			pop[i] = new Demes(input[i]);
		}
		return pop;
	}
	
		
	public void writeDemeToDeme2 () throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		if(this.deme != null){
			int deme_size = this.deme.length;
			HashMap<Integer, HashMap<Integer, Demes>> map = new HashMap<Integer,HashMap<Integer,Demes>> (deme_size*2);
			for(int i = 0; i < deme_size; i++){
				int gen_size = this.deme[i].length;
				Integer deme_index = new Integer (i);
				HashMap<Integer,Demes> submap = new HashMap<Integer,Demes> (gen_size);
				for(int j = 0; j < gen_size; j++){
					Integer gen_index = new Integer (j);
					Demes spec = new Demes (this.deme[i][j],j);
					submap.put(gen_index, spec);
				}
				map.put(deme_index, submap);
			}
			this.deme2 = new HashMap<Integer,HashMap<Integer,Demes>> (map);
		}
	}
	
	public void writeDeme2ToDeme () throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		if(this.deme2 != null){
			int deme_size = this.demenumb;
			int gen_size = this.gen;
			this.deme = new Demes [deme_size][gen_size];
            for(int i = 0; i < deme_size; i++){
            	for(int j = 0; j < gen_size; j++){
            		this.deme[i][j] = this.getDeme2(i, j);
            	}
            }
		}
	}
	
	public int getLargestSubHashMap (){
		HashMap<Integer,HashMap<Integer,Demes>> map = this.getDeme2();
		Set<Integer> keyset = map.keySet();
		Iterator<Integer> it = keyset.iterator();
		int size = 0;
		while(it.hasNext()){
			HashMap<Integer,Demes> submap = map.get(it.next());
			if(submap.size() >= size){
				size = submap.size();
			}
		}
		return size;
	}
	
	public HashMap<Integer,Demes> getDemesOfSameGen (int gen_index) {//throws ArrayIndexOutOfBoundsException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		HashMap<Integer,HashMap<Integer,Demes>> deme_map = this.getDeme2();
		HashMap<Integer,Demes> map = new HashMap<Integer,Demes> ();
		Set<Integer> keyset = deme_map.keySet();
		Iterator<Integer> it = keyset.iterator();
		while(it.hasNext()){
			Integer deme_index = it.next();
			Demes spec = this.getDeme2(deme_index.intValue(), gen_index);
			map.put(deme_index, spec);
		}
		return map;
	}
	
	public double[] getErrorDev (int gen_index, boolean CMDM) throws ArrayIndexOutOfBoundsException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		HashMap<Integer,Demes> map = this.getDemesOfSameGen(gen_index);
		Set<Integer> keyset = map.keySet();
		double[] stdev_array = new double[keyset.size()];
		Iterator<Integer> it = keyset.iterator();
		if(CMDM){
			while(it.hasNext()){
				Integer index = it.next();
				stdev_array[index.intValue()] = map.get(index).getCMstdev();
			}
		}
		else{
			while(it.hasNext()){
				Integer index = it.next();
				stdev_array[index.intValue()] = map.get(index).getDMstdev();
			}
		}
		return stdev_array;
	}
	
	public double getMeanStDev(int gen_index, boolean CMDM) throws ArrayIndexOutOfBoundsException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		double[] st_dev = this.getErrorDev(gen_index, CMDM);
		int dim =  st_dev.length;
		double sum = 0.0;
		for(int i = 0; i < dim; i++){
			sum += st_dev[i];
		}
		return sum/((double) dim);
	}
	
	public boolean isLessThanThreshold (int gen_index, boolean CMDM, double threshold) throws ArrayIndexOutOfBoundsException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		double average = this.getMeanStDev(gen_index, CMDM);
		if(average <= threshold){
			return true;
		}
		else{
			return false;
		}
	}
	
	/*public void isPopHomogenous (int gen) throws ArrayIndexOutOfBoundsException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		HashMap<Integer,Species> pop_map = this.getDemesOfSameGen(gen);
		Set<Integer> keyset = pop_map.keySet();
		if(this.is_homogen == null){
			this.is_homogen = new HashMap<Integer,Integer> ();
		}
		Iterator<Integer> it = keyset.iterator();
		Integer gen_index = new Integer (gen);
		double metric = 0.0, metric2 = 0.0;
		while(it.hasNext()){
			Integer deme_ind = it.next();
			metric += pop_map.get(deme_ind).getAverageMetric();
			metric2+= pop_map.get(deme_ind).getAverageMetric2();
		}
		if(metric == 0.0 || metric2 == 0.0){
			this.is_homogen.put(gen_index,new Integer (1));
		}
		else{
			this.is_homogen.put(gen_index, new Integer (0));
		}
	}*/
	
	
	
	public boolean isGenXHomogenous (int gen) throws ArrayIndexOutOfBoundsException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		if(this.deme2 != null){
			HashMap<Integer,Demes> deme_map = this.getDemesOfSameGen(gen);
			Integer deme0 = new Integer (0);
			HashSet<Pair<Integer>> pairs = new HashSet<Pair<Integer>> (deme_map.get(deme0).getSubHash());
			int sum = pairs.size();
			deme_map.remove(deme0);
			Set<Integer> keyset = deme_map.keySet();
			Iterator<Integer> it = keyset.iterator();
			while(it.hasNext()){
				Integer deme = it.next();
				HashSet<Pair<Integer>> subhash = new HashSet<Pair<Integer>> (deme_map.get(deme).getSubHash());
				if(subhash.size() != 0){
					sum += subhash.size();
				}
			}
			int benchmark = (int) (((double) this.numofConts)*0.8);
			int av_frequency = (int) (((double) sum)/((double) this.demenumb));
			if(av_frequency >= benchmark){
				System.out.println("Population is homogenous...");
				return true;
			}
			else{
				System.out.println("Population is inhomogenous...");
				return false;
			}
		}
		else{
			if(this.deme != null){
				HashSet<Pair<Integer>> pairs = new HashSet<Pair<Integer>> (this.deme[0][gen].getSubHash());
				int sum_size_of_common_contacts = pairs.size();
				for(int i = 1; i < this.demenumb; i++){
					sum_size_of_common_contacts += this.deme[i][gen].getSubHash().size();
					if(pairs.size() == 0){
						break;
					}
				}
				double ratio = (double)sum_size_of_common_contacts/((double) this.demenumb); 
				int benchmark = (int) (((double) this.numofConts)*0.8);
				if(ratio >= benchmark){
					System.out.println("Population is homogenous...");
					return true;
				}
				else{
					System.out.println("Population is inhomogenous...");
					return false;
				}
			}
			else{
				throw new NullPointerException ("Neither of the fields 'deme' nor 'deme2' was initialized before calling this method!!");
			}
		}
	}
	
	public void swapSpecies (int gen_index) throws ArrayIndexOutOfBoundsException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		if(this.deme2 != null){
		HashMap<Integer,HashMap<Integer,Demes>> population = this.getDeme2();			//getting the whole population
		
		Integer gen_num = new Integer (gen_index);										//generation index
		
		HashMap<Integer,Demes> submap1 = this.getDemesOfSameGen(gen_index);			//getting the all species at the specified generation
		
		HashSet<Integer> surveillieur = new HashSet<Integer> ();						//checker, checks that no index is taken twice
		
		Set<Integer> keyset = submap1.keySet();											//keyset to iterate over the species of the specified generation
		
		Iterator<Integer> it = keyset.iterator();
		Individuals[] array = convertHashMapToIndiArray(submap1, this.pop_size);		//converting the species hashmap to an individuals array
		
		int size = array.length;//, deme_size = this.demenumb;
		/*int[] deme_index = new int[size];
		for(int i = 0; i < size; i++){
			Random rand = new Random ();
			deme_index[i] = rand.nextInt(deme_size);
		}
		for(int i = 0; i < deme_size; i++){
			Integer deme_val = new Integer (i);
			HashMap<Integer,Species> submap = population.get(deme_val);
			Individuals[] subpop = new Individuals[this.pop_size];
			int counter = 0;
			for(int j = 0; j < size; j++){
				if(deme_index[j] == i){
					subpop[counter] = new Individuals(array[j]);
					counter++;
				}
			}
			Species spec = new Species (subpop,gen_index);
			submap.put(gen_num, spec);
			population.put(deme_val, submap);
		}*/
		while(it.hasNext()){															//iterating over the species hashmap
			
			Integer deme_num = it.next();												//deme index
			
			HashMap<Integer,Demes> submap2 = population.get(deme_num);				//getting the deme specified by the deme index
			
			Individuals[] subarray = new Individuals[this.pop_size];					//a subarray used to stored transferred individuals
			
			int counter = 0;
			while(counter < this.pop_size){												//
				Random rand = new Random ();
				int int_val = rand.nextInt(size);
				Integer i_val = new Integer (int_val);
				if(!surveillieur.contains(i_val)){
					subarray[counter] = new Individuals (array[int_val]);
					surveillieur.add(i_val);
					counter++;
				}
			}
			Demes deme = new Demes (subarray,gen_index);
			submap2.put(gen_num, deme);
			population.put(deme_num, submap2);
		}
		this.deme2.putAll(population);
		System.out.println("All Individuals of generation " + gen_index + " were successfully transferred to new demes.");
		//}
		//else{
		//	System.out.println("Transfer of individuals was aborted due to lack of variation.");
		//}
		}
		else{
			if(this.deme != null){
				int deme_size=this.demenumb, popsize = this.pop_size, generation = this.gen;				
				Individuals[] indiarray = new Individuals[deme_size * popsize];
				for(int i = 0; i < deme_size; i++){
					for(int j = 0; j < popsize; j++){
						indiarray[popsize*i+j] = new Individuals(deme[i][generation].getPop(j));
					}
				}
				int deme_num = 0;
				int[][] surveilleur = generatePairSet(deme_size,popsize);
				while(deme_num < deme_size){
					int counter = 0;
					Individuals[] newarray = new Individuals[popsize];
					while(counter < popsize){
						Random ran = new Random ();
						int rand_val = ran.nextInt(surveilleur.length);
						int f_val = surveilleur[rand_val][0], s_val = surveilleur[rand_val][1];
						newarray[counter] = new Individuals (indiarray[popsize*f_val+s_val]);
						surveilleur = removeEntrie(surveilleur,rand_val);
						counter++;
					}
					this.setDeme(deme_num, gen_index, new Demes (newarray,gen_index));
					deme_num++;
					
				}
				/*
				while(deme_num < deme_size){
					int counter = 0;
					Individuals[] newarray = new Individuals[popsize];
					while(counter < popsize){
						Random rand1 = new Random ();
						Random rand2 = new Random ();
						int f_val = rand1.nextInt(deme_size),sec_val=rand2.nextInt(popsize);
						Integer first = new Integer (f_val), secon = new Integer (sec_val);
						Pair<Integer> pair = new Pair<Integer> (first,secon);
						if(!surveilleur.contains(pair)){
							newarray[counter] = new Individuals(indiarray[popsize*f_val+sec_val]);
							counter++;
							surveilleur.add(pair);
						}
					}
				}*/
			}
			else{
				throw new NullPointerException ("Neither of the two fields 'deme' nor deme2' was initialized before calling this method!");
			}
		}
	}
	
	/**
	 * method, returning the deme of this <code> Species </code> as an array of an array of Population
	 * instances.
	 * @return this Population instance's demes as a matrix copy of Population instances
	 * @throws PdbLoadError 
	 * @throws PdbCodeNotFoundError 
	 * @throws SQLException 
	 */
	public Demes[][] getDeme() throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim1 = this.deme.length, dim2 = this.deme[0].length;
		Demes[][] thisdemes = new Demes[dim1][dim2];
		for(int i = 0; i < dim1; i++){
			for(int j = 0; j < dim2; j++){
				thisdemes[i][j] = new Demes(this.deme[i][j]);
			}
		}
		return thisdemes;
	}
	
	/**
	 * method, returning the i-th deme at the j-th generation of this <code> Species </code>.
	 * @param i
	 * @param j
	 * @return
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public Demes getDeme(int i, int j) {//throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		if(this.deme[i][j] != null){
			return new Demes(this.deme[i][j]);
		}
		else{
			String nullp = "The Species at position i = " + i + ", j = " + j + " was not initialized.";
			throw new NullPointerException(nullp);
		}
	}
	
	public HashMap<Integer,HashMap<Integer,Demes>> getDeme2 (){
		return new HashMap<Integer,HashMap<Integer,Demes>> (this.deme2);
	}
	
	public HashMap<Integer,Demes> getDeme2 (int i) throws ArrayIndexOutOfBoundsException {
		int deme_size = this.deme2.size();
		if(i <  deme_size && 0 <= i){
			Integer deme_index = new Integer (i);
			return new HashMap<Integer,Demes> (this.deme2.get(deme_index));
		}
		else{
			if(i >= deme_size){
				String arrayind = "The index i = " + i + " is greater than the size of this deme s = " + deme_size;
				throw new ArrayIndexOutOfBoundsException (arrayind);
			}
			else{
				String arrayind = "The index i = " + i + " must be greater than zero";
				throw new ArrayIndexOutOfBoundsException (arrayind);
			}
		}
	}
	
	/**
	 * Getter: gets the i-th deme and the j-th generation of the field <code>deme2</code>.
	 * @param i deme index
	 * @param j generation index
	 * @return the Species in the i-th deme and the j-th generation
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 * @throws ArrayIndexOutOfBoundsException
	 */
	public Demes getDeme2 (int i, int j) {// throws SQLException, PdbCodeNotFoundError, PdbLoadError, ArrayIndexOutOfBoundsException {
		int deme_size = this.deme2.size();
		if((0 <= i) && (i < deme_size)){
			Integer deme_index = new Integer (i);
			HashMap<Integer,Demes> subset = this.getDeme2(deme_index);
			int gen_size = subset.size();
			if((0 <= j) && (j <= gen_size)){
				Integer gen_index  = new Integer (j);
				//Pair<Integer> pair = new Pair<Integer> (deme_index, gen_index);
				return new Demes (subset.get(gen_index),j);
			}
			else{
				if(j > gen_size){
					String arrayind = "The index j = " + j + " is greater than the number of generations s = " + deme_size;
					throw new ArrayIndexOutOfBoundsException (arrayind);
				}
				else{
					String arrayind = "The index j = " + j + " must be greater than zero";
					throw new ArrayIndexOutOfBoundsException (arrayind);
				}
			}
		}
		else{
			if(i >= deme_size){
				String arrayind = "The index i = " + i + " is greater than the size of this deme s = " + deme_size;
				throw new ArrayIndexOutOfBoundsException (arrayind);
			}
			else{
				String arrayind = "The index i = " + i + " must be greater than zero";
				throw new ArrayIndexOutOfBoundsException (arrayind);
			}
		}
	}
	
	/**
	 * method returing the field <code> gen </code> of this <code> Species </code>. Since this field
	 * is not used at all, it might as well be deleted.
	 * @return
	 */
	public int getGen(){
		return this.gen;
	}
	
	/**
	 * method returing the field <code> demnumb </code> of this <code> Species </code>.
	 * @return
	 */
	public int getDemeSize(){
		return this.demenumb;
	}
	
	/**
	 * method returing the field <code> evosteps </code> of this <code> Species </code>.
	 * @return
	 */
	public int getEvoSteps(){
		return this.evosteps;
	}
	
	/**
	 * method returing the field <code> evotype </code> of this <code> Species </code>.
	 * @return
	 */
	public boolean getEvoType(){
		return this.evoType;
	}
	
	/**
	 * method writing the currently evolved deme to a temporary directory in case of a
	 * severe error occurring, interrupting the run.
	 * @param k - the k-th deme to be written to a contact map file
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public void printDemesToTempFiles (int deme_num, int gen_num) throws SQLException, PdbCodeNotFoundError, PdbLoadError, FileNotFoundException, IOException{
		String pdbcode = "";
		if(this.deme[deme_num][gen_num] != null){
			pdbcode = pdbcode + this.deme[deme_num][gen_num].getPop(0).getName();
			File dir_test = new File(path + pdbcode + "/tempfiles/temp" + deme_num +"/");
			if(!dir_test.exists()){
				dir_test.mkdirs();
			}
			Demes buffer = new Demes(this.deme[deme_num][gen_num],gen_num);
			buffer.printAllIndisToTempDir(path,pdbcode + "/tempfiles/temp" + deme_num +"/");
			buffer.printToFile(path,pdbcode + "/temp" + deme_num +"/" + pdbcode + "deme" + deme_num);
		}
		else{
			if(this.getDeme2(deme_num,gen_num) != null){
				pdbcode = pdbcode + this.deme[deme_num][gen_num].getPop(0).getName();
				File dir_test = new File(path + pdbcode + "/tempfiles/temp" + deme_num +"/");
				if(!dir_test.exists()){
					dir_test.mkdirs();
				}
				Demes buffer = new Demes(this.deme[deme_num][gen_num],gen_num);
				buffer.printAllIndisToTempDir(path,pdbcode + "/tempfiles/temp" + deme_num +"/");
				buffer.printToFile(path,pdbcode + "/temp" + deme_num +"/" + pdbcode + "deme" + deme_num);
			}
			else{
				String nullp = "Neither the field 'deme' nor the field 'deme2' are intialized.";
				throw new NullPointerException(nullp);
			}
		}
	}
	
	public static void initializeIndiArray (Individuals[] array){
		if(Demes.hasNull(array)){
			for(int i = 0; i < array.length; i++){
				array[i] = new Individuals();
			}
		}
	}
	
	public static int[][] removeEntrie (int[][] array, int entry_index){
		int dim1 = array.length, dim2 = array[0].length;
		int[][] copy = new int[dim1-1][dim2];
		for(int i = 0; i < dim1; i++){
			if(i < entry_index){
				System.arraycopy(array[i], 0, copy[i], 0, dim2);
			}
			if(i > entry_index){
				System.arraycopy(array[i], 0, copy[i-1], 0, dim2);
			}
		}
		return copy;
	}
	
	/**
	 * method returing the field <code> path </code> of this <code> Species </code>.
	 * @return
	 */
	public static String getPath(){
		return new String(path);
	}
	
	public static final void setPath(String dir){
		path = new String(dir);
	}
	
	public static Integer[] convertHashMap (HashMap<Integer,Demes> map){
		int size = map.size(), counter = 0;
		Integer[] array = new Integer[size];
		Set<Integer> keyset = map.keySet();
		Iterator<Integer> it = keyset.iterator();
		while(it.hasNext()){
			Integer int_val = it.next();
			array[counter] = int_val;
			counter++;
		}
		Arrays.sort(array);
		return array;
	}
	
	public static Individuals[] convertHashMapToIndiArray (HashMap<Integer, Demes> map, int popsize){
		int size = map.size(), counter = 0;
		Individuals[] array = new Individuals[size*popsize];
		Set<Integer> keyset = map.keySet();
		Iterator<Integer> it = keyset.iterator();
		while(it.hasNext()){
			Integer int_val = it.next();
			Individuals[] subarray = map.get(int_val).getPop();
			int length = subarray.length;
			System.arraycopy(subarray, 0, array, counter, length);
			counter += length;
		}
		return array;
	}
	
	public static int[][] generatePairSet(int seed1, int seed2){
		if( seed1 > 0 && seed2 > 0){
			int[][] pairset = new int[seed1 * seed2][2]; 
			for(int i = 0; i < seed1; i++){
				for(int j = 0; j < seed2; j++){
					pairset[i * seed1 + seed2][0] = i;
					pairset[i * seed1 + seed2][1] = j;
				}
			}
			return pairset;
		}
		else{
			throw new IllegalArgumentException ("Neither of the seed parameter must be less or equal zero!");
		}
	}
	
	public static HashMap<Integer,HashSet<Pair<Integer>>> generateRandomSubSets (int seed1, int seed2){
		if(seed1 > 0 && seed2 > 0){
	//		HashSet<Integer> pairset = new HashSet<Integer> (seed2 * 2);
			HashMap<Integer,HashSet<Pair<Integer>>> randSets = new HashMap<Integer,HashSet<Pair<Integer>>> (2*seed1);
			return randSets;
		}
		else{
			throw new IllegalArgumentException ("Neither of the seed parameter must be less or equal zero!");
		}
	}

	public static void main (String[] args) throws SQLException, IOException, PdbCodeNotFoundError, PdbLoadError{
		//Population pop1 = new Population (path + "1e0l/");
		/*int[] array = {2};//{0,2,3,7,8,9,10,11,12,13,14};
		String helper = "1,4,5,6";
		//int[] num_conts = {2,5};
		for(int j = 0; j < array.length; j++){
			String str = (new Integer (array[j])).toString();
			if(helper.contains(str)){
				Population ec = new Population(array[j]);//,num_conts[0],"dummy");		
				ec.evolve("dummy");
			}
			else{
				Population ec = new Population(array[j]);//,num_conts[1],"dummy");		
				ec.evolve("dummy");
			}
			//for(int i = 0; i < ec.getDeme().length; i++){
			//ContactMaps.printToFile(ec.getDeme()[i],getPath());
			//}
		}
		/*String dir = "/home/gmueller/Arbeiten/ContactMaps/tests/main_run/1bkr/tempfiles/";
		Population p = new Population(dir);
		p.evolve();
		for(int i = 0; i < p.getPop().length; i++){
			ContactMaps.printToFile(p.getPop()[i],path);
		}*/
		//String dir = "";
		//Individuals starter = new Individuals (dir);
		Species sp = new Species (2,25,0.0003,true);
		sp.evolve2(20);
	}
}