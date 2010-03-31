package owl.embed.contactmaps;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import owl.core.structure.PdbCodeNotFoundError;
import owl.core.structure.PdbLoadError;
import owl.core.util.RegexFileFilter;

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
	
	private static final String dm = "deme";
	private static final String gn  = "gen";
	
	/*-----------------------------------fields-----------------------------------------*/
	
	/**
	 * field - deme, a matrix of {@link Demes} instances
	 */
	//public Demes[][] deme;
	
	/**
	 * int fields: demenumb - number of demes; evosteps - number of evolution steps;
	 * gen - number of generation
	 */
	private int demenumb, evosteps, gen, pop_size, numofConts;
	
	/**
	 * field: a threshold value used to discriminate fast convergence
	 */
	private static double threshold;
	
	/**
	 * field: determines whether 'CMError' (evoType = true) or 'DMError' (evoType = false)
	 * are used for evolution
	 */
	private boolean evoType;
	
	private boolean homogenous;
	
	/**
	 * field: all contacts with a frequency higher than a given threshold are saved here
	 */
	private HashMap<Pair<Integer>,Double> contacts_hfreq;
	//public String character;
	
	private HashMap<Pair<String>,Demes> deme2;
	
	/**
	 * static field: specifies the directory to which all files will be written
	 */
	private static String path = "/project/StruPPi/gabriel/Arbeiten/";
	
	//private String temp_path;	
	
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
		setDeme2(5, protentry, 14, 5, dummy,true);
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
		evosteps = 30;
		setDeme2(5, protentry, 20, numofconts, "dummy",evo);
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
		setDeme2(5,protentry, 20, num_of_conts, dummy,evo);
	}	
	
	public Species (String dir) throws FileNotFoundException, IOException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		readStarterFromFile(dir);
	}
	
	public Species (Demes[] deme, int steps, String dir, String temp_dir, double thres){
		setDeme2(deme,steps,dir,temp_dir,thres);
	}
	
	public Species (String pdb, int pop, int deme, int steps, double cont_percent, String output) throws IllegalArgumentException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		this(pdb,pop,deme,steps,cont_percent,output,true);
	}
	public Species (String pdb, int pop, int deme, int steps, double cont_percent, String output, boolean evoType) throws IllegalArgumentException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		setDeme2(pdb,pop,deme,steps,cont_percent,output);
		this.evoType = evoType;
	}
	
	/*------------------------------------------------------setter--------------------------------------------------------------*/
	
		
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
			deme2 = new HashMap<Pair<String>,Demes> (deme_size_mult);
			Demes[] starter_pop = new Demes[demesize];
			for(int i = 0; i < demesize; i++){
				Integer deme_index = new Integer (i);
				Integer gen_index = new Integer (0);
				starter_pop[i] = new Demes(popsize, protentry, num_of_conts, 0, "deme" + i);
				Pair<String> pair = new Pair<String> (dm+deme_index,gn+gen_index);
				deme2.put(pair,starter_pop[i]);
				System.out.println(starter_pop[i].getPop(0).getName() + ", deme " + i + " initialized...");
			}
			gen = 0;
			pop_size = popsize;
			demenumb = demesize;
			numofConts = starter_pop[0].getNumOfContacts();
			evoType = evo;
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
		demenumb = demesize;
		evosteps = evsteps;
		deme2 = new HashMap<Pair<String>,Demes> (2*demesize);
		for(int i = 0; i < demesize; i++){
			Integer gen_index  = new Integer (0);
			Integer deme_index = new Integer (i);
			Demes deme = new Demes (popsize, protentry, 5, 0, "deme" + i);
			Pair<String> pair = new Pair<String>(dm+deme_index,gn+gen_index);
			deme2.put(pair, deme);
		}
		pop_size = popsize;
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
			evosteps = 20 - gen;															//number of evolution steps are set to default: 20 - gen (generation, where stopped)
			deme2 = new HashMap<Pair<String>,Demes>();
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
									gen = (int) dnum;															//field: gen - number of generation at which the evolution
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
					Pair<String> pair = new Pair<String> (dm+k,gn+0);
					deme2.put(pair,new Demes(starterindis,0));															//starter are initialized

				}
				numofConts = deme2.get(new Pair<String>(dm+0,gn+0)).getPop(0).getNumOfContacts();
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
	
	public void setDeme2 (Demes[] deme, int steps, String dir, String temp_dir, double thres){
		int length = deme.length;
		deme2 = new HashMap<Pair<String>,Demes> ();
		for(int i = 0; i < length; i++){
			Pair<String> pair = new Pair<String> (dm+i,gn+0);
			deme2.put(pair, deme[i]);
		}
		demenumb   = length;
		evosteps   = steps;
		evoType    = true;
		gen        = 0;
		pop_size   = deme[0].getSize();
		numofConts = deme[0].getNumOfContacts();
		path = new String (dir);
		//temp_path = new String(temp_dir);
		threshold = thres;		
	}
	
	public void setDeme2 (String pdb, int pop, int deme, int steps, double num_conts, String out_dir) throws IllegalArgumentException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		demenumb = deme;
		pop_size = pop;
		evoType = true;
		gen = 0;
		path = new String (out_dir);
		threshold = 0.003;
		deme2 = new HashMap<Pair<String>,Demes>();
		for(int i = 0; i < deme; i++){
			Demes dms = new Demes (pdb,pop,num_conts);
			if(i==0) numofConts = dms.getPop(0).getNumOfContacts();
			Pair<String> pair = new Pair<String>(dm+i,gn+0);
			deme2.put(pair,dms);
			System.out.println(pdb + " deme "+i+" initialized...");
		}
	}
	
	public void clear (){
		deme2 = null;
		pop_size = 0; numofConts = 0; gen = 0;
		contacts_hfreq = null;
		Individuals.clearFullCMandDM();
	}

	
	/*--------------------------------------------------evolution-----------------------------------------------------------*/

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
		HashMap<Pair<String>,Demes> map = getDeme2();
		Set<Pair<String>> demeset = map.keySet();
		Iterator<Pair<String>> it = demeset.iterator();
		if(threshold == 0.0){
			threshold = 0.003;
		}
		int gennum = getGen();
		//if(gennum >= 10){
			boolean isless = isGenXHomogenous(gennum);//isLessThanThreshold(gennum, false, threshold);
			if(isless){
				swapSpecies(gennum);
				map = getDeme2();
				demeset = map.keySet();
				it = demeset.iterator();
			}
		//}
		while(it.hasNext()){
			Pair<String> pair = it.next();
			int deme_nu = (int) Double.parseDouble(pair.getFirst().replaceAll(dm, ""));
			evolve2(deme_nu,gennum);
		}
		gen += 1;
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
		int currentgen = gen;
		evosteps = final_generation + 1;
		//boolean homogenity = false;
		//if(currentgen >= 20){
			//homogenity = isGenXHomogenous(currentgen);
		//}
		boolean homogeneity = isGenXHomogenous(currentgen);
		if(currentgen >= 10){
			while(currentgen < final_generation && !homogeneity){
				homogeneity = isGenXHomogenous(currentgen);
				evolve2();
				currentgen = gen;
			}
		}
		else{
			if(!homogeneity){
				if(currentgen == 0){
					for(int i = 0; i < demenumb; i++){
						Demes spec = getDeme2(i, 0);
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
					checkForHomogeneity(currentgen,0.8);
					if(homogenous){
						break;
					}
					evolve2();
					currentgen = gen;
				}
				if(homogenous){
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
		HashMap<Integer,Demes> submap    = getDeme2(deme_index.intValue());
		Pair<String> npair = new Pair<String> (dm+deme_index,gn+gen_plus_one);
		if(!submap.containsKey(gen_plus_one)){
			boolean evo = evoType;
			Demes off = new Demes (submap.get(new Integer(current_gen)).evolve(evo, "dummy"), gen_plus_one);
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
			deme2.put(npair, off);
			//gen = current_gen + 1;
		}
	}
	
	/**
	 * This auxiliary method chekcs whether all contact are already present in the Population and therefore
	 * any evolution will not yield any better result than the current one.
	 * @param gen_index
	 */
	public void checkForHomogeneity(int gen_index, double threshold_val){
		if(deme2 != null){
			HashMap<String,Demes> submap = getDemesOfSameGen(gen_index);
			Integer zero = new Integer (0);
			HashSet<Pair<Integer>> common_contacts = new HashSet<Pair<Integer>> (submap.get(dm+zero).getSubHash());
			submap.remove(dm+zero);
			Set<String> keyset = submap.keySet();
			Iterator<String> it = keyset.iterator();
			while(it.hasNext()){
				String deme_index = it.next();
				common_contacts.retainAll(submap.get(deme_index).getSubHash());
			}
			int frequency = (int) (((double) common_contacts.size())/threshold_val);
			if(frequency >= numofConts){
				homogenous = true;
			}
		}
		else throw new NullPointerException ("Deme field not initialized!");
	}
	
	public void printStarterGenerationToFile (String dir, String addname) throws FileNotFoundException, IOException{
		if(deme2 != null){
			HashMap<String,Demes> starterset = getDemesOfSameGen(0);
			Set<String> keyset = starterset.keySet();
			Iterator<String> it = keyset.iterator();
			int counter = 0;
			while(it.hasNext()){
				String index = it.next();
				Demes spec = starterset.get(index);
				spec.printAllIndisToTempDir(dir, addname + counter);
				counter++;
			}
		}
		else throw new NullPointerException ("Deme field not initialized!");
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
		Individuals ne = new Individuals(getDeme2(0,0).getPop(0));
		String cont = "#CMVIEW GRAPH FILE ver: 1.0\n#SEQUENCE: "+ne.getSequence()+"\n"+
		"#PDB: "+ne.getName()+ "\n#PDB CHAIN CODE: "+ne.getChainCode()+"\n#CT: "+Individuals.getContactT()+ "\n#CUTOFF: "+Individuals.getContactDist()+"\n"+ 
		"#GENERATION: "+0 + "\n#Species SIZE: " + getDeme2(0,0).getSize() + "\n#CMError: "+ 0 + "\n#CMstDev: " + 0 + "\n#DMError: "+ 0 + 
		"\n#DMstDev: " + 0 + "\n#NUMB. CONTACTS: " + ne.getNumOfContacts() +"\n" + "#NUMB. OF ALL CONTACTS: " + ne.getFullContact() + "\n";
		HashSet<Pair<Integer>> keyset = new HashSet<Pair<Integer>> (contacts_hfreq.keySet());
		int[][] index_array = SortIntArray.converter(keyset);
		int length = index_array[0].length;
		//TreeSet<CompPair> trset = new TreeSet<CompPair> (trmap.keySet());
		HashMap<Pair<Integer>,Double> hashmap = new HashMap<Pair<Integer>,Double> (contacts_hfreq);
		//HashSet<Pair<Integer>> hashset = getKey();
		//Iterator<CompPair> it = trset.iterator();
		cont = cont + "\n";
		//while(it.hasNext()){
		for(int i = 0; i < length; i++){
			int index1 = index_array[0][i], index2 = index_array[1][i];
			cont = cont + index1 + "\t" + index2 + "\t" + hashmap.get(new Pair<Integer>(new Integer(index1),new Integer(index2))) + "\n"; 
		}
		FileOutputStream file = new FileOutputStream(path+name+"-cons"+getDeme2(0,0).getNumOfContacts()+".cmap");
		PrintStream printa = new PrintStream(file);
		printa.print(cont);
		printa.close();
		file.close();
		System.out.println(ne.getName()+" at generation "+gen+" written to file...");
	}
	
	public Demes[] copyPops(Demes[] input) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim1 = input.length;
		Demes[] pop = new Demes[dim1];
		for(int i = 0; i < dim1; i++){
			pop[i] = new Demes(input[i]);
		}
		return pop;
	}
	
	public int getLargestSubHashMap (){
		int size = 0;
		for(int i = 0; i < demenumb; i++){
			HashMap<Integer,Demes> map = getDeme2(i);
			if(map.size() >= size) size = map.size();
		}
		return size;
	}
	
	public HashMap<String,Demes> getDemesOfSameGen (int gen_index) {//throws ArrayIndexOutOfBoundsException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		HashMap<Pair<String>,Demes> deme_map = getDeme2();
		HashMap<String,Demes> map = new HashMap<String,Demes> ();
		Set<Pair<String>> keyset = deme_map.keySet();
		Iterator<Pair<String>> it = keyset.iterator();
		while(it.hasNext()){
			Pair<String> pair = it.next();
			if(pair.getSecond().matches(gn+gen_index)){
				String deme_index = pair.getFirst();
				Demes spec = deme_map.get(pair);
				map.put(deme_index, spec);
			}
		}
		return map;
	}
	
	public double[] getErrorDev (int gen_index, boolean CMDM) throws ArrayIndexOutOfBoundsException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		HashMap<String,Demes> map = getDemesOfSameGen(gen_index);
		Set<String> keyset = map.keySet();
		double[] stdev_array = new double[keyset.size()];
		Iterator<String> it = keyset.iterator();
		if(CMDM){
			while(it.hasNext()){
				String index = it.next();
				int index_val = (int) Double.parseDouble(index.replaceAll(dm, ""));
				stdev_array[index_val] = map.get(index).getCMstdev();
			}
		}
		else{
			while(it.hasNext()){
				String index = it.next();
				int index_val = (int) Double.parseDouble(index.replaceAll(dm, ""));
				stdev_array[index_val] = map.get(index).getDMstdev();
			}
		}
		return stdev_array;
	}
	
	public double getMeanStDev(int gen_index, boolean CMDM) throws ArrayIndexOutOfBoundsException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		double[] st_dev = getErrorDev(gen_index, CMDM);
		int dim =  st_dev.length;
		double sum = 0.0;
		for(int i = 0; i < dim; i++){
			sum += st_dev[i];
		}
		return sum/((double) dim);
	}
	
	public boolean isLessThanThreshold (int gen_index, boolean CMDM, double threshold) throws ArrayIndexOutOfBoundsException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		double average = getMeanStDev(gen_index, CMDM);
		if(average <= threshold){
			return true;
		}
		else{
			return false;
		}
	}
	
	/*public void isPopHomogenous (int gen) throws ArrayIndexOutOfBoundsException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		HashMap<Integer,Species> pop_map = getDemesOfSameGen(gen);
		Set<Integer> keyset = pop_map.keySet();
		if(is_homogen == null){
			is_homogen = new HashMap<Integer,Integer> ();
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
			is_homogen.put(gen_index,new Integer (1));
		}
		else{
			is_homogen.put(gen_index, new Integer (0));
		}
	}*/
	
	
	
	public boolean isGenXHomogenous (int gen) throws ArrayIndexOutOfBoundsException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		if(deme2 != null){
			HashMap<String,Demes> deme_map = getDemesOfSameGen(gen);
			Integer deme0 = new Integer (0);
			HashSet<Pair<Integer>> pairs = new HashSet<Pair<Integer>> (deme_map.get(dm+deme0).getSubHash());
			int sum = pairs.size();
			deme_map.remove(dm+deme0);
			Set<String> keyset = deme_map.keySet();
			Iterator<String> it = keyset.iterator();
			while(it.hasNext()){
				String deme = it.next();
				HashSet<Pair<Integer>> subhash = new HashSet<Pair<Integer>> (deme_map.get(deme).getSubHash());
				if(subhash.size() != 0){
					sum += subhash.size();
				}
			}
			int benchmark = (int) (((double) numofConts)*0.8);
			int av_frequency = (int) (((double) sum)/((double) demenumb));
			if(av_frequency >= benchmark){
				System.out.println("Population is homogenous...");
				return true;
			}
			else{
				System.out.println("Population is inhomogenous...");
				return false;
			}
		}
		else throw new NullPointerException ("Deme field not initialized!");
	}
	
	public void swapSpecies (int gen_index) throws ArrayIndexOutOfBoundsException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		if(deme2 != null){
		HashMap<Pair<String>,Demes> population = getDeme2();			//getting the whole population
		
		Integer gen_num = new Integer (gen_index);										//generation index
		
		HashMap<String,Demes> submap1 = getDemesOfSameGen(gen_index);			//getting the all species at the specified generation
		
		HashSet<Integer> surveillieur = new HashSet<Integer> ();						//checker, checks that no index is taken twice
		
		Set<String> keyset = submap1.keySet();											//keyset to iterate over the species of the specified generation
		
		Iterator<String> it = keyset.iterator();
		Individuals[] array = convertHashMapToIndiArray(submap1, pop_size);		//converting the species hashmap to an individuals array
		
		int size = array.length;//, deme_size = demenumb;
		/*int[] deme_index = new int[size];
		for(int i = 0; i < size; i++){
			Random rand = new Random ();
			deme_index[i] = rand.nextInt(deme_size);
		}
		for(int i = 0; i < deme_size; i++){
			Integer deme_val = new Integer (i);
			HashMap<Integer,Species> submap = population.get(deme_val);
			Individuals[] subpop = new Individuals[pop_size];
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
			
			String deme_num = it.next();												//deme index
			
			Pair<String> pair = new Pair<String>(deme_num,gn+gen_num);
			
			//Demes submap2 = submap1.get(deme_num);				//getting the deme specified by the deme index
			
			Individuals[] subarray = new Individuals[pop_size];					//a subarray used to store transferred individuals
			
			int counter = 0;
			while(counter < pop_size){												//
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
			population.put(pair, deme);
		}
		deme2.putAll(population);
		System.out.println("All Individuals of generation " + gen_index + " were successfully transferred to new demes.");
		//}
		//else{
		//	System.out.println("Transfer of individuals was aborted due to lack of variation.");
		//}
		}
		else throw new NullPointerException ("Deme field not initialized!");
	}
	
	public HashMap<Pair<String>,Demes> getDeme2 (){
		return new HashMap<Pair<String>,Demes> (deme2);
	}
	
	public HashMap<Integer,Demes> getDeme2 (int i) throws ArrayIndexOutOfBoundsException {
		int deme_size = deme2.size();
		if(i <  deme_size && 0 <= i){
			Set<Pair<String>> keys = getDeme2().keySet();
			Iterator<Pair<String>> it = keys.iterator();
			HashMap<Integer,Demes> gen_map = new HashMap<Integer,Demes>();
			while(it.hasNext()){
				Pair<String> pair = it.next();
				int val1 = (int) Double.parseDouble(pair.getFirst().replaceAll(dm, ""));
				int val2 = (int) Double.parseDouble(pair.getSecond().replaceAll(gn, ""));
				if(val1 == i) gen_map.put(new Integer(val2), getDeme2().get(pair));
			}
			return gen_map;
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
		Pair<String> pair = new Pair<String> (dm+i,gn+j);
		if(deme2.containsKey(pair)) return new Demes (deme2.get(pair));
		else throw new IllegalArgumentException ("No such entry!");
	}
	
	/**
	 * method returing the field <code> gen </code> of this <code> Species </code>. Since this field
	 * is not used at all, it might as well be deleted.
	 * @return
	 */
	public int getGen(){
		return gen;
	}
	
	/**
	 * method returing the field <code> demnumb </code> of this <code> Species </code>.
	 * @return
	 */
	public int getDemeSize(){
		return demenumb;
	}
	
	/**
	 * method returing the field <code> evosteps </code> of this <code> Species </code>.
	 * @return
	 */
	public int getEvoSteps(){
		return evosteps;
	}
	
	/**
	 * method returing the field <code> evotype </code> of this <code> Species </code>.
	 * @return
	 */
	public boolean getEvoType(){
		return evoType;
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
			if(getDeme2(deme_num,gen_num) != null){
				Demes buffer = getDeme2(deme_num,gen_num);
				pdbcode = pdbcode + buffer.getPop(0).getName();
				File dir_test = new File(path + pdbcode + "/tempfiles/temp" + deme_num +"/");
				if(!dir_test.exists()){
					dir_test.mkdirs();
				}
				//buffer.printAllIndisToTempDir(path,pdbcode + "/tempfiles/temp" + deme_num +"/");
				buffer.printToFile(path,pdbcode + "/temp" + deme_num +"/" + pdbcode + "deme" + deme_num);
			}
			else{
				String nullp = "Neither the field 'deme' nor the field 'deme2' are intialized.";
				throw new NullPointerException(nullp);
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
	
	public static Individuals[] convertHashMapToIndiArray (HashMap<String, Demes> map, int popsize){
		int size = map.size(), counter = 0;
		Individuals[] array = new Individuals[size*popsize];
		Set<String> keyset = map.keySet();
		Iterator<String> it = keyset.iterator();
		while(it.hasNext()){
			String int_val = it.next();
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
		String[] prots = {"1bkr","1e0l","1e6k","1o8w","1odd","1onl","1pzc","1r9h","1sha","1ugm","2ci2"};
		for(int i = 0; i < prots.length; i++){
			Species sp = new Species(prots[i+1],20,5,20,10.0,"/project/StruPPi/gabriel/Arbeiten/");
			sp.evolve2(5);
		}
	}
}