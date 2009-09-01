import gnu.getopt.Getopt;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;

import proteinstructure.DecoyScoreSet;
import proteinstructure.DecoyScoreSetsGroup;
import proteinstructure.FileFormatError;
import proteinstructure.DecoyScore;
import proteinstructure.Pdb;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbfilePdb;
import proteinstructure.Scorer;
import tools.MySQLConnection;
import tools.RegexFileFilter;


public class scoreDecoys {
	
	public static final String DECOYS_BASEDIR = "/scratch/local/graphletDecoys/decoys/Scratch/dd/multiple";
	public static final String[] DECOYSETS = 
	{"4state_reduced", "fisa", "fisa_casp3", "hg_structal", "ig_structal", "ig_structal_hires", "lattice_ssfit", "lmds", "vhp_mcmd"};
	
	private static final int DEFAULT_MIN_SEQ_SEP = 3;
	
	private static final String SCORES_TABLE_NAME = "scores";
	private static final boolean DEBUG = false;
	
	public static void main(String[] args) throws IOException, FileFormatError, SQLException {

		String help = 
			"\nScores a set of decoy sets \n" +
			"Usage:\n" +
			"scoreDecoys -o <out_dir> -a <file> -r <file> [-m <min_seq_sep>]\n"+
			"  -o <dir>      : output directory where all output files will be written\n" +
			"  -s <file>     : file with scoring matrix\n" +
			"  -m <int>      : minimum sequence separation to consider a contact. Default: "+DEFAULT_MIN_SEQ_SEP+"\n" +
			"  -w <db_name>  : writes results also to given database name, table "+SCORES_TABLE_NAME+"\n";

		
		File outDir = null;
		File scMatFile = null;
		int minSeqSep = DEFAULT_MIN_SEQ_SEP;
		String dbName = null;
		MySQLConnection conn = null;

		Getopt g = new Getopt("scoreDecoys", args, "o:s:m:w:h?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'o':
				outDir = new File(g.getOptarg());
				break;
			case 's':
				scMatFile = new File(g.getOptarg());
				break;
			case 'm':
				minSeqSep = Integer.parseInt(g.getOptarg());
				break;				
			case 'w':
				dbName = g.getOptarg();
				conn = new MySQLConnection();
				break;
			case 'h':
			case '?':
				System.out.println(help);
				System.exit(0);
				break; // getopt() already printed an error
			}
		}
		
		if (outDir==null) {
			System.err.println("Must specify an output directory (-o)");
			System.err.println(help);
			System.exit(1);
		}

		if (scMatFile==null) {
			System.err.println("A scoring matrix file must be specified (-s)");
			System.err.println(help);
			System.exit(1);			
		}
		

		
		
		Scorer scorer = Scorer.readScoreMatFromFile(scMatFile);
		
		File decoysSetDir = new File(DECOYS_BASEDIR);
		int countAllDecoys = 0;
		int countValDecoys = 0;
		for (String decoySet: DECOYSETS) {
			System.out.println(decoySet);
			File decoysDir = new File(decoysSetDir,decoySet);
			File listFile = new File(decoysDir,"list");
			if (!listFile.canRead()) {
				System.err.println("Can't find list file "+listFile);
				continue;
			}
			
			DecoyScoreSetsGroup setsGroup = new DecoyScoreSetsGroup();
			
			BufferedReader br = new BufferedReader(new FileReader(listFile));
			String line;
			while ((line=br.readLine())!=null) {
				if (line.isEmpty() || line.startsWith("#")) continue;
				String decoy = line.trim();
				countAllDecoys++;
				System.out.println(decoy);
				File dir = new File(decoysDir,decoy);
				File[] files = dir.listFiles(new RegexFileFilter("^.*\\.pdb"));			

				DecoyScoreSet decoyScoreSet = new DecoyScoreSet(decoy, decoySet, scorer);
				
				File rmsdFile = new File(dir,"rmsds");
				HashMap<String,Double> rmsds = null;
				try {
					rmsds = DecoyScoreSet.readRmsds(rmsdFile);
				} catch (IOException e) {
					System.err.println("Couldn't read rmsds file "+rmsdFile+". Won't score this decoy");
					continue;
				}
								
				for (File file:files) {
					Pdb pdb = new PdbfilePdb(file.getAbsolutePath());

					try {
						String[] chains = pdb.getChains();
						if (chains==null) {
							System.out.print("x");
							if (DEBUG) System.out.println("\nFile "+file+" doesn't contain any chain");
							continue;
						}
						if (chains.length>1) System.out.println("\nWarning. More than one chain in "+file);
						pdb.load(chains[0]);
						System.out.print(".");
						
					} catch (PdbLoadError e) {
						System.out.print("x");
						if (DEBUG) System.err.println("\nCouldn't load "+file+". Error: "+e.getMessage());
						continue;
					}
					
					if (pdb.isAllAtom()) {
						double score = scorer.scoreIt(pdb,minSeqSep);

						if (rmsds.containsKey(file.getName())) {
							decoyScoreSet.addDecoyScore(new DecoyScore(file,score,rmsds.get(file.getName())));
						} else {
							System.err.println("Couldn't find rmsd value for "+file+". Skipping it.");
						}
					} else {
						System.out.println("Skipping "+file+"because it's not an all-atom structure");
					}

				}
				System.out.println();
			
				if (decoyScoreSet.containsNative()) {
					setsGroup.addDecoyScoreSet(decoyScoreSet);
					countValDecoys++;
				} else {
					System.err.println("\nCouldn't find native structure from decoy "+decoy+" of set "+decoySet+". Will exclude this decoy from stats.");
				}

				// writing scores files
				File scoreFile = new File(outDir,decoySet+"_"+decoy+".scores");
				decoyScoreSet.writeToFile(scoreFile);
				
			}
			
			br.close();
			
			// writing stats file
			File statsFile = new File(outDir,decoySet+".stats");
			setsGroup.writeStatsToFile(statsFile);
			
			// writing all scores to db if switch -w specified
			if (dbName!=null) setsGroup.writeToDb(conn, dbName, SCORES_TABLE_NAME);
			
		}
		System.out.println("Done. Total decoys "+countAllDecoys+", scored including native: "+countValDecoys);

	}
	
	
		
}
