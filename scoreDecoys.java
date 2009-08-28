import gnu.getopt.Getopt;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import proteinstructure.AtomCountScorer;
import proteinstructure.AtomTypeScorer;
import proteinstructure.DecoyScoreSet;
import proteinstructure.DecoySetStats;
import proteinstructure.FileFormatError;
import proteinstructure.DecoyScore;
import proteinstructure.Pdb;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbfilePdb;
import proteinstructure.ResCountScorer;
import proteinstructure.ResTypeScorer;
import proteinstructure.Scorer;
import tools.RegexFileFilter;


public class scoreDecoys {
	
	private static final String DECOYS_BASEDIR = "/scratch/local/decoys/ioannis/dd/multiple";
	private static final String[] DECOYSETS = 
	{"4state_reduced", "fisa", "fisa_casp3", "hg_structal", "ig_structal", "ig_structal_hires", "lattice_ssfit", "lmds", "vhp_mcmd"};
	
	private static final int DEFAULT_MIN_SEQ_SEP = 3;
	

	
	public static void main(String[] args) throws IOException, FileFormatError {

		String help = 
			"\nScores a set of decoy sets \n" +
			"Usage:\n" +
			"scoreDecoys -o <out_dir> -a <file> -r <file> [-m <min_seq_sep>]\n"+
			"  -o <dir>      : output directory where all output files will be written\n" +
			"  -a <file>     : file with atom scoring matrix\n" +
			"  -r <file>     : file with residue scoring matrix\n" +
			"  -m <int>      : minimum sequence separation to consider a contact. Default: "+DEFAULT_MIN_SEQ_SEP+"\n\n" +
			"Either an atom scoring matrix file or a residue scoring matrix file or both can be provided. If both\n" +
			"are provided then the stats file will contain the statistics for both side by side.";

		
		File outDir = null;
		File atomScMatFile = null;
		File resScMatFile = null;
		int minSeqSep = DEFAULT_MIN_SEQ_SEP;

		Getopt g = new Getopt("scoreDecoys", args, "o:a:r:m:h?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'o':
				outDir = new File(g.getOptarg());
				break;
			case 'a':
				atomScMatFile = new File(g.getOptarg());
				break;
			case 'r':
				resScMatFile = new File(g.getOptarg());
				break;				
			case 'm':
				minSeqSep = Integer.parseInt(g.getOptarg());
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

		if (atomScMatFile==null && resScMatFile==null) {
			System.err.println("At least an atom or a residue scoring matrix file must be specified");
			System.err.println(help);
			System.exit(1);			
		}
		

		
		
		Scorer atomScorer = null;
		if (atomScMatFile!=null) {
			atomScorer = Scorer.readScoreMatFromFile(atomScMatFile);
			if (!(atomScorer instanceof AtomTypeScorer) && !(atomScorer instanceof AtomCountScorer)) {
				System.err.println("Wrong score matrix file "+atomScMatFile+", was expecting an atom scoring matrix file");
				System.exit(1);
			}
		}
		Scorer resScorer = null;
		if (resScMatFile!=null) {
			resScorer = Scorer.readScoreMatFromFile(resScMatFile);
			if (!(resScorer instanceof ResTypeScorer) && !(resScorer instanceof ResCountScorer)) {
				System.err.println("Wrong score matrix file "+resScMatFile+", was expecting a residue scoring matrix file");
				System.exit(1);
			}
		}
		
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
			HashMap<String, DecoySetStats> resStats = null;
			if (resScMatFile!=null) 
				resStats = new HashMap<String, DecoySetStats>();
			HashMap<String, DecoySetStats> atomStats = null;
			if (atomScMatFile!=null)
				atomStats = new HashMap<String, DecoySetStats>();
			
			BufferedReader br = new BufferedReader(new FileReader(listFile));
			String line;
			while ((line=br.readLine())!=null) {
				if (line.isEmpty() || line.startsWith("#")) continue;
				String decoy = line.trim();
				countAllDecoys++;
				System.out.println(decoy);
				File dir = new File(decoysDir,decoy);
				File[] files = dir.listFiles(new RegexFileFilter("^.*\\.pdb"));			

				DecoyScoreSet allAtomScores = null;
				if (atomScMatFile!=null) 
					allAtomScores = new DecoyScoreSet(decoy);
				DecoyScoreSet allResScores = null;
				if (resScMatFile!=null)
					allResScores = new DecoyScoreSet(decoy);
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
					System.out.print(".");
					try {
						String[] chains = pdb.getChains();
						if (chains==null) {
							System.out.println("\nFile "+file+" doesn't contain any chain");
							continue;
						}
						if (chains.length>1) System.out.println("\nWarning. More than one chain in "+file);
						pdb.load(chains[0]);
					} catch (PdbLoadError e) {
						System.err.println("\nCouldn't load "+file+". Error: "+e.getMessage());
						continue;
					}
					
					if (pdb.isAllAtom()) {
						double atomScore = Double.NaN;
						if (atomScorer!=null) 
							atomScore = atomScorer.scoreIt(pdb,minSeqSep);
						double resScore = Double.NaN;
						if (resScorer!=null)
							resScore = resScorer.scoreIt(pdb,minSeqSep);

						if (rmsds.containsKey(file.getName())) {
							if (allAtomScores!=null) allAtomScores.addDecoyScore(new DecoyScore(file,atomScore,rmsds.get(file.getName())));
							if (allResScores!=null)  allResScores.addDecoyScore(new DecoyScore(file,resScore,rmsds.get(file.getName())));
						} else {
							System.err.println("Couldn't find rmsd value for "+file+". Skipping it.");
						}
					} else {
						System.out.println("Skipping "+file+"because it's not an all-atom structure");
					}

				}
				System.out.println();
			
				if ((allResScores!=null && allResScores.containsDecoyScore(decoy+".pdb")) ||
						(allAtomScores!=null && allAtomScores.containsDecoyScore(decoy+".pdb"))) {
					String nativeFileName = decoy+".pdb";
					if (resStats!=null)  resStats.put(decoy, new DecoySetStats(allResScores,nativeFileName));
					if (atomStats!=null) atomStats.put(decoy, new DecoySetStats(allAtomScores,nativeFileName));
					countValDecoys++;
				} else {
					System.err.println("\nCouldn't find native structure from decoy "+decoy+" of set "+decoySet+". Will exclude this decoy from stats.");
				}

				// writing scores files
				File resScoreFile  = new File(outDir,decoySet+"_"+decoy+".res.scores");
				File atomScoreFile = new File(outDir,decoySet+"_"+decoy+".atom.scores");
				if (allResScores!=null)  allResScores.writeToFile(resScoreFile);
				if (allAtomScores!=null) allAtomScores.writeToFile(atomScoreFile);
				
			}
			
			br.close();
			
			// writing stats file
			File statsFile = new File(outDir,decoySet+".stats");
			if (resStats!=null && atomStats!=null){
				DecoySetStats.writeStats(statsFile, resStats, atomStats);
			} else if (resStats!=null){
				DecoySetStats.writeStats(statsFile, resStats);
			} else if (atomStats!=null) {
				DecoySetStats.writeStats(statsFile, atomStats);
			}
			
		}
		System.out.println("Done. Total decoys "+countAllDecoys+", scored including native: "+countValDecoys);

	}
		
}
