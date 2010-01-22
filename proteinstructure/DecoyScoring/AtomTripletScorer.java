package proteinstructure.DecoyScoring;

import gnu.getopt.Getopt;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;

import proteinstructure.AIGNode;
import proteinstructure.AIGraph;
import proteinstructure.FileFormatError;
import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbasePdb;
import proteinstructure.TemplateList;

import tools.MySQLConnection;
import tools.Triplet;

public class AtomTripletScorer extends TripletScorer {
	
	private static final File DEFAULT_LIST_FILE = new File("/project/StruPPi/jose/emp_potential/cullpdb_pc20_res1.6_R0.25_d090728_chains1627.list");
	private static final double DEFAULT_CUTOFF = 4.0;
	private static final int DEFAULT_MIN_SEQ_SEP = 3;

	private static final int NUM_ATOM_TYPES = 167;
	

	/**
	 * Constructs a AtomTripletScorer by taking a list of structure ids (pdbCodes+pdbChainCodes),
	 * and parameters contact type, distance cutoff and minimum sequence separation used for 
	 * calculating the scoring matrix. 
	 * @param listFile the file with a list of PDB ids (pdbCode+pdbChainCode)
	 * @param cutoff the distance cutoff to be used as definition of contacts
	 * @param minSeqSep the minimum sequence separation for a contact to be considered
	 * @throws SQLException if can't establish connection to db server
	 */
	public AtomTripletScorer(File listFile, double cutoff, int minSeqSep) throws SQLException {
		
		this.scoringMethod = ScoringMethod.ATOMTRIPLET;
		
		this.types2indices = new HashMap<String, Integer>();
		this.indices2types = new HashMap<Integer, String>();

		this.structureIds = new ArrayList<String>();
		this.listFile = listFile;
		this.cutoff = cutoff;
		this.minSeqSep = minSeqSep;
		
		this.numEntities = NUM_ATOM_TYPES;
		entityCounts = new int[numEntities];
		tripletCounts = new int[numEntities][numEntities][numEntities];
		
		this.conn = new MySQLConnection();
		
	}

	/**
	 * Constructs a AtomTripletScorer given a file with a scoring matrix. All values of the matrix 
	 * and some other necessary fields will be read from the file.
	 * @param scMatFile
	 * @throws IOException
	 * @throws FileFormatError
	 */
	public AtomTripletScorer(File scMatFile) throws IOException, FileFormatError  {
		this.scoringMethod = ScoringMethod.ATOMTRIPLET;
		readScMatFromFile(scMatFile);
	}

	
	@Override
	public void countTriplets() throws SQLException, IOException {
		Scorer.initAtomMap(types2indices, indices2types);

		for (String id:TemplateList.readIdsListFile(listFile)) {
			String pdbCode = id.substring(0,4);
			String pdbChainCode = id.substring(4,5);
			Pdb pdb = null;
			try {
				pdb = new PdbasePdb(pdbCode,DB,conn);
				pdb.load(pdbChainCode);
				if (!isValidPdb(pdb)) {
					System.err.println(id+" didn't pass the quality checks to be included in training set");
					continue;
				}
				System.out.println(id);
				this.structureIds.add(id);
			} catch (PdbCodeNotFoundError e) {
				System.err.println("Couldn't find pdb "+pdbCode);
				continue;
			} catch (PdbLoadError e) {
				System.err.println("Couldn't load pdb "+pdbCode);
				continue;
			}

			
			AIGraph graph = pdb.getAllAtomGraph(cutoff);
			graph.restrictContactsToMinRange(minSeqSep);
			for (AIGNode node:graph.getVertices()) {
				// we have to avoid OXT atoms, they are not in the map
				if (node.getAtomName().equals("OXT")) continue;
				countEntity(types2indices.get(node.getParent().getResidueType()+node.getAtomName()));
			}
			//for (Triplet<AIGNode> trip:graph.getTriangles()) {
			for (Triplet<AIGNode> trip:graph.getTriplets()) {
				AIGNode inode = trip.getFirst();
				AIGNode jnode = trip.getSecond();
				AIGNode knode = trip.getThird();
				if (inode.getAtomName().equals("OXT") || 
						jnode.getAtomName().equals("OXT") || 
						knode.getAtomName().equals("OXT")) continue;
				countTriplet(types2indices.get(inode.getParent().getResidueType()+inode.getAtomName()),
						types2indices.get(jnode.getParent().getResidueType()+jnode.getAtomName()),
						types2indices.get(knode.getParent().getResidueType()+knode.getAtomName()));
			}
		}
		this.totalStructures = structureIds.size();
	}
	
	@Override
	public double scoreIt(Pdb pdb) {
		AIGraph graph = pdb.getAllAtomGraph(this.cutoff);
		graph.restrictContactsToMinRange(minSeqSep);
		return scoreIt(graph);
	}

	protected double scoreIt(AIGraph graph) {
		double totalScore = 0;
		//Collection<Triplet<AIGNode>> triplets = graph.getTriangles();
		Collection<Triplet<AIGNode>> triplets = graph.getTriplets();
 		for (Triplet<AIGNode> trip:triplets) {
			AIGNode inode = trip.getFirst();
			AIGNode jnode = trip.getSecond();
			AIGNode knode = trip.getSecond();
			if (inode.getAtomName().equals("OXT") || jnode.getAtomName().equals("OXT")) continue;
			int[] ind = {types2indices.get(inode.getParent().getResidueType()+inode.getAtomName()),
					types2indices.get(jnode.getParent().getResidueType()+jnode.getAtomName()),
					types2indices.get(knode.getParent().getResidueType()+knode.getAtomName())};
			Arrays.sort(ind);
			totalScore+= scoringMat[ind[0]][ind[1]][ind[2]];
		}
		
		return (totalScore/triplets.size());
	}
	
	
	public static void main(String[] args) throws Exception {
		String help = 
			"\nCompiles a scoring matrix based on atom types from a given file with a list of \n" +
			"pdb structures.\n" +
			"Usage:\n" +
			"AtomTripletScorer -o <output_matrix_file> [-l <list_file> -c <cutoff> -m <min_seq_sep>]\n"+
			"  -o <file>     : file to write the scoring matrix to\n" +
			"  -l <file>     : file with list of pdbCodes+pdbChainCodes to use as training set \n" +
			"                  the scoring matrix. Default is "+DEFAULT_LIST_FILE+"\n" +
			"  -c <float>    : distance cutoff for the contacts. Default: "+DEFAULT_CUTOFF+"\n" +
			"  -m <int>      : minimum sequence separation to consider a contact. Default: "+DEFAULT_MIN_SEQ_SEP+"\n" +
			"  -w            : write counts matrix as well as scoring matrix \n";

			File listFile = DEFAULT_LIST_FILE;
			File scMatFile = null;
			double cutoff = DEFAULT_CUTOFF;
			int minSeqSep = DEFAULT_MIN_SEQ_SEP;
			boolean writeCounts = false;

			Getopt g = new Getopt("AtomTripletScorer", args, "l:c:o:m:wh?");
			int c;
			while ((c = g.getopt()) != -1) {
				switch(c){
				case 'l':
					listFile = new File(g.getOptarg());
					break;
				case 'c':
					cutoff = Double.parseDouble(g.getOptarg());
					break;				
				case 'o':
					scMatFile = new File(g.getOptarg());
					break;
				case 'm':
					minSeqSep = Integer.parseInt(g.getOptarg());
					break;
				case 'w':
					writeCounts = true;
					break;									
				case 'h':
				case '?':
					System.out.println(help);
					System.exit(0);
					break; // getopt() already printed an error
				}
			}

			if (scMatFile==null) {
				System.err.println("Missing argument output matrix file (-o)");
				System.exit(1);
			}


		AtomTripletScorer sc = new AtomTripletScorer(listFile,cutoff,minSeqSep);
		sc.countTriplets();
		sc.calcScoringMat();
		sc.writeScMatToFile(scMatFile,writeCounts);

	}

}
