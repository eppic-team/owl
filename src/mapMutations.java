import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.*;

import owl.core.connections.CSAConnection;
import owl.core.runners.DsspRunner;
import owl.core.runners.NaccessRunner;
import owl.core.structure.AminoAcid;
import owl.core.structure.Pdb;
import owl.core.structure.PdbCodeNotFoundError;
import owl.core.structure.PdbLoadError;
import owl.core.structure.PdbasePdb;
import owl.core.structure.PdbfilePdb;
import owl.core.structure.features.CatalSiteSet;
import owl.core.structure.features.CatalyticSite;
import owl.core.structure.features.SecStrucElement;
import owl.core.structure.graphs.RIGEdge;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;

import edu.uci.ics.jung.algorithms.shortestpath.Distance;
import edu.uci.ics.jung.algorithms.shortestpath.UnweightedShortestPath;


/**
 * Given a structure and a list of mutations, visualize the mutated sites on the structure
 * and calculate their surface accessibility and secondary structure state.
 * 
 * Todo:
 * Evaluate graph properties of the mutated site: distance from center, interface or active site
 * @author stehr
 */
public class mapMutations {

	private static final String NACCESS_EXECUTABLE = "/project/StruPPi/bin/naccess";
	private static final String NACCESS_PARAMETERS = "";
	private static final String DSSP_EXECUTABLE = "/project/StruPPi/Software/dssp/dsspcmbi";
	private static final String DSSP_PARAMETERS = "--";
	private static final double EXPOSURE_CUTOFF = 5.0; // everything above this cutoff is considered exposed
	private static final String PYMOL_EXECUTABLE = "/project/StruPPi/bin/pymol -c";
	
	private static final boolean keepPymolScript = true;	// if false, script will be deleted on exit
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if(args.length < 2) {
			System.out.println(mapMutations.class.getName() + " <pdbcode or filename> <mutated positions>");
			System.exit(1);
		}
		
		// read structure
		Pdb pdb = readStructureOrExit(args[0]);
		File outFile = new File(args[0] + ".png");
		
		// read mutated positions
		int[] mutations = new int[args.length-1];
		for (int i = 1; i < args.length; i++) {
			int pos = Integer.parseInt(args[i]);
			if(pos > pdb.getObsLength()) {
				System.err.println("Error: Position " + pos + " is bigger than length of protein (" + pdb.getObsLength() + "). Skipping.");
			} else {
				mutations[i-1] = pos;
			}
		}
		AminoAcid newAA = AminoAcid.ALA;
		
		
		// calculate surface accessibilites
		try {
			NaccessRunner naccRunner = new NaccessRunner(new File(NACCESS_EXECUTABLE), NACCESS_PARAMETERS);
			naccRunner.runNaccess(pdb);
		} catch (IOException e) {
			System.err.println("Error running NACCESS: " + e.getMessage());
			System.exit(1);
		}
		
		// calculate secondary structure
		try {
			pdb.setSecondaryStructure(DsspRunner.runDssp(pdb,DSSP_EXECUTABLE, DSSP_PARAMETERS));
		} catch (IOException e) {
			System.err.println("Error running DSSP: " + e.getMessage());
			System.exit(1);
		}
		
		// read catalytic site annotation
		try {
			CSAConnection.parseCSA(pdb,CatalSiteSet.LATEST_VERSION, false);
		} catch (IOException e) {
			System.err.println("Error acessing CSA: " + e.getMessage());
			System.exit(1);
		}

		// get set of catalytic residues
		System.out.println();
		CatalSiteSet css = pdb.getCSA();
		HashSet<Integer> catalyticResidues = new HashSet<Integer>();
		if(!css.isEmpty()) {
			System.out.println("Catalytic site annotation:");		
			System.out.println(css);
			Iterator<CatalyticSite> it = css.getIterator();
			while(it.hasNext()) {
				CatalyticSite cs = it.next();
				catalyticResidues.addAll(cs.getRes());
			}
		}
		System.out.println("Catalytic residues: " + catalyticResidues);

		// get set of surface residues
		HashSet<Integer> surfaceResidues = new HashSet<Integer>();
		for(int r:pdb.getAllSortedResSerials()) {
			if(pdb.getAllRsaFromResSerial(r) > EXPOSURE_CUTOFF) {
				surfaceResidues.add(r);
			}
		}
		
		// graph analysis
		RIGraph rig = pdb.getRIGraph("Ca", 8.0);
		// find center, find surface nodes (find catalytic nodes, find interface nodes, if pdb code known or close homolog)
		// calculate shortest paths from mutated node to all other nodes
		// report distance to center, distance to surface (distance to active site, distance to interface, if known)
		// normalize somehow
		// do statistics over cancer mutations, compare to random mutations
		UnweightedShortestPath<RIGNode, RIGEdge> sp = new UnweightedShortestPath<RIGNode, RIGEdge>(rig);
		int center = getCentralNode(rig, sp);
		System.out.println("Most central node: " + center + "\n");
		RIGNode c = rig.getNodeFromSerial(center);
		
		// print information about mutated residues (residue type, exposure, secondary structure)
		for (int i = 0; i < mutations.length; i++) {
			int pos = mutations[i];
			boolean exposed = pdb.getAllRsaFromResSerial(pos) > EXPOSURE_CUTOFF;
			char ssState = pdb.getSecondaryStructure().getSecStrucElement(pos).getType();
			System.out.printf("Position %d :                    %s %s\n", pos, pdb.getResTypeFromResSerial(pos), getMutChemPropStr(AminoAcid.getByThreeLetterCode(pdb.getResTypeFromResSerial(pos)), newAA));
			System.out.printf("Relative surface accessibility :  %2.0f%% (%s)\n", pdb.getAllRsaFromResSerial(pos), exposed?"exposed":"buried");
			System.out.printf("Secondary structure state:        %3s (%s)\n", ssState, SecStrucElement.getTypeDescription(ssState));
			System.out.printf("Distance from center:             %3d\n", sp.getDistance(c, rig.getNodeFromSerial(pos)));
			System.out.printf("Distance from surface:            %3d\n", minDistFromSet(rig, pos, sp, surfaceResidues));	
			System.out.println();
		}
		
		createImage(pdb, mutations, catalyticResidues, center, outFile);
		
	}

	/**
	 * Returns a description of the chemical properties which have changed upon mutation of an amino acid
	 */
	public static String getMutChemPropStr(AminoAcid from, AminoAcid to) {
		String str = "( ";
		if(from.isSmall() && !to.isSmall()) {
			str += "Small->Large ";
		}
		if(!from.isSmall() && to.isSmall()) {
			str += "Large->Small ";
		}
		if(from.isHydrophobic() && !to.isHydrophobic()) {
			str += "Hydrophobic->Hydrophilic ";
		}
		if(!from.isHydrophobic() && to.isHydrophobic()) {
			str += "Hydrophilic->Hydrophobic ";
		}
		if(from.isPolar() && !to.isPolar()) {
			str += "Polar->Nonpolar ";
		}
		if(!from.isPolar() && to.isPolar()) {
			str += "Nonpolar->Polar ";
		}
		str += ")";
		return str;
	}
	
	/**
	 * Returns the minimum distance (=shortest path length) between a residue and a set of residues
	 * @param rig
	 * @param distances
	 * @param set
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public static int minDistFromSet(RIGraph rig, int source, Distance distances, Set<Integer> set) {
		int min = Integer.MAX_VALUE;
		for(int setRes:set) {
			int dist = (Integer) distances.getDistance(rig.getNodeFromSerial(source), rig.getNodeFromSerial(setRes));
			if(dist < min) {
				min = dist;
			}
		}
		return min;
	}
	
	/**
	 * Returns the number a node that is most central according to closeness centrality
	 * @param rig
	 * @param distances
	 * @return the number of a most central node or -1 if no node was found
	 */
	@SuppressWarnings("unchecked")
	public static int getCentralNode(RIGraph rig, Distance distances) {
		int center = -1;
		double min = Double.MAX_VALUE;
		for(RIGNode n:rig.getVertices()) {
			int sum = 0;
			int cnt = 0;
			for(RIGNode n2:rig.getVertices()) {
				Integer d = (Integer) distances.getDistance(n, n2);
				sum += d;
				cnt++;
			}
			double avg = 1.0*sum/cnt;
			if(avg < min) {
				center = n.getResidueSerial();
				min = avg;
			}
		}
		return center;
	}
	
	public static void createImage(Pdb pdb, int[] positions, Collection<Integer> catalyticResidues, int center, File outFile) {
				
		try {
			File scriptFile = File.createTempFile("temp", ".pml");
			if(!keepPymolScript) scriptFile.deleteOnExit();
			File pdbFile = File.createTempFile("temp", ".pdb");
			pdbFile.deleteOnExit();
			pdb.writeToPDBFile(pdbFile.toString());
			PrintWriter out = new PrintWriter(scriptFile);
			out.println("Pymol script:");
			out.println("load " + pdbFile);
			out.println("hide all");
			out.println("show cartoon");
			out.println("orient");
			out.println("set sphere_scale, 4");
			out.println("set sphere_transparency, 0.6");
			out.println("set label_color, yellow");
			out.println("set label_size, 20");
			// mutations
			for(int pos:positions) {
				out.println("select resi " + pos);
				out.println("show sticks, sele");
				out.println("color red, sele");
				out.println("show spheres, sele and name ca");
				//out.println("label sele and name ca, ' %s %s' % (resi, resn)");
			}
			// catalytic residues
			for(int pos:catalyticResidues) {
				out.println("select resi " + pos);
				out.println("color yellow, sele");
				out.println("show spheres, sele and name ca");
			}
			// show center
			out.println("select resi " + center);
			out.println("color blue, sele");
			out.println("show spheres, sele and name ca");			
			
			// draw
			out.println("ray");
			out.println("png " + outFile);
			out.println("quit");
			out.close();
			
			String cmdLine = PYMOL_EXECUTABLE + " " + scriptFile;
			Process p = Runtime.getRuntime().exec(cmdLine);
			p.waitFor();
			
		} catch (IOException e) {
			System.err.println("Error while creating image file:" + e.getMessage());
			System.exit(1);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
	}
	
	/**
	 * Loads a pdb structure where arg can be a pdbcode+chaincode or a pdb file name.
	 * If something goes wrong, prints an error message and exits.
	 * @param arg a pdbcode+chaincode (e.g. 1tdrB) or a pdb file name
	 * @return the structure object
	 */
	public static Pdb readStructureOrExit(String arg) {
		Pdb pdb = null;
		
		// check if argument is a filename
		File inFile = new File(arg);
		if(inFile.canRead()) {
			System.out.println("Reading file " + arg);
			pdb = new PdbfilePdb(arg);
			try {
				String[] chains = pdb.getChains();
				System.out.println("Loading chain " + chains[0]);
				pdb.load(chains[0]);
			} catch (PdbLoadError e) {
				System.err.println("Error loading file " + arg + ":" + e.getMessage());
			}
		} else {
			// check if argument is a pdb code
			if(arg.length() < 4 || arg.length() > 5) {
				System.err.println(arg + "is neither a valid file name nor a valid pdb code");
				System.exit(1);
			} else {
				String pdbCode = arg.substring(0,4);
				String chainCode = arg.substring(4,5);
				try {
					System.out.println("Loading pdb code " + pdbCode);
					pdb = new PdbasePdb(pdbCode);
					if(chainCode.length() == 0) {
						try {
							chainCode = pdb.getChains()[0];
						} catch (PdbLoadError e) {
							System.err.println("Error loading pdb structure:" + e.getMessage());
							System.exit(1);
						}
					}
					try {
						System.out.println("Loading chain " + chainCode);
						pdb.load(pdb.getChains()[0]);
					} catch (PdbLoadError e) {
						System.err.println("Error loading pdb structure:" + e.getMessage());
						System.exit(1);
					}
				} catch (SQLException e) {
					System.err.println("Database error: " + e.getMessage());
					System.exit(1);
				} catch (PdbCodeNotFoundError e) {
					System.err.println("Pdb code " + pdbCode + " not found in database.");
					System.exit(1);
				}

			}
		}
		return pdb;
	}
	
}
