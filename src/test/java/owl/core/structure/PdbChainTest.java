package owl.core.structure;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Properties;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.Point3d;
import javax.vecmath.Point3i;
import javax.vecmath.Tuple3d;
import javax.vecmath.Vector3d;






//import junit.framework.Assert;
import org.junit.Assert;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import owl.core.runners.DsspRunner;
import owl.core.runners.MaxClusterRunner;
import owl.core.runners.NaccessRunner;
import owl.core.sequence.alignment.AlignmentConstructionException;
import owl.core.sequence.alignment.MultipleSequenceAlignment;
import owl.core.structure.Atom;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbLoadException;
import owl.core.structure.Residue; 
import owl.core.structure.AaResidue;
import owl.core.structure.TemplateList;
import owl.core.structure.features.SecondaryStructure;
import owl.core.structure.graphs.FileRIGraph;
import owl.core.structure.graphs.RIGEdge;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.FileFormatException;
import owl.core.util.GeometryTools;
import owl.core.util.Interval;
import owl.core.util.IntervalSet;
import owl.core.util.OptSuperposition;
import owl.tests.TestsSetup;
import edu.uci.ics.jung.graph.util.Pair;


public class PdbChainTest {

	private static String MAX_CLUSTER_EXEC; 
	private static String NACCESS_EXEC; 
	private static String DSSP_EXEC;
	
	private static final String CIFDIR="/nfs/data/dbs/pdb/data/structures/all/mmCIF/";
	
	private static final String TESTDATADIR = "/owl/core/structure";
	private static final String TEST_PDB_FILE_1_PATH = TESTDATADIR+"/1tdrA.pdb";
	private static final String TEST_PDB_FILE_2_PATH = TESTDATADIR+"/1tdrB.pdb";
	private static final String TEST_CHAIN_1 = "A";
	private static final String TEST_CHAIN_2 = "B";
	private static final InputStream NACCESS_OUTPUT_REF = PdbChainTest.class.getResourceAsStream(TESTDATADIR+"/1tdrA.rsa");
	private static final String TESTSET10_LIST_PATH = TESTDATADIR+"/testset10.list";
	private static final String TESTSETRMSD_LIST_PATH = TESTDATADIR+"/testset_rmsd.list";
	
	private static final String TEST_PDB_3 = "12as";
	private static final String TEST_CHAIN_3 = "A";
	private static final String TEST_PDB_4 = "12as";
	private static final String TEST_CHAIN_4 = "B";

	private static final int NTHREADS = Runtime.getRuntime().availableProcessors(); // number of threads for ASA calculation
	
	
	private static final String[] CMFILES = { 
		"1d2sA_ALL_5.0.cm", "1d2sA_BB_5.0.cm", "1d2sA_BB-SC_5.0.cm", "1d2sA_Ca_8.0.cm", "1d2sA_Cb_8.0.cm", "1d2sA_SC_5.0.cm", 
		"1od3A_ALL_5.0.cm", "1od3A_BB_5.0.cm", "1od3A_BB-SC_5.0.cm", "1od3A_Ca_8.0.cm", "1od3A_Cb_8.0.cm", "1od3A_SC_5.0.cm", 
		"1oewA_ALL_5.0.cm", "1oewA_BB_5.0.cm", "1oewA_BB-SC_5.0.cm", "1oewA_Ca_8.0.cm", "1oewA_Cb_8.0.cm", "1oewA_SC_5.0.cm", 
		"1w0nA_ALL_5.0.cm", "1w0nA_BB_5.0.cm", "1w0nA_BB-SC_5.0.cm", "1w0nA_Ca_8.0.cm", "1w0nA_Cb_8.0.cm", "1w0nA_SC_5.0.cm", 
		"1wb4A_ALL_5.0.cm", "1wb4A_BB_5.0.cm", "1wb4A_BB-SC_5.0.cm", "1wb4A_Ca_8.0.cm", "1wb4A_Cb_8.0.cm", "1wb4A_SC_5.0.cm", 
		"1wvfA_ALL_5.0.cm", "1wvfA_BB_5.0.cm", "1wvfA_BB-SC_5.0.cm", "1wvfA_Ca_8.0.cm", "1wvfA_Cb_8.0.cm", "1wvfA_SC_5.0.cm", 
		"2gpiA_ALL_5.0.cm", "2gpiA_BB_5.0.cm", "2gpiA_BB-SC_5.0.cm", "2gpiA_Ca_8.0.cm", "2gpiA_Cb_8.0.cm", "2gpiA_SC_5.0.cm", 
		"2h6fB_ALL_5.0.cm", "2h6fB_BB_5.0.cm", "2h6fB_BB-SC_5.0.cm", "2h6fB_Ca_8.0.cm", "2h6fB_Cb_8.0.cm", "2h6fB_SC_5.0.cm", 
		"3nulA_ALL_5.0.cm", "3nulA_BB_5.0.cm", "3nulA_BB-SC_5.0.cm", "3nulA_Ca_8.0.cm", "3nulA_Cb_8.0.cm", "3nulA_SC_5.0.cm", 
		"7odcA_ALL_5.0.cm", "7odcA_BB_5.0.cm", "7odcA_BB-SC_5.0.cm", "7odcA_Ca_8.0.cm", "7odcA_Cb_8.0.cm", "7odcA_SC_5.0.cm"};
	
	
	private static File TEST_PDB_FILE_1;
	private static File TEST_PDB_FILE_2;
	private static File TESTSET10_LIST;
	private static File TESTSETRMSD_LIST;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		Properties p = TestsSetup.readPaths();
		NACCESS_EXEC = p.getProperty("NACCESS_EXEC");
		MAX_CLUSTER_EXEC = p.getProperty("MAX_CLUSTER_EXEC");
		DSSP_EXEC = p.getProperty("DSSP_EXEC");
		
		TEST_PDB_FILE_1 = TestsSetup.inputStreamToTempFile(PdbAsymUnitTest.class.getResourceAsStream(TEST_PDB_FILE_1_PATH), "PdbChainTest", ".pdb");
		TEST_PDB_FILE_2 = TestsSetup.inputStreamToTempFile(PdbAsymUnitTest.class.getResourceAsStream(TEST_PDB_FILE_2_PATH), "PdbChainTest", ".pdb");
		TESTSET10_LIST = TestsSetup.inputStreamToTempFile(PdbAsymUnitTest.class.getResourceAsStream(TESTSET10_LIST_PATH), "PdbChainTest", ".list");
		TESTSETRMSD_LIST = TestsSetup.inputStreamToTempFile(PdbAsymUnitTest.class.getResourceAsStream(TESTSETRMSD_LIST_PATH), "PdbChainTest", ".list");
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testRunNaccess() throws IOException, PdbLoadException, FileFormatException {
		
		PdbAsymUnit fullpdb = new PdbAsymUnit(TEST_PDB_FILE_1);
		PdbChain pdb = fullpdb.getChain(TEST_CHAIN_1);
		Assert.assertFalse(pdb.hasASA());
		NaccessRunner naccRunner = new NaccessRunner(new File(NACCESS_EXEC), "");
		naccRunner.runNaccess(pdb);
		Assert.assertTrue(pdb.hasASA());
		
		// open our naccess reference result (rsa) file		
		BufferedReader br = new BufferedReader(new InputStreamReader(NACCESS_OUTPUT_REF));
		String line;
		int resser = 0;
		while ((line=br.readLine())!=null) {
			if (!line.startsWith("RES")) continue;
			resser++;
			double allAtomAsa = Double.parseDouble(line.split("\\s+")[5]);
			double scAsa = Double.parseDouble(line.split("\\s+")[7]);
			Assert.assertEquals(allAtomAsa, pdb.getAllRsaFromResSerial(resser), 0);
			Assert.assertEquals(scAsa, pdb.getScRsaFromResSerial(resser), 0);
		}
		br.close();
		Assert.assertEquals(resser,pdb.getStdAaObsLength());
		
		
	}
	
	@Test
	public void testASAcalcVsNaccess() throws IOException, PdbLoadException, FileFormatException {
		
		System.out.println("Using "+NTHREADS+" CPUs for ASA calculations");
		System.out.println("Matching of our ASA values against NACCESS's. Areas with disagreement >20% reported");
		String[] pdbIds = TemplateList.readIdsListFile(TESTSET10_LIST);
		
		for (String pdbId:pdbIds) {
			System.out.println(pdbId);
			String pdbCode = pdbId.substring(0,4);
			String pdbChainCode = pdbId.substring(4,5);

			PdbChain theirs = null;
			PdbChain ours = null;
			

			PdbAsymUnit theirsFull = new PdbAsymUnit(TestsSetup.gunzipFile(new File(CIFDIR,pdbCode+".cif.gz")));
			theirs = theirsFull.getChain(pdbChainCode);
			// note we run nacccess with -h to also include the het residues
			NaccessRunner naccRunner = new NaccessRunner(new File(NACCESS_EXEC),"-h");
			naccRunner.runNaccess(theirs);

			PdbAsymUnit oursFull = new PdbAsymUnit(TestsSetup.gunzipFile(new File(CIFDIR,pdbCode+".cif.gz")));
			ours = oursFull.getChain(pdbChainCode);
			long start = System.currentTimeMillis();
			ours.calcASAs(960, NTHREADS, true);
			long end = System.currentTimeMillis();
			System.out.printf("Time: %4.1fs\n",((end-start)/1000.0));


			
			checkASAsMatch(theirs, ours, 0.20, 0.11);
			
		}

		System.out.println("Matching of our ASAs from a chain to itself rotated. Areas with disagreement >15% reported");
		for (String pdbId:pdbIds) {
			// test whether values calculated for a rotated molecule match to those of the original molecule. 
			// The values don't match exactly at all! (NACCESS had this problem too). It must be something inherent to the algorithm,
			// I don't think it can be just rounding, I don't really understand it!
			System.out.println(pdbId);
			String pdbCode = pdbId.substring(0,4);
			String pdbChainCode = pdbId.substring(4,5);
			PdbChain ours = null;
			try {
				PdbAsymUnit oursFull = new PdbAsymUnit(TestsSetup.gunzipFile(new File(CIFDIR,pdbCode+".cif.gz")));
				ours = oursFull.getChain(pdbChainCode);

				
				PdbChain rotated = ours.copy(oursFull);
				PdbChain translated = ours.copy(oursFull);

				rotated.rotate(new Vector3d(0,0,1), Math.PI/4.0);
				translated.doCrystalTranslation(new Point3i(1,1,1));
				
				System.out.println("960 sphere points");
				ours.calcASAs(960,NTHREADS,true);
				rotated.calcASAs(960,NTHREADS,true);
				translated.calcASAs(960,NTHREADS,true);
				checkASAsMatch(ours,rotated, 0.15, 0.10);
				checkASAsMatch(ours,translated, 0.01, 0.01);

				System.out.println("2000 sphere points");
				ours.calcASAs(2000,NTHREADS,true);
				rotated.calcASAs(2000,NTHREADS,true);
				checkASAsMatch(ours,rotated, 0.15, 0.10);
				
				System.out.println("3000 sphere points");
				ours.calcASAs(3000,NTHREADS,true);
				rotated.calcASAs(3000,NTHREADS,true);
				checkASAsMatch(ours,rotated, 0.15, 0.10);
				
				System.out.println("9600 sphere points");
				ours.calcASAs(9600,NTHREADS,true);
				rotated.calcASAs(9600,NTHREADS,true);
				checkASAsMatch(ours,rotated, 0.15, 0.10);
				

			} catch (IOException e) {
				System.err.println("Could not find pdb code "+pdbCode+". Error: "+e.getMessage());
				continue;
			} catch (FileFormatException e) {
				System.err.println("Could not find pdb code "+pdbCode+". Error: "+e.getMessage());
				continue;				
			}
		}
		
		
	}
	
	private void checkASAsMatch(PdbChain theirs, PdbChain ours, double toleranceSingleVal, double toleranceGlobal) {
		int misMatches = 0;
		for (int resser:ours.getAllResSerials()) {
			Residue tRes = theirs.getResidue(resser);
			Residue oRes = ours.getResidue(resser);
			//String notWithin = "";
			// we allow for a max 20% discrepancy
			if (Math.abs(tRes.getAsa()-oRes.getAsa())>tRes.getAsa()*toleranceSingleVal){
				//notWithin = "x";
				misMatches++;
			}
			//System.out.printf("%6.2f\t%6.2f\t"+notWithin+"\n",tRes.getAsa(),oRes.getAsa());
			
		}
		System.out.printf("%4.1f%% mismatches (%d)\n",100.0*(double)misMatches/(double)ours.getObsLength(),misMatches);
		// we require that 90% of them must be within the 20% tolerance
		Assert.assertTrue(misMatches<(double)ours.getObsLength()*toleranceGlobal);

	}
	
	// disabling test for the moment because we only have a 32-bit version of dssp and it does not run on Ubuntu anymore...
	//@Test
	public void testRunDssp() throws PdbLoadException, IOException, FileFormatException {
		PdbAsymUnit fullpdb = new PdbAsymUnit(TEST_PDB_FILE_1);
		PdbChain pdb = fullpdb.getChain(TEST_CHAIN_1);

		PdbAsymUnit dsspAssignedFullpdb = new PdbAsymUnit(TEST_PDB_FILE_1);
		PdbChain pdbDsspAssigned = dsspAssignedFullpdb.getChain(TEST_CHAIN_1);

		SecondaryStructure secondaryStructure = DsspRunner.runDssp(pdb, DSSP_EXEC, "--");
		pdbDsspAssigned.setSecondaryStructure(secondaryStructure);
		
		// TODO the TEST_PDB_FILE_1 should contain the secondary structure so that we could perform this comparison 
		//for (int resser:pdb.getAllSortedResSerials()) {
		//	Assert.assertEquals(pdb.getResidue(resser).getSsElem(),pdbDsspAssigned.getResidue(resser).getSsElem());
		//}
		// all we can do right now is check for nulls
		Assert.assertNotNull(pdbDsspAssigned.getSecondaryStructure());
		for (int resser:pdb.getAllStdAaResSerials()) {
			Assert.assertNotNull(pdbDsspAssigned.getResidue(resser).getSsElem());
		}
	}
	
	@Test
	public void testRmsd() throws PdbLoadException, IOException, FileFormatException {
		
		PdbAsymUnit fullpdb = new PdbAsymUnit(TEST_PDB_FILE_1);
		PdbChain pdb1 = fullpdb.getChain(TEST_CHAIN_1);

		PdbAsymUnit fullpdb1p = new PdbAsymUnit(TEST_PDB_FILE_1);
		PdbChain pdb1p = fullpdb1p.getChain(TEST_CHAIN_1);

		PdbAsymUnit fullpdb2 = new PdbAsymUnit(TEST_PDB_FILE_2);
		PdbChain pdb2 = fullpdb2.getChain(TEST_CHAIN_2);

		
		String[] cts = {"Ca", "Cb", "ALL", "BB", "SC"};
		for (String ct:cts) {
			Assert.assertEquals(0.0,pdb1.rmsd(pdb1p, ct),0.0001);
			Assert.assertEquals(0.0,pdb1p.rmsd(pdb1, ct),0.0001);			
		}
		// rmsd on intervals check
		IntervalSet intervSet = new IntervalSet();
		intervSet.add(new Interval(pdb1.getFirstResidue().getSerial()+10,pdb1.getLastResidue().getSerial()-10));
		// check that rmsd is 0 on an interval
		for (String ct:cts) {
			Assert.assertEquals(0.0,pdb1.rmsd(pdb1p, ct, intervSet),0.0001);
			Assert.assertEquals(0.0,pdb1p.rmsd(pdb1, ct, intervSet),0.0001);			
		}

		MaxClusterRunner mcr = new MaxClusterRunner(MAX_CLUSTER_EXEC);
		double mcrmsd = mcr.calculatePairwiseScoreSeqDependent(TEST_PDB_FILE_1.getAbsolutePath(), TEST_PDB_FILE_2.getAbsolutePath(), MaxClusterRunner.ScoreType.RMSD);
		double ourrmsd = pdb1.rmsd(pdb2, "Ca");
		Assert.assertEquals(mcrmsd, ourrmsd, 0.001);
		ourrmsd = pdb2.rmsd(pdb1, "Ca");
		Assert.assertEquals(mcrmsd, ourrmsd, 0.001);
		
		
		// a more exhaustive test on a few test entries
		
		String[] pdbIds = TemplateList.readIdsListFile(TESTSETRMSD_LIST);
		
		for (String pdbId:pdbIds) {
			
			PdbAsymUnit pdbi = new PdbAsymUnit(TestsSetup.gunzipFile(new File(CIFDIR,pdbId.substring(0,4)+".cif.gz")));
			PdbAsymUnit pdbip = new PdbAsymUnit(TestsSetup.gunzipFile(new File(CIFDIR,pdbId.substring(0,4)+".cif.gz")));
		
			// first we write pdb files for every chain
			ArrayList<File> pdbFiles = new ArrayList<File>();
			for (PdbChain pdbChaini:pdbi.getPolyChains()) {
				
				File pdbFile = File.createTempFile("rmsdtest_"+pdbi.getPdbCode()+pdbChaini.getPdbChainCode()+"_",".pdb");
				pdbFile.deleteOnExit();
				pdbChaini.writeToPDBFileWithPdbChainCodes(pdbFile, false); 
				pdbFiles.add(pdbFile);
			}	
			
		
			int i = 0;
			for (PdbChain pdbChaini:pdbi.getPolyChains()) {
				
				PdbChain pdbpChaini = pdbip.getChain(pdbChaini.getPdbChainCode());
				
				for (String ct:cts) {
					Assert.assertEquals(0.0,pdbChaini.rmsd(pdbpChaini, ct),0.0001);
					Assert.assertEquals(0.0,pdbpChaini.rmsd(pdbChaini, ct),0.0001);			
				}
				
				// rmsd on intervals check
				intervSet = new IntervalSet();
				intervSet.add(new Interval(pdbChaini.getFirstResidue().getSerial()+10,
						pdbChaini.getLastResidue().getSerial()-10));
				// check that rmsd is 0 on an interval
				try {
					for (String ct:cts) {
						Assert.assertEquals(0.0,pdbChaini.rmsd(pdbpChaini, ct, intervSet),0.0001);
						Assert.assertEquals(0.0,pdbpChaini.rmsd(pdbChaini, ct, intervSet),0.0001);			
					}
				} catch (IllegalArgumentException e) {
					// whenever a residue from the interval is unobserved, rmsd throws an IllegalArgumentException
				}

				i++;
				int j = 0;
				for (PdbChain pdbChainj:pdbi.getPolyChains()) {
					
					if (!pdbi.areChainsInSameCluster(pdbChaini.getPdbChainCode(), pdbChainj.getPdbChainCode())) {
						j++;
						continue;						
					}
						
					// in order to make sure that the rmsd call is not altering the coordinates we copy here and compare in the end
					PdbChain newPdbChaini = pdbChaini.copy(null);
					PdbChain newPdbChainj = pdbChainj.copy(null);
					
					mcrmsd = mcr.calculatePairwiseScoreSeqDependent(pdbFiles.get(i-1).getAbsolutePath(), pdbFiles.get(j).getAbsolutePath(), MaxClusterRunner.ScoreType.RMSD);
					ourrmsd = pdbChaini.rmsd(pdbChainj, "Ca");
					System.out.println(pdbi.getPdbCode()+": "+pdbChaini.getPdbChainCode()+" vs "+pdbChainj.getPdbChainCode()+": "
							+String.format("%.3f", ourrmsd)+" (ours) "+String.format("%.3f", mcrmsd)+" (maxcluster)" );
					
					Assert.assertEquals(mcrmsd, ourrmsd, 0.001);
					ourrmsd = pdbChainj.rmsd(pdbChaini, "Ca");
					Assert.assertEquals(mcrmsd, ourrmsd, 0.001);
					
					// finally checking that the rmsd call didn't alter the coordinates
					Assert.assertEquals(0.0,pdbChaini.rmsd(newPdbChaini, "Ca"),0.000001);
					Assert.assertEquals(0.0,pdbChainj.rmsd(newPdbChainj, "Ca"),0.000001);
					
					j++;
				}
			}
			
			
			// test for rmsd call with transform=true
			i = 0;
			for (PdbChain pdbChaini:pdbi.getPolyChains()) {

				i++;
				
				// when transforming the coordinates we can't use a chain twice (because the coordinates are altered)				
				if (i>1) continue;
				// In any case beware this test fails if the i chain is first used with n atoms and then later used 
				// again with m (>n) atoms. The left-over m-n atoms had not been transformed and the rmsd will not coincide
				// since n atoms will be in the transformed positions and m-n will be in their original positions
				// e.g. this happens for 4fc4
								
				int j = 0;
				for (PdbChain pdbChainj:pdbi.getPolyChains()) {
					
					if (!pdbi.areChainsInSameCluster(pdbChaini.getPdbChainCode(), pdbChainj.getPdbChainCode())) {
						j++;
						continue;						
					}					
					
					ArrayList<Tuple3d> conf1AL = new ArrayList<Tuple3d>();
					ArrayList<Tuple3d> conf2AL = new ArrayList<Tuple3d>();	
					// there might be unobserved residues or some missing atoms for a residue
					// here we get the ones that are in common
					for (int iresser:pdbChaini.getAllResSerials()) {
						Residue iRes = pdbChaini.getResidue(iresser);
						if (!(iRes instanceof AaResidue)) continue;
						if (pdbChainj.containsResidue(iresser)) { 
							Residue jRes = pdbChainj.getResidue(iresser);
							if (!(jRes instanceof AaResidue)) continue;
							if (iRes.containsAtom("CA") && jRes.containsAtom("CA")) {
								conf1AL.add(iRes.getAtom("CA").getCoords());
								conf2AL.add(jRes.getAtom("CA").getCoords());
							}
						}
					}

					// converting the ArrayLists to arrays
					Tuple3d[] conformation1 = new Tuple3d[conf1AL.size()]; 
					Tuple3d[] conformation2 = new Tuple3d[conf2AL.size()];
					conf1AL.toArray(conformation1);
					conf2AL.toArray(conformation2);
					
					System.out.println("rmsd with transform for "+pdbChaini.getPdbChainCode()+" vs "+pdbChainj.getPdbChainCode()+" on "+conformation1.length+" common CA atoms");
					
					OptSuperposition os = GeometryTools.calcOptimalSuperposition(conformation1, conformation2, true);
					double rmsd = os.getRmsd();

					mcrmsd = mcr.calculatePairwiseScoreSeqDependent(pdbFiles.get(i-1).getAbsolutePath(), pdbFiles.get(j).getAbsolutePath(), MaxClusterRunner.ScoreType.RMSD);
					
					Assert.assertEquals(mcrmsd, rmsd, 0.001);
					
					// finally checking that the rmsd call did alter the coordinates
					Point3d icentroid = GeometryTools.getCentroid(conformation1);
					Point3d jcentroid = GeometryTools.getCentroid(conformation2);
					
					Assert.assertTrue(icentroid.epsilonEquals(jcentroid, 0.00001));
					
					Assert.assertEquals(rmsd, GeometryTools.getCoordinatesRmsd(conformation1, conformation2),0.00001);

					j++;
				}
			}
		}
	}
	
	@Test
	public void testGetGraph() throws IOException, PdbLoadException, FileFormatException {
		// NOTE: to test we compare to previously calculated graphs with an older version of our software
		//		In principle edges should coincide 100%. However there can be rounding problems that make one
		//		version consider a border-case to be an edge or not. Indeed I've seen there is a difference 
		//		with graphs calculated with our genGraph program compiled with java 1.5 vs to the one compiled 
		//		with java 1.6.
		
		
		String[] cts = 		{"Ca", "Cb", "ALL", "BB", "SC", "BB/SC"};
		double[] cutoffs = 	{ 8.0,  8.0,   5.0,  5.0,  5.0,     5.0};
		
		
		// first we read files to which we will compare to
		HashMap<String, RIGraph> refGraphs = new HashMap<String, RIGraph>();
		for (String fileName:CMFILES) {
			Pattern p = Pattern.compile("^(\\d\\w\\w\\w\\w)_([A-Za-z\\-]+)_(\\d\\.\\d)\\.cm$");
			Matcher m = p.matcher(fileName);
			if (m.matches()) {
				String pdbId = m.group(1);
				String ct = m.group(2);
				ct = ct.replace("-", "/");
				String cutoff = m.group(3);
				File file = TestsSetup.inputStreamToTempFile(
						PdbChainTest.class.getResourceAsStream(TESTDATADIR+"/"+fileName), fileName, "");
				refGraphs.put(pdbId+ct+cutoff,new FileRIGraph(file.getAbsolutePath()));
			}
		}
		
		String[] pdbIds = TemplateList.readIdsListFile(TESTSET10_LIST);

		for (String pdbId:pdbIds) {
			System.out.print(pdbId+"\t");
			String pdbCode = pdbId.substring(0,4);
			String pdbChainCode = pdbId.substring(4,5);

			PdbAsymUnit fullpdb = new PdbAsymUnit(TestsSetup.gunzipFile(new File(CIFDIR,pdbCode+".cif.gz")));
			PdbChain pdb = fullpdb.getChain(pdbChainCode);
			//System.out.println(pdb.get_length());
			
			for (int i=0;i<cts.length;i++) {
				System.out.print(cts[i]+"\t"+cutoffs[i]+"\t");
				RIGraph graph = pdb.getRIGraph(cts[i], cutoffs[i]);
				
				// getting the refGraph from pre-read reference files
				RIGraph refGraph = refGraphs.get(pdbId+cts[i]+cutoffs[i]);
				// asserting
				Assert.assertEquals(refGraph.getPdbCode(), graph.getPdbCode());
				Assert.assertEquals(refGraph.getPdbChainCode(), graph.getPdbChainCode());
				Assert.assertEquals(refGraph.getChainCode(), graph.getChainCode());
				Assert.assertEquals(refGraph.getContactType(), graph.getContactType());
				Assert.assertEquals(refGraph.getCutoff(), graph.getCutoff(),0.0001);
				Assert.assertEquals(refGraph.getSequence(), graph.getSequence());
				Assert.assertEquals(refGraph.getFullLength(), graph.getFullLength());
				// can't check these two. When reading from file present nodes are stricyly those with edges,
				// however when reading from db and calculating graph present nodes are all observed nodes, even 
				// if they don't have any edges
				//Assert.assertEquals(refGraph.getObsLength(), graph.getObsLength());
				//Assert.assertEquals(refGraph.getVertexCount(), graph.getVertexCount());
				Assert.assertEquals(refGraph.getEdgeCount(), graph.getEdgeCount());
				for (RIGNode refNode:refGraph.getVertices()) {
					RIGNode node = graph.getNodeFromSerial(refNode.getResidueSerial());
					Assert.assertEquals(refNode.getResidueType(),node.getResidueType());
				}
				int nodesPresentInEdgeList = 0;
				int nodesMissingInEdgeList = 0;
				for (RIGNode node:graph.getVertices()) {
					// we have to check if it exists in the refGraph read from file, because there can be missing nodes
					// when observed residues have no contacts
					if (refGraph.containsVertexI(node.getResidueSerial())) {
						RIGNode refNode = refGraph.getNodeFromSerial(node.getResidueSerial());
						Assert.assertEquals(refNode.getResidueType(),node.getResidueType());
						nodesPresentInEdgeList++;
					} else {
						nodesMissingInEdgeList++;
					}	
				}
				// we at least require that not more than 5% weren't present in the edge list
				Assert.assertTrue(((double)nodesMissingInEdgeList/(double)nodesPresentInEdgeList)<0.05);
				
				for (RIGEdge edge:refGraph.getEdges()) {
					Pair<RIGNode> pair = refGraph.getEndpoints(edge);
					int iresser = pair.getFirst().getResidueSerial();
					int jresser = pair.getSecond().getResidueSerial();
					Assert.assertTrue(graph.containsEdgeIJ(iresser, jresser));
				}
				// and we check the opposite to be sure
				for (RIGEdge edge:graph.getEdges()) {
					Pair<RIGNode> pair = graph.getEndpoints(edge);
					int iresser = pair.getFirst().getResidueSerial();
					int jresser = pair.getSecond().getResidueSerial();
					Assert.assertTrue(refGraph.containsEdgeIJ(iresser, jresser));
				}
			}
			System.out.println();
		}
		
	}
	
	@Test
	public void testConstructors() throws PdbLoadException, IOException, FileFormatException {
		PdbAsymUnit fullpdb = new PdbAsymUnit(TEST_PDB_FILE_1);
		PdbChain pdb = fullpdb.getChain(TEST_CHAIN_1);
		Vector3d[] conformation = new Vector3d[pdb.getObsLength()];
		int i=0;
		for (int resser:pdb.getAllResSerials()) {
			Atom atom = pdb.getResidue(resser).getAtom("CA");
			conformation[i]=new Vector3d(atom.getCoords());
			i++;
		}
		PdbChain model = new PdbChain(pdb.getSequence(), conformation, "CA");
		Assert.assertEquals(PdbAsymUnit.NO_PDB_CODE, model.getPdbCode());
		Assert.assertEquals(PdbChain.DEFAULT_CHAIN, model.getChainCode());
		Assert.assertEquals(PdbChain.DEFAULT_CHAIN, model.getPdbChainCode());
		Assert.assertEquals(pdb.getObsLength(), model.getObsLength());
		Assert.assertEquals(model.getObsLength(),model.getNumAtoms());
		Assert.assertEquals(pdb.getFullLength(), model.getFullLength());
		Assert.assertEquals(pdb.getSequence(), model.getSequence());
		for (int resser:model.getAllResSerials()) {
			Residue residue = model.getResidue(resser);
			Residue refRes = pdb.getResidue(resser);
			Assert.assertEquals(((AaResidue)refRes).getAaType(),((AaResidue)residue).getAaType());
			Assert.assertEquals(refRes.getAtom("CA").getCoords(),residue.getAtom("CA").getCoords());
		}
		
		Assert.assertEquals(0.0,model.rmsd(pdb, "Ca"),0.001);
		// this also tests whether rmsd method properly finds the common atoms among the 2 structures
		Assert.assertEquals(0.0,model.rmsd(pdb, "ALL"),0.001); 
		
	}
	
	
	//@Test
	public void testCalculateAtomDistMatrix() {
		//TODO get the atom distance matrix and at least compare to a contact map
	}
	
	@Test
	public void testGetDiffDistMap () throws PdbLoadException, AlignmentConstructionException, FileNotFoundException, IOException, FileFormatException {
		
		System.out.println("Loading pdb objects...");
		PdbAsymUnit pdbfull1 = new PdbAsymUnit(TestsSetup.gunzipFile(new File(CIFDIR,TEST_PDB_3+".cif.gz")));
		PdbChain pdb1 = pdbfull1.getChain(TEST_CHAIN_3);
		Assert.assertNotNull(pdb1);
		PdbAsymUnit pdbfull2 = new PdbAsymUnit(TestsSetup.gunzipFile(new File(CIFDIR,TEST_PDB_4+".cif.gz")));
		PdbChain pdb2 = pdbfull2.getChain(TEST_CHAIN_4); 

		Assert.assertNotNull(pdb2);
		
		System.out.println("Calculating distance maps...");
		HashMap<Pair<Integer>,Double> distMap1 = pdb1.calcDistMatrix("Ca");
		Assert.assertNotNull(distMap1);
		HashMap<Pair<Integer>,Double> distMap2 = pdb2.calcDistMatrix("Ca");
		Assert.assertNotNull(distMap2);
		
		// create an alignment
		String name1 = TEST_PDB_3+TEST_CHAIN_3;
		String name2 = TEST_PDB_4+TEST_CHAIN_4;
		String[] tags = {name1, name2};
		String[] seqs = {pdb1.getSequence().getSeq(), pdb2.getSequence().getSeq()};
		MultipleSequenceAlignment ali = new MultipleSequenceAlignment(tags, seqs);
		
		System.out.println("Calculating difference distance map...");
		HashMap<Pair<Integer>,Double> diffDistMap = pdb1.getDiffDistMap("Ca", pdb2, "Ca", ali, name1, name2);
		Assert.assertNotNull(diffDistMap);
		Assert.assertEquals(distMap1.size(), diffDistMap.size());
		Assert.assertEquals(distMap2.size(), diffDistMap.size());
		
//		for(int i=1; i<=pdb1.getFullLength();i++) {
//			if (!pdb1.hasCoordinates(i)) System.out.print(i+" ");	
//		}
//		System.out.println();
//
//		for(int i=1; i<=pdb2.getFullLength();i++) {
//			if (!pdb2.hasCoordinates(i)) System.out.print(i+" ");
//		}
//		System.out.println();

		System.out.print("Missing matrix values (dim=" + pdb1.getFullLength() + "): ");
		int mis = 0;
		for(int i=1; i<=pdb1.getFullLength();i++) {
			for(int j=1; j<i;j++) {
				if(!diffDistMap.containsKey(new Pair<Integer>(j,i))) {
					//System.out.print("(" + i + "," + j + ") ");
					mis++;
				} else {
					Pair<Integer> e = new Pair<Integer>(j,i);
					double dist1 = distMap1.get(e);
					double dist2 = distMap2.get(e);
					double diffDist = diffDistMap.get(e);
					Assert.assertEquals(Math.abs(dist1 - dist2), diffDist, 0.00001);					
				}
			}
		}
		System.out.println(mis);
		
		// assuming trivial alignment
		int expMissing = pdb1.getFullLength()*(pdb1.getFullLength()-1)/2 - pdb1.getStdAaObsLength()*(pdb1.getStdAaObsLength()-1)/2;
		Assert.assertEquals(expMissing, mis);
		
		double min = Collections.min(diffDistMap.values());
		double max = Collections.max(diffDistMap.values());
		
		System.out.println("Checking difference distance map for " + TEST_PDB_3+TEST_CHAIN_3 + " and " + TEST_PDB_4+TEST_CHAIN_4);
		System.out.println("size=" + diffDistMap.size() + " min=" + min + " max= " + max);

	}
	
	@Test
	public void testGetAllPhiPsi() throws IOException, PdbLoadException, FileFormatException {
		String[] pdbIds = TemplateList.readIdsListFile(TESTSET10_LIST);

		for (String pdbId:pdbIds) {
			System.out.print(pdbId+"\t");
			String pdbCode = pdbId.substring(0,4);
			String pdbChainCode = pdbId.substring(4,5);

			PdbAsymUnit pdbfull = new PdbAsymUnit(TestsSetup.gunzipFile(new File(CIFDIR,pdbCode+".cif.gz"))); 
			PdbChain pdb = pdbfull.getChain(pdbChainCode);

			TreeMap<Integer,double[]> phipsi = pdb.getAllPhiPsi();
			Assert.assertEquals(pdb.getObsLength(), phipsi.size());
			// Can't do much asserting here, ideas?
			// At least we check if there aren't exceptions in running getAllPhiPsi, which 
			// in turn is a good test for hasCoordinates(int,String)
		}
	}
	
	// to debug the testing code (run as java program so that we can use normal debugger)
	public static void main(String[] args) throws Exception {
		PdbChainTest pdbTest = new PdbChainTest();
		setUpBeforeClass();
		//pdbTest.testRunDssp();
		pdbTest.testRmsd();
		//pdbTest.testRunNaccess();
		//pdbTest.testASAcalcVsNaccess();
	}
}
