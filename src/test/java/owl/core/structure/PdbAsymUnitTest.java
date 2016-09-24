package owl.core.structure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.TreeMap;

import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

import org.junit.Assert;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.xml.sax.SAXException;

import owl.core.structure.AsaCalculator;
import owl.core.structure.ChainCluster;
import owl.core.structure.ChainInterface;
import owl.core.structure.ChainInterfaceList;
import owl.core.structure.CrystalCell;
import owl.core.structure.InterfacesFinder;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbLoadException;
import owl.core.structure.SpaceGroup;
import owl.core.util.FileFormatException;
import owl.tests.TestsSetup;


public class PdbAsymUnitTest {

	@SuppressWarnings("unused") // we keep it here since we might still want to use naccess for area calculations
	private static String NACCESS_EXEC; 
	
	private static String LOCAL_CIF_DIR;
	
	private static final String TESTDATADIR = "/owl/core/structure";
	private static final String LISTFILEPATH = TESTDATADIR+"/testset_interfaces.txt";
	private static final String LISTFILE2PATH = TESTDATADIR+"/testset_interfaces2.txt";
	private static final String CULLPDB20FILEPATH = TESTDATADIR+"/cullpdb_20";
	
	private static File LISTFILE;
	private static File LISTFILE2;
	private static File CULLPDB20FILE;
	
	private static final double CUTOFF = 5.9;
	

	private static final boolean CONSIDER_HETATOMS = true;
	private static final boolean CONSIDER_NONPOLY = false;
	private static final int     CONSIDER_COFACTORS = -1;
	
	private static final int NTHREADS = Runtime.getRuntime().availableProcessors(); // number of threads for ASA calculation
	
	// for tests that don't need to compare to pisa we shave off the small interfaces
	private static final double MIN_AREA_TO_KEEP = 35;

	private static List<File> listFileCifFiles;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		Properties p = TestsSetup.readPaths();
		NACCESS_EXEC = p.getProperty("NACCESS_EXEC");
		LOCAL_CIF_DIR = p.getProperty("LOCAL_CIF_DIR");
		
		LISTFILE = TestsSetup.inputStreamToTempFile(PdbAsymUnitTest.class.getResourceAsStream(LISTFILEPATH), "PdbAsymUnitTest", ".list");
		LISTFILE2 = TestsSetup.inputStreamToTempFile(PdbAsymUnitTest.class.getResourceAsStream(LISTFILE2PATH), "PdbAsymUnitTest", ".list");
		CULLPDB20FILE = TestsSetup.inputStreamToTempFile(PdbAsymUnitTest.class.getResourceAsStream(CULLPDB20FILEPATH), "PdbAsymUnitTest", ".list");
		
		// for pdbs in LISTFILE (used for several tests) we 
		listFileCifFiles = new ArrayList<File>();
		List<String> pdbCodes = readListFile(LISTFILE);
		System.out.println("Grabbing cif files in list "+LISTFILEPATH);
		for (String pdbCode: pdbCodes) {
			System.out.print(".");
			File cifFile = new File(System.getProperty("java.io.tmpdir"),pdbCode+".pdbasymunittest.cif");
			cifFile.deleteOnExit();
			listFileCifFiles.add(cifFile);
			PdbAsymUnit.grabCifFile(LOCAL_CIF_DIR, null, pdbCode, cifFile, false);
		}
		System.out.println();
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
	public void testChainClusters() throws IOException {
		
		System.out.println("Checking unique and representative chains");
		for (File cifFile:listFileCifFiles) {
			
			System.out.println(cifFile.getName());

			PdbAsymUnit pdb = null;
			try {
				pdb = new PdbAsymUnit(cifFile);
			} catch (PdbLoadException e) {
				System.err.println("PDB load error, cause: "+e.getMessage());
				continue;
			} catch (FileFormatException e) {
				System.err.println("PDB load error, cause: "+e.getMessage());
				continue;
			}
			
			List<ChainCluster> clusters = pdb.getProtChainClusters();
			Assert.assertTrue(clusters.size()<=pdb.getNumPolyChains() && clusters.size()>0);
			List<String> allchains = new ArrayList<String>();
			for (ChainCluster cluster:clusters) {
				for (PdbChain member:cluster.getMembers()) {
					allchains.add(member.getPdbChainCode());
				}
			}
			Assert.assertTrue(allchains.size()==pdb.getNumPolyChains());
			for (String pdbChainCode:pdb.getPdbChainCodes()) {
				Assert.assertNotNull(pdb.getProtChainCluster(pdbChainCode));
			}
			
		}
	}
	
	@Test
	public void testInterfaceClusters() throws IOException {

		List<String> pdbCodes = readListFile(LISTFILE2);

		System.out.println("Interface clusters calculation test");
		System.out.println("Will use "+NTHREADS+" CPUs for ASA calculations");

		
		for (String pdbCode: pdbCodes) {

			System.out.println("\n##"+pdbCode);
			File cifFile = new File(System.getProperty("java.io.tmpdir"),pdbCode+".pdbasymunittest.cif");
			cifFile.deleteOnExit();
			PdbAsymUnit.grabCifFile(LOCAL_CIF_DIR, null, pdbCode, cifFile, false);

			PdbAsymUnit pdb = null;
			try {
				pdb = new PdbAsymUnit(cifFile);
			} catch (PdbLoadException e) {
				System.err.println("PDB load error, cause: "+e.getMessage());
				continue;
			} catch (FileFormatException e) {
				System.err.println("PDB load error, cause: "+e.getMessage());
				continue;
			}

			pdb.removeHatoms();

			System.out.println(pdb.getSpaceGroup().getShortSymbol()+" ("+pdb.getSpaceGroup().getId()+")");

			long start = System.currentTimeMillis();

			ChainInterfaceList interfaces = 
					pdb.getAllInterfaces(CUTOFF, AsaCalculator.DEFAULT_N_SPHERE_POINTS, NTHREADS, CONSIDER_HETATOMS, CONSIDER_NONPOLY, CONSIDER_COFACTORS, 35);

			long end = System.currentTimeMillis();
			
			System.out.println("Time: "+((end-start)/1000)+"s");
			System.out.println("Total number of interfaces found (above 35A2 area): "+interfaces.getNumInterfaces());


			// testing the clusters
			
			
			System.out.print("Calculating clusters... ");
			
			start = System.currentTimeMillis();
			interfaces.initialiseClusters(pdb, 2.0, 10, "CA");
			end = System.currentTimeMillis();
			
			System.out.printf("Time: %6.3f\n",((end-start)/1000.0));
			
			List<InterfaceCluster> clusters = interfaces.getClusters();

			System.out.println("Total number of interface clusters found (above 35A2 area): "+clusters.size());
			
			// number of clusters must be smaller than number of interfaces
			Assert.assertTrue(interfaces.getNumInterfaces()>=clusters.size());
			
			HashSet<Integer> interfaceIds = new HashSet<Integer>();
			HashSet<Integer> clusterIds = new HashSet<Integer>();
			for (InterfaceCluster cluster:clusters) {
				
				// just a sanity check: clusters ids are not duplicated
				Assert.assertTrue(clusterIds.add(cluster.getId()));
				
				for (ChainInterface interf:cluster.getMembers()) {					
					// main test: no duplicates should exist 
					Assert.assertTrue("There should not be any duplicate interface ids in the interface clusters",interfaceIds.add(interf.getId()));
				}
			}
			
		}
	}

	@Test
	public void testInterfacesFinderNumCells() throws IOException, SAXException {
		
		
		List<String> pdbCodes = readListFile(CULLPDB20FILE);
		
		System.out.println("Interface calculation - number of neighboring cells test ("+pdbCodes.size()+" structures to test)");
		System.out.println("Will use "+NTHREADS+" CPUs for ASA calculations");
		
		Map<String, List<Integer>> pdbsNeedingExtraNeighbors = new TreeMap<String, List<Integer>>();
		
		for (String pdbCode: pdbCodes) {
					
			System.out.println("\n##"+pdbCode);
			File cifFile = new File(System.getProperty("java.io.tmpdir"),pdbCode+".pdbasymunittest.cif");
			cifFile.deleteOnExit();
			try {
				PdbAsymUnit.grabCifFile(LOCAL_CIF_DIR, null, pdbCode, cifFile, false);
			} catch (IOException e) {
				System.err.println("Something went wrong while grabbing cif file: "+e.getMessage());
				continue;
			}

			
			PdbAsymUnit pdb = null;
			try {
				pdb = new PdbAsymUnit(cifFile);				
			} catch (PdbLoadException e) {
				System.err.println("PDB load error, cause: "+e.getMessage());
				continue;
			} catch (FileFormatException e) {
				System.err.println("PDB load error, cause: "+e.getMessage());
				continue;
			}

			pdb.removeHatoms();

			//List<Integer> numInterfacesFound = new ArrayList<Integer>();			
			
			int numInterfaces = -1;
			for (int numCells=1;numCells<=7;numCells++) {
				long start = 0, end =0;

				InterfacesFinder interfFinder = new InterfacesFinder(pdb.copy());
				interfFinder.setNumCells(numCells);
				
				start = System.currentTimeMillis();
				
				// we use a small value for number of sphere points because we are not interested in the ASAs at all here 
				ChainInterfaceList interfaces = 
						interfFinder.getAllInterfaces(CUTOFF, 10, NTHREADS, CONSIDER_HETATOMS, CONSIDER_NONPOLY, CONSIDER_COFACTORS, MIN_AREA_TO_KEEP);
				
				end = System.currentTimeMillis();
				long total = (end-start)/100; // i.e. in deciseconds


				System.out.println("Tried "+numCells+" cells. Total number of interfaces found: "+interfaces.size());
				// time in deciseconds
				System.out.println("Time: "+total+"ds");

				if (numInterfaces!=-1 && numInterfaces!=interfaces.size()) {
					if (pdbsNeedingExtraNeighbors.containsKey(pdbCode)) {
						pdbsNeedingExtraNeighbors.get(pdbCode).add(numCells);
					} else {
						List<Integer> list = new ArrayList<Integer>();
						list.add(numCells);
						pdbsNeedingExtraNeighbors.put(pdbCode, list);
					}
				}
				
				numInterfaces = interfaces.size();
				//numInterfacesFound.add(interfaces.size());

			}
			
			// we don't actually assert anything here, because the num of interfaces will change in some cases when changing the numCells
			//for (int i=0;i<numInterfacesFound.size();i++){
			//	if (i>0)
			//		Assert.assertEquals(numInterfacesFound.get(i-1),numInterfacesFound.get(i));
			//}
			
			
		}
		
		System.out.println("PDBs that needed more than 1 cell neighbours to find all interfaces:");
		for (String pdbCode:pdbsNeedingExtraNeighbors.keySet()) {
			System.out.print(pdbCode+":");
			for (int numCells:pdbsNeedingExtraNeighbors.get(pdbCode)) {
				System.out.print(" "+numCells);
			}
			System.out.println();
		}

	}
	
	@Test
	public void testCrystalCellTransformations() throws IOException {
		
		for (File cifFile:listFileCifFiles) {
			
			System.out.println("\n##"+cifFile.getName());

			PdbAsymUnit pdb = null;
			try {
				pdb = new PdbAsymUnit(cifFile);
			} catch (PdbLoadException e) {
				System.err.println("PDB load error, cause: "+e.getMessage());
				continue;
			} catch (FileFormatException e) {
				System.err.println("PDB load error, cause: "+e.getMessage());
				continue;
			}
			
			pdb.removeHatoms();
			
			SpaceGroup sg = pdb.getSpaceGroup();
			CrystalCell cell = pdb.getCrystalCell();
			
			for (int op=0;op<sg.getNumOperators();op++) {
				Matrix4d transfXtal = new Matrix4d(sg.getTransformation(op));				
				
				Matrix4d transfOrthon = cell.transfToOrthonormal(transfXtal);
				
				Matrix4d transfBackToXtal = cell.transfToCrystal(transfOrthon);
				
				Matrix4d transfBackToOrthon = cell.transfToOrthonormal(transfBackToXtal);
				
				Assert.assertTrue(transfXtal.epsilonEquals(transfBackToXtal, 0.00001));
				
				Assert.assertTrue(transfOrthon.epsilonEquals(transfBackToOrthon, 0.00001));
				
				
				// transforms on translations done with Matrix4d or Tuple3d transformations must coincide
				Vector3d xtalTransl = new Vector3d(2,5,1);
				transfXtal.setTranslation(xtalTransl); 
				
				transfOrthon = cell.transfToOrthonormal(transfXtal);
				
				Vector3d orthonTransl = new Vector3d(xtalTransl);
				cell.transfToOrthonormal(orthonTransl); 
				
				Vector3d translFromMatrix4dTransform = new Vector3d(transfOrthon.m03, transfOrthon.m13, transfOrthon.m23);
				Assert.assertTrue(translFromMatrix4dTransform.epsilonEquals(orthonTransl, 0.0001)); 
				
			}
			
			int numCells = 3;
			for (int i=-numCells;i<=numCells;i++) {
				for (int j=-numCells;j<=numCells;j++) {
					for (int k=-numCells;k<=numCells;k++) {
						if (i==0 && j==0 && k==0) continue; 
						
						Vector3d translXtal = new Vector3d (i,j,k);
						Vector3d translOrthon = new Vector3d(translXtal);
						cell.transfToOrthonormal(translOrthon);
						Vector3d translBackToXtal = new Vector3d(translOrthon);
						cell.transfToCrystal(translBackToXtal);

						Assert.assertTrue(translXtal.epsilonEquals(translBackToXtal, 0.00001));
					}
				}
			}

		}
		
	}
	
	private static List<String> readListFile(File file) throws IOException {
		List<String> pdbCodes = new ArrayList<String>();
		BufferedReader flist = new BufferedReader(new FileReader(file));
		String line;
		while ((line = flist.readLine() ) != null ) {
			if (line.startsWith("#")) continue;
			if (line.isEmpty()) break;
			String pdbCode = line.split("\\s+")[0].toLowerCase();
			pdbCodes.add(pdbCode);
		}
		flist.close();
		return pdbCodes;
	}
	

}
