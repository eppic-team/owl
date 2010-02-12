package tests.proteinstructure;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.Vector3d;

import junit.framework.Assert;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import edu.uci.ics.jung.graph.util.Pair;

import proteinstructure.Alignment;
import proteinstructure.AlignmentConstructionError;
import proteinstructure.Atom;
import proteinstructure.ConformationsNotSameSizeError;
import proteinstructure.FileRIGraph;
import proteinstructure.FileFormatError;
import proteinstructure.MaxClusterRunner;
import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbasePdb;
import proteinstructure.PdbfilePdb;
import proteinstructure.RIGEdge;
import proteinstructure.RIGNode;
import proteinstructure.RIGraph;
import proteinstructure.Residue;
import proteinstructure.TemplateList;
import tools.Interval;
import tools.IntervalSet;
import tools.MySQLConnection;
import tools.RegexFileFilter;

public class PdbTest {

	private static final String TESTDATADIR = "tests/proteinstructure/data";
	private static final String TEST_PDB_FILE_1 = "tests/proteinstructure/data/1tdrA.pdb";
	private static final String TEST_PDB_FILE_2 = "tests/proteinstructure/data/1tdrB.pdb";
	private static final String TEST_CHAIN_1 = "A";
	private static final String TEST_CHAIN_2 = "B";
	private static final String NACCESS_EXEC = "/project/StruPPi/Software/naccess2.1.1/naccess";
	private static final String NACCESS_OUTPUT_REF = "tests/proteinstructure/data/1tdrA.rsa";
	private static final String TESTSET10_LIST = "tests/proteinstructure/data/testset10.list";
	
	private static final String TEST_PDB_3 = "12as";
	private static final String TEST_CHAIN_3 = "A";
	private static final String TEST_PDB_4 = "12as";
	private static final String TEST_CHAIN_4 = "B";

	
	private static final String MAX_CLUSTER_EXE = "/project/StruPPi/bin/maxcluster";
	private static final String PDBASE_DB = "pdbase";
	private static final String MYSQLSERVER = "talyn";

	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
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
	public void testRunNaccess() throws IOException, PdbLoadError {
		
		Pdb pdb = new PdbfilePdb(TEST_PDB_FILE_1);
		pdb.load(TEST_CHAIN_1);
		Assert.assertFalse(pdb.hasASA());
		pdb.runNaccess(NACCESS_EXEC, "");
		Assert.assertTrue(pdb.hasASA());
		
		// open our naccess reference result (rsa) file
		File outRef = new File(NACCESS_OUTPUT_REF);
		BufferedReader br = new BufferedReader(new FileReader(outRef));
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
		Assert.assertEquals(resser,pdb.getObsLength());
		
		
	}
	
	@Test
	public void testCheckConsurfHssp() throws IOException, PdbLoadError {
		Pdb pdb = new PdbfilePdb(TEST_PDB_FILE_1);
		pdb.load(TEST_CHAIN_1);
		pdb.checkConsurfHssp(false);
		for (int resser:pdb.getAllSortedResSerials()){
			// very basic test, not sure how to test this better 
			Assert.assertNotNull(pdb.getConsurfhsspColorFromResSerial(resser));
			Assert.assertNotNull(pdb.getConsurfhsspScoreFromResSerial(resser));
		}
 
	}
	
	@Test
	public void testRmsd() throws PdbLoadError, ConformationsNotSameSizeError, IOException {
		Pdb pdb1 = new PdbfilePdb(TEST_PDB_FILE_1);
		pdb1.load(TEST_CHAIN_1);
		Pdb pdb1p = new PdbfilePdb(TEST_PDB_FILE_1);
		pdb1p.load(TEST_CHAIN_1);
		Pdb pdb2 = new PdbfilePdb(TEST_PDB_FILE_2);
		pdb2.load(TEST_CHAIN_2);
		
		String[] cts = {"Ca", "Cb", "ALL", "BB", "SC", "Cg"};
		for (String ct:cts) {
			Assert.assertEquals(0.0,pdb1.rmsd(pdb1p, ct),0.0001);
			Assert.assertEquals(0.0,pdb1p.rmsd(pdb1, ct),0.0001);			
		}
		// rmsd on intervals check
		IntervalSet intervSet = new IntervalSet();
		intervSet.add(new Interval(pdb1.getMinObsResSerial()+10,pdb1.getMaxObsResSerial()-10));
		// check that rmsd is 0 on an interval
		for (String ct:cts) {
			Assert.assertEquals(0.0,pdb1.rmsd(pdb1p, ct, intervSet),0.0001);
			Assert.assertEquals(0.0,pdb1p.rmsd(pdb1, ct, intervSet),0.0001);			
		}

		MaxClusterRunner mcr = new MaxClusterRunner(MAX_CLUSTER_EXE);
		double mcrmsd = mcr.calculatePairwiseScore(TEST_PDB_FILE_1, TEST_PDB_FILE_2, MaxClusterRunner.ScoreType.RMSD);
		double ourrmsd = pdb1.rmsd(pdb2, "Ca");
		Assert.assertEquals(mcrmsd, ourrmsd, 0.001);
		ourrmsd = pdb2.rmsd(pdb1, "Ca");
		Assert.assertEquals(mcrmsd, ourrmsd, 0.001);
	}
	
	@Test
	public void testGetGraph() throws IOException, PdbCodeNotFoundError, SQLException, PdbLoadError, FileFormatError {
		// NOTE: to test we compare to previously calculated graphs with an older version of our software
		//		In principle edges should coincide 100%. However there can be rounding problems that make one
		//		version consider a border-case to be an edge or not. Indeed I've seen there is a difference 
		//		with graphs calculated with our genGraph program compiled with java 1.5 vs to the one compiled 
		//		with java 1.6.
		
		
		String[] cts = 		{"Ca", "Cb", "ALL", "BB", "SC", "BB/SC"};
		double[] cutoffs = 	{ 8.0,  8.0,   5.0,  5.0,  5.0,     5.0};
		
		// first we read files to which we will compare to
		File dir = new File(TESTDATADIR);
		File[] files = dir.listFiles(new RegexFileFilter("^\\d\\w\\w\\w\\w_.*\\.cm$"));
		HashMap<String, RIGraph> refGraphs = new HashMap<String, RIGraph>();
		for (File file:files) {
			Pattern p = Pattern.compile("^(\\d\\w\\w\\w\\w)_([A-Za-z:]+)_(\\d\\.\\d)\\.cm$");
			Matcher m = p.matcher(file.getName());
			if (m.matches()) {
				String pdbId = m.group(1);
				String ct = m.group(2);
				ct = ct.replace(":", "/");
				String cutoff = m.group(3);
				refGraphs.put(pdbId+ct+cutoff,new FileRIGraph(file.getAbsolutePath()));
			}
		}
		
		String[] pdbIds = TemplateList.readIdsListFile(new File(TESTSET10_LIST));
		MySQLConnection conn = new MySQLConnection(MYSQLSERVER,PDBASE_DB);
		for (String pdbId:pdbIds) {
			System.out.print(pdbId+"\t");
			String pdbCode = pdbId.substring(0,4);
			String pdbChainCode = pdbId.substring(4,5);

			Pdb pdb = new PdbasePdb(pdbCode,PDBASE_DB,conn);
			pdb.load(pdbChainCode);
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
	public void testConstructors() throws PdbLoadError, ConformationsNotSameSizeError {
		Pdb pdb = new PdbfilePdb(TEST_PDB_FILE_1);
		pdb.load(TEST_CHAIN_1);
		Vector3d[] conformation = new Vector3d[pdb.getObsLength()];
		int i=0;
		for (int resser:pdb.getAllSortedResSerials()) {
			Atom atom = pdb.getResidue(resser).getAtom("CA");
			conformation[i]=new Vector3d(atom.getCoords());
			i++;
		}
		Pdb model = new Pdb(pdb.getSequence(), conformation, "CA");
		Assert.assertTrue(pdb.isDataLoaded()); // doesn't make much sense here but still good to test for sanity
		Assert.assertEquals(Pdb.NO_PDB_CODE, model.getPdbCode());
		Assert.assertEquals(Pdb.DEFAULT_CHAIN, model.getChainCode());
		Assert.assertEquals(Pdb.DEFAULT_CHAIN, model.getPdbChainCode());
		Assert.assertEquals(pdb.getObsLength(), model.getObsLength());
		Assert.assertEquals(model.getObsLength(),model.getNumAtoms());
		Assert.assertEquals(pdb.getFullLength(), model.getFullLength());
		Assert.assertEquals(pdb.getSequence(), model.getSequence());
		for (int resser:model.getAllSortedResSerials()) {
			Residue residue = model.getResidue(resser);
			Residue refRes = pdb.getResidue(resser);
			Assert.assertEquals(refRes.getAaType(),residue.getAaType());
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
	public void testGetDiffDistMap () throws SQLException, PdbCodeNotFoundError, PdbLoadError, AlignmentConstructionError {
		
		MySQLConnection conn = new MySQLConnection(MYSQLSERVER,PDBASE_DB);
		
		System.out.println("Loading pdb objects...");
		Pdb pdb1 = new PdbasePdb(TEST_PDB_3,PDBASE_DB,conn);
		pdb1.load(TEST_CHAIN_3);
		Assert.assertNotNull(pdb1);
		Pdb pdb2 = new PdbasePdb(TEST_PDB_4,PDBASE_DB,conn);
		pdb2.load(TEST_CHAIN_4);
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
		String[] seqs = {pdb1.getSequence(), pdb2.getSequence()};
		Alignment ali = new Alignment(tags, seqs);
		
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
					Assert.assertEquals(Math.abs(dist1 - dist2), diffDist);					
				}
			}
		}
		System.out.println(mis);
		
		// assuming trivial alignment
		int expMissing = pdb1.getFullLength()*(pdb1.getFullLength()-1)/2 - pdb1.getObsLength()*(pdb1.getObsLength()-1)/2;
		Assert.assertEquals(expMissing, mis);
		
		double min = Collections.min(diffDistMap.values());
		double max = Collections.max(diffDistMap.values());
		
		System.out.println("Checking difference distance map for " + TEST_PDB_3+TEST_CHAIN_3 + " and " + TEST_PDB_4+TEST_CHAIN_4);
		System.out.println("size=" + diffDistMap.size() + " min=" + min + " max= " + max);

	}
	
	@Test
	public void testGetAllPhiPsi() throws IOException, PdbCodeNotFoundError, SQLException, PdbLoadError {
		String[] pdbIds = TemplateList.readIdsListFile(new File(TESTSET10_LIST));
		MySQLConnection conn = new MySQLConnection(MYSQLSERVER,PDBASE_DB);
		for (String pdbId:pdbIds) {
			System.out.print(pdbId+"\t");
			String pdbCode = pdbId.substring(0,4);
			String pdbChainCode = pdbId.substring(4,5);

			Pdb pdb = new PdbasePdb(pdbCode,PDBASE_DB,conn);
			pdb.load(pdbChainCode);
			
			TreeMap<Integer,double[]> phipsi = pdb.getAllPhiPsi();
			Assert.assertEquals(pdb.getObsLength(), phipsi.size());
			// Can't do much asserting here, ideas?
			// At least we check if there aren't exceptions in running getAllPhiPsi, which 
			// in turn is a good test for hasCoordinates(int,String)
		}
	}
	
	// to debug the testing code (run as java program so that we can use normal debugger)
	public static void main(String[] args) throws Exception {
		PdbTest pdbTest = new PdbTest();
		pdbTest.testGetAllPhiPsi();
	}
}
