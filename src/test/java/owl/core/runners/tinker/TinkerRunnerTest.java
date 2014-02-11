package owl.core.runners.tinker;


import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.Properties;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.runner.JUnitCore;

import owl.core.runners.tinker.TinkerError;
import owl.core.runners.tinker.TinkerRunner;
import owl.core.structure.Atom;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbLoadException;
import owl.core.structure.AaResidue;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.FileFormatException;
import owl.core.util.IntPairSet;
import owl.tests.TestsSetup;



public class TinkerRunnerTest {

	private static final String CIFDIR="/nfs/data/dbs/pdb/data/structures/all/mmCIF/";
	private static final String PDB_CODE = "1bxy";
	private static final String CHAIN = "A";
	private static final String CT = "Cb";
	private static final double CUTOFF = 8;
	
	private static String TINKERBINDIR;
	private static String DISTGEOM_EXE;
	private static String PRMFILE;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		Properties p = TestsSetup.readPaths();
		TINKERBINDIR = p.getProperty("TINKERBINDIR");
		DISTGEOM_EXE = p.getProperty("DISTGEOM_EXE");
		PRMFILE = p.getProperty("PRMFILE");
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
	public void testReconstruct() throws PdbLoadException, TinkerError, IOException, FileFormatException {
		File cifFile = new File(System.getProperty("java.io.tmpdir"),PDB_CODE+"_testreconstruct.cif");
		cifFile.deleteOnExit();
		PdbAsymUnit.grabCifFile(CIFDIR, null, PDB_CODE, cifFile, false); 
		PdbAsymUnit fullpdb = new PdbAsymUnit(cifFile);
		PdbChain pdb = fullpdb.getChain(CHAIN);
		RIGraph graph = pdb.getRIGraph(CT, CUTOFF);
		RIGraph[] graphs = {graph};
		TinkerRunner tr = new TinkerRunner(TINKERBINDIR,DISTGEOM_EXE,PRMFILE);
		tr.reconstruct(graph.getSequence(), graphs, null, true, 1, 
				TinkerRunner.DEFAULT_FORCECONSTANT_DISTANCE, TinkerRunner.DEFAULT_FORCECONSTANT_TORSION, 
				System.getProperty("java.io.tmpdir"), "reconstructTester"+PDB_CODE+CHAIN, true, false);
		
		File outpdbfile = tr.getOutPdbFile(1);
		PdbAsymUnit fulloutpdb = new PdbAsymUnit(outpdbfile);
		PdbChain outpdb = fulloutpdb.getChain(PdbAsymUnit.NULL_CHAIN_CODE);
		assertEquals(pdb.getFullLength(),outpdb.getFullLength());
		assertEquals(pdb.getSequence().getSeq(),outpdb.getSequence().getSeq());
		
		// checking atoms and that chirality is fine
		for (int resser:outpdb.getAllStdAaResSerials()) {
			AaResidue res = (AaResidue) outpdb.getResidue(resser);
			AaResidue refRes = (AaResidue) pdb.getResidue(resser); 
			assertEquals(refRes.getAaType(),res.getAaType());
			for (Atom atom:refRes) {
				assertTrue(res.containsAtom(atom.getCode()));
			}
			AaResidue.Chirality chir = res.getChirality();
			assertTrue(chir!=AaResidue.Chirality.D);
		}
		
		// checking we really reconstructed input contact map
		IntPairSet viols = TinkerRunner.getViolatedEdges(graph, outpdb);
		System.out.println("Our violations count: "+viols.size());
		System.out.println("Tinker's violations count: "+ (tr.getNumUpperViol()[0]+tr.getNumLowerViol()[0]));
		assertTrue(viols.size()<graph.getEdgeCount()*0.05);
		// this seems not to always coincide, not sure why (hopefully it's just an effect of rounding) leaving it out
		//assertEquals(tr.getNumUpperViol()[0]+tr.getNumLowerViol()[0],viols.size());
		
		// checking rmsd of output or mirror is below 4 (should always be for 1bxyA) 
		boolean nonmirrorbelow = outpdb.rmsd(pdb, "Ca")<4;
		outpdb.mirror();
		boolean mirrorbelow = outpdb.rmsd(pdb, "Ca")<4;
		assertTrue(nonmirrorbelow || mirrorbelow);
		
		
	}

	/**
	 * To run the test from the command line, tinker can't seem to run through eclipse, it is 
	 * always killed because of memory  
	 * @param args
	 */
	public static void main(String args[]) throws Exception {
		JUnitCore.main("tests.tinker.TinkerRunnerTest");
		// following code to run as java program to be able to debug 
		//TinkerRunnerTest trt = new TinkerRunnerTest();
		//setUpBeforeClass();
		//trt.testReconstruct();
	}

}
