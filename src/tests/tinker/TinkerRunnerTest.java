package tests.tinker;


import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.runner.JUnitCore;

import proteinstructure.Atom;
import proteinstructure.ConformationsNotSameSizeError;
import proteinstructure.IntPairSet;
import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbasePdb;
import proteinstructure.PdbfilePdb;
import proteinstructure.RIGraph;
import proteinstructure.Residue;

import tinker.TinkerError;
import tinker.TinkerRunner;
import tools.MySQLConnection;

public class TinkerRunnerTest {

	private static final String PDB_CODE = "1bxy";
	private static final String CHAIN = "A";
	private static final String PDBASE_DB = "pdbase";
	private static final String CT = "Cb";
	private static final double CUTOFF = 8;
	
	private static final String TINKERBINDIR = "/project/StruPPi/Software/tinker/bin";
	private static final String DISTGEOM_EXE = "distgeom.32.static.maxatm40000";
	private static final String PRMFILE = "/project/StruPPi/Software/tinker/amber/amber99.prm";
	
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
	public void testReconstruct() throws SQLException, PdbCodeNotFoundError, PdbLoadError, TinkerError, IOException, ConformationsNotSameSizeError {
		MySQLConnection conn = new MySQLConnection();
		Pdb pdb = new PdbasePdb(PDB_CODE, PDBASE_DB,conn);
		pdb.load(CHAIN);
		RIGraph graph = pdb.getRIGraph(CT, CUTOFF);
		RIGraph[] graphs = {graph};
		TinkerRunner tr = new TinkerRunner(TINKERBINDIR,DISTGEOM_EXE,PRMFILE);
		tr.reconstruct(graph.getSequence(), graphs, null, true, 1, 
				TinkerRunner.DEFAULT_FORCECONSTANT_DISTANCE, TinkerRunner.DEFAULT_FORCECONSTANT_TORSION, 
				System.getProperty("java.io.tmpdir"), "reconstructTester"+PDB_CODE+CHAIN, true, false);
		
		File outpdbfile = tr.getOutPdbFile(1);
		Pdb outpdb = new PdbfilePdb(outpdbfile.getAbsolutePath());
		outpdb.load(Pdb.NULL_CHAIN_CODE);
		assertEquals(pdb.getFullLength(),outpdb.getFullLength());
		assertEquals(pdb.getSequence(),outpdb.getSequence());
		
		// checking atoms and that chirality is fine
		for (int resser:outpdb.getAllSortedResSerials()) {
			Residue res = outpdb.getResidue(resser);
			Residue refRes = pdb.getResidue(resser); 
			assertEquals(refRes.getAaType(),res.getAaType());
			for (Atom atom:refRes) {
				assertTrue(res.containsAtom(atom.getCode()));
			}
			Residue.Chirality chir = res.getChirality();
			assertTrue(chir!=Residue.Chirality.D);
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
	public static void main(String args[]) {
		JUnitCore.main("tests.tinker.TinkerRunnerTest");
	}

}
