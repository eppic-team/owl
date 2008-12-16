package tests.proteinstructure;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import junit.framework.Assert;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import proteinstructure.Pdb;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbfilePdb;

public class PdbTest {

	private static final String TEST_PDB_FILE = "tests/proteinstructure/data/1tdrA.pdb";
	private static final String TEST_CHAIN = "A";
	private static final String NACCESS_EXEC = "/project/StruPPi/Software/naccess2.1.1/naccess";
	private static final String NACCESS_OUTPUT_REF = "tests/proteinstructure/data/1tdrA.rsa";
	
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
		
		Pdb pdb = new PdbfilePdb(TEST_PDB_FILE);
		pdb.load(TEST_CHAIN);
		pdb.runNaccess(NACCESS_EXEC, "");
		
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
		Assert.assertEquals(resser,pdb.get_length());
		
		
	}
	
}
