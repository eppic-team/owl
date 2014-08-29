package owl.core.structure;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Properties;
import java.util.zip.GZIPInputStream;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import owl.core.util.FileFormatException;
import owl.tests.TestsSetup;

public class HbplusTest {
	
	private static String LOCAL_PDB_DIR;
	private static File HBPLUS_EXE;
	private static final String[] PDBIDS_TO_TEST = {"1d2s", "1lom", "1oew", "1smt", "1tdr", "2gpi", "2h6f", "3fsl", "3nul", "7odc"};

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		System.out.println("...SETTING UP HBPLUS TESTS...");
		Properties p = TestsSetup.readPaths();
		LOCAL_PDB_DIR = p.getProperty("LOCAL_PDB_DIR");
		HBPLUS_EXE = new File(p.getProperty("HBPLUS_EXE"));
		System.out.println("...TESTING HBPLUS...");
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		for (int i = 0; i < PDBIDS_TO_TEST.length; i++) {
			String pdbidToTest = PDBIDS_TO_TEST[i];
			File testPdbFileGz = new File(LOCAL_PDB_DIR + "/pdb" + pdbidToTest + ".ent.gz");
			byte[] buffer = new byte[1024];
			GZIPInputStream gZIPInputStream = new GZIPInputStream(new FileInputStream(testPdbFileGz.getPath()));
			FileOutputStream fileOutputStream = new FileOutputStream(pdbidToTest + ".pdb");
			int bytes_read;
			while ((bytes_read = gZIPInputStream.read(buffer)) > 0) {
				fileOutputStream.write(buffer, 0, bytes_read);
			}
			gZIPInputStream.close();
			fileOutputStream.close();
		}
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testHbplus() throws IOException, InterruptedException, FileFormatException, PdbLoadException {
		for (int i = 0; i < PDBIDS_TO_TEST.length; i++) {
			String pdbidToTest = PDBIDS_TO_TEST[i];
			PdbAsymUnit pdb = new PdbAsymUnit(new File(pdbidToTest + ".pdb"));
			pdb.removeHatoms();
			ChainInterfaceList interfaces = pdb.getAllInterfaces(5.5, 3000, 1, false, false, -1, 35);
			System.out.println("=====" + pdbidToTest.toUpperCase() + "=====");
			for (ChainInterface interf : interfaces) {
				if (interf.isFirstProtein() && interf.isSecondProtein()) {
					File pdbFile = new File(".", pdbidToTest + "." + interf.getId() + ".pdb.gz");
					interf.writeToPdbFile(pdbFile, true, true);
					if (HBPLUS_EXE != null && HBPLUS_EXE.exists()) {
						interf.runHBPlus(HBPLUS_EXE, pdbFile);
						System.out.println("HBPlus successfully run on Interface " + interf.getId() + ".");
					}
					else {
						System.out.println("Interface " + interf.getId() + ": Couldn't run HBPlus. Non-fatal error: an "
								+ "alternate method for hydrogen-bond prediction will be used.");
					}
				}
			}
		}
	}
}