package owl.tests.core.runners.blast;


import java.io.File;
import java.util.Iterator;
import java.util.Properties;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import owl.core.runners.blast.BlastHit;
import owl.core.runners.blast.BlastHitList;
import owl.core.runners.blast.BlastHsp;
import owl.core.runners.blast.BlastRunner;
import owl.core.runners.blast.BlastTabularParser;
import owl.core.runners.blast.BlastXMLParser;
import owl.tests.TestsSetup;



/**
 * JUnit4 test for Blast parsers: will test our blast tab parser and blast xml parser comparing results.
 * 
 * @author duarte
 *
 */
public class BlastParsersTest {

	private static final String[] testSetCasp7 = 
	{"T0346", "T0290", "T0340", "T0345", "T0315", "T0305", "T0366", "T0295", "T0317", "T0288"};
	private static final String CASP7_TARGETS_DIR = "src/owl/tests/core/structure/data/casp7";
	
	private static String BLAST_BIN_DIR;
	private static String BLAST_DB_DIR;
	private static String BLAST_DB;
	
	private static final String TMP_DIR = System.getProperty("java.io.tmpdir");
	
	int numThreads = Runtime.getRuntime().availableProcessors(); // we run blast with as many CPUs we have available in executing host
	
	File[] tabFiles;
	File[] xmlFiles;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		Properties p = TestsSetup.readPaths();
		BLAST_BIN_DIR = p.getProperty("BLAST_BIN_DIR");
		BLAST_DB_DIR = p.getProperty("BLAST_DB_DIR");
		BLAST_DB = p.getProperty("BLAST_DB");
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		System.out.println("Blasting");
		tabFiles = new File[testSetCasp7.length];
		xmlFiles = new File[testSetCasp7.length];
		int i=0;
		for (String target: testSetCasp7) {
			File queryFile = new File(CASP7_TARGETS_DIR,target+".fa");
			File outFileXML = new File(TMP_DIR,target+".blast.xml");
			File outFileTAB = new File(TMP_DIR,target+".blast.tab");
			tabFiles[i] = outFileTAB;
			xmlFiles[i] = outFileXML;
			BlastRunner br = new BlastRunner(BLAST_BIN_DIR,BLAST_DB_DIR);
			br.runBlastp(queryFile, BLAST_DB, outFileXML, 7, true, numThreads);
			br.runBlastp(queryFile, BLAST_DB, outFileTAB, 8, true, numThreads);
			i++;
		}
	}

	@After
	public void tearDown() throws Exception {
		System.out.println("Removing files");
		for (int i=0;i<tabFiles.length;i++) {
			tabFiles[i].delete();
			xmlFiles[i].delete();
		}
	}

	@Test
	public void testXMLAgainstTAB() throws Exception {

		for (int i=0;i<tabFiles.length;i++) {
			System.out.println(tabFiles[i]+" vs "+xmlFiles[i]);
			
			// comparing XML against TAB parser
			BlastTabularParser btp = new BlastTabularParser(tabFiles[i]);
			BlastHitList hitListTAB = btp.getHits();
			
			BlastXMLParser bxm = new BlastXMLParser(xmlFiles[i]);
			BlastHitList hitListXML = bxm.getHits();
			
			// it seems that the tabular output lists a maximum of 250 hsp lines, thus the following
			if (hitListXML.size()>250) {
				Assert.assertTrue(hitListTAB.size()==250);
			} else {
				Assert.assertEquals(hitListTAB.size(), hitListXML.size());
			}
			
			int hspCount = 0;
			Iterator<BlastHit> it = hitListTAB.iterator();
			for (BlastHit hitXML: hitListXML) {
				System.out.println("  "+hitXML.getSubjectId());
				// it seems that the tabular output lists a maximum of 250 hsp lines, thus the following
				if (hspCount>=250) break;
				BlastHit hitTAB = it.next();
				Iterator<BlastHsp> hspIt= hitTAB.iterator();
				Assert.assertEquals(hitTAB.getSubjectId(),hitXML.getSubjectId());
				for (BlastHsp hspXML:hitXML) {
					hspCount++;
					BlastHsp hspTAB = hspIt.next();
					Assert.assertEquals(hspTAB.getAliLength(), hspXML.getAliLength());
					assertEvalues(hspTAB.getEValue(),hspXML.getEValue());
					Assert.assertEquals(hspTAB.getPercentIdentity(),hspXML.getPercentIdentity(),0.01);
					Assert.assertEquals(hspTAB.getQueryEnd(),hspXML.getQueryEnd());
					Assert.assertEquals(hspTAB.getQueryStart(),hspXML.getQueryStart());
					Assert.assertEquals(hspTAB.getScore(),hspXML.getScore(),1);
					Assert.assertEquals(hspTAB.getSubjectEnd(),hspXML.getSubjectEnd());
					Assert.assertEquals(hspTAB.getSubjectStart(),hspXML.getSubjectStart());

					// alignment
					Assert.assertNull(hspTAB.getAlignment());
					Assert.assertEquals(2,hspXML.getAlignment().getNumberOfSequences());

					String queryId = hitXML.getQueryId();
					String subjectId = hitXML.getSubjectId();
					int n = 1;
					for (String tag:hspXML.getAlignment().getTags()) {
						if (n==1) Assert.assertEquals(queryId, tag);
						if (n==2) Assert.assertEquals(subjectId, tag);
						Assert.assertNotNull(hspXML.getAlignment().getAlignedSequence(tag));
						n++;
					}

					// sanity checks
					Assert.assertTrue(hitXML.getQueryLength()>=hspXML.getQueryEnd());
					Assert.assertTrue(hitXML.getQueryLength()>=hspXML.getQueryStart());
					Assert.assertTrue(hitXML.getSubjectLength()>=hspXML.getSubjectEnd());
					Assert.assertTrue(hitXML.getSubjectLength()>=hspXML.getSubjectStart());
					Assert.assertTrue(hitXML.getQueryLength()>=hspXML.getQueryEnd()-hspXML.getQueryStart());
					Assert.assertTrue(hitXML.getSubjectLength()>=hspXML.getSubjectEnd()-hspXML.getSubjectStart());
					Assert.assertTrue(hspXML.getAliLength()<hitXML.getQueryLength()+hitXML.getSubjectLength());
					Assert.assertEquals(hspXML.getAlignment().getAlignmentLength(), hspXML.getAliLength());
				}
			}
		}
	}
	
	private void assertEvalues(double d1, double d2) {
		String d1str = String.valueOf(d1);
		String d2str = String.valueOf(d2);
		if (d1str.contains("E")) {
			int exp1=Integer.parseInt(d1str.substring(d1str.indexOf("E")+1, d1str.length()));
			int exp2=Integer.parseInt(d2str.substring(d2str.indexOf("E")+1, d2str.length()));
			if (exp1==exp2) {
				double b1 = Double.parseDouble(d1str.substring(0, d1str.indexOf("E")));
				double b2 = Double.parseDouble(d2str.substring(0, d2str.indexOf("E")));
				Assert.assertEquals(b1,b2,0.5);
			} else {
				Assert.assertEquals(exp1,exp2,1);
			}
		} else {
			Assert.assertEquals(d1,d2,0.5);
		}
	}
	
	// to debug the testing code (run as java program so that we can use normal debugger)
	public static void main(String[] args) throws Exception {
		BlastParsersTest tester = new BlastParsersTest();
		setUpBeforeClass();
		tester.setUp();
		tester.testXMLAgainstTAB();
	}

}
