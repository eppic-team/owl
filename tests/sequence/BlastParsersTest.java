package tests.sequence;


import java.io.File;
import java.util.Iterator;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import sequence.BlastHit;
import sequence.BlastHitList;
import sequence.BlastRunner;
import sequence.BlastTabularParser;
import sequence.BlastXMLParser;


/**
 * JUnit4 test for Blast parsers: will test our blast tab parser and blast xml parser comparing results.
 * 
 * @author duarte
 *
 */
public class BlastParsersTest {

	int numThreads = Runtime.getRuntime().availableProcessors(); // we run blast with as many CPUs we have available in executing host
	
	String[] testSetCasp7 = {"T0346", "T0290", "T0340", "T0345", "T0315", "T0305", "T0366", "T0295", "T0317", "T0288"};
	String targetDir = "/project/StruPPi/CASP7/targets";
	String outDir = "/scratch/local/temp";
	String db = "seqs_pdbase_20080903.fix.reps.fa";
	
	File[] tabFiles;
	File[] xmlFiles;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
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
			File queryFile = new File(targetDir,target+".fa");
			File outFileXML = new File(outDir,target+".blast.xml");
			File outFileTAB = new File(outDir,target+".blast.tab");
			tabFiles[i] = outFileTAB;
			xmlFiles[i] = outFileXML;
			BlastRunner br = new BlastRunner("/project/StruPPi/bin","/project/StruPPi/CASP8/blast_dbs");
			br.runBlastp(queryFile, db, outFileXML, 7, true, numThreads);
			br.runBlastp(queryFile, db, outFileTAB, 8, true, numThreads);
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
			
			Assert.assertEquals(hitListTAB.size(), hitListXML.size());
			
			Iterator<BlastHit> it = hitListTAB.iterator();
			for (BlastHit hitXML: hitListXML) {
				BlastHit hitTAB = it.next();
				System.out.println("  "+hitXML.getSubjectId());
				Assert.assertEquals(hitTAB.getSubjectId(),hitXML.getSubjectId());
				Assert.assertEquals(hitTAB.getAliLength(), hitXML.getAliLength());
				assertEvalues(hitTAB.getEValue(),hitXML.getEValue());
				Assert.assertEquals(hitTAB.getPercentIdentity(),hitXML.getPercentIdentity(),0.01);
				Assert.assertEquals(hitTAB.getQueryEnd(),hitXML.getQueryEnd());
				Assert.assertEquals(hitTAB.getQueryStart(),hitXML.getQueryStart());
				Assert.assertEquals(hitTAB.getScore(),hitXML.getScore(),1);
				Assert.assertEquals(hitTAB.getSubjectEnd(),hitXML.getSubjectEnd());
				Assert.assertEquals(hitTAB.getSubjectStart(),hitXML.getSubjectStart());
				
				// alignment
				Assert.assertNull(hitTAB.getAlignment());
				Assert.assertEquals(2,hitXML.getAlignment().getNumberOfSequences());
				
				String queryId = hitXML.getQueryId();
				String subjectId = hitXML.getSubjectId();
				int n = 1;
				for (String tag:hitXML.getAlignment().getTags()) {
					if (n==1) Assert.assertEquals(queryId, tag);
					if (n==2) Assert.assertEquals(subjectId, tag);
					Assert.assertNotNull(hitXML.getAlignment().getAlignedSequence(tag));
					n++;
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
}
