package owl.tests.core.structure;


import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import owl.core.structure.AAinfo;
import owl.core.structure.AminoAcid;


public class AminoAcidTest {

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
	public void testAminoAcid() {
		// test copied from main method in AminoAcid
    	// iterate over all 20 amino acids
		// perform circular conversion and verify result
		// make sure that wrong input will lead to desired result

		String three; char one; int num; AminoAcid result;
		
		//Testing conversion functions of class AminoAcid

		// Step 1
		for (AminoAcid aa : AminoAcid.values()) {
			num = aa.getNumber();
			result = AminoAcid.getByNumber(AminoAcid.three2num(AminoAcid.one2three(AminoAcid.num2one(num))));
			assertTrue(aa.equals(result));
		}
		
		// Step 2
		for (AminoAcid aa : AminoAcid.values()) {
			one = aa.getOneLetterCode();
			result = AminoAcid.getByOneLetterCode(AminoAcid.three2one(AminoAcid.num2three(AminoAcid.one2num(one))));
			assertTrue(aa.equals(result));
		}
		
		// Step 3
		for (AminoAcid aa : AminoAcid.values()) {	
			three = aa.getThreeLetterCode();
			result = AminoAcid.getByThreeLetterCode(AminoAcid.one2three(AminoAcid.num2one(AminoAcid.three2num(three))));
			assertTrue(aa.equals(result));
		}
		
		// Step 4
		assertNull(AminoAcid.getByNumber(-2)); 
		assertEquals('?',AminoAcid.num2one(-2));
		assertNull(AminoAcid.num2three(-2));

		assertNull(AminoAcid.getByOneLetterCode('?'));
		assertEquals(-2,AminoAcid.one2num('?'));
		assertNull(AminoAcid.one2three('?'));
		
		assertNull(AminoAcid.getByThreeLetterCode(""));
		assertEquals('?',AminoAcid.three2one(""));
		assertEquals(-2,AminoAcid.three2num(""));	
		
		// Step 5 : compare number of atoms with data in AAinfo
		for(AminoAcid aa : AminoAcid.values()) {
			// for unknown amino acids, the reported number of atoms should be -1
			if(aa == AminoAcid.XXX || aa == AminoAcid.STP) {
				assertEquals(-1,aa.getNumberOfAtoms());
			} else {	
				assertEquals(AAinfo.getNumberAtoms(aa.getThreeLetterCode()),aa.getNumberOfAtoms() + 4);
			}
			if (aa!=AminoAcid.XXX && aa!=AminoAcid.STP) {
				// check validation methods
				assertTrue(AminoAcid.isStandardAA(aa.getOneLetterCode()));
				assertTrue(AminoAcid.isStandardAA(aa.getThreeLetterCode()));
			}
		}
	}

}
