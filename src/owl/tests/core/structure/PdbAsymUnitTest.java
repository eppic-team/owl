package owl.tests.core.structure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Properties;

//import javax.vecmath.Matrix4d;

import junit.framework.Assert;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.xml.sax.SAXException;

import owl.core.connections.pisa.PisaConnection;
import owl.core.structure.ChainInterface;
import owl.core.structure.ChainInterfaceList;
import owl.core.structure.Pdb;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbCodeNotFoundError;
import owl.core.structure.PdbLoadError;
import owl.core.structure.Residue;
import owl.core.util.MySQLConnection;
import owl.tests.TestsSetup;


public class PdbAsymUnitTest {

	private static String NACCESS_EXEC; 
	
	private static final String TESTDATADIR = "src/owl/tests/core/structure/data";
	private static final String LISTFILE = TESTDATADIR+"/testset_interfaces.txt";
	
	private static final String PDBASE_DB = "pdbase";
	
	private static final String PISA_INTERFACES_URL = "http://www.ebi.ac.uk/msd-srv/pisa/cgi-bin/interfaces.pisa?";
	
	private static final double CUTOFF = 5.9;

	// we allow for a 20% discrepancy from PISA in area values (we calculate with NACCESS and results will disagree always)
	private static final double TOLERANCE_ASA = 0.30;
	private static final double TOLERANCE_BSA = 0.20;
	// at least so many residues have to be in agreement within TOLERANCE above
	private static final double TOLERANCE_RESIDUE_AGREEMENT = 0.90;

	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		Properties p = TestsSetup.readPaths();
		NACCESS_EXEC = p.getProperty("NACCESS_EXEC");
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
	public void testGetAllInterfaces() throws IOException, SQLException, SAXException {
		MySQLConnection conn = new MySQLConnection();

		List<String> pdbCodes = new ArrayList<String>();
		BufferedReader flist = new BufferedReader(new FileReader(LISTFILE));
		String line;
		while ((line = flist.readLine() ) != null ) {
			if (line.startsWith("#")) continue;
			if (line.isEmpty()) break;
			String pdbCode = line.split("\\s+")[0].toLowerCase();
			pdbCodes.add(pdbCode);
		}
		flist.close();

		// getting PISA interfaces
		PisaConnection pc = new PisaConnection(PISA_INTERFACES_URL, null, null);
		System.out.println("Downloading PISA interfaces");
		Map<String, ChainInterfaceList> all = pc.getInterfacesDescription(pdbCodes);

		for (String pdbCode: pdbCodes) {
		
			ChainInterfaceList pisaInterfaces = all.get(pdbCode);
			// we sort them on interface area because pisa doesn't always sort them like that (it does some kind of grouping)
			pisaInterfaces.sort();
			
			System.out.println("\n##"+pdbCode);

			PdbAsymUnit pdb = null;
			try {
				pdb = new PdbAsymUnit(pdbCode, conn, PDBASE_DB);
			} catch (PdbCodeNotFoundError e) {
				System.err.println("Missing PDB code "+pdbCode);
				continue;
			} catch (PdbLoadError e) {
				System.err.println("PDB load error, cause: "+e.getMessage());
				continue;
			}

			System.out.println(pdb.getSpaceGroup().getShortSymbol()+" ("+pdb.getSpaceGroup().getId()+")");
			
			long start = System.currentTimeMillis();
			ChainInterfaceList interfaces = pdb.getAllInterfaces(CUTOFF, new File(NACCESS_EXEC));
			long end = System.currentTimeMillis();
			System.out.println("Time: "+((end-start)/1000)+"s");
			System.out.println("Total number of interfaces found: "+interfaces.size());
			
			int pisaCount = pisaInterfaces.getNumProtProtInterfaces();
			System.out.println("PISA interface count: "+pisaCount);
			Assert.assertEquals(pisaCount, interfaces.size());

			int i = 0;
			for (int p=0;p<pisaInterfaces.size();p++) {
				ChainInterface pisaInterf = pisaInterfaces.get(p);
				if (!pisaInterf.isProtein()) continue;
				System.out.println("\nInterface "+(i+1));
				ChainInterface myInterf = interfaces.get(i);
				
				
				System.out.printf("Areas, pisa: %8.2f\tmy: %8.2f\tdiff: %4.1f%%\n",pisaInterf.getInterfaceArea(),myInterf.getInterfaceArea(),
						(pisaInterf.getInterfaceArea()-myInterf.getInterfaceArea())*100.0/pisaInterf.getInterfaceArea());
				Assert.assertEquals(pisaInterf.getInterfaceArea(), myInterf.getInterfaceArea(), pisaInterf.getInterfaceArea()*0.10);
				
				// asa/bsas of individual residues, we allow for some discrepancy from PISA
				
				Pdb myFirstMol = myInterf.getFirstMolecule();
				if (!myFirstMol.getPdbChainCode().equals(pisaInterf.getFirstMolecule().getPdbChainCode())) {
					myFirstMol = myInterf.getSecondMolecule();
				}
				Pdb mySecondMol = myInterf.getSecondMolecule();
				if (!mySecondMol.getPdbChainCode().equals(pisaInterf.getSecondMolecule().getPdbChainCode())) {
					mySecondMol = myInterf.getFirstMolecule();
				}

				
				System.out.println("Chain 1");
				int[] counts1 = checkResidues(pisaInterf.getFirstMolecule(), myFirstMol);
				int[] counts2 = null;
				
				if (checkCounts(counts1)) {
					System.out.println("Chain 2");
					counts2 = checkResidues(pisaInterf.getSecondMolecule(), mySecondMol);

				} else {
					System.out.println("Counts of first PISA chain to our first didn't match. Trying swapping chains.");
					System.out.println("Chain 1");
					counts1 = checkResidues(pisaInterf.getSecondMolecule(), myFirstMol);
					System.out.println("Chain 2");
					counts2 = checkResidues(pisaInterf.getFirstMolecule(), mySecondMol);
				}
				
				if (!checkCounts(counts1) || !checkCounts(counts2)) {
					System.out.println("Failure for "+pdbCode+", interface "+(i+1));
				}
				Assert.assertTrue(checkCounts(counts1));
				Assert.assertTrue(checkCounts(counts2));

				i++;
			}
			
		}
	}
	
	private static boolean deltaComp(double a, double b, double delta) {
		boolean within = false;
		if (delta<0.2) { // for small values we have to have a bigger margin, chose 0.50 more or less arbitrarily
			within = (Math.abs(a-b)<0.50); 
		} else {
			within = (Math.abs(a-b)<=delta);
		}
		
		return within;
	}
	
	/**
	 * Checks individual residues' ASAs and BSAs values between pisa and our values
	 * @param pisaMolecule
	 * @param myMolecule
	 * @return
	 */
	private static int[] checkResidues(Pdb pisaMolecule, Pdb myMolecule) {
		int a = 0;
		int b = 0;
		int t = 0;
		for (Residue residue:pisaMolecule.getResidues().values()) {
			int resser = myMolecule.getResSerFromPdbResSer(residue.getPdbSerial());
			double pisaAsa = residue.getAsa();
			double pisaBsa = residue.getBsa();
			Residue myRes = myMolecule.getResidue(resser);
			if (myRes!=null) {
				double myAsa = myRes.getAsa();
				double myBsa = myRes.getBsa();
				System.out.printf("%s\t%s\t%d\t%s\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t",
						residue.getPdbSerial(),residue.getAaType().getThreeLetterCode(),resser,myRes.getAaType().getThreeLetterCode(),pisaAsa,pisaBsa,myAsa,myBsa);
				Assert.assertEquals(residue.getAaType(),myRes.getAaType());
				if (deltaComp(pisaAsa, myAsa, pisaAsa*TOLERANCE_ASA)) {
					System.out.print(" ");
					a++;
				} else {
					System.out.print("x");
				}
				if (deltaComp(pisaBsa, myBsa, pisaBsa*TOLERANCE_BSA)) {
					System.out.print(" ");
					b++;
				} else {
					System.out.print("x");
				}
				t++;
				System.out.println();
			}
		}
		System.out.println("Total: "+t+". Agreements within "+String.format("%4.2f(ASA) %4.2f(BSA)",TOLERANCE_ASA,TOLERANCE_BSA)+" tolerance: ASA "+a+" BSA "+b);
		int[] counts = {a,b,t};
		return counts;
	}

	/**
	 * Returns true if both asa and bsa counts are above the predefined tolerance
	 * threshold {@value #TOLERANCE_RESIDUE_AGREEMENT} for number of residues in agreement with PISA
	 * @param counts
	 * @return
	 */
	private static boolean checkCounts(int[] counts) {
		return (counts[0]>(TOLERANCE_RESIDUE_AGREEMENT*(double)counts[2]) && counts[1]>(TOLERANCE_RESIDUE_AGREEMENT*(double)counts[2]));
	}
}
