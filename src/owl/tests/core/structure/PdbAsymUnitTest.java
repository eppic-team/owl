package owl.tests.core.structure;

import java.io.BufferedReader;
//import java.io.File;
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
import owl.core.connections.pisa.PisaInterfaceList;
import owl.core.structure.Asa;
import owl.core.structure.ChainInterface;
import owl.core.structure.ChainInterfaceList;
import owl.core.structure.HetResidue;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbLoadException;
import owl.core.structure.Residue;
import owl.core.util.FileFormatException;
import owl.tests.TestsSetup;


public class PdbAsymUnitTest {

	@SuppressWarnings("unused") // we keep it here since we might still want to use naccess for area calculations
	private static String NACCESS_EXEC; 
	
	private static String LOCAL_CIF_DIR;
	
	private static final String TESTDATADIR = "src/owl/tests/core/structure/data";
	private static final String LISTFILE = TESTDATADIR+"/testset_interfaces.txt";
	
	private static final String PISA_INTERFACES_URL = "http://www.ebi.ac.uk/msd-srv/pisa/cgi-bin/interfaces.pisa?";
	
	private static final double CUTOFF = 5.9;
	
	private static final double CLASH_DISTANCE = 1.5;

	// we allow for a 20% discrepancy from PISA in area values (we calculate with NACCESS/our own implementation and results will disagree always)
	private static final double TOLERANCE_ASA = 0.30;
	private static final double TOLERANCE_BSA = 0.20;
	// at least so many residues have to be in agreement within TOLERANCE above
	private static final double TOLERANCE_RESIDUE_AGREEMENT = 0.90;

	private static final boolean CONSIDER_HETATOMS = true;
	private static final boolean CONSIDER_NONPOLY = false;
	private static final boolean PRINT_PER_RES = false; // whether to print areas agreement per residue or not
	
	private static final int NTHREADS = Runtime.getRuntime().availableProcessors(); // number of threads for ASA calculation

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		Properties p = TestsSetup.readPaths();
		NACCESS_EXEC = p.getProperty("NACCESS_EXEC");
		LOCAL_CIF_DIR = p.getProperty("LOCAL_CIF_DIR");
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
	public void testGetAllRepChains() throws IOException {
		List<String> pdbCodes = readListFile();
		System.out.println("Checking unique and representative chains");
		for (String pdbCode: pdbCodes) {
			
			System.out.println(pdbCode);
			File cifFile = new File(System.getProperty("java.io.tmpdir"),pdbCode+".pdbasymunittest.cif");
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
			List<String> repChains = pdb.getAllRepChains();
			Assert.assertTrue(repChains.size()<=pdb.getNumPolyChains() && repChains.size()>0);
			List<String> allchains = new ArrayList<String>();
			for (String repChain:repChains) {
				List<String> chains = pdb.getSeqIdenticalGroup(repChain);
				allchains.addAll(chains);
			}
			Assert.assertTrue(allchains.size()==pdb.getNumPolyChains());
			for (String pdbChainCode:pdb.getPdbChainCodes()) {
				Assert.assertTrue(repChains.contains(pdb.getRepChain(pdbChainCode)));
			}
			
		}
	}
	
	@Test
	public void testGetAllInterfaces() throws IOException, SQLException, SAXException {

		List<String> pdbCodes = readListFile();

		System.out.println("Will use "+NTHREADS+" CPUs for ASA calculations");
		
		// getting PISA interfaces
		PisaConnection pc = new PisaConnection(PISA_INTERFACES_URL, null, null);
		System.out.println("Downloading PISA interfaces");
		Map<String, PisaInterfaceList> all = pc.getInterfacesDescription(pdbCodes);

		for (String pdbCode: pdbCodes) {
					
			System.out.println("\n##"+pdbCode);
			File cifFile = new File(System.getProperty("java.io.tmpdir"),pdbCode+".pdbasymunittest.cif");
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

			ChainInterfaceList pisaInterfaces = all.get(pdbCode).convertToChainInterfaceList(pdb);
			// we sort them on interface area because pisa doesn't always sort them like that (it does some kind of grouping)
			pisaInterfaces.sort();
			
			System.out.println(pdb.getSpaceGroup().getShortSymbol()+" ("+pdb.getSpaceGroup().getId()+")");
			
			long start = System.currentTimeMillis();
			//ChainInterfaceList interfaces = pdb.getAllInterfaces(CUTOFF, new File(NACCESS_EXEC));
			ChainInterfaceList interfaces = pdb.getAllInterfaces(CUTOFF, null, Asa.DEFAULT_N_SPHERE_POINTS, NTHREADS, CONSIDER_HETATOMS, CONSIDER_NONPOLY, false);
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
				
				// make sure there are no clashes
				Assert.assertFalse(myInterf.hasClashes(CLASH_DISTANCE));
				
				// asa/bsas of individual residues, we allow for some discrepancy from PISA
				
				PdbChain myFirstMol = myInterf.getFirstMolecule();
				if (!myFirstMol.getPdbChainCode().equals(pisaInterf.getFirstMolecule().getPdbChainCode())) {
					myFirstMol = myInterf.getSecondMolecule();
				}
				PdbChain mySecondMol = myInterf.getSecondMolecule();
				if (!mySecondMol.getPdbChainCode().equals(pisaInterf.getSecondMolecule().getPdbChainCode())) {
					mySecondMol = myInterf.getFirstMolecule();
				}

				
				System.out.println("Chain 1");
				int[] counts1 = checkResidues(pisaInterf.getFirstMolecule(), myFirstMol, PRINT_PER_RES, CONSIDER_HETATOMS);
				int[] counts2 = null;
				
				if (checkCounts(counts1)) {
					System.out.println("Chain 2");
					counts2 = checkResidues(pisaInterf.getSecondMolecule(), mySecondMol, PRINT_PER_RES, CONSIDER_HETATOMS);

				} else {
					System.out.println("Counts of first PISA chain to our first didn't match. Trying swapping chains.");
					System.out.println("Chain 1");
					counts1 = checkResidues(pisaInterf.getSecondMolecule(), myFirstMol, PRINT_PER_RES, CONSIDER_HETATOMS);
					System.out.println("Chain 2");
					counts2 = checkResidues(pisaInterf.getFirstMolecule(), mySecondMol, PRINT_PER_RES, CONSIDER_HETATOMS);
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
	private static int[] checkResidues(PdbChain pisaMolecule, PdbChain myMolecule, boolean printPerRes, boolean hetAtom) {
		int a = 0;
		int b = 0;
		int t = 0;
		for (Residue residue:pisaMolecule) {
			int resser = myMolecule.getResSerFromPdbResSer(residue.getPdbSerial());
			double pisaAsa = residue.getAsa();
			double pisaBsa = residue.getBsa();
			Residue myRes = myMolecule.getResidue(resser);
			if (myRes!=null) {
				if (!hetAtom && (myRes instanceof HetResidue)) continue;
				double myAsa = myRes.getAsa();
				double myBsa = myRes.getBsa();
				if (printPerRes) {
					System.out.printf("%s\t%s\t%d\t%s\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t",
						residue.getPdbSerial(),residue.getLongCode(),resser,myRes.getLongCode(),pisaAsa,pisaBsa,myAsa,myBsa);
				}
				Assert.assertEquals(residue.getLongCode(),myRes.getLongCode());
				if (deltaComp(pisaAsa, myAsa, pisaAsa*TOLERANCE_ASA)) {
					if (printPerRes) System.out.print(" ");
					a++;
				} else {
					if (printPerRes) System.out.print("x");
				}
				if (deltaComp(pisaBsa, myBsa, pisaBsa*TOLERANCE_BSA)) {
					if (printPerRes) System.out.print(" ");
					b++;
				} else {
					if (printPerRes) System.out.print("x");
				}
				t++;
				if (printPerRes) System.out.println();
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
	
	private static List<String> readListFile() throws IOException {
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
		return pdbCodes;
	}
	
	// to debug the testing code (run as java program so that we can use normal debugger)
	public static void main(String[] args) throws Exception {
		PdbAsymUnitTest pdbAUTest = new PdbAsymUnitTest();
		setUpBeforeClass();
		pdbAUTest.testGetAllInterfaces();
	}

}
