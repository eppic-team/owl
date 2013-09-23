package owl.tests.core.structure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Properties;

import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

import org.junit.Assert;
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
import owl.core.structure.CrystalCell;
import owl.core.structure.HetResidue;
import owl.core.structure.InterfacesFinder;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbLoadException;
import owl.core.structure.Residue;
import owl.core.structure.SpaceGroup;
import owl.core.util.FileFormatException;
import owl.tests.TestsSetup;


public class PdbAsymUnitTest {

	@SuppressWarnings("unused") // we keep it here since we might still want to use naccess for area calculations
	private static String NACCESS_EXEC; 
	
	private static String LOCAL_CIF_DIR;
	
	private static final String TESTDATADIR = "src/owl/tests/core/structure/data";
	private static final String LISTFILE = TESTDATADIR+"/testset_interfaces.txt";
	private static final String LISTFILE2 = TESTDATADIR+"/testset_interfaces2.txt";
	private static final String CULLPDB20FILE = TESTDATADIR+"/cullpdb_20";
	
	private static final double CUTOFF = 5.9;
	
	// we allow for a 20% discrepancy from PISA in area values (we calculate with NACCESS/our own implementation and results will disagree always)
	private static final double TOLERANCE_ASA = 0.30;
	private static final double TOLERANCE_BSA = 0.20;
	// at least so many residues have to be in agreement within TOLERANCE above
	private static final double TOLERANCE_RESIDUE_AGREEMENT = 0.90;

	private static final boolean CONSIDER_HETATOMS = true;
	private static final boolean CONSIDER_NONPOLY = false;
	private static final int     CONSIDER_COFACTORS = -1;
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
		List<String> pdbCodes = readListFile(new File(LISTFILE));
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
	public void testInterfacesVsPisa() throws IOException, SQLException, SAXException {

		List<String> pdbCodes = readListFile(new File(LISTFILE));

		System.out.println("Interface calculation vs PISA test ("+pdbCodes.size()+" structures to test)");
		System.out.println("Will use "+NTHREADS+" CPUs for ASA calculations");
		
		// getting PISA interfaces
		PisaConnection pc = new PisaConnection();
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
			
			pdb.removeHatoms();

			ChainInterfaceList pisaInterfaces = all.get(pdbCode).convertToChainInterfaceList(pdb);
			// we sort them on interface area because pisa doesn't always sort them like that (it does some kind of grouping)
			pisaInterfaces.sort();
			
			System.out.println(pdb.getSpaceGroup().getShortSymbol()+" ("+pdb.getSpaceGroup().getId()+")");
			
			long start = System.currentTimeMillis();
			//ChainInterfaceList interfaces = pdb.getAllInterfaces(CUTOFF, new File(NACCESS_EXEC));
			ChainInterfaceList interfaces = pdb.getAllInterfaces(CUTOFF, Asa.DEFAULT_N_SPHERE_POINTS, NTHREADS, CONSIDER_HETATOMS, CONSIDER_NONPOLY, CONSIDER_COFACTORS);
			long end = System.currentTimeMillis();
			System.out.println("Time: "+((end-start)/1000)+"s");
			System.out.println("Total number of interfaces found: "+interfaces.size());
			
			int pisaCount = pisaInterfaces.getNumProtProtInterfaces();
			System.out.println("PISA interface count: "+pisaCount);
			Assert.assertEquals(pisaCount, interfaces.size());

			int i = 0;
			for (int p=0;p<pisaInterfaces.size();p++) {
				ChainInterface pisaInterf = pisaInterfaces.get(p+1);
				if (!pisaInterf.isProtein()) continue;
				System.out.println("\nInterface "+(i+1));
				ChainInterface myInterf = interfaces.get(i+1);
				
				
				System.out.printf("Areas, pisa: %8.2f\tmy: %8.2f\tdiff: %4.1f%%\n",pisaInterf.getInterfaceArea(),myInterf.getInterfaceArea(),
						(pisaInterf.getInterfaceArea()-myInterf.getInterfaceArea())*100.0/pisaInterf.getInterfaceArea());
				Assert.assertEquals(pisaInterf.getInterfaceArea(), myInterf.getInterfaceArea(), pisaInterf.getInterfaceArea()*0.10);
				
				// make sure there are no clashes
				Assert.assertFalse(myInterf.hasClashes());
				
				// asa/bsas of individual residues, we allow for some discrepancy from PISA
				
				PdbChain myFirstMol = myInterf.getFirstMolecule();
				if (!myFirstMol.getPdbChainCode().equals(pisaInterf.getFirstMolecule().getPdbChainCode())) {
					myFirstMol = myInterf.getSecondMolecule();
				}
				PdbChain mySecondMol = myInterf.getSecondMolecule();
				if (!mySecondMol.getPdbChainCode().equals(pisaInterf.getSecondMolecule().getPdbChainCode())) {
					mySecondMol = myInterf.getFirstMolecule();
				}

				//Assert.assertEquals(pisaInterf.getSecondTransf().getTransformId(),myInterf.getSecondTransf().getTransformId());
				//Assert.assertEquals(pisaInterf.getSecondTransf().getCrystalTranslation(),myInterf.getSecondTransf().getCrystalTranslation());
				
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
	
	@Test
	public void testInterfacesVsPisaCountsOnly() throws IOException, SQLException, SAXException {

		List<String> pdbCodes = readListFile(new File(LISTFILE2));

		System.out.println("Interface calculation vs PISA test ("+pdbCodes.size()+" structures to test). Only checking total interface counts");
		System.out.println("Will use "+NTHREADS+" CPUs for ASA calculations");
		
		// getting PISA interfaces
		PisaConnection pc = new PisaConnection();
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
			
			pdb.removeHatoms();

			ChainInterfaceList pisaInterfaces = all.get(pdbCode).convertToChainInterfaceList(pdb);
			// we sort them on interface area because pisa doesn't always sort them like that (it does some kind of grouping)
			pisaInterfaces.sort();
			
			System.out.println(pdb.getSpaceGroup().getShortSymbol()+" ("+pdb.getSpaceGroup().getId()+")");
			
			long start = System.currentTimeMillis();
			ChainInterfaceList interfaces = pdb.getAllInterfaces(CUTOFF, Asa.DEFAULT_N_SPHERE_POINTS, NTHREADS, CONSIDER_HETATOMS, CONSIDER_NONPOLY, CONSIDER_COFACTORS);
			long end = System.currentTimeMillis();
			System.out.println("Time: "+((end-start)/1000)+"s");
			System.out.println("Total number of interfaces found above 30A2 area: "+interfaces.getNumInterfacesAboveArea(30));
			
			int pisaCount = pisaInterfaces.getNumProtProtInterfacesAboveArea(30);
			System.out.println("PISA interface count above 30A2 area: "+pisaCount);
			Assert.assertEquals(pisaCount, interfaces.getNumInterfacesAboveArea(30));

			
			
		}
	}
	
	@Test
	public void testInterfacesFinderRedundancyElimination() throws IOException, SQLException, SAXException {
		
		// A exclude list for entries we don't want to test here
		// At the moment entries in exclude list are the ones for which a different number of interfaces is found
		// because slightly different number of contacts are found for symmetry-equivalent interfaces due to rounding
		// leading thus to apparent different interfaces since we base uniqueness of interfaces in same exact number of contacts between same atoms
		// e.g. in 1vyi there are 9 interfaces found if no symmetry redundancy elimination is used (they are actually only 6 unique ones following pisa)
		// The repeated ones are the pairs 1-2, 3-4 and 6-7. For instance in 6-7 atoms CD (residue 7) and O (residue 43) are at cutoff limit 5.9
		// This depends of course on cut-off values used. 
		String[] excludeCodesCutoffIssues = {"1vyi"};
		
		
		List<String> pdbCodes = readListFile(new File(CULLPDB20FILE));
		
		for (String excludeCode:excludeCodesCutoffIssues) {
			if (pdbCodes.remove(excludeCode)) {
				System.out.println("Removing code "+excludeCode+" because it is in the exclude list");
			}
		}
		
		
		System.out.println("Interface calculation - redundancy elimination test ("+pdbCodes.size()+" structures to test)");
		System.out.println("Will use "+NTHREADS+" CPUs for ASA calculations");
		
		ArrayList<String> withDuplicatesInRedundElim = new ArrayList<String>();
		
		for (String pdbCode: pdbCodes) {
					
			System.out.println("\n##"+pdbCode);
			File cifFile = new File(System.getProperty("java.io.tmpdir"),pdbCode+".pdbasymunittest.cif");
			try {
				PdbAsymUnit.grabCifFile(LOCAL_CIF_DIR, null, pdbCode, cifFile, false);
			} catch (IOException e) {
				System.err.println("Something went wrong while grabbing cif file: "+e.getMessage());
				continue;
			}

			
			// in some cases like 2pyq if one reuses the pdb for calculating with/without redund elim, the result is different (in 2pyq 21 vs 22 interfaces)
			// I guess there is some issue with the caching of bounding boxes
			// anyway to solve this we can simply clone the pdb

			PdbAsymUnit pdb = null;
			PdbAsymUnit pdb2 = null;
			try {
				pdb = new PdbAsymUnit(cifFile);				
			} catch (PdbLoadException e) {
				System.err.println("PDB load error, cause: "+e.getMessage());
				continue;
			} catch (FileFormatException e) {
				System.err.println("PDB load error, cause: "+e.getMessage());
				continue;
			}

			pdb.removeHatoms();
			pdb2 = pdb.copy();
			
			SpaceGroup sg = pdb.getSpaceGroup();
			if (sg==null) System.out.println("No space group");
			else System.out.println(sg.getShortSymbol()+" ("+sg.getId()+")");
			System.out.println(pdb.getNumPolyChains()+" polymer chains in AU");
			
			long start = 0, end =0;
			
			start = System.currentTimeMillis();
			
			InterfacesFinder interfFinder = new InterfacesFinder(pdb);
			interfFinder.setWithRedundancyElimination(false);
			ChainInterfaceList interfacesWithRedundancy = 
					interfFinder.getAllInterfaces(CUTOFF, Asa.DEFAULT_N_SPHERE_POINTS, NTHREADS, CONSIDER_HETATOMS, CONSIDER_NONPOLY, CONSIDER_COFACTORS);
			end = System.currentTimeMillis();
			long totalWithRedundancy = (end-start)/1000;
			
			InterfacesFinder interfFinder2 = new InterfacesFinder(pdb2);
			start = System.currentTimeMillis();
			interfFinder2.setWithRedundancyElimination(true);
			ChainInterfaceList interfacesRedundancyElim = 
					interfFinder2.getAllInterfaces(CUTOFF,Asa.DEFAULT_N_SPHERE_POINTS, NTHREADS, CONSIDER_HETATOMS, CONSIDER_NONPOLY, CONSIDER_COFACTORS);
			end = System.currentTimeMillis();
			long totalRedundancyElim = (end-start)/1000;			
			
			System.out.println("Time: ");
			System.out.println(" with redundancy: "+totalWithRedundancy+"s");
			System.out.println(" redundancy elimination: "+totalRedundancyElim+"s");
			
			System.out.println("Total number of interfaces found: "+interfacesWithRedundancy.size()+" "+interfacesRedundancyElim.size());
						
			System.out.println("Duplicates: "+
			 interfFinder.getDuplicatesCount1()+","+interfFinder.getDuplicatesCount2()+","+interfFinder.getDuplicatesCount3()+
			 "  -  "+
			 interfFinder2.getDuplicatesCount1()+","+interfFinder2.getDuplicatesCount2()+","+interfFinder2.getDuplicatesCount3());
			if (interfFinder2.getDuplicatesCount1()+interfFinder2.getDuplicatesCount2()+interfFinder2.getDuplicatesCount3()>0) {
				System.out.println(" warning! duplicates in redundancy elimination");
				withDuplicatesInRedundElim.add(pdbCode);
			}

			
			Assert.assertEquals(interfacesWithRedundancy.size(), interfacesRedundancyElim.size());

			System.out.println(("Areas: "));
			int countDiscAbove2pc = 0;
			for (int i=0;i<interfacesRedundancyElim.size();i++) {
				ChainInterface wr = interfacesWithRedundancy.get(i+1);
				ChainInterface re = interfacesRedundancyElim.get(i+1);
				
				//System.out.printf("%7.2f %7.2f - %5d %5d\n",wr.getInterfaceArea(),re.getInterfaceArea(),wr.getNumContacts(),re.getNumContacts());

				// because same interfaces can be found in different places in redundancy/non-redundancy
				// it can happen that the area values are slightly different (due to asa calc algorithm giving different results for same molec in different orientations)
				// then the sorting based on areas can be different... but at most there would be a different in 1 position, 
				// e.g. 1hbn
				// actually some weird cases do have differences in 2 positions! e.g. 1idp (because of 3 interfaces with very similar area)
				// in any case, instead of checking number of contacts, we check areas and only allow discrepancies within 10%
				// and additionally that there are no more than 2 values that differ in more than 2% (because the majority should be 0% differ)
				
				System.out.printf(" wr: %8.2f\tre: %8.2f\tdiff: %4.1f%%\n",
						wr.getInterfaceArea(),re.getInterfaceArea(),
						(wr.getInterfaceArea()-re.getInterfaceArea())*100.0/wr.getInterfaceArea());

				// we don't check the too small ones because they are problematic (and anyway not interesting)
				if (wr.getInterfaceArea()>10) {
					Assert.assertEquals(wr.getInterfaceArea(), re.getInterfaceArea(), wr.getInterfaceArea()*0.10); // 10% discrepancy permitted
					// additionally we count the ones with discrepancies above 3%, we won't allow more than 2
					if (Math.abs(wr.getInterfaceArea() - re.getInterfaceArea())>0.03*wr.getInterfaceArea()) {
						countDiscAbove2pc++;
					}
				}
			}
			Assert.assertTrue(countDiscAbove2pc<=2);
		
			// checking that we have the same set of number of contacts in both (it has to match 100%, there's no ambiguities like with areas)
			ArrayList<Integer> contactsWR = new ArrayList<Integer>();
			ArrayList<Integer> contactsRE = new ArrayList<Integer>();
			for (int i=0;i<interfacesRedundancyElim.size();i++) {
				ChainInterface wr = interfacesWithRedundancy.get(i+1);
				ChainInterface re = interfacesRedundancyElim.get(i+1);
				contactsWR.add(wr.getNumContacts());
				contactsRE.add(re.getNumContacts());
			} 
			Collections.sort(contactsWR);
			Collections.sort(contactsRE);
			
			for (int i=0;i<contactsWR.size();i++) {
				Assert.assertEquals(contactsWR.get(i),contactsRE.get(i));
			}
			
		}
		if (!withDuplicatesInRedundElim.isEmpty()) {
			System.out.println("Warning! "+withDuplicatesInRedundElim.size()+" structures had duplicate interfaces in redundancy elimination calculation: ");
			for (String pdbCode:withDuplicatesInRedundElim) {
				System.out.println(pdbCode);
			}
			System.out.println("\nDuplicates in redundancy elimination procedure are unusual but they can still occur (due mainly to NCS). ");			
		}

	}
	
	@Test
	public void testCrystalCellTransformations() throws IOException {
		List<String> pdbCodes = readListFile(new File(LISTFILE));

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
			
			pdb.removeHatoms();
			
			SpaceGroup sg = pdb.getSpaceGroup();
			CrystalCell cell = pdb.getCrystalCell();
			
			for (int op=0;op<sg.getNumOperators();op++) {
				Matrix4d transfXtal = sg.getTransformation(op);				
				
				Matrix4d transfOrthon = cell.transfToOrthonormal(transfXtal);
				
				Matrix4d transfBackToXtal = cell.transfToCrystal(transfOrthon);
				
				Matrix4d transfBackToOrthon = cell.transfToOrthonormal(transfBackToXtal);
				
				Assert.assertTrue(transfXtal.epsilonEquals(transfBackToXtal, 0.00001));
				
				Assert.assertTrue(transfOrthon.epsilonEquals(transfBackToOrthon, 0.00001));
				
				
				// transforms on translations done with Matrix4d or Tuple3d transformations must coincide
				Vector3d xtalTransl = new Vector3d(2,5,1);
				transfXtal.setTranslation(xtalTransl); 
				
				transfOrthon = cell.transfToOrthonormal(transfXtal);
				
				Vector3d orthonTransl = new Vector3d(xtalTransl);
				cell.transfToOrthonormal(orthonTransl); 
				
				Vector3d translFromMatrix4dTransform = new Vector3d(transfOrthon.m03, transfOrthon.m13, transfOrthon.m23);
				Assert.assertTrue(translFromMatrix4dTransform.epsilonEquals(orthonTransl, 0.0001)); 
				
			}
			
			int numCells = 3;
			for (int i=-numCells;i<=numCells;i++) {
				for (int j=-numCells;j<=numCells;j++) {
					for (int k=-numCells;k<=numCells;k++) {
						if (i==0 && j==0 && k==0) continue; 
						
						Vector3d translXtal = new Vector3d (i,j,k);
						Vector3d translOrthon = new Vector3d(translXtal);
						cell.transfToOrthonormal(translOrthon);
						Vector3d translBackToXtal = new Vector3d(translOrthon);
						cell.transfToCrystal(translBackToXtal);

						Assert.assertTrue(translXtal.epsilonEquals(translBackToXtal, 0.00001));
					}
				}
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
	
	private static List<String> readListFile(File file) throws IOException {
		List<String> pdbCodes = new ArrayList<String>();
		BufferedReader flist = new BufferedReader(new FileReader(file));
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
		pdbAUTest.testInterfacesVsPisa();
	}

}
