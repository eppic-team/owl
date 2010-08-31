package owl.tests.core.structure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

//import javax.vecmath.Matrix4d;

import junit.framework.Assert;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.xml.sax.SAXException;

import owl.core.connections.pisa.PisaConnection;
import owl.core.connections.pisa.PisaInterface;
import owl.core.runners.NaccessRunner;
import owl.core.structure.ChainInterface;
import owl.core.structure.Pdb;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbCodeNotFoundError;
import owl.core.structure.PdbLoadError;
import owl.core.util.MySQLConnection;
import owl.tests.TestsSetup;


public class PdbAsymUnitTest {

	private static String NACCESS_EXEC; 
	
	private static final String TESTDATADIR = "src/owl/tests/core/structure/data";
	private static final String LISTFILE = TESTDATADIR+"/testset_interfaces.txt";
	
	private static final String PDBASE_DB = "pdbase";
	
	private static final String PISA_INTERFACES_URL = "http://www.ebi.ac.uk/msd-srv/pisa/cgi-bin/interfaces.pisa?";
	
	private static final double CUTOFF = 6.0;

	
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
			String pdbCode = line.split("\\s+")[0].toLowerCase();
			pdbCodes.add(pdbCode);
		}
		flist.close();

		// getting PISA interfaces
		PisaConnection pc = new PisaConnection(PISA_INTERFACES_URL, null, null);
		System.out.println("Downloading PISA interfaces");
		Map<String, List<PisaInterface>> all = pc.getInterfacesDescription(pdbCodes);

		for (String pdbCode: pdbCodes) {
		
			List<PisaInterface> pisaInterfaces = all.get(pdbCode);
			// we sort them on interface area because pisa doesn't always sort them like that (it does some kind of grouping)
			List<PisaInterface> pisaSortedInterfaces = new ArrayList<PisaInterface>();
			pisaSortedInterfaces.addAll(pisaInterfaces);
			Collections.sort(pisaSortedInterfaces);
			
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
			Set<ChainInterface> interfaces = pdb.getAllInterfaces(CUTOFF);
			long end = System.currentTimeMillis();
			System.out.println("Time: "+((end-start)/1000)+"s");
			System.out.println("Total number of interfaces found: "+interfaces.size());
			System.out.println("Calculating BSAs with NACCESS");

			for (Pdb chain:pdb.getAllChains()) {
				NaccessRunner nar = new NaccessRunner(new File(NACCESS_EXEC), ""); 
				nar.runNaccess(chain);
			}
			HashMap<String, HashMap<Integer,Double>> asas = pdb.getAbsSurfaceAccessibilities();

			for (ChainInterface interf:interfaces) {
				System.out.print(".");
				interf.calcBSAnaccess(new File(NACCESS_EXEC),asas);
			}
			System.out.println();

			List<ChainInterface> sortedInterfaces = new ArrayList<ChainInterface>();
			sortedInterfaces.addAll(interfaces);
			Collections.sort(sortedInterfaces);
			
			int pisaCount = getNumProtProtInterfaces(pisaInterfaces);
			System.out.println("PISA interface count: "+pisaCount);
			Assert.assertEquals(pisaCount, sortedInterfaces.size());

			int i = sortedInterfaces.size()-1;
			for (int p=pisaSortedInterfaces.size()-1;p>=0;p--) {
				PisaInterface pisaInterf = pisaSortedInterfaces.get(p);
				if (!pisaInterf.isProtein()) continue;
				System.out.println("\nInterface "+(sortedInterfaces.size()-i));
				ChainInterface myInterf = sortedInterfaces.get(i);
				
				
				//ChainInterface interf = findInterface(sortedInterfaces, pisaInterf);
				//Assert.assertNotNull(interf);
				
				
				// we allow for a 5% discrepancy from PISA
				System.out.printf("Areas, pisa: %8.2f my: %8.2f\n",pisaInterf.getInterfaceArea(),myInterf.getInterfaceArea()/2.0);
				Assert.assertEquals(pisaInterf.getInterfaceArea(), myInterf.getInterfaceArea()/2.0, pisaInterf.getInterfaceArea()*0.10);
				i--;
			}
			
		}
	}

	private static int getNumProtProtInterfaces (List<PisaInterface> pisaInterfs) {
		int count = 0;
		for (PisaInterface interf:pisaInterfs) {
			if (interf.isProtein()) {
				count++;
			}
		}
		return count;
	}
	
//	private static ChainInterface findInterface(List<ChainInterface> list, PisaInterface pisaInterf) {
//		String pdbChainCode1 = pisaInterf.getFirstMolecule().getChainId();
//		String pdbChainCode2 = pisaInterf.getSecondMolecule().getChainId();
//		for (ChainInterface interf:list) {
//			Matrix4d firstTransf = interf.getFirstTransf();
//			Matrix4d secondTransf = interf.getSecondTransf();
//			if (interf.getFirstMolecule().getPdbChainCode().equals(pdbChainCode1) && interf.getSecondMolecule().getPdbChainCode().equals(pdbChainCode2)){
//				
//			} else if (interf.getSecondMolecule().getPdbChainCode().equals(pdbChainCode1) && interf.getFirstMolecule().getPdbChainCode().equals(pdbChainCode2)) {
//				firstTransf = interf.getSecondTransf();
//				secondTransf = interf.getFirstTransf();
//			} else {
//				return null;
//			}
//			if (deltaComp(firstTransf,pisaInterf.getFirstMolecule().getSymOp()) && deltaComp(secondTransf,pisaInterf.getSecondMolecule().getSymOp())) {
//				return interf;
//			}
//		}
//		return null;
//
//	}
	
//	private static boolean deltaComp(Matrix4d m1, Matrix4d m2) {
//		double delta = 0.000001;
//		return (deltaComp(m1.m00,m2.m00,delta) &&
//				deltaComp(m1.m01,m2.m01,delta) &&
//				deltaComp(m1.m02,m2.m02,delta) &&
//				deltaComp(m1.m03,m2.m03,delta) &&
//				deltaComp(m1.m10,m2.m10,delta) &&
//				deltaComp(m1.m11,m2.m11,delta) &&
//				deltaComp(m1.m12,m2.m12,delta) &&
//				deltaComp(m1.m13,m2.m13,delta) &&
//				deltaComp(m1.m20,m2.m20,delta) &&
//				deltaComp(m1.m21,m2.m21,delta) &&
//				deltaComp(m1.m22,m2.m22,delta) &&
//				deltaComp(m1.m23,m2.m23,delta) &&
//				deltaComp(m1.m30,m2.m30,delta) &&
//				deltaComp(m1.m31,m2.m31,delta) &&
//				deltaComp(m1.m32,m2.m32,delta) &&
//				deltaComp(m1.m33,m2.m33,delta));
//		
//	}
//	
//	private static boolean deltaComp(double d1, double d2, double delta) {
//		return Math.abs(d1-d2)<delta;
//	}

}
