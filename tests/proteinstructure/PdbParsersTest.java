package tests.proteinstructure;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;

import junit.framework.Assert;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import proteinstructure.CiffileFormatError;
import proteinstructure.CiffilePdb;
import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbaseInconsistencyError;
import proteinstructure.PdbasePdb;
import proteinstructure.SecondaryStructure;
import tools.MySQLConnection;

public class PdbParsersTest {
	
	
	String cifdir="/project/StruPPi/BiO/DBd/PDB-REMEDIATED/data/structures/unzipped/all/mmCIF/";
	String listFile = "/project/StruPPi/michael/cullpdb/cullpdb_20";
	String pdbaseDB = "pdbase";
	String mysqlServer = "talyn";
	
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
	public void testCIFagainstPDBASE() throws CiffileFormatError, PdbaseInconsistencyError, SQLException, IOException {

		MySQLConnection conn = new MySQLConnection(mysqlServer, pdbaseDB);
		
		BufferedReader flist = new BufferedReader(new FileReader(listFile));
		String line;
		while ((line = flist.readLine() ) != null ) {
			String pdbCode = line.split("\\s+")[0].toLowerCase();
			String pdbChainCode = line.split("\\s+")[1];
			
			System.out.println(pdbCode+" "+pdbChainCode);
			
			Pdb ciffilePdb = null;
			Pdb pdbasePdb = null;

			try {
				
				ciffilePdb = new CiffilePdb(new File(cifdir,pdbCode+".cif"));
				ciffilePdb.load(pdbChainCode);
				
				pdbasePdb = new PdbasePdb(pdbCode, pdbaseDB, conn);
				pdbasePdb.load(pdbChainCode);
				
				// asserting
				
				// identifiers
				Assert.assertEquals(pdbasePdb.getPdbCode(), ciffilePdb.getPdbCode());
				Assert.assertEquals(pdbasePdb.getChainCode(), ciffilePdb.getChainCode());
				Assert.assertEquals(pdbasePdb.getPdbChainCode(), ciffilePdb.getPdbChainCode());
				Assert.assertEquals(pdbasePdb.getModel(), ciffilePdb.getModel());

				// sequences
				Assert.assertEquals(pdbasePdb.getSequence(), ciffilePdb.getSequence());
				Assert.assertEquals(pdbasePdb.getObservedSequence(), ciffilePdb.getObservedSequence());
				
				// lengths
				Assert.assertEquals(pdbasePdb.getFullLength(), ciffilePdb.getFullLength());
				Assert.assertEquals(pdbasePdb.get_length(), ciffilePdb.get_length());
				Assert.assertEquals(pdbasePdb.getNumAtoms(), ciffilePdb.getNumAtoms());
				
				// info from atom serials 
				for (int atomser:pdbasePdb.getAllAtomSerials()) {
					Assert.assertEquals(pdbasePdb.getAtomCoord(atomser), ciffilePdb.getAtomCoord(atomser));
					Assert.assertEquals(pdbasePdb.get_resser_from_atomser(atomser), ciffilePdb.get_resser_from_atomser(atomser));
					Assert.assertEquals(pdbasePdb.getAtomNameFromAtomSer(atomser), ciffilePdb.getAtomNameFromAtomSer(atomser));
				}
				
				for (int atomser:ciffilePdb.getAllAtomSerials()) {
					Assert.assertEquals(ciffilePdb.getAtomCoord(atomser), pdbasePdb.getAtomCoord(atomser));
					Assert.assertEquals(ciffilePdb.get_resser_from_atomser(atomser), pdbasePdb.get_resser_from_atomser(atomser));
				}

				// info from residue serials
				for (int resser:pdbasePdb.getAllSortedResSerials()) {
					Assert.assertEquals(pdbasePdb.get_pdbresser_from_resser(resser), ciffilePdb.get_pdbresser_from_resser(resser));
					String pdbresser = pdbasePdb.get_pdbresser_from_resser(resser);
					Assert.assertEquals(pdbasePdb.get_resser_from_pdbresser(pdbresser), ciffilePdb.get_resser_from_pdbresser(pdbresser));
					
					Assert.assertEquals(pdbasePdb.getResTypeFromResSerial(resser), ciffilePdb.getResTypeFromResSerial(resser));
				}

				for (int resser:ciffilePdb.getAllSortedResSerials()) {
					Assert.assertEquals(ciffilePdb.get_pdbresser_from_resser(resser), pdbasePdb.get_pdbresser_from_resser(resser));
					String pdbresser = ciffilePdb.get_pdbresser_from_resser(resser);
					Assert.assertEquals(ciffilePdb.get_resser_from_pdbresser(pdbresser), pdbasePdb.get_resser_from_pdbresser(pdbresser));
					
					Assert.assertEquals(ciffilePdb.getResTypeFromResSerial(resser), pdbasePdb.getResTypeFromResSerial(resser));
				}
				
				// secondary structure
				SecondaryStructure pdbaseSS = pdbasePdb.getSecondaryStructure();
				SecondaryStructure ciffileSS = ciffilePdb.getSecondaryStructure();
				Assert.assertEquals(pdbaseSS.getNumElements(), ciffileSS.getNumElements());
				
				for (int resser:pdbasePdb.getAllSortedResSerials()) {
					Assert.assertEquals("Failed for resser "+resser,pdbaseSS.getSecStrucElement(resser), ciffileSS.getSecStrucElement(resser));
				}
				
				
			}
			catch (PdbCodeNotFoundError e) {
				System.err.println("pdb code not found in pdbase");
			}
			catch (PdbLoadError e) {
				System.err.println("pdb load error, cause: "+e.getMessage());
			}

			
		}
		
		flist.close();

	}
	
	@Test
	public void testPdbfileParser() {
		//TODO implement it!, interesting things to test: 
		// lengths, observed, unobserved... basically every special case that we try to catch in the parser we should test here
	}
}
