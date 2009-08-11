package tests.proteinstructure;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import junit.framework.Assert;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import proteinstructure.AAinfo;
import proteinstructure.CiffileFormatError;
import proteinstructure.CiffilePdb;
import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbaseInconsistencyError;
import proteinstructure.PdbasePdb;
import proteinstructure.PdbfilePdb;
import proteinstructure.SecondaryStructure;
import tools.MySQLConnection;

public class PdbParsersTest {
	
	
	String cifdir="/project/StruPPi/BiO/DBd/PDB-REMEDIATED/data/structures/unzipped/all/mmCIF/";
	String pdbdir="/project/StruPPi/BiO/DBd/PDB-REMEDIATED/data/structures/unzipped/all/pdb/";
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

//	@Test
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
	public void testPdbfileParser() throws IOException, SQLException {
		BufferedReader flist = new BufferedReader(new FileReader(listFile));
		String line;
		while ((line = flist.readLine() ) != null ) {
			String pdbCode = line.split("\\s+")[0].toLowerCase();
			//String pdbChainCode = line.split("\\s+")[1];
			
			System.out.println(pdbCode);
			
			Pdb pdbfilePdb = null;
			Pdb pdbasePdb = null;
			
			try {
				pdbfilePdb = new PdbfilePdb(new File(pdbdir,"pdb"+pdbCode+".ent").getAbsolutePath());
				pdbasePdb = new PdbasePdb(pdbCode);
				String[] chains = pdbfilePdb.getChains();
				// test getChains/getModels
				Assert.assertTrue(chains.length>0);
				Integer[] models = pdbfilePdb.getModels();
				Assert.assertTrue(models.length>0);
				Assert.assertFalse(pdbfilePdb.isDataLoaded());
				
				// loading all chains and all models
				for (int model:models) {
					System.out.println(model);
					for (String chain:chains) {
						pdbasePdb.load(chain);
						System.out.println(chain);
						pdbfilePdb.load(chain,model);
						Assert.assertTrue(pdbfilePdb.isDataLoaded());
						
						// pdbCode properly read
						Pattern p = Pattern.compile("^\\d\\w\\w\\w$");
						Matcher m = p.matcher(pdbfilePdb.getPdbCode());
						Assert.assertTrue(m.matches());
						Assert.assertEquals(pdbasePdb.getPdbCode(), pdbfilePdb.getPdbCode());
						
						// chain codes coincide
						Assert.assertEquals(chain,pdbfilePdb.getChainCode());
						Assert.assertEquals(pdbfilePdb.getChainCode(), pdbfilePdb.getPdbChainCode());
						Assert.assertEquals(pdbasePdb.getPdbChainCode(), pdbfilePdb.getPdbChainCode());
						// model as input
						Assert.assertEquals(model,pdbfilePdb.getModel());
						
						// sequence
						Assert.assertTrue(pdbfilePdb.getObservedSequence().length()<=pdbfilePdb.getSequence().length());
						Assert.assertTrue(pdbfilePdb.get_length()==pdbfilePdb.getObservedSequence().length());
						Assert.assertTrue(pdbfilePdb.getFullLength()==pdbfilePdb.getSequence().length());
						Assert.assertEquals(pdbasePdb.getSequence(), pdbfilePdb.getSequence());
						String seq = pdbfilePdb.getSequence();
						for (int resser:pdbfilePdb.getAllSortedResSerials()) {
							Assert.assertEquals(pdbfilePdb.getResTypeFromResSerial(resser), 
												AAinfo.oneletter2threeletter(String.valueOf(seq.charAt(resser-1))));
						}

						// atom number: at least 1 atom per observed residue
						Assert.assertTrue(pdbfilePdb.getNumAtoms()>=pdbfilePdb.get_length());
						
						// info from atom serials
						for (int atomser:pdbfilePdb.getAllAtomSerials()) {
							Assert.assertNotNull(pdbfilePdb.getAtomCoord(atomser));
							Assert.assertNotNull(pdbfilePdb.get_resser_from_atomser(atomser));
						}

						// sec structure
						SecondaryStructure ss = pdbfilePdb.getSecondaryStructure();
						Assert.assertNotNull(ss);
						Assert.assertEquals(pdbfilePdb.getSequence(), ss.getSequence());

					}
				}		
				
			} catch (PdbLoadError e) {
				System.err.println("pdb load error, cause: "+e.getMessage());
			} catch (PdbCodeNotFoundError e) {
				System.err.println("pdb code not found in pdbase");			
			}
		}
		flist.close();
	}
}
