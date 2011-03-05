package owl.tests.core.structure;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import junit.framework.Assert;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import owl.core.structure.AminoAcid;
import owl.core.structure.CiffilePdb;
import owl.core.structure.CrystalCell;
import owl.core.structure.Pdb;
import owl.core.structure.PdbCodeNotFoundException;
import owl.core.structure.PdbLoadException;
import owl.core.structure.PdbasePdb;
import owl.core.structure.PdbfilePdb;
import owl.core.structure.Residue;
import owl.core.structure.features.SecondaryStructure;
import owl.core.util.FileFormatException;
import owl.core.util.MySQLConnection;


public class PdbParsersTest {
	
	// paths to the mirror of the PDB ftp repository with all pdb/mmCIF compressed files
	private static final String CIFDIR="/nfs/data/dbs/pdb/data/structures/all/mmCIF/";
	private static final String PDBDIR="/nfs/data/dbs/pdb/data/structures/all/pdb/";
	
	private static final String LISTFILE = "src/owl/tests/core/structure/data/cullpdb_20";
	private static final String PDBASE_DB = "pdbase";
	
	private static final String TINKERPDBFILE = "src/owl/tests/core/structure/data/1bxyA_tinker.pdb";
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		// to throw away jaligner logging
		System.setProperty("java.util.logging.config.file","/dev/null");
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testSpecialPDBFiles() throws PdbLoadException {
		// testing some special PDB files
		
		// tinker file (no occupancy or bfactors)
		Pdb pdb = new PdbfilePdb(TINKERPDBFILE);
		pdb.load(Pdb.NULL_CHAIN_CODE);		
	}
	
	@Test
	public void testCIFagainstPDBASE() throws FileFormatException, SQLException, IOException {

		MySQLConnection conn = new MySQLConnection();
		
		BufferedReader flist = new BufferedReader(new FileReader(LISTFILE));
		String line;
		while ((line = flist.readLine() ) != null ) {
			if (line.startsWith("#")) continue;
			String pdbCode = line.split("\\s+")[0].toLowerCase();
			String pdbChainCode = line.split("\\s+")[1];
			
			String message = "Failed for "+pdbCode+pdbChainCode;
			
			System.out.println(pdbCode+" "+pdbChainCode);
			
			Pdb ciffilePdb = null;
			Pdb pdbasePdb = null;

			try {
				File cifFile = unzipFile(new File(CIFDIR,pdbCode+".cif.gz"));
				ciffilePdb = new CiffilePdb(cifFile);
				ciffilePdb.load(pdbChainCode);
				
				pdbasePdb = new PdbasePdb(pdbCode, PDBASE_DB, conn);
				pdbasePdb.load(pdbChainCode);
				
				// asserting
				
				// getChains/getModels
				String[] ciffileChains = ciffilePdb.getChains(); 
				String[] pdbaseChains = pdbasePdb.getChains();
				Integer[] ciffileModels = ciffilePdb.getModels(); 
				Integer[] pdbaseModels = pdbasePdb.getModels();
				Assert.assertTrue(ciffileChains.length==pdbaseChains.length);
				Assert.assertTrue(ciffileModels.length==pdbaseModels.length);
				HashSet<String> ciffileChainsSet = new HashSet<String>(Arrays.asList(ciffileChains));
				HashSet<String> pdbaseChainsSet = new HashSet<String>(Arrays.asList(pdbaseChains));
				for (String chain:pdbaseChains) {
					Assert.assertTrue(ciffileChainsSet.contains(chain)); 
				}
				for (String chain:ciffileChains) {
					Assert.assertTrue(pdbaseChainsSet.contains(chain)); 
				}
				HashSet<Integer> ciffileModelsSet = new HashSet<Integer>(Arrays.asList(ciffileModels));
				HashSet<Integer> pdbaseModelsSet = new HashSet<Integer>(Arrays.asList(pdbaseModels));
				for (int model:pdbaseModels) {
					Assert.assertTrue(ciffileModelsSet.contains(model));
				}
				for (int model:ciffileModels) {
					Assert.assertTrue(pdbaseModelsSet.contains(model));
				}
				// crystal data
				Assert.assertEquals(pdbasePdb.getSpaceGroup(), ciffilePdb.getSpaceGroup());
				CrystalCell pdbaseCell = pdbasePdb.getCrystalCell();
				CrystalCell ciffileCell = ciffilePdb.getCrystalCell();
				if (pdbaseCell!=null && ciffileCell!=null) { // nulls will happen for NMR entries
					Assert.assertEquals(pdbaseCell.getA(), ciffileCell.getA(), 0.001);
					Assert.assertEquals(pdbaseCell.getB(), ciffileCell.getB(), 0.001);
					Assert.assertEquals(pdbaseCell.getC(), ciffileCell.getC(), 0.001);
					Assert.assertEquals(pdbaseCell.getAlpha(), ciffileCell.getAlpha(), 0.001);
					Assert.assertEquals(pdbaseCell.getBeta(), ciffileCell.getBeta(), 0.001);
					Assert.assertEquals(pdbaseCell.getGamma(), ciffileCell.getGamma(), 0.001);
				}
				
				// exp data and quality parameters
				Assert.assertEquals(pdbasePdb.getExpMethod(),ciffilePdb.getExpMethod());
				Assert.assertEquals(pdbasePdb.getResolution(),ciffilePdb.getResolution(),0.0001);
				Assert.assertEquals(pdbasePdb.getRfree(),ciffilePdb.getRfree(),0.0001);
				Assert.assertEquals(pdbasePdb.getRsym(),ciffilePdb.getRsym(),0.0001);
				
				// identifiers
				Assert.assertEquals(pdbasePdb.getPdbCode(), ciffilePdb.getPdbCode());
				Assert.assertEquals(pdbasePdb.getChainCode(), ciffilePdb.getChainCode());
				Assert.assertEquals(pdbasePdb.getPdbChainCode(), ciffilePdb.getPdbChainCode());
				Assert.assertEquals(pdbasePdb.getModel(), ciffilePdb.getModel());

				// sequences
				Assert.assertEquals(pdbasePdb.getSequence(), ciffilePdb.getSequence());
				Assert.assertEquals(pdbasePdb.getObsSequence(), ciffilePdb.getObsSequence());
				
				// lengths
				Assert.assertEquals(pdbasePdb.getFullLength(), ciffilePdb.getFullLength());
				Assert.assertEquals(pdbasePdb.getObsLength(), ciffilePdb.getObsLength());
				Assert.assertEquals(pdbasePdb.getNumAtoms(), ciffilePdb.getNumAtoms());
				
				// info from atom serials 
				for (int atomser:pdbasePdb.getAllAtomSerials()) {
					Assert.assertEquals(pdbasePdb.getAtomCoord(atomser), ciffilePdb.getAtomCoord(atomser));
					Assert.assertEquals(pdbasePdb.getAtom(atomser).getOccupancy(), ciffilePdb.getAtom(atomser).getOccupancy(),0.001);
					Assert.assertEquals(pdbasePdb.getAtom(atomser).getBfactor(), ciffilePdb.getAtom(atomser).getBfactor(),0.001);
					Assert.assertEquals(pdbasePdb.getResSerFromAtomSer(atomser), ciffilePdb.getResSerFromAtomSer(atomser));
					Assert.assertEquals(pdbasePdb.getAtomNameFromAtomSer(atomser), ciffilePdb.getAtomNameFromAtomSer(atomser));
				}
				
				for (int atomser:ciffilePdb.getAllAtomSerials()) {
					Assert.assertEquals(ciffilePdb.getAtomCoord(atomser), pdbasePdb.getAtomCoord(atomser));
					Assert.assertEquals(ciffilePdb.getResSerFromAtomSer(atomser), pdbasePdb.getResSerFromAtomSer(atomser));
				}

				// info from residue serials
				for (int resser:pdbasePdb.getAllSortedResSerials()) {
					Assert.assertEquals(pdbasePdb.getPdbResSerFromResSer(resser), ciffilePdb.getPdbResSerFromResSer(resser));
					String pdbresser = pdbasePdb.getPdbResSerFromResSer(resser);
					Assert.assertEquals(pdbasePdb.getResSerFromPdbResSer(pdbresser), ciffilePdb.getResSerFromPdbResSer(pdbresser));
								
					Assert.assertEquals(pdbasePdb.getResidue(resser).getAaType(),ciffilePdb.getResidue(resser).getAaType());
				}

				for (int resser:ciffilePdb.getAllSortedResSerials()) {
					Assert.assertEquals(ciffilePdb.getPdbResSerFromResSer(resser), pdbasePdb.getPdbResSerFromResSer(resser));
					String pdbresser = ciffilePdb.getPdbResSerFromResSer(resser);
					Assert.assertEquals(ciffilePdb.getResSerFromPdbResSer(pdbresser), pdbasePdb.getResSerFromPdbResSer(pdbresser));
					
				}
				
				// secondary structure
				SecondaryStructure pdbaseSS = pdbasePdb.getSecondaryStructure();
				SecondaryStructure ciffileSS = ciffilePdb.getSecondaryStructure();
				Assert.assertEquals(message,pdbaseSS.getNumElements(), ciffileSS.getNumElements());

				for (int resser:pdbasePdb.getAllSortedResSerials()) {
					String resserMsg = "Failed for "+pdbCode+pdbChainCode+" and resser "+resser; 
					Assert.assertEquals(resserMsg,pdbaseSS.getSecStrucElement(resser), ciffileSS.getSecStrucElement(resser));
					// checking that the 2 ways of accessing the sec struct element coincide
					Assert.assertSame(pdbaseSS.getSecStrucElement(resser), pdbasePdb.getResidue(resser).getSsElem());
					Assert.assertEquals(resserMsg,pdbaseSS.getSecStrucElement(resser),pdbasePdb.getResidue(resser).getSsElem());
					Assert.assertSame(ciffileSS.getSecStrucElement(resser), ciffilePdb.getResidue(resser).getSsElem());
					Assert.assertEquals(resserMsg,ciffileSS.getSecStrucElement(resser),ciffilePdb.getResidue(resser).getSsElem());
				}
				
				
			} catch (FileNotFoundException e){
				System.err.println("File missing. "+e.getMessage());
			} catch (PdbCodeNotFoundException e) {
				System.err.println("pdb code not found in pdbase");
			} catch (PdbLoadException e) {
				System.err.println("pdb load error, cause: "+e.getMessage());
			}

			
		}
		
		flist.close();

	}
	
	@Test
	public void testPdbfileParser() throws IOException, SQLException {
		
		ArrayList<String> warnings = new ArrayList<String>();
		// testing a list of PDB files from PDB
		BufferedReader flist = new BufferedReader(new FileReader(LISTFILE));
		String line;
		while ((line = flist.readLine() ) != null ) {
			if (line.startsWith("#")) continue;
			String pdbCode = line.split("\\s+")[0].toLowerCase();
			//String pdbChainCode = line.split("\\s+")[1];
			
			System.out.println(pdbCode);
			
			Pdb pdbfilePdb = null;
			Pdb pdbasePdb = null;
			
			try {
				File pdbFile = unzipFile(new File(PDBDIR,"pdb"+pdbCode+".ent.gz"));
				pdbfilePdb = new PdbfilePdb(pdbFile.getAbsolutePath());
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
						
						// crystal data
						if (pdbasePdb.getSpaceGroup()!=null) {
							Assert.assertEquals(pdbasePdb.getSpaceGroup().getId(),pdbfilePdb.getSpaceGroup().getId());
						}
						Assert.assertNotNull(pdbfilePdb.getCrystalCell());
						
						// exp data and quality parameters
						Assert.assertEquals(pdbasePdb.getExpMethod(),pdbfilePdb.getExpMethod());
						Assert.assertEquals(pdbasePdb.getResolution(),pdbfilePdb.getResolution(),0.01);
						Assert.assertEquals(pdbasePdb.getRfree(),pdbfilePdb.getRfree(),0.01);
						Assert.assertEquals(pdbasePdb.getRsym(),pdbfilePdb.getRsym(),0.001);

						// sequence
						Assert.assertTrue(pdbfilePdb.getObsSequence().length()<=pdbfilePdb.getSequence().length());
						Assert.assertTrue(pdbfilePdb.getObsLength()==pdbfilePdb.getObsSequence().length());
						Assert.assertTrue(pdbfilePdb.getFullLength()==pdbfilePdb.getSequence().length());
						// we can't assert the sequences against pdbase because of some very weird entries. See http://pdbwiki.org/index.php/1ejg
						//Assert.assertEquals(pdbasePdb.getSequence(), pdbfilePdb.getSequence());
						Assert.assertEquals(pdbasePdb.getObsSequence(), pdbfilePdb.getObsSequence());
						String seq = pdbfilePdb.getSequence();
						for (int resser:pdbfilePdb.getAllSortedResSerials()) {
							Assert.assertEquals(pdbfilePdb.getResidue(resser).getAaType().getThreeLetterCode(), 
												AminoAcid.one2three(seq.charAt(resser-1)));
							// at least 1 atom per observed residue
							Assert.assertTrue(pdbfilePdb.getResidue(resser).getNumAtoms()>0);
							// at least the CA atom must be present, doesn't work for all, e.g. 2cioB residue 80
							//Assert.assertNotNull(pdbfilePdb.getResidue(resser).getAtom("CA"));
							
							if (pdbasePdb.containsResidue(resser)) {
								Residue pdbaseRes = pdbasePdb.getResidue(resser);
								Residue pdbfileRes = pdbfilePdb.getResidue(resser);
								Assert.assertEquals(pdbaseRes.getAaType(), pdbfileRes.getAaType());
								Assert.assertEquals(pdbaseRes.getNumAtoms(), pdbfileRes.getNumAtoms());
								Assert.assertEquals(pdbaseRes.getNumHeavyAtoms(), pdbfileRes.getNumHeavyAtoms());
							} else { 							
								// In some cases the alignment is ambiguous, see 2nwr at residues 191 onwards:
								// SEQRES           151 LTERGTTFGYNNLVVDFRSLPIMKQWAKVIYDATHSVQLPGGLGDKSGGM    200
								//                      |||||||||||||||||||||||||||||||||||||||||      |||
								// ATOM             150 LTERGTTFGYNNLVVDFRSLPIMKQWAKVIYDATHSVQLPG------GGM    193
								// the cif file places the G on the left of the gap rather than on the right
								// as we do (and they are right looking at the 3D coords)
								// We could guess it from distances in 3D or so (at least for the 2nwr case) 
								// but in general it's quite difficult to solve
								// Anyway this is a special case (and doesn't have many important consequences). It 
								// happens often enough (1 in a 100 or so)
								// we still want to test for all other cases to make sure we are aligning well, 
								// that's why here we only print a warning
								System.err.println("Residue "+resser+" wrongly mapped in pdbfile vs pdbase (not observed in pdbase and observed in pdbfile)");
								warnings.add(pdbCode+chain+": wronly mapped residue "+resser+", not observed in pdbase but observed in pdbfile");								
							}

						}

						// info from atom serials
						for (int atomser:pdbfilePdb.getAllAtomSerials()) {
							Assert.assertNotNull(pdbfilePdb.getAtomCoord(atomser));
							Assert.assertNotNull(pdbfilePdb.getResSerFromAtomSer(atomser));
						}

						// sec structure
						SecondaryStructure ss = pdbfilePdb.getSecondaryStructure();
						Assert.assertNotNull(ss);
						Assert.assertEquals(pdbfilePdb.getSequence(), ss.getSequence());

						for (int resser:pdbfilePdb.getAllSortedResSerials()) {
							String resserMsg = "Failed for "+pdbCode+chain+" and resser "+resser;
							// checking that the 2 ways of accessing the sec struct element coincide
							Assert.assertSame(ss.getSecStrucElement(resser), pdbfilePdb.getResidue(resser).getSsElem());
							Assert.assertEquals(resserMsg,ss.getSecStrucElement(resser),pdbfilePdb.getResidue(resser).getSsElem());
						}

					}
				}		
				
			} catch (FileNotFoundException e){
				System.err.println("File missing. "+e.getMessage());
			}catch (PdbLoadException e) {
				System.err.println("pdb load error, cause: "+e.getMessage());
			} catch (PdbCodeNotFoundException e) {
				System.err.println("pdb code not found in pdbase");			
			}
		}
		flist.close();
		System.err.println("Warnings for: ");
		for (String warning:warnings) {
			System.err.println(warning);
		}
	}
	
	private static File unzipFile(File repoGzFile) throws FileNotFoundException {
		if (!repoGzFile.exists()) {
			throw new FileNotFoundException("PDB repository file "+repoGzFile+" could not be found.");
		}
		File unzippedFile = null;
		try {
			String prefix = repoGzFile.getName().substring(0,repoGzFile.getName().lastIndexOf(".gz"));
			unzippedFile = File.createTempFile(prefix,"");
			unzippedFile.deleteOnExit();

			GZIPInputStream zis = new GZIPInputStream(new FileInputStream(repoGzFile));
			FileOutputStream os = new FileOutputStream(unzippedFile);
			int b;
			while ( (b=zis.read())!=-1) {
				os.write(b);
			}
			zis.close();
			os.close();
		} catch (IOException e) {
			System.err.println("Couldn't uncompress "+repoGzFile+" file into "+unzippedFile);
			System.err.println(e.getMessage());
			System.exit(1);
		}
		return unzippedFile;
	}

	// to debug the testing code (run as java program so that we can use normal debugger)
	public static void main(String[] args) throws Exception {
		PdbParsersTest pdbTest = new PdbParsersTest();
		pdbTest.setUp();
		//pdbTest.testSpecialPDBFiles();
		//pdbTest.testCIFagainstPDBASE();
		pdbTest.testPdbfileParser();
	}

}
