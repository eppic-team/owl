import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;

import proteinstructure.CiffilePdb;
//import proteinstructure.MsdsdPdb;
import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbasePdb;
import proteinstructure.PdbfilePdb;
import tools.MySQLConnection;


public class testGetChains {

	public static void main(String[] args) throws IOException, SQLException{
		
		String pdbFilesDir = "/project/StruPPi/BiO/DBd/PDB-REMEDIATED/data/structures/unzipped/all/pdb";
		String cifFilesDir = "/project/StruPPi/BiO/DBd/PDB-REMEDIATED/data/structures/unzipped/all/mmCIF";
		MySQLConnection conn = new MySQLConnection("white","duarte","nieve");
		
		//String[] pdbCodes = {"12as", "3eca", "1bxy"} ;
		
		String listFile = "/project/StruPPi/jose/workspace/aglappe-jung/pdbcodes.testset";
		BufferedReader fpdb = new BufferedReader(new FileReader(listFile));
		String  line;
		while ((line=fpdb.readLine())!=null) {
			String pdbCode = line.split("\t")[0].toLowerCase();
		
			System.out.println("pdb "+pdbCode);
			
			try {
				Pdb pdbasepdb = new PdbasePdb(pdbCode, "pdbase", conn);
				String[] pdbasechains = pdbasepdb.getChains();
				Integer[] pdbasemodels = pdbasepdb.getModels();
				pdbasepdb.load(pdbasechains[0]);


//				Pdb msdsdpdb = new MsdsdPdb(pdbCode,"msdsd_00_07_a",conn);
//				String[] msdsdchains = msdsdpdb.getChains();
//				Integer[] msdsdmodels = msdsdpdb.getModels();
//				msdsdpdb.load(msdsdchains[0]);

				Pdb ciffilepdb = new CiffilePdb(new File(cifFilesDir,pdbCode+".cif"));
				String[] ciffilechains = ciffilepdb.getChains();
				Integer[] ciffilemodels = ciffilepdb.getModels();
				ciffilepdb.load(ciffilechains[0]);


				Pdb pdbfilepdb = new PdbfilePdb(pdbFilesDir+"/pdb"+pdbCode+".ent");
				String[] pdbfilechains = pdbfilepdb.getChains();
				Integer[] pdbfilemodels = pdbfilepdb.getModels();
				pdbfilepdb.load(pdbfilechains[0]);
				

//			System.out.println("online cif file...");
//			pdb = new CiffilePdb(pdbCode);
//			allChains = pdb.getChains();
//			System.out.print("pdb: "+pdbCode+", chains: ");
//			for (String chain:allChains) {
//				System.out.print(chain+ " ");
//			}
//			System.out.println();
//			pdb.load("A");
//			System.out.println("length:" + pdb.getFullLength());

				boolean chainsmatch = true;
				boolean modelsmatch = true;
				if (pdbasechains.length != ciffilechains.length || pdbasechains.length!= pdbfilechains.length ) { //|| pdbasechains.length!=msdsdchains.length) {
					System.out.println("different number of chains");
					chainsmatch=false;
				}
				else {
					for (int i=0; i<pdbasechains.length;i++) {
						if (!pdbasechains[i].equals(ciffilechains[i]) || !pdbasechains[i].equals(pdbfilechains[i]) ){ // || !pdbasechains[i].equals(msdsdchains[i])){
							System.out.println("chain differs: pdbase "+pdbasechains[i]+", ciffile "+ciffilechains[i]+", pdbfile "+pdbfilechains[i]);//+", msdsd "+msdsdchains[i]);
							chainsmatch=false;
						}
					}
				}
				
				if (chainsmatch) {
					for (String chain:pdbfilechains) {
						System.out.print(chain+ " ");
					}
					System.out.println();

					if (pdbasepdb.getFullLength()!=ciffilepdb.getFullLength() || pdbasepdb.getFullLength()!=pdbfilepdb.getFullLength()){ // || pdbasepdb.getFullLength()!=msdsdpdb.getFullLength()) {
						System.out.println("length for chain "+pdbfilechains[0]+" differs");
					}
				}
				
				if (pdbasemodels.length != ciffilemodels.length || pdbasemodels.length != pdbfilemodels.length ) {//|| pdbasemodels.length != msdsdmodels.length) {
					System.out.println("different number of models");
					modelsmatch = false;
				} else {
					for (int i=0;i<pdbasemodels.length;i++){
						if (pdbasemodels[i]!=ciffilemodels[i] || pdbasemodels[i]!=pdbfilemodels[i] ){ //|| pdbasemodels[i]!=msdsdmodels[i]){
							System.out.println("model differs: pdbase "+pdbasemodels[i]+", ciffile "+ciffilemodels[i]+", pdbfile "+pdbfilemodels[i]);//+", msdsd "+msdsdmodels[i]);
							modelsmatch = false;
						}
					}
				}
				if (modelsmatch) {
					for (int model:pdbasemodels) {
						System.out.print(model+" ");
					}
					System.out.println();
				}
				
				
			} catch (PdbCodeNotFoundError e) {
				System.err.println("pdb code was not found in pdbase");
			} catch (PdbLoadError e) {
				System.err.println("couldn't load pdb data, specific error: "+e.getMessage());
			}
		}
		fpdb.close();
		
		

		


	}

}
