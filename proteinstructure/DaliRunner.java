package proteinstructure;

import java.io.File;
import java.io.IOException;

public class DaliRunner {	
	
	public static String doDALI(Pdb first, Pdb second,String dali_executable) throws IOException, InterruptedException {
		
		File workdir = createTempDirectory();
		first.dump2pdbfile(workdir.getAbsolutePath()+"/mod1.pdb", true);
		second.dump2pdbfile(workdir.getAbsolutePath()+"/mod2.pdb", true);
		
		String daliOutputFilename = workdir.getAbsolutePath()+"/"+first.pdbCode+first.pdbChainCode+
									"-"+second.pdbCode+second.pdbChainCode+".html";
		
		Process dali = Runtime.getRuntime().exec(new String[]{dali_executable,
				"-pairwise",
				workdir.getAbsolutePath()+"/mod1.pdb",
				workdir.getAbsolutePath()+"/mod2.pdb"}, 
				new String[]{"",""}, workdir);
		dali.waitFor();		
		new File(workdir.getAbsolutePath()+"/aln.html").renameTo(new File(daliOutputFilename));
		return daliOutputFilename;
	}
	
	
	
	private static File createTempDirectory()
    throws IOException
{
    final File temp;

    temp = File.createTempFile("temp", Long.toString(System.nanoTime()));

    if(!(temp.delete()))
    {
        throw new IOException("Could not delete temp file: " + temp.getAbsolutePath());
    }

    if(!(temp.mkdir()))
    {
        throw new IOException("Could not create temp directory: " + temp.getAbsolutePath());
    }

    return (temp);
}

	/**
	 * Test method
	 * @throws IOException 
	 * @throws PdbLoadError 
	 * @throws InterruptedException 
	 * @throws FileFormatError 
	 */
	
	public static void main(String[] args) throws IOException, PdbLoadError, InterruptedException {
		
	
		Pdb pdb1 = new PdbfilePdb("/project/StruPPi/matthias/test/7ADH.pdb");
		Pdb pdb2 = new PdbfilePdb("/project/StruPPi/matthias/test/2OUI.pdb");
		pdb1.load(pdb1.getChains()[0]);
		pdb2.load(pdb2.getChains()[2]);
		String fileName = DaliRunner.doDALI(pdb1,pdb2,"/project/StruPPi/matthias/DaliLite");
		

	}

}
