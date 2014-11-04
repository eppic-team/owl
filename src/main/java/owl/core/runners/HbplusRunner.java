package owl.core.runners;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;

import owl.core.structure.Atom;
import owl.core.structure.ChainInterface;
import owl.core.structure.PdbChain;
import edu.uci.ics.jung.graph.util.Pair;

public class HbplusRunner {

	private File pathToHBPlus;
	private File pdbFile;
	private HashSet<Pair<Atom>> hBondPairs;

	public HbplusRunner(File pathToHBPlus, File pdbFile) {
		this.pathToHBPlus = pathToHBPlus;
		this.pdbFile = pdbFile;
		this.hBondPairs = new HashSet<Pair<Atom>>();
	}

	public File getPdbFile() {
		return pdbFile;
	}

	public String runHBPlus() throws IOException, InterruptedException {
		byte[] buffer = new byte[1024];
		GZIPInputStream gZIPInputStream = new GZIPInputStream(new FileInputStream(pdbFile.getPath()));
		FileOutputStream fileOutputStream = new FileOutputStream(pdbFile.getPath().substring(0, pdbFile.getPath().lastIndexOf(".")));
		int bytes_read;
		while ((bytes_read = gZIPInputStream.read(buffer)) > 0) {
			fileOutputStream.write(buffer, 0, bytes_read);
		}
		gZIPInputStream.close();
		fileOutputStream.close();
		String command = pathToHBPlus.toString() + " " + pdbFile.getPath().substring(0, pdbFile.getPath().lastIndexOf("."));
		StringBuffer output = new StringBuffer();
		Process p;
		p = Runtime.getRuntime().exec(command);
		p.waitFor();
		BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
		String line = "";
		while ((line = reader.readLine()) != null) {
			output.append(line + "\n");
		}
		return output.toString();
	}

	public void parseOutput(ChainInterface chainInterface) throws IOException {
		BufferedReader reader = null;
		String hb2fileName = pdbFile.getName().substring(0, pdbFile.getName().lastIndexOf("."));
		hb2fileName = hb2fileName.substring(0, hb2fileName.lastIndexOf(".")) + ".hb2";
		// note HBPlus outputs to current directory (whatever that is), we move the file to our dir and use it from there
		File orighb2File = new File(".",hb2fileName);
		File hb2File = new File(pdbFile.getParent(), hb2fileName);
		hb2File.deleteOnExit();
		Files.move(orighb2File.toPath(), hb2File.toPath(), StandardCopyOption.REPLACE_EXISTING);
		reader = new BufferedReader(new FileReader(hb2File));
		String line;
		while ((line = reader.readLine()) != null) {
			if (Character.isDigit(line.charAt(1)) && Character.isDigit(line.charAt(2))) {
				String[] columns = line.split("\\s+");
				String donorChain = columns[0].substring(0, 1);
				String donorResID = String.valueOf(Integer.parseInt(columns[0].substring(1, 5)));
				String donorCode = columns[1];
				String acceptorChain = columns[2].substring(0, 1);
				String acceptorResID = String.valueOf(Integer.parseInt(columns[2].substring(1, 5)));
				String acceptorCode = columns[3];

				PdbChain firstMolecule = chainInterface.getFirstMolecule();
				PdbChain secondMolecule = chainInterface.getSecondMolecule();

				if ((firstMolecule.getChainCode().equals(donorChain)) && (chainInterface.getSecondPdbChainCodeForOutput().equals(acceptorChain))) {
					Atom donor = firstMolecule.getResidue(firstMolecule.getResSerFromPdbResSer(donorResID)).getAtom(donorCode);
					Atom acceptor = secondMolecule.getResidue(secondMolecule.getResSerFromPdbResSer(acceptorResID)).getAtom(acceptorCode);
					Pair<Atom> pair = new Pair<Atom>(donor, acceptor);
					hBondPairs.add(pair);
				}
				else if ((chainInterface.getSecondPdbChainCodeForOutput().equals(donorChain)) && (firstMolecule.getChainCode().equals(acceptorChain))) {
					Atom donor = secondMolecule.getResidue(secondMolecule.getResSerFromPdbResSer(donorResID)).getAtom(donorCode);
					Atom acceptor = firstMolecule.getResidue(firstMolecule.getResSerFromPdbResSer(acceptorResID)).getAtom(acceptorCode);
					Pair<Atom> pair = new Pair<Atom>(donor, acceptor);
					hBondPairs.add(pair);
				}
			}
		}
		reader.close();
	}

	public HashSet<Pair<Atom>> getPairs() {
		return hBondPairs;
	}
}