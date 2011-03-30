package owl.core.runners;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.TreeMap;

import owl.core.structure.PdbChain;
import owl.core.structure.features.SecStrucElement;
import owl.core.structure.features.SecondaryStructure;


/**
 * Class to run the DSSP program (assignment of secondary structure)
 * As of September 2007, a DSSP executable can be downloaded from http://swift.cmbi.ru.nl/gv/dssp/
 * after filling in a license agreement.
 * 
 * @author duarte_j
 *
 */
public class DsspRunner {
	

	/** 
	 * Runs an external DSSP executable and (re)assigns the secondary structure annotation from the parsed output.
	 * The resulting secondary structure information will have 4 states.
	 * @param dsspExecutable
	 * @param dsspParameters for current version of DSSP set this to "--" (two hyphens)
	 * @return the secondary structure annotation object
	 */
	public static SecondaryStructure runDssp(PdbChain pdb, String dsspExecutable, String dsspParameters) throws IOException {
		return runDssp(pdb, dsspExecutable, dsspParameters, SecStrucElement.ReducedState.FOURSTATE, SecStrucElement.ReducedState.FOURSTATE);
	}

	/** 
	 * Runs an external DSSP executable and (re)assigns the secondary structure annotation from the parsed output.
	 * @param dsspExecutable
	 * @param dsspParameters for current version of DSSP set this to "--" (two hyphens)
	 * @param state4Type
	 * @param state4Id
	 * @return the secondary structure annotation object 
	 * @return
	 */
	public static SecondaryStructure runDssp(PdbChain pdb, String dsspExecutable, String dsspParameters, SecStrucElement.ReducedState state4Type, SecStrucElement.ReducedState state4Id) throws IOException {
		String pdbCode = pdb.getPdbCode();
		String chainCode = pdb.getChainCode();
		
		String startLine = "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC";
		String line;
		int lineCount = 0;
		char ssType, sheetLabel;
		TreeMap<Integer, Character> ssTypes;
		TreeMap<Integer, Character> sheetLabels;
		int resNum;
		String resNumStr;
		File test = new File(dsspExecutable);
		if(!test.canRead()) throw new IOException("DSSP Executable is not readable");
		Process myDssp = Runtime.getRuntime().exec(dsspExecutable + " " + dsspParameters);
		PrintStream dsspInput = new PrintStream(myDssp.getOutputStream());
		BufferedReader dsspOutput = new BufferedReader(new InputStreamReader(myDssp.getInputStream()));
		BufferedReader dsspError = new BufferedReader(new InputStreamReader(myDssp.getErrorStream()));
		pdb.writeAtomLines(dsspInput);	// pipe atom lines to dssp
		dsspInput.close();
		ssTypes = new TreeMap<Integer,Character>();
		sheetLabels = new TreeMap<Integer,Character>();
		while((line = dsspOutput.readLine()) != null) {
			lineCount++;
			if(line.startsWith(startLine)) {
				//System.out.println("Dssp Output: ");
				break;
			}
		}
		while((line = dsspOutput.readLine()) != null) {
			lineCount++;
			resNumStr = line.substring(5,10).trim();
			ssType = line.charAt(16);			
			sheetLabel = line.charAt(33);
			if (state4Id == SecStrucElement.ReducedState.FOURSTATE && SecStrucElement.getReducedStateTypeFromDsspType(ssType, state4Id) == SecStrucElement.OTHER) {
				sheetLabel = ' ';
			}
			if(!resNumStr.equals("")) {		// this should only happen if dssp inserts a line indicating a chain break
				try {
					resNum = Integer.valueOf(resNumStr);
					ssTypes.put(resNum, ssType);
					sheetLabels.put(resNum, sheetLabel);
				} catch (NumberFormatException e) {
					System.err.println("Error while parsing DSSP output for "+pdbCode+"_"+chainCode+". Expected residue number, found '" + resNumStr + "' in line " + lineCount);
				}
			}
		}
		//for(char c:ssTypes) {System.out.print(c);}; System.out.println(".");
		dsspOutput.close();
		dsspError.close();

		if(ssTypes.size() == 0) {
			throw new IOException("No DSSP output found.");
		}

		if(ssTypes.size() != pdb.getObsLength()) {	// compare with number of observed residues
			System.err.println("Error: DSSP output size (" + ssTypes.size() + ") for "+pdbCode+"_"+chainCode+" does not match number of observed residues in structure (" + pdb.getObsLength() + ").");
		}

		// assign secondary structure
		SecondaryStructure secondaryStructure = new SecondaryStructure(pdb.getSequence().getSeq());
		char lastType = SecStrucElement.getReducedStateTypeFromDsspType(ssTypes.get(ssTypes.firstKey()), state4Id);
		int lastResSer = ssTypes.firstKey();
		char lastSheet = sheetLabels.get(lastResSer);
		char thisType, thisSheet, reducedType;
		int start = 1;
		int elementCount = 0;
		SecStrucElement ssElem;
		String ssId;
		for(int resSer:ssTypes.keySet()) {
			thisType = SecStrucElement.getReducedStateTypeFromDsspType(ssTypes.get(resSer), state4Id);
			thisSheet = sheetLabels.get(resSer);
			if(thisType != lastType || thisSheet != lastSheet || resSer > lastResSer+1) {
				// finish previous element, start new one
				elementCount++;
				reducedType = SecStrucElement.getReducedStateTypeFromDsspType(ssTypes.get(lastResSer), state4Type);
				ssId = new Character(lastType).toString() + 
						(lastSheet==' '?"":String.valueOf(lastSheet)) + new Integer(elementCount).toString();
				ssElem = new SecStrucElement(reducedType,start,lastResSer,ssId);
				secondaryStructure.add(ssElem);
				start = resSer;
				lastType = thisType;
				lastSheet = thisSheet;
			}
			lastResSer = resSer;
		}
		// finish last element
		elementCount++;
		reducedType = SecStrucElement.getReducedStateTypeFromDsspType(ssTypes.get(ssTypes.lastKey()), state4Type);
		ssId = new Character(lastType).toString() + 
				(lastSheet==' '?"":String.valueOf(lastSheet)) + new Integer(elementCount).toString();
		ssElem = new SecStrucElement(reducedType, start,ssTypes.lastKey(),ssId);
		secondaryStructure.add(ssElem);

		secondaryStructure.setComment("DSSP");
		return secondaryStructure;
	}

}
