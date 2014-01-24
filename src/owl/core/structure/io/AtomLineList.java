package owl.core.structure.io;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;

import owl.core.structure.AminoAcid;
import owl.core.structure.Nucleotide;
import owl.core.util.FileFormatException;


public class AtomLineList implements Iterable<AtomLine> {
	
	private ArrayList<AtomLine> list;
	private TreeSet<String> allAtomAltLocs;
	
	private HashMap<String,Integer> insCodesPerChain; // pdb chain codes to number of insertion codes per chain
	
	private TreeMap<String,String> pdbchaincode2chaincode;
	private TreeMap<String,ArrayList<AtomLine>> atomLineGroups;
	
	public AtomLineList() {
		list = new ArrayList<AtomLine>();
		allAtomAltLocs = new TreeSet<String>(); 
	}
	
	public void addAtomLine(AtomLine atomLine) {
		list.add(atomLine);
		if (!atomLine.labelAltId.equals(".")) {
			allAtomAltLocs.add(atomLine.labelAltId);
		}
	}

	public boolean isEmpty() {
		return list.isEmpty();
	}
	
	public int size() {
		return list.size();
	}
	
	public AtomLine get(int i) {
		return list.get(i);
	}
	
	@Override
	public Iterator<AtomLine> iterator() {
		return list.iterator();
	}
	
	/**
	 * Returns the first (in alphabetical order) alt loc id found for this list of AtomLines (excluding the default ".")
	 * or null if no alt locs exist for this list.
	 * @return
	 */
	public String getAtomAltLoc() {
		if (allAtomAltLocs.isEmpty()) {
			return null;
		} else {
			return allAtomAltLocs.first();
		}
	}
	
	/**
	 * Returns the map with mapping of pdbChainCodes to our assigned chainCodes (CIF-like codes)
	 * Calculated after call to {@link #sortIntoChains()}
	 * To be used in PDB file parsing only.
	 * @return
	 */
	public TreeMap<String,String> getPdbChainCode2chainCode() {
		return pdbchaincode2chaincode;
	}
	
	/**
	 * Returns a map of our assigned chainCodes (CIF-like codes) to lists of AtomLines belonging to them 
	 * Calculated after call to {@link #sortIntoChains()}
	 * To be used in PDB file parsing only. 
	 * @return
	 */
	public TreeMap<String,ArrayList<AtomLine>> getAtomLineGroups() {
		return atomLineGroups;
	}
	
	/**
	 * Returns the number of insertion codes for given pdb chain code (auth asym id)
	 * To be used in PDB file parsing only. 
	 * @param pdbChainCode
	 * @return
	 */
	public int getNumInsCodeForChain(String pdbChainCode) {
		if (insCodesPerChain==null) {
			insCodesPerChain = new HashMap<String, Integer>();
			for (AtomLine atomLine:this) {
				if (!atomLine.insCode.equals(".")) {
					if (!insCodesPerChain.containsKey(atomLine.authAsymId)) {
						insCodesPerChain.put(atomLine.authAsymId,1);
					} else {
						int n = insCodesPerChain.get(atomLine.authAsymId)+1;
						insCodesPerChain.put(atomLine.authAsymId,n);
					}
				}
			}
		} 
		if (!insCodesPerChain.containsKey(pdbChainCode)) {
			return 0; 
		}
		return insCodesPerChain.get(pdbChainCode);
	}
	
	/**
	 * Sorts the atom lines into groups of chains assigning CIF codes (labelAsymIds or chainCodes) to them
	 * Also assigns the isNonPoly field of all AtomLines in this list and finds out the pdbchaincode2chaincode mapping.
	 * Retrieve the data subsequently with {@link #getAtomLineGroups()} and {@link #getPdbChainCode2chainCode()}
	 * To be used in PDB file parsing only. 
	 * @param terRecordSeen whether at least one TER record is present in PDB file or not
	 * @throws FileFormatException 
	 */
	public void sortIntoChainsPDBFormat(boolean terRecordSeen) throws FileFormatException {
		atomLineGroups = new TreeMap<String, ArrayList<AtomLine>>();
		
		String lastPdbChainCode = null;
		boolean lastOutOfPolyChain = false;
		boolean lastIsHetAtom = false;
		 
		String chainCode = "A";
		String currentChainCode = chainCode;
		
		// two strategies to follow here: 
		//  1) when at least a TER record is present we can rely on chain codes and outOfPolyChain fields exclusively
		if (terRecordSeen) {
			for (AtomLine atomLine:this) {
				if (	lastPdbChainCode==null || 
						!lastPdbChainCode.equals(atomLine.authAsymId) || 
						(lastOutOfPolyChain==false && atomLine.outOfPolyChain==true)) {

					ArrayList<AtomLine> list = new ArrayList<AtomLine>();
					list.add(atomLine);
					atomLineGroups.put(chainCode,list);
					currentChainCode = chainCode;
					chainCode = getNextChainCode(chainCode);
				} else {
					atomLineGroups.get(currentChainCode).add(atomLine);
				}

				atomLine.labelAsymId = currentChainCode;

				lastPdbChainCode = atomLine.authAsymId;
				lastOutOfPolyChain = atomLine.outOfPolyChain;
			}
		// 2) no TER records at all: we must also try to guess if a HETATOM is not peptide-linked 
		//    (otherwise lacking the TER records, we would extend the poly chains into non-poly parts, e.g. 1c52.pdb without TER records)
		//    This strategy is dangerous because our knowledge of HETATOMs that are peptide-linked is not 
		//    very comprehensive (see HetAtom.isPeptideLinked()). We will sometimes think wrongly that a certain hetatom line
		//    belongs to a not-peptide-linked residue and then cut wrongly the chain into a poly and a non-poly part
		} else {
			for (AtomLine atomLine:this) {
				if (	lastPdbChainCode==null || 
						!lastPdbChainCode.equals(atomLine.authAsymId) || 
						(lastOutOfPolyChain==false && atomLine.outOfPolyChain==true) ||
						// this is the new condition is strategy 2 vs 1)
						(lastIsHetAtom==false && atomLine.lineIsHetAtm==true && !atomLine.isPeptideLinked()) ) {

					ArrayList<AtomLine> list = new ArrayList<AtomLine>();
					list.add(atomLine);
					atomLineGroups.put(chainCode,list);
					currentChainCode = chainCode;
					chainCode = getNextChainCode(chainCode);
				} else {
					atomLineGroups.get(currentChainCode).add(atomLine);
				}

				atomLine.labelAsymId = currentChainCode;

				lastPdbChainCode = atomLine.authAsymId;
				lastOutOfPolyChain = atomLine.outOfPolyChain;
				lastIsHetAtom = atomLine.lineIsHetAtm;
			}
		}
		
		pdbchaincode2chaincode = new TreeMap<String, String>();
		
		// we now check for consistency in atom numbering
		for (ArrayList<AtomLine> group:atomLineGroups.values()) {
			int lastAtomSerial = -1;	
			for (AtomLine line:group) {
				if(line.atomserial <= lastAtomSerial) {
					throw new FileFormatException("Atom serials do not occur in ascending order in PDB file for atom " + line.atomserial + ")");
				}
				lastAtomSerial = line.atomserial;
			}
		}
		
		// now chains are all assigned to chainCodes
		// we can now fill the isNonPoly field indicating whether a chain is non-polymer or polymer (it's polymer when it has at least one ATOM line)
		// we have to do this for the cases when a polymer chain starts with a HETATM, otherwise we could rely directly on the outOfPolyChain field
		for (ArrayList<AtomLine> group:atomLineGroups.values()) {
			
			boolean hasOnePolyAtom = false;
			for (AtomLine atomLine:group) {
				if (!atomLine.lineIsHetAtm && 
					// some bad PDB files (e.g. output of old phenix versions) don't use HETATM to tag a het atom but simply ATOM, 
					// that's why we also check this (UNK we check for the rare cases when a deposited PDB file has a poly chain with only unknown residues):
					(AminoAcid.isStandardAA(atomLine.res_type) || Nucleotide.isStandardNuc(atomLine.res_type) || atomLine.res_type.equals("UNK"))) 
						hasOnePolyAtom = true;
			}
			// if hasOnePolyAtom is true that means that the chain is poly, we assign isNonPoly accordingly:
			if (hasOnePolyAtom) {
				String previous = pdbchaincode2chaincode.put(group.get(0).authAsymId,group.get(0).labelAsymId);
				
				// We check that we are not reassigning the poly chain: if we are that means there's something wrong in the file
				// For instance a badly placed TER record before a HETATM line (followed by ATOM lines) can cause this problem
				if (previous!=null) 
					throw new FileFormatException("Polymer PDB chain code "+group.get(0).authAsymId+" assigned twice, to chain codes: "+
							group.get(0).labelAsymId+" and "+previous+
							". Most likely there is something wrong with this PDB file: check that TER records are correctly placed and that all chains have a unique chain code");
				
				for (AtomLine atomLine:group) {
					atomLine.isNonPoly = false;
				}
			} else {
				for (AtomLine atomLine:group) {
					atomLine.isNonPoly = true;
				}
			}			

		}
		
	}
	
	/**
	 * Sorts the atom lines into groups of chains assigning poly/non-poly chains
	 * Also assigns the pdbchaincode2chaincode mapping.
	 * Retrieve the data subsequently with {@link #getAtomLineGroups()} and {@link #getPdbChainCode2chainCode()}
	 * To be used in mmCIF file parsing only when no pdbx_poly_seq_scheme is available  
	 * @throws FileFormatException
	 */
	public void sortIntoChainsCIFFormat() throws FileFormatException {
		atomLineGroups = new TreeMap<String, ArrayList<AtomLine>>();
		
		String lastChainCode = null;
		
		for (AtomLine atomLine:this) {
			if (	lastChainCode==null || 
					!lastChainCode.equals(atomLine.labelAsymId) 
					) {

				ArrayList<AtomLine> list = new ArrayList<AtomLine>();
				list.add(atomLine);
				atomLineGroups.put(atomLine.labelAsymId,list);
			} else {
				atomLineGroups.get(atomLine.labelAsymId).add(atomLine);
			}

			lastChainCode = atomLine.labelAsymId;
		}
		
		pdbchaincode2chaincode = new TreeMap<String, String>();
		
		// we now check for consistency in atom numbering
		for (ArrayList<AtomLine> group:atomLineGroups.values()) {
			int lastAtomSerial = -1;	
			for (AtomLine line:group) {
				if(line.atomserial <= lastAtomSerial) {
					throw new FileFormatException("Atom serials do not occur in ascending order in CIF file for atom " + line.atomserial + ")");
				}
				lastAtomSerial = line.atomserial;
			}
		}
		
		// now chains are all assigned to chainCodes
		// we can now fill the isNonPoly field indicating whether a chain is non-polymer or polymer (it's polymer when it has at least one ATOM line)
		for (ArrayList<AtomLine> group:atomLineGroups.values()) {
			
			boolean hasOnePolyAtom = false;
			for (AtomLine atomLine:group) {
				if (!atomLine.lineIsHetAtm && 
					// some bad PDB files (e.g. output of old phenix versions) don't use HETATM to tag a het atom but simply ATOM, 
					// that's why we also check this (UNK we check for the rare cases when a deposited PDB file has a poly chain with only unknown residues):
					// (not sure if this applies to CIF files too, but keeping it in case)
					(AminoAcid.isStandardAA(atomLine.res_type) || Nucleotide.isStandardNuc(atomLine.res_type) || atomLine.res_type.equals("UNK"))) 
						hasOnePolyAtom = true;
			}
			// if hasOnePolyAtom is true that means that the chain is poly, we assign isNonPoly accordingly:
			if (hasOnePolyAtom) {
				String previous = pdbchaincode2chaincode.put(group.get(0).authAsymId,group.get(0).labelAsymId);
				
				// We check that we are not reassigning the poly chain: if we are that means there's something wrong in the file
				// Not sure if this can happen inf CIF file but we'll keep the check here anyway
				if (previous!=null) 
					throw new FileFormatException("Polymer PDB chain code "+group.get(0).authAsymId+" assigned twice, to chain codes: "+
							group.get(0).labelAsymId+" and "+previous+
							". Most likely there is something wrong with this CIF file: check that all chains have a unique chain code");
				
				for (AtomLine atomLine:group) {
					atomLine.isNonPoly = false;
				}
			} else {
				for (AtomLine atomLine:group) {
					atomLine.isNonPoly = true;
				}
			}			

		}

		
	}
	
	/**
	 * Gets the next CIF chain code given a CIF chain code according to the convention followed by the PDB
	 * i.e.: A,B,...,Z,AA,BA,CA,...,ZA,AB,BB,CB,...,ZB,.......,ZZ,AAA,BAA,CAA,...
	 * @param chainCode
	 * @return
	 */
	private String getNextChainCode(String chainCode) {
		if (chainCode.length()==1) {
			if (!chainCode.equals("Z")) {
				return Character.toString(getNextChar(chainCode.charAt(0)));
			} else {
				return "AA";
			}
		} else if (chainCode.length()==2) {
			if (chainCode.equals("ZZ")) {
				return "AAA";
			}
			char[] c = new char[2];
			chainCode.getChars(0, 2, c, 0);
			c[0] = getNextChar(c[0]);
			if (c[0]=='A') {
				c[1] = getNextChar(c[1]);
			} 
			return new String(c);
		} else if (chainCode.length()==3) {
			char[] c = new char[3];
			chainCode.getChars(0, 3, c, 0);
			c[0] = getNextChar(c[0]);
			if (c[0]=='A') {
				c[1] = getNextChar(c[1]);
				if (c[1]=='A') {
					c[2] = getNextChar(c[2]);
				}
			}
			return new String(c);
		}
		return null;
	}
	
	private char getNextChar(char c) {
		if (c!='Z') {
			return ((char)(c+1));
		} else {
			return 'A';
		}
	}
}
