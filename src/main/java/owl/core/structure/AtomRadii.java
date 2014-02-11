package owl.core.structure;

/**
 * Class containing methods to return the van der Waals atom radii 
 * for standard amino acids and nucleotides as defined by
 * Chothia (1976) J.Mol.Biol.105,1-14
 * 
 * @author duarte_j
 *
 */
public class AtomRadii {

	// Chothia's amino acid atoms vdw radii
	public static final double TRIGONAL_CARBON_VDW = 1.76;
	public static final double TETRAHEDRAL_CARBON_VDW = 1.87;
	public static final double TRIGONAL_NITROGEN_VDW = 1.65;
	public static final double TETRAHEDRAL_NITROGEN_VDW = 1.50;
	public static final double SULFUR_VDW = 1.85;
	public static final double OXIGEN_VDW = 1.40;
	
	// Chothia's nucleotide atoms vdw radii
	public static final double NUC_CARBON_VDW = 1.80;
	public static final double NUC_NITROGEN_VDW = 1.60;
	public static final double PHOSPHOROUS_VDW = 1.90;
	

	/**
	 * Gets the radius for given amino acid and atom 
	 * @param aa
	 * @param atom
	 * @return
	 */	
	private static double getRadius(AminoAcid aa, Atom atom) {
		
		if (atom.getType().equals(AtomType.H)) return AtomType.H.getRadius();
		// some unusual entries (e.g. 1tes) contain Deuterium atoms in standard aminoacids
		if (atom.getType().equals(AtomType.D)) return AtomType.D.getRadius();

		String atomCode = atom.getCode();
		
		// here we use the values that Chothia gives in his paper (as NACCESS does)
		if (atom.getType()==AtomType.O) {
			return OXIGEN_VDW;
		} 
		else if (atom.getType()==AtomType.S) {
			return SULFUR_VDW;
		}
		else if (atom.getType()==AtomType.N) {
			if (atomCode.equals("NZ")) return TETRAHEDRAL_NITROGEN_VDW; // tetrahedral Nitrogen
			return TRIGONAL_NITROGEN_VDW;								// trigonal Nitrogen
		}
		else if (atom.getType()==AtomType.C) { // it must be a carbon
			if (atomCode.equals("C") || 
					atomCode.equals("CE1") || atomCode.equals("CE2") || atomCode.equals("CE3") ||
					atomCode.equals("CH2") || 
					atomCode.equals("CZ") || atomCode.equals("CZ2") || atomCode.equals("CZ3")) {
				return TRIGONAL_CARBON_VDW; 							// trigonal Carbon
			}
			else if (atomCode.equals("CA") || atomCode.equals("CB") || 
					atomCode.equals("CE") ||
					atomCode.equals("CG1") || atomCode.equals("CG2")) {
				return TETRAHEDRAL_CARBON_VDW;							// tetrahedral Carbon
			}
			// left cases depend on amino acid: CD, CD1, CD2, CG
			else {
				switch (aa) {
					case PHE:						
					case TRP:
					case TYR:
					case HIS:
					case ASP:
					case ASN:
						return TRIGONAL_CARBON_VDW;
						
					case PRO:
					case LYS:
					case ARG:
					case MET:
					case ILE:
					case LEU:
						return TETRAHEDRAL_CARBON_VDW;
						
					case GLN:
					case GLU:
						if (atomCode.equals("CD")) return TRIGONAL_CARBON_VDW;
						else if (atomCode.equals("CG")) return TETRAHEDRAL_CARBON_VDW;
						
					default:
						System.err.println("Warning: unexpected carbon atom "+atomCode+" for aminoacid "+aa.getThreeLetterCode()+", assigning its standard vdw radius");
						return AtomType.C.getRadius();
				}
			}
			
		// not any of the expected atoms
		} else {
			System.err.println("Warning: unexpected atom "+atomCode+" for aminoacid "+aa.getThreeLetterCode()+", assigning its standard vdw radius");
			return atom.getType().getRadius();
		}
	}

	
	/**
	 * Gets the radius for given nucleotide and atom 
	 * @param nuc
	 * @param atom
	 * @return
	 */
	private static double getRadius(Nucleotide nuc, Atom atom) {
		
		if (atom.getType().equals(AtomType.H)) return AtomType.H.getRadius();
		if (atom.getType().equals(AtomType.D)) return AtomType.D.getRadius();
		
		if (atom.getType()==AtomType.C) return NUC_CARBON_VDW;
		
		if (atom.getType()==AtomType.N) return NUC_NITROGEN_VDW;
		
		if (atom.getType()==AtomType.P) return PHOSPHOROUS_VDW;
		
		if (atom.getType()==AtomType.O) return OXIGEN_VDW;
		
		System.err.println("Warning: unexpected atom "+atom.getCode()+" for nucleotide "+nuc.getTwoLetterCode()+", assigning its standard vdw radius");
		return atom.getType().getRadius();
	}
	
	
	/**
	 * Gets the van der Waals radius of the given atom following the values defined by 
	 * Chothia (1976) J.Mol.Biol.105,1-14
	 * NOTE: the vdw values defined by the paper assume no Hydrogens and thus "inflates" 
	 * slightly the heavy atoms to account for Hydrogens. Thus this method cannot be used
	 * in a structure that contains Hydrogens!
	 * 
	 * If atom is neither part of a nucleotide nor of a standard aminoacid,
	 * the default vdw radius for the element is returned. If atom is of
	 * unknown type (element) the vdw radius of {@link #AtomType().X} is returned
	 * 
	 * @param atom
	 * @return
	 */
	public static double getRadius(Atom atom) {
		
		if (atom.getType()==null) {
			System.err.println("Warning: unrecognised atom "+atom.getCode()+" with serial "+atom.getSerial()+
					", assigning the vdw radius of the default unknown atom (Nitrogen).");
			return AtomType.X.getRadius();
		}
		
		Residue res = atom.getParentResidue();
		
		if (res==null) {
			System.err.println("Warning: unknown parent residue for atom "+atom.getCode()+" with serial "+
					atom.getSerial()+", assigning its default vdw radius");
			return atom.getType().getRadius();
		}
		
		if (res instanceof AaResidue) return getRadius(((AaResidue)res).getAaType(), atom);
		
		if (res instanceof NucResidue) return getRadius(((NucResidue)res).getNucType(), atom);
		
		
		return atom.getType().getRadius();
	}
	
}
