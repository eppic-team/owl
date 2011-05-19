package owl.core.structure;

import javax.vecmath.Point3d;

/**
 * A class to store all data that we read from the atom lines of cif/pdbase
 */
public class AtomLine {

	public String labelAsymId; // chainCode
	public String labelAltId;
	public int atomserial;
	public String atom;
	public String element;
	public String res_type;
	public int resSerial;
	public int pdbResSerial;
	public String insCode;
	public Point3d coords;
	public double occupancy;
	public double bfactor;
	public String authAsymId; // pdbChainCode
	public boolean outOfPolyChain;
	public boolean lineIsHetAtm; // true whenever an atom comes from a line starting with HETATM
	
	public boolean isNonPoly;

	public AtomLine(String labelAsymId, String labelAltId, int atomserial, String atom, String element, 
			String res_type, int resSerial, int pdbResSerial, String insCode, Point3d coords, double occupancy, double bfactor, String authAsymId,
			boolean outOfPolyChain, boolean lineIsHetAtm) {
		this.labelAsymId = labelAsymId;
		this.labelAltId = labelAltId;
		this.atomserial = atomserial;
		this.atom = atom;
		this.element = element;
		this.res_type = res_type;
		this.resSerial = resSerial;
		this.pdbResSerial = pdbResSerial;
		this.insCode = insCode;
		this.coords = coords;
		this.occupancy = occupancy;
		this.bfactor = bfactor;
		// this is to support legacy PDB files that can contain a blank in the chain code field, we convert them to "A"
		if (authAsymId.equals(" ")) authAsymId = PdbfileParser.NULL_chainCode;
		this.authAsymId = authAsymId;
		this.outOfPolyChain = outOfPolyChain;
		this.lineIsHetAtm = lineIsHetAtm;
	}

	public String getPdbResSerialWithInsCode() {
		return pdbResSerial+(insCode.equals(".")?"":insCode);
	}
	
	public String toString() {
		return (lineIsHetAtm?"HETATM ":"ATOM ")+labelAsymId+" "+labelAltId+" "+atomserial+" "+atom+" "+element+" "+res_type+" "+
		resSerial+" "+pdbResSerial+" "+(insCode==null?".":insCode)+" "+
		coords+" "+
		String.format("%4.2f",occupancy)+" "+
		String.format("%4.2f",bfactor)+" "+
		authAsymId;
	}
}
