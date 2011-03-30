package owl.mutanom.core;

import java.io.File;
import java.io.IOException;

import owl.core.structure.AminoAcid;
import owl.core.structure.PdbChain;
import owl.mutanom.output.PyMolScriptMaker;
import owl.mutanom.output.PyMolScriptMaker.Color;



/**
 * A mutation in a protein sequence. Originally only for missense mutations, now extended
 * to also hold basic data about frameshifts and early-stop mutations.
 * @author stehr
 */
public class Mutation {

	/*--------------------------- type definitions --------------------------*/
	public enum MutType {MISSENSE, FRAMESHIFT, EARLYSTOP, INDEL, SILENT};
	
	/*--------------------------- member variables --------------------------*/
	
	MutType type;
	AminoAcid before;
	AminoAcid after;
	int position;
	Gene parent;
	String dnaMutStr;			// use to store the original genomic mutation
	String comment;
	String commentDetails;
	int	mutId;	// e.g. from COSMIC
	
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Initialize a new missense or silent mutation.
	 * TODO: Check that parameters are legal.
	 */
	public Mutation(AminoAcid before, AminoAcid after, int pos) {
		this.before = before;
		this.after = after;
		this.position = pos;
		if(before.equals(after)) {
			this.type = MutType.SILENT;
		} else {
			this.type = MutType.MISSENSE;
		}
	}
	
	/**
	 * Initialize a new mutation of the given type.
	 * This is private because it is possible to pass a wrong type here.
	 */
	private Mutation(AminoAcid before, AminoAcid after, int pos, MutType type) {
		this.before = before;
		this.after = after;
		this.position = pos;
		this.type = type;
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * @returns true if this mutations is classified as a missense mutation, false otherwise
	 */
	public boolean isMissense() {
		return this.type == MutType.MISSENSE;
	}
	
	/**
	 * Sets the parent gene for this mutation
	 */
	public void setParent(Gene g) {
		this.parent = g;
	}

	/**
	 * Sets an optional genomic mutation string for this mutation (for frameshifts etc.)
	 */
	public void setDnaMutStr(String c) {
		this.dnaMutStr = c;
	}
	
	/**
	 * Sets an optional mutation id (e.g. COSMIC identifier)
	 * @param id the mutation id
	 */
	public void setMutId(int id) {
		this.mutId = id;
	}
	
	/**
	 * Sets the comment for this mutation
	 */
	public void setComment(String c) {
		this.comment = c;
	}
	
	/**
	 * Sets the comment details for this mutation
	 */
	public void setCommentDetails(String c) {
		this.commentDetails = c;
	}
	
	/**
	 * @return The comment set for this mutation (can be null)
	 */
	public String getComment() {
		return this.comment;
	}

	/**
	 * @return The detailed comment set for this mutation (can be null)
	 */
	public String getCommentDetails() {
		return this.commentDetails;
	}
	
	/**
	 * @return the (optional) genomic mutation string for this mutation set by the user (can be null)
	 */
	public String getDnaMutStr() {
		return this.dnaMutStr;
	}
	
	/**
	 * @return the mutation type (missesen, early-stop, etc.)
	 */
	public MutType getType() {
		return this.type;
	}
	
	/**
	 * Returns this mutation's position.
	 * @return the position
	 */
	public int getPos() {
		return this.position;
	}
	
	/**
	 * Returns the amino acid before mutation.
	 * @return
	 */
	public AminoAcid getWtAA() {
		return this.before;
	}
	
	/**
	 * Returns the amino acid after mutation.
	 * @return
	 */
	public AminoAcid getMutAA() {
		return this.after;
	}
	
	/**
	 * @return a string representation of this mutation
	 */
	public String toString() {
		return String.format("%s%d%s", this.before.getOneLetterCode(), this.position, this.after.getOneLetterCode());
	}
	
	/**
	 * Returns true iff this mutation is in a known or mutated substructure of its parent gene. 
	 */
	public boolean inKnownOrPredictedStructure() {
		return false;
	}
	
	// Note: Use Gene.getSubstructure(position) instead.
//	/**
//	 * Return the substructure this mutation occurs in or null if no such substructure exists or
//	 * the parent gene of this mutation is not known.
//	 * @return the first substructure of the parent gene covering the position of this mutation
//	 */
//	public Substructure getSubStructure() {
//		if(this.parent != null) {
//			return this.parent.getSubstructure(this.position);
//		} else {
//			return null;
//		}
//	}
	
	/**
	 * Visualize this mutation using PyMol and write a session file and preview image.
	 * The files are written to basename.png (image) and basename.pse (session).
	 * @param baseName the basename of the files to be written.
	 * @param outDir directory where files are being written
	 */
	public void writePngAndPyMolSession(File outDir, String baseName, File pdbFileDir) throws IOException {
		File sessionFile = new File(outDir, baseName + ".pse");
		File pngFile = new File(outDir, baseName + ".png");
		String geneName = this.parent.getGeneName();
		Substructure ss = this.parent.getSubstructure(this.position);
		if(ss == null) {
			System.out.println("No substructure found for " + geneName + " " + this.toString());
		} else {
			String pdbFileName = ss.getDefaultPdbFileName(geneName);
			String objName = geneName + "_" + ss.getRange();
			String selName = String.format("%s%d", this.before.getOneLetterCode(), this.position);
			//String mutSelName = String.format("%d%s", this.position, this.after.getOneLetterCode());			
			//String mutObjName = String.format("%s_%s%d%s", objName, this.before.getOneLetterCode(), this.position, this.after.getOneLetterCode());
			String mutObjName = String.format("%s%d%s", this.before.getOneLetterCode(), this.position, this.after.getOneLetterCode());
			
			PyMolScriptMaker psm = new PyMolScriptMaker(true);
			psm.loadDefaultSettings();
			psm.load(new File(pdbFileDir, pdbFileName), objName);
			if(ss.getChainCode() == null || ss.getChainCode().length() == 0) {
				System.err.println("Error: No chain code for mutation " + this.parent==null?this:this.parent.getGeneName()+this);
			}
			psm.highlightResidue(objName, this.position, ss.getChainCode().charAt(0), Color.RED, selName, null);				
			psm.createObject(mutObjName, selName);
			//psm.copyObject(objName, mutObjName);
			//psm.highlightResidue(mutObjName, this.position, ss.getChainCode().charAt(0), Color.CYAN, mutSelName);
			psm.mutateResidue(mutObjName, this.position, this.after.getThreeLetterCode());
			psm.writePng(pngFile);
			psm.saveSession(sessionFile);
			psm.executeAndClose();
		}
	}
	
	/**
	 * Returns the default basename for files associated with this mutation, e.g. FBK7_R465H or null if parent is unknown.
	 * Warning: Currently only works for point mutations.
	 * @return
	 */
	public String getFileBasename() {
		if(this.parent != null) {
			return this.parent.getGeneName() + "_" + this.before.getOneLetterCode() + this.position + this.after.getOneLetterCode();
		} else return null;
	}
	
	/*---------------------------- static methods ---------------------------*/
	
	/**
	 * Returns a mutation object given the textual notation of the mutation, e.g. K23G
	 * or null if an error occured.
	 */
	public static Mutation parseMutation(String s) {
		String s2 = s.trim();
		AminoAcid before;
		AminoAcid after;
		MutType type = MutType.MISSENSE;
		int pos;
		// workaround for frameshifts, early stops, etc: if mutation string starts with a number, take that number as the position of the stop/fs/etc.
		if(s2.indexOf("f") > 0) {
			// frame shift
			pos = Integer.parseInt(s2.substring(0, s2.indexOf("f")));
			before = AminoAcid.XXX;
			after = AminoAcid.STP;
			type = MutType.FRAMESHIFT;
		} else
			if(s2.indexOf("*") > 0) {
				// early stop
				pos = Integer.parseInt(s2.substring(0, s2.indexOf("*")));
				before = AminoAcid.XXX;
				after = AminoAcid.XXX;
				type = MutType.EARLYSTOP;
			} else		
				if(s2.indexOf("_") > 0) {
					// insertion/deletion
					pos = Integer.parseInt(s2.substring(0, s2.indexOf("_")));
					before = AminoAcid.XXX;
					after = AminoAcid.XXX;
					type = MutType.INDEL;
				} else {
					if(s2.length() < 3) return null;
					before = AminoAcid.getByOneLetterCode(s2.charAt(0));
					after = AminoAcid.getByOneLetterCode(s2.charAt(s2.length()-1));
					pos = Integer.parseInt(s2.substring(1, s2.length()-1));
				}
		return new Mutation(before, after, pos, type);
	}
	
	/**
	 * Returns a mutation object given a pdb object and a mutated position,
	 * or null if the given position can not be determined.
	 * Target amino acid will be marked as unknown ('X');
	 * @param pos
	 * @param pdb
	 * @return
	 */
	public static Mutation fromPos(int pos, PdbChain pdb) {
		if(pos < 0 || pos > pdb.getFullLength()) return null;
		return new Mutation(pdb.getResidue(pos).getAaType(),AminoAcid.XXX,pos);
	}
	
	
}
