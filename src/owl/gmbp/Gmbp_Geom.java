package owl.gmbp;

import java.util.HashMap;

import javax.vecmath.Point3d;

import owl.core.structure.Residue;

public class Gmbp_Geom {
	
//	private Residue iRes;
//	private Residue jRes;
	
	/**
	 * Extracts and returns the coordinates of certain atom (CA,C,N,...) from residue
	 * @param residue 
	 * @param atom type
	 * @return coordinates certain atom
	 */
	public static Point3d getRelativeCoord (Residue iRes, String atomType){
		Point3d relativeCoord = new Point3d(0,0,0);
		//TODO: compute and return the relative coordinate of the atom of type atomType
		// coord of CAatom of iRes should be (0,0,0) ...
		
		//Test if all necessary atoms CA, C and N available
		if (iRes.containsAtom("CA") && iRes.containsAtom("C") && iRes.containsAtom("N")){
			
		}
		return relativeCoord;
	}
	/**
	 * Extracts and returns the coordinates of all atoms (CA,C,N,...) from residue
	 * @param residue 
	 * @return coordinates all atoms of residue
	 */
	public static HashMap<String,Point3d> getTheResidueCoord (Residue iRes) {
		HashMap<String,Point3d> atom_coord = new HashMap<String, Point3d>();
		//TODO: compute and return the relative coordinate of all available atoms of the residue
		// coord of CAatom of iRes should be (0,0,0) ...

		//Test if all necessary atoms CA, C and N available
		if (iRes.containsAtom("CA") && iRes.containsAtom("C") && iRes.containsAtom("N")){
			
		}
		return atom_coord;
	}
	/**
	 * Extracts and returns the relative coordinates of certain atom (CA,C,N,...) 
	 * for residue jRes with respect to iRes
	 * @param central residue
	 * @param residue to rotate (with respect to iRes)
	 * @param atom type
	 * @return coordinates certain atom
	 */
	public static Point3d getRelativeCoord (Residue iRes, Residue jRes, String atomType){
		Point3d relativeCoord = new Point3d(0,0,0);
		//TODO: compute and return the relative coordinate of the atom of type atomType
		
		//Test if all necessary atoms CA, C and N available
		if (iRes.containsAtom("CA") && iRes.containsAtom("C") && iRes.containsAtom("N")){
			
		}
		return relativeCoord;
	}
	/**
	 * Extracts and returns the relative coordinates of all atoms (CA,C,N,...) 
	 * for residue jRes with respect to iRes
	 * @param central residue
	 * @param residue to rotate (with respect to iRes)
	 * @return coordinates all atoms of residue
	 */
	public static HashMap<String,Point3d> getTheResidueCoord (Residue iRes, Residue jRes) {
		HashMap<String,Point3d> atom_coord = new HashMap<String, Point3d>();
		//TODO: compute and return the relative coordinate of all available atoms of the residue jRes

		//Test if all necessary atoms CA, C and N available
		if (iRes.containsAtom("CA") && iRes.containsAtom("C") && iRes.containsAtom("N")){
			
		}
		return atom_coord;
	}
	/**
	 * Extracts and returns the relative coordinates of the coordinates with respect to iRes
	 * @param central residue
	 * @param coordinates to transform with respect to iRes
	 * @return coordinates all atoms of residue
	 */
	public static HashMap<String,Point3d> getTheResidueCoord (Residue iRes, Point3d coord) {
		HashMap<String,Point3d> atom_coord = new HashMap<String, Point3d>();
		//TODO: compute and return the relative coordinate (based on coord) with respect to iRes

		//Test if all necessary atoms CA, C and N available
		if (iRes.containsAtom("CA") && iRes.containsAtom("C") && iRes.containsAtom("N")){
			
		}
		return atom_coord;
	}
	/**
	 * Extracts and returns the relative coordinates of the coordinates (coord) with respect to iRes
	 * @param CA coordinate of the central residue
	 * @param C coordinate of the central residue
	 * @param N coordinate of the central residue
	 * @param coordinates to transform with respect to iRes
	 * @return coordinates all atoms of residue
	 */
	public static HashMap<String,Point3d> getTheResidueCoord (Point3d coord_CA, Point3d coord_C, Point3d coord_N, Point3d coord) {
		HashMap<String,Point3d> atom_coord = new HashMap<String, Point3d>();
		//TODO: compute and return the relative coordinate (based on coord) with respect to iRes

		
		return atom_coord;
	}
	

}
