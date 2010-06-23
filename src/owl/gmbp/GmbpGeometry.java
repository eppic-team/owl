package owl.gmbp;

import java.util.HashMap;
import java.util.Set;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import owl.core.structure.AAinfo;
import owl.core.structure.Atom;
import owl.core.structure.Pdb;
import owl.core.structure.Residue;

public class GmbpGeometry {
	
	/* Class that implements rotation and translation invariant framework
	 * All information about the central residue I and the contacting residue J 
	 * have to be handed over to the method getNeighborsTransRotatedCoord.
	 * 
	 * */ 
	 
	
	
	public GmbpGeometry(){}
	
	/**
	 * Extracts and returns the necessary coordinates (of the coordinate system defining atoms: C_alpha, C and N)
	 * @param residue
	 * @return coordinates for all atoms
	 */
	public HashMap<String,Point3d> getTheResidueCoord(Residue res){
		//Method to return the coordinates of specific residue. 
		HashMap<String,Point3d> atom_coord = new HashMap<String, Point3d>();		

		//Test if all necessary atoms CA, C and N available
		if (res.containsAtom("CA") && res.containsAtom("C") && res.containsAtom("N")){
			for(Atom atom:res.getAtoms()){
				atom_coord.put(atom.getCode(), atom.getCoords());
			}			
		}
		
		return atom_coord;
	}
	/**
	 * Extracts and returns the necessary coordinates (of the coordinate system defining atoms: C_alpha, C and N)
	 * @param residue serial (residue number in sequence)
	 * @param residue type
	 * @param pdb
	 * @return coordinates for all atoms
	 */
	public static HashMap<String,Point3d> getTheResidueCoord (int resser, String res_type, Pdb pdb) {
		//Method to return the coordinates of specific residue. 
		
		HashMap<String,Point3d> atom_coord = new HashMap<String, Point3d>();
		
		// Setting up the HashMap "atomcoord" if all the necessary coordinates exist. Declaring as NULL otherwise.  
		if (pdb.hasCoordinates(resser, "CA") && pdb.hasCoordinates(resser, "C") && pdb.hasCoordinates(resser, "N")) {
			
			// Generic code to extract the coordinates from the present amino acid 
			Set<String> atoms = AAinfo.getAtomsForCTAndRes("ALL", res_type);
			for (String atom:atoms) {
				if (pdb.hasCoordinates(resser, atom)) {// Making sure the text-book atom indeed exists! 	
					atom_coord.put(atom, pdb.getAtomCoord(resser, atom));
				}
			}
		}
		else { 
			atom_coord = null;
		}
		
		return atom_coord;  
	}
	
	
	public HashMap<String,Vector3d> getNeighborsTransRotatedCoord (HashMap<String,Point3d> iCoord, HashMap<String,Point3d> jCoord, 
			String iResType, String jResType, Residue iRes, Residue jRes, Boolean toRotate) {
		// METHOD for finding the transformed coordinates of the neighboring residue of a 'node'
		// Doing it according to the EDGE-LIST

		// ONE-TIME Definitions used for the Central Residue as well as the Neighborhood Residues
		HashMap<String,Vector3d> atom_coord_rotated = new HashMap<String, Vector3d>();
		Matrix4d RotMatrix01 = new Matrix4d(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);
		Matrix4d RotMatrix02 = new Matrix4d(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);

		// Get the NEW ORIGIN: c-alpha Coordinates of the central residue
		Vector3d new_origin = new Vector3d (iCoord.get("CA"));//The Origin of the New Neighborhood

		// Get the translated (to c-alpha) coordinates
		atom_coord_rotated = getTranslatedCoord(iRes, iResType, iCoord);
		// CENTRAL RESIDUE: TRANSLATION TO ORIGIN OVER // 		

		if (toRotate) {//Performing the rotations only if the  "toRotate=true"
			RotMatrix01 = getRotationMatrix01 (atom_coord_rotated);
			atom_coord_rotated = getCoordAfterRotation (iRes, iResType, atom_coord_rotated, RotMatrix01);
			// CENTRAL RESIDUE: FIRST ROTATION OVER // 		
	
			RotMatrix02 = getRotationMatrix02 (atom_coord_rotated);
			atom_coord_rotated = getCoordAfterRotation (iRes, iResType, atom_coord_rotated, RotMatrix02);
			// CENTRAL RESIDUE: SECOND ROTATION OVER //
		}

		// FOR EXTRACTING & WRITING THE TRANS-ROTATED COORDINATE INFORMATION OF A *SPECIFIC NEIGHBOR* OF THE *PRESENT NODE*. 
		atom_coord_rotated = getTranslatedCoordNbr (jRes, jResType, jCoord, new_origin); 
		if (toRotate) {//Obtaining the rotated coordinates only if the  "toRotate=true"
			atom_coord_rotated = getCoordAfterRotation(jRes, jResType, atom_coord_rotated, RotMatrix01); 
			atom_coord_rotated = getCoordAfterRotation (jRes, jResType, atom_coord_rotated, RotMatrix02);
		}
		return atom_coord_rotated;
	}
	@SuppressWarnings("unused")
	private static HashMap<String,Vector3d> getNeighborsTransRotatedCoord (int resser_i, String res_type_i, int resser_j, String res_type_j, Pdb pdb, Boolean toRotate) {
		// METHOD for finding the transformed coordinates of the neighboring residue of a 'node'
		// Doing it according to the EDGE-LIST

		// ONE-TIME Definitions used for the Central Residue as well as the Neighborhood Residues
		HashMap<String,Point3d> atom_coord = new HashMap<String, Point3d>(); 
		HashMap<String,Vector3d> atom_coord_rotated = new HashMap<String, Vector3d>();
		Matrix4d RotMatrix01 = new Matrix4d(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);
		Matrix4d RotMatrix02 = new Matrix4d(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);

		// Get the NEW ORIGIN: c-alpha Coordinates of the central residue
		atom_coord = getTheResidueCoord(resser_i, res_type_i, pdb);
		Vector3d new_origin = new Vector3d (atom_coord.get("CA"));//The Origin of the New Neighborhood

		// Get the translated (to c-alpha) coordinates
		atom_coord_rotated = getTranslatedCoord (resser_i, res_type_i, pdb, atom_coord);
		// CENTRAL RESIDUE: TRANSLATION TO ORIGIN OVER // 		

		if (toRotate) {//Performing the rotations only if the  "toRotate=true"
			RotMatrix01 = getRotationMatrix01 (atom_coord_rotated);
			atom_coord_rotated = getCoordAfterRotation (resser_i, res_type_i, pdb, atom_coord_rotated, RotMatrix01);
			// CENTRAL RESIDUE: FIRST ROTATION OVER // 		
	
			RotMatrix02 = getRotationMatrix02 (atom_coord_rotated);
			atom_coord_rotated = getCoordAfterRotation (resser_i, res_type_i, pdb, atom_coord_rotated, RotMatrix02);
			// CENTRAL RESIDUE: SECOND ROTATION OVER //
		}

		// FOR EXTRACTING & WRITING THE TRANS-ROTATED COORDINATE INFORMATION OF A *SPECIFIC NEIGHBOR* OF THE *PRESENT NODE*. 
		atom_coord = getTheResidueCoord(resser_j, res_type_j, pdb);
		atom_coord_rotated = getTranslatedCoordNbr (resser_j, res_type_j, pdb, atom_coord, new_origin); 
		if (toRotate) {//Obtaining the rotated coordinates only if the  "toRotate=true"
		atom_coord_rotated = getCoordAfterRotation (resser_j, res_type_j, pdb, atom_coord_rotated, RotMatrix01);
		atom_coord_rotated = getCoordAfterRotation (resser_j, res_type_j, pdb, atom_coord_rotated, RotMatrix02);
		}
		return atom_coord_rotated;
	}
	
	private static HashMap<String,Vector3d> getTranslatedCoord (Residue res, String res_type, HashMap<String,Point3d> atom_coord) {	
		// Method for getting the coordinates "translated (shifted)" to the new-origin
		// NOTE: The origin is extracted from the available coordinates. Origin defined as the coordinates of CA atom. 
		
		Vector3d new_origin = new Vector3d (atom_coord.get("CA")); //Defining the ORIGIN

		Set<String> atoms = AAinfo.getAtomsForCTAndRes("ALL", res_type);

		// Recomputing the coordinates/vectors with the new origin
		HashMap<String,Vector3d> atom_coord_rotated = new HashMap<String, Vector3d>();
		for (String atom:atoms) {
			if (res.containsAtom(atom)) // Making sure the text-book atom indeed exists! 
			{ 	
				atom_coord_rotated.put(atom, new Vector3d(atom_coord.get(atom)));
				atom_coord_rotated.get(atom).sub(new_origin);
			}
		}

		return atom_coord_rotated;
	}
	private static HashMap<String,Vector3d> getTranslatedCoord (int resser, String res_type, Pdb pdb, HashMap<String,Point3d> atom_coord) {	
		// Method for getting the coordinates "translated (shifted)" to the new-origin
		// NOTE: The origin is extracted from the available coordinates. Origin defined as the coordinates of CA atom. 
		
		Vector3d new_origin = new Vector3d (atom_coord.get("CA")); //Defining the ORIGIN

		Set<String> atoms = AAinfo.getAtomsForCTAndRes("ALL", res_type);

		// Recomputing the coordinates/vectors with the new origin
		HashMap<String,Vector3d> atom_coord_rotated = new HashMap<String, Vector3d>();
		for (String atom:atoms) {
			if (pdb.hasCoordinates(resser, atom)) {// Making sure the text-book atom indeed exists! 	
				atom_coord_rotated.put(atom, new Vector3d(atom_coord.get(atom)));
				atom_coord_rotated.get(atom).sub(new_origin);
			}
		}

		return atom_coord_rotated;
	}	
	
	private static HashMap<String,Vector3d> getTranslatedCoordNbr (Residue res_nbr, String res_type_nbr, HashMap<String,Point3d> atom_coord, Vector3d new_origin) {	
		// Method for getting the coordinates of the neighboring residue, "translated (shifted)" to wrt the central residue (origin).
		// NOTE: The reference origin is passed on as a parameter. 

		Set<String> atoms = AAinfo.getAtomsForCTAndRes("ALL", res_type_nbr);

		// Recomputing the coordinates/vectors with the new origin
		HashMap<String,Vector3d> atom_coord_rotated = new HashMap<String, Vector3d>();
		for (String atom:atoms) {
			if (res_nbr.containsAtom(atom)) {// Making sure the text-book atom indeed exists! 
				atom_coord_rotated.put(atom, new Vector3d(atom_coord.get(atom)));
				atom_coord_rotated.get(atom).sub(new_origin);//Shifting each coordinate of neighbor wrt "C-Alpha of Central Residue"   
			}
		}

		return atom_coord_rotated;
	}	
	private static HashMap<String,Vector3d> getTranslatedCoordNbr (int resser_nbr, String res_type_nbr, Pdb pdb, HashMap<String,Point3d> atom_coord, Vector3d new_origin) {	
		// Method for getting the coordinates of the neighboring residue, "translated (shifted)" to wrt the central residue (origin).
		// NOTE: The reference origin is passed on as a parameter. 

		Set<String> atoms = AAinfo.getAtomsForCTAndRes("ALL", res_type_nbr);

		// Recomputing the coordinates/vectors with the new origin
		HashMap<String,Vector3d> atom_coord_rotated = new HashMap<String, Vector3d>();
		for (String atom:atoms) {
			if (pdb.hasCoordinates(resser_nbr, atom)) {// Making sure the text-book atom indeed exists! 	
				atom_coord_rotated.put(atom, new Vector3d(atom_coord.get(atom)));
				atom_coord_rotated.get(atom).sub(new_origin);//Shifting each coordinate of neighbor wrt "C-Alpha of Central Residue"   
			}
		}

		return atom_coord_rotated;
	}	

	private static HashMap<String,Vector3d> getCoordAfterRotation (Residue res, String res_type, HashMap<String,Vector3d> atom_coord_rotated, Matrix4d RotMatrix_4d) {
		//Method to return coordinates after a certain transformation (as specified by a Matrix4d operator).
		
		Set<String> atoms = AAinfo.getAtomsForCTAndRes("ALL", res_type);	
		//=========== PERFORMING THE FIRST ROTATION WITH theta01 ANGLE ========//
		// Computing the vectors in the rotated frame of reference
		for (String atom:atoms) {
//			if (pdb.hasCoordinates(resser, atom)) {// Making sure the text-book atom indeed exists!
			if (res.containsAtom(atom)) {// Making sure the text-book atom indeed exists! 
				RotMatrix_4d.transform(atom_coord_rotated.get(atom));
			}
		}

		return atom_coord_rotated;
	}
	private static HashMap<String,Vector3d> getCoordAfterRotation (int resser, String res_type, Pdb pdb, HashMap<String,Vector3d> atom_coord_rotated, Matrix4d RotMatrix_4d) {
		//Method to return coordinates after a certain transformation (as specified by a Matrix4d operator).
		
		Set<String> atoms = AAinfo.getAtomsForCTAndRes("ALL", res_type);	
		//=========== PERFORMING THE FIRST ROTATION WITH theta01 ANGLE ========//
		// Computing the vectors in the rotated frame of reference
		for (String atom:atoms) {
			if (pdb.hasCoordinates(resser, atom)) {// Making sure the text-book atom indeed exists! 	
				RotMatrix_4d.transform(atom_coord_rotated.get(atom));
			}
		}

		return atom_coord_rotated;
	}
	
	private static Matrix4d getRotationMatrix01 (HashMap<String,Vector3d> atom_coord_rotated) {

		//Defining the unit vectors for X-, Y- and Z-axis.
		Vector3d iVec = new Vector3d(1,0,0); 

		//================== STEP 01: THE FIRST ROTATION =========================//
		// Finding out the arbitrary vector 'r' around which rotation is to be made. 
		Vector3d rVec = new Vector3d(0,0,0);
		Vector3d cVec = new Vector3d (atom_coord_rotated.get("C"));//Defining ca-c vector
		rVec.cross(cVec, iVec); // Setting r = v (ca-->c) x i
		rVec.normalize(rVec); // Normalizing 'rVec'	

		// Finding the angle by which the FIRST rotation is to be done
		double theta01 = cVec.angle(iVec); // Angle between ca-c vector and the x-axis

		//======== Generating Rotation Matrix Automatically ==========//
		//Setting up the Rotation Axis (rVec)  and the Rotation Angle (theta01)
		AxisAngle4d axisangle01 = new AxisAngle4d();
		axisangle01 = new AxisAngle4d(rVec, theta01);
		//axisangle01 = new AxisAngle4d(0.0, 1.0, 0.0, 0.5236); //Testing (y-axis/30 Degrees)

		// Creating a Matrix4d Element for Rotation Matrix. 
		Matrix4d RotMatrix01_4d = new Matrix4d(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);

		// Generating the ROTATION MATRIX using SetRotation
		RotMatrix01_4d.setRotation(axisangle01);
		//----------- Generating Rotation Matrix Automatically ------------//

		return RotMatrix01_4d; 

	}   

	private static Matrix4d getRotationMatrix02 (HashMap<String,Vector3d> atom_coord_rotated) {

		//Defining the unit vectors for X-, Y- and Z-axis.
		Vector3d iVec = new Vector3d(1,0,0); 
		Vector3d kVec = new Vector3d(0,0,1); 

		//================== STEP 02: THE SECOND ROTATION =========================//				
		Vector3d NVec = new Vector3d (atom_coord_rotated.get("N"));//Defining ca-N vector
		Vector3d wp = new Vector3d(0.0, NVec.y, NVec.z);//The projection of ca-N vector on YZ-plane 
		// NOTE: The projection (wp) needs to be taken on the plane that is perpendicular to the axis 
		// around which one intends to rotate the space eventually.  

		// Finding the angle by which the SECOND rotation is to be done
		// depending on the position of the XYprojection(ca-N) OR wp vector.
		// WARNING: The angle calculation is presently SCREWED UP. Need to fix it!!
		double theta02 = 0.0;
		if (NVec.y > 0.0) {
			theta02 = wp.angle(kVec); // Positive angle made by wp vector with the Z-axis
		} else {
			theta02 = -wp.angle(kVec); // Negative angle made by wp vector with the Z-axis
		}

		//======== Generating Rotation Matrix Automatically ==========//
		//Setting up the Rotation Axis (X-axis, i) and the Rotation Angle (theta02)
		AxisAngle4d axisangle02 = new AxisAngle4d();
		axisangle02 = new AxisAngle4d(iVec, theta02);//Rotation around the X-axis

		// Creating a Matrix4d Element for Rotation Matrix. 
		Matrix4d RotMatrix02_4d = new Matrix4d(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);

		// Generating the ROTATION MATRIX using SetRotation
		RotMatrix02_4d.setRotation(axisangle02);
		//----------- Generating Rotation Matrix Automatically ------------//

		return RotMatrix02_4d;      
	}        
	
	public Vector3d getSphericalFromCartesian (Vector3d Atom_Coord) {
		// METHOD for conversion from Cartesian to Spherical Coordinates
		Vector3d Atom_Coord_sph = new Vector3d(0,0,0);
		
		// Sph(r,theta,phi) or Sph(r,phi,lambda)
		Atom_Coord_sph.x = Math.sqrt( Math.pow(Atom_Coord.x, 2) + Math.pow(Atom_Coord.y, 2) + Math.pow(Atom_Coord.z, 2) ); //r
		//Atom_Coord_sph.y = Math.atan2(Math.sqrt(Math.pow(Atom_Coord.x, 2)+Math.pow(Atom_Coord.y, 2)), Atom_Coord.z); //theta
		Atom_Coord_sph.y = Math.acos( Atom_Coord.z/Atom_Coord_sph.x); //theta
		Atom_Coord_sph.z = Math.atan2(Atom_Coord.y, Atom_Coord.x); //phi

		return Atom_Coord_sph;
	}
	
	public Vector3d getCartesianFromSpherical(Vector3d Atom_Coord){
		// METHOD for conversion from Cartesian to Spherical Coordinates
		Vector3d Atom_Coord_car = new Vector3d(0,0,0);
		
		Atom_Coord_car.x = Atom_Coord.x * Math.sin(Atom_Coord.y) * Math.cos(Atom_Coord.z);
		Atom_Coord_car.y = Atom_Coord.x * Math.sin(Atom_Coord.y) * Math.sin(Atom_Coord.z);
		Atom_Coord_car.x = Atom_Coord.x * Math.cos(Atom_Coord.y);
		
		return Atom_Coord_car;
	}

}
