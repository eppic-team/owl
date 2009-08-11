package proteinstructure;

import java.io.*;
import java.util.*;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * Some tools for molecular geometry calculations.
 * @author stehr
 *
 */
public class GeometryTools {
	
	/**
	 * Returns the angle between two plains defined by six points.
	 * Plain 1 is defined by points a1,a2,a3 and plain 2 is defined by points b1,b2,b3.
	 * @return the torsion angle in radians
	 */
	public static double torsionAngle(Point3d a1, Point3d a2, Point3d a3, Point3d b1, Point3d b2, Point3d b3) {
		Vector3d vA1 = new Vector3d(a2);
		vA1.sub(a1);
		
		Vector3d vA2 = new Vector3d(a3);
		vA2.sub(a2);
		
		Vector3d vB1 = new Vector3d(b2);
		vB1.sub(b1);
		
		Vector3d vB2 = new Vector3d(b3);
		vB2.sub(b2);
		
		Vector3d normA = new Vector3d();	// norm vector for plane A
		normA.cross(vA1, vA2);
		
		Vector3d normB = new Vector3d();	// norm vector for plane B
		normB.cross(vB1, vB2);
		
		double cos = normA.dot(normB) / (normA.length()*normB.length());	// angle between planes
		
		double sign = Math.signum(normA.dot(vB2)); 				// sign of the angle

		return sign * Math.acos(cos);
	}
	
	/**
	 * Returns the dihedral angle in radians given the coordinates of four consecutive atoms in a polymer.
	 * e.g. for amino acid i
	 * phi: C(i-1), N, C-alpha, C
	 * psi: N, C-alpha, C, N(i+1)
	 * omega: C-alpha, C, N(i+1), C-alpha(i+1)
	 * @return the dihedral angle in radians
	 */
	public static double dihedralAngle(Point3d p1, Point3d p2, Point3d p3, Point3d p4) {
		return torsionAngle(p1,p2,p3,p2,p3,p4);
	}
	
	/**
	 * Read backbone coordinates in pdb format (see http://www.wwpdb.org/documentation/format32/v3.2.html)
	 * @param pdbFile the input file
	 * @param chain the chain for which coordinates are read
	 * @return a map from residue number to Backbone objects holding the atom coordinates
	 * @throws IOException
	 */
	public static TreeMap<Integer, BackboneCoords> readBackbone(File pdbFile, char chain) throws IOException {
		TreeMap<Integer, BackboneCoords> backbone = new TreeMap<Integer, BackboneCoords>();
		BufferedReader in = new BufferedReader(new FileReader(pdbFile));
		String line;
		while((line=in.readLine()) != null) {
			if(!line.startsWith("ATOM")) continue;		// process only atom lines
			if(line.charAt(21) != chain) continue;		// process only given chain
			Integer resNum = new Integer(line.substring(22,26).trim());
			String atomType = line.substring(13, 15).trim();
			double x = Double.parseDouble(line.substring(31, 38).trim());
			double y = Double.parseDouble(line.substring(38, 46).trim());
			double z = Double.parseDouble(line.substring(46, 54).trim());			
			BackboneCoords bb;
			if(backbone.containsKey(resNum)) {
				bb = backbone.get(resNum);
			} else {
				bb = new BackboneCoords();
			}
			if(atomType.equals("CA")) bb.Ca = new Point3d(x,y,z);
                        if(atomType.equals("C")) bb.C = new Point3d(x,y,z);
                        if(atomType.equals("N")) bb.N = new Point3d(x,y,z);
                        if(atomType.equals("O")) bb.O = new Point3d(x,y,z);
			backbone.put(resNum, bb);
		}
		in.close();
		return backbone;
	}
	
	/**
	 * Returns a predicted secondary structure state based on dihedral angles
	 * @param phi the phi angle
	 * @param psi the psi angle
	 * @return one of the states 'H' (helix), 'E' (exteded strand), 'O' (other)
	 */
	public static char classifyByAngle(double phi, double psi) {
		return 'O';
	}
	
	/**
	 * Returns a predicted secondary structure state based on Ca_i-Ca_iplus3 distance
	 * @param distCa the distance between the C-alpha atoms at position i and i+3
	 * @return one of the states 'H' (helix), 'E' (exteded strand), 'O' (other)
	 */
	public static char classifyByDistance(double distCa) {
		return 'O';
	}
	
	/**
	 * Takes the name of a pdb file, reads in the coordinates and prints out the dihedral angles.
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		if(args.length < 1) {
			System.out.println("Usage: " + GeometryTools.class.getName() + " <pdb_file_name>");
			System.exit(1);
		}

		File pdbFile = new File(args[0]);
		char chain = 'A';	// define the chain to be processed
		
		// read coordinates
		TreeMap<Integer,BackboneCoords> myBackbone = readBackbone(pdbFile, chain);

		// calculate angles and predict secondary structure state
		System.out.println("resNum\tphi\tpsi\tomega\tpredSS\tCaDist\tpredSS");
		int maxResNum = myBackbone.lastKey();
		for (int i = 2; i <= maxResNum - 1; i++) {
			BackboneCoords bb_i = myBackbone.get(i);
			BackboneCoords bb_iplus1 = myBackbone.get(i+1);
			BackboneCoords bb_iminus1 = myBackbone.get(i-1);
			double phi = Math.toDegrees(dihedralAngle(bb_iminus1.C, bb_i.N, bb_i.Ca, bb_i.C)); // phi: C(i-1), N, C-alpha, C
			double psi = Math.toDegrees(dihedralAngle(bb_i.N, bb_i.Ca, bb_i.C, bb_iplus1.N));  // psi: N, C-alpha, C, N(i+1)
			double omega = Math.toDegrees(dihedralAngle(bb_i.Ca, bb_i.C, bb_iplus1.N, bb_iplus1.Ca)); // omega: C-alpha, C, N(i+1), C-alpha(i+1)
			char ssPred = classifyByAngle(phi, psi);
			System.out.printf("%3d\t%4.0f\t%4.0f\t%4.0f\t%s", i, phi, psi, omega, ssPred);
			if(myBackbone.containsKey(i+3)) {
				BackboneCoords bb_iplus3 = myBackbone.get(i+3);
				double distCa = bb_i.Ca.distance(bb_iplus3.Ca);
				char ssPred2 = classifyByDistance(distCa);
				System.out.printf("\t%4.1f\t%s", distCa, ssPred2);
			}
			System.out.println();
		}
	}
}
