package proteinstructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Locale;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

/**
 * Class representing a polymer's 3-dimensional structure
 * The polymer can be composed by multiple chains of polypeptide/DNA/RNA...
 * @author duarte
 *
 */
public class Polymer implements Iterable<PolResidue> {

	private static final double DEFAULT_B_FACTOR = 0.00;		// default value if no b-factor is given
	private static final double DEFAULT_OCCUPANCY = 1.00;		// default value if no occupancy is given

	/**
	 * A class to represent a moment of inertia together with its associated axis.
	 * Implements Comparable, being the comparison based on the values of the moment of 
	 * inertia.
	 * @author duarte
	 *
	 */
	private class InertiaMomentAndAxis implements Comparable<InertiaMomentAndAxis> {
		public double moment;
		public Vector3d axis;
		public InertiaMomentAndAxis(double moment, Vector3d axis) {
			this.moment = moment;
			this.axis = axis;
		}
		public int compareTo(InertiaMomentAndAxis o) {
			if (this.moment<o.moment) return -1;
			if (this.moment>o.moment) return 1;
			return 0;
		}
	}
	
	/**
	 * The Map containing all residues. The keys are of the form 
	 * chainCode+residueSerial, e.g. A123 for residue 123 of chain A
	 */
	private TreeMap<String,PolResidue> residues;


	/**
	 * Creates a new Polymer with no residues. Use {@link #addResidue(PolResidue)}
	 * to add residues subsequently.
	 */
	public Polymer() {
		residues = new TreeMap<String, PolResidue>();
	}

	/**
	 * Reads a multi chain molecule from a file in PDB format.
	 * @param pdbFile
	 * @param chains chains that will be read from the PDB file
	 * @throws IOException
	 */
	public void readFromPdbFile(File pdbFile, String[] chains) throws IOException {

		String chainCodeRegex = "[";
		for (int i=0;i<chains.length;i++) {
			chainCodeRegex+=chains[i];
		}
		chainCodeRegex+="]";
		
		String atomTypeRegex = "CA |P  "; // for protein we take the CA, for DNA the P
		
		BufferedReader fpdb = new BufferedReader(new FileReader(pdbFile));
		int linecount=0;
		String line;

		while((line = fpdb.readLine()) != null ) {
			linecount++;
			if (line.startsWith("ATOM")) {
				Pattern pl = Pattern.compile("^.{6}(.....).{2}("+atomTypeRegex+").{1}(...).{1}("+chainCodeRegex+")(.{4})(.).{3}(.{8})(.{8})(.{8})",Pattern.CASE_INSENSITIVE);
				Matcher ml = pl.matcher(line);
				if (ml.find()) {
					String atomType = ml.group(2).trim();
					String resType = ml.group(3).trim();
					String chainCode = ml.group(4);
					int resSerial = Integer.parseInt(ml.group(5).trim());
					double x = Double.parseDouble(ml.group(7).trim());
					double y = Double.parseDouble(ml.group(8).trim());
					double z = Double.parseDouble(ml.group(9).trim());
					Point3d coords = new Point3d(x,y,z);
					this.addResidue(new PolResidue(resType,resSerial,chainCode,chainCode,atomType,coords));
				}
			}		
		}
		fpdb.close();
	}
	
	/**
	 * Adds a residue to this Polymer
	 * @param residue
	 */
	public void addResidue(PolResidue residue) {
		residues.put(residue.getChainCode()+residue.getSerial(),residue);
	}
	
	/**
	 * Gets the Residue given a residue serial and chain code 
	 * @param serial
	 * @param chainCode
	 * @return
	 */
	public PolResidue getResidue(int serial, String chainCode) {
		return residues.get(chainCode+serial);
	}

	public Iterator<PolResidue> iterator() {
		return residues.values().iterator();
	}

	/**
	 * Gets the center of mass 
	 * @return
	 */
	public Point3d getCenterOfMass() {
		double massSum = 0;
		Point3d sum = new Point3d();
		for (PolResidue residue:this) {
			Point3d coord = new Point3d(residue.getCoords());
			coord.scale(residue.getMass());
			sum.add(coord);
			massSum+=residue.getMass();
		}
		sum.scale(1.0/massSum);
		return sum;
	}
	
	/**
	 * Removes given residues from this Polymer 
	 * @param residues
	 */
	public void removeResidues(Set<String> residues) {
		Iterator<PolResidue> it = this.iterator();
		while (it.hasNext()){
			PolResidue residue = it.next();
			if (residues.contains(residue.getChainCode()+residue.getSerial())) {
				it.remove();
			}
		}
	}
	
	/**
	 * Gets all residues of given type.
	 * @param type
	 * @return
	 */
	public ArrayList<PolResidue> getResiduesOfType(String type) {
		ArrayList<PolResidue> residues = new ArrayList<PolResidue>();
		for (PolResidue residue:this) {
			if (residue.getType().equals(type)) {
				residues.add(residue);
			}
		}
		return residues;
	}
	
	/**
	 * Returns the number of residues in this polymer
	 * @return
	 */
	public int size() {
		return residues.size();
	}
	
	/**
	 * Gets all coordinates of all atoms of this polymer in an array.
	 * The order of the atoms in the array is arbitrary. 
	 * @return
	 */
	private Point3d[] getCoordsArray() {
		Point3d[] coordsArray = new Point3d[this.size()];
		int i=0;
		for (PolResidue residue:this) {
			coordsArray[i] = residue.getCoords();
			i++;
		}
		return coordsArray;
	}
	
	/**
	 * Gets the distance matrix (only upper half filled) with pairwise distances
	 * between all atoms of this polymer.
	 * The indices of the atoms are as in {@link #getCoordsArray()} 
	 * @return
	 */
	public double[][] getDistanceMatrix() {
		Point3d[] coordsArray = getCoordsArray();
		double[][] distMat = new double[this.size()][this.size()];
		for (int i=0;i<this.size();i++) {
			for (int j=i+1;j<this.size();j++) {
				distMat[i][j] = coordsArray[i].distance(coordsArray[j]);
			}
		}
		return distMat;
	}
	
	/**
	 * Gets the maximum distance of any two atoms in this polymer, i.e.
	 * some sort of diameter of gyration.
	 * @return
	 */
	public double getDiameter() {
		double[][] distMat = getDistanceMatrix();
		ArrayList<Double> distances = new ArrayList<Double>();
		for (int i=0;i<this.size();i++) {
			for (int j=i+1;j<this.size();j++) {
				distances.add(distMat[i][j]);
			}
		}
		return Collections.max(distances);
	}
	
	/**
	 * Gets the momemt of inertia tensor for the given center point.
	 * @param center
	 * @return a 3x3 array with containing the moments of inertia tensor
	 */
	public double[][] getMomentInertiaTensor(Point3d center) {
		double[][] tensor = new double[3][3];
		for (PolResidue residue:this) {
			Point3d coords = new Point3d(residue.getCoords());
			coords.sub(center);
			double m = residue.getMass();
			tensor[0][0] += m*(coords.y*coords.y+coords.z*coords.z); 
			tensor[1][1] += m*(coords.x*coords.x+coords.z*coords.z);
			tensor[2][2] += m*(coords.x*coords.x+coords.y*coords.y);
			tensor[0][1] += -m*coords.x*coords.y;
			tensor[0][2] += -m*coords.x*coords.z;
			tensor[1][2] += -m*coords.y*coords.z;
		}
		tensor[1][0] = tensor[0][1];
		tensor[2][0] = tensor[0][2];
		tensor[2][1] = tensor[1][2];
		return tensor;
	}
	
	/**
	 * Gets the principal moments of inertia, i.e. the result of diagonalizing the
	 * moments of inertia tensor.
	 * @param center
	 * @return
	 */
	public ArrayList<InertiaMomentAndAxis> getPrincipalMomentsInertia(Point3d center) {
		Matrix tensor = new Matrix(getMomentInertiaTensor(center));
		EigenvalueDecomposition eig = tensor.eig();
		Matrix eigVals = eig.getD();
		Matrix eigVecs = eig.getV();
		//eigVals.print(12, 3);
		//eigVecs.print(12, 3);
		ArrayList<InertiaMomentAndAxis> momentsAndAxes = new ArrayList<InertiaMomentAndAxis>();
		for (int i=0;i<3;i++) {
			momentsAndAxes.add(new InertiaMomentAndAxis(eigVals.get(i, i), new Vector3d(eigVecs.get(0, i), eigVecs.get(1, i), eigVecs.get(2, i))));
		}
		return momentsAndAxes;
	}
	
	/**
	 * Gets the inertia axis corresponding to the maximum principal moment of 
	 * inertia.
	 * @param center
	 * @return
	 */
	public Vector3d getBiggestInertiaAxis(Point3d center) {
		ArrayList<InertiaMomentAndAxis> momentsAndAxes = getPrincipalMomentsInertia(center);
		Collections.sort(momentsAndAxes); // sorts in ascending order: biggest value is last
		return momentsAndAxes.get(2).axis;
	}
	
	/**
	 * Transform the coordinates of this 3-D structure translating them to the given center
	 * and rotating them so that the given axis aligns with the z-axis
	 * @param center
	 * @param axis
	 */
	public void transformRefFrameToCenterAndAxis(Point3d center, Vector3d axis) {
		// finding the rotation matrix to align z axis to the given inertia axis
		Vector3d r = new Vector3d();
		Vector3d k = new Vector3d(0,0,1);
		r.cross(axis, k); // this is the axis of rotation
		double alpha = axis.angle(k); // this is the angle to rotate
		AxisAngle4d axisAngle = new AxisAngle4d(r, alpha);
		// note that the matrix needs to be initialised to the unit matrix otherwise setRotation() doesn't work properly
		Matrix4d rot = new Matrix4d(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1); 
		rot.setRotation(axisAngle);
		for (PolResidue residue:this) { 
			Point3d coords = residue.getCoords();
			// translate to new origin
			coords.sub(center);
			// rotate so that z axis is the given axis
			rot.transform(coords);			
		}	
	}
	
	/**
	 * Perform a rotate of this polymer around rotAxis with the given rotAngle 
	 * @param rotAxis the vector around which the rotation will be performed
	 * @param rotAngle the rotation angle in radians
	 */
	public void rotate(Vector3d rotAxis, double rotAngle) {
		AxisAngle4d axisAngle = new AxisAngle4d(rotAxis, rotAngle);
		// note that the matrix needs to be initialised to the unit matrix otherwise setRotation() doesn't work properly
		Matrix4d rot = new Matrix4d(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1); 
		rot.setRotation(axisAngle);
		for (PolResidue residue:this) { 
			Point3d coords = residue.getCoords();
			// rotate
			rot.transform(coords);			
		}			
	}
	
	/**
	 * Translate this Polymer to given center
	 * @param center
	 */
	public void translate(Point3d center) {
		for (PolResidue residue:this) { 
			Point3d coords = residue.getCoords();
			// translate to new origin
			coords.sub(center);
		}	
	}

	/**
	 * Write this polymer in PDB file format.
	 * The order of atoms is arbitrary!
	 * @param file
	 * @throws FileNotFoundException
	 */
	public void writeToPdbFile(File file) throws FileNotFoundException {
		PrintStream out = new PrintStream(new FileOutputStream(file));
		writeAtomLines(out);
		out.close();
	}
	
	public void writeAtomLines(PrintStream out) {
		int atomser = 1;
		for (PolResidue residue:this) {
			String atom = residue.getAtomType();
			String res = residue.getType();
			String chainCodeStr = residue.getChainCode();
			int resser = residue.getSerial();
			Point3d coords = residue.getCoords();
			double occupancy = DEFAULT_OCCUPANCY;
			double bFactor = DEFAULT_B_FACTOR;
			String atomType = "C";
			
			// Local.US is necessary, otherwise java prints the doubles locale-dependant (i.e. with ',' for some locales)
			out.printf(Locale.US,"ATOM  %5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %s\n",
						atomser, atom, res, chainCodeStr, resser, coords.x, coords.y, coords.z, occupancy, bFactor, atomType);
			atomser++;
		}

	}
	
	/**
	 * Deep copies this Polymer
	 * @return
	 */
	public Polymer copy() {
		Polymer mol = new Polymer();
		for (PolResidue residue:this) {
			mol.addResidue(residue.copy());
		}
		return mol;
	}
	
}
