package owl.core.structure;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3d;

/**
 * Class containing static methods to calculate Accessible Surface Areas based on
 * the rolling ball algorithm by Shrake and Rupley.
 * 
 * The code is taken from a python implementation at http://boscoh.com/protein/asapy
 * and adapted to java.
 * Thanks again to Bosco K. Ho for another great piece of code (our rmsd calculation is also his)
 * 
 * See Shrake, A., and J. A. Rupley. "Environment and Exposure to Solvent of Protein Atoms. 
 * Lysozyme and Insulin." JMB (1973) 79:351-371.
 * 
 * @author duarte_j
 *
 */
public class Asa {

	// Bosco uses as default 960, Shrake and Rupley seem to use in their paper 92
	// to get a bit faster computations I've reduced it to 2/3 of 960, still gets acceptable agreement with default NACCESS accuracy
	private static final int DEFAULT_N_SPHERE_POINTS = 640;
	private static final double DEFAULT_PROBE_SIZE = 1.4;
	
	/**
	 * Returns list of 3d coordinates of points on a sphere using the
	 * Golden Section Spiral algorithm.
	 * @param n the number of points to be used in generating the spherical dot-density
	 * @return
	 */
	private static List<Point3d> generateSpherePoints(int n) {
	    List<Point3d> points = new ArrayList<Point3d>();
	    double inc = Math.PI * (3.0 - Math.sqrt(5.0));
	    double offset = 2.0 / (double)n; 
	    for (int k=0;k<n;k++) {
	        double y = k * offset - 1.0 + (offset / 2.0);
	        double r = Math.sqrt(1.0 - y*y);
	        double phi = k * inc;
	        points.add(new Point3d(Math.cos(phi)*r, y, Math.sin(phi)*r));
	    }
	    return points;
	}

	/**
	 * Returns list of indices of atoms within probe distance to atom k.
	 * @param atoms an array of atoms
	 * @param probe the probe size
	 * @param k index of atom for which we want neighbor indices
	 */
	private static ArrayList<Integer> findNeighborIndices(Atom[] atoms, double probe, int k) {
	    ArrayList<Integer> neighbor_indices = new ArrayList<Integer>();
	    double radius = atoms[k].getRadius() + probe + probe;
	    for (int i=0;i<atoms.length;i++) {
	    	if (i==k) continue;
	        double dist = atoms[i].getCoords().distance(atoms[k].getCoords()); 
	        if (dist < radius + atoms[i].getType().getRadius()) {
	            neighbor_indices.add(i);
	        }
	    }
	    return neighbor_indices;
	}

	/**
	 * Calculates the Accessible Surface Areas of the given atoms
	 * Probe size is default value {@value #DEFAULT_PROBE_SIZE} and number of sphere points is also default value
	 * {@value #DEFAULT_N_SPHERE_POINTS}
	 * @param atoms
	 * @return an array with asa values matching the input atoms array
	 */
	public static double[] calculateAsa(Atom[] atoms) {
		return calculateAsa(atoms, DEFAULT_PROBE_SIZE, DEFAULT_N_SPHERE_POINTS);
	}
	
	/**
	 * Calculates the Accessible Surface Areas of the given atoms, using given probe size.
	 * @param atoms
	 * @param probe the probe size
	 * @param nSpherePoints the number of points to be used in generating the spherical 
	 * dot-density, the more points the more accurate (and slower) calculation.
	 * @return an array with asa values matching the input atoms array
	 */
	public static double[] calculateAsa(Atom[] atoms, double probe, int nSpherePoints) { 
		double[] asas = new double[atoms.length];
	    List<Point3d> sphere_points = generateSpherePoints(nSpherePoints);

	    double cons = 4.0 * Math.PI / (double)sphere_points.size(); 
	    Point3d test_point = new Point3d();

	    for (int i=0;i<atoms.length;i++) {
	    	Atom atom_i = atoms[i];
	    	ArrayList<Integer> neighbor_indices = findNeighborIndices(atoms, probe, i);
	        int n_neighbor = neighbor_indices.size();
	        int j_closest_neighbor = 0;
	        double radius = probe + atom_i.getRadius();

	        int n_accessible_point = 0;
	        
	        for (Point3d point: sphere_points){
	            boolean is_accessible = true;

	            test_point.x = point.x*radius + atom_i.getCoords().x;
	            test_point.y = point.y*radius + atom_i.getCoords().y;
	            test_point.z = point.z*radius + atom_i.getCoords().z;

	            int[] cycled_indices = new int[n_neighbor];
	            int arind = 0;
	            for (int ind=j_closest_neighbor;ind<n_neighbor;ind++) {
	            	cycled_indices[arind] = ind;
	            	arind++;
	            }
	            for (int ind=0;ind<j_closest_neighbor;ind++){
	            	cycled_indices[arind] = ind;
	            	arind++;
	            }

	            for (int j: cycled_indices) {
	                Atom atom_j = atoms[neighbor_indices.get(j)];
	                double r = atom_j.getRadius() + probe;
	                double diff_sq = atom_j.getCoords().distanceSquared(test_point);
	                if (diff_sq < r*r) {
	                    j_closest_neighbor = j;
	                    is_accessible = false;
	                    break;
	                }
	            }
	            if (is_accessible) {
	                n_accessible_point++;
	            }
	        }
	        double area = cons*n_accessible_point*radius*radius ;
	        asas[i] = area;
	        //atom_i.setAsa(area);
	    }
	    return asas;
	}

	/**
	 * To test the class
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		Pdb pdb = new PdbfilePdb(args[0]);
		pdb.load(pdb.getChains()[0]);
		pdb.calcASAs();
		
		double tot = 0;
		for (int resser:pdb.getAllSortedResSerials()) {
			Residue res = pdb.getResidue(resser);
			
			System.out.printf("%3d\t%s\t%6.2f\n",res.getSerial(),res.getAaType().getThreeLetterCode(),res.getAsa());
			tot+=res.getAsa();
		}
		System.out.printf("Total area: %9.2f\n",tot);
		
	}
}
