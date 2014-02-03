package owl.core.structure;

import java.io.File;
import java.util.ArrayList;

import javax.vecmath.Point3d;


/**
 * Class containing static methods to calculate Accessible Surface Areas based on
 * the rolling ball algorithm by Shrake and Rupley.
 * 
 * The code is taken from a python implementation at http://boscoh.com/protein/asapy
 * and adapted to java.
 * Thanks again to Bosco K. Ho for another great piece of code (our rmsd calculation is also his)
 * 
 * See 
 * Shrake, A., and J. A. Rupley. "Environment and Exposure to Solvent of Protein Atoms. 
 * Lysozyme and Insulin." JMB (1973) 79:351-371.
 * Lee, B., and Richards, F.M. "The interpretation of Protein Structures: Estimation of
 * Static Accessibility" JMB (1971) 55:379-400
 * @author duarte_j
 *
 */
public class Asa {

	// Bosco uses as default 960, Shrake and Rupley seem to use in their paper 92 (not sure if this is actually the same parameter)
	public static final int DEFAULT_N_SPHERE_POINTS = 960;
	public static final double DEFAULT_PROBE_SIZE = 1.4;
	public static final int DEFAULT_NTHREADS = 1;
	
	private class GroupASACalcThread extends Thread {
		
		int start;
		int end;
		Atom[] atoms;
		Point3d[] sphere_points;
		double[] asas;
		double probe;
		double cons;
		
		public GroupASACalcThread(int start, int end, Atom[] atoms, Point3d[] sphere_points, double[] asas, double probe, double cons) {
			this.start = start;
			this.end = end;
			this.atoms = atoms;
			this.sphere_points = sphere_points;
			this.asas = asas;
			this.probe = probe;
			this.cons = cons;
		}

		public void run() {
			calcGroupOfAsas(start, end, atoms, sphere_points, asas, probe, cons);
		}
	}
	
	
	/**
	 * Returns list of 3d coordinates of points on a sphere using the
	 * Golden Section Spiral algorithm.
	 * @param n the number of points to be used in generating the spherical dot-density
	 * @return
	 */
	private static Point3d[] generateSpherePoints(int n) {
	    Point3d[] points = new Point3d[n];
	    double inc = Math.PI * (3.0 - Math.sqrt(5.0));
	    double offset = 2.0 / (double)n; 
	    for (int k=0;k<n;k++) {
	        double y = k * offset - 1.0 + (offset / 2.0);
	        double r = Math.sqrt(1.0 - y*y);
	        double phi = k * inc;
	        points[k] = new Point3d(Math.cos(phi)*r, y, Math.sin(phi)*r);
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
		// looking at a typical protein case, number of neighbours are from ~10 to ~50, with an average of ~30
		// Thus 40 seems to be a good compromise for the starting capacity
	    ArrayList<Integer> neighbor_indices = new ArrayList<Integer>(40);
	    
	    double radius = AtomRadii.getRadius(atoms[k]) + probe + probe;
	    
	    for (int i=0;i<atoms.length;i++) {
	    	if (i==k) continue;
	    	
	        double dist = atoms[i].getCoords().distance(atoms[k].getCoords());
	    	
	        if (dist < radius + AtomRadii.getRadius(atoms[i])) {
	            neighbor_indices.add(i);
	        }
	    }
	    
	    return neighbor_indices;
	}

	/**
	 * Calculates the Accessible Surface Areas of the given atoms, using given probe size.
	 * @param atoms
	 * @param probe the probe size
	 * @param nSpherePoints the number of points to be used in generating the spherical 
	 * dot-density, the more points the more accurate (and slower) calculation
	 * @param nThreads the number of parallel threads to use for the calculation
	 * @return an array with asa values matching the input atoms array
	 */
	public static double[] calculateAsa(Atom[] atoms, double probe, int nSpherePoints, int nThreads) { 
		
		double[] asas = new double[atoms.length];
	    Point3d[] sphere_points = generateSpherePoints(nSpherePoints);

	    double cons = 4.0 * Math.PI / (double)nSpherePoints; 

	    if (nThreads==1) {
		    for (int i=0;i<atoms.length;i++) {	    	
		        asas[i] = calcSingleAsa(atoms, sphere_points, i, probe, cons); 
		        //atom_i.setAsa(area);
		    }
	    } else {
	    	// NOTE the multithreaded calculation does not scale up well (4 CPUs ~ x2.8, 8CPUs ~ x2.9)
	    	// tried copying the arrays (atoms and sphere_points) as new arrays but the scaling behaves the same
	    	// also tried dividing the asas array in parts and them joining all together but same scaling behaviour
	    	// I guess it's simply memory bottlenecks of the architecture :(
	    	GroupASACalcThread[] threads = new GroupASACalcThread[nThreads];
	    	
	    	int[] startIndices = getStartingIdxForGroups(atoms.length, nThreads);


		    for (int k=0;k<nThreads;k++) {
		    	threads[k] = new Asa().new GroupASACalcThread(startIndices[k], startIndices[k+1], atoms, sphere_points, asas, probe, cons);
		    	threads[k].start();
		    }
	    	
		    for (int k=0;k<nThreads;k++) {
		    	try {
		    		threads[k].join();
		    	} catch (InterruptedException e) {
		    		System.err.println("Unexpected error while running multi-threaded ASA calculation. Exiting.");
		    		e.printStackTrace();
		    		System.exit(1);
		    	}
		    }
	    }
	    
	    return asas;
	}
	
	private static void calcGroupOfAsas(int startIdx, int endIdx, Atom[] atoms, Point3d[] sphere_points, double[] asas, double probe, double cons) {
		for (int i=startIdx;i<endIdx;i++) {
			asas[i] = calcSingleAsa(atoms, sphere_points, i, probe, cons);
		}
	}

	private static double calcSingleAsa(Atom[] atoms, Point3d[] sphere_points, int i, double probe, double cons) {
    	Atom atom_i = atoms[i];
    	ArrayList<Integer> neighbor_indices = findNeighborIndices(atoms, probe, i);
        int n_neighbor = neighbor_indices.size();
        int j_closest_neighbor = 0;
        double radius = probe + AtomRadii.getRadius(atom_i);

        int n_accessible_point = 0;
        
        for (Point3d point: sphere_points){
            boolean is_accessible = true;
            Point3d test_point = new Point3d(point.x*radius + atom_i.getCoords().x,
            								point.y*radius + atom_i.getCoords().y,
            								point.z*radius + atom_i.getCoords().z);

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
                double r = AtomRadii.getRadius(atom_j) + probe;
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
        return cons*n_accessible_point*radius*radius;
	}
	
	private static int[] getStartingIdxForGroups(int n, int nGroups) {
		int[] indices = new int[nGroups+1];

		int baseSize = n/nGroups;
		int remainder = n%nGroups;

		indices[0] = 0;
		for (int k=1;k<nGroups;k++){
			indices[k] = indices[k-1]+baseSize;
			if (k<remainder) {
				indices[k]+=1;
			}
		}
		indices[nGroups] = n;
		return indices;
	}

	/**
	 * To test the class
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		PdbAsymUnit pdb = new PdbAsymUnit(new File(args[0]));
		int nThreads = Integer.parseInt(args[1]);
		
		pdb.removeHatoms();
		
		long start = System.currentTimeMillis();
		pdb.calcASAs(1000,nThreads,false);
		long end = System.currentTimeMillis();

		
		double tot = 0;
		
		for (PdbChain chain:pdb.getAllChains()) {
			for (Residue res:chain) {

				System.out.printf("%3d\t%s\t%6.2f\n",res.getSerial(),res.getLongCode(),res.getAsa());
				tot+=res.getAsa();
			}
		}
		System.out.printf("Total area: %9.2f\n",tot);
		System.out.printf("Time: %4.1fs\n",((end-start)/1000.0));
		

	}
}
