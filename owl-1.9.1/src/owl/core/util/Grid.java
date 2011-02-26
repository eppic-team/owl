package owl.core.util;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import javax.vecmath.Point3d;
import javax.vecmath.Point3i;

import owl.core.structure.Atom;

/**
 * A grid to be used for calculating atom contacts through our geometric hashing algorithm.
 * 
 * The grid is composed of cells of size of the cutoff so that the distances that need to be calculated
 * are reduced to those within each cell and to the neighbouring cells.
 * 
 * @author duarte_j
 *
 */
public class Grid {
	
	private static final int SCALE=100; // i.e. we use units of hundredths of Amstrongs (thus cutoffs can be specified with a maximum precission of 0.01A)
	
	private HashMap<Point3i,GridCell> cells; // the grid
	
	private double cutoff;
	private int cellSize;
	
	private int[] iAtomSerials; // map from matrix indices to i atom serials
	private int[] jAtomSerials; // map from matrix indices to j atom serials

	public Grid(double cutoff) {
		this.cutoff = cutoff;
		this.cells = new HashMap<Point3i,GridCell>();
		this.cellSize = (int) Math.floor(cutoff*SCALE);
	}
	
	public void putIatoms(Collection<Atom> atoms) {
		this.iAtomSerials = new int[atoms.size()];
		int i = 0;
		for (Atom atom:atoms) {
			Point3d coord = atom.getCoords();
			int floorX = cellSize*((int)Math.floor(coord.x*SCALE/cellSize));
			int floorY = cellSize*((int)Math.floor(coord.y*SCALE/cellSize));
			int floorZ = cellSize*((int)Math.floor(coord.z*SCALE/cellSize));
			Point3i floor = new Point3i(floorX,floorY,floorZ);
			if (cells.containsKey(floor)){
				// we put the coords for atom i in its corresponding box (identified by floor)
				cells.get(floor).put_i_Point(i, coord);
			} else {
				GridCell box = new GridCell(floor);
				box.put_i_Point(i, coord);
				cells.put(floor,box);
			}
			iAtomSerials[i] = atom.getSerial();
			i++;
		}

	}
	
	public void putJatoms(Collection<Atom> atoms) {
		this.jAtomSerials = new int[atoms.size()];
		int j = 0;
		for (Atom atom:atoms) {
			Point3d coord = atom.getCoords();
			int floorX = cellSize*((int)Math.floor(coord.x*SCALE/cellSize));
			int floorY = cellSize*((int)Math.floor(coord.y*SCALE/cellSize));
			int floorZ = cellSize*((int)Math.floor(coord.z*SCALE/cellSize));
			Point3i floor = new Point3i(floorX,floorY,floorZ);
			if (cells.containsKey(floor)){
				// we put the coords for atom i in its corresponding box (identified by floor)
				cells.get(floor).put_j_Point(j, coord);
			} else {
				GridCell box = new GridCell(floor);
				box.put_j_Point(j, coord);
				cells.put(floor,box);
			}
			jAtomSerials[j] = atom.getSerial();
			j++;
		}
		
	}
	
	/**
	 * Calculates a distance matrix for i to j atoms. The distance of any 2 atoms that 
	 * are more than 2 cellSizes apart need not be calculated, in that case we set the 
	 * distance value in the matrix to 0.0f
	 * 
	 * The procedure is first calculate all pairwise distances within the same cell and then
	 * distances of points of each cell to all its neighbouring cells.
	 * @param crossed
	 * @return
	 */
	public float[][] getDistMatrix(boolean crossed) {
		// to minimise memory footprint we use floats
		float[][]distMatrix = new float[iAtomSerials.length][jAtomSerials.length];

		for (Point3i floor:cells.keySet()){ // for each cell
			// distances of points within this cell
			cells.get(floor).getDistancesWithinCell(distMatrix,crossed);

			//TODO should iterate only through half of the neighbours here 
			// distances of points from this box to all neighbouring boxes: 26 iterations (26 neighbouring boxes)
			for (int x=floor.x-cellSize;x<=floor.x+cellSize;x+=cellSize){
				for (int y=floor.y-cellSize;y<=floor.y+cellSize;y+=cellSize){
					for (int z=floor.z-cellSize;z<=floor.z+cellSize;z+=cellSize){
						if (!((x==floor.x)&&(y==floor.y)&&(z==floor.z))) { // skip this box
							Point3i neighbor = new Point3i(x,y,z);
							if (cells.containsKey(neighbor)){
								cells.get(floor).getDistancesToNeighborCell(cells.get(neighbor),distMatrix,crossed);
							}
						}
					}
				}
			} 
		}
		return distMatrix;
	}
	
	public void countDensity(Map<Integer,Integer> densityCount) {
		// count density
		for(Point3i floor:cells.keySet()) {
			//int size = boxes.get(floor).size();
			int size = getNumGridNbs(floor, cellSize);	// count number of neighbouring grid cells with points in them
			if(densityCount.containsKey(size)) {
				int old = densityCount.get(size);
				densityCount.put(size, ++old);
			} else {
				densityCount.put(size, 1);
			}
		}

	}
	
	/** 
	 * Returns the number of neighbours of given grid cell
	 * @param floor the cell's floor
	 * @param boxSize
	 * @return 
	 */
	private int getNumGridNbs(Point3i floor, int boxSize) {
		Point3i neighbor;
		int nbs = 0;
		for (int x=floor.x-boxSize;x<=floor.x+boxSize;x+=boxSize){
			for (int y=floor.y-boxSize;y<=floor.y+boxSize;y+=boxSize){
				for (int z=floor.z-boxSize;z<=floor.z+boxSize;z+=boxSize){
					neighbor = new Point3i(x,y,z);
					if (cells.containsKey(neighbor)) nbs++;
				}
			}
		} 
		// compensate for counting myself as a neighbour
		if(cells.containsKey(floor)) nbs--;
		return nbs;
	}

	
	public int getAtomSerFromIidx(int i) {
		return iAtomSerials[i];
	}
	
	public int getAtomSerFromJidx(int j) {
		return jAtomSerials[j];
	}
	
	public double getCutoff() {
		return cutoff;
	}
}
