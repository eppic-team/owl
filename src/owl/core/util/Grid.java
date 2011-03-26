package owl.core.util;

import java.util.Arrays;
import java.util.Map;

import javax.vecmath.Point3d;

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
	
	private GridCell[][][] cells;
	
	private double cutoff;
	private int cellSize;
	
	private Atom[] iAtoms;
	private Atom[] jAtoms;
	
	private int[] bounds;
	
	private boolean noOverlap; // if the 2 sets of atoms are found not to overlap then this is set to true
	
	public Grid(double cutoff) {
		this.cutoff = cutoff;
		this.cellSize = (int) Math.floor(cutoff*SCALE);
		this.noOverlap = false;
	}
	
	private int getFloor(double number) {
		return (cellSize*((int)Math.floor(number*SCALE/cellSize)));
	}
	
	private int xintgrid2xgridindex(int xgridDim) {
		return (xgridDim-bounds[0])/cellSize;
	}
	
	private int yintgrid2ygridindex(int ygridDim) {
		return (ygridDim-bounds[1])/cellSize;
	}
	
	private int zintgrid2zgridindex(int zgridDim) {
		return (zgridDim-bounds[2])/cellSize;
	}
	
	/**
	 * Creates the grid based on the boundaries defined by all atoms given (iAtoms and jAtoms)
	 * and places the atoms in their corresponding grid cells.
	 * @param iAtoms
	 * @param jAtoms
	 */
	public void fillGrid(Atom[] iAtoms, Atom[] jAtoms) {
		this.iAtoms = iAtoms;
		this.jAtoms = jAtoms;
		
		if (!findGridBounds()) {
			noOverlap = true;
			return;
		}
		
		cells = new GridCell[1+(bounds[3]-bounds[0])/cellSize]
		                    [1+(bounds[4]-bounds[1])/cellSize]
		                    [1+(bounds[5]-bounds[2])/cellSize];
		
		int i = 0;
		for (Atom atom:iAtoms) {
			Point3d coord = atom.getCoords();
			int xind = xintgrid2xgridindex(getFloor(coord.x));
			int yind = yintgrid2ygridindex(getFloor(coord.y));
			int zind = zintgrid2zgridindex(getFloor(coord.z));
			if (cells[xind][yind][zind]==null) {
				cells[xind][yind][zind] = new GridCell();
			}
			cells[xind][yind][zind].addIindex(i);
			i++;
		}
		
		int j = 0;
		for (Atom atom:jAtoms) {
			Point3d coord = atom.getCoords();
			int xind = xintgrid2xgridindex(getFloor(coord.x));
			int yind = yintgrid2ygridindex(getFloor(coord.y));
			int zind = zintgrid2zgridindex(getFloor(coord.z));
			if (cells[xind][yind][zind]==null) {
				cells[xind][yind][zind] = new GridCell();
			}
			cells[xind][yind][zind].addJindex(j);
			j++;
		}
		
	}
	
	/**
	 * Calculates an int array of size 6 into member variable bounds:
	 * - elements 0,1,2: minimum x,y,z of the iAtoms and jAtoms
	 * - elements 3,4,5: maximum x,y,z of the iAtoms and jAtoms
	 * Checks also if the i and j grid overlap, i.e. the enclosing bounds
	 * are no more than one cell size apart.
	 * @return true if i and j grids overlap, false if they don't
	 */
	private boolean findGridBounds() {
		bounds = new int[6];
		int[] ibounds = new int[6];
		int[] jbounds = ibounds;
		
		// i
		double[] xs = new double[iAtoms.length];
		double[] ys = new double[iAtoms.length];
		double[] zs = new double[iAtoms.length];
		int c = 0;
		for (Atom atom:iAtoms) {
			xs[c] = atom.getCoords().x;
			ys[c] = atom.getCoords().y;
			zs[c] = atom.getCoords().z;
			c++;
		}
		
		double[] ixminmax = getMinMax(xs);
		double[] iyminmax = getMinMax(ys);
		double[] izminmax = getMinMax(zs);
		ibounds[0] = getFloor(ixminmax[0]);
		ibounds[1] = getFloor(iyminmax[0]);
		ibounds[2] = getFloor(izminmax[0]);
		ibounds[3] = getFloor(ixminmax[1]);
		ibounds[4] = getFloor(iyminmax[1]);
		ibounds[5] = getFloor(izminmax[1]);

		double[] jxminmax = ixminmax;
		double[] jyminmax = iyminmax;
		double[] jzminmax = izminmax;
		
		if (jAtoms!=iAtoms) {
			xs = new double[jAtoms.length];
			ys = new double[jAtoms.length];
			zs = new double[jAtoms.length];
			c = 0;
			for (Atom atom:jAtoms) {
				xs[c] = atom.getCoords().x;
				ys[c] = atom.getCoords().y;
				zs[c] = atom.getCoords().z;
				c++;
			}
			jxminmax = getMinMax(xs);
			jyminmax = getMinMax(ys);
			jzminmax = getMinMax(zs);
			jbounds = new int[6];
			jbounds[0] = getFloor(jxminmax[0]);
			jbounds[1] = getFloor(jyminmax[0]);
			jbounds[2] = getFloor(jzminmax[0]);
			jbounds[3] = getFloor(jxminmax[1]);
			jbounds[4] = getFloor(jyminmax[1]);
			jbounds[5] = getFloor(jzminmax[1]);
		}

		if (!checkGridsOverlap(ibounds,jbounds)) {
			return false;
		}
		
		bounds[0] = Math.min(ibounds[0],jbounds[0]);
		bounds[1] = Math.min(ibounds[1],jbounds[1]);
		bounds[2] = Math.min(ibounds[2],jbounds[2]);
		bounds[3] = Math.max(ibounds[3],jbounds[3]);
		bounds[4] = Math.max(ibounds[4],jbounds[4]);
		bounds[5] = Math.max(ibounds[5],jbounds[5]);
		return true;

	}
	
	/**
	 * Gets an array of size 2 with min and max values contained in given array
	 * @param array
	 * @return
	 */
	private double[] getMinMax(double[] array) {
		double[] minmax = new double[2];
		
		double max = Double.MIN_VALUE;
		double min = Double.MAX_VALUE;

		for(double value : array) {
			if(value > max) max = value;
			if(value < min) min = value; 
		}

		minmax[0] = min;
		minmax[1] = max;
		return minmax;
	}
	
	private class Bound implements Comparable<Bound> {
		int cardinal;
		int value;
		public Bound(int cardinal,int value) {
			this.cardinal = cardinal;
			this.value = value;
		}
		@Override
		public int compareTo(Bound o) {
			return new Integer(this.value).compareTo(o.value);
		}
		public String toString() {
			return "["+cardinal+","+value+"]";
		}
	}
	
	private boolean checkGridsOverlap(int[] ibounds,int[] jbounds) {
		if (ibounds==jbounds) return true;
		// x dimension
		if (!areOverlapping(ibounds[0],ibounds[3],jbounds[0],jbounds[3])) {
			return false;
		}		
		// y dimension
		if (!areOverlapping(ibounds[1],ibounds[4],jbounds[1],jbounds[4])) {
			return false;
		}		
		// z dimension
		if (!areOverlapping(ibounds[2],ibounds[5],jbounds[2],jbounds[5])) {
			return false;
		}		
		return true;
	}
	
	private boolean areOverlapping(int imin, int imax, int jmin, int jmax) {
		Bound[] bounds = {new Bound(0,imin), new Bound(1,imax),
				   		   new Bound(2,jmin), new Bound(3,jmax)};
		Arrays.sort(bounds);
		if ((bounds[0].cardinal==0 && bounds[1].cardinal==1)) {
			if ((bounds[2].value-bounds[1].value)>cellSize) return false;
		} else if (bounds[0].cardinal==2 && bounds[1].cardinal==3) {
			if ((bounds[0].value-bounds[3].value)>cellSize) return false;
		}
		return true;
			
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
		float[][]distMatrix = new float[iAtoms.length][jAtoms.length];
		
		// if the 2 sets of atoms are not overlapping they are too far away and no need to calculate anything
		if (noOverlap) return distMatrix;
		
		for (int xind=0;xind<cells.length;xind++) {
			for (int yind=0;yind<cells[xind].length;yind++) {
				for (int zind=0;zind<cells[xind][yind].length;zind++) {
					// distances of points within this cell
					GridCell thisCell = cells[xind][yind][zind];
					if (thisCell==null) continue;
					thisCell.getDistancesWithinCell(distMatrix,iAtoms,jAtoms,crossed);
					
					// distances of points from this box to all neighbouring boxes: 26 iterations (26 neighbouring boxes)
					for (int x=xind-1;x<=xind+1;x++) {
						for (int y=yind-1;y<=yind+1;y++) {
							for (int z=zind-1;z<=zind+1;z++) {
								if (x==xind && y==yind && z==zind) continue;
								if (x>=0 && x<cells.length && y>=0 && y<cells[x].length && z>=0 && z<cells[x][y].length) {
									if (cells[x][y][z] == null) continue;
									thisCell.getDistancesToNeighborCell(cells[x][y][z],distMatrix,iAtoms,jAtoms,crossed);
								}
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
		
		for (int xind=0;xind<cells.length;xind++) {
			for (int yind=0;yind<cells[xind].length;yind++) {
				for (int zind=0;zind<cells[xind][yind].length;zind++) {
					if (cells[xind][yind][zind]==null) continue;
					int size = getNumGridNbs(xind,yind,zind);	// count number of neighbouring grid cells with points in them
					
					if(densityCount.containsKey(size)) {
						int old = densityCount.get(size);
						densityCount.put(size, ++old);
					} else {
						densityCount.put(size, 1);
					}
				}
			}
		}
		
	}
	
	/** 
	 * Returns the number of neighbours of given grid cell (cells with points in them)
	 * @param xind x index of cell
	 * @param yind y index of cell
	 * @param zind z index of cell
	 * @return 
	 */
	private int getNumGridNbs(int xind,int yind,int zind) {
		int nbs = 0;
		
		for (int x=xind-1;x<=xind+1;x++) {
			for (int y=yind-1;y<=yind+1;y++) {
				for (int z=zind-1;z<=zind+1;z++) {
					if (x==xind && y==yind && z==zind) continue;
					if (x>=0 && x<cells.length && y>=0 && y<cells[x].length && z>=0 && z<cells[x][y].length) {
						if (cells[x][y][z] != null) nbs++;
					}
				}
			}
		}
		return nbs;
	}
	
	public double getCutoff() {
		return cutoff;
	}
}
