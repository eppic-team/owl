package owl.core.util;

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
	
	public Grid(double cutoff) {
		this.cutoff = cutoff;
		this.cellSize = (int) Math.floor(cutoff*SCALE);
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
		
		bounds = findGridBounds();
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
	 * Returns an int array of size 6:
	 * - elements 0,1,2: minimum x,y,z of the iAtoms and jAtoms
	 * - elements 3,4,5: maximum x,y,z of the iAtoms and jAtoms
	 * @return
	 */
	private int[] findGridBounds() {
		int[] bounds = new int[6];
		int size = iAtoms.length;
		if (jAtoms!=iAtoms) {
			size += jAtoms.length;
		}
		double[] xs = new double[size];
		double[] ys = new double[size];
		double[] zs = new double[size];
		int c = 0;
		for (Atom atom:iAtoms) {
			xs[c] = atom.getCoords().x;
			ys[c] = atom.getCoords().y;
			zs[c] = atom.getCoords().z;
			c++;
		}
		if (jAtoms!=iAtoms) {
			for (Atom atom:jAtoms) {
				xs[c] = atom.getCoords().x;
				ys[c] = atom.getCoords().y;
				zs[c] = atom.getCoords().z;
				c++;
			}
		}
		double[] xminmax = getMinMax(xs);
		double[] yminmax = getMinMax(ys);
		double[] zminmax = getMinMax(zs);
		bounds[0] = getFloor(xminmax[0]);
		bounds[1] = getFloor(yminmax[0]);
		bounds[2] = getFloor(zminmax[0]);
		bounds[3] = getFloor(xminmax[1]);
		bounds[4] = getFloor(yminmax[1]);
		bounds[5] = getFloor(zminmax[1]);
		return bounds;
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
