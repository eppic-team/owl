package owl.core.util;

import java.util.ArrayList;
import javax.vecmath.Point3i;

import owl.core.structure.Atom;

/**
 * A grid cell to be used in contact calculation via our geometric hashing algorithm.
 * 
 * @author duarte_j
 *
 */
public class GridCell {
	
	
	Point3i floor; // the floor of the cell, which identifies it
	
	
	private ArrayList<Integer> iIndices;
	private ArrayList<Integer> jIndices;
	
	public GridCell(Point3i floor){
		this.floor=floor;
		iIndices = new ArrayList<Integer>();
		jIndices = new ArrayList<Integer>();
	}
	
	public void addIindex(int serial){
		iIndices.add(serial);
	}

	public void addJindex(int serial){
		jIndices.add(serial);
	}

	public void getDistancesWithinCell(float[][] distMatrix, Atom[] iAtoms, Atom[] jAtoms, boolean crossed){
		for (int i:iIndices) {
			for (int j:jIndices) {
				if (!crossed) {
					if (j>i) { 
						// this only works if previously we have made sure that atom serials are sequential from 0 to MAXATOMSERIAL
						distMatrix[i][j] = (float)iAtoms[i].getCoords().distance(jAtoms[j].getCoords());
					} 
				} else {
					// We could check if two atoms in the same cell have the same serial.
					// This could happen when the 2 contact types (i/j) have overlapping atoms, e.g. ALL/BB
					// The check is not strictly necessary, because distance in case i=j would be 0 (atom to itself). 
					// It's just to make sure that there wouldn't be rounding problems in comparing to 0.0 in getAIgraph in Pdb
					distMatrix[i][j] = (float)iAtoms[i].getCoords().distance(jAtoms[j].getCoords());
				}
			}
		}

	}
	
	public void getDistancesToNeighborCell(GridCell nbBox ,float[][] distMatrix, Atom[] iAtoms, Atom[] jAtoms, boolean crossed){
		for (int i:iIndices) {
			for (int j:nbBox.jIndices) {
				if (!crossed) {
					if (j>i) {
						// this only works if previously we have made sure that atom serials are sequential from 0 to MAXATOMSERIAL
						if (distMatrix[i][j]==0.0f){ // i.e. if we haven't passed through this cell yet
							distMatrix[i][j] = (float)iAtoms[i].getCoords().distance(jAtoms[j].getCoords());
						}
					}
				} else {
					distMatrix[i][j] = (float)iAtoms[i].getCoords().distance(jAtoms[j].getCoords());
				}
			}
		}
	}
	
}
