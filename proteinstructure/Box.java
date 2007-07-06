package proteinstructure;

import java.util.TreeMap;
import javax.vecmath.Point3d;
import javax.vecmath.Point3i;

public class Box {
	
	Point3i floor;
	
	private TreeMap<Integer,Point3d> pointsInBox;
	
	public Box(Point3i floor){
		this.floor=floor;
		pointsInBox = new TreeMap<Integer, Point3d>(); 
	}
	
	public void putPoint(int serial, Point3d point){
		pointsInBox.put(serial, point);
	}

	public void getDistancesWithinBox(double[][] distMatrix){ //we pass a reference to the distMatrix that we alter within this method
		for (int i_serial:pointsInBox.keySet()) {
			for (int j_serial:pointsInBox.keySet()) {
				if (j_serial>i_serial){
					// this only works if previously we have made sure that atom serials are sequential from 0 to MAXATOMSERIAL
					distMatrix[i_serial][j_serial] = pointsInBox.get(i_serial).distance(pointsInBox.get(j_serial));
				}
			}
		}

	}
	
	public void getDistancesToNeighborBox(Box nbBox ,double[][] distMatrix){
		for (int i_serial:pointsInBox.keySet()){
			for (int j_serial:nbBox.pointsInBox.keySet()){
				if (j_serial>i_serial){
					// this only works if previously we have made sure that atom serials are sequential from 0 to MAXATOMSERIAL
					if (distMatrix[i_serial][j_serial]==0.0){ // i.e. if we haven't passed through this cell yet
						distMatrix[i_serial][j_serial] = pointsInBox.get(i_serial).distance(nbBox.pointsInBox.get(j_serial));
					}
				}
			}
		}
	}
	
}
